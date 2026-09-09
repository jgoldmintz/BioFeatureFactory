"""Isolated adabmDCA subprocess attempts and verified completion records."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from biofeaturefactory.mutation_effects.bin.adabmdca_task import (
    THREAD_VARIABLES,
    _gpu_free_memory,
    lease_gpu,
    run_task,
    verify_completion,
)


FAKE_BACKEND = '''
import argparse
import json
import os
import sys
import time
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--gene")
parser.add_argument("--fasta")
parser.add_argument("--mutations")
parser.add_argument("--msa")
parser.add_argument("--codon-msa")
parser.add_argument("--protein-params")
parser.add_argument("--codon-params")
parser.add_argument("--output")
parser.add_argument("--adabmdca-device")
parser.add_argument("--mode", default="success")
arguments = parser.parse_args()
side = "protein" if arguments.protein_params else "codon"
params = Path(arguments.protein_params or arguments.codon_params)
if params.is_dir():
    params = params / f"{arguments.gene}.{side}_adabm_params"
reused = params.exists()
if not reused:
    params.write_text("partial checkpoint")
Path("observed.json").write_text(json.dumps({
    "cwd": os.getcwd(), "reused": reused, "pid": os.getpid(),
    "environment": {name: os.environ.get(name) for name in (
        "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "BFF_TASK_THREADS", "CUDA_VISIBLE_DEVICES",
    )},
    "device": arguments.adabmdca_device,
}))
if "BFF_TEST_RESOURCES" in os.environ:
    Path(os.environ["BFF_TASK_RESOURCE_REPORT"]).write_text(os.environ["BFF_TEST_RESOURCES"])
if arguments.mode == "wait":
    while not Path(os.environ["BFF_TEST_GATE"]).exists():
        time.sleep(0.02)
if arguments.codon_msa:
    Path(arguments.codon_msa + ".encoded.fasta").write_text("encoded")
if arguments.mode == "cuda_oom":
    print("torch.OutOfMemoryError: CUDA out of memory. Tried to allocate 9.04 GiB", file=sys.stderr)
    sys.exit(1)
if arguments.mode == "cpu_oom":
    print("RuntimeError: DefaultCPUAllocator: can't allocate memory", file=sys.stderr)
    sys.exit(1)
if arguments.mode == "ordinary_error":
    print("ValueError: invalid alignment", file=sys.stderr)
    sys.exit(7)
if arguments.mode == "oom_hint":
    print("Hint: CUDA out of memory can happen on small GPUs", file=sys.stderr)
    sys.exit(2)
if arguments.mode == "missing":
    sys.exit(0)
if not reused:
    params.write_text("complete parameters")
if arguments.mode == "empty_params":
    params.write_text("")
output = Path(arguments.output) / arguments.gene / "adabmDCA"
output.mkdir(parents=True, exist_ok=True)
tsv = output / f"{arguments.gene}.{side}.tsv"
tsv.write_text("pkey\\tscore\\nGENE:A1C\\t0.5\\n" if arguments.mode != "empty_tsv" else "")
'''


@pytest.fixture
def task(tmp_path, monkeypatch):
    backend = tmp_path / "backend.py"
    backend.write_text(FAKE_BACKEND)
    fasta = tmp_path / "gene.fasta"
    fasta.write_text(">GENE\nAAA\n")
    msa = tmp_path / "gene.msa.fasta"
    msa.write_text(">GENE\nAAA\n")
    executable_directory = tmp_path / "fake-bin"
    executable_directory.mkdir()
    nvidia_smi = executable_directory / "nvidia-smi"
    nvidia_smi.write_text(f"#!{sys.executable}\nprint('GPU-first, 81920\\nGPU-second, 81920')\n")
    nvidia_smi.chmod(0o755)
    monkeypatch.setenv("PATH", str(executable_directory) + os.pathsep + os.environ["PATH"])
    plan = {
        "gene": "GENE", "side": "protein", "fingerprint": "input-and-model-hash",
        "cpu_memory_gib": 128, "gpu_memory_gib": 60, "params": None,
        "eligible_gpu_uuids": ["GPU-first"], "lease_dir": str(tmp_path / "leases"),
        "gpu_wait_timeout": 0,
    }
    return tmp_path, backend, plan, ["--fasta", str(fasta), "--msa", str(msa)]


def _attempts(directory):
    return sorted((directory / ".adabmdca-attempts").iterdir())


def _observed(directory):
    return [json.loads((attempt / "observed.json").read_text()) for attempt in _attempts(directory)]


def test_success_publishes_only_complete_hashed_artifacts(task):
    directory, backend, plan, arguments = task
    assert run_task(plan, "cpu", backend, 3, arguments, directory) == 0
    manifest = directory / "GENE.protein.complete.json"
    assert verify_completion(manifest, plan["fingerprint"])
    assert not verify_completion(manifest, "different-inputs")
    assert (directory / "GENE.protein_adabm_params").read_text() == "complete parameters"
    (directory / "GENE.protein.tsv").write_text("changed")
    assert not verify_completion(manifest, plan["fingerprint"])


def test_completion_verifies_separately_published_params_without_copying(task):
    directory, backend, plan, arguments = task
    assert run_task(plan, "cpu", backend, 3, arguments, directory) == 0
    params = directory / "GENE.protein_adabm_params"
    published = directory / "published_params"
    published.mkdir()
    relocated = params.rename(published / params.name)
    manifest = directory / "GENE.protein.complete.json"
    paths = {params.name: relocated, "GENE.protein.tsv": directory / "GENE.protein.tsv"}
    assert verify_completion(manifest, plan["fingerprint"], paths)
    assert not verify_completion(manifest, plan["fingerprint"])
    wrong_name = published / "OTHER.protein_adabm_params"
    relocated.rename(wrong_name)
    assert not verify_completion(manifest, plan["fingerprint"], {**paths, params.name: wrong_name})
    wrong_name.rename(relocated)
    relocated.write_text("modified parameters")
    assert not verify_completion(manifest, plan["fingerprint"], paths)


@pytest.mark.parametrize("mutation", ["filename", "success", "params_hash", "extra_file", "missing_file", "empty_params"])
def test_completion_rejects_inconsistent_artifacts(task, mutation):
    directory, backend, plan, arguments = task
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    manifest = directory / "GENE.protein.complete.json"
    contents = json.loads(manifest.read_text())
    if mutation == "filename":
        manifest = directory / "OTHER.protein.complete.json"
    elif mutation == "success":
        contents["success"] = False
    elif mutation == "params_hash":
        contents["files"]["GENE.protein_adabm_params"]["sha256"] = "wrong"
    elif mutation == "extra_file":
        contents["files"]["../unrelated"] = {}
    elif mutation == "missing_file":
        (directory / "GENE.protein.tsv").unlink()
    elif mutation == "empty_params":
        (directory / "GENE.protein_adabm_params").write_text("")
    manifest.write_text(json.dumps(contents))
    assert not verify_completion(manifest, plan["fingerprint"])


def test_gpu_oom_leaves_partial_checkpoint_private_and_cpu_retry_is_fresh(task):
    directory, backend, plan, arguments = task
    assert run_task(plan, "cuda", backend, 2, [*arguments, "--mode", "cuda_oom"], directory) == 0
    oom = json.loads((directory / "GENE.protein.gpu_oom.json").read_text())
    assert oom["reason"] == "cuda_oom"
    assert oom["fingerprint"] == plan["fingerprint"]
    assert not (directory / "GENE.protein.tsv").exists()
    assert not (directory / "GENE.protein_adabm_params").exists()
    first_attempt = _attempts(directory)[0]
    assert (first_attempt / "params" / "GENE.protein_adabm_params").read_text() == "partial checkpoint"
    assert run_task(plan, "cpu", backend, 4, arguments, directory) == 0
    assert len(_attempts(directory)) == 2
    assert all(not observed["reused"] for observed in _observed(directory))
    assert not (directory / "GENE.protein.gpu_oom.json").exists()
    assert verify_completion(directory / "GENE.protein.complete.json", plan["fingerprint"])


@pytest.mark.parametrize("device,mode,expected", [
    ("cuda", "ordinary_error", 7), ("cuda", "cpu_oom", 1),
    ("cpu", "cpu_oom", 1), ("cpu", "cuda_oom", 1), ("cuda", "oom_hint", 2),
    ("cuda", "missing", 1), ("cpu", "empty_tsv", 1), ("cpu", "empty_params", 1),
])
def test_failures_do_not_publish_completion_or_trigger_wrong_fallback(task, device, mode, expected):
    directory, backend, plan, arguments = task
    assert run_task(plan, device, backend, 1, [*arguments, "--mode", mode], directory) == expected
    assert not (directory / "GENE.protein.complete.json").exists()
    assert not (directory / "GENE.protein.gpu_oom.json").exists()
    assert not (directory / "GENE.protein.tsv").exists()
    assert not (directory / "GENE.protein_adabm_params").exists()


@pytest.mark.parametrize("device,expected", [
    ("cuda", "GPU-first"), ("cpu", ""),
])
def test_thread_limits_and_leased_gpu_binding(task, monkeypatch, device, expected):
    directory, backend, plan, arguments = task
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "GPU-existing")
    for variable in THREAD_VARIABLES:
        monkeypatch.setenv(variable, "99")
    assert run_task(plan, device, backend, 6, arguments, directory) == 0
    observed = _observed(directory)[0]
    assert observed["device"] == device
    assert observed["environment"]["CUDA_VISIBLE_DEVICES"] == expected
    assert all(observed["environment"][variable] == "6" for variable in THREAD_VARIABLES)
    assert os.environ["CUDA_VISIBLE_DEVICES"] == "GPU-existing"


def test_gpu_lease_requires_explicit_identity_and_current_free_memory(task):
    directory, backend, plan, arguments = task
    plan["eligible_gpu_uuids"] = []
    with pytest.raises(ValueError, match="explicit eligible"):
        run_task(plan, "cuda", backend, 1, arguments, directory)
    plan["eligible_gpu_uuids"] = ["GPU-first"]
    with pytest.raises(ValueError, match="timed out"):
        with lease_gpu(plan, probe=lambda: {"GPU-first": 59}):
            pytest.fail("Insufficient free memory must not launch the backend")
    assert not (directory / "GENE.protein.gpu_oom.json").exists()
    with lease_gpu(plan, probe=lambda: {"GPU-first": 80}) as (uuid, descriptor):
        assert uuid == "GPU-first"
        assert descriptor >= 0


def test_gpu_probe_rejects_malformed_memory(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(args, 0, "GPU-first, N/A"))
    with pytest.raises(ValueError, match="Invalid GPU resource probe response"):
        _gpu_free_memory()


def test_gpu_lease_rechecks_external_memory_pressure(task):
    directory, backend, plan, arguments = task
    observations = iter([{"GPU-first": 20}, {"GPU-first": 80}])
    with lease_gpu({**plan, "gpu_wait_timeout": 1}, probe=lambda: next(observations), poll_interval=0.001) as (uuid, descriptor):
        assert uuid == "GPU-first"
    assert list(observations) == []


@pytest.mark.parametrize("report", ['{"cuda_peak_allocated_gib": 3.25}', "unfinished {", "[]"])
@pytest.mark.parametrize("mode", ["success", "cuda_oom"])
def test_resource_reports_do_not_mask_backend_result(task, monkeypatch, report, mode):
    directory, backend, plan, arguments = task
    monkeypatch.setenv("BFF_TEST_RESOURCES", report)
    assert run_task(plan, "cuda", backend, 1, [*arguments, "--mode", mode], directory) == 0
    suffix = "complete" if mode == "success" else "gpu_oom"
    manifest = json.loads((directory / f"GENE.protein.{suffix}.json").read_text())
    assert manifest["subprocess_peak_rss_gib"] >= 0
    assert manifest["backend_returncode"] == (0 if mode == "success" else 1)
    if report.startswith('{"'):
        assert manifest["backend_resources"] == {"cuda_peak_allocated_gib": 3.25}
    else:
        assert "backend_resources" not in manifest
    attempt = json.loads((_attempts(directory)[0] / "attempt.json").read_text())
    assert attempt["subprocess_peak_rss_gib"] == manifest["subprocess_peak_rss_gib"]


def _wrapper_command(directory, backend, plan, arguments, device, mode="success"):
    directory.mkdir(exist_ok=True)
    plan_path = directory / "plan.json"
    plan_path.write_text(json.dumps(plan))
    wrapper = Path(__file__).resolve().parents[2] / "biofeaturefactory" / "mutation_effects" / "bin" / "adabmdca_task.py"
    return [
        sys.executable, str(wrapper), "--plan", str(plan_path), "--device", device,
        "--backend-script", str(backend), "--threads", "1", "--", *arguments,
        "--mode", mode,
    ]


def _wait_for_observation(directory, timeout=10):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        candidates = list((directory / ".adabmdca-attempts").glob("*/observed.json"))
        if candidates:
            try:
                return json.loads(candidates[0].read_text())
            except json.JSONDecodeError:
                pass
        time.sleep(0.01)
    pytest.fail(f"Backend did not start in {directory}")


def test_two_gpu_processes_run_exclusively_while_cpu_bypasses_gpu_queue(task, monkeypatch):
    directory, backend, plan, arguments = task
    plan.update(eligible_gpu_uuids=["GPU-first", "GPU-second"], gpu_wait_timeout=10)
    gate = directory / "release-backends"
    monkeypatch.setenv("BFF_TEST_GATE", str(gate))
    processes = []
    try:
        for name in ("gpu-one", "gpu-two"):
            work = directory / name
            processes.append(subprocess.Popen(_wrapper_command(work, backend, plan, arguments, "cuda", "wait"), cwd=work))
        first = _wait_for_observation(directory / "gpu-one")
        second = _wait_for_observation(directory / "gpu-two")
        assert {first["environment"]["CUDA_VISIBLE_DEVICES"], second["environment"]["CUDA_VISIBLE_DEVICES"]} == {"GPU-first", "GPU-second"}
        queued_work = directory / "gpu-queued"
        processes.append(subprocess.Popen(_wrapper_command(queued_work, backend, plan, arguments, "cuda"), cwd=queued_work))
        cpu_work = directory / "cpu"
        cpu = subprocess.Popen(_wrapper_command(cpu_work, backend, plan, arguments, "cpu"), cwd=cpu_work)
        processes.append(cpu)
        assert cpu.wait(timeout=10) == 0
        assert _wait_for_observation(cpu_work)["environment"]["CUDA_VISIBLE_DEVICES"] == ""
        assert not list((queued_work / ".adabmdca-attempts").glob("*/observed.json"))
        gate.touch()
        assert all(process.wait(timeout=10) == 0 for process in processes)
        assert _wait_for_observation(queued_work)["environment"]["CUDA_VISIBLE_DEVICES"] in {"GPU-first", "GPU-second"}
    finally:
        gate.touch()
        for process in processes:
            if process.poll() is None:
                process.kill()
            process.wait(timeout=10)


def test_os_releases_gpu_lock_after_holder_crashes(task):
    directory, backend, plan, arguments = task
    code = (
        "import json, os, sys; "
        "from biofeaturefactory.mutation_effects.bin.adabmdca_task import lease_gpu; "
        "allocation = lease_gpu(json.loads(sys.argv[1])); "
        "allocation.__enter__(); os._exit(29)"
    )
    result = subprocess.run([sys.executable, "-c", code, json.dumps(plan)], timeout=10)
    assert result.returncode == 29
    with lease_gpu(plan) as (uuid, descriptor):
        assert uuid == "GPU-first"


def test_crashed_wrapper_keeps_gpu_locked_until_inherited_backend_exits(task, monkeypatch):
    directory, backend, plan, arguments = task
    gate = directory / "release-backend"
    monkeypatch.setenv("BFF_TEST_GATE", str(gate))
    work = directory / "crashed-wrapper"
    process = subprocess.Popen(_wrapper_command(work, backend, plan, arguments, "cuda", "wait"), cwd=work)
    try:
        _wait_for_observation(work)
        process.kill()
        process.wait(timeout=10)
        with pytest.raises(ValueError, match="timed out"):
            with lease_gpu(plan):
                pytest.fail("The orphaned backend must retain its inherited GPU lock")
        gate.touch()
        with lease_gpu({**plan, "gpu_wait_timeout": 10}) as (uuid, descriptor):
            assert uuid == "GPU-first"
    finally:
        gate.touch()
        if process.poll() is None:
            process.kill()
        process.wait(timeout=10)


def test_explicit_params_are_snapshotted_and_published(task):
    directory, backend, plan, arguments = task
    source = directory / "explicit.params"
    source.write_text("user supplied parameters")
    plan["params"] = str(source)
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    assert _observed(directory)[0]["reused"]
    assert (directory / "GENE.protein_adabm_params").read_text() == source.read_text()
    assert verify_completion(directory / "GENE.protein.complete.json", plan["fingerprint"])


def test_explicit_cli_params_directory_resolves_matching_gene(task):
    directory, backend, plan, arguments = task
    del plan["params"]
    supplied = directory / "supplied"
    supplied.mkdir()
    (supplied / "GENE.protein_adabm_params").write_text("explicit matching model")
    (supplied / "OTHER.protein_adabm_params").write_text("unrelated model")
    arguments += ["-pp", str(supplied)]
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    assert _observed(directory)[0]["reused"]
    assert (directory / "GENE.protein_adabm_params").read_text() == "explicit matching model"


def test_plan_null_params_prevents_implicit_cached_params_reuse(task):
    directory, backend, plan, arguments = task
    stale = directory / "GENE.protein_adabm_params"
    stale.write_text("partial checkpoint from an old run")
    arguments += ["--protein-params", str(stale)]
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    assert not _observed(directory)[0]["reused"]
    assert stale.read_text() == "complete parameters"


def test_failed_rerun_does_not_expose_stale_success(task):
    directory, backend, plan, arguments = task
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    plan["fingerprint"] = "changed-inputs"
    assert run_task(plan, "cuda", backend, 1, [*arguments, "--mode", "cuda_oom"], directory) == 0
    assert not (directory / "GENE.protein.complete.json").exists()
    assert not (directory / "GENE.protein.tsv").exists()
    previous = list((directory / ".adabmdca-attempts").glob("*/previous/GENE.protein.complete.json"))
    assert len(previous) == 1


def test_codon_encoding_stays_inside_attempt(task):
    directory, backend, plan, arguments = task
    plan["side"] = "codon"
    arguments[2] = "--codon-msa"
    original_msa = Path(arguments[3])
    assert run_task(plan, "cpu", backend, 1, arguments, directory) == 0
    assert not Path(str(original_msa) + ".encoded.fasta").exists()
    assert list(_attempts(directory)[0].glob("inputs/*.encoded.fasta"))
    assert verify_completion(directory / "GENE.codon.complete.json", plan["fingerprint"])


def test_cli_forces_attempt_output_and_planned_device(task):
    directory, backend, plan, arguments = task
    plan_path = directory / "plan.json"
    plan_path.write_text(json.dumps(plan))
    wrapper = Path(__file__).resolve().parents[2] / "biofeaturefactory" / "mutation_effects" / "bin" / "adabmdca_task.py"
    result = subprocess.run([
        sys.executable, str(wrapper), "--plan", str(plan_path), "--device", "cpu",
        "--backend-script", str(backend), "--threads", "2", "--", *arguments,
        "--output", str(directory / "wrong-output"), "--adabmdca-device", "cuda",
    ], cwd=directory, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert not (directory / "wrong-output").exists()
    assert verify_completion(directory / "GENE.protein.complete.json", plan["fingerprint"])
