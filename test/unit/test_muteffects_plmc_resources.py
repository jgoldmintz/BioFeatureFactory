"""Controller integration for model-sized EVmutation RAM admission."""

import json
import struct
import sys

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller
from biofeaturefactory.mutation_effects.bin import evmutation_cache, plmc_resources
from biofeaturefactory.mutation_effects.bin.codon_encoding import CODON_ALPHABET


def _write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


@pytest.fixture
def make_panel(tmp_path):
    def create_panel(label, sites=100, gene_count=1, memory_gib=64, cpus=16,
                     msa_sides=("protein", "codon"), gpus=()):
        directory = tmp_path / label
        source = directory / "source"
        genes = [f"GENE{index + 1}" for index in range(gene_count)]
        files = {}
        for gene in genes:
            orf = "ATG" + "GCT" * (sites - 1)
            protein = "M" + "A" * (sites - 1)
            files[gene] = {
                "fasta": _write(source / gene / "fastas" / f"{gene}.fasta", f">ORF\n{orf}\n"),
                "mutations": _write(source / gene / "mappings" / "mutations" / f"{gene}_mutations.csv", "mutant\nG4A\n"),
            }
            if "protein" in msa_sides:
                files[gene]["protein"] = _write(
                    source / gene / "MSA" / f"{gene}.msa.a2m",
                    f">{gene}\n{protein}\n>homolog\n{protein}\n",
                )
            if "codon" in msa_sides:
                files[gene]["codon"] = _write(
                    source / gene / "CodonMSA" / f"{gene}.codon.msa.fasta",
                    f">ORF\n{orf}\n>homolog\n{orf}\n",
                )
        hardware = _write(directory / "hardware.json", json.dumps({
            "cpus": cpus, "memory_gib": memory_gib, "gpus": list(gpus),
        }))
        return {
            "source": source, "output": directory / "output", "genes": genes,
            "hardware": hardware, "files": files,
        }

    return create_panel


def _parse(monkeypatch, panel, *options, both_backends=False):
    arguments = [
        "mutEffects_controller", "--fasta", str(panel["source"]),
        "--output", str(panel["output"]),
        "--resource-hardware", str(panel["hardware"]),
    ]
    if not both_backends:
        arguments.append("--evmutation-only")
    monkeypatch.setattr(sys, "argv", [*arguments, *map(str, options)])
    return controller.parse_args()


def _prepare(args, panel):
    manifest = controller.build_manifest(panel["genes"], args)
    config, adabm_plans = controller.prepare_resource_plans(args, manifest, panel["genes"])
    return config, adabm_plans, manifest


def _native_params(panel, gene, side, sites=3):
    alphabet = CODON_ALPHABET[1:] if side == "codon" else "ACDEFGHIKLMNPQRSTVWY"
    states = len(alphabet)
    suffix = "codon_model_params" if side == "codon" else "model_params"
    path = panel["output"] / suffix / f"{gene}.{suffix}"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = struct.pack("=5i5f", sites, states, 2, 0, 1, 0.2, 0.01, 16.2, 0.0, 2.0)
    payload += alphabet.encode() + struct.pack("=2f", 1.0, 1.0)
    payload += (alphabet[1] * sites).encode()
    payload += struct.pack(f"={sites}i", *range(1, sites + 1))
    payload += struct.pack(f"={sites * states}f", *([1 / states] * (sites * states)))
    payload += bytes(4 * sites * states)
    payload += bytes(8 * (sites * (sites - 1) // 2) * states * states)
    path.write_bytes(payload)
    return path


def _deny_external_process(*arguments, **keywords):
    pytest.fail("Resource planning unexpectedly launched an external process")


def test_model_length_changes_ram_concurrency_not_mutation_count(make_panel, monkeypatch):
    small_panel = make_panel("small", sites=100, gene_count=4)
    large_panel = make_panel("large", sites=1200, gene_count=4)
    small_args = _parse(monkeypatch, small_panel)
    large_args = _parse(monkeypatch, large_panel)
    small, small_adabm, small_manifest = _prepare(small_args, small_panel)
    large, large_adabm, large_manifest = _prepare(large_args, large_panel)

    assert small_adabm == large_adabm == []
    assert small["hardware"] == large["hardware"]
    assert small["evmutation"]["memory_gib"] is None
    assert small["cpu_allocation"] == {"threads": 4, "concurrent_jobs": 4}
    assert large["cpu_allocation"] == {"threads": 16, "concurrent_jobs": 1}
    assert len(small["evmutation_plans"]) == len(large["evmutation_plans"]) == 4
    assert max(plan["memory_gib"] for plan in small["evmutation_plans"]) < 8
    assert min(plan["memory_gib"] for plan in large["evmutation_plans"]) > 32
    for gene in small_panel["genes"]:
        assert small_panel["files"][gene]["mutations"].read_bytes() == large_panel["files"][gene]["mutations"].read_bytes()
    for config, args, panel, manifest in (
        (small, small_args, small_panel, small_manifest),
        (large, large_args, large_panel, large_manifest),
    ):
        tasks = controller.pending_resource_tasks(args, manifest, panel["genes"], config, [])
        assert [task["memory_gib"] for task in tasks] == [plan["memory_gib"] for plan in config["evmutation_plans"]]
        assert all(plan["threads"] == config["threads"] for plan in config["evmutation_plans"])


def test_model_side_changes_estimated_state_count_and_ram(make_panel, monkeypatch):
    panel = make_panel("sides", sites=300)
    protein_args = _parse(monkeypatch, panel, "--msa", panel["source"])
    codon_args = _parse(monkeypatch, panel, "--codon-msa", panel["source"])
    protein_config, _, _ = _prepare(protein_args, panel)
    codon_config, _, _ = _prepare(codon_args, panel)
    protein, = protein_config["evmutation_plans"]
    codon, = codon_config["evmutation_plans"]

    assert protein["side"] == "protein"
    assert codon["side"] == "codon"
    assert protein["sites"] == codon["sites"] == 300
    assert protein["states"] == 20
    assert codon["states"] == 64
    assert codon["memory_gib"] > 5 * protein["memory_gib"]


@pytest.mark.parametrize("side", ["protein", "codon"])
def test_native_prebuilt_models_use_header_scoring_shape_without_plmc(
    make_panel, monkeypatch, side,
):
    panel = make_panel(f"prebuilt_{side}", sites=300)
    option = "--msa" if side == "protein" else "--codon-msa"
    args = _parse(monkeypatch, panel, option, panel["source"])
    training_config, _, _ = _prepare(args, panel)
    params = _native_params(panel, "GENE1", side)
    config, adabm_plans, manifest = _prepare(args, panel)
    plan, = config["evmutation_plans"]

    assert not adabm_plans
    assert args.plmc_binary is None
    assert plan["prebuilt"]
    assert plan["params"] == str(params.resolve())
    assert plan["sites"] == 3
    assert plan["states"] == (64 if side == "codon" else 20)
    assert plan["training_memory_gib"] == 0
    assert plan["memory_gib"] == plan["scoring_memory_gib"]
    assert plan["memory_gib"] < training_config["evmutation_plans"][0]["memory_gib"]
    monkeypatch.setattr(controller.importlib.util, "find_spec", _deny_external_process)
    controller.validate_backend_tools(args, manifest, panel["genes"])
    controller.validate_db_coverage(panel["genes"], manifest, args)


def test_verified_ev_result_is_excluded_before_resource_estimation(make_panel, monkeypatch):
    panel = make_panel("completed")
    args = _parse(monkeypatch, panel)
    params = _native_params(panel, "GENE1", "protein")
    before = controller.build_manifest(panel["genes"], args)
    directory = panel["output"] / "GENE1" / "EVmutation"
    tsv = _write(directory / "GENE1.protein.tsv", "pkey\tscore\nGENE1:G4A\t0.25\n")
    evmutation_cache.write_completion(
        directory / "GENE1.protein.routing.json",
        before["ev_fingerprints"]["GENE1"]["protein"], tsv, params,
        msa=panel["files"]["GENE1"]["protein"],
    )
    monkeypatch.setattr(plmc_resources, "plan_evmutation_task", _deny_external_process)

    config, adabm_plans, manifest = _prepare(args, panel)

    assert manifest["EVmutation"] == ["GENE1"]
    assert config["evmutation_plans"] == adabm_plans == []
    assert controller.pending_resource_tasks(args, manifest, panel["genes"], config, []) == []
    assert config["cpu_allocation"]["concurrent_jobs"] == 0


@pytest.mark.parametrize("memory_floor", [8, 60])
def test_explicit_memory_is_a_floor_not_an_unsafe_estimate_override(
    make_panel, monkeypatch, memory_floor,
):
    panel = make_panel(f"floor_{memory_floor}", sites=1200)
    automatic, _, _ = _prepare(_parse(monkeypatch, panel), panel)
    explicit, _, _ = _prepare(
        _parse(monkeypatch, panel, "--evmutation-memory-gib", memory_floor), panel,
    )
    estimated = automatic["evmutation_plans"][0]["memory_gib"]
    plan, = explicit["evmutation_plans"]

    assert estimated > 8
    assert plan["memory_gib"] == max(memory_floor, estimated)
    assert explicit["evmutation"]["memory_gib"] == memory_floor
    assert plan["memory_floor_gib"] == memory_floor


def test_missing_msa_defers_ev_estimate_and_counts_one_shared_generation_task(
    make_panel, monkeypatch,
):
    panel = make_panel("missing_msa", msa_sides=())
    args = _parse(monkeypatch, panel, "--msa-memory-gib", 12, both_backends=True)
    config, adabm_plans, manifest = _prepare(args, panel)
    tasks = controller.pending_resource_tasks(args, manifest, panel["genes"], config, adabm_plans)

    assert config["evmutation_plans"] == adabm_plans == []
    assert tasks == [{
        "id": "GENE1.protein.msa", "memory_gib": 12,
        "eligible_gpu_uuids": [], "threads": None,
    }]
    assert config["cpu_allocation"] == {"threads": 16, "concurrent_jobs": 1}


def test_ev_and_gpu_host_memory_share_the_same_admission_budget(make_panel, monkeypatch):
    panel = make_panel("shared_ram", sites=1200, memory_gib=256, gpus=[{
        "index": "0", "uuid": "GPU-allocated", "memory_gib": 80,
    }])
    args = _parse(monkeypatch, panel, both_backends=True)
    initial, initial_adabm, _ = _prepare(args, panel)
    ev_memory = initial["evmutation_plans"][0]["memory_gib"]
    gpu_host_memory = initial_adabm[0]["gpu_host_memory_gib"]
    hardware = json.loads(panel["hardware"].read_text())
    hardware["memory_gib"] = max(ev_memory, gpu_host_memory) + min(ev_memory, gpu_host_memory) / 2
    panel["hardware"].write_text(json.dumps(hardware))

    config, adabm_plans, manifest = _prepare(args, panel)
    tasks = controller.pending_resource_tasks(args, manifest, panel["genes"], config, adabm_plans)

    assert adabm_plans[0]["device"] == "cuda"
    assert len(tasks) == 2
    assert {task["id"] for task in tasks} == {"GENE1.protein.evmutation", "GENE1.protein.adabmdca"}
    assert sorted(task["memory_gib"] for task in tasks) == sorted([ev_memory, gpu_host_memory])
    assert sum(task["memory_gib"] for task in tasks) > hardware["memory_gib"]
    assert config["cpu_allocation"] == {"threads": 16, "concurrent_jobs": 1}


@pytest.mark.parametrize("memory_floor", [None, 12])
def test_nextflow_command_forwards_ev_memory_only_when_explicit(
    make_panel, monkeypatch, memory_floor,
):
    panel = make_panel(f"command_{memory_floor}")
    options = [] if memory_floor is None else ["--evmutation-memory-gib", memory_floor]
    args = _parse(monkeypatch, panel, *options)
    config, _, _ = _prepare(args, panel)
    args.resource_config = config
    args.resource_config_path = panel["output"] / "resources.json"

    command = controller.build_nextflow_cmd(args, "manifest.json")

    if memory_floor is None:
        assert "--evmutation_memory" not in command
    else:
        assert command[command.index("--evmutation_memory") + 1] == "12.0 GB"
    assert command[command.index("--evmutation_cpus") + 1] == str(config["threads"])


def test_plan_only_reports_ev_estimates_without_outputs_or_launch(make_panel, monkeypatch, capsys):
    panel = make_panel("plan_only", sites=1200)
    args = _parse(monkeypatch, panel, "--resource-plan-only")
    capsys.readouterr()
    monkeypatch.setattr(controller.subprocess, "Popen", _deny_external_process)
    monkeypatch.setattr(controller, "write_resource_snapshot", _deny_external_process)

    controller.run_controller(args)

    report = json.loads(capsys.readouterr().out)
    task, = report["tasks"]
    assert task["backend"] == "evmutation"
    assert task["sites"] == 1200
    assert task["memory_gib"] > 32
    assert report["cpu_allocation"]["concurrent_jobs"] == 1
    assert not panel["output"].exists()


@pytest.mark.parametrize("value", ["0", "-1", "nan", "inf"])
def test_invalid_explicit_ev_memory_is_rejected(make_panel, monkeypatch, value):
    panel = make_panel(f"invalid_{value}")

    with pytest.raises(SystemExit) as exit_info:
        _parse(monkeypatch, panel, "--evmutation-memory-gib", value)

    assert exit_info.value.code == 2


@pytest.mark.parametrize("options", [(), ("--evmutation-memory-gib", 128)])
def test_unschedulable_ev_request_fails_before_any_launch(make_panel, monkeypatch, capsys, options):
    panel = make_panel(f"unschedulable_{len(options)}", sites=1200, memory_gib=8)
    args = _parse(monkeypatch, panel, *options)
    monkeypatch.setattr(controller.subprocess, "Popen", _deny_external_process)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    assert "EVmutation unschedulable" in capsys.readouterr().err
    assert not panel["output"].exists()
