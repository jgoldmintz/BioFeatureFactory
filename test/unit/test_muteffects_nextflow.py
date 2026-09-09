"""Real controller/Nextflow input routing; scoring is replaced by a local stub.

Run with BFF_TEST_NEXTFLOW=1 to enable the Nextflow integration check.
"""

import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import textwrap

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller


ADABMDCA_TEST_SETTINGS = {
    "model": "bmDCA", "dtype": "float32", "nchains": 64, "nepochs": 1,
    "tol": 0.001, "patience": 3, "check_every": 10, "target": 0.95,
    "lr": 0.01, "nsweeps": 10, "seed": 0, "skip_codon": False,
    "validation_log": None,
}


def ev_params_bytes(states):
    alphabet = b"-ACDEFGHIKLMNPQRSTVWY" if states == 21 else bytes(range(33, 98))
    return (
        struct.pack("=5i5f", 2, states, 1, 0, 500, 0.2, 0.01, 0.01, 0.0, 1.0)
        + alphabet + struct.pack("=f", 1.0) + b"AA" + struct.pack("=2i", 1, 2)
        + bytes((4 * states + 2 * states * states) * 4)
    )


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test routing",
)
@pytest.mark.parametrize("backend,prebuilt_ev,automatic_threads,generate_protein,blocked_codon", [
    ("adabmDCA", False, False, False, False), ("EVmutation", False, False, False, False),
    ("EVmutation", True, False, False, False), ("adabmDCA", False, True, False, False),
    ("EVmutation", False, True, False, False), ("EVmutation", False, False, True, False),
    ("EVmutation", False, False, False, True),
])
def test_controller_routes_per_gene_inputs_through_nextflow(tmp_path, backend, prebuilt_ev, automatic_threads, generate_protein, blocked_codon):
    source_root = tmp_path / "source"
    output_root = tmp_path / "results"
    runtime = tmp_path / "runtime"
    (runtime / "bin").mkdir(parents=True)
    shutil.copy2(controller.__file__, runtime / "mutEffects_controller.py")
    shutil.copy2(controller.NEXTFLOW_SCRIPT, runtime / "bin" / "main.nf")
    for name in ("resource_planner.py", "plmc_resources.py", "adabmdca_task.py", "evmutation_cache.py", "codon_encoding.py", "__init__.py"):
        source = Path(controller.NEXTFLOW_SCRIPT).parent / name
        shutil.copy2(source, runtime / "bin" / name)
    (runtime / "nextflow.config").write_text(
        "params.msa_cpus = 1\n"
        "params.msa_memory = '512 MB'\n"
        "params.evmutation_cpus = 1\n"
        "params.evmutation_memory = '512 MB'\n"
        "trace { enabled = true; file = 'trace.tsv'; overwrite = true; fields = 'name,cpus' }\n"
    )
    stub = textwrap.dedent('''\
        import argparse
        import json
        import os
        from pathlib import Path
        import struct

        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', type=Path)
        parser.add_argument('--mutations', type=Path)
        parser.add_argument('--msa', type=Path)
        parser.add_argument('--codon-msa', type=Path)
        parser.add_argument('--output', type=Path)
        parser.add_argument('--protein-params', type=Path)
        parser.add_argument('--codon-params', type=Path)
        parser.add_argument('--model-params', type=Path)
        parser.add_argument('--codon-model-params', type=Path)
        parser.add_argument('--plmc-binary')
        parser.add_argument('--score-missense-codon', action='store_true')
        parser.add_argument('--adabmdca-model')
        args, unused = parser.parse_known_args()
        gene = args.fasta.read_text().splitlines()[0][1:]
        side = 'protein' if args.msa else 'codon'
        msa = args.msa or args.codon_msa
        if side == 'codon' and os.environ.get('BFF_READY_SIDE_MARKER'):
            Path(os.environ['BFF_READY_SIDE_MARKER']).touch()
        if os.environ.get('BFF_EXPECT_FORCED_CODON'):
            assert side == 'codon' and args.score_missense_codon
        assert msa.read_text().splitlines()[0] == '>' + gene
        assert args.mutations.read_text().splitlines()[0] == 'mutant'
        for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            assert os.environ[variable] == os.environ['BFF_EXPECT_THREADS']
        observed = {
            'fasta': str(args.fasta.resolve()),
            'mutations': str(args.mutations.resolve()),
            'msa' if args.msa else 'codon_msa': msa.read_text(),
        }
        if Path(__file__).stem == 'adabmdca_pipeline':
            assert args.adabmdca_model == 'pseudoDCA'
            assert os.environ['BFF_TASK_THREADS'] == os.environ['BFF_EXPECT_THREADS']
            suffix = f'{side}_adabm_params'
            output = args.output / gene / 'adabmDCA'
            output.mkdir(parents=True)
            params_path = args.protein_params or args.codon_params
            if params_path.is_dir():
                params_path = params_path / f'{gene}.{suffix}'
        else:
            suffix = 'model_params' if side == 'protein' else 'codon_model_params'
            output = Path('.')
            params_path = Path(f'{gene}.{suffix}')
            if os.environ.get('BFF_EXPECT_EV_PARAMS'):
                supplied = args.model_params or args.codon_model_params
                assert supplied.is_file()
                assert not args.plmc_binary
        (output / f'{gene}.{side}.tsv').write_text(json.dumps(observed))
        if not os.environ.get('BFF_EXPECT_EV_PARAMS'):
            if Path(__file__).stem == 'adabmdca_pipeline':
                params_path.write_text('routing-test stub; not model parameters')
            else:
                states = 21 if side == 'protein' else 65
                alphabet = b'-ACDEFGHIKLMNPQRSTVWY' if states == 21 else bytes(range(33, 98))
                params_path.write_bytes(
                    struct.pack('=5i5f', 2, states, 1, 0, 500, 0.2, 0.01, 0.01, 0.0, 1.0)
                    + alphabet + struct.pack('=f', 1.0) + b'AA' + struct.pack('=2i', 1, 2)
                    + bytes((4 * states + 2 * states * states) * 4)
                )
    ''')
    for name in ("adabmdca_pipeline.py", "evmutation_pipeline.py"):
        (runtime / name).write_text(stub)

    if generate_protein:
        package = runtime / "biofeaturefactory"
        core = package / "core"
        core.mkdir(parents=True)
        (package / "__init__.py").write_text("from pkgutil import extend_path\n__path__ = extend_path(__path__, __name__)\n")
        (core / "__init__.py").write_text("")
        (core / "msa_generation_pipeline.py").write_text(textwrap.dedent('''\
            import argparse
            import os
            from pathlib import Path
            import time

            parser = argparse.ArgumentParser()
            parser.add_argument('--fasta', type=Path)
            parser.add_argument('--threads', type=int)
            args, unused = parser.parse_known_args()
            assert args.threads == int(os.environ['BFF_EXPECT_THREADS'])
            marker = Path(os.environ['BFF_READY_SIDE_MARKER'])
            deadline = time.monotonic() + 20
            while not marker.exists() and time.monotonic() < deadline:
                time.sleep(0.05)
            assert marker.exists(), 'Ready codon EV scoring waited for protein MSA generation'
            gene = args.fasta.read_text().splitlines()[0][1:]
            Path(f'{gene}.msa.a2m').write_text(f'>{gene}\\nMA\\n')
            Path(f'{gene}.msa.stats.json').write_text('{}')
        '''))

    expected = {}
    for gene in (("NPM1",) if automatic_threads or generate_protein or blocked_codon else ("NPM1", "PAM")):
        files = {
            "fasta": (f"fastas/{gene}.fasta", f">{gene}\nATGGCTTAA\n"),
            "mutations": (f"mappings/mutations/{gene}_mutations.csv", "mutant\nG4A\nT6C\n"),
            "msa": (f"MSA/{gene}.msa.a2m", f">{gene}\nMA\n"),
            "codon_msa": (f"CodonMSA/{gene}.codon.msa.fasta", f">{gene}\nATGGCT\n"),
        }
        expected[gene] = {}
        for artifact, (relative, content) in files.items():
            if generate_protein and artifact == "msa":
                continue
            if blocked_codon and artifact == "codon_msa":
                content = f">{gene}\n" + "ATG" * 1000 + "\n"
            target = source_root / gene / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content)
            expected[gene][artifact] = str(target.resolve())
        if prebuilt_ev:
            for artifact in ("model_params", "codon_model_params"):
                model = output_root / artifact / f"{gene}.{artifact}"
                model.parent.mkdir(parents=True, exist_ok=True)
                model.write_bytes(ev_params_bytes(21 if artifact == "model_params" else 65))

    params_file = tmp_path / "prebuilt-params.dat"
    params_file.write_text("routing-test stub; not model parameters")
    hardware_file = tmp_path / "hardware.json"
    hardware_file.write_text(json.dumps({"cpus": 20 if automatic_threads else 2, "memory_gib": 8, "gpus": []}))
    command = [
        sys.executable, str(runtime / "mutEffects_controller.py"),
        "--fasta", str(source_root), "--output", str(output_root),
        "--adabmdca-protein-params", str(params_file),
        "--adabmdca-codon-params", str(params_file),
        "--adabmdca-nchains", "64",
        "--adabmdca-device", "cpu",
        "--resource-hardware", str(hardware_file),
        "--msa-memory-gib", "0.5", "--evmutation-memory-gib", "0.5",
    ]
    if not automatic_threads:
        command.extend(["--threads", "1", "--adabmdca-nepochs", "1"])
    if generate_protein:
        database = tmp_path / "database"
        database.mkdir()
        (database / "uniref90.fasta").write_text(">database\nMA\n")
        command.extend(["--db-root", str(database)])
    if backend == "adabmDCA":
        command.append("--adabmdca-only")
    else:
        command.append("--evmutation-only")
        if not prebuilt_ev:
            command.extend(["--plmc-binary", "/bin/true"])
    environment = dict(os.environ)
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": str(Path(controller.__file__).resolve().parents[2]),
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
        "NXF_ANSI_LOG": "false",
        "BFF_EXPECT_THREADS": "10" if automatic_threads else "1",
    })
    if prebuilt_ev:
        environment["BFF_EXPECT_EV_PARAMS"] = "1"
    if generate_protein:
        environment["PYTHONPATH"] = str(runtime) + os.pathsep + environment["PYTHONPATH"]
        environment["BFF_READY_SIDE_MARKER"] = str(tmp_path / "codon-scoring-started")
    result = subprocess.run(
        command, cwd=runtime, env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120,
    )
    assert result.returncode == int(blocked_codon), result.stdout
    manifest = json.loads((output_root / ".evmutation_manifest.json").read_text())
    assert manifest["input_files"] == expected
    for gene, inputs in expected.items():
        for side, artifact in (("protein", "msa"), ("codon", "codon_msa")):
            if blocked_codon and side == "codon":
                assert not (output_root / gene / backend / f"{gene}.{side}.tsv").exists()
                continue
            observed = json.loads((output_root / gene / backend / f"{gene}.{side}.tsv").read_text())
            assert observed == {
                "fasta": inputs["fasta"], "mutations": inputs["mutations"],
                artifact: Path(inputs[artifact]).read_text() if inputs.get(artifact) else f">{gene}\nMA\n",
            }
            if backend == "adabmDCA":
                complete = json.loads(
                    (output_root / gene / backend / f"{gene}.{side}.complete.json").read_text()
                )
                assert complete["device"] == "cpu"
                assert complete["threads"] == (10 if automatic_threads else 1)
                assert complete["settings"]["model"] == "pseudoDCA"
                assert complete["settings"]["nepochs"] == (None if automatic_threads else 1)
            else:
                marker = output_root / gene / backend / f"{gene}.{side}.routing.json"
                assert json.loads(marker.read_text())["fingerprint"] == manifest["ev_fingerprints"][gene][side]
                plan = json.loads((output_root / "resource_plans" / f"{gene}.{side}.evmutation.plan.json").read_text())
                assert plan["threads"] == (10 if automatic_threads else 1)
                assert plan["memory_gib"] > 0
        if prebuilt_ev:
            for artifact in ("model_params", "codon_model_params"):
                assert (output_root / artifact / f"{gene}.{artifact}").read_bytes() == ev_params_bytes(21 if artifact == "model_params" else 65)
    tasks = list((runtime / "work").glob("*/*/.command.sh"))
    assert len(tasks) == len(expected) * (2 if blocked_codon else 4 + int(generate_protein))
    for row in (runtime / "trace.tsv").read_text().splitlines()[1:]:
        task_name, cpus = row.split("\t")
        if task_name.startswith("run_"):
            assert int(cpus) == (10 if automatic_threads else 1)
    assert sum("msa_generation_pipeline" in task.read_text() for task in tasks) == int(generate_protein)
    assert all("codon_msa_pipeline.py" not in task.read_text() for task in tasks)
    if blocked_codon:
        failures = manifest["resource_errors"]
        assert len(failures) == 1
        failure = failures[0]
        assert (failure["gene"], failure["side"], failure["backend"]) == ("NPM1", "codon", "evmutation")
        assert result.stdout.count(failure["message"]) == 1
        assert result.stdout.index(failure["message"]) > result.stdout.rindex("Submitted process >")
        assert not (output_root / "resource_plans/NPM1.codon.evmutation.plan.json").exists()
        assert "Closing CacheDB done" in (runtime / ".nextflow.log").read_text()
        return
    if backend == "EVmutation":
        cached_files = {
            path: path.stat().st_mtime_ns
            for path in output_root.glob("*/EVmutation/*")
            if path.is_file()
        }
        repeated = subprocess.run(
            command, cwd=runtime, env=environment, capture_output=True, text=True, timeout=120,
        )
        assert repeated.returncode == 0, repeated.stdout + repeated.stderr
        assert "Submitted process >" not in repeated.stdout
        assert set((runtime / "work").glob("*/*/.command.sh")) == set(tasks)
        assert all(path.stat().st_mtime_ns == modified for path, modified in cached_files.items())
        if prebuilt_ev:
            forced_environment = {**environment, "BFF_EXPECT_FORCED_CODON": "1"}
            forced = subprocess.run(
                [*command, "--codon-msa", str(source_root)], cwd=runtime,
                env=forced_environment, capture_output=True, text=True, timeout=120,
            )
            assert forced.returncode == 0, forced.stdout + forced.stderr
            assert len(list((runtime / "work").glob("*/*/.command.sh"))) == len(tasks) + 2 * len(expected)
            for gene in expected:
                marker = output_root / gene / "EVmutation" / f"{gene}.codon.routing.json"
                assert json.loads(marker.read_text())["fingerprint"] != manifest["ev_fingerprints"][gene]["codon"]
            restored = subprocess.run(
                command, cwd=runtime, env=environment, capture_output=True, text=True, timeout=120,
            )
            assert restored.returncode == 0, restored.stdout + restored.stderr
            assert len(list((runtime / "work").glob("*/*/.command.sh"))) == len(tasks) + 4 * len(expected)
            for gene in expected:
                marker = output_root / gene / "EVmutation" / f"{gene}.codon.routing.json"
                assert json.loads(marker.read_text())["fingerprint"] == manifest["ev_fingerprints"][gene]["codon"]
                protein = output_root / gene / "EVmutation" / f"{gene}.protein.tsv"
                assert protein.stat().st_mtime_ns == cached_files[protein]


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test routing",
)
@pytest.mark.parametrize("generate_msa,gpu_oom,scenario", [
    (False, False, None), (True, False, None), (False, True, None),
    (False, False, "protein"), (False, False, "codon"),
    (True, False, "stream"), (False, False, "completed"), (False, False, "skip"),
    (False, False, "resource_blocked"),
])
def test_resource_plans_route_each_side_and_retry_oom(tmp_path, generate_msa, gpu_oom, scenario):
    worker_threads = 2 if generate_msa and scenario is None else 1
    runtime = tmp_path / "runtime"
    binary_dir = runtime / "bin"
    binary_dir.mkdir(parents=True)
    for name in ("main.nf", "resource_planner.py", "adabmdca_task.py"):
        shutil.copy2(Path(controller.NEXTFLOW_SCRIPT).parent / name, binary_dir / name)
    (runtime / "adabmdca_pipeline.py").write_text(textwrap.dedent('''\
        import argparse
        import json
        import os
        from pathlib import Path
        import sys

        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', type=Path)
        parser.add_argument('--msa', type=Path)
        parser.add_argument('--codon-msa', type=Path)
        parser.add_argument('--adabmdca-device')
        parser.add_argument('--adabmdca-nchains', type=int)
        parser.add_argument('--adabmdca-nepochs', type=int)
        parser.add_argument('--adabmdca-model')
        parser.add_argument('--output', type=Path)
        parser.add_argument('--protein-params', type=Path)
        parser.add_argument('--codon-params', type=Path)
        parser.add_argument('--skip-codon', action='store_true')
        parser.add_argument('--score-missense-codon', action='store_true')
        args, unused = parser.parse_known_args()
        gene = args.fasta.read_text().splitlines()[0][1:]
        side = 'protein' if args.msa else 'codon'
        msa = args.msa or args.codon_msa
        assert msa.exists()
        assert args.adabmdca_nchains == 64
        assert args.adabmdca_nepochs == 1
        assert args.adabmdca_model == 'bmDCA'
        for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BFF_TASK_THREADS'):
            assert os.environ[variable] == os.environ['BFF_EXPECT_THREADS']
        if os.environ.get('BFF_EXPECT_PROTEIN_ONLY'):
            assert side == 'protein' and args.skip_codon
        if os.environ.get('BFF_EXPECT_CODON_ONLY'):
            assert side == 'codon' and args.score_missense_codon
        if side == 'codon' and os.environ.get('BFF_READY_SIDE_MARKER'):
            Path(os.environ['BFF_READY_SIDE_MARKER']).touch()
        if args.adabmdca_device == 'cuda':
            print('torch.OutOfMemoryError: CUDA out of memory.', file=sys.stderr)
            sys.exit(1)
        output = args.output / gene / 'adabmDCA'
        output.mkdir(parents=True)
        (output / f'{gene}.{side}.tsv').write_text('mutant\\tscore\\nG4A\\t0.25\\n')
        params_path = args.protein_params or args.codon_params
        if params_path.is_dir():
            params_path = params_path / f'{gene}.{side}_adabm_params'
        params_path.write_text('stub parameters')
    '''))
    source = tmp_path / "inputs"
    source.mkdir()
    fasta = source / "NPM1.fasta"
    fasta.write_text(">NPM1\nATGGCTTAA\n")
    mutations = source / "NPM1_mutations.csv"
    mutations.write_text("mutant\nG4A\n")
    protein_msa = source / "NPM1.msa.a2m"
    protein_msa.write_text(">NPM1\nMA\n")
    codon_msa = source / "NPM1.codon.msa.fasta"
    codon_msa.write_text(">NPM1\nATGGCT\n")
    inputs = {"fasta": str(fasta), "mutations": str(mutations)}
    manifest = {"input_files": {"NPM1": inputs}}
    available_sides = set() if generate_msa else {"protein", "codon"}
    needed_sides = {"protein", "codon"}
    if scenario in {"protein", "skip"}:
        available_sides = needed_sides = {"protein"}
    elif scenario in {"codon", "completed", "resource_blocked"}:
        available_sides = needed_sides = {"codon"}
    elif scenario == "stream":
        available_sides = {"codon"}
    for side, artifact, path in (("protein", "msa", protein_msa), ("codon", "codon_msa", codon_msa)):
        if side in available_sides:
            inputs[artifact] = str(path)
            manifest[artifact] = ["NPM1"]
    if scenario:
        forced_codon = scenario in {"codon", "skip"}
        manifest["routing"] = {"NPM1": {
            "protein": not forced_codon, "codon": scenario != "protein",
            "mode": "codon" if forced_codon else "protein" if scenario == "protein" else "auto",
            "score_missense_codon": forced_codon, "classes": ["missense"],
        }}
    if scenario == "completed":
        manifest["adabmdca_protein"] = ["NPM1"]
    if scenario == "resource_blocked":
        manifest["resource_errors"] = [{
            "gene": "NPM1", "side": "protein", "backend": "adabmdca",
            "message": "NPM1/protein: prebuilt resource plan blocked", "stage": "resource_planning",
        }]
    manifest_file = tmp_path / "manifest.json"
    manifest_file.write_text(json.dumps(manifest))
    package = runtime / "biofeaturefactory"
    core = package / "core"
    core.mkdir(parents=True)
    for target in (package / "__init__.py", core / "__init__.py"):
        target.write_text("")
    for module, side in (("msa_generation_pipeline", "protein"), ("codon_msa_pipeline", "codon")):
        stub = textwrap.dedent('''\
            import argparse
            from pathlib import Path

            parser = argparse.ArgumentParser()
            parser.add_argument('--fasta', type=Path)
            parser.add_argument('--threads', type=int)
            args, unused = parser.parse_known_args()
            import os
            assert args.threads == int(os.environ['BFF_EXPECT_THREADS'])
            gene = args.fasta.read_text().splitlines()[0][1:]
        ''')
        if side == "protein":
            if scenario == "stream":
                stub += textwrap.dedent('''\
                    import os
                    import time

                    marker = Path(os.environ['BFF_READY_SIDE_MARKER'])
                    deadline = time.monotonic() + 20
                    while not marker.exists() and time.monotonic() < deadline:
                        time.sleep(0.05)
                    assert marker.exists(), 'Ready codon MSA scoring waited for protein MSA generation'
                ''')
            stub += "Path(f'{gene}.msa.a2m').write_text(f'>{gene}\\nMA\\n')\n"
            stub += "Path(f'{gene}.msa.stats.json').write_text('{}')\n"
        else:
            stub += "Path(f'{gene}.codon.msa.fasta').write_text(f'>{gene}\\nATGGCT\\n')\n"
            stub += "Path(f'{gene}.codon.msa.manifest.tsv').write_text('stub')\n"
            stub += "Path(f'{gene}.codon.msa.stats.json').write_text('{}')\n"
        (core / f"{module}.py").write_text(stub)
    config_file = tmp_path / "resources.json"
    hardware = {"cpus": 2, "memory_gib": 8, "gpus": []}
    if gpu_oom:
        hardware["gpus"] = [{"uuid": "GPU-test", "memory_gib": 80}]
    config_file.write_text(json.dumps({
        "hardware": hardware, "threads": worker_threads, "device": "auto",
        "lease_dir": str(tmp_path / "gpu_leases"), "gpu_wait_timeout": 10,
        "settings": {**ADABMDCA_TEST_SETTINGS, "skip_codon": scenario == "skip"},
        "routing": manifest.get("routing", {}),
        "overrides": {
            f"NPM1.{side}": {
                "cpu_memory_gib": 1, "gpu_host_memory_gib": 0.25, "gpu_memory_gib": 1,
            }
            for side in ("protein", "codon")
        },
    }))
    executor_config = runtime / "nextflow.config"
    executor_config.write_text(
        f"params.msa_cpus = {worker_threads}\nparams.msa_memory = '512 MB'\n"
        "params.gpu_slots = 1\n"
        "trace { enabled = true; file = 'trace.tsv'; fields = 'name,cpus,memory,workdir,status' }\n"
    )
    output = tmp_path / "results"
    command = [
        "nextflow", "run", str(binary_dir / "main.nf"),
        "--fasta", str(fasta), "--mutations", str(mutations),
        "--manifest", str(manifest_file), "--skip_evmutation", "true",
        "--resource_config", str(config_file),
        "--resource_executor", "local",
        "--output_dir", str(output),
    ]
    if "protein" in needed_sides - available_sides:
        command.extend(["--uniref90_db", str(source)])
    if "codon" in needed_sides - available_sides:
        command.extend(["--db_root", str(source)])
    if scenario == "skip":
        command.extend(["--skip_codon_adabmdca", "true"])
    if scenario == "resource_blocked":
        command.extend(["--resource_errors", str(tmp_path / "resource-errors.json")])
    environment = dict(os.environ)
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": str(runtime),
        "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_ANSI_LOG": "false",
        "BFF_EXPECT_THREADS": str(worker_threads),
    })
    if scenario in {"protein", "skip"}:
        environment["BFF_EXPECT_PROTEIN_ONLY"] = "1"
    if scenario == "codon":
        environment["BFF_EXPECT_CODON_ONLY"] = "1"
    if scenario == "stream":
        environment["BFF_READY_SIDE_MARKER"] = str(tmp_path / "codon-scoring-started")
    if gpu_oom:
        inventory = binary_dir / "nvidia-smi"
        inventory.write_text("#!/bin/sh\nprintf 'GPU-test, 81920\\n'\n")
        inventory.chmod(0o755)
        environment["PATH"] = str(binary_dir) + os.pathsep + environment["PATH"]
    result = subprocess.run(
        command, cwd=runtime, env=environment, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    if scenario == "resource_blocked":
        assert json.loads((tmp_path / "resource-errors.json").read_text()) == []
        assert not (output / "resource_plans/NPM1.protein.plan.json").exists()
        assert not (output / "NPM1/adabmDCA/NPM1.protein.tsv").exists()
    tasks = list((runtime / "work").glob("*/*/.command.sh"))
    gpu_tasks = [task for task in tasks if "'--device' 'cuda'" in task.read_text()]
    cpu_tasks = [task for task in tasks if "'--device' 'cpu'" in task.read_text()]
    assert len(gpu_tasks) == (2 if gpu_oom else 0)
    assert len(cpu_tasks) == len(needed_sides)
    assert {task.parent for task in gpu_tasks}.isdisjoint({task.parent for task in cpu_tasks})
    for side in needed_sides:
        plan = json.loads((output / "resource_plans" / f"NPM1.{side}.plan.json").read_text())
        assert plan["device"] == ("cuda" if gpu_oom else "cpu")
        complete = json.loads((output / "NPM1/adabmDCA" / f"NPM1.{side}.complete.json").read_text())
        assert complete["device"] == "cpu"
        assert (output / f"adabmdca_{side}_params" / f"NPM1.{side}_adabm_params").exists()
    trace = (runtime / "trace.tsv").read_text()
    for row in trace.splitlines()[1:]:
        task_name, cpus, unused = row.split("\t", 2)
        if not task_name.startswith("plan_adabmdca"):
            assert int(cpus) == worker_threads
    assert trace.count("run_adabmdca_cpu (") == len(needed_sides)
    assert trace.count("generate_protein_msa (") == int("protein" in needed_sides - available_sides)
    assert trace.count("generate_codon_msa (") == int("codon" in needed_sides - available_sides)
    if gpu_oom:
        assert "256 MB" in trace
        assert "1 GB" in trace
        assert trace.count("run_adabmdca_gpu (") == 2
    if scenario == "resource_blocked":
        strict_runtime = tmp_path / "strict-runtime"
        strict_runtime.mkdir()
        strict = subprocess.run(
            command[:-2], cwd=strict_runtime, env=environment, capture_output=True, text=True, timeout=90,
        )
        assert strict.returncode != 0
        assert "Manifest contains resource planning errors" in strict.stdout + strict.stderr
        assert "Submitted process >" not in strict.stdout + strict.stderr


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test routing",
)
def test_nextflow_queues_seven_jobs_with_two_gpu_and_three_cpu_slots(tmp_path):
    runtime = tmp_path / "runtime"
    binary_dir = runtime / "bin"
    binary_dir.mkdir(parents=True)
    for name in ("main.nf", "resource_planner.py", "adabmdca_task.py"):
        shutil.copy2(Path(controller.NEXTFLOW_SCRIPT).parent / name, binary_dir / name)
    (runtime / "adabmdca_pipeline.py").write_text(textwrap.dedent('''\
        import argparse
        import fcntl
        import json
        import os
        from pathlib import Path
        import time

        parser = argparse.ArgumentParser()
        parser.add_argument('--gene')
        parser.add_argument('--msa')
        parser.add_argument('--output', type=Path)
        parser.add_argument('--protein-params', type=Path)
        parser.add_argument('--codon-params', type=Path)
        parser.add_argument('--adabmdca-device')
        args, unused = parser.parse_known_args()
        events = Path(os.environ['BFF_TEST_EVENTS'])
        release = events.with_suffix('.release')
        gpu_ready = events.with_suffix('.gpu-ready')
        side = 'protein' if args.msa else 'codon'
        record = {
            'gene': args.gene, 'device': args.adabmdca_device,
            'gpu': os.environ.get('CUDA_VISIBLE_DEVICES'), 'side': side,
        }
        with events.open('a+') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            handle.seek(0)
            previous = [json.loads(line) for line in handle]
            handle.write(json.dumps({**record, 'event': 'start'}) + '\\n')
            handle.flush()
            if args.adabmdca_device == 'cuda' and sum(
                event['event'] == 'start' and event['device'] == 'cuda' for event in previous
            ) == 1:
                gpu_ready.touch()
            if sum(event['event'] == 'start' for event in previous) == 4:
                release.touch()
            fcntl.flock(handle, fcntl.LOCK_UN)
        deadline = time.monotonic() + 30
        while not release.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        assert release.exists(), 'Five runnable jobs did not start together'
        time.sleep(0.3)
        output = args.output / args.gene / 'adabmDCA'
        output.mkdir(parents=True)
        (output / f'{args.gene}.{side}.tsv').write_text('mutant\\tscore\\nG4A\\t0.25\\n')
        params = args.protein_params or args.codon_params
        (params / f'{args.gene}.{side}_adabm_params').write_text('stub parameters')
        with events.open('a') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            handle.write(json.dumps({**record, 'event': 'finish'}) + '\\n')
    '''))
    inventory = binary_dir / "nvidia-smi"
    inventory.write_text("#!/bin/sh\nprintf 'GPU-first, 81920\\nGPU-second, 81920\\n'\n")
    inventory.chmod(0o755)
    jobs = []
    for index in range(7):
        gene = f"GENE{index}"
        side = "protein" if index % 2 == 0 else "codon"
        device = "cuda" if index < 3 else "cpu"
        fasta = tmp_path / f"{gene}.fasta"
        fasta.write_text(f">{gene}\nATGGCTTAA\n")
        msa = tmp_path / f"{gene}.msa"
        msa.write_text(f">{gene}\n" + ("MA\n" if side == "protein" else "ATGGCT\n"))
        mutations = tmp_path / f"{gene}.csv"
        mutations.write_text("mutant\nG4A\n")
        plan = {
            "gene": gene, "side": side, "device": device, "threads": 1,
            "cpu_memory_gib": 1, "gpu_host_memory_gib": 0.25,
            "gpu_memory_gib": 1, "params": None, "fingerprint": gene + "-test",
            "settings": ADABMDCA_TEST_SETTINGS,
            "eligible_gpu_uuids": ["GPU-first", "GPU-second"],
            "lease_dir": str(tmp_path / "leases"), "gpu_wait_timeout": 20,
        }
        plan_file = tmp_path / f"{gene}.plan.json"
        plan_file.write_text(json.dumps(plan))
        jobs.append({
            "gene": gene, "side": side, "fasta": str(fasta), "msa": str(msa),
            "mutations": str(mutations), "plan_file": str(plan_file), "plan": plan,
        })
    tasks_file = tmp_path / "tasks.json"
    tasks_file.write_text(json.dumps(jobs))
    harness = binary_dir / "schedule.nf"
    harness.write_text(textwrap.dedent('''\
        nextflow.enable.dsl = 2
        include { run_adabmdca_gpu; run_adabmdca_cpu } from './main.nf'

        workflow {
            def gpu_ready = Channel.watchPath(params.gpu_ready, 'create').first()
            def definitions = new groovy.json.JsonSlurper().parse(file(params.tasks).toFile())
            def routes = Channel.fromList(definitions)
                .map { entry ->
                    tuple(entry.gene, entry.side, file(entry.fasta), file(entry.msa),
                          file(entry.mutations), file(entry.plan_file), entry.plan)
                }
                .branch { gene, side, fasta, msa, mutations, plan_file, plan ->
                    gpu: plan.device == 'cuda'
                    cpu: plan.device == 'cpu'
                }
            run_adabmdca_gpu(routes.gpu)
            def cpu_ready = routes.cpu.combine(gpu_ready)
                .map { gene, side, fasta, msa, mutations, plan_file, plan, ready_file ->
                    tuple(gene, side, fasta, msa, mutations, plan_file, plan)
                }
            run_adabmdca_cpu(cpu_ready)
        }
    '''))
    (runtime / "nextflow.config").write_text(
        "params.gpu_slots = 2\n"
        "executor { cpus = 5; memory = '3.5 GB' }\n"
    )
    events_file = tmp_path / "events.jsonl"
    environment = dict(os.environ)
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1", "BFF_TEST_EVENTS": str(events_file),
        "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_ANSI_LOG": "false",
        "PATH": str(binary_dir) + os.pathsep + environment["PATH"],
    })
    result = subprocess.run(
        ["nextflow", "run", str(harness), "--tasks", str(tasks_file),
         "--gpu_ready", str(events_file.with_suffix(".gpu-ready")),
         "--output_dir", str(tmp_path / "results")],
        cwd=runtime, env=environment, capture_output=True, text=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    events = [json.loads(line) for line in events_file.read_text().splitlines()]
    initial = events[:5]
    assert all(event["event"] == "start" for event in initial)
    assert [event["device"] for event in initial].count("cuda") == 2
    assert [event["device"] for event in initial].count("cpu") == 3
    assert len({event["gpu"] for event in initial if event["device"] == "cuda"}) == 2
    assert all(event["gpu"] == "" for event in events if event["device"] == "cpu")
    active_gpus = set()
    active_cpus = set()
    completed = set()
    for event in events:
        if event["device"] == "cuda":
            if event["event"] == "start":
                assert event["gpu"] not in active_gpus
                active_gpus.add(event["gpu"])
                assert len(active_gpus) <= 2
            else:
                active_gpus.remove(event["gpu"])
        elif event["event"] == "start":
            active_cpus.add(event["gene"])
            assert len(active_cpus) <= 3
        else:
            active_cpus.remove(event["gene"])
        assert 0.25 * len(active_gpus) + len(active_cpus) <= 3.5
        if event["event"] == "finish":
            completed.add(event["gene"])
    assert len(completed) == 7
    assert not active_gpus


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test admission",
)
@pytest.mark.parametrize("configured", [True, False])
def test_nextflow_queues_heavy_evmutation_until_shared_ram_is_free(tmp_path, configured):
    from biofeaturefactory.mutation_effects.bin.plmc_resources import plan_evmutation_task

    runtime = tmp_path / "runtime"
    binary_dir = runtime / "bin"
    binary_dir.mkdir(parents=True)
    for name in ("main.nf", "plmc_resources.py", "resource_planner.py", "codon_encoding.py"):
        shutil.copy2(Path(controller.NEXTFLOW_SCRIPT).parent / name, binary_dir / name)
    jobs = []
    config = {
        "hardware": {"cpus": 8, "memory_gib": 16, "gpus": []},
        "threads": 1, "evmutation": {"memory_gib": None},
    }
    for gene, side, sequence in (("LIGHT1", "protein", "MA"), ("LIGHT2", "protein", "MA"), ("HEAVY", "codon", "ATG" * 120)):
        fasta = tmp_path / f"{gene}.fasta"
        fasta.write_text(f">{gene}\nATGGCTTAA\n")
        msa = tmp_path / f"{gene}.msa.fasta"
        msa.write_text(f">{gene}\n{sequence}\n")
        mutations = tmp_path / f"{gene}.csv"
        mutations.write_text("mutant\nG4A\n")
        plan_config = config if configured else {"threads": 1, "evmutation": {"memory_gib": None}}
        plan = plan_evmutation_task(gene, side, msa, plan_config)
        jobs.append({
            "gene": gene, "side": side, "fasta": str(fasta), "msa": str(msa),
            "mutations": str(mutations), "memory_gib": plan["memory_gib"],
        })
    budget = jobs[-1]["memory_gib"] + 0.6
    assert sum(job["memory_gib"] for job in jobs[:2]) + 0.5 < budget
    assert jobs[0]["memory_gib"] > 0.6
    tasks_file = tmp_path / "tasks.json"
    tasks_file.write_text(json.dumps(jobs))
    config_file = tmp_path / "resources.json"
    config_file.write_text(json.dumps(config))
    events_file = tmp_path / "events.jsonl"
    ready = tmp_path / "light-ready"
    planned = tmp_path / "heavy-planned"
    (runtime / "evmutation_pipeline.py").write_text(textwrap.dedent('''\
        import argparse
        import fcntl
        import json
        import os
        from pathlib import Path
        import time

        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', type=Path)
        parser.add_argument('--msa')
        args, unused = parser.parse_known_args()
        gene = args.fasta.read_text().splitlines()[0][1:]
        side = 'protein' if args.msa else 'codon'
        for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            assert os.environ[variable] == '1'
        events = Path(os.environ['BFF_TEST_EVENTS'])
        with events.open('a+') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            handle.seek(0)
            previous = [json.loads(line) for line in handle]
            handle.write(json.dumps({'gene': gene, 'event': 'start'}) + '\\n')
            handle.flush()
            if side == 'protein' and any(event['event'] == 'start' for event in previous):
                Path(os.environ['BFF_LIGHT_READY']).touch()
        if side == 'protein':
            planned = Path(os.environ['BFF_HEAVY_PLANNED'])
            deadline = time.monotonic() + 30
            while not planned.exists() and time.monotonic() < deadline:
                time.sleep(0.05)
            assert planned.exists(), 'Heavy EV plan did not finish while light jobs held RAM'
            time.sleep(1)
        Path(f'{gene}.{side}.tsv').write_text('mutant\\tscore\\nG4A\\t0.25\\n')
        suffix = 'model_params' if side == 'protein' else 'codon_model_params'
        Path(f'{gene}.{suffix}').write_text('stub parameters')
        with events.open('a') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            handle.write(json.dumps({'gene': gene, 'event': 'finish'}) + '\\n')
    '''))
    harness = binary_dir / "schedule_ev.nf"
    harness.write_text(textwrap.dedent('''\
        nextflow.enable.dsl = 2
        include { plan_evmutation; run_protein_evmutation; run_codon_evmutation } from './main.nf'

        workflow {
            def ready = Channel.watchPath(params.light_ready, 'create').first()
            def definitions = new groovy.json.JsonSlurperClassic().parse(file(params.tasks).toFile())
            def inputs = Channel.fromList(definitions)
                .map { entry ->
                    tuple(entry.gene, entry.side, file(entry.fasta), file(entry.msa),
                          file(entry.mutations), [skip_codon: true, model_params: null, fingerprint: null])
                }
                .branch { gene, side, fasta, msa, mutations, options ->
                    light: side == 'protein'
                    heavy: side == 'codon'
                }
            def heavy_ready = inputs.heavy.combine(ready)
                .map { gene, side, fasta, msa, mutations, options, ready_file ->
                    tuple(gene, side, fasta, msa, mutations, options)
                }
            def config = params.resource_config ? file(params.resource_config) : []
            def planned = plan_evmutation(inputs.light.mix(heavy_ready), config)
                .map { gene, side, fasta, msa, mutations, options, plan_file ->
                    def plan = new groovy.json.JsonSlurperClassic().parse(plan_file.toFile())
                    if (side == 'codon')
                        new File(params.heavy_planned).text = 'ready'
                    tuple(gene, side, fasta, msa, mutations, options, plan_file, plan)
                }
                .branch { gene, side, fasta, msa, mutations, options, plan_file, plan ->
                    protein: side == 'protein'
                    codon: side == 'codon'
                }
            run_protein_evmutation(planned.protein)
            run_codon_evmutation(planned.codon)
        }
    '''))
    (runtime / "nextflow.config").write_text(
        "params.evmutation_cpus = 1\n"
        f"executor {{ cpus = 8; memory = '{budget:.6f} GB' }}\n"
        "trace { enabled = true; file = 'trace.tsv'; raw = true; fields = 'name,cpus,memory' }\n"
    )
    command = [
        "nextflow", "run", str(harness), "--tasks", str(tasks_file),
        "--light_ready", str(ready), "--heavy_planned", str(planned),
        "--output_dir", str(tmp_path / "results"),
    ]
    if configured:
        command.extend(["--resource_config", str(config_file)])
    environment = {
        **os.environ, "PYTHONDONTWRITEBYTECODE": "1", "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true", "NXF_ANSI_LOG": "false",
        "BFF_TEST_EVENTS": str(events_file), "BFF_LIGHT_READY": str(ready),
        "BFF_HEAVY_PLANNED": str(planned),
        "PYTHONPATH": str(Path(controller.__file__).resolve().parents[2]),
    }
    result = subprocess.run(command, cwd=runtime, env=environment, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    events = [json.loads(line) for line in events_file.read_text().splitlines()]
    assert {event["gene"] for event in events[:2]} == {"LIGHT1", "LIGHT2"}
    assert all(event["event"] == "start" for event in events[:2])
    heavy_start = events.index({"gene": "HEAVY", "event": "start"})
    assert all(events.index({"gene": gene, "event": "finish"}) < heavy_start for gene in ("LIGHT1", "LIGHT2"))
    active = set()
    memory = {job["gene"]: job["memory_gib"] for job in jobs}
    for event in events:
        if event["event"] == "start":
            active.add(event["gene"])
        else:
            active.remove(event["gene"])
        assert sum(memory[gene] for gene in active) <= budget
    assert not active and len(events) == 6
    for row in (runtime / "trace.tsv").read_text().splitlines()[1:]:
        name, cpus, requested = row.split("\t")
        assert int(cpus) == 1
        if name.startswith("run_"):
            gene = name.split("(")[1].split()[0]
            assert int(requested) / 1024 ** 3 == pytest.approx(memory[gene], abs=1e-6)


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test deferred errors",
)
@pytest.mark.parametrize("backend", ["evmutation", "adabmdca"])
def test_nextflow_defers_generated_msa_resource_errors(tmp_path, backend):
    runtime = tmp_path / "runtime"
    binary_dir = runtime / "bin"
    binary_dir.mkdir(parents=True)
    for name in ("main.nf", "plmc_resources.py", "resource_planner.py", "codon_encoding.py", "adabmdca_task.py"):
        shutil.copy2(Path(controller.NEXTFLOW_SCRIPT).parent / name, binary_dir / name)
    (runtime / f"{backend}_pipeline.py").write_text(textwrap.dedent('''\
        import argparse
        import os
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument('--fasta', type=Path)
        parser.add_argument('--msa', type=Path)
        parser.add_argument('--codon-msa', type=Path)
        parser.add_argument('--output', type=Path)
        parser.add_argument('--protein-params', type=Path)
        args, unused = parser.parse_known_args()
        assert args.msa and not args.codon_msa, 'Unschedulable codon task launched'
        gene = args.fasta.read_text().splitlines()[0][1:]
        if args.protein_params:
            output = args.output / gene / 'adabmDCA'
            output.mkdir(parents=True)
            params_path = args.protein_params
            if params_path.is_dir():
                params_path = params_path / f'{gene}.protein_adabm_params'
        else:
            output = args.output
            params_path = output / f'{gene}.model_params'
        (output / f'{gene}.protein.tsv').write_text('mutant\\tscore\\nG4A\\t0.25\\n')
        params_path.write_text('stub parameters')
        Path(os.environ['BFF_HEALTHY_FINISHED']).touch()
    '''))
    core = runtime / "biofeaturefactory" / "core"
    core.mkdir(parents=True)
    (core.parent / "__init__.py").write_text("")
    (core / "__init__.py").write_text("")
    (core / "codon_msa_pipeline.py").write_text(textwrap.dedent('''\
        import os
        from pathlib import Path
        import time

        ready = Path(os.environ['BFF_HEALTHY_FINISHED'])
        deadline = time.monotonic() + 20
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        assert ready.exists(), 'Ready protein scoring did not run during MSA generation'
        Path('GENE.codon.msa.fasta').write_text('>ORF\\n' + 'ATG' * 1000 + '\\n')
        Path('GENE.codon.msa.manifest.tsv').write_text('stub')
        Path('GENE.codon.msa.stats.json').write_text('{}')
    '''))
    fasta = tmp_path / "GENE.fasta"
    fasta.write_text(">GENE\nATGGCTTAA\n")
    mutations = tmp_path / "GENE.csv"
    mutations.write_text("mutant\nG4A\nT6C\n")
    msa = tmp_path / "GENE.msa.a2m"
    msa.write_text(">GENE\nMA\n")
    manifest_file = tmp_path / "manifest.json"
    manifest_file.write_text(json.dumps({
        "input_files": {"GENE": {"fasta": str(fasta), "mutations": str(mutations), "msa": str(msa)}},
    }))
    config_file = tmp_path / "resources.json"
    config_file.write_text(json.dumps({
        "hardware": {"cpus": 4, "memory_gib": 8, "gpus": []},
        "threads": 1, "device": "auto", "settings": ADABMDCA_TEST_SETTINGS,
        "lease_dir": str(tmp_path / "leases"), "gpu_wait_timeout": 10,
        "evmutation": {"memory_gib": None},
    }))
    (runtime / "nextflow.config").write_text(
        "params.msa_cpus = 1\nparams.msa_memory = '512 MB'\nparams.evmutation_cpus = 1\n"
        "executor { cpus = 4; memory = '8 GB' }\n"
        "trace { enabled = true; file = 'trace.tsv'; fields = 'name,status' }\n"
    )
    output = tmp_path / "results"
    report = tmp_path / "resource-errors.json"
    command = [
        "nextflow", "run", str(binary_dir / "main.nf"),
        "--fasta", str(fasta), "--mutations", str(mutations),
        "--manifest", str(manifest_file), "--resource_config", str(config_file),
        "--output_dir", str(output), "--db_root", str(tmp_path),
        "--skip_adabmdca" if backend == "evmutation" else "--skip_evmutation", "true",
    ]
    if backend == "evmutation":
        command.extend(["--plmc_binary", "/bin/true"])
    report.write_text("[]")
    command.extend(["--resource_errors", str(report)])
    environment = {
        **os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": str(runtime),
        "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_ANSI_LOG": "false",
        "BFF_HEALTHY_FINISHED": str(tmp_path / "healthy-finished"),
    }
    result = subprocess.run(command, cwd=runtime, env=environment, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    directory = "EVmutation" if backend == "evmutation" else "adabmDCA"
    assert (output / "GENE" / directory / "GENE.protein.tsv").is_file()
    assert not (output / "GENE" / directory / "GENE.codon.tsv").exists()
    plan_name = "GENE.codon.evmutation.plan.json" if backend == "evmutation" else "GENE.codon.plan.json"
    blocked = json.loads((output / "resource_plans" / plan_name).read_text())
    assert blocked["device"] == "blocked"
    failure = blocked["resource_error"]
    assert (failure["gene"], failure["side"], failure["backend"], failure["stage"]) == (
        "GENE", "codon", backend, "resource_planning",
    )
    assert json.loads(report.read_text()) == [failure]
    assert failure["message"] not in result.stdout + result.stderr
    assert "Failed to invoke" not in result.stdout + result.stderr
    assert "Closing CacheDB done" in (runtime / ".nextflow.log").read_text()
    rows = (runtime / "trace.tsv").read_text().splitlines()[1:]
    assert len(rows) == 4
    assert all(row.endswith("\tCOMPLETED") for row in rows)
