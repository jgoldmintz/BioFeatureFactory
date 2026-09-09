"""Keep schedulable models running and summarize resource failures afterward."""

import json
from pathlib import Path
import sys

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller


def _write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


@pytest.fixture
def panel(tmp_path):
    def make(genes, memory_gib=32):
        source = tmp_path / "input"
        for gene, (sites, mutations) in genes.items():
            orf = "ATG" + "GCT" * (sites - 1)
            protein = "M" + "A" * (sites - 1)
            _write(source / gene / "fastas" / f"{gene}.fasta", f">ORF\n{orf}\n")
            _write(source / gene / "mappings" / "mutations" / f"{gene}_mutations.csv", "mutant\n" + mutations)
            _write(source / gene / "MSA" / f"{gene}.msa.a2m", f">{gene}\n{protein}\n>homolog\n{protein}\n")
            _write(source / gene / "CodonMSA" / f"{gene}.codon.msa.fasta", f">ORF\n{orf}\n>homolog\n{orf}\n")
        hardware = _write(tmp_path / "hardware.json", json.dumps({
            "cpus": 12, "memory_gib": memory_gib, "gpus": [],
        }))
        return {"source": source, "hardware": hardware, "output": tmp_path / "output", "genes": list(genes)}

    return make


def _parse(monkeypatch, inputs, *options, both_backends=False):
    arguments = [
        "mutEffects_controller", "--fasta", str(inputs["source"]),
        "--output", str(inputs["output"]), "--resource-hardware", str(inputs["hardware"]),
    ]
    if not both_backends:
        arguments.append("--evmutation-only")
    monkeypatch.setattr(sys, "argv", [*arguments, *map(str, options)])
    return controller.parse_args()


def _prepare(args, inputs):
    manifest = controller.build_manifest(inputs["genes"], args)
    config, plans = controller.prepare_resource_plans(args, manifest, inputs["genes"])
    return manifest, config, plans


def _deny(*arguments, **keywords):
    pytest.fail("Blocked resource task unexpectedly required a tool or launched Nextflow")


def _capture_nextflow(monkeypatch, capsys, returncode=0, dynamic_errors=()):
    launches = []

    class Process:
        def __init__(self, command, **options):
            launches.append(command)
            manifest_path = Path(command[command.index("--manifest") + 1])
            self.manifest = json.loads(manifest_path.read_text())
            self.report = Path(command[command.index("--resource_errors") + 1])
            assert json.loads(self.report.read_text()) is None

        def wait(self):
            before_completion = capsys.readouterr()
            assert "WARNING:" not in before_completion.err
            assert "ERROR:" not in before_completion.err
            assert "unschedulable" not in before_completion.out
            self.report.write_text(json.dumps(list(dynamic_errors)))
            print("healthy tasks completed")
            return returncode

    monkeypatch.setattr(controller.subprocess, "Popen", Process)
    return launches


def test_blocked_codon_preserves_healthy_protein_other_gene_and_routing(panel, monkeypatch):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n"), "LIGHT": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs)

    manifest, config, plans = _prepare(args, inputs)

    error, = manifest["resource_errors"]
    assert {key: error[key] for key in ("gene", "side", "backend", "stage")} == {
        "gene": "HEAVY", "side": "codon", "backend": "evmutation", "stage": "resource_planning",
    }
    assert "unschedulable" in error["message"]
    assert controller.side_enabled(args, manifest, "HEAVY", "evmutation", "codon")
    assert not controller.side_pending(args, manifest, "HEAVY", "evmutation", "codon")
    assert {(plan["gene"], plan["side"]) for plan in config["evmutation_plans"]} == {
        ("HEAVY", "protein"), ("LIGHT", "protein"),
    }
    tasks = controller.pending_resource_tasks(args, manifest, inputs["genes"], config, plans)
    assert {task["id"] for task in tasks} == {"HEAVY.protein.evmutation", "LIGHT.protein.evmutation"}
    assert config["cpu_allocation"] == {"threads": 6, "concurrent_jobs": 2}


def test_ev_failures_do_not_block_adabmdca_for_same_gene_side(panel, monkeypatch):
    inputs = panel({"GENE": (30, "G4A\nT6C\n")})
    args = _parse(monkeypatch, inputs, "--evmutation-memory-gib", 64, both_backends=True)

    manifest, config, plans = _prepare(args, inputs)

    assert {(error["backend"], error["side"]) for error in manifest["resource_errors"]} == {
        ("evmutation", "protein"), ("evmutation", "codon"),
    }
    assert config["evmutation_plans"] == []
    assert {plan["side"] for plan in plans} == {"protein", "codon"}
    assert all(controller.side_pending(args, manifest, "GENE", "adabmdca", side) for side in ("protein", "codon"))
    assert all(not controller.side_pending(args, manifest, "GENE", "evmutation", side) for side in ("protein", "codon"))
    monkeypatch.setattr(controller.importlib.util, "find_spec", lambda name: object())
    controller.validate_backend_tools(args, manifest, inputs["genes"])


def test_adabmdca_failure_does_not_block_evmutation(panel, monkeypatch):
    inputs = panel({"GENE": (30, "G4A\n")})
    overrides = _write(inputs["source"].parent / "overrides.json", json.dumps({
        "GENE.protein": {"cpu_memory_gib": 100},
    }))
    args = _parse(monkeypatch, inputs, "--resource-overrides", overrides, both_backends=True)

    manifest, config, plans = _prepare(args, inputs)

    error, = manifest["resource_errors"]
    assert error["backend"] == "adabmdca"
    assert error["gene"] == "GENE"
    assert error["side"] == "protein"
    assert len(config["evmutation_plans"]) == 1
    assert plans == []
    assert controller.side_pending(args, manifest, "GENE", "evmutation", "protein")
    assert not controller.side_pending(args, manifest, "GENE", "adabmdca", "protein")


@pytest.mark.parametrize("exception", [ValueError("invalid alignment shape"), OSError("alignment unavailable")])
def test_planning_input_error_is_collected_without_hiding_other_gene(panel, monkeypatch, exception):
    inputs = panel({"BAD": (30, "G4A\n"), "GOOD": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs)
    original = controller.plmc_resources.plan_evmutation_task

    def selectively_fail(gene, *arguments, **keywords):
        if gene == "BAD":
            raise exception
        return original(gene, *arguments, **keywords)

    monkeypatch.setattr(controller.plmc_resources, "plan_evmutation_task", selectively_fail)

    manifest, config, plans = _prepare(args, inputs)

    error, = manifest["resource_errors"]
    assert error["gene"] == "BAD"
    assert str(exception) in error["message"]
    assert [plan["gene"] for plan in config["evmutation_plans"]] == ["GOOD"]


@pytest.mark.parametrize("returncode", [0, 17])
def test_healthy_jobs_finish_before_resource_error_summary(panel, monkeypatch, capsys, returncode):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n"), "LIGHT": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    capsys.readouterr()
    launches = _capture_nextflow(monkeypatch, capsys, returncode=returncode)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == (returncode or 1)
    assert len(launches) == 1
    after_completion = capsys.readouterr()
    assert "healthy tasks completed" in after_completion.out
    assert "HEAVY/codon" in after_completion.err
    assert "unschedulable" in after_completion.err
    command = launches[0]
    manifest = json.loads(Path(command[command.index("--manifest") + 1]).read_text())
    assert manifest["resource_errors"][0]["gene"] == "HEAVY"


def test_routing_warnings_are_emitted_only_after_completion(panel, monkeypatch, capsys):
    inputs = panel({"GENE": (30, "T6C\n")})
    args = _parse(monkeypatch, inputs, "--msa", inputs["source"], "--plmc-binary", "/unused/plmc")
    capsys.readouterr()
    _capture_nextflow(monkeypatch, capsys)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 0
    after_completion = capsys.readouterr()
    assert "healthy tasks completed" in after_completion.out
    assert "this will not produce biologically accurate results" in after_completion.err


def test_all_blocked_collects_every_error_without_tools_or_output(panel, monkeypatch, capsys):
    inputs = panel({"FIRST": (600, "T6C\n"), "SECOND": (600, "T6C\n")})
    args = _parse(monkeypatch, inputs)
    capsys.readouterr()
    monkeypatch.setattr(controller.subprocess, "Popen", _deny)
    monkeypatch.setattr(controller, "validate_backend_tools", _deny)
    monkeypatch.setattr(controller, "write_resource_snapshot", _deny)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    printed = capsys.readouterr()
    assert "FIRST/codon" in printed.err
    assert "SECOND/codon" in printed.err
    assert printed.err.count("unschedulable") == 2
    assert not inputs["output"].exists()


def test_plan_only_includes_errors_and_healthy_tasks_without_launch_or_writes(panel, monkeypatch, capsys):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n"), "LIGHT": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs, "--resource-plan-only")
    capsys.readouterr()
    monkeypatch.setattr(controller.subprocess, "Popen", _deny)
    monkeypatch.setattr(controller, "validate_backend_tools", _deny)
    monkeypatch.setattr(controller, "write_resource_snapshot", _deny)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    report = json.loads(capsys.readouterr().out)
    assert {(task["gene"], task["side"]) for task in report["tasks"]} == {
        ("HEAVY", "protein"), ("LIGHT", "protein"),
    }
    assert report["resource_errors"][0]["gene"] == "HEAVY"
    assert "warnings" in report
    assert not inputs["output"].exists()


def test_current_dynamic_failure_report_is_merged_at_completion(panel, monkeypatch, capsys):
    inputs = panel({"GENE": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    stale = _write(inputs["output"] / ".bff-resources" / "old.errors.json", json.dumps([{
        "gene": "STALE", "side": "codon", "backend": "evmutation", "message": "previous failure",
        "stage": "resource_planning",
    }]))
    error = {
        "gene": "GENE", "side": "protein", "backend": "evmutation", "message": "generated MSA too large",
        "stage": "resource_planning",
    }
    capsys.readouterr()
    launches = _capture_nextflow(monkeypatch, capsys, dynamic_errors=[error])

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    printed = capsys.readouterr()
    assert "generated MSA too large" in printed.err
    assert "STALE" not in printed.err
    report = Path(launches[0][launches[0].index("--resource_errors") + 1])
    assert report != stale
    assert json.loads(stale.read_text())[0]["gene"] == "STALE"


def test_report_path_is_unique_without_invalidating_resource_snapshot(panel, monkeypatch, capsys):
    inputs = panel({"GENE": (30, "G4A\n")})
    launches = _capture_nextflow(monkeypatch, capsys)
    for attempt in range(2):
        args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
        capsys.readouterr()
        with pytest.raises(SystemExit) as stopped:
            controller.run_controller(args)
        assert stopped.value.code == 0

    first, second = launches
    assert first[first.index("--resource_errors") + 1] != second[second.index("--resource_errors") + 1]
    assert first[first.index("--resource_config") + 1] == second[second.index("--resource_config") + 1]
    assert first[first.index("--manifest") + 1] == second[second.index("--manifest") + 1]


@pytest.mark.parametrize("contents", ["not json", "{}", '[{"gene":"GENE"}]', None])
def test_unreadable_dynamic_report_cannot_silently_succeed(panel, monkeypatch, capsys, contents):
    inputs = panel({"GENE": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    capsys.readouterr()

    class Process:
        def __init__(self, command, **options):
            self.report = Path(command[command.index("--resource_errors") + 1])

        def wait(self):
            if contents is None:
                self.report.unlink()
            else:
                self.report.write_text(contents)
            print("healthy tasks completed")
            return 0

    monkeypatch.setattr(controller.subprocess, "Popen", Process)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    printed = capsys.readouterr()
    assert "healthy tasks completed" in printed.out
    assert "Cannot read resource-error report" in printed.err


def test_nextflow_launch_failure_keeps_collected_resource_errors(panel, monkeypatch, capsys):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    capsys.readouterr()

    def unavailable(*arguments, **keywords):
        raise FileNotFoundError("nextflow unavailable")

    monkeypatch.setattr(controller.subprocess, "Popen", unavailable)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    printed = capsys.readouterr()
    assert "HEAVY/codon" in printed.err
    assert "unschedulable" in printed.err
    assert "nextflow unavailable" in printed.err


@pytest.mark.parametrize("stage", ["snapshot", "run_report"])
def test_resource_setup_write_failure_keeps_collected_errors(panel, monkeypatch, capsys, stage):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    capsys.readouterr()
    monkeypatch.setattr(controller.subprocess, "Popen", _deny)
    message = f"{stage} directory is not writable"

    if stage == "snapshot":
        def failed_snapshot(*arguments, **keywords):
            raise OSError(message)

        monkeypatch.setattr(controller, "write_resource_snapshot", failed_snapshot)
    else:
        original = controller.tempfile.mkstemp

        def failed_report(*arguments, **keywords):
            if keywords.get("prefix") == "run-":
                raise OSError(message)
            return original(*arguments, **keywords)

        monkeypatch.setattr(controller.tempfile, "mkstemp", failed_report)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    printed = capsys.readouterr()
    assert "HEAVY/codon" in printed.err
    assert "unschedulable" in printed.err
    assert message in printed.err


def test_missing_backend_dependency_preserves_collected_resource_errors(panel, monkeypatch, capsys):
    inputs = panel({"HEAVY": (600, "G4A\nT6C\n")})
    args = _parse(monkeypatch, inputs)
    capsys.readouterr()
    monkeypatch.setattr(controller.subprocess, "Popen", _deny)

    with pytest.raises(SystemExit, match="--plmc-binary"):
        controller.run_controller(args)

    printed = capsys.readouterr()
    assert "HEAVY/codon" in printed.err
    assert "unschedulable" in printed.err
    assert not inputs["output"].exists()


def test_nextflow_zero_exit_without_final_report_is_not_success(panel, monkeypatch, capsys):
    inputs = panel({"GENE": (30, "G4A\n")})
    args = _parse(monkeypatch, inputs, "--plmc-binary", "/unused/plmc")
    capsys.readouterr()
    reports = []

    class Process:
        def __init__(self, command, **options):
            report = Path(command[command.index("--resource_errors") + 1])
            assert json.loads(report.read_text()) is None
            reports.append(report)

        def wait(self):
            print("Nextflow exited without publishing its final report")
            return 0

    monkeypatch.setattr(controller.subprocess, "Popen", Process)

    with pytest.raises(SystemExit) as stopped:
        controller.run_controller(args)

    assert stopped.value.code == 1
    assert len(reports) == 1
    assert json.loads(reports[0].read_text()) is None
    printed = capsys.readouterr()
    assert "Cannot read resource-error report" in printed.err
