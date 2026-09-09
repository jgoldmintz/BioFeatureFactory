"""CLI resource planning and provenance-aware mutation-effects resume behavior."""

import json
import sys
from pathlib import Path

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller
from biofeaturefactory.mutation_effects.bin.adabmdca_task import _file_record


def _write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


@pytest.fixture
def inputs(tmp_path):
    source = tmp_path / "source"
    files = {
        "fasta": _write(source / "NPM1" / "fastas" / "NPM1.fasta", ">ORF\nATGGCT\n"),
        "mutations": _write(source / "NPM1" / "mappings" / "mutations" / "NPM1_mutations.csv", "mutant\nG4A\nT6C\n"),
        "msa": _write(source / "NPM1" / "MSA" / "NPM1.msa.a2m", ">ORF\nMA\n>homolog\nM-\n"),
        "codon_msa": _write(source / "NPM1" / "CodonMSA" / "NPM1.codon.msa.fasta", ">ORF\nATGGCT\n>homolog\nATG---\n"),
    }
    hardware = _write(tmp_path / "hardware.json", json.dumps({
        "cpus": 8, "memory_gib": 64,
        "gpus": [{"index": "0", "uuid": "GPU-allocated", "memory_gib": 80}],
    }))
    return {
        "source": source, "output": tmp_path / "output", "files": files,
        "hardware": hardware, "lease_dir": tmp_path / "leases",
    }


def _parse(monkeypatch, inputs, *extra):
    monkeypatch.setattr(sys, "argv", [
        "mutEffects_controller", "--fasta", str(inputs["source"]),
        "--output", str(inputs["output"]), "--adabmdca-only",
        "--resource-hardware", str(inputs["hardware"]),
        "--gpu-lease-dir", str(inputs["lease_dir"]), *map(str, extra),
    ])
    return controller.parse_args()


def _prepare(args):
    genes = controller._resolve_genes(args.fasta)
    manifest = controller.build_manifest(genes, args)
    config, plans = controller.prepare_resource_plans(args, manifest, genes)
    return config, plans, manifest, genes


def _publish_completion(output, plan):
    gene, side = plan["gene"], plan["side"]
    artifact_directory = output / gene / "adabmDCA"
    tsv = _write(artifact_directory / f"{gene}.{side}.tsv", "pkey\tscore\nNPM1:G4A\t0.5\n")
    params = _write(output / f"adabmdca_{side}_params" / f"{gene}.{side}_adabm_params", "h 0 A 0.5\n")
    marker = _write(artifact_directory / f"{gene}.{side}.complete.json", json.dumps({
        "gene": gene, "side": side, "fingerprint": plan["fingerprint"],
        "device": plan["device"], "success": True,
        "files": {artifact.name: _file_record(artifact) for artifact in (tsv, params)},
    }))
    return marker, tsv, params


def _deny_launch(*args, **kwargs):
    pytest.fail("Resource planning unexpectedly launched an external process")


def test_verified_completion_skips_training_plan_and_dependency_checks(inputs, monkeypatch):
    args = _parse(monkeypatch, inputs)
    config, plans, manifest, genes = _prepare(args)
    assert {plan["side"] for plan in plans} == {"protein", "codon"}
    for plan in plans:
        _publish_completion(args.output, plan)
    monkeypatch.setattr(controller.resource_planner, "plan_task", lambda *args, **kwargs: pytest.fail("Completed artifacts were replanned for training"))
    config, completed, manifest, genes = _prepare(args)
    assert all(plan["verified_complete"] and plan["device"] == "cached" for plan in completed)
    for side in ("protein", "codon"):
        assert manifest[f"adabmdca_{side}"] == ["NPM1"]
        assert manifest[f"adabmdca_{side}_params"] == ["NPM1"]
    monkeypatch.setattr(controller.importlib.util, "find_spec", lambda name: pytest.fail(f"Unneeded dependency lookup: {name}"))
    controller.validate_backend_tools(args, manifest, genes)


@pytest.mark.parametrize("option,value", [
    ("--adabmdca-seed", "1"), ("--adabmdca-nepochs", "3"),
    ("--adabmdca-model", "bmDCA"), ("--adabmdca-dtype", "float64"),
])
def test_model_setting_changes_invalidate_completed_results(inputs, monkeypatch, option, value):
    args = _parse(monkeypatch, inputs)
    config, plans, manifest, genes = _prepare(args)
    for plan in plans:
        _publish_completion(args.output, plan)
    changed_args = _parse(monkeypatch, inputs, option, value)
    config, changed, manifest, genes = _prepare(changed_args)
    assert all(not plan["verified_complete"] and plan["params"] is None for plan in changed)
    original = {plan["side"]: plan["fingerprint"] for plan in plans}
    assert all(plan["fingerprint"] != original[plan["side"]] for plan in changed)
    assert manifest["adabmdca_protein"] == manifest["adabmdca_codon"] == []


@pytest.mark.parametrize("kind,contents,invalidated", [
    ("fasta", ">ORF\nATGGCC\n", {"protein", "codon"}),
    ("mutations", "mutant\nG4T\nT6C\n", {"protein", "codon"}),
    ("msa", ">ORF\nMA\n>homolog\nMT\n", {"protein"}),
    ("codon_msa", ">ORF\nATGGCT\n>homolog\nATGGCC\n", {"codon"}),
])
def test_input_changes_invalidate_only_affected_completed_sides(inputs, monkeypatch, kind, contents, invalidated):
    args = _parse(monkeypatch, inputs)
    config, plans, manifest, genes = _prepare(args)
    for plan in plans:
        _publish_completion(args.output, plan)
    inputs["files"][kind].write_text(contents)
    config, changed, manifest, genes = _prepare(args)
    assert {plan["side"] for plan in changed if not plan["verified_complete"]} == invalidated


@pytest.mark.parametrize("problem", ["missing_marker", "false_success", "wrong_fingerprint", "partial_params", "changed_tsv"])
def test_partial_or_unverified_old_artifacts_are_not_reused(inputs, monkeypatch, problem):
    args = _parse(monkeypatch, inputs, "--skip-codon")
    config, plans, manifest, genes = _prepare(args)
    marker, tsv, params = _publish_completion(args.output, plans[0])
    if problem == "missing_marker":
        marker.unlink()
    elif problem == "partial_params":
        params.write_text("partial checkpoint")
    elif problem == "changed_tsv":
        tsv.write_text("partial scores")
    else:
        contents = json.loads(marker.read_text())
        contents["success" if problem == "false_success" else "fingerprint"] = False if problem == "false_success" else "unrelated run"
        marker.write_text(json.dumps(contents))
    config, changed, manifest, genes = _prepare(args)
    assert not changed[0]["verified_complete"]
    assert changed[0]["params"] is None
    assert manifest["adabmdca_protein"] == manifest["adabmdca_protein_params"] == []


def test_completed_tasks_resume_when_current_ram_cannot_fit_training(inputs, monkeypatch):
    args = _parse(monkeypatch, inputs)
    config, plans, manifest, genes = _prepare(args)
    for plan in plans:
        _publish_completion(args.output, plan)
    inputs["hardware"].write_text(json.dumps({"cpus": 1, "memory_gib": 0.01, "gpus": []}))
    low_resource_args = _parse(monkeypatch, inputs, "--resume")
    config, cached, manifest, genes = _prepare(low_resource_args)
    assert config["hardware"]["memory_gib"] == 0.01
    assert all(plan["verified_complete"] and plan["device"] == "cached" for plan in cached)


@pytest.mark.parametrize("directory_mode", [False, True])
def test_explicit_prebuilt_params_are_inventoried_and_skip_training_dependencies(inputs, monkeypatch, directory_mode):
    options = []
    supplied = {}
    for side in ("protein", "codon"):
        params = _write(inputs["source"].parent / f"explicit-{side}" / f"NPM1.{side}_adabm_params", "h 0 A 0.25\n")
        supplied[side] = params
        options.extend([f"--adabmdca-{side}-params", params.parent if directory_mode else params])
    args = _parse(monkeypatch, inputs, *options)
    config, plans, manifest, genes = _prepare(args)
    for plan in plans:
        assert plan["device"] == "cpu"
        assert plan["gpu_memory_gib"] == 0
        assert Path(plan["params"]) == supplied[plan["side"]]
        assert manifest[f"adabmdca_{plan['side']}_params"] == ["NPM1"]
    monkeypatch.setattr(controller.importlib.util, "find_spec", lambda name: pytest.fail(f"Prebuilt params triggered training dependency lookup: {name}"))
    controller.validate_backend_tools(args, manifest, genes)


def test_resource_plan_only_writes_nothing_and_launches_nothing(inputs, monkeypatch, capsys):
    _write(inputs["output"] / ".evmutation_manifest.json", "existing inventory\n")
    args = _parse(monkeypatch, inputs, "--resource-plan-only")
    capsys.readouterr()
    root = inputs["source"].parent
    before = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
    monkeypatch.setattr(controller.subprocess, "Popen", _deny_launch)
    monkeypatch.setattr(controller.resource_planner, "detect_hardware", _deny_launch)
    monkeypatch.setattr(controller, "validate_backend_tools", _deny_launch)
    controller.run_controller(args)
    after = {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}
    assert after == before
    assert not inputs["lease_dir"].exists()
    printed = json.loads(capsys.readouterr().out)
    assert printed["hardware"]["cpus"] == 8
    assert {plan["side"] for plan in printed["tasks"]} == {"protein", "codon"}


@pytest.mark.parametrize("option,value", [
    ("--threads", "0"), ("--resource-cpus", "0"), ("--resource-memory-gib", "0"),
    ("--resource-memory-gib", "nan"), ("--resource-headroom", "0"),
    ("--resource-headroom", "1"), ("--resource-memory-margin", "0.9"),
    ("--resource-memory-margin", "nan"), ("--gpu-wait-timeout", "-1"),
    ("--msa-memory-gib", "0"), ("--evmutation-memory-gib", "inf"),
    ("--adabmdca-device", "cuda:invalid"),
])
def test_invalid_resource_cli_values_are_rejected_before_launch(inputs, monkeypatch, option, value):
    monkeypatch.setattr(controller.subprocess, "Popen", _deny_launch)
    with pytest.raises(SystemExit):
        args = _parse(monkeypatch, inputs, option, value, "--resource-plan-only")
        controller.run_controller(args)
    assert not inputs["output"].exists()


@pytest.mark.parametrize("override", [
    {"NPM1.protein": {"threads": 1.5}},
    {"NPM1.protein": {"gpu_memory_gib": -1}},
    {"NPM1.protein": {"unknown_resource": 1}},
])
def test_invalid_resource_overrides_report_controller_error(inputs, monkeypatch, override):
    overrides = _write(inputs["source"].parent / "overrides.json", json.dumps(override))
    args = _parse(monkeypatch, inputs, "--resource-overrides", overrides, "--resource-plan-only")
    monkeypatch.setattr(controller.subprocess, "Popen", _deny_launch)
    with pytest.raises(SystemExit, match="ERROR:"):
        controller.run_controller(args)
    assert not inputs["output"].exists()
