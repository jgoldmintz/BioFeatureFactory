"""EVmutation cache validity across scoring routes and changed artifacts."""

import json
import sys

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller
from biofeaturefactory.mutation_effects.bin import evmutation_cache


def _write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


@pytest.fixture
def inputs(tmp_path):
    source = tmp_path / "source"
    output = tmp_path / "output"
    files = {
        "fasta": _write(source / "GENE" / "fastas" / "GENE.fasta", ">ORF\nATGGCTTAA\n"),
        "mutations": _write(source / "GENE" / "mappings" / "mutations" / "GENE_mutations.csv", "mutant\nG4A\nT6C\n"),
        "msa": _write(source / "GENE" / "MSA" / "GENE.msa.a2m", ">ORF\nMA\n>homolog\nMT\n"),
        "codon_msa": _write(source / "GENE" / "CodonMSA" / "GENE.codon.msa.fasta", ">ORF\nATGGCTTAA\n>homolog\nATGGCCTAA\n"),
        "protein_params": _write(output / "model_params" / "GENE.model_params", "protein fixture bytes\n"),
        "codon_params": _write(output / "codon_model_params" / "GENE.codon_model_params", "codon fixture bytes\n"),
    }
    hardware = _write(tmp_path / "hardware.json", json.dumps({"cpus": 8, "memory_gib": 32, "gpus": []}))
    return source, output, files, hardware


def _parse(monkeypatch, inputs, *extra):
    source, output, files, hardware = inputs
    monkeypatch.setattr(sys, "argv", [
        "mutEffects_controller", "--fasta", str(source), "--output", str(output),
        "--evmutation-only", "--resource-hardware", str(hardware), *map(str, extra),
    ])
    return controller.parse_args()


def _manifest(args):
    return controller.build_manifest(["GENE"], args)


def _publish(inputs, manifest, side):
    source, output, files, hardware = inputs
    directory = output / "GENE" / "EVmutation"
    tsv = _write(directory / f"GENE.{side}.tsv", "pkey\tscore\nGENE:G4A\t0.25\n")
    marker = directory / f"GENE.{side}.routing.json"
    msa = files["msa" if side == "protein" else "codon_msa"]
    evmutation_cache.write_completion(
        marker, manifest["ev_fingerprints"]["GENE"][side], tsv,
        files[f"{side}_params"], msa=msa,
    )
    return marker, tsv


def test_forced_codon_then_automatic_routing_invalidates_codon_cache(inputs, monkeypatch):
    source, output, files, hardware = inputs
    forced = _parse(monkeypatch, inputs, "-cm", files["codon_msa"])
    before = _manifest(forced)
    assert before["routing"]["GENE"]["score_missense_codon"]
    _publish(inputs, before, "codon")
    assert _manifest(forced)["codon_EVmutation"] == ["GENE"]
    automatic = _parse(monkeypatch, inputs)
    after = _manifest(automatic)
    assert not after["routing"]["GENE"]["score_missense_codon"]
    assert before["ev_fingerprints"]["GENE"]["codon"] != after["ev_fingerprints"]["GENE"]["codon"]
    assert after["codon_EVmutation"] == []
    assert controller.side_pending(automatic, after, "GENE", "evmutation", "codon")


def test_skip_codon_change_under_automatic_routing_invalidates_protein_cache(inputs, monkeypatch):
    skipped = _parse(monkeypatch, inputs, "--skip-codon", "evmutation")
    before = _manifest(skipped)
    assert before["routing"]["GENE"]["mode"] == "auto"
    _publish(inputs, before, "protein")
    assert _manifest(skipped)["EVmutation"] == ["GENE"]
    automatic = _parse(monkeypatch, inputs)
    after = _manifest(automatic)
    assert before["ev_fingerprints"]["GENE"]["protein"] != after["ev_fingerprints"]["GENE"]["protein"]
    assert after["EVmutation"] == []
    assert controller.side_pending(automatic, after, "GENE", "evmutation", "protein")


@pytest.mark.parametrize("mode", ["auto", "protein", "codon", "both"])
def test_unchanged_routes_reuse_verified_results(inputs, monkeypatch, mode):
    source, output, files, hardware = inputs
    options = []
    if mode in {"protein", "both"}:
        options += ["--msa", files["msa"]]
    if mode in {"codon", "both"}:
        options += ["-cm", files["codon_msa"]]
    args = _parse(monkeypatch, inputs, *options)
    before = _manifest(args)
    for side in before["ev_fingerprints"]["GENE"]:
        _publish(inputs, before, side)
    after = _manifest(args)
    for side in before["ev_fingerprints"]["GENE"]:
        artifact = "EVmutation" if side == "protein" else "codon_EVmutation"
        assert after[artifact] == ["GENE"]
        assert not controller.side_pending(args, after, "GENE", "evmutation", side)


@pytest.mark.parametrize("side,changed", [
    ("protein", "msa"), ("codon", "codon_msa"),
    ("protein", "protein_params"), ("codon", "codon_params"),
])
def test_changed_msa_or_params_invalidates_only_its_side(inputs, monkeypatch, side, changed):
    source, output, files, hardware = inputs
    args = _parse(monkeypatch, inputs)
    before = _manifest(args)
    for selected in ("protein", "codon"):
        _publish(inputs, before, selected)
    files[changed].write_text(files[changed].read_text() + "changed input bytes\n")
    after = _manifest(args)
    artifact = "EVmutation" if side == "protein" else "codon_EVmutation"
    other = "codon_EVmutation" if side == "protein" else "EVmutation"
    assert after[artifact] == []
    assert after[other] == ["GENE"]
    assert controller.side_pending(args, after, "GENE", "evmutation", side)


@pytest.mark.parametrize("problem", ["missing_marker", "changed_tsv"])
def test_unverified_or_changed_tsv_is_not_a_completed_task(inputs, monkeypatch, problem):
    args = _parse(monkeypatch, inputs)
    manifest = _manifest(args)
    marker, tsv = _publish(inputs, manifest, "protein")
    if problem == "missing_marker":
        marker.unlink()
    else:
        tsv.write_text("partial output")
    assert _manifest(args)["EVmutation"] == []


@pytest.mark.parametrize("codon_only", [False, True])
def test_prebuilt_models_still_need_no_plmc_binary(inputs, monkeypatch, codon_only):
    source, output, files, hardware = inputs
    options = []
    if codon_only:
        files["protein_params"].unlink()
        options += ["-cm", files["codon_msa"]]
    args = _parse(monkeypatch, inputs, *options)
    assert args.plmc_binary is None
    manifest = _manifest(args)
    assert manifest["codon_model_params"] == ["GENE"]
    assert manifest["model_params"] == ([] if codon_only else ["GENE"])
    controller.validate_backend_tools(args, manifest, ["GENE"])
    controller.validate_db_coverage(["GENE"], manifest, args)
