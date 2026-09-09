"""EVmutation routing provenance and atomic completion marker regressions."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from biofeaturefactory.mutation_effects.bin import evmutation_cache as cache


@pytest.fixture
def artifacts(tmp_path):
    files = {
        "fasta": tmp_path / "GENE.fasta",
        "mutations": tmp_path / "GENE_mutations.csv",
        "validation_log": tmp_path / "validation.log",
        "tsv": tmp_path / "GENE.protein.tsv",
        "params": tmp_path / "external-prebuilt.model_params",
        "msa": tmp_path / "GENE.msa.a2m",
        "marker": tmp_path / "GENE.protein.routing.json",
    }
    files["fasta"].write_text(">ORF\nATGGCT\n")
    files["mutations"].write_text("mutant\nG4A\nT6C\n")
    files["validation_log"].write_text("GENE: validation complete\n")
    files["tsv"].write_text("nt_mutant\tscore\nG4A\t0.5\n")
    files["params"].write_bytes(b"fixture binary params\x00\x01")
    files["msa"].write_text(">ORF\nMA\n>homolog\nM-\n")
    return files


def fingerprint(artifacts, **options):
    return cache.routing_fingerprint(
        "GENE", "protein", artifacts["fasta"], artifacts["mutations"],
        artifacts["validation_log"], **options,
    )


def write_marker(artifacts, expected=None):
    expected = expected or fingerprint(artifacts)
    cache.write_completion(artifacts["marker"], expected, artifacts["tsv"],
                           artifacts["params"], artifacts["msa"])
    return expected


@pytest.mark.parametrize("setting", ["skip_codon", "score_missense_codon"])
def test_routing_mode_changes_invalidate_fingerprint(artifacts, setting):
    assert fingerprint(artifacts) != fingerprint(artifacts, **{setting: True})


@pytest.mark.parametrize("input_name", ["fasta", "mutations", "validation_log"])
def test_input_content_changes_invalidate_fingerprint(artifacts, input_name):
    before = fingerprint(artifacts)
    path = artifacts[input_name]
    path.write_text(path.read_text() + "changed\n")
    assert fingerprint(artifacts) != before


@pytest.mark.parametrize("filename", ["evmutation_pipeline.py", "codon_encoding.py"])
def test_backend_implementation_changes_invalidate_fingerprint(artifacts, monkeypatch, filename):
    before = fingerprint(artifacts)
    original_digest = cache.file_digest

    def changed_digest(path):
        return "changed-implementation" if Path(path).name == filename else original_digest(path)

    monkeypatch.setattr(cache, "file_digest", changed_digest)
    assert fingerprint(artifacts) != before


def test_fingerprint_ignores_input_locations_and_tracks_gene_side(artifacts, tmp_path):
    before = fingerprint(artifacts)
    relocated = dict(artifacts)
    for name in ("fasta", "mutations", "validation_log"):
        relocated[name] = tmp_path / f"relocated_{artifacts[name].name}"
        relocated[name].write_bytes(artifacts[name].read_bytes())
    assert fingerprint(relocated) == before
    assert cache.routing_fingerprint(
        "OTHER", "protein", artifacts["fasta"], artifacts["mutations"], artifacts["validation_log"],
    ) != before
    assert cache.routing_fingerprint(
        "GENE", "codon", artifacts["fasta"], artifacts["mutations"], artifacts["validation_log"],
    ) != before


def test_verified_completion_reuses_prebuilt_params_without_copying(artifacts):
    params_contents = artifacts["params"].read_bytes()
    expected = write_marker(artifacts)
    assert cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"],
                                   artifacts["params"], artifacts["msa"])
    assert artifacts["params"].read_bytes() == params_contents
    assert not (artifacts["params"].parent / "GENE.model_params").exists()


@pytest.mark.parametrize("artifact", ["tsv", "params", "msa"])
@pytest.mark.parametrize("change", ["same_size", "empty", "missing"])
def test_changed_missing_or_partial_artifacts_are_not_reused(artifacts, artifact, change):
    expected = write_marker(artifacts)
    path = artifacts[artifact]
    if change == "same_size":
        contents = path.read_bytes()
        path.write_bytes(bytes([contents[0] ^ 1]) + contents[1:])
    elif change == "empty":
        path.write_bytes(b"")
    else:
        path.unlink()
    assert not cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"],
                                       artifacts["params"], artifacts["msa"])


@pytest.mark.parametrize("contents", [None, "{unfinished", "[]", "{}", '{"success": true}'])
def test_missing_legacy_or_malformed_markers_are_not_reused(artifacts, contents):
    if contents is not None:
        artifacts["marker"].write_text(contents)
    assert not cache.verify_completion(artifacts["marker"], fingerprint(artifacts), artifacts["tsv"],
                                       artifacts["params"], artifacts["msa"])


def test_changed_fingerprint_and_unverified_msa_are_rejected(artifacts):
    expected = write_marker(artifacts)
    assert not cache.verify_completion(artifacts["marker"], "other-fingerprint", artifacts["tsv"],
                                       artifacts["params"], artifacts["msa"])
    assert not cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"], artifacts["params"])
    cache.write_completion(artifacts["marker"], expected, artifacts["tsv"], artifacts["params"])
    assert cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"], artifacts["params"])
    assert not cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"],
                                       artifacts["params"], artifacts["msa"])


@pytest.mark.parametrize("artifact", ["tsv", "params", "msa"])
def test_marker_is_not_replaced_when_an_artifact_is_incomplete(artifacts, artifact):
    expected = write_marker(artifacts)
    previous = artifacts["marker"].read_bytes()
    artifacts[artifact].write_bytes(b"")
    with pytest.raises(ValueError, match="Missing or empty"):
        write_marker(artifacts, expected)
    assert artifacts["marker"].read_bytes() == previous


def test_atomic_replace_failure_preserves_previous_marker_and_cleans_temporary(artifacts, monkeypatch):
    write_marker(artifacts)
    previous = artifacts["marker"].read_bytes()
    previous_names = set(artifacts["marker"].parent.iterdir())

    def failed_replace(source, destination):
        assert json.loads(Path(source).read_text())["success"] is True
        assert Path(destination).read_bytes() == previous
        raise OSError("simulated rename failure")

    monkeypatch.setattr(cache.os, "replace", failed_replace)
    with pytest.raises(OSError, match="simulated rename failure"):
        write_marker(artifacts)
    assert artifacts["marker"].read_bytes() == previous
    assert set(artifacts["marker"].parent.iterdir()) == previous_names


def test_cli_writes_complete_marker_with_selected_msa(artifacts):
    expected = fingerprint(artifacts)
    result = subprocess.run([
        sys.executable, cache.__file__, "write",
        "--marker", str(artifacts["marker"]), "--fingerprint", expected,
        "--tsv", str(artifacts["tsv"]), "--params", str(artifacts["params"]),
        "--msa", str(artifacts["msa"]),
    ], capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    assert result.returncode == 0, result.stderr
    assert cache.verify_completion(artifacts["marker"], expected, artifacts["tsv"],
                                   artifacts["params"], artifacts["msa"])
