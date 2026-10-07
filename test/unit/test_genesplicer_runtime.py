"""Exercise GeneSplicer command launching with a real executable fixture."""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from biofeaturefactory.genesplicer import genesplicer_ensemble as pipeline


def make_runtime(tmp_path, monkeypatch, spaced=False):
    binary_dir = tmp_path / ("installation with spaces" if spaced else "installation")
    model_dir = tmp_path / ("human models" if spaced else "models")
    temporary_dir = tmp_path / ("sequence files" if spaced else "sequences")
    for directory in (binary_dir, model_dir, temporary_dir):
        directory.mkdir()
    (model_dir / "config_file").write_text("fixture")
    executable = binary_dir / "genesplicer"
    executable.write_text(f"#!{sys.executable}\n" + '''import json
import os
import sys
from pathlib import Path

assert len(sys.argv) == 3, sys.argv
fasta = Path(sys.argv[1])
model = Path(sys.argv[2])
assert (model / "config_file").is_file()
(model / "invocation.json").write_text(json.dumps({
    "arguments": sys.argv[1:],
    "cwd": os.getcwd(),
    "fasta_path": str(fasta.resolve()),
    "fasta": fasta.read_text(),
}))
if (model / "fail").exists():
    print("fixture prediction failed", file=sys.stderr)
    raise SystemExit(7)
if not (model / "empty").exists():
    print("4 5 8.25 High donor")
''')
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(binary_dir))
    monkeypatch.setattr(pipeline.tempfile, "tempdir", str(temporary_dir))
    return executable, model_dir, temporary_dir


@pytest.mark.parametrize("spaced", (False, True))
@pytest.mark.parametrize("use_install_dir", (False, True))
def test_paths_remain_single_arguments(tmp_path, monkeypatch, spaced, use_install_dir):
    executable, model_dir, temporary_dir = make_runtime(tmp_path, monkeypatch, spaced)
    original_cwd = Path.cwd()
    install_dir = str(executable.parent) if use_install_dir else None

    result = pipeline._run_genesplicer_on_seq("TEST_WT", "ACGTAGT", install_dir, str(model_dir))

    assert result.to_dict("records") == [{
        "End5": 4, "End3": 5, "Score": 8.25,
        "confidence": "High", "splice_site_type": "donor",
    }]
    invocation = json.loads((model_dir / "invocation.json").read_text())
    assert invocation["arguments"][1] == str(model_dir)
    assert invocation["fasta"] == ">TEST_WT\nACGTAGT"
    assert Path(invocation["cwd"]) == (executable.parent if use_install_dir else original_cwd)
    fasta_path = Path(invocation["fasta_path"])
    assert fasta_path.parent == (executable.parent if use_install_dir else temporary_dir)
    assert invocation["arguments"][0] == (fasta_path.name if use_install_dir else str(fasta_path))
    assert not fasta_path.exists()
    assert Path.cwd() == original_cwd


@pytest.mark.parametrize("mapping_has_pkey", (False, True))
@pytest.mark.parametrize(
    "ref_nt,alt_nt,variant_class",
    [("A", "G", "snv"), ("ACGT" * 750, "TGCA" * 750, "mnv")],
    ids=("snv", "multikilobase_replacement"),
)
def test_process_gene_uses_pkey_for_mutant_fasta_header(
    tmp_path, monkeypatch, mapping_has_pkey, ref_nt, alt_nt, variant_class,
):
    executable, model_dir, temporary_dir = make_runtime(tmp_path, monkeypatch)
    gene_name = "GENE1"
    wt_sequence = f"CCGA{ref_nt}TTAG"
    mutant_sequence = f"CCGA{alt_nt}TTAG"
    mutant_token = f"{ref_nt}5{alt_nt}"
    fasta_path = tmp_path / f"{gene_name}.fasta"
    fasta_path.write_text(f">genomic\n{wt_sequence}\n")
    mapping_row = {"mutant": mutant_token, "genomic": mutant_token}
    minted_pkey = pipeline.mint_pkey(gene_name, mutant_token)
    expected_pkey = minted_pkey
    if mapping_has_pkey:
        expected_pkey = "GENE1-0123456789ab"
        assert expected_pkey != minted_pkey
        mapping_row["pkey"] = expected_pkey

    events, sites, variants, stats = pipeline._process_gene(
        fasta_path, pd.DataFrame([mapping_row]), str(executable.parent), str(model_dir),
        window=151, report_radius=151, visibility_threshold=1.0, high_cutoff=5.0,
        shift_bp=3, distance_k=75, cluster_radius=3, max_cluster_span=0, failure_map={},
    )

    invocation = json.loads((model_dir / "invocation.json").read_text())
    assert invocation["fasta"] == f">{expected_pkey}\n{mutant_sequence}"
    assert variants == [{
        "pkey": expected_pkey, "variant_class": variant_class, "length_delta": 0,
    }]
    assert len(events) == len(sites) == 1
    assert set(events[0]["pkey"]) == {expected_pkey}
    assert set(sites[0]["pkey"]) == {expected_pkey}
    assert stats["total_rows"] == stats["processed"] == 1
    assert stats["gene_skipped"] is None
    assert stats["skipped_detail"] == []
    assert not Path(invocation["fasta_path"]).exists()
    assert not list(executable.parent.glob("*.fasta"))
    assert not list(temporary_dir.glob("*.fasta"))


@pytest.mark.parametrize("use_install_dir", (False, True))
def test_nonzero_exit_remains_failure_and_cleans_fasta(tmp_path, monkeypatch, use_install_dir):
    executable, model_dir, temporary_dir = make_runtime(tmp_path, monkeypatch)
    (model_dir / "fail").touch()
    install_dir = str(executable.parent) if use_install_dir else None

    with pytest.raises(RuntimeError, match=r"GeneSplicer failed for TEST_WT \(rc=7\): fixture prediction failed"):
        pipeline._run_genesplicer_on_seq("TEST_WT", "ACGTAGT", install_dir, str(model_dir))

    invocation = json.loads((model_dir / "invocation.json").read_text())
    assert not Path(invocation["fasta_path"]).exists()
    assert not list(temporary_dir.glob("*.fasta"))


@pytest.mark.parametrize("use_install_dir", (False, True))
@pytest.mark.parametrize("failure", ("missing", "not_executable"))
def test_launch_errors_remain_runtime_failures_and_clean_fasta(
    tmp_path, monkeypatch, use_install_dir, failure,
):
    executable, model_dir, temporary_dir = make_runtime(tmp_path, monkeypatch)
    if failure == "missing":
        executable.unlink()
    else:
        executable.chmod(0o644)
    install_dir = str(executable.parent) if use_install_dir else None

    with pytest.raises(RuntimeError, match="GeneSplicer failed for TEST_WT") as caught:
        pipeline._run_genesplicer_on_seq("TEST_WT", "ACGTAGT", install_dir, str(model_dir))

    assert isinstance(caught.value.__cause__, OSError)
    assert not list(executable.parent.glob("*.fasta"))
    assert not list(temporary_dir.glob("*.fasta"))


def test_successful_empty_prediction_remains_empty(tmp_path, monkeypatch):
    executable, model_dir, _ = make_runtime(tmp_path, monkeypatch)
    (model_dir / "empty").touch()

    result = pipeline._run_genesplicer_on_seq(
        "TEST_WT", "ACGTAGT", str(executable.parent), str(model_dir))

    assert result.empty
    assert list(result.columns) == ["End5", "End3", "Score", "confidence", "splice_site_type"]
    invocation = json.loads((model_dir / "invocation.json").read_text())
    assert not Path(invocation["fasta_path"]).exists()
