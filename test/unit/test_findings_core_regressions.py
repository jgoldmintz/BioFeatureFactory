import json
import os
from pathlib import Path
import sys

import pysam
import pytest

from biofeaturefactory.core import codon_msa_pipeline as codon
from biofeaturefactory.core import vcf_converter as vcf
from biofeaturefactory.lib.utility import load_validation_failures, trim_muts


@pytest.fixture
def vcf_inputs(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">1\nAACCGGTTAACCGGTT\n")
    pysam.faidx(str(reference))
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text("GENE\t1\t+\t1\t16\t1\t16\n")
    mutations = tmp_path / "GENE_mutations.csv"
    mutations.write_text("mutant\ngd.A1G\n")
    return reference, annotation, mutations, tmp_path / "results"


@pytest.mark.parametrize("change", ["touch", "content_same_stat"])
def test_annotation_changes_invalidate_vcf_cli_cache(vcf_inputs, monkeypatch, capsys, change):
    reference, annotation, mutations, output = vcf_inputs
    monkeypatch.setattr(sys, "argv", [
        "vcf_converter", "-m", str(mutations), "-r", str(reference),
        "-a", str(annotation), "-o", str(output), "--chromosome-format", "simple",
    ])
    vcf.main()
    capsys.readouterr()
    vcf.main()
    assert "cached 1 variants" in capsys.readouterr().out
    previous = annotation.stat()
    if change == "touch":
        os.utime(annotation, ns=(previous.st_atime_ns, previous.st_mtime_ns + 1000000))
    else:
        annotation.write_text("GENE\t1\t+\t2\t16\t2\t16\n")
        os.utime(annotation, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    vcf.main()
    assert "cached" not in capsys.readouterr().out
    rows = [line.split("\t") for line in (output / "GENE/vcf/GENE.vcf").read_text().splitlines()
            if not line.startswith("#")]
    assert rows[0][1] == ("1" if change == "touch" else "2")


def test_legacy_vcf_cache_is_not_reused(vcf_inputs, monkeypatch, capsys):
    reference, annotation, mutations, output = vcf_inputs
    monkeypatch.setattr(sys, "argv", ["vcf_converter", "-m", str(mutations),
        "-r", str(reference), "-a", str(annotation), "-o", str(output)])
    vcf.main()
    cache_path = output / ".vcf_converter_cache.json"
    cache = json.loads(cache_path.read_text())
    for entry in cache.values():
        entry.pop("annotation_state", None)
    cache_path.write_text(json.dumps(cache))
    capsys.readouterr()
    vcf.main()
    assert "cached" not in capsys.readouterr().out


@pytest.mark.parametrize("mapped,accepted", [
    ("AAC1A", True), ("AAT1A", False), ("TTT15T", False), ("A99G", False),
    ("A1G", True), ("A1AC", True), ("AAC1GGT", True),
])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_mapping_validation_checks_entire_ref(vcf_inputs, capsys, mapped, accepted, strand):
    reference, annotation, mutations, output = vcf_inputs
    mutations.write_text("mutant\nAAC1A\n")
    annotation.write_text(f"GENE\t1\t{strand}\t1\t16\t1\t16\n")
    mapping = reference.parent / "mapping.csv"
    mapping.write_text(f"mutant,chromosome\nAAC1A,{mapped}\n")
    success, count, error, path = vcf.process_single_file(
        str(mutations), str(output), reference_fasta=str(reference),
        annotation_file=str(annotation), chromosome_mapping_input=str(mapping),
        validate_mapping=True,
    )
    assert success, error
    assert count == int(accepted)
    rows = [line for line in Path(path).read_text().splitlines() if not line.startswith("#")]
    assert len(rows) == int(accepted)
    if not accepted:
        assert "REF_MISMATCH" in capsys.readouterr().err


@pytest.mark.parametrize("token", ["ACAA12A", "A12ACT", "AC12GT", "U12A", "a12act", "gd.AC12A", "ch.A12G"])
def test_multibase_validation_log_reaches_filter(tmp_path, token):
    log = tmp_path / "validation.log"
    log.write_text(f"GENE: mutation {token} expects reference bases\n")
    mutations = tmp_path / "GENE_mutations.csv"
    mutations.write_text(f"mutant\n{token}\nG8A\n")
    assert load_validation_failures(log) == {"GENE": {token}}
    assert trim_muts(mutations, log=log, gene_name="GENE") == ["G8A"]


@pytest.mark.parametrize("token", ["AC0A", "A1Gjunk", "A1G/garbage", "R3W", "A1N"])
def test_validation_log_does_not_accept_partial_or_invalid_tokens(tmp_path, token):
    log = tmp_path / "validation.log"
    log.write_text(f"GENE: mutation {token} expects reference bases\n")
    assert load_validation_failures(log) == {}


def test_refseq_source_is_independent_of_index_location(tmp_path, monkeypatch):
    database = tmp_path / "Bio_DBs"
    database.mkdir()
    source = database / "refseq_proteins_merged.faa"
    source.write_text(">NP_1.1\nMA\n")
    index = tmp_path / "index" / "refseq"
    index.parent.mkdir()
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        if command[1] == "convertalis":
            Path(command[5]).write_text("focus\tNP_1.1\t1\t2\t1\t1\n")

    monkeypatch.setattr(codon.subprocess, "run", run)
    hits = codon._search_focus_protein_against_refseq(
        "MA", "focus", db_root=str(database), target_db_base=str(index),
    )
    assert hits == ["NP_1.1"]
    assert calls[0] == ["mmseqs", "createdb", str(source), str(index)]
    assert codon._default_mmseqs_target_paths(str(index), None) == (
        str(index), str(index.parent / "refseq_proteins_merged.faa"),
    )
