"""Mutation-directed model selection using real FASTA and mutation CSV inputs."""

import pytest

from biofeaturefactory.mutation_effects.bin.mutation_routing import classify_gene, choose_route


ORF = "ATGAAAGAATGGCTGACCTGTGATTAA"


@pytest.fixture
def inputs(tmp_path):
    fasta = tmp_path / "GENE.fasta"
    mutations = tmp_path / "GENE_mutations.csv"
    fasta.write_text(f">ORF\n{ORF}\n")
    mutations.write_text("mutant\n")
    return fasta, mutations


@pytest.mark.parametrize("token,expected", [
    ("A4G", "MISSENSE"),
    ("A6G", "SYNONYMOUS"),
    ("G11A", "STOP_GAIN"),
    ("T25C", "STOP_LOSS"),
    ("AAA4AAG", "SYNONYMOUS"),
    ("AAAGAA4GGATTT", "MNV"),
    ("TGG10TAG", "STOP_GAIN"),
    ("TAA25CAA", "STOP_LOSS"),
    ("TAA25TAG", "SYNONYMOUS"),
    ("AAAG4A", "INFRAME_DEL"),
    ("A4AGAA", "INFRAME_INS"),
    ("AAA4A", "FRAMESHIFT"),
    ("a6g", "SYNONYMOUS"),
    ("U2C", "MISSENSE"),
    ("invalid", "UNKNOWN"),
    ("A100G", "UNKNOWN"),
    ("G6A", "UNKNOWN"),
    ("A0G", "UNKNOWN"),
    ("gd.T5000C", "UNKNOWN"),
    ("ch.T70050267A", "UNKNOWN"),
])
def test_shared_core_classification(inputs, token, expected):
    fasta, mutations = inputs
    mutations.write_text(f"mutant\n{token}\n")
    assert classify_gene(fasta, mutations, "GENE") == [expected]


def test_real_mapping_csv_preserves_classes_and_removes_markers(inputs):
    fasta, mutations = inputs
    mutations.write_text("mutant\nA6G*\n A4G \nG11A\nA6G\n\n")
    assert classify_gene(fasta, mutations, "GENE") == ["MISSENSE", "STOP_GAIN", "SYNONYMOUS"]


def test_orf_header_is_preferred_over_other_fasta_records(inputs):
    fasta, mutations = inputs
    fasta.write_text(f">GENE genomic\nNNNNNN\n>orf transcript\n{ORF.lower()}\n")
    mutations.write_text("mutant\nA6G\n")
    assert classify_gene(fasta, mutations, "GENE") == ["SYNONYMOUS"]


def test_first_fasta_record_is_fallback(inputs):
    fasta, mutations = inputs
    fasta.write_text(f">GENE\n{ORF}\n>alternate\nNNNN\n")
    mutations.write_text("mutant\nA4G\n")
    assert classify_gene(fasta, mutations, "GENE") == ["MISSENSE"]


@pytest.mark.parametrize("sequence,token", [
    ("ANN", "A1G"), ("ATGA", "A4G"), ("", "A1G"), ("ATGAAA", "AAAG4A"),
])
def test_ambiguous_partial_empty_and_overrunning_reference_is_unknown(inputs, sequence, token):
    fasta, mutations = inputs
    fasta.write_text(f">ORF\n{sequence}\n")
    mutations.write_text(f"mutant\n{token}\n")
    assert classify_gene(fasta, mutations, "GENE") == ["UNKNOWN"]


def test_validation_filter_matches_shared_case_handling(inputs, tmp_path):
    fasta, mutations = inputs
    mutations.write_text("mutant\na4g\nA6G\n")
    validation_log = tmp_path / "validation.log"
    validation_log.write_text("GENE: mutation A4G expects A at position 4\n")
    assert classify_gene(fasta, mutations, "gene", validation_log) == ["SYNONYMOUS"]


def test_empty_or_fully_filtered_mutations_have_no_classes(inputs, tmp_path):
    fasta, mutations = inputs
    assert classify_gene(fasta, mutations, "GENE") == []
    mutations.write_text("mutant\nA4G\n")
    validation_log = tmp_path / "all_failed.log"
    validation_log.write_text("GENE: mutation A4G expects A at position 4\n")
    assert classify_gene(fasta, mutations, "GENE", validation_log) == []


@pytest.mark.parametrize("classes,protein,codon", [
    ([], False, False),
    (["MISSENSE"], True, False),
    (["SYNONYMOUS"], False, True),
    (["STOP_GAIN", "STOP_LOSS"], False, True),
    (["MNV", "INFRAME_DEL", "INFRAME_INS", "FRAMESHIFT"], True, False),
    (["MISSENSE", "SYNONYMOUS"], True, True),
    (["UNKNOWN"], True, True),
    (["synonymous", "SYNONYMOUS"], False, True),
])
def test_automatic_routing(classes, protein, codon):
    route = choose_route(classes)
    assert route == {
        "protein": protein, "codon": codon,
        "classes": sorted({label.upper() for label in classes}),
        "mode": "auto", "score_missense_codon": False, "warnings": [],
    }


@pytest.mark.parametrize("classes", [["SYNONYMOUS"], ["STOP_GAIN"], ["STOP_LOSS"]])
def test_forced_protein_warns_about_absent_codon_effects(classes, capsys):
    route = choose_route(classes, protein_explicit=True)
    assert route["protein"] is True
    assert route["codon"] is False
    assert route["mode"] == "protein"
    assert route["score_missense_codon"] is False
    assert "this will not produce biologically accurate results" in route["warnings"][0]
    assert "codon effects" in route["warnings"][0]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_forced_codon_scores_missense_and_explains_cost():
    route = choose_route(["MISSENSE"], codon_explicit=True)
    assert route["protein"] is False
    assert route["codon"] is True
    assert route["mode"] == "codon"
    assert route["score_missense_codon"] is True
    warning = route["warnings"][0]
    assert "higher memory" in warning
    assert "per-codon" in warning and "per-amino-acid" in warning
    assert "exponential" not in warning


def test_both_explicit_inputs_force_both_routes():
    route = choose_route(["MISSENSE"], protein_explicit=True, codon_explicit=True)
    assert route["protein"] is True
    assert route["codon"] is True
    assert route["mode"] == "both"
    assert route["score_missense_codon"] is False
    assert route["warnings"] == []


def test_both_explicit_inputs_preserve_normal_mixed_variant_distribution():
    route = choose_route(["MISSENSE", "SYNONYMOUS", "STOP_GAIN"], protein_explicit=True, codon_explicit=True)
    assert route["protein"] is True
    assert route["codon"] is True
    assert route["mode"] == "both"
    assert route["score_missense_codon"] is False
    assert route["warnings"] == []


def test_both_explicit_inputs_include_protein_for_synonymous_variants():
    route = choose_route(["SYNONYMOUS"], protein_explicit=True, codon_explicit=True)
    assert route["protein"] is True
    assert route["codon"] is True
    assert route["mode"] == "both"
    assert route["score_missense_codon"] is False
    assert route["warnings"] == []


@pytest.mark.parametrize("classes,protein_explicit,codon_explicit", [
    (["MISSENSE"], True, False), (["SYNONYMOUS"], False, True),
])
def test_matching_explicit_route_needs_no_warning(classes, protein_explicit, codon_explicit):
    route = choose_route(classes, protein_explicit=protein_explicit, codon_explicit=codon_explicit)
    assert route["protein"] is protein_explicit
    assert route["codon"] is codon_explicit
    assert route["warnings"] == []


@pytest.mark.parametrize("protein_explicit,codon_explicit", [(True, False), (False, True), (True, True)])
def test_explicit_alignment_does_not_train_without_remaining_mutations(protein_explicit, codon_explicit):
    route = choose_route([], protein_explicit=protein_explicit, codon_explicit=codon_explicit)
    assert route["protein"] is False
    assert route["codon"] is False
    assert route["score_missense_codon"] is False
    assert route["warnings"] == []
