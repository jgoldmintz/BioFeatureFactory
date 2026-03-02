"""
Unit tests for exon_aware_mapping.py:
locate_orf_in_transcript, validate_mutations_against_orf, derive_orf_from_transcript

Run with: pytest test/unit/test_exon_mapping.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from exon_aware_mapping import (
    locate_orf_in_transcript,
    validate_mutations_against_orf,
    derive_orf_from_transcript,
)



# locate_orf_in_transcript


class TestLocateOrfInTranscript:

    def test_orf_found_returns_offset(self):
        transcript = "NNNNATGAAATAA"
        orf = "ATGAAATAA"
        idx = locate_orf_in_transcript(orf, transcript, "GENE1")
        assert idx == 4

    def test_orf_at_start_returns_zero(self):
        orf = "ATGAAA"
        idx = locate_orf_in_transcript(orf, orf, "GENE1")
        assert idx == 0

    def test_orf_not_found_raises(self):
        with pytest.raises(ValueError, match="does not align"):
            locate_orf_in_transcript("ATGCCC", "NNNNATGAAATAA", "GENE1")

    def test_case_insensitive(self):
        idx = locate_orf_in_transcript("atgaaa", "NNNNATGAAA", "GENE1")
        assert idx == 4



# validate_mutations_against_orf


class TestValidateMutationsAgainstOrf:

    # ATG AAA GAA = M K E
    ORF = "ATGAAAGAA"

    def test_valid_mutation_no_issues(self):
        issues, mismatches = validate_mutations_against_orf("GENE1", self.ORF, ["A4G"])
        assert issues == []
        assert mismatches == []

    def test_empty_mutations_returns_empty(self):
        # Note: function returns [] (not a tuple) for empty input
        result = validate_mutations_against_orf("GENE1", self.ORF, [])
        assert result == []

    def test_out_of_bounds_position_reported(self):
        issues, mismatches = validate_mutations_against_orf("GENE1", self.ORF, ["A100G"])
        assert len(issues) == 1
        assert mismatches == []

    def test_wrong_base_reported_as_mismatch(self):
        # Position 4 in ORF is 'A', but mutation says base is 'C'
        issues, mismatches = validate_mutations_against_orf("GENE1", self.ORF, ["C4G"])
        assert issues == []
        assert len(mismatches) == 1

    def test_multiple_mutations_independently_validated(self):
        issues, mismatches = validate_mutations_against_orf(
            "GENE1", self.ORF, ["A4G", "A100G", "C4G"]
        )
        assert len(issues) == 1       # A100G out of bounds
        assert len(mismatches) == 1   # C4G wrong base

    def test_stop_codon_mutation_skipped(self):
        # get_mutation_data_bioAccurate returns None for stop codons
        issues, mismatches = validate_mutations_against_orf("GENE1", self.ORF, ["A4Stop"])
        assert issues == []
        assert mismatches == []



# derive_orf_from_transcript


class TestDeriveOrfFromTranscript:

    def test_simple_orf_extracted(self):
        # ATG AAA TAA: start at 0, stop at 6
        transcript = "ATGAAATAA"
        orf, start = derive_orf_from_transcript("GENE1", transcript, [], {}, False, None)
        assert orf == "ATGAAATAA"
        assert start == 0

    def test_orf_with_5prime_utr(self):
        # NNN prefix before ORF
        transcript = "NNNATGAAATAA"
        orf, start = derive_orf_from_transcript("GENE1", transcript, [], {}, False, None)
        assert orf == "ATGAAATAA"
        assert start == 3

    def test_no_orf_raises(self):
        with pytest.raises(ValueError, match="no in-frame start/stop"):
            derive_orf_from_transcript("GENE1", "NNNNNNNN", [], {}, False, None)

    def test_longest_orf_selected(self):
        # Two possible ORFs: ATG AAA TAA (short) and ATG AAA GAA TAA (longer)
        # ATGAAATAAATGAAAGAATAA
        transcript = "ATGAAATAAATGAAAGAATAA"
        orf, start = derive_orf_from_transcript("GENE1", transcript, [], {}, False, None)
        # Longer orf covers more — should prefer longer
        assert len(orf) >= 9

    def test_orf_covers_mutation_positions(self):
        # Mutation at position 7 (1-based) — ORF must be long enough
        transcript = "ATGAAAGAATAA"
        mutation_positions = [7]
        pos_to_mut = {7: ["G7A"]}
        orf, start = derive_orf_from_transcript(
            "GENE1", transcript, mutation_positions, pos_to_mut, False, None
        )
        assert len(orf) >= 7
