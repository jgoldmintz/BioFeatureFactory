"""
Unit tests for codon utility functions in utility.py:
compute_cai, compute_tai, extract_codon_with_bicodons, detect_alphabet

Run with: pytest test/unit/test_utils_codon.py -v
"""

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import compute_cai, compute_tai, extract_codon_with_bicodons, detect_alphabet



# compute_cai


class TestComputeCAI:

    # Controlled weights so tests don't depend on human reference table
    W = {"ATG": 1.0, "AAA": 0.5, "GAA": 0.8}

    def test_geometric_mean_of_weights(self):
        # ATG (1.0) + AAA (0.5) → exp((ln1.0 + ln0.5) / 2) = sqrt(0.5)
        result = compute_cai("ATGAAA", w_values=self.W)
        expected = math.exp((math.log(1.0) + math.log(0.5)) / 2)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_single_codon(self):
        result = compute_cai("ATG", w_values={"ATG": 0.8})
        assert result == pytest.approx(0.8, rel=1e-6)

    def test_stop_codon_skipped(self):
        # TAA is stop — should be skipped, only ATG counted
        result = compute_cai("ATGTAA", w_values={"ATG": 0.8})
        assert result == pytest.approx(0.8, rel=1e-6)

    def test_unknown_codon_skipped(self):
        # NNN has no weight entry — skipped, only ATG counted
        result = compute_cai("ATGNNN", w_values={"ATG": 0.6})
        assert result == pytest.approx(0.6, rel=1e-6)

    def test_empty_or_all_stops_returns_none(self):
        assert compute_cai("TAATAG", w_values=self.W) is None

    def test_lowercase_input_handled(self):
        result = compute_cai("atgaaa", w_values=self.W)
        assert result is not None

    def test_rna_u_converted(self):
        result_dna = compute_cai("ATGAAA", w_values=self.W)
        result_rna = compute_cai("AUGAAA", w_values=self.W)
        assert result_dna == pytest.approx(result_rna, rel=1e-6)

    def test_result_between_zero_and_one(self):
        result = compute_cai("ATGAAAGAA", w_values=self.W)
        assert 0.0 < result <= 1.0



# compute_tai


class TestComputeTAI:

    TAI = {"ATG": 0.9, "AAA": 0.3, "GAA": 0.6}

    def test_geometric_mean(self):
        result = compute_tai("ATGAAA", tai_weights=self.TAI)
        expected = math.exp((math.log(0.9) + math.log(0.3)) / 2)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_stop_codon_skipped(self):
        result = compute_tai("ATGTAA", tai_weights={"ATG": 0.9})
        assert result == pytest.approx(0.9, rel=1e-6)

    def test_zero_weight_skipped(self):
        result = compute_tai("ATGAAA", tai_weights={"ATG": 0.9, "AAA": 0.0})
        assert result == pytest.approx(0.9, rel=1e-6)

    def test_all_skipped_returns_none(self):
        assert compute_tai("TAATAG", tai_weights=self.TAI) is None

    def test_result_between_zero_and_one(self):
        result = compute_tai("ATGAAAGAA", tai_weights=self.TAI)
        assert 0.0 < result <= 1.0



# extract_codon_with_bicodons


class TestExtractCodonWithBicodons:

    # ATG AAA GAA TGG = M K E W  (12 nt, 4 codons)
    SEQ = "ATGAAAGAATGG"

    def test_middle_codon_has_both_bicodons(self):
        # codon 2 = AAA (positions 4-6), middle codon
        codon, fwd, rev, pos_in_codon, pos, codon_num = extract_codon_with_bicodons("A4G", self.SEQ)
        assert codon == "AAA"
        assert fwd != ""   # forward bicodon available
        assert rev != ""   # reverse bicodon available
        assert codon_num == 2

    def test_first_codon_has_forward_only(self):
        # codon 1 = ATG (positions 1-3)
        codon, fwd, rev, pos_in_codon, pos, codon_num = extract_codon_with_bicodons("A1T", self.SEQ)
        assert codon == "ATG"
        assert fwd != ""   # forward bicodon (ATG + AAA)
        assert rev == ""   # no reverse (no preceding codon)
        assert codon_num == 1

    def test_last_codon_has_reverse_only(self):
        # codon 4 = TGG (positions 10-12)
        codon, fwd, rev, pos_in_codon, pos, codon_num = extract_codon_with_bicodons("T10A", self.SEQ)
        assert codon == "TGG"
        assert fwd == ""   # no forward (last codon)
        assert rev != ""   # reverse bicodon (GAA + TGG)
        assert codon_num == 4

    def test_forward_bicodon_correct_sequence(self):
        codon, fwd, rev, _, _, _ = extract_codon_with_bicodons("A4G", self.SEQ)
        assert fwd == "AAAGAA"   # AAA + GAA

    def test_reverse_bicodon_correct_sequence(self):
        codon, fwd, rev, _, _, _ = extract_codon_with_bicodons("A4G", self.SEQ)
        assert rev == "ATGAAA"   # ATG + AAA

    def test_pos_in_codon_first_position(self):
        _, _, _, pos_in_codon, _, _ = extract_codon_with_bicodons("A4G", self.SEQ)
        assert pos_in_codon == 0   # 1st base of codon

    def test_pos_in_codon_third_position(self):
        _, _, _, pos_in_codon, _, _ = extract_codon_with_bicodons("A6G", self.SEQ)
        assert pos_in_codon == 2   # 3rd base of codon

    def test_stop_codon_mutation_returns_none(self):
        codon, fwd, rev, _, _, _ = extract_codon_with_bicodons("A4Stop", self.SEQ)
        assert codon is None



# detect_alphabet


class TestDetectAlphabet:

    def test_pure_dna_is_nucleotide(self):
        assert detect_alphabet("ACGTACGTACGT") == "nucleotide"

    def test_pure_rna_is_nucleotide(self):
        assert detect_alphabet("ACGUACGUACGU") == "nucleotide"

    def test_protein_sequence_is_protein(self):
        # MKEFW has F and W which are not IUPAC nucleotide codes
        assert detect_alphabet("MKEFWYLP") == "protein"

    def test_mostly_nucleotide_is_nucleotide(self):
        # 11/12 = 91.7% nucleotide chars → nucleotide
        assert detect_alphabet("ACGTACGTACGM") == "nucleotide"

    def test_below_threshold_is_protein(self):
        # ACGT + MKEFW: 4 nt / 9 total = 44% → protein
        assert detect_alphabet("ACGTMKEFW") == "protein"

    def test_gaps_excluded_from_ratio(self):
        # Gaps stripped before ratio computed; ACGT only → nucleotide
        assert detect_alphabet("A-C-G-T") == "nucleotide"

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            detect_alphabet("")
