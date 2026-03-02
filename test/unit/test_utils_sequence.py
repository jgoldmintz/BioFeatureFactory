"""
Unit tests for sequence and coordinate utility functions in utility.py.

Run with: pytest test/unit/test_utils_sequence.py -v
"""

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    update_str,
    subseq,
    get_mutation_data,
    get_mutation_data_bioAccurate,
    convert_position,
    translate_orf_sequence,
    infer_aamutation_from_nt,
)



# update_str


class TestUpdateStr:

    def test_replace_first_char(self):
        assert update_str("ACGT", "T", 0) == "TCGT"

    def test_replace_last_char(self):
        assert update_str("ACGT", "A", 3) == "ACGA"

    def test_replace_middle(self):
        assert update_str("ACGT", "X", 2) == "ACXT"

    def test_original_unchanged(self):
        s = "ACGT"
        update_str(s, "T", 0)
        assert s == "ACGT"



# subseq


class TestSubseq:

    def test_centered_window(self):
        # pos=4, l=3 → half=1, start=3, end=6
        assert subseq("ACGTACGT", 4, 3) == "TAC"

    def test_truncates_at_left_boundary(self):
        # pos=0, l=3 → half=1, start=0, end=2 (clipped)
        assert subseq("ACGT", 0, 3) == "AC"

    def test_truncates_at_right_boundary(self):
        # pos=3, l=3 → half=1, start=2, end=4 (len=4, clipped)
        assert subseq("ACGT", 3, 3) == "GT"

    def test_window_length_one(self):
        assert subseq("ACGT", 2, 1) == "G"

    def test_full_window_returned_when_centered(self):
        # pos=3, l=3 in seq of length 7 → no clipping
        result = subseq("ACGTACG", 3, 3)
        assert len(result) == 3

    def test_invalid_even_length_raises(self):
        with pytest.raises(AssertionError):
            subseq("ACGT", 2, 4)

    def test_out_of_range_pos_raises(self):
        with pytest.raises(AssertionError):
            subseq("ACGT", 10, 3)



# get_mutation_data


class TestGetMutationData:

    def test_parses_ref_and_alt(self):
        pos, (ref, alt) = get_mutation_data("G123A")
        assert ref == "G"
        assert alt == "A"

    def test_position_above_one_is_zero_based(self):
        pos, _ = get_mutation_data("G2A")
        assert pos == 1

    def test_position_100(self):
        pos, _ = get_mutation_data("A100T")
        assert pos == 99

    def test_position_one_behavior(self):
        # position 1 returns 1, not 0 — documents current behavior
        pos, _ = get_mutation_data("G1A")
        assert pos == 1

    def test_multi_digit_position(self):
        pos, (ref, alt) = get_mutation_data("C1234G")
        assert pos == 1233
        assert ref == "C"
        assert alt == "G"



# get_mutation_data_bioAccurate


class TestGetMutationDataBioAccurate:

    def test_parses_normally(self):
        pos, (ref, alt) = get_mutation_data_bioAccurate("G123A")
        assert pos == 123
        assert ref == "G"
        assert alt == "A"

    def test_position_one_returns_one(self):
        pos, _ = get_mutation_data_bioAccurate("G1A")
        assert pos == 1

    def test_stop_codon_returns_none(self):
        pos, nts = get_mutation_data_bioAccurate("G123Stop")
        assert pos is None
        assert nts is None

    def test_sto_variant_returns_none(self):
        pos, nts = get_mutation_data_bioAccurate("A45Sto")
        assert pos is None
        assert nts is None



# convert_position


class TestConvertPosition:

    def test_identity_no_gaps(self):
        pos, err = convert_position("ACGT", "ACGT", 1)
        assert pos == 1
        assert err is None

    def test_position_2_no_gaps(self):
        pos, err = convert_position("ACGT", "ACGT", 2)
        assert pos == 2
        assert err is None

    def test_gap_in_seq1_shifts_seq2_position(self):
        # seq1: A - C G T  (C is 2nd non-gap residue)
        # seq2: A A C G T  (C is 3rd non-gap residue)
        pos, err = convert_position("A-CGT", "AACGT", 2)
        assert pos == 3
        assert err is None

    def test_seq1_position_aligns_with_gap_in_seq2(self):
        # seq1: A C G T
        # seq2: A C - T  (seq1 pos 3 = G aligns with gap in seq2)
        pos, err = convert_position("ACGT", "AC-T", 3)
        assert err is not None
        assert "gap" in err.lower()

    def test_zero_position_sets_error(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, err = convert_position("ACGT", "ACGT", 0)
        assert err is not None



# translate_orf_sequence


class TestTranslateOrfSequence:

    def test_basic_translation(self):
        # ATG (M) + AAA (K) + TAA (stop)
        result = translate_orf_sequence("ATGAAATAA")
        assert result == "MK"

    def test_stop_codon_stripped(self):
        result = translate_orf_sequence("ATGAAATAA")
        assert "*" not in result

    def test_empty_input_returns_empty(self):
        assert translate_orf_sequence("") == ""

    def test_rna_u_converted(self):
        # U should be treated as T
        result = translate_orf_sequence("AUGAAAUAA")
        assert result == "MK"

    def test_lowercase_input_handled(self):
        result = translate_orf_sequence("atgaaataa")
        assert result == "MK"

    def test_incomplete_trailing_codon_handled(self):
        # 10 nt: ATG AAA TAA + G (trailing incomplete codon)
        result = translate_orf_sequence("ATGAAATAAG")
        assert result == "MK"



# infer_aamutation_from_nt


class TestInferAamutationFromNt:

    # ATG AAA GAA TGG = M  K  E  W
    SEQ = "ATGAAAGAATGG"

    def test_synonymous_mut_same_aa(self):
        # codon 2 = AAA (K); 3rd base A6→G gives AAG (still K)
        result = infer_aamutation_from_nt("A6G", self.SEQ)
        assert result is not None
        aa_pos, wt_aa, mut_aa = result
        assert wt_aa == mut_aa == "K"

    def test_missense_mut_changes_aa(self):
        # codon 2 = AAA (K); 1st base A4→G gives GAA (E)
        result = infer_aamutation_from_nt("A4G", self.SEQ)
        assert result is not None
        aa_pos, wt_aa, mut_aa = result
        assert wt_aa == "K"
        assert mut_aa == "E"

    def test_stop_codon_returns_none(self):
        result = infer_aamutation_from_nt("A4Stop", self.SEQ)
        assert result is None
