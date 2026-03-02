"""
Unit tests for MSA diversity metrics in utility.py:
compute_sequence_weights, compute_neff

Run with: pytest test/unit/test_utils_neff.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import compute_sequence_weights, compute_neff



# compute_sequence_weights


class TestComputeSequenceWeights:

    def test_single_sequence_weight_is_one(self):
        weights = compute_sequence_weights({"seq1": "ACGT"})
        assert weights["seq1"] == pytest.approx(1.0)

    def test_identical_sequences_each_half_weight(self):
        msa = {"seq1": "ACGT", "seq2": "ACGT"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8)
        assert weights["seq1"] == pytest.approx(0.5)
        assert weights["seq2"] == pytest.approx(0.5)

    def test_completely_different_sequences_each_full_weight(self):
        # AAAA vs CCCC: 0% identity → not neighbors
        msa = {"seq1": "AAAA", "seq2": "CCCC"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8)
        assert weights["seq1"] == pytest.approx(1.0)
        assert weights["seq2"] == pytest.approx(1.0)

    def test_empty_msa_returns_empty(self):
        assert compute_sequence_weights({}) == {}

    def test_weights_sum_equals_neff(self):
        msa = {"s1": "ACGT", "s2": "ACGT", "s3": "TTTT"}
        weights = compute_sequence_weights(msa)
        neff = compute_neff(msa)
        assert sum(weights.values()) == pytest.approx(neff)

    def test_all_keys_present_in_output(self):
        msa = {"s1": "ACGT", "s2": "ACGC", "s3": "TTTT"}
        weights = compute_sequence_weights(msa)
        assert set(weights.keys()) == {"s1", "s2", "s3"}

    def test_gap_positions_excluded_from_identity(self):
        # Both sequences have gaps at the same positions — gaps excluded
        # s1: A-GT, s2: A-GT → aligned non-gap positions: A,G,T (3 matches / 3) = 100%
        msa = {"s1": "A-GT", "s2": "A-GT"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8)
        assert weights["s1"] == pytest.approx(0.5)
        assert weights["s2"] == pytest.approx(0.5)



# compute_neff


class TestComputeNeff:

    def test_single_sequence_neff_is_one(self):
        assert compute_neff({"seq1": "ACGT"}) == pytest.approx(1.0)

    def test_two_identical_sequences_neff_is_one(self):
        # Each has weight 0.5 → neff = 1.0
        msa = {"s1": "ACGT", "s2": "ACGT"}
        assert compute_neff(msa) == pytest.approx(1.0)

    def test_two_diverse_sequences_neff_is_two(self):
        msa = {"s1": "AAAA", "s2": "CCCC"}
        assert compute_neff(msa) == pytest.approx(2.0)

    def test_neff_increases_with_diversity(self):
        identical = {"s1": "ACGT", "s2": "ACGT", "s3": "ACGT"}
        diverse = {"s1": "AAAA", "s2": "CCCC", "s3": "GGGG"}
        assert compute_neff(diverse) > compute_neff(identical)

    def test_neff_never_exceeds_sequence_count(self):
        msa = {"s1": "AAAA", "s2": "CCCC", "s3": "GGGG"}
        assert compute_neff(msa) <= len(msa)

    def test_empty_msa_returns_zero(self):
        assert compute_neff({}) == pytest.approx(0.0)
