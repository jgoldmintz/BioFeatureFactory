"""
Unit tests for MSA diversity metrics in utility.py:
compute_sequence_weights, compute_neff

Run with: pytest test/unit/test_utils_neff.py -v
"""

import sys
from pathlib import Path

import pytest

from biofeaturefactory.lib.utility import compute_sequence_weights, compute_neff



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



# codon_mode behavior


class TestCodonMode:
    """Tests for codon_mode=True in compute_sequence_weights / compute_neff.

    Codon-level identity is the right comparison unit for q=64 codon Potts
    fits. The codon_mode=True path chunks the sequence into 3-nt triplets
    and treats any triplet containing a gap character as a gap codon ('---').
    """

    def test_codon_mode_single_codon_weight_is_one(self):
        weights = compute_sequence_weights({"seq1": "AAA"}, codon_mode=True)
        assert weights["seq1"] == pytest.approx(1.0)

    def test_codon_mode_identical_codon_msas_cluster(self):
        msa = {"s1": "AAACCC", "s2": "AAACCC"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8, codon_mode=True)
        assert weights["s1"] == pytest.approx(0.5)
        assert weights["s2"] == pytest.approx(0.5)

    def test_codon_mode_one_codon_diff_below_threshold(self):
        # 4 codons, 1 differs (AAA vs AAC) -> 75% codon identity < 0.8 -> not neighbors.
        # Note: at the nt level, 11/12 = 91.7% identity, which WOULD cluster.
        msa = {"s1": "AAACCCGGGTTT", "s2": "AACCCCGGGTTT"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8, codon_mode=True)
        assert weights["s1"] == pytest.approx(1.0)
        assert weights["s2"] == pytest.approx(1.0)

    def test_codon_mode_normalizes_partial_gap_codons(self):
        # 'A-G' contains a gap: it must be treated as a gap codon ('---'),
        # not as an informative residue. Before the fix, identical partial-
        # gap codons would falsely cluster these sequences.
        msa = {"s1": "A-GA-G", "s2": "A-GA-G"}
        weights = compute_sequence_weights(msa, identity_threshold=0.8, codon_mode=True)
        # All codons are gap-containing -> 0 informative aligned positions ->
        # no clustering -> each gets full weight.
        assert weights["s1"] == pytest.approx(1.0)
        assert weights["s2"] == pytest.approx(1.0)

    def test_codon_mode_neff_diverges_from_nt_mode_for_synonymous_diff(self):
        # Two seqs differing by one synonymous nucleotide change inside one codon:
        #   nt-mode  : 5/6 = 83.3% identity > 0.8  -> cluster -> N_eff = 1.0
        #   codon-mode: 1/2 = 50%  identity < 0.8  -> no cluster -> N_eff = 2.0
        # The two modes give different answers; the codon-mode one is the
        # right one for q=64 codon Potts fits.
        msa = {"s1": "AAAGGG", "s2": "AAGGGG"}
        nt_neff = compute_neff(msa, identity_threshold=0.8, codon_mode=False)
        codon_neff = compute_neff(msa, identity_threshold=0.8, codon_mode=True)
        assert nt_neff == pytest.approx(1.0)
        assert codon_neff == pytest.approx(2.0)

    def test_codon_mode_empty_msa_returns_zero(self):
        assert compute_neff({}, codon_mode=True) == pytest.approx(0.0)

    def test_codon_mode_default_is_false(self):
        # Backward compatibility: callers that don't pass codon_mode get
        # the original nt-level behavior.
        msa = {"s1": "AAAGGG", "s2": "AAGGGG"}
        default_neff = compute_neff(msa, identity_threshold=0.8)
        nt_neff = compute_neff(msa, identity_threshold=0.8, codon_mode=False)
        assert default_neff == pytest.approx(nt_neff)
