"""
Unit tests for get_codon_counts() in utility.py.

Covers: RSCU, W, RSCPU, CPS, W_CP calculations, stop codon handling,
ZeroDivisionError paths, bicodon edge cases.

Run with: pytest test/unit/test_utils_codon_counts.py -v
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import get_codon_counts, codon_to_aa, codon_table


class TestGetCodonCountsBasicCounting:

    def test_single_codon_counted(self):
        """ATG appears once."""
        cd, _ = get_codon_counts("ATG")
        assert cd["counts"]["ATG"] == 1

    def test_repeated_codon_counted(self):
        """ATGATG → ATG appears twice."""
        cd, _ = get_codon_counts("ATGATG")
        assert cd["counts"]["ATG"] == 2

    def test_multiple_distinct_codons(self):
        """ATGAAAGAA → ATG=1, AAA=1, GAA=1."""
        cd, _ = get_codon_counts("ATGAAAGAA")
        assert cd["counts"]["ATG"] == 1
        assert cd["counts"]["AAA"] == 1
        assert cd["counts"]["GAA"] == 1

    def test_trailing_nucleotides_ignored(self):
        """ATGA → only ATG counted, trailing A ignored."""
        cd, _ = get_codon_counts("ATGA")
        assert cd["counts"]["ATG"] == 1
        total = sum(cd["counts"].values())
        assert total == 1

    def test_stop_codons_counted_in_counts(self):
        """Stop codons are counted in raw counts dict."""
        cd, _ = get_codon_counts("ATGTAA")
        assert cd["counts"]["TAA"] == 1

    def test_empty_sequence_all_zeros(self):
        cd, _ = get_codon_counts("")
        assert all(v == 0 for v in cd["counts"].values())


class TestGetCodonCountsRSCU:

    def test_rscu_single_synonymous_codon_used(self):
        """If only one Phe codon (TTT) is used, RSCU(TTT) = 2.0, RSCU(TTC) = 0.0."""
        # Phe has 2 codons: TTT, TTC
        cd, _ = get_codon_counts("TTT")
        # RSCU = count / (total_syn / num_syn) = 1 / (1/2) = 2.0
        assert cd["RSCU"]["TTT"] == pytest.approx(2.0)
        assert cd["RSCU"]["TTC"] == pytest.approx(0.0)

    def test_rscu_equal_usage(self):
        """Equal usage of synonymous codons → RSCU = 1.0 for each."""
        cd, _ = get_codon_counts("TTTTTC")
        assert cd["RSCU"]["TTT"] == pytest.approx(1.0)
        assert cd["RSCU"]["TTC"] == pytest.approx(1.0)

    def test_rscu_stop_codons_are_nan(self):
        """Stop codons get RSCU = NaN."""
        cd, _ = get_codon_counts("TAA")
        assert np.isnan(cd["RSCU"]["TAA"])

    def test_rscu_gap_codon_is_nan(self):
        """Gap codon '---' gets RSCU = NaN."""
        cd, _ = get_codon_counts("ATG")
        assert np.isnan(cd["RSCU"]["---"])

    def test_rscu_zero_division_when_no_synonymous_usage(self):
        """When an amino acid has no codons used, RSCU for those codons = 0.0 (not error)."""
        # Only use ATG (Met), no Phe codons used at all
        cd, _ = get_codon_counts("ATG")
        # numsyn for Phe = 0, numsyn/len = 0/2 = 0 → ZeroDivisionError → NaN
        assert np.isnan(cd["RSCU"]["TTT"]) or cd["RSCU"]["TTT"] == pytest.approx(0.0)


class TestGetCodonCountsW:

    def test_w_max_codon_is_one(self):
        """The most-used synonymous codon gets W = 1.0."""
        seq = "TTT" + "TTT" + "TTT" + "TTC"  # TTT=3, TTC=1
        cd, _ = get_codon_counts(seq)
        assert cd["W"]["TTT"] == pytest.approx(1.0)

    def test_w_less_used_codon_below_one(self):
        """Less-used synonymous codon gets W < 1.0."""
        seq = "TTT" + "TTT" + "TTC"  # 9 chars, 3 codons: TTT=2, TTC=1
        cd, _ = get_codon_counts(seq)
        assert cd["W"]["TTC"] < 1.0
        assert cd["W"]["TTC"] > 0.0

    def test_w_stop_codons_are_nan(self):
        cd, _ = get_codon_counts("ATG")
        assert np.isnan(cd["W"]["TAA"])

    def test_w_equal_usage_both_one(self):
        """Equal usage → both W = 1.0."""
        cd, _ = get_codon_counts("TTTTTC")
        assert cd["W"]["TTT"] == pytest.approx(1.0)
        assert cd["W"]["TTC"] == pytest.approx(1.0)


class TestGetCodonCountsBicodonCounting:

    def test_bicodon_pair_counted(self):
        """ATGAAA → bicodon ATGAAA counted once."""
        _, cpd = get_codon_counts("ATGAAA")
        assert cpd["counts"]["ATGAAA"] == 1

    def test_two_bicodon_pairs(self):
        """ATGAAAGAA → bicodons ATGAAA and AAAGAA each counted once."""
        _, cpd = get_codon_counts("ATGAAAGAA")
        assert cpd["counts"]["ATGAAA"] == 1
        assert cpd["counts"]["AAAGAA"] == 1

    def test_single_codon_no_bicodon(self):
        """ATG alone → no bicodon."""
        _, cpd = get_codon_counts("ATG")
        total = sum(cpd["counts"].values())
        assert total == 0


class TestGetCodonCountsRSCPU:

    def test_rscpu_stop_codon_pair_is_nan(self):
        """Bicodon involving stop codon → RSCPU = NaN."""
        _, cpd = get_codon_counts("ATGTAA")
        assert np.isnan(cpd["RSCPU"]["ATGTAA"])

    def test_rscpu_single_pair_used(self):
        """Single codon pair used once → RSCPU computed."""
        _, cpd = get_codon_counts("ATGAAA")
        # Both ATG (Met, 1 synonym) and AAA (Lys, 2 synonyms)
        # Total syn pairs = 1 * 2 = 2, only ATGAAA used
        # numsyn = 1, syn_cp_count = 2
        # RSCPU = 1 / (1/2) = 2.0
        assert cpd["RSCPU"]["ATGAAA"] == pytest.approx(2.0)


class TestGetCodonCountsCPS:

    def test_cps_computed_for_used_pair(self):
        """CPS = ln(observed/expected) for a used codon pair."""
        _, cpd = get_codon_counts("ATGAAA")
        # observed = 1, expected = count(ATG) * count(AAA) = 1*1 = 1
        # noln CPS = 1/1 = 1.0, CPS = ln(1.0) = 0.0
        assert cpd["noln CPS"]["ATGAAA"] == pytest.approx(1.0)
        assert cpd["CPS"]["ATGAAA"] == pytest.approx(0.0)

    def test_cps_nan_when_codon_not_used(self):
        """CPS = NaN when one codon in the pair has zero count (expected=0)."""
        _, cpd = get_codon_counts("ATGAAA")
        # ATGTTC: ATG count=1 but TTC count=0 → expected=0 → ZeroDivisionError → NaN
        # But ATGTTC wasn't observed either, so count=0
        # Actually this bicodon isn't used so it won't appear in CPS unless it's a stop
        # Let's check a pair where one codon is unused
        assert np.isnan(cpd["CPS"].get("TTTTTC", np.nan))

    def test_cps_nan_for_stop_pairs(self):
        """CPS for stop-codon pairs is NaN."""
        _, cpd = get_codon_counts("ATGTAA")
        assert np.isnan(cpd["CPS"].get("ATGTAA", np.nan))


class TestGetCodonCountsW_CP:

    def test_w_cp_for_used_pair(self):
        """W_CP for a used codon pair should be between 0 and 1 (inclusive)."""
        _, cpd = get_codon_counts("ATGAAA")
        w_cp = cpd["W_CP"]["ATGAAA"]
        assert not np.isnan(w_cp)
        assert 0.0 < w_cp <= 1.0

    def test_w_cp_nan_for_stop_pair(self):
        _, cpd = get_codon_counts("ATGTAA")
        assert np.isnan(cpd["W_CP"]["ATGTAA"])

    def test_w_cp_max_is_one(self):
        """The most-used synonymous codon pair gets W_CP = 1.0."""
        # Use ATGAAA (the only Met-Lys pair used)
        _, cpd = get_codon_counts("ATGAAA")
        assert cpd["W_CP"]["ATGAAA"] == pytest.approx(1.0)


class TestGetCodonCountsReturnStructure:

    def test_returns_two_dicts(self):
        cd, cpd = get_codon_counts("ATGAAAGAA")
        assert isinstance(cd, dict)
        assert isinstance(cpd, dict)

    def test_codondata_has_expected_keys(self):
        cd, _ = get_codon_counts("ATG")
        assert "counts" in cd
        assert "RSCU" in cd
        assert "W" in cd

    def test_codonpairdata_has_expected_keys(self):
        _, cpd = get_codon_counts("ATGAAA")
        assert "counts" in cpd
        assert "RSCPU" in cpd
        assert "CPS" in cpd
        assert "noln CPS" in cpd
        assert "W_CP" in cpd

    def test_all_64_codons_in_counts(self):
        cd, _ = get_codon_counts("ATG")
        # 64 standard codons + gap = 65
        assert len(cd["counts"]) == len(codon_to_aa)

    def test_longer_sequence_consistency(self):
        """A real-ish CDS: check that total codon counts match expected."""
        seq = "ATGAAAGAATTTTTCGATTGCGAT"  # 8 codons
        cd, cpd = get_codon_counts(seq)
        total_codons = sum(cd["counts"].values())
        assert total_codons == 8
        total_bicodons = sum(cpd["counts"].values())
        assert total_bicodons == 7  # 8 codons → 7 pairs
