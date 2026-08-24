# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Tests for the additive length-aware variant core.

The load-bearing test in this file is TestLegacyEquivalence: it asserts the new
parser agrees with the existing SNV parsers on every legacy token, which is what
makes the addition safe to land without touching them.
"""

import pytest

from biofeaturefactory.lib.utility import (
    Variant,
    _translate_codons,
    protein_consequence,
    canonical_token,
    get_mutation_data,
    get_mutation_data_bioAccurate,
    parse_variant,
    rc_variant,
    revcomp_seq,
    splice_seq,
    update_str,
)

# Real 27 bp window carrying the CAACAA tandem, from the Enh13 case. Index map:
# 0:G 1:C 2:A 3:A 4:G 5:C 6:A 7:A 8:A 9:C 10:C 11:A 12:C 13:A 14:A 15:C 16:A
# 17:A 18:T 19:G 20:G 21:T 22:C 23:A 24:G 25:A 26:C
TANDEM = "GCAAGCAAACCACAACAATGGTCAGAC"


class TestVariantRecord:

    def test_pos_is_one_based_and_pos0_slices(self):
        v = Variant(pos=1, ref="A", alt="G")
        assert v.pos0 == 0

    def test_rejects_zero_position(self):
        with pytest.raises(ValueError, match="1-based"):
            Variant(pos=0, ref="A", alt="G")

    def test_rejects_empty_allele(self):
        with pytest.raises(ValueError, match="non-empty"):
            Variant(pos=5, ref="ACAA", alt="")

    def test_rejects_unknown_orientation(self):
        with pytest.raises(ValueError, match="orientation"):
            Variant(pos=5, ref="A", alt="G", orientation="plus")

    def test_hashable_so_it_can_key_a_dict(self):
        v = Variant(pos=5, ref="A", alt="G")
        assert {v: "x"}[Variant(pos=5, ref="A", alt="G")] == "x"

    @pytest.mark.parametrize("ref,alt,kind,delta", [
        ("A", "G", "snv", 0),
        ("AC", "GT", "mnv", 0),
        ("A", "ACAA", "insertion", 3),
        ("ACAA", "A", "deletion", -3),
        ("ACAA", "GT", "delins", -2),
    ])
    def test_kind_and_length_delta(self, ref, alt, kind, delta):
        v = Variant(pos=10, ref=ref, alt=alt)
        assert v.kind == kind
        assert v.length_delta == delta
        assert v.is_snv is (kind == "snv")


class TestLegacyEquivalence:
    """Every legacy SNV token must parse to the identical (pos, ref, alt)."""

    LEGACY = ["G1A", "A2T", "C123G", "T3435A", "G1234A", "A999C"]

    @pytest.mark.parametrize("token", LEGACY)
    def test_agrees_with_bioaccurate_one_based(self, token):
        pos, (ref, alt) = get_mutation_data_bioAccurate(token, is_nt=True)
        v = parse_variant(token, is_nt=True)
        assert (v.pos, v.ref, v.alt) == (pos, ref, alt)

    @pytest.mark.parametrize("token", LEGACY)
    def test_agrees_with_get_mutation_data_zero_based(self, token):
        pos0, (ref, alt) = get_mutation_data(token)
        v = parse_variant(token, is_nt=True)
        assert (v.pos0, v.ref, v.alt) == (pos0, ref, alt)

    @pytest.mark.parametrize("token", LEGACY)
    def test_splice_seq_matches_update_str_on_snvs(self, token):
        seq = "ACGT" * 1000
        v = parse_variant(token, is_nt=True)
        assert splice_seq(seq, v.pos0, v.ref, v.alt, validate=False) == \
            update_str(seq, v.alt, v.pos0)


class TestParseVariant:

    @pytest.mark.parametrize("token,pos,ref,alt", [
        ("ACAA112217430A", 112217430, "ACAA", "A"),
        ("A112217430ACAA", 112217430, "A", "ACAA"),
        ("AT100G", 100, "AT", "G"),
        ("g123a", 123, "G", "A"),
    ])
    def test_parses_non_snv_tokens(self, token, pos, ref, alt):
        v = parse_variant(token, is_nt=True)
        assert (v.pos, v.ref, v.alt) == (pos, ref, alt)

    @pytest.mark.parametrize("token", [
        "R213W",        # amino acid on the nt path
        "175Stop",      # stop codon
        "A5X",          # off-alphabet alt
        "X5A",          # off-alphabet ref
        "123",          # no alleles
        "AG",           # no position
        "A0G",          # position below 1
        "",
        None,
    ])
    def test_returns_none_never_raises(self, token):
        assert parse_variant(token, is_nt=True) is None

    def test_amino_acid_path_accepts_aa_tokens(self):
        v = parse_variant("R213W", is_nt=False)
        assert (v.pos, v.ref, v.alt) == (213, "R", "W")

    def test_carries_gene_and_orientation(self):
        v = parse_variant("G1A", is_nt=True, gene="TP53", orientation="transcript")
        assert v.gene == "TP53" and v.orientation == "transcript"


class TestSpliceSeq:

    def test_deletion_shortens(self):
        assert splice_seq("ACAATGG", 0, "ACAA", "A") == "ATGG"

    def test_insertion_lengthens(self):
        assert splice_seq("ATGG", 0, "A", "ACAA") == "ACAATGG"

    def test_rejects_ref_mismatch_by_default(self):
        with pytest.raises(ValueError, match="REF mismatch"):
            splice_seq("ACGT", 0, "TT", "T")

    def test_validate_false_skips_the_check(self):
        assert splice_seq("ACGT", 0, "TT", "T", validate=False) == "TGT"

    def test_rejects_span_past_end(self):
        with pytest.raises(ValueError, match="only"):
            splice_seq("ACGT", 3, "TAAA", "T")

    def test_rejects_negative_position(self):
        with pytest.raises(ValueError, match=">= 0"):
            splice_seq("ACGT", -1, "A", "G")

    def test_u_and_t_compare_equal(self):
        assert splice_seq("ACGU", 3, "T", "A") == "ACGA"

    def test_round_trips_through_the_record(self):
        v = parse_variant("ACAA12A", is_nt=True)
        assert splice_seq(TANDEM, v.pos0, v.ref, v.alt) == \
            "GCAAGCAAACCACAATGGTCAGAC"


class TestCanonicalToken:
    """The pkey-collision guard: one edit must yield one key."""

    def test_two_spellings_of_one_deletion_collapse(self):
        # Same 3 bp deletion inside the CAACAA tandem, written at two different
        # anchors. Both produce the identical sequence.
        left = Variant(pos=12, ref="ACAA", alt="A")
        right = Variant(pos=15, ref="ACAA", alt="A")
        assert splice_seq(TANDEM, left.pos0, left.ref, left.alt) == \
            splice_seq(TANDEM, right.pos0, right.ref, right.alt)
        # Without the reference they stay distinct -- representation only.
        assert canonical_token(left) != canonical_token(right)
        # With it they collapse to one key.
        assert canonical_token(left, seq=TANDEM) == \
            canonical_token(right, seq=TANDEM) == "CACA11C"

    def test_snv_token_is_unchanged(self):
        assert canonical_token(Variant(pos=123, ref="G", alt="A")) == "G123A"

    def test_uppercases(self):
        assert canonical_token(Variant(pos=1, ref="g", alt="a")) == "g1a".upper()

    def test_trims_shared_flanks(self):
        # AGT/AGA at 10 is really T/A at 12.
        assert canonical_token(Variant(pos=10, ref="AGT", alt="AGA")) == "T12A"

    def test_keeps_the_vcf_anchor_base(self):
        # ACAA/A must NOT collapse to an empty allele.
        assert canonical_token(Variant(pos=12, ref="ACAA", alt="A")) == "ACAA12A"

    def test_out_of_range_seq_leaves_representation_alone(self):
        v = Variant(pos=900, ref="ACAA", alt="A")
        assert canonical_token(v, seq=TANDEM) == "ACAA900A"


class TestProteinConsequence:
    """Codon-aware protein effect. ORF translates to M K E W L T C D *."""

    ORF = "ATGAAAGAATGGCTGACCTGTGATTAA"

    def _c(self, token):
        v = parse_variant(token, is_nt=True)
        assert v is not None, token
        return protein_consequence(v, self.ORF)

    def test_orf_translates_with_single_char_stop(self):
        # codon_to_aa maps a stop to the 4-char 'Stop'; per-residue indexing and
        # len() require one character, so it must normalize to '*'.
        assert _translate_codons(self.ORF) == "MKEWLTCD*"

    def test_snv_output_is_byte_identical_to_single_residue(self):
        r = self._c("A4G")
        assert (r["wt_aa"], r["mut_aa"]) == ("K", "E")
        assert (r["n_aa_wt"], r["n_aa_mut"]) == (1, 1)
        assert r["aa_consequence"] == "snv"
        assert r["aa_pos"] == 2

    def test_inframe_deletion_empties_the_mutant_allele(self):
        r = self._c("AAAG4A")
        assert (r["wt_aa"], r["mut_aa"]) == ("E", "")
        assert (r["n_aa_wt"], r["n_aa_mut"]) == (1, 0)
        assert r["aa_consequence"] == "inframe_del"

    def test_inframe_insertion_empties_the_wt_allele(self):
        r = self._c("A4AGAA")
        assert (r["wt_aa"], r["mut_aa"]) == ("", "R")
        assert r["aa_consequence"] == "inframe_ins"

    def test_multi_residue_span_is_a_string_not_a_list(self):
        r = self._c("AAAGAAT4A")
        assert r["wt_aa"] == "KEW" and isinstance(r["wt_aa"], str)
        assert r["n_aa_wt"] == 3
        assert r["aa_consequence"] == "inframe_delins"

    def test_frameshift_leaves_both_alleles_empty(self):
        r = self._c("AAG4A")
        assert (r["wt_aa"], r["mut_aa"]) == ("", "")
        assert (r["n_aa_wt"], r["n_aa_mut"]) == (0, 0)
        assert r["aa_consequence"] == "frameshift"
        assert r["new_stop_aa_pos"] is not None

    def test_stop_gained_detected(self):
        r = self._c("TGG10TAG")
        assert r["aa_consequence"] == "stop_gained"
        assert r["mut_aa"] == "*"
        assert r["new_stop_aa_pos"] == 4

    def test_synonymous_trims_to_empty(self):
        assert self._c("AAA4AAG")["aa_consequence"] == "synonymous"

    def test_returns_none_outside_the_orf(self):
        assert protein_consequence(Variant(pos=900, ref="A", alt="G"), self.ORF) is None

    def test_every_allele_field_is_a_string(self):
        """The whole point of the string convention: one dtype in every row."""
        for token in ["A4G", "AAAG4A", "A4AGAA", "AAG4A", "AAAGAAT4A"]:
            r = self._c(token)
            assert isinstance(r["wt_aa"], str) and isinstance(r["mut_aa"], str), token
            assert isinstance(r["n_aa_wt"], int) and isinstance(r["n_aa_mut"], int), token


class TestRevcomp:

    def test_multibase_is_reversed_and_complemented(self):
        assert revcomp_seq("ACAA") == "TTGT"

    def test_single_base(self):
        assert revcomp_seq("A") == "T"

    def test_raises_on_unknown_base(self):
        with pytest.raises(ValueError, match="cannot complement"):
            revcomp_seq("ACXA")

    def test_rc_variant_flips_alleles_and_orientation(self):
        v = Variant(pos=100, ref="ACAA", alt="A", orientation="transcript")
        r = rc_variant(v)
        assert (r.ref, r.alt, r.orientation) == ("TTGT", "T", "genomic")

    def test_rc_variant_corrects_multibase_start(self):
        # genomic_pos is the genomic coordinate of the transcript-FIRST ref base,
        # which on the minus strand is the rightmost of a 4 bp span.
        v = Variant(pos=50, ref="ACAA", alt="A", orientation="transcript")
        assert rc_variant(v, genomic_pos=1000).pos == 997

    def test_rc_variant_snv_needs_no_correction(self):
        v = Variant(pos=50, ref="A", alt="G", orientation="transcript")
        assert rc_variant(v, genomic_pos=1000).pos == 1000
