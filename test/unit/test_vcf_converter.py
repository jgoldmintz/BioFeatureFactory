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

"""vcf_converter: chromosome naming, ORF->genomic projection, VCF row shape."""

import pytest

from biofeaturefactory.core.vcf_converter import (
    chromosome_to_refseq,
    format_chromosome,
    reverse_complement,
)


class TestChromosomeNaming:
    def test_refseq_for_autosome(self):
        assert format_chromosome("5", "refseq") == "NC_000005.10"

    def test_ucsc_prefix(self):
        assert format_chromosome("5", "ucsc") == "chr5"

    def test_chromosome_to_refseq_matches_format(self):
        assert chromosome_to_refseq("5") == format_chromosome("5", "refseq")

    def test_refseq_is_zero_padded_to_six_digits(self):
        # NC_0000<NN>.<ver>; a naive f"NC_{chrom}" would give NC_5
        assert format_chromosome("1", "refseq").startswith("NC_000001.")
        assert format_chromosome("22", "refseq").startswith("NC_000022.")


class TestReverseComplement:
    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("g") == "C"

    def test_multibase_is_returned_unchanged(self):
        """Documents a LIMIT, not desired behaviour.

        The implementation is comp.get(nucleotide.upper(), nucleotide) -- a
        single-base lookup. Any multi-base string misses the dict and comes back
        untouched, so this is a base complement, not a reverse complement:
        'ATGC' should be 'GCAT'. Harmless today only because impact analysis
        reports ZERO callers. If a caller ever appears -- non-SNV REF/ALT on a
        minus-strand gene is the obvious one -- it needs a real implementation
        first, and this test should start asserting 'GCAT'.
        """
        assert reverse_complement("ATGC") == "ATGC"
        assert reverse_complement("AAAA") == "AAAA"

    def test_unknown_character_passes_through(self):
        assert reverse_complement("N") == "N"
