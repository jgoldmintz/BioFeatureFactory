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
"""codon_usage: codon context extraction and bicodon edges."""

from biofeaturefactory.codon_usage.codon_usage_pipeline import _codon_context


class TestCodonContext:
    SEQ = "ATGAAAGAATGG"  # ATG AAA GAA TGG

    def test_middle_codon_has_both_bicodons(self):
        codon, fwd, rev = _codon_context(self.SEQ, 1)
        assert codon == "AAA"
        assert fwd == "AAAGAA"   # this codon + next
        assert rev == "ATGAAA"   # previous + this codon

    def test_first_codon_has_no_reverse_bicodon(self):
        codon, fwd, rev = _codon_context(self.SEQ, 0)
        assert codon == "ATG"
        assert fwd == "ATGAAA"
        assert rev is None       # None, not "" -- there is no preceding codon

    def test_last_codon_has_no_forward_bicodon(self):
        codon, fwd, rev = _codon_context(self.SEQ, 3)
        assert codon == "TGG"
        assert fwd is None
        assert rev == "GAATGG"
