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
"""EVmutation: codon alphabet and row-cell formatting.

Importing this module also pins the sys.path shim: the EVmutation clone lives at
<mutation_effects>/EVmutation/, one level above bin/, and anchoring the shim only
at bin/ made `from EVmutation.model import CouplingsModel` raise
ModuleNotFoundError on every invocation.
"""

from biofeaturefactory.mutation_effects.evmutation_pipeline import (
    CODON_ALPHABET,
    _join,
    aa_symbol,
)


class TestModuleImports:
    def test_evmutation_clone_is_on_the_path(self):
        from EVmutation.model import CouplingsModel  # noqa: F401


class TestCodonAlphabet:
    def test_has_sixty_five_symbols(self):
        # 64 codons + one gap symbol; q=65 is what the codon Potts model uses and
        # what the codon capacity tier's memory arithmetic assumes.
        assert len(CODON_ALPHABET) == 65

    def test_symbols_are_unique(self):
        assert len(set(CODON_ALPHABET)) == len(CODON_ALPHABET)

    def test_gap_is_first(self):
        # plmc treats the FIRST alphabet symbol as the gap for -g/--gapignore.
        assert CODON_ALPHABET[0] == "-"

    def test_one_character_per_codon(self):
        assert all(len(c) == 1 for c in CODON_ALPHABET)


class TestRowCells:
    def test_aa_symbol_passthrough(self):
        assert aa_symbol("A") == "A"

    def test_join_uses_comma(self):
        assert _join([1, 2]) == "1,2"

    def test_join_empty_is_empty_string(self):
        # Empty must not render as "0" -- an unmeasurable cell reads as absent.
        assert _join([]) == ""
