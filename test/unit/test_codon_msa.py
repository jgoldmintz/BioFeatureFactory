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
"""codon_msa: codon translation and protein-id extraction."""

from biofeaturefactory.core.codon_msa_pipeline import (
    _extract_protein_ids,
    _translate_codon,
)


class TestTranslateCodon:
    def test_start_codon(self):
        assert _translate_codon("ATG") == "M"

    def test_stop_codon_is_star(self):
        assert _translate_codon("TGA") == "*"

    def test_ambiguous_codon_is_x(self):
        assert _translate_codon("NNN") == "X"

    def test_truncated_codon_is_x_not_an_exception(self):
        assert _translate_codon("at") == "X"


class TestExtractProteinIds:
    def test_finds_multiple_accessions(self):
        assert _extract_protein_ids("XP_1.1 and NP_2.2") == ["XP_1.1", "NP_2.2"]

    def test_falls_back_to_first_token_when_no_accession_matches(self):
        """A header with no RefSeq-style accession still yields an id.

        The fallback is text.split()[0], so a non-RefSeq header is keyed on its
        first whitespace field rather than dropped. Only an EMPTY string yields
        no id at all.
        """
        assert _extract_protein_ids("no accessions here") == ["no"]
        assert _extract_protein_ids("") == []

    def test_uniref90_prefix_is_stripped(self):
        assert _extract_protein_ids("UniRef90_P12345") == ["P12345"]
