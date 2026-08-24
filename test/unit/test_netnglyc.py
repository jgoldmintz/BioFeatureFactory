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
"""netNglyc: gene/mutation id normalisation, AA consequence classification."""

from biofeaturefactory.netNglyc.netnglyc_pipeline import (
    _aa_consequence_from_alleles,
    normalize_gene_symbol,
    normalize_mutation_id,
)


class TestNormalizeGeneSymbol:
    def test_uppercases_and_strips(self):
        assert normalize_gene_symbol(" brca1 ") == "BRCA1"

    def test_already_normal_is_unchanged(self):
        assert normalize_gene_symbol("BRCA1") == "BRCA1"


class TestNormalizeMutationId:
    def test_plain_token_is_preserved(self):
        assert normalize_mutation_id("A123G") == "A123G"


class TestAaConsequence:
    def test_substitution(self):
        assert _aa_consequence_from_alleles("K", "R") == "snv"

    def test_stop_gained(self):
        assert _aa_consequence_from_alleles("K", "*") == "stop_gained"

    def test_synonymous_is_not_reported_as_snv(self):
        # Same residue in and out is not a protein-level substitution.
        assert _aa_consequence_from_alleles("K", "K") != "snv"
