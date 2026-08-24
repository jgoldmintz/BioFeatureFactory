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
"""miranda: gene~substrate keys, distance weighting, safe division."""

import pytest

from biofeaturefactory.miranda.miranda_ensemble import (
    _gene_key,
    _split_gene_key,
    _substrate_of,
    distance_to_variant,
    exp_weight,
    safe_div,
)


class TestGeneKey:
    def test_round_trips(self):
        key = _gene_key("BRCA1", "utr3")
        assert _split_gene_key(key) == ("BRCA1", "utr3")
        assert _substrate_of(key) == "utr3"

    def test_substrate_is_uppercased_in_the_key(self):
        # The KEY carries an upper-cased substrate; _split_gene_key returns the
        # lower-cased form, so the two are not string-identical by design.
        assert _gene_key("BRCA1", "utr3") == "BRCA1~UTR3"

    def test_separator_is_tilde_not_dash(self):
        # A '-' would collide with pkey ({GENE}-{sha}); '~' cannot appear in either.
        assert "~" in _gene_key("BRCA1", "utr3")
        assert "-" not in _gene_key("BRCA1", "utr3")


class TestSafeDiv:
    def test_normal_division(self):
        assert safe_div(3.0, 2.0) == 1.5

    def test_zero_denominator_returns_zero_not_raises(self):
        assert safe_div(1.0, 0.0) == 0.0


class TestDistanceToVariant:
    def test_distance_is_absolute(self):
        assert distance_to_variant(100, 120, 1) == 20

    def test_inside_the_span_is_zero(self):
        assert distance_to_variant(10, 10, 3) == 0

    def test_missing_site_position_returns_none(self):
        # None means "not measurable" and must not be read as distance 0.
        assert distance_to_variant(None, 5, 1) is None


class TestExpWeight:
    def test_zero_distance_is_full_weight(self):
        assert exp_weight(0, 0.1) == 1.0

    def test_missing_distance_is_full_weight(self):
        assert exp_weight(None, 0.1) == 1.0

    def test_weight_decays_with_distance(self):
        assert exp_weight(50, 0.1) < exp_weight(5, 0.1) <= 1.0
