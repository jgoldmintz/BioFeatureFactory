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
"""adabmDCA: alphabet symbols, row joining, zero-sum gauge."""

import numpy as np

from biofeaturefactory.mutation_effects.adabmdca_pipeline import (
    _join,
    aa_symbol,
    apply_zero_sum_gauge,
)


class TestAaSymbol:
    def test_standard_residue(self):
        assert aa_symbol("A") == "A"

    def test_stop_and_gap_pass_through(self):
        assert aa_symbol("*") == "*"
        assert aa_symbol("-") == "-"


class TestJoin:
    def test_joins_with_comma(self):
        assert _join([1, 2]) == "1,2"

    def test_empty_is_empty_string_not_none(self):
        # An empty cell must be "" so a consumer reads it as absent, not as 0.
        assert _join([]) == ""


class TestZeroSumGauge:
    def test_h_rows_sum_to_zero(self):
        rng = np.random.default_rng(0)
        L, q = 5, 4
        h = rng.normal(size=(L, q))
        J = rng.normal(size=(L, q, L, q))
        h_g, _ = apply_zero_sum_gauge(h, J)
        assert np.allclose(h_g.sum(axis=1), 0.0, atol=1e-12)

    def test_coupling_rows_and_columns_sum_to_zero(self):
        rng = np.random.default_rng(1)
        L, q = 4, 3
        h = rng.normal(size=(L, q))
        J = rng.normal(size=(L, q, L, q))
        _, J_g = apply_zero_sum_gauge(h, J)
        # dim 1 = a (left symbol), dim 3 = b (right symbol)
        assert np.allclose(J_g.sum(axis=1), 0.0, atol=1e-12)
        assert np.allclose(J_g.sum(axis=3), 0.0, atol=1e-12)

    def test_gauge_is_idempotent(self):
        rng = np.random.default_rng(2)
        L, q = 4, 3
        h = rng.normal(size=(L, q))
        J = rng.normal(size=(L, q, L, q))
        h1, J1 = apply_zero_sum_gauge(h, J)
        h2, J2 = apply_zero_sum_gauge(h1, J1)
        assert np.allclose(h1, h2, atol=1e-12)
        assert np.allclose(J1, J2, atol=1e-12)
