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
"""rare_codon: Poisson-binomial tail probability.

rare_codon_pipeline.py imports its cg_cotrans siblings by bare name, so it is
loaded from its own directory the way the pipeline itself is invoked rather than
as a package module.
"""

import importlib.util
import itertools
import sys
from pathlib import Path

import pytest

_RC_DIR = Path(__file__).resolve().parents[2] / "biofeaturefactory" / "rare_codon"
if str(_RC_DIR) not in sys.path:
    sys.path.insert(0, str(_RC_DIR))
_spec = importlib.util.spec_from_file_location("rc_pipeline", _RC_DIR / "rare_codon_pipeline.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

prob_ntuple = rc._prob_ntuple_poisson_binomial


def _brute_force(p, n):
    """P(X >= n) by explicit enumeration -- the 2**N form this DP replaced."""
    total = 0.0
    for mask in itertools.product([0, 1], repeat=len(p)):
        if sum(mask) < n:
            continue
        w = 1.0
        for bit, pi in zip(mask, p):
            w *= pi if bit else (1.0 - pi)
        total += w
    return total


class TestProbNtupleEdges:
    def test_n_zero_is_certain(self):
        assert prob_ntuple([0.5] * 4, 0) == 1.0

    def test_n_greater_than_population_is_impossible(self):
        assert prob_ntuple([0.5] * 4, 5) == 0.0

    def test_certain_success(self):
        assert prob_ntuple([1.0, 0.0, 0.0], 1) == pytest.approx(1.0)

    def test_two_fair_coins(self):
        assert prob_ntuple([0.5, 0.5], 2) == pytest.approx(0.25)
        assert prob_ntuple([0.5, 0.5], 1) == pytest.approx(0.75)


class TestProbNtupleMatchesEnumeration:
    """The DP replaced a 2**N enumeration that took 19h21m at N=35.

    Correctness is pinned against that enumeration at sizes where brute force is
    still cheap; the DP is O(N**2) and returns the same numbers.
    """

    @pytest.mark.parametrize(
        "p",
        [
            [0.5, 0.5, 0.5],
            [0.1, 0.9, 0.4, 0.6],
            [0.05, 0.05, 0.05, 0.05, 0.05],
            [0.2, 0.8, 0.5, 0.3, 0.7, 0.9],
        ],
    )
    def test_all_thresholds_agree(self, p):
        for n in range(len(p) + 2):
            assert prob_ntuple(p, n) == pytest.approx(_brute_force(p, n), abs=1e-12)

    def test_tail_is_monotone_non_increasing(self):
        p = [0.3, 0.6, 0.2, 0.8, 0.5]
        tails = [prob_ntuple(p, n) for n in range(len(p) + 1)]
        assert all(a >= b - 1e-15 for a, b in zip(tails, tails[1:]))
