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
"""RNAfold: folding window extraction, JSD, worker autodetect."""

import pytest

from biofeaturefactory.RNAfold.run_viennaRNA_pipeline import (
    _autodetect_workers,
    jsd_unpaired,
    window_around,
)


class TestWindowAround:
    SEQ = "ACGTACGTACGT"  # 12 nt

    def test_centred_window(self):
        win, off = window_around(self.SEQ, 6, 5)
        assert win == "ACGTA"
        assert len(win) == 5
        # off is the 1-based position of the variant WITHIN the window
        assert self.SEQ[6 - 1] == win[off - 1]

    def test_window_clamps_at_sequence_start(self):
        win, off = window_around(self.SEQ, 1, 5)
        assert off == 1
        assert self.SEQ[0] == win[off - 1]

    def test_window_longer_than_sequence_returns_whole_sequence(self):
        win, off = window_around("ACGT", 2, 101)
        assert win == "ACGT"
        assert off == 2

    def test_even_length_is_rejected(self):
        # A centred window needs an odd length or the variant has no centre.
        with pytest.raises(AssertionError):
            window_around(self.SEQ, 6, 4)


class TestJsdUnpaired:
    def test_identical_distributions_are_zero(self):
        assert jsd_unpaired([0.5, 0.5], [0.5, 0.5]) == 0.0

    def test_disjoint_distributions_approach_one_bit(self):
        # JSD is bounded at 1 bit for base-2 logs; disjoint support hits the bound.
        assert jsd_unpaired([1.0, 0.0], [0.0, 1.0]) == pytest.approx(1.0, abs=1e-9)

    def test_is_symmetric(self):
        a, b = [0.9, 0.1], [0.2, 0.8]
        assert jsd_unpaired(a, b) == pytest.approx(jsd_unpaired(b, a))


class TestAutodetectWorkers:
    def test_never_exceeds_task_count(self):
        assert _autodetect_workers(1) == 1

    def test_respects_cap(self):
        assert _autodetect_workers(100, cap=5) == 5

    def test_default_cap_applies(self):
        assert _autodetect_workers(1000) <= 8
