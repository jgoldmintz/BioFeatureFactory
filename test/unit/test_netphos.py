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
"""netphos: phosphosite gain/loss classification."""

import pytest

from biofeaturefactory.netphos.netphos_pipeline import _classify_netphos_event


class TestClassifyNetphosEvent:
    def test_site_gained(self):
        label, code, delta = _classify_netphos_event(0.2, 0.9, False, True)
        assert label == "gained"
        assert code == 2
        assert delta == pytest.approx(0.7)

    def test_site_lost(self):
        label, code, delta = _classify_netphos_event(0.9, 0.2, True, False)
        assert label == "lost"
        assert code == -2
        assert delta == pytest.approx(-0.7)

    def test_gain_and_loss_codes_are_antisymmetric(self):
        _, gained, _ = _classify_netphos_event(0.2, 0.9, False, True)
        _, lost, _ = _classify_netphos_event(0.9, 0.2, True, False)
        assert gained == -lost

    def test_both_above_threshold_is_stable(self):
        label, code, delta = _classify_netphos_event(0.9, 0.9, True, True)
        assert label == "stable"
        assert code == 0

    def test_neither_above_threshold_is_subthreshold_not_stable(self):
        """A site absent in BOTH is distinct from one present in both.

        Collapsing them onto 'stable' would report an unphosphorylated residue as
        an unchanged phosphosite; the -3 code keeps them separable downstream.
        """
        label, code, _ = _classify_netphos_event(0.2, 0.2, False, False)
        assert label == "subthreshold"
        assert code == -3

    def test_small_delta_above_threshold_stays_stable(self):
        label, _, delta = _classify_netphos_event(0.90, 0.92, True, True)
        assert label == "stable"
        assert delta == pytest.approx(0.02)
