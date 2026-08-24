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
"""genesplicer: confidence-to-weight mapping."""

from biofeaturefactory.genesplicer.genesplicer_ensemble import _confidence_to_weight


class TestConfidenceToWeight:
    def test_known_levels_are_ordered(self):
        hi = _confidence_to_weight("High")
        med = _confidence_to_weight("Medium")
        lo = _confidence_to_weight("Low")
        assert hi == 1.0
        assert med == 0.75
        assert lo == 0.5
        assert hi > med > lo

    def test_case_insensitive(self):
        assert _confidence_to_weight("high") == _confidence_to_weight("High")

    def test_unknown_label_falls_back_to_medium(self):
        """Documents the fallback: an unrecognised confidence is NOT dropped.

        GeneSplicer emits its confidence as free text, so an unexpected label
        would otherwise silently weight a real site at zero.
        """
        assert _confidence_to_weight("bogus") == _confidence_to_weight("Medium")
