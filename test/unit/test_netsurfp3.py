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
"""NetSurfP3: burial and disorder change classification."""

from biofeaturefactory.NetSurfP3.netsurfp3_pipeline import (
    classify_burial_change,
    classify_disorder_change,
)


class TestClassifyBurialChange:
    def test_no_change_within_a_class(self):
        assert classify_burial_change(0.05, 0.05) == 0
        assert classify_burial_change(0.5, 0.6) == 0

    def test_buried_to_exposed_is_positive(self):
        assert classify_burial_change(0.05, 0.5) == 2

    def test_exposed_to_buried_is_negative(self):
        assert classify_burial_change(0.5, 0.05) == -2

    def test_sign_is_antisymmetric(self):
        assert classify_burial_change(0.05, 0.5) == -classify_burial_change(0.5, 0.05)


class TestClassifyDisorderChange:
    def test_ordered_to_disordered_is_positive(self):
        assert classify_disorder_change(0.1, 0.9) == 2

    def test_no_change_is_zero(self):
        assert classify_disorder_change(0.1, 0.1) == 0

    def test_sign_is_antisymmetric(self):
        assert classify_disorder_change(0.1, 0.9) == -classify_disorder_change(0.9, 0.1)
