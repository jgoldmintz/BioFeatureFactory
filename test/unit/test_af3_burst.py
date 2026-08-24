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
"""AlphaFold3 burst: failure flag encoding."""

from biofeaturefactory.alphafold3.burst import _failed_flag


class TestFailedFlag:
    def test_encodes_the_reason(self):
        assert _failed_flag("boom") == "FAILED:boom"

    def test_prefix_is_machine_greppable(self):
        # Ingest distinguishes a failed shard from a pending one by this prefix,
        # so it must stay a fixed literal rather than free text.
        assert _failed_flag("anything").startswith("FAILED:")
