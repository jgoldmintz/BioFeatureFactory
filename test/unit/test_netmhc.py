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
"""netMHC: working-directory stem construction."""

from biofeaturefactory.netMHC.netmhc_pipeline import workdir_stem


class TestWorkdirStem:
    def test_gene_and_mutation_joined(self):
        assert workdir_stem("BRCA1", "A123G") == "BRCA1_A123G"

    def test_long_stem_is_truncated_and_hash_suffixed(self):
        """A filesystem-safe stem, kept unique.

        Non-SNV tokens can be hundreds of characters (a knockout carries the whole
        deleted span), which overruns path limits. Truncating alone would collide
        two different long tokens onto one directory, so a hash of the full token
        is appended.
        """
        stem = workdir_stem("G", "X" * 200)
        assert len(stem) < 200
        assert stem.startswith("G_XXXX")

    def test_two_long_tokens_do_not_collide(self):
        a = workdir_stem("G", "A" * 200)
        b = workdir_stem("G", "A" * 199 + "T")
        assert a != b

    def test_limit_is_respected(self):
        assert len(workdir_stem("G", "X" * 200, limit=40)) <= 40 + 16
