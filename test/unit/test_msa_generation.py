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
"""msa_generation / alphafold3: query length, focus id, VCF contig."""

from biofeaturefactory.alphafold3.alphafold3_pipeline import parse_vcf_chrom
from biofeaturefactory.core.msa_generation_pipeline import get_focus_id, get_query_length


def _write(p, text):
    p.write_text(text)
    return str(p)


class TestQueryFasta:
    FASTA = ">ORF\nATGAAATAG\n>transcript\nAAATGAAATAGCC\n"

    def test_query_length_is_the_first_record(self):
        import tempfile, pathlib
        p = pathlib.Path(tempfile.mkdtemp()) / "G.fasta"
        assert get_query_length(_write(p, self.FASTA)) == 9

    def test_focus_id_is_the_first_header(self):
        import tempfile, pathlib
        p = pathlib.Path(tempfile.mkdtemp()) / "G.fasta"
        assert get_focus_id(_write(p, self.FASTA)) == "ORF"


class TestParseVcfChrom:
    def test_reads_contig_from_first_data_row(self, tmp_path):
        v = tmp_path / "x.vcf"
        v.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "NC_000005.10\t100\t.\tA\tG\t.\t.\t.\n"
        )
        assert parse_vcf_chrom(str(v)) == "NC_000005.10"

    def test_header_only_vcf_returns_none(self, tmp_path):
        """None, not "" -- an absent contig must not be read as a real one."""
        v = tmp_path / "e.vcf"
        v.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\n")
        assert parse_vcf_chrom(str(v)) is None
