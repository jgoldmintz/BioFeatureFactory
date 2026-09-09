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
"""msa_generation / alphafold3: query length, VCF routing, POSTAR contigs."""

from biofeaturefactory.alphafold3.alphafold3_pipeline import (
    _resolve_gene_vcf,
    parse_vcf_chrom,
)
from biofeaturefactory.alphafold3.bin.rbp_database import POSTAR3Database
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


class TestResolveGeneVcf:
    def test_resolves_variant_mapping_layout_and_ignores_spliceai_vcf(self, tmp_path):
        vcf_dir = tmp_path / "PAM" / "vcf"
        vcf_dir.mkdir(parents=True)
        canonical = vcf_dir / "PAM.vcf"
        canonical.write_text("canonical\n")
        (vcf_dir / "PAM.spliceai.vcf").write_text("annotated\n")

        assert _resolve_gene_vcf(tmp_path, "PAM") == canonical

    def test_preserves_flat_vcf_directory_support(self, tmp_path):
        canonical = tmp_path / "PAM.vcf"
        canonical.write_text("canonical\n")

        assert _resolve_gene_vcf(tmp_path, "PAM") == canonical

    def test_preserves_explicit_vcf_file_support(self, tmp_path):
        explicit = tmp_path / "custom.vcf"
        explicit.write_text("explicit\n")

        assert _resolve_gene_vcf(explicit, "PAM") == explicit


class TestPostarChromosomeAliases:
    def test_refseq_autosome_matches_chr_prefixed_postar_record(self, tmp_path):
        postar = tmp_path / "postar.tsv"
        postar.write_text(
            "chr5\t99\t101\tid\t+\tRBP5\tCLIP\tcell\tacc\t1.0\n"
        )
        database = POSTAR3Database(str(postar), use_tabix=False)

        hits = database.query("NC_000005.10", 99, 101)

        assert [site.rbp_name for site in hits] == ["RBP5"]

    def test_refseq_sex_chromosomes_match_chr_prefixed_records(self, tmp_path):
        postar = tmp_path / "postar.tsv"
        postar.write_text(
            "chrX\t99\t101\tx\t+\tRBPX\tCLIP\tcell\tacc\t1.0\n"
            "chrY\t199\t201\ty\t+\tRBPY\tCLIP\tcell\tacc\t1.0\n"
        )
        database = POSTAR3Database(str(postar), use_tabix=False)

        assert database.query("NC_000023.11", 99, 101)[0].rbp_name == "RBPX"
        assert database.query("NC_000024.10", 199, 201)[0].rbp_name == "RBPY"
