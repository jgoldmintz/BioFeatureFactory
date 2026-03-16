"""
Unit tests for build_mutant_sequences_for_gene in utility.py.

Run with: pytest test/unit/test_build_mutants.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import build_mutant_sequences_for_gene


def _write(path, text):
    Path(path).write_text(text)
    return str(path)


# ===========================================================================
# build_mutant_sequences_for_gene
# ===========================================================================


class TestBuildMutantSequencesForGene:
    """Test mutant AA sequence generation from mutation files."""

    # ORF: ATG AAA TGC GAT -> M K C D
    WT_NT = "ATGAAATGCGAT"
    WT_AA = "MKCD"

    def test_single_column_nt_mutation(self, tmp_path):
        # ORF: ATG AAA TGC GAT -> M K C D
        # Positions (1-based): A1 T2 G3 A4 A5 A6 T7 G8 C9 G10 A11 T12
        # Mutation A4T: codon AAA (pos 4,5,6) -> TAA => K->Stop, skipped
        # Mutation G8A: codon TGC (pos 7,8,9) -> TAC => C->Y at AA pos 3
        mut_file = _write(tmp_path / "muts.csv", "mutant\nG8A\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="nt",
        )
        # G8A: codon TGC -> TAC => C->Y at AA pos 3
        assert "TEST-G8A" in result
        assert result["TEST-G8A"][2] == "Y"  # position 2 (0-indexed) changed from C to Y

    def test_single_column_aa_mutation(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv", "mutant\nK2R\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="aa",
        )
        assert "TEST-K2R" in result
        assert result["TEST-K2R"] == "MRCD"

    def test_csv_format_with_headers(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv",
            "mutant,aamutant\nG7T,C3F\n"
        )
        result = build_mutant_sequences_for_gene(
            gene_name="GENE",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="nt",
        )
        assert "GENE-G7T" in result
        assert result["GENE-G7T"][2] == "F"

    def test_skip_mutation_via_failure_map(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv", "mutant\nK2R\n")
        failure_map = {"TEST": {"K2R"}}
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=failure_map,
            input_type="aa",
        )
        assert result == {}

    def test_wt_mismatch_skipped(self, tmp_path):
        # Claim C is at pos 2, but WT has K at pos 2
        mut_file = _write(tmp_path / "muts.csv", "mutant\nC2F\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="aa",
        )
        assert result == {}

    def test_out_of_range_skipped(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv", "mutant\nK999R\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="aa",
        )
        assert result == {}

    def test_missing_file_returns_empty(self):
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence="ATGAAA",
            aa_sequence="MK",
            mapping_file="/nonexistent/path.csv",
            log_path=None,
            failure_map=None,
        )
        assert result == {}

    def test_none_mapping_file(self):
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence="ATGAAA",
            aa_sequence="MK",
            mapping_file=None,
            log_path=None,
            failure_map=None,
        )
        assert result == {}

    def test_multiple_mutations(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv", "mutant\nK2R\nC3F\nD4E\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="aa",
        )
        assert "TEST-K2R" in result
        assert "TEST-C3F" in result
        assert "TEST-D4E" in result
        assert result["TEST-K2R"] == "MRCD"
        assert result["TEST-C3F"] == "MKFD"
        assert result["TEST-D4E"] == "MKCE"

    def test_header_line_skipped_single_column(self, tmp_path):
        mut_file = _write(tmp_path / "muts.csv", "mutant\nK2R\n")
        result = build_mutant_sequences_for_gene(
            gene_name="TEST",
            nt_sequence=self.WT_NT,
            aa_sequence=self.WT_AA,
            mapping_file=mut_file,
            log_path=None,
            failure_map=None,
            input_type="aa",
        )
        # "mutant" header should not be treated as a mutation
        assert len(result) == 1
        assert "TEST-K2R" in result
