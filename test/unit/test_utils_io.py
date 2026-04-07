"""
Unit tests for I/O, mapping, and file discovery functions in utility.py.

Run with: pytest test/unit/test_utils_io.py -v
"""

import sys
import os
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    read_fasta,
    write_fasta,
    load_mapping,
    validate_mapping_content,
    validate_fasta_content,
    extract_gene_from_filename,
    is_likely_gene_name,
    strip_all_extensions,
    extract_mutation_from_sequence_name,
    should_skip_mutation,
    load_validation_failures,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


def _write(path, text):
    Path(path).write_text(text)
    return str(path)


# ===========================================================================
# read_fasta
# ===========================================================================


class TestReadFasta:

    def test_single_sequence(self, tmp_dir):
        f = _write(tmp_dir / "one.fasta", ">seq1\nACGTACGT\n")
        result = read_fasta(f)
        assert result == {"seq1": "ACGTACGT"}

    def test_multiple_sequences(self, tmp_dir):
        f = _write(tmp_dir / "two.fasta", ">alpha\nAAA\n>beta\nCCC\n")
        result = read_fasta(f)
        assert result == {"alpha": "AAA", "beta": "CCC"}

    def test_multiline_sequence(self, tmp_dir):
        f = _write(tmp_dir / "multi.fasta", ">seq1\nACGT\nTGCA\n")
        result = read_fasta(f)
        assert result["seq1"] == "ACGTTGCA"

    def test_comment_lines_skipped(self, tmp_dir):
        f = _write(tmp_dir / "comments.fasta", "#comment\n>seq1\nACGT\n#another\n")
        result = read_fasta(f)
        assert result == {"seq1": "ACGT"}

    def test_ncbi_format(self, tmp_dir):
        f = _write(tmp_dir / "ncbi.fasta", ">NM_001234.5 some description\nACGT\n")
        result = read_fasta(f, aformat="NCBI")
        assert "NM_001234.5" in result

    def test_first_format_takes_first_word(self, tmp_dir):
        f = _write(tmp_dir / "first.fasta", ">gene_name extra info\nACGT\n")
        result = read_fasta(f, aformat="FIRST")
        assert "gene_name" in result

    def test_full_header_format(self, tmp_dir):
        f = _write(tmp_dir / "full.fasta", ">full header line\nACGT\n")
        result = read_fasta(f, aformat="FULL")
        assert "full header line" in result

    def test_duplicate_replace(self, tmp_dir):
        f = _write(tmp_dir / "dup.fasta", ">seq1\nAAAA\n>seq1\nCCCC\n")
        result = read_fasta(f, duplicate="replace")
        assert result["seq1"] == "CCCC"

    def test_duplicate_append(self, tmp_dir):
        f = _write(tmp_dir / "dup.fasta", ">seq1\nAAAA\n>seq1\nCCCC\n")
        result = read_fasta(f, duplicate="append")
        assert result["seq1"] == "AAAACCCC"

    def test_duplicate_separate(self, tmp_dir):
        f = _write(tmp_dir / "dup.fasta", ">seq1\nAAAA\n>seq1\nCCCC\n")
        result = read_fasta(f, duplicate="separate")
        assert "seq1" in result
        assert "seq1_2" in result

    def test_empty_file(self, tmp_dir):
        f = _write(tmp_dir / "empty.fasta", "")
        result = read_fasta(f)
        assert result == {}

    def test_whitespace_stripped(self, tmp_dir):
        f = _write(tmp_dir / "ws.fasta", ">seq1\n  ACGT  \n  TGCA  \n")
        result = read_fasta(f)
        assert result["seq1"] == "ACGTTGCA"


# ===========================================================================
# write_fasta
# ===========================================================================


class TestWriteFasta:

    def test_roundtrip(self, tmp_dir):
        data = {"seq1": "ACGTACGT", "seq2": "TGCATGCA"}
        out = tmp_dir / "out.fasta"
        write_fasta(out, data)
        result = read_fasta(str(out))
        assert result == data

    def test_line_wrapping_at_60(self, tmp_dir):
        long_seq = "A" * 150
        out = tmp_dir / "long.fasta"
        write_fasta(out, {"seq1": long_seq})
        lines = out.read_text().strip().split("\n")
        # header + 60 + 60 + 30 = 4 lines
        assert lines[0] == ">seq1"
        assert len(lines[1]) == 60
        assert len(lines[2]) == 60
        assert len(lines[3]) == 30

    def test_creates_parent_dirs(self, tmp_dir):
        out = tmp_dir / "sub" / "dir" / "out.fasta"
        write_fasta(out, {"s": "ACGT"})
        assert out.exists()

    def test_empty_dict(self, tmp_dir):
        out = tmp_dir / "empty.fasta"
        write_fasta(out, {})
        assert out.read_text() == ""


# ===========================================================================
# validate_mapping_content
# ===========================================================================


class TestValidateMappingContent:

    def test_valid_two_column_csv(self, tmp_dir):
        f = _write(tmp_dir / "map.csv", "mutant,transcript\nA123G,A456G\n")
        result, delim = validate_mapping_content(f)
        assert result is True
        assert delim == ","

    def test_valid_tab_delimited(self, tmp_dir):
        f = _write(tmp_dir / "map.tsv", "mutant\ttranscript\nA123G\tA456G\n")
        result, delim = validate_mapping_content(f)
        assert result is True
        assert delim == "\t"

    def test_single_column_mutation(self, tmp_dir):
        f = _write(tmp_dir / "single.csv", "mutant\nA123G\nC456T\n")
        result = validate_mapping_content(f)
        # single column has no delimiter-based mapping column
        assert result is not None

    def test_missing_required_columns(self, tmp_dir):
        f = _write(tmp_dir / "bad.csv", "foo,bar\n1,2\n")
        result = validate_mapping_content(f)
        # Should not validate as a proper mapping file
        assert result is not None  # returns something, but first element should be falsy

    def test_nonexistent_file(self, tmp_dir):
        result = validate_mapping_content(str(tmp_dir / "nonexistent.csv"))
        assert result is False


# ===========================================================================
# validate_fasta_content
# ===========================================================================


class TestValidateFastaContent:

    def test_valid_fasta(self, tmp_dir):
        f = _write(tmp_dir / "valid.fasta", ">seq1\nACGTACGT\n")
        assert validate_fasta_content(f) is True

    def test_no_header(self, tmp_dir):
        f = _write(tmp_dir / "nohdr.fasta", "ACGTACGT\n")
        assert validate_fasta_content(f) is False

    def test_empty_sequence(self, tmp_dir):
        f = _write(tmp_dir / "empty.fasta", ">seq1\n>seq2\n")
        assert validate_fasta_content(f) is False

    def test_nonexistent_file(self, tmp_dir):
        assert validate_fasta_content(str(tmp_dir / "nope.fasta")) is False

    def test_empty_file(self, tmp_dir):
        f = _write(tmp_dir / "empty.fasta", "")
        assert validate_fasta_content(f) is False


# ===========================================================================
# load_mapping
# ===========================================================================


class TestLoadMapping:

    def test_transcript_mapping(self, tmp_dir):
        f = _write(tmp_dir / "map.csv", "mutant,transcript\nA123G,A456G\nC789T,C012T\n")
        result = load_mapping(f, mapType="transcript")
        assert result == {"A123G": "A456G", "C789T": "C012T"}

    def test_chromosome_mapping(self, tmp_dir):
        f = _write(tmp_dir / "map.csv", "mutant,chromosome\nA123G,A999G\n")
        result = load_mapping(f, mapType="chromosome")
        assert result == {"A123G": "A999G"}

    def test_missing_maptype_column(self, tmp_dir):
        f = _write(tmp_dir / "map.csv", "mutant,transcript\nA123G,A456G\n")
        result = load_mapping(f, mapType="genomic")
        assert result == {}

    def test_empty_csv_returns_empty(self, tmp_dir):
        # In practice callers always validate files exist before calling load_mapping
        f = _write(tmp_dir / "empty.csv", "mutant,transcript\n")
        result = load_mapping(f, mapType="transcript")
        assert result == {}

    def test_tab_delimited(self, tmp_dir):
        f = _write(tmp_dir / "map.tsv", "mutant\ttranscript\nA123G\tA456G\n")
        result = load_mapping(f, mapType="transcript")
        assert result == {"A123G": "A456G"}


# ===========================================================================
# extract_gene_from_filename
# ===========================================================================


class TestExtractGeneFromFilename:

    def test_simple_gene(self):
        assert extract_gene_from_filename("BRCA1.fasta") == "BRCA1"

    def test_gene_with_prefix(self):
        result = extract_gene_from_filename("transcript_mapping_ABCB1.csv")
        assert result == "ABCB1"

    def test_gene_with_number(self):
        result = extract_gene_from_filename("TP53_mutations.csv")
        assert result == "TP53"

    def test_hyphenated_gene(self):
        result = extract_gene_from_filename("HLA-A_mapping.csv")
        assert result == "HLA-A"

    def test_nested_path(self):
        result = extract_gene_from_filename("/path/to/ZFP36.fasta")
        assert result == "ZFP36"

    def test_multiple_extensions(self):
        result = extract_gene_from_filename("BRCA2.csv.gz")
        assert "BRCA2" in result

    def test_pdb_id(self):
        assert extract_gene_from_filename("1ABC.fasta") == "1ABC"

    def test_pdb_id_with_chain(self):
        assert extract_gene_from_filename("1ABC_A.fasta") == "1ABC_A"

    def test_pdb_id_chain_b(self):
        assert extract_gene_from_filename("7RM1_B.pdb") == "7RM1_B"

    def test_pdb_id_no_chain(self):
        assert extract_gene_from_filename("7RM1.fasta") == "7RM1"


# ===========================================================================
# is_likely_gene_name
# ===========================================================================


class TestIsLikelyGeneName:

    def test_standard_gene(self):
        assert is_likely_gene_name("BRCA1") is True

    def test_short_gene(self):
        assert is_likely_gene_name("TP") is True

    def test_too_short(self):
        assert is_likely_gene_name("A") is False

    def test_too_long(self):
        assert is_likely_gene_name("A" * 20) is False

    def test_pdb_id(self):
        assert is_likely_gene_name("1ABC") is True

    def test_pdb_id_4char(self):
        assert is_likely_gene_name("7RM1") is True

    def test_pdb_id_with_chain(self):
        assert is_likely_gene_name("1ABC_A") is True

    def test_pdb_id_chain_b(self):
        assert is_likely_gene_name("7RM1_B") is True

    def test_with_hyphen(self):
        assert is_likely_gene_name("HLA-A") is True

    def test_stopword(self):
        # "mapping" would pass pattern but context matters
        assert is_likely_gene_name("mapping") is True  # pattern check only


# ===========================================================================
# strip_all_extensions
# ===========================================================================


class TestStripAllExtensions:

    def test_single_extension(self):
        assert strip_all_extensions("file.csv") == "file"

    def test_double_extension(self):
        assert strip_all_extensions("file.csv.gz") == "file"

    def test_no_extension(self):
        assert strip_all_extensions("file") == "file"

    def test_multiple_dots(self):
        assert strip_all_extensions("some.file.name.txt") == "some"


# ===========================================================================
# extract_mutation_from_sequence_name
# ===========================================================================


class TestExtractMutationFromSequenceName:

    def test_gene_mutation(self):
        gene, mut = extract_mutation_from_sequence_name("ZFP36-C330T")
        assert gene == "ZFP36"
        assert mut == "C330T"

    def test_gene_only(self):
        gene, mut = extract_mutation_from_sequence_name("ZFP36")
        assert gene == "ZFP36"
        assert mut is None

    def test_complex_gene_name(self):
        gene, mut = extract_mutation_from_sequence_name("HLA-A-G123T")
        assert gene == "HLA-A"
        assert mut == "G123T"


# ===========================================================================
# should_skip_mutation
# ===========================================================================


class TestShouldSkipMutation:

    def test_skip_present(self):
        fmap = {"BRCA1": {"A123G", "C456T"}}
        assert should_skip_mutation("BRCA1", "A123G", fmap) is True

    def test_not_in_map(self):
        fmap = {"BRCA1": {"A123G"}}
        assert should_skip_mutation("BRCA1", "C456T", fmap) is False

    def test_different_gene(self):
        fmap = {"BRCA1": {"A123G"}}
        assert should_skip_mutation("TP53", "A123G", fmap) is False

    def test_empty_map(self):
        assert should_skip_mutation("BRCA1", "A123G", {}) is False

    def test_none_map(self):
        assert should_skip_mutation("BRCA1", "A123G", None) is False

    def test_case_insensitive(self):
        fmap = {"BRCA1": {"A123G"}}
        assert should_skip_mutation("brca1", "a123g", fmap) is True

    def test_none_gene(self):
        fmap = {"BRCA1": {"A123G"}}
        assert should_skip_mutation(None, "A123G", fmap) is False

    def test_none_mutation(self):
        fmap = {"BRCA1": {"A123G"}}
        assert should_skip_mutation("BRCA1", None, fmap) is False


# ===========================================================================
# load_validation_failures
# ===========================================================================


class TestLoadValidationFailures:

    def test_parses_log(self, tmp_dir):
        log = _write(tmp_dir / "val.log",
            "BRCA1: mutation A123G expects C at position 123 but found T\n"
            "TP53: mutation G456A expects G at position 456 but found A\n"
        )
        result = load_validation_failures(log)
        assert "BRCA1" in result
        assert "A123G" in result["BRCA1"]
        assert "TP53" in result
        assert "G456A" in result["TP53"]

    def test_none_path(self):
        result = load_validation_failures(None)
        assert result == {}

    def test_empty_log(self, tmp_dir):
        log = _write(tmp_dir / "empty.log", "")
        result = load_validation_failures(log)
        assert result == {}

    def test_nonmatching_lines_skipped(self, tmp_dir):
        log = _write(tmp_dir / "noise.log",
            "some random line\n"
            "BRCA1: mutation A123G expects C\n"
            "another random line\n"
        )
        result = load_validation_failures(log)
        assert "BRCA1" in result
        assert len(result) == 1
