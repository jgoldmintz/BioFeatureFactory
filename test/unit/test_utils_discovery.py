"""
Unit tests for file discovery, batch processing, and helper functions in utility.py.

Run with: pytest test/unit/test_utils_discovery.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    discover_fasta_files,
    discover_mapping_files,
    validate_fasta_content,
    validate_mapping_content,
    split_fasta_into_batches,
    read_fasta,
    ExtractGeneFromFASTA,
    extract_mutation_from_sequence_name,
    _normalize_logs,
    _collect_failures_from_logs,
    _FILTER_LOG_CACHE,
    trim_muts,
    should_skip_mutation,
    load_wt_sequences,
    resolve_output_base,
)


def _write(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)
    return str(path)


# ===========================================================================
# discover_fasta_files
# ===========================================================================


class TestDiscoverFastaFiles:

    def test_finds_fasta(self, tmp_path):
        _write(tmp_path / "BRCA1.fasta", ">seq\nACGT\n")
        result = discover_fasta_files(str(tmp_path))
        assert "BRCA1" in result

    def test_finds_fa_extension(self, tmp_path):
        _write(tmp_path / "TP53.fa", ">seq\nACGT\n")
        result = discover_fasta_files(str(tmp_path))
        assert "TP53" in result

    def test_ignores_invalid_fasta(self, tmp_path):
        _write(tmp_path / "bad.fasta", "not a fasta\n")
        result = discover_fasta_files(str(tmp_path))
        assert result == {}

    def test_empty_dir(self, tmp_path):
        result = discover_fasta_files(str(tmp_path))
        assert result == {}

    def test_nonexistent_dir(self, tmp_path):
        result = discover_fasta_files(str(tmp_path / "nope"))
        assert result == {}

    def test_none_input(self):
        result = discover_fasta_files(None)
        assert result == {}

    def test_multiple_genes(self, tmp_path):
        _write(tmp_path / "BRCA1.fasta", ">seq\nACGT\n")
        _write(tmp_path / "TP53.fasta", ">seq\nTGCA\n")
        result = discover_fasta_files(str(tmp_path))
        assert len(result) == 2


# ===========================================================================
# discover_mapping_files
# ===========================================================================


class TestDiscoverMappingFiles:

    def test_finds_csv(self, tmp_path):
        _write(tmp_path / "transcript_mapping_BRCA1.csv", "mutant,transcript\nA1G,A2G\n")
        result = discover_mapping_files(str(tmp_path))
        assert "BRCA1" in result

    def test_single_file(self, tmp_path):
        f = _write(tmp_path / "BRCA1.csv", "mutant,transcript\nA1G,A2G\n")
        result = discover_mapping_files(f)
        assert "BRCA1" in result

    def test_nonexistent(self, tmp_path):
        result = discover_mapping_files(str(tmp_path / "nope"))
        assert result == {}

    def test_none_input(self):
        result = discover_mapping_files(None)
        assert result == {}


# ===========================================================================
# split_fasta_into_batches
# ===========================================================================


class TestSplitFastaIntoBatches:

    def test_splits_correctly(self, tmp_path):
        content = "".join(f">seq{i}\nACGT\n" for i in range(5))
        f = _write(tmp_path / "input.fasta", content)
        batches = split_fasta_into_batches(f, batch_size=2, temp_dir=str(tmp_path))
        assert len(batches) == 3  # 2+2+1
        # First batch should have 2 sequences
        b1 = read_fasta(batches[0])
        assert len(b1) == 2

    def test_single_batch_when_small(self, tmp_path):
        content = ">seq1\nACGT\n>seq2\nTGCA\n"
        f = _write(tmp_path / "input.fasta", content)
        batches = split_fasta_into_batches(f, batch_size=10, temp_dir=str(tmp_path))
        assert len(batches) == 1

    def test_empty_fasta(self, tmp_path):
        f = _write(tmp_path / "empty.fasta", "")
        batches = split_fasta_into_batches(f, temp_dir=str(tmp_path))
        assert batches == []


# ===========================================================================
# ExtractGeneFromFASTA
# ===========================================================================


class TestExtractGeneFromFASTA:

    def test_basic_extraction(self, tmp_path):
        f = _write(tmp_path / "test.fasta", ">BRCA1-A123G\nACGT\n")
        assert ExtractGeneFromFASTA(f) == "BRCA1"

    def test_with_count(self, tmp_path):
        f = _write(tmp_path / "test.fasta", ">BRCA1-A123G\nACGT\n>BRCA1-C456T\nTGCA\n")
        gene, count = ExtractGeneFromFASTA(f, count=True)
        assert gene == "BRCA1"
        assert count == 2

    def test_no_separator(self, tmp_path):
        f = _write(tmp_path / "test.fasta", ">BRCA1\nACGT\n")
        assert ExtractGeneFromFASTA(f) == "BRCA1"

    def test_empty_file(self, tmp_path):
        f = _write(tmp_path / "empty.fasta", "")
        assert ExtractGeneFromFASTA(f) is None


# ===========================================================================
# _normalize_logs
# ===========================================================================


class TestNormalizeLogs:

    def test_single_file(self, tmp_path):
        f = _write(tmp_path / "test.log", "content")
        result = _normalize_logs(f)
        assert len(result) == 1

    def test_directory_gathers_logs(self, tmp_path):
        _write(tmp_path / "a.log", "content")
        _write(tmp_path / "b.log", "content")
        result = _normalize_logs(str(tmp_path))
        assert len(result) == 2

    def test_none_returns_empty(self):
        assert _normalize_logs(None) == []

    def test_list_of_paths(self, tmp_path):
        f1 = _write(tmp_path / "a.log", "content")
        f2 = _write(tmp_path / "b.log", "content")
        result = _normalize_logs([f1, f2])
        assert len(result) == 2


# ===========================================================================
# trim_muts
# ===========================================================================


class TestTrimMuts:

    def test_reads_mutations(self, tmp_path):
        f = _write(tmp_path / "muts.csv", "mutant\nA123G\nC456T\n")
        result = trim_muts(f)
        assert "A123G" in result
        assert "C456T" in result

    def test_skips_header(self, tmp_path):
        f = _write(tmp_path / "muts.csv", "mutant\nA123G\n")
        result = trim_muts(f)
        assert "mutant" not in result

    def test_strips_asterisks(self, tmp_path):
        f = _write(tmp_path / "muts.csv", "mutant\nA123G*\n")
        result = trim_muts(f)
        assert "A123G" in result

    def test_filters_with_log(self, tmp_path):
        # Clear cache to avoid stale results
        _FILTER_LOG_CACHE.clear()
        muts = _write(tmp_path / "muts.csv", "mutant\nA123G\nC456T\n")
        log = _write(tmp_path / "val.log", "TESTGENE: mutation A123G expects C\n")
        result = trim_muts(muts, log=log, gene_name="TESTGENE")
        assert "C456T" in result
        assert "A123G" not in result


# ===========================================================================
# load_wt_sequences
# ===========================================================================


class TestLoadWtSequences:

    def test_loads_from_dir(self, tmp_path):
        _write(tmp_path / "BRCA1.fasta", ">transcript\nACGTACGT\n>ORF\nATGAAA\n")
        result = load_wt_sequences(str(tmp_path), wt_header="transcript")
        assert "BRCA1" in result
        assert result["BRCA1"] == "ACGTACGT"

    def test_prefers_specified_header(self, tmp_path):
        _write(tmp_path / "TP53.fasta", ">transcript\nAAAA\n>ORF\nCCCC\n")
        result = load_wt_sequences(str(tmp_path), wt_header="ORF")
        assert result["TP53"] == "CCCC"

    def test_fallback_to_longest(self, tmp_path):
        _write(tmp_path / "GENE.fasta", ">short\nAA\n>long\nACGTACGTACGT\n")
        result = load_wt_sequences(str(tmp_path), wt_header="nonexistent")
        assert "GENE" in result
        assert result["GENE"] == "ACGTACGTACGT"

    def test_single_file(self, tmp_path):
        f = _write(tmp_path / "BRCA1.fasta", ">transcript\nACGT\n")
        result = load_wt_sequences(f, wt_header="transcript")
        assert "BRCA1" in result

    def test_nonexistent_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_wt_sequences(str(tmp_path / "nope"))


# ===========================================================================
# resolve_output_base
# ===========================================================================


class TestResolveOutputBase:

    def test_creates_nested_dir(self, tmp_path):
        result = resolve_output_base(str(tmp_path), "BRCA1.fasta", "NetPhos")
        assert "BRCA1" in result
        assert "NetPhos" in result
        assert Path(result).parent.exists()

    def test_returns_gene_as_basename(self, tmp_path):
        result = resolve_output_base(str(tmp_path), "TP53_mutations.csv", "SpliceAI")
        assert result.endswith("TP53")
