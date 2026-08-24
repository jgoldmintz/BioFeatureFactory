"""
Unit tests for spliceai-parser.py:
parse_spliceai_entries

Run with: pytest test/unit/test_spliceai_parser.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "spliceai" / "bin"))

# spliceai-parser.py has a hyphen so import via importlib
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "spliceai_parser",
    Path(__file__).parent.parent.parent / "biofeaturefactory" / "spliceai" / "bin" / "spliceai-parser.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
parse_spliceai_entries = _mod.parse_spliceai_entries



# parse_spliceai_entries


class TestParseSpliceAIEntries:

    VALID = "SpliceAI=T|BRCA1|0.01|0.00|0.02|0.00|5|-7|23|-42"

    def test_valid_entry_parsed(self):
        results = parse_spliceai_entries(self.VALID)
        assert len(results) == 1
        e = results[0]
        assert e["allele"] == "T"
        assert e["symbol"] == "BRCA1"
        assert e["ds_ag"] == pytest.approx(0.01)
        assert e["ds_al"] == pytest.approx(0.00)
        assert e["ds_dg"] == pytest.approx(0.02)
        assert e["ds_dl"] == pytest.approx(0.00)
        assert e["dp_ag"] == 5
        assert e["dp_al"] == -7
        assert e["dp_dg"] == 23
        assert e["dp_dl"] == -42

    def test_non_spliceai_field_returns_empty(self):
        assert parse_spliceai_entries("DP=50;AF=0.3") == []

    def test_malformed_entry_skipped(self):
        # Only 5 pipe-separated parts instead of 10
        assert parse_spliceai_entries("SpliceAI=T|BRCA1|0.01|0.00|0.02") == []

    def test_multiple_entries_parsed(self):
        info = "SpliceAI=T|BRCA1|0.01|0.00|0.02|0.00|5|-7|23|-42,T|TP53|0.50|0.00|0.00|0.00|1|2|3|4"
        results = parse_spliceai_entries(info)
        assert len(results) == 2
        assert results[0]["symbol"] == "BRCA1"
        assert results[1]["symbol"] == "TP53"

    def test_duplicate_entries_deduplicated(self):
        entry = "T|BRCA1|0.01|0.00|0.02|0.00|5|-7|23|-42"
        info = f"SpliceAI={entry},{entry}"
        results = parse_spliceai_entries(info)
        assert len(results) == 1

    def test_block_labels_assigned(self):
        entries = ",".join([
            "T|GENE1|0.01|0.00|0.02|0.00|5|-7|23|-42",
            "T|GENE2|0.10|0.00|0.00|0.00|1|2|3|4",
            "T|GENE3|0.20|0.00|0.00|0.00|1|2|3|4",
            "T|GENE4|0.30|0.00|0.00|0.00|1|2|3|4",
        ])
        results = parse_spliceai_entries(f"SpliceAI={entries}")
        labels = [r["block_label"] for r in results]
        assert labels == ["A", "B", "C", "D"]

    def test_scores_are_floats(self):
        results = parse_spliceai_entries(self.VALID)
        for key in ("ds_ag", "ds_al", "ds_dg", "ds_dl"):
            assert isinstance(results[0][key], float)

    def test_positions_are_ints(self):
        results = parse_spliceai_entries(self.VALID)
        for key in ("dp_ag", "dp_al", "dp_dg", "dp_dl"):
            assert isinstance(results[0][key], int)

    # --- INFO is a ';'-delimited field list; SpliceAI need not come first -----

    def test_spliceai_field_found_after_other_info_keys(self):
        """A real annotated VCF carries other INFO keys before SpliceAI=.

        Matching only info.startswith('SpliceAI=') discarded every such record
        with no diagnostic, so an externally annotated VCF (accepted by main.nf
        via --input_vcf_file) produced an empty output file.
        """
        info = "AC=1;AN=2;SpliceAI=T|BRCA1|0.01|0.00|0.02|0.00|5|-7|23|-42;DP=30"
        results = parse_spliceai_entries(info)
        assert len(results) == 1
        assert results[0]["symbol"] == "BRCA1"
        assert results[0]["dp_dl"] == -42

    def test_spliceai_first_with_trailing_keys(self):
        info = "SpliceAI=T|BRCA1|0.01|0.00|0.02|0.00|5|-7|23|-42;DP=30"
        results = parse_spliceai_entries(info)
        assert len(results) == 1
        assert results[0]["dp_dl"] == -42

    def test_info_value_helper(self):
        assert _mod.spliceai_info_value("SpliceAI=A|B") == "A|B"
        assert _mod.spliceai_info_value("AC=1;SpliceAI=A|B") == "A|B"
        assert _mod.spliceai_info_value("AC=1;AF=0.3") is None
        assert _mod.spliceai_info_value(".") is None

    # --- refused blocks and unannotated records are COUNTED, never silent ----

    def test_malformed_block_counted_in_stats(self):
        stats = {}
        assert parse_spliceai_entries("SpliceAI=T|BRCA1|0.01|0.00|0.02", stats) == []
        assert stats["malformed_blocks"] == 1

    def test_unparseable_score_counted_in_stats(self):
        stats = {}
        assert parse_spliceai_entries(
            "SpliceAI=T|BRCA1|notafloat|0.00|0.02|0.00|5|-7|23|-42", stats) == []
        assert stats["malformed_blocks"] == 1

    def test_record_without_spliceai_does_not_count_as_malformed(self):
        """No SpliceAI field is a different fact from a corrupt one."""
        stats = {}
        assert parse_spliceai_entries(".", stats) == []
        assert "malformed_blocks" not in stats
