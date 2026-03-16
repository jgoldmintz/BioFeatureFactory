"""
Unit tests for write_tsv and mutation_class in utility.py.

Run with: pytest test/unit/test_write_tsv.py -v
"""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import write_tsv, mutation_class


# ===========================================================================
# write_tsv
# ===========================================================================


class TestWriteTsv:

    def test_basic_write(self, tmp_path):
        rows = [{"a": 1, "b": 2}, {"a": 3, "b": 4}]
        out = tmp_path / "out.tsv"
        write_tsv(rows, out, ["a", "b"])
        lines = out.read_text().strip().split("\n")
        assert lines[0] == "a\tb"
        assert lines[1] == "1\t2"
        assert lines[2] == "3\t4"

    def test_fieldnames_inferred_from_rows(self, tmp_path):
        rows = [{"x": 10, "y": 20}]
        out = tmp_path / "out.tsv"
        write_tsv(rows, out)
        lines = out.read_text().strip().split("\n")
        assert "x" in lines[0]
        assert "y" in lines[0]

    def test_empty_rows_with_fieldnames_writes_header_only(self, tmp_path):
        out = tmp_path / "out.tsv"
        write_tsv([], out, ["col1", "col2"])
        content = out.read_text().strip()
        assert content == "col1\tcol2"

    def test_empty_rows_no_fieldnames_no_file(self, tmp_path):
        out = tmp_path / "out.tsv"
        write_tsv([], out)
        assert not out.exists()

    def test_none_rows_no_fieldnames_no_file(self, tmp_path):
        out = tmp_path / "out.tsv"
        write_tsv(None, out)
        assert not out.exists()

    def test_mkdir_creates_parents(self, tmp_path):
        out = tmp_path / "sub" / "dir" / "out.tsv"
        write_tsv([{"a": 1}], out)
        assert out.exists()

    def test_mkdir_false_no_parent_creation(self, tmp_path):
        out = tmp_path / "nonexistent" / "out.tsv"
        with pytest.raises(FileNotFoundError):
            write_tsv([{"a": 1}], out, mkdir=False)

    def test_extrasaction_ignore(self, tmp_path):
        rows = [{"a": 1, "b": 2, "extra": 99}]
        out = tmp_path / "out.tsv"
        write_tsv(rows, out, ["a", "b"], extrasaction='ignore')
        content = out.read_text()
        assert "extra" not in content
        assert "99" not in content

    def test_extrasaction_raise(self, tmp_path):
        rows = [{"a": 1, "extra": 99}]
        out = tmp_path / "out.tsv"
        with pytest.raises(ValueError):
            write_tsv(rows, out, ["a"], extrasaction='raise')

    def test_tab_delimiter(self, tmp_path):
        rows = [{"col1": "val1", "col2": "val2"}]
        out = tmp_path / "out.tsv"
        write_tsv(rows, out, ["col1", "col2"])
        lines = out.read_text().strip().split("\n")
        assert "\t" in lines[0]
        assert "\t" in lines[1]

    def test_string_path(self, tmp_path):
        rows = [{"a": 1}]
        out = str(tmp_path / "out.tsv")
        write_tsv(rows, out)
        assert Path(out).exists()

    def test_roundtrip_with_csv_reader(self, tmp_path):
        rows = [{"name": "BRCA1", "score": "0.95"}, {"name": "TP53", "score": "0.87"}]
        out = tmp_path / "out.tsv"
        write_tsv(rows, out, ["name", "score"])
        with open(out) as f:
            reader = csv.DictReader(f, delimiter='\t')
            read_rows = list(reader)
        assert len(read_rows) == 2
        assert read_rows[0]["name"] == "BRCA1"
        assert read_rows[1]["score"] == "0.87"

    def test_multiple_writes_overwrite(self, tmp_path):
        out = tmp_path / "out.tsv"
        write_tsv([{"a": 1}], out, ["a"])
        write_tsv([{"a": 2}, {"a": 3}], out, ["a"])
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 rows


# ===========================================================================
# mutation_class
# ===========================================================================


class TestMutationClass:

    def test_synonymous(self):
        assert mutation_class("A", "A") == "SYNONYMOUS"

    def test_missense(self):
        assert mutation_class("A", "V") == "MISSENSE"

    def test_stop_gain(self):
        assert mutation_class("A", "Stop") == "STOP_GAIN"

    def test_stop_loss(self):
        assert mutation_class("Stop", "A") == "STOP_LOSS"

    def test_unknown_wt(self):
        assert mutation_class("X", "A") == "UNKNOWN"

    def test_unknown_mut(self):
        assert mutation_class("A", "X") == "UNKNOWN"

    def test_both_unknown(self):
        assert mutation_class("X", "X") == "UNKNOWN"

    def test_stop_to_stop(self):
        assert mutation_class("Stop", "Stop") == "SYNONYMOUS"
