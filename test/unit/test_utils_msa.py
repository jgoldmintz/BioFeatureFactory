"""
Unit tests for MSA parsing and conversion functions in utility.py.

Run with: pytest test/unit/test_utils_msa.py -v
"""

import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import parse_stockholm, stockholm_to_a2m, filter_msa_by_gaps



# helpers


def write_sto(tmp_path, content):
    p = tmp_path / "test.sto"
    p.write_text(textwrap.dedent(content))
    return str(p)



# parse_stockholm


class TestParseStockholm:

    def test_sequences_parsed(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            seq1  ACGT
            seq2  AC-T
            #=GC RF xxxx
            //
        """)
        msa, _ = parse_stockholm(sto)
        assert msa["seq1"] == "ACGT"
        assert msa["seq2"] == "AC-T"

    def test_rf_annotation_captured(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            seq1  ACGT-A
            seq2  AC-TAA
            #=GC RF xxxx.x
            //
        """)
        _, rf = parse_stockholm(sto)
        assert rf == "xxxx.x"

    def test_rf_none_when_absent(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            seq1  ACGT
            //
        """)
        _, rf = parse_stockholm(sto)
        assert rf is None

    def test_multiblock_sequences_concatenated(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            seq1  ACGT
            seq2  AC-T
            #=GC RF xxxx

            seq1  GGCC
            seq2  GG-C
            #=GC RF xxxx
            //
        """)
        msa, rf = parse_stockholm(sto)
        assert msa["seq1"] == "ACGTGGCC"
        assert msa["seq2"] == "AC-TGG-C"
        assert rf == "xxxxxxxx"

    def test_other_gc_annotations_not_in_msa(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            seq1  ACGT
            #=GC SS_cons ....
            #=GC RF xxxx
            //
        """)
        msa, rf = parse_stockholm(sto)
        assert "#=GC" not in msa
        assert "SS_cons" not in msa
        assert rf == "xxxx"

    def test_empty_file_returns_empty(self, tmp_path):
        sto = write_sto(tmp_path, """\
            # STOCKHOLM 1.0
            //
        """)
        msa, rf = parse_stockholm(sto)
        assert msa == {}
        assert rf is None



# stockholm_to_a2m


class TestStockholmToA2m:
    # Alignment:  col 0-3 are match (x), col 4 is insert (.)
    # seq1: A C G T - A
    # seq2: A C - T A A
    # RF:   x x x x . x

    MSA = {"seq1": "ACGT-A", "seq2": "AC-TAA"}
    RF = "xxxx.x"

    def test_match_columns_uppercased(self):
        a2m = stockholm_to_a2m(self.MSA, "seq1", self.RF)
        for i in [0, 1, 2, 3, 5]:
            assert a2m["seq1"][i].isupper() or a2m["seq1"][i] == "-"

    def test_match_column_gap_becomes_dash(self):
        a2m = stockholm_to_a2m(self.MSA, "seq1", self.RF)
        # seq2 col2 = '-' in a match column
        assert a2m["seq2"][2] == "-"

    def test_insert_column_gap_becomes_dot(self):
        a2m = stockholm_to_a2m(self.MSA, "seq1", self.RF)
        # seq1 col4 = '-' in insert column
        assert a2m["seq1"][4] == "."

    def test_insert_column_residue_becomes_lowercase(self):
        a2m = stockholm_to_a2m(self.MSA, "seq1", self.RF)
        # seq2 col4 = 'A' in insert column
        assert a2m["seq2"][4] == "a"

    def test_all_sequences_same_length(self):
        a2m = stockholm_to_a2m(self.MSA, "seq1", self.RF)
        lengths = {len(v) for v in a2m.values()}
        assert len(lengths) == 1

    def test_rf_none_fallback_uses_focus_gaps(self):
        # No RF: match cols = non-gap positions of focus sequence
        # focus="ACGT" has no gaps → all 4 cols are match → all uppercase
        msa = {"focus": "ACGT", "other": "AC-T"}
        a2m = stockholm_to_a2m(msa, "focus", rf_annotation=None)
        assert a2m["focus"] == "ACGT"
        assert a2m["other"] == "AC-T"

    def test_focus_not_found_raises(self):
        with pytest.raises(ValueError, match="not found"):
            stockholm_to_a2m({"seq1": "ACGT"}, "missing", rf_annotation=None)

    def test_focus_partial_match_resolves(self):
        msa = {"ORF|NM_000123": "ACGT"}
        a2m = stockholm_to_a2m(msa, "ORF", rf_annotation=None)
        assert "ORF|NM_000123" in a2m



# filter_msa_by_gaps


class TestFilterMsaByGaps:

    def test_gappy_sequence_removed(self):
        msa = {"good": "ACGT", "gappy": "----"}
        filtered = filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=1.0)
        assert "good" in filtered
        assert "gappy" not in filtered

    def test_gappy_column_removed(self):
        msa = {"seq1": "-CGT", "seq2": "-CGT", "seq3": "-CGT"}
        filtered = filter_msa_by_gaps(msa, max_seq_gaps=1.0, max_col_gaps=0.6)
        for seq in filtered.values():
            assert len(seq) == 3

    def test_empty_msa_returns_empty(self):
        assert filter_msa_by_gaps({}) == {}

    def test_a2m_insert_columns_not_counted_as_gaps(self):
        # lowercase and '.' are insert states; must not count as gaps in match-col accounting
        msa = {
            "seq1": "ACgt",   # match cols: A,C  insert cols: g,t
            "seq2": "AC..",   # match cols: A,C  insert cols: . (insert gap)
        }
        filtered = filter_msa_by_gaps(msa, max_seq_gaps=0.0, max_col_gaps=1.0, a2m_format=True)
        assert "seq1" in filtered
        assert "seq2" in filtered

    def test_a2m_match_column_gaps_trigger_removal(self):
        msa = {
            "seq1": "ACgt",
            "seq2": "--gt",   # 100% gaps in match columns
        }
        filtered = filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=1.0, a2m_format=True)
        assert "seq1" in filtered
        assert "seq2" not in filtered
