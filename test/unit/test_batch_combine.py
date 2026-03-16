"""
Unit tests for batch output combining functions in utility.py:
combine_batch_outputs, _combine_glycosylation_outputs, _combine_phosphorylation_outputs

Run with: pytest test/unit/test_batch_combine.py -v
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    combine_batch_outputs,
    _combine_glycosylation_outputs,
    _combine_phosphorylation_outputs,
)


# ── Fixtures ─────────────────────────────────────────────────────────────

SAMPLE_NETNGLYC_BATCH = """\
>batch1-netnglyc\t3 amino acids

SeqName                 Position  Potential  N-Glyc result  Comment
======================================================================
GENE1-A100G              50       0.7432     +              N-glyc site
GENE1-A100G              120      0.5100     +              N-glyc site

    GENE1-A100G  MNTPKQLSYFLLLSGALLAAAPQEFQN
"""

SAMPLE_NETNGLYC_BATCH_2 = """\
>batch2-netnglyc\t2 amino acids

SeqName                 Position  Potential  N-Glyc result  Comment
======================================================================
GENE1-B200C              80       0.8100     +              N-glyc site

    GENE1-B200C  MNTPKQLSYFLLLSGALLAAAPQEFQN
"""

SAMPLE_NETPHOS_BATCH = """\
>GENE1\t5 amino acids
#
#  prediction results
#
# Sequence\t\t   # x   Context     Score   Kinase    Answer
# -------------------------------------------------------------------
# GENE1-A100G     Sequence   50  S  MQLSYFL  0.950  PKC       YES
# GENE1-A100G     Sequence  120  T  PQEFQNT  0.830  CKI       YES
"""

SAMPLE_NETPHOS_BATCH_2 = """\
>GENE1\t3 amino acids
#
#  prediction results
#
# Sequence\t\t   # x   Context     Score   Kinase    Answer
# -------------------------------------------------------------------
# GENE1-B200C     Sequence   80  Y  MQLSYFL  0.710  EGFR      YES
"""


@pytest.fixture
def tmp_dir():
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


def _write_batch(tmp_dir, name, content):
    p = tmp_dir / name
    p.write_text(content)
    return str(p)


# ── combine_batch_outputs dispatcher ─────────────────────────────────────

class TestCombineBatchOutputsDispatcher:

    def test_empty_list_returns_false(self, tmp_dir):
        out = str(tmp_dir / "out.out")
        assert combine_batch_outputs([], out) is False

    def test_unsupported_format_returns_false(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", "content")
        out = str(tmp_dir / "out.out")
        assert combine_batch_outputs([f], out, format_type="bogus") is False

    def test_dispatches_netnglyc(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined.out")
        result = combine_batch_outputs([f], out, format_type="netnglyc")
        assert result is True
        assert os.path.exists(out)

    def test_dispatches_netphos(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        out = str(tmp_dir / "combined.out")
        result = combine_batch_outputs([f], out, format_type="netphos")
        assert result is True
        assert os.path.exists(out)


# ── _combine_glycosylation_outputs ───────────────────────────────────────

class TestCombineGlycosylationOutputs:

    def test_single_batch_combined(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined-netnglyc.out")
        result = _combine_glycosylation_outputs([f], out)
        assert result is True
        content = Path(out).read_text()
        assert "SeqName" in content
        assert "GENE1-A100G" in content

    def test_two_batches_merged(self, tmp_dir):
        f1 = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        f2 = _write_batch(tmp_dir, "b2.out", SAMPLE_NETNGLYC_BATCH_2)
        out = str(tmp_dir / "combined-netnglyc.out")
        result = _combine_glycosylation_outputs([f1, f2], out)
        assert result is True
        content = Path(out).read_text()
        # Both genes present
        assert "GENE1-A100G" in content
        assert "GENE1-B200C" in content

    def test_prediction_lines_preserved(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined-netnglyc.out")
        _combine_glycosylation_outputs([f], out)
        content = Path(out).read_text()
        # Prediction data should be in output
        assert "0.7432" in content or "0.5100" in content

    def test_header_has_sequence_count(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined-netnglyc.out")
        _combine_glycosylation_outputs([f], out)
        content = Path(out).read_text()
        assert "amino acids" in content

    def test_sequence_sections_preserved(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined-netnglyc.out")
        _combine_glycosylation_outputs([f], out)
        content = Path(out).read_text()
        assert "MNTPKQLSYFLLLSGALLAAAPQEFQN" in content

    def test_bad_batch_file_skipped(self, tmp_dir):
        good = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        out = str(tmp_dir / "combined-netnglyc.out")
        result = _combine_glycosylation_outputs([good, "/nonexistent/file.out"], out)
        assert result is True
        content = Path(out).read_text()
        assert "GENE1-A100G" in content

    def test_with_original_fasta(self, tmp_dir):
        """When original_fasta_file is provided and exists, uses it for count."""
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETNGLYC_BATCH)
        fasta = tmp_dir / "gene.fasta"
        fasta.write_text(">GENE1-wt\nMNTPKQLSYFLLLSGALLAAAP\n>GENE1-A100G\nMNTPKQLSYFLLLSGALLAAAP\n")
        out = str(tmp_dir / "combined-netnglyc.out")
        result = _combine_glycosylation_outputs([f], out, original_fasta_file=str(fasta))
        assert result is True

    def test_fallback_name_count(self, tmp_dir):
        """Batch with Name: lines uses name count fallback."""
        content = "Name: seq1\nName: seq2\nsome prediction data\n"
        f = _write_batch(tmp_dir, "b1.out", content)
        out = str(tmp_dir / "combined-netnglyc.out")
        result = _combine_glycosylation_outputs([f], out)
        assert result is True


# ── _combine_phosphorylation_outputs ─────────────────────────────────────

class TestCombinePhosphorylationOutputs:

    def test_single_batch_combined(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        out = str(tmp_dir / "combined-netphos.out")
        result = _combine_phosphorylation_outputs([f], out)
        assert result is True
        content = Path(out).read_text()
        assert "GENE1-A100G" in content

    def test_two_batches_merged(self, tmp_dir):
        f1 = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        f2 = _write_batch(tmp_dir, "b2.out", SAMPLE_NETPHOS_BATCH_2)
        out = str(tmp_dir / "combined-netphos.out")
        result = _combine_phosphorylation_outputs([f1, f2], out)
        assert result is True
        content = Path(out).read_text()
        assert "GENE1-A100G" in content
        assert "GENE1-B200C" in content

    def test_prediction_lines_extracted(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        out = str(tmp_dir / "combined-netphos.out")
        _combine_phosphorylation_outputs([f], out)
        content = Path(out).read_text()
        # Prediction lines start with "# " and have 7+ fields
        lines = [l for l in content.split("\n") if "PKC" in l or "CKI" in l]
        assert len(lines) >= 1

    def test_header_written(self, tmp_dir):
        f = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        out = str(tmp_dir / "combined-netphos.out")
        _combine_phosphorylation_outputs([f], out)
        content = Path(out).read_text()
        assert "prediction results" in content
        assert "Kinase" in content

    def test_sequence_count_summed(self, tmp_dir):
        f1 = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        f2 = _write_batch(tmp_dir, "b2.out", SAMPLE_NETPHOS_BATCH_2)
        out = str(tmp_dir / "combined-netphos.out")
        _combine_phosphorylation_outputs([f1, f2], out)
        content = Path(out).read_text()
        # 5 + 3 = 8
        assert "8 amino acids" in content

    def test_bad_batch_file_skipped(self, tmp_dir):
        good = _write_batch(tmp_dir, "b1.out", SAMPLE_NETPHOS_BATCH)
        out = str(tmp_dir / "combined-netphos.out")
        result = _combine_phosphorylation_outputs([good, "/nonexistent/file.out"], out)
        assert result is True

    def test_header_fallback_no_seq_count(self, tmp_dir):
        """Batch without standard header line → fallback count."""
        content = "# GENE1-X     Sequence   50  S  MQLSYFL  0.950  PKC       YES\n"
        f = _write_batch(tmp_dir, "b1.out", content)
        out = str(tmp_dir / "combined-netphos.out")
        result = _combine_phosphorylation_outputs([f], out)
        assert result is True
