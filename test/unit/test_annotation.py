"""
Unit tests for annotation parsing and get_genome_loc in utility.py.

Run with: pytest test/unit/test_annotation.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    get_genome_loc,
    _detect_annotation_format,
    _parse_attributes,
    _normalize_chrom_name,
    _infer_transcript_priority,
    _split_multi_value,
    _prepare_custom_annotation,
)


def _write(path, text):
    Path(path).write_text(text)
    return str(path)


# ===========================================================================
# _detect_annotation_format
# ===========================================================================


class TestDetectAnnotationFormat:

    def test_gtf_by_extension(self, tmp_path):
        f = _write(tmp_path / "anno.gtf", "")
        assert _detect_annotation_format(f) == "gtf"

    def test_gff3_by_extension(self, tmp_path):
        f = _write(tmp_path / "anno.gff3", "")
        assert _detect_annotation_format(f) == "gff3"

    def test_gff_by_extension(self, tmp_path):
        f = _write(tmp_path / "anno.gff", "")
        assert _detect_annotation_format(f) == "gff3"

    def test_gtf_by_content(self, tmp_path):
        line = 'chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG001"; gene_name "TP53";\n'
        f = _write(tmp_path / "anno.txt", line)
        assert _detect_annotation_format(f) == "gtf"

    def test_gff3_by_content(self, tmp_path):
        line = "chr1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-TP53;Name=TP53\n"
        f = _write(tmp_path / "anno.txt", line)
        assert _detect_annotation_format(f) == "gff3"

    def test_custom_format(self, tmp_path):
        line = "TP53\tchr17\t+\t100\t200\t100,150\t130,200\n"
        f = _write(tmp_path / "anno.txt", line)
        assert _detect_annotation_format(f) == "custom"

    def test_nonexistent_file(self, tmp_path):
        assert _detect_annotation_format(str(tmp_path / "nope.txt")) == "custom"


# ===========================================================================
# _parse_attributes
# ===========================================================================


class TestParseAttributes:

    def test_gtf_format(self):
        attr = 'gene_id "ENSG001"; gene_name "TP53";'
        result = _parse_attributes(attr, "gtf")
        assert result["gene_id"] == "ENSG001"
        assert result["gene_name"] == "TP53"

    def test_gff3_format(self):
        attr = "ID=gene-TP53;Name=TP53;gene_id=GeneID:7157"
        result = _parse_attributes(attr, "gff3")
        assert result["ID"] == "gene-TP53"
        assert result["Name"] == "TP53"
        assert result["gene_id"] == "GeneID:7157"

    def test_empty_string(self):
        assert _parse_attributes("", "gtf") == {}
        assert _parse_attributes("", "gff3") == {}

    def test_gtf_no_value(self):
        attr = "some_flag;"
        result = _parse_attributes(attr, "gtf")
        assert result["some_flag"] == ""


# ===========================================================================
# _normalize_chrom_name
# ===========================================================================


class TestNormalizeChromName:

    def test_chr_prefix_stripped(self):
        result = _normalize_chrom_name("chr17", "GRCh38")
        assert result == "17"

    def test_bare_number_unchanged(self):
        result = _normalize_chrom_name("17", "GRCh38")
        assert result == "17"

    def test_refseq_accession_resolved(self):
        result = _normalize_chrom_name("NC_000017.11", "GRCh38")
        assert result == "17"

    def test_none_input(self):
        assert _normalize_chrom_name(None, "GRCh38") is None

    def test_empty_string(self):
        assert _normalize_chrom_name("", "GRCh38") == ""


# ===========================================================================
# _split_multi_value
# ===========================================================================


class TestSplitMultiValue:

    def test_semicolon_separated(self):
        assert _split_multi_value("a;b;c") == ["a", "b", "c"]

    def test_comma_separated(self):
        assert _split_multi_value("a,b,c") == ["a", "b", "c"]

    def test_mixed(self):
        assert _split_multi_value("a;b,c") == ["a", "b", "c"]

    def test_empty(self):
        assert _split_multi_value("") == []

    def test_none(self):
        assert _split_multi_value(None) == []

    def test_whitespace_trimmed(self):
        assert _split_multi_value(" a ; b ") == ["a", "b"]


# ===========================================================================
# _infer_transcript_priority
# ===========================================================================


class TestInferTranscriptPriority:

    def test_mane_select_highest(self):
        attrs = {"tag": "MANE Select"}
        assert _infer_transcript_priority(attrs) == 4

    def test_refseq_select(self):
        attrs = {"tag": "RefSeq Select"}
        assert _infer_transcript_priority(attrs) == 3

    def test_protein_coding_biotype(self):
        attrs = {"transcript_biotype": "protein_coding"}
        assert _infer_transcript_priority(attrs) == 2

    def test_mrna_biotype(self):
        attrs = {"biotype": "mRNA"}
        assert _infer_transcript_priority(attrs) == 1

    def test_no_priority_signals(self):
        attrs = {"some_other": "value"}
        assert _infer_transcript_priority(attrs) == 0

    def test_ccds_tag(self):
        attrs = {"tag": "CCDS"}
        assert _infer_transcript_priority(attrs) == 2

    def test_product_protein(self):
        attrs = {"product": "protein phosphatase"}
        assert _infer_transcript_priority(attrs) == 1


# ===========================================================================
# _prepare_custom_annotation
# ===========================================================================


class TestPrepareCustomAnnotation:

    def test_basic_lookup(self, tmp_path):
        content = "TP53\tchr17\t-\t7571720\t7590868\t7571720,7578176\t7573008,7580603\n"
        f = _write(tmp_path / "anno.txt", content)
        result = _prepare_custom_annotation("TP53", f)
        assert result is not None
        assert result["chrom"] == "17"
        assert result["strand"] == "-"
        assert result["tx_start"] == 7571720
        assert result["tx_end"] == 7590868
        assert len(result["exons"]) == 2

    def test_case_insensitive(self, tmp_path):
        content = "TP53\tchr17\t-\t100\t200\t100\t200\n"
        f = _write(tmp_path / "anno.txt", content)
        result = _prepare_custom_annotation("tp53", f)
        assert result is not None

    def test_gene_not_found(self, tmp_path):
        content = "BRCA1\tchr17\t+\t100\t200\t100\t200\n"
        f = _write(tmp_path / "anno.txt", content)
        result = _prepare_custom_annotation("TP53", f)
        assert result is None

    def test_comments_skipped(self, tmp_path):
        content = "#header\nTP53\tchr17\t-\t100\t200\t100\t200\n"
        f = _write(tmp_path / "anno.txt", content)
        result = _prepare_custom_annotation("TP53", f)
        assert result is not None


# ===========================================================================
# get_genome_loc (integration of the above)
# ===========================================================================


class TestGetGenomeLoc:

    def test_custom_annotation(self, tmp_path):
        content = "BRCA1\tchr17\t-\t43044295\t43170245\t43044295,43090943\t43047642,43091032\n"
        f = _write(tmp_path / "anno.txt", content)
        result = get_genome_loc("BRCA1", f)
        assert result is not None
        assert result["strand"] == "-"
        assert result["chrom"] == "17"
        assert len(result["exons"]) == 2

    def test_gtf_annotation(self, tmp_path):
        lines = [
            'chr17\tENSEMBL\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\texon\t100\t200\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\texon\t300\t500\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "TP53";\n',
        ]
        f = _write(tmp_path / "anno.gtf", "".join(lines))
        result = get_genome_loc("TP53", f, assembly="GRCh38")
        assert result is not None
        assert result["strand"] == "+"
        assert len(result["exons"]) == 2
        assert result["tx_start"] == 100
        assert result["tx_end"] == 500
        assert result["transcript_id"] == "ENST001"

    def test_gene_not_found(self, tmp_path):
        lines = [
            'chr17\tENSEMBL\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG001"; gene_name "TP53";\n',
        ]
        f = _write(tmp_path / "anno.gtf", "".join(lines))
        result = get_genome_loc("BRCA2", f)
        assert result is None

    def test_nonexistent_file(self, tmp_path):
        result = get_genome_loc("TP53", str(tmp_path / "nope.gtf"))
        assert result is None

    def test_specific_transcript_id(self, tmp_path):
        lines = [
            'chr17\tENSEMBL\tgene\t100\t800\t.\t+\t.\tgene_id "ENSG001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\texon\t100\t200\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "TP53";\n',
            'chr17\tENSEMBL\ttranscript\t100\t800\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST002"; gene_name "TP53";\n',
            'chr17\tENSEMBL\texon\t100\t300\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST002"; gene_name "TP53";\n',
            'chr17\tENSEMBL\texon\t600\t800\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST002"; gene_name "TP53";\n',
        ]
        f = _write(tmp_path / "anno.gtf", "".join(lines))
        result = get_genome_loc("TP53", f, transcript_id="ENST002")
        assert result is not None
        assert result["transcript_id"] == "ENST002"
        assert len(result["exons"]) == 2

    def test_gff3_annotation(self, tmp_path):
        lines = [
            "chr1\tRefSeq\tgene\t100\t500\t.\t+\t.\tID=gene-ABCB1;Name=ABCB1\n",
            "chr1\tRefSeq\tmRNA\t100\t500\t.\t+\t.\tID=rna-NM_001;Parent=gene-ABCB1;transcript_id=NM_001\n",
            "chr1\tRefSeq\texon\t100\t250\t.\t+\t.\tParent=rna-NM_001;transcript_id=NM_001\n",
            "chr1\tRefSeq\texon\t350\t500\t.\t+\t.\tParent=rna-NM_001;transcript_id=NM_001\n",
        ]
        f = _write(tmp_path / "anno.gff3", "".join(lines))
        result = get_genome_loc("ABCB1", f, assembly="GRCh38")
        assert result is not None
        assert len(result["exons"]) == 2
        assert result["tx_start"] == 100
        assert result["tx_end"] == 500
