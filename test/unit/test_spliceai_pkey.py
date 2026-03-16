"""
Unit tests for generate_pkey_with_mapping in spliceai-parser.py.

Run with: pytest test/unit/test_spliceai_pkey.py -v
"""

import sys
import importlib.util
from pathlib import Path

import pytest

# spliceai-parser.py has a hyphen so import via importlib
_spec = importlib.util.spec_from_file_location(
    "spliceai_parser",
    Path(__file__).parent.parent.parent / "biofeaturefactory" / "spliceai" / "bin" / "spliceai-parser.py"
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
generate_pkey_with_mapping = _mod.generate_pkey_with_mapping


# ===========================================================================
# generate_pkey_with_mapping
# ===========================================================================


class TestGeneratePkeyWithMapping:

    # Typical dual mapping scenario:
    #   chromosome_mapping: {"A123G": "A87504250G"} (ORF mutation -> genomic coordinate)
    #   transcript_mapping: {"A123G": "A123G"}      (ORF mutation -> ORF mutation or mapped form)
    # VCF line has pos=87504250, ref=A, alt=G => chromosome_notation = "A87504250G"
    # Reverse lookup in chromosome_mapping finds "A123G"
    # Then transcript_mapping["A123G"] gives the ORF mutation ID

    CHROM_MAP = {
        "A123G": "A87504250G",
        "C456T": "C87505000T",
        "G789A": "G87506000A",
    }
    TRANS_MAP = {
        "A123G": "A123G",
        "C456T": "C456T",
        "G789A": "G789A",
    }

    def test_basic_mapping(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
        )
        assert result == "BRCA1-A123G"

    def test_second_mutation(self):
        result = generate_pkey_with_mapping(
            pos="87505000", ref="C", alt="T",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
        )
        assert result == "BRCA1-C456T"

    def test_no_match_in_chromosome_map(self):
        result = generate_pkey_with_mapping(
            pos="99999999", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
        )
        assert result is None

    def test_match_in_chrom_but_not_transcript(self):
        partial_trans = {"A123G": "A123G"}  # missing C456T
        result = generate_pkey_with_mapping(
            pos="87505000", ref="C", alt="T",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=partial_trans,
        )
        assert result is None

    def test_empty_gene_context(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
        )
        assert result is None

    def test_empty_mappings(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping={},
            transcript_mapping={},
        )
        assert result is None

    def test_skip_mutations(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
            skip_mutations={"A123G"},
        )
        assert result is None

    def test_skip_mutations_case_insensitive(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
            skip_mutations={"a123g"},
        )
        assert result is None

    def test_skip_does_not_affect_other_mutations(self):
        result = generate_pkey_with_mapping(
            pos="87505000", ref="C", alt="T",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
            skip_mutations={"A123G"},
        )
        assert result == "BRCA1-C456T"

    def test_none_skip_mutations(self):
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=self.TRANS_MAP,
            skip_mutations=None,
        )
        assert result == "BRCA1-A123G"

    def test_different_transcript_mapping_value(self):
        """Transcript mapping can remap to a different mutation ID."""
        trans_map = {"A123G": "A100G"}  # remapped
        result = generate_pkey_with_mapping(
            pos="87504250", ref="A", alt="G",
            gene_context="BRCA1",
            chromosome_mapping=self.CHROM_MAP,
            transcript_mapping=trans_map,
        )
        assert result == "BRCA1-A100G"
