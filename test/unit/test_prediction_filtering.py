"""
Unit tests for prediction filtering functions in utility.py:
process_single_mutation_for_sequence, parse_predictions_with_mutation_filtering

Run with: pytest test/unit/test_prediction_filtering.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "utils"))
from utility import (
    process_single_mutation_for_sequence,
    parse_predictions_with_mutation_filtering,
    extract_mutation_from_sequence_name,
    should_skip_mutation,
    get_mutation_data_bioAccurate,
)


# ── Helper factories ──────────────────────────────────────────────────────

def _netphos_pred(gene, pos, aa, score=0.9, answer="YES", kinase="PKC"):
    return {
        "Gene": gene,
        "pos": pos,
        "amino_acid": aa,
        "score": score,
        "answer": answer,
        "kinase": kinase,
    }


def _netnglyc_pred(seq_name, position, sequon, potential=0.7):
    return {
        "seq_name": seq_name,
        "position": position,
        "sequon": sequon,
        "potential": potential,
    }


# ── extract_mutation_from_sequence_name ──────────────────────────────────

class TestExtractMutationFromSequenceName:

    def test_gene_and_mutation(self):
        gene, mut = extract_mutation_from_sequence_name("ZFP36-C330T")
        assert gene == "ZFP36"
        assert mut == "C330T"

    def test_no_mutation(self):
        gene, mut = extract_mutation_from_sequence_name("ZFP36")
        assert gene == "ZFP36"
        assert mut is None

    def test_multi_hyphen_takes_last(self):
        gene, mut = extract_mutation_from_sequence_name("HLA-B-Y110F")
        assert gene == "HLA-B"
        assert mut == "Y110F"


# ── should_skip_mutation ─────────────────────────────────────────────────

class TestShouldSkipMutation:

    def test_skip_when_in_failure_map(self):
        fm = {"GENE1": {"C330T"}}
        assert should_skip_mutation("gene1", "c330t", fm) is True

    def test_no_skip_when_absent(self):
        fm = {"GENE1": {"C330T"}}
        assert should_skip_mutation("gene1", "A100G", fm) is False

    def test_no_skip_empty_map(self):
        assert should_skip_mutation("g", "m", {}) is False

    def test_no_skip_none_values(self):
        assert should_skip_mutation(None, None, {"X": set()}) is False


# ── get_mutation_data_bioAccurate ────────────────────────────────────────

class TestGetMutationDataBioAccurate:

    def test_normal_mutation(self):
        pos, (orig, mut) = get_mutation_data_bioAccurate("Y110F")
        assert pos == 110
        assert orig == "Y"
        assert mut == "F"

    def test_stop_codon_returns_none(self):
        pos, _ = get_mutation_data_bioAccurate("Y110Stop")
        assert pos is None

    def test_position_one(self):
        pos, (orig, mut) = get_mutation_data_bioAccurate("M1V")
        assert pos == 1


# ── process_single_mutation_for_sequence (netphos) ───────────────────────

class TestProcessSingleMutationNetphos:

    MAPPING = {"C330T": "Y110F"}

    def test_matching_prediction_returned(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F")]
        results = process_single_mutation_for_sequence(
            "ZFP36-C330T", preds, self.MAPPING, tool_type="netphos"
        )
        assert len(results) == 1
        assert results[0]["pkey"] == "ZFP36-C330T"
        assert results[0]["Gene"] == "ZFP36"

    def test_wrong_position_filtered(self):
        preds = [_netphos_pred("ZFP36-C330T", 50, "F")]
        results = process_single_mutation_for_sequence(
            "ZFP36-C330T", preds, self.MAPPING, tool_type="netphos"
        )
        assert results == []

    def test_wrong_amino_acid_filtered(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "S")]
        results = process_single_mutation_for_sequence(
            "ZFP36-C330T", preds, self.MAPPING, tool_type="netphos"
        )
        assert results == []

    def test_no_mutation_in_name_returns_empty(self):
        preds = [_netphos_pred("ZFP36", 110, "F")]
        results = process_single_mutation_for_sequence(
            "ZFP36", preds, self.MAPPING, tool_type="netphos"
        )
        assert results == []

    def test_mutation_not_in_mapping_returns_empty(self):
        preds = [_netphos_pred("ZFP36-A100G", 50, "V")]
        results = process_single_mutation_for_sequence(
            "ZFP36-A100G", preds, {"C330T": "Y110F"}, tool_type="netphos"
        )
        assert results == []

    def test_is_mutant_false_raises(self):
        with pytest.raises(ValueError, match="should only be used for mutant"):
            process_single_mutation_for_sequence(
                "ZFP36-C330T", [], self.MAPPING, is_mutant=False
            )

    def test_unsupported_tool_type_raises(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F")]
        with pytest.raises(ValueError, match="Unsupported tool_type"):
            process_single_mutation_for_sequence(
                "ZFP36-C330T", preds, self.MAPPING, tool_type="bogus"
            )

    def test_failure_map_skips_mutation(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F")]
        fm = {"ZFP36": {"C330T"}}
        results = process_single_mutation_for_sequence(
            "ZFP36-C330T", preds, self.MAPPING, tool_type="netphos", failure_map=fm
        )
        assert results == []


# ── process_single_mutation_for_sequence (netnglyc) ──────────────────────

class TestProcessSingleMutationNetnglyc:

    # A300G → N100S means target_aa = 'S' (the mutant aa)
    # netnglyc pred_aa = sequon[0], so sequon must start with 'S' to match
    MAPPING = {"A300G": "N100S"}

    def test_matching_prediction_returned(self):
        # sequon starts with 'S' (target_aa), position matches
        preds = [_netnglyc_pred("GENE1-A300G", 100, "SSTV")]
        results = process_single_mutation_for_sequence(
            "GENE1-A300G", preds, self.MAPPING, tool_type="netnglyc"
        )
        assert len(results) == 1
        assert results[0]["pkey"] == "GENE1-A300G"
        # Field renaming: position → pos, sequon → Sequon
        assert "pos" in results[0]
        assert "Sequon" in results[0]
        assert "seq_name" not in results[0]

    def test_wrong_position_filtered(self):
        preds = [_netnglyc_pred("GENE1-A300G", 50, "SSTV")]
        results = process_single_mutation_for_sequence(
            "GENE1-A300G", preds, self.MAPPING, tool_type="netnglyc"
        )
        assert results == []

    def test_wrong_sequon_aa_filtered(self):
        """sequon[0] = 'N' doesn't match target_aa = 'S'."""
        preds = [_netnglyc_pred("GENE1-A300G", 100, "NSTV")]
        results = process_single_mutation_for_sequence(
            "GENE1-A300G", preds, self.MAPPING, tool_type="netnglyc"
        )
        assert results == []

    def test_empty_sequon_filtered(self):
        preds = [_netnglyc_pred("GENE1-A300G", 100, "")]
        results = process_single_mutation_for_sequence(
            "GENE1-A300G", preds, self.MAPPING, tool_type="netnglyc"
        )
        assert results == []


# ── parse_predictions_with_mutation_filtering ────────────────────────────

class TestParsePredictionsMutationFiltering:

    MAPPING = {"C330T": "Y110F"}

    def test_basic_mutant_filtering(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.95)]
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, tool_type="netphos"
        )
        assert len(results) == 1

    def test_threshold_filters_low_score(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.3)]
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, threshold=0.5, tool_type="netphos"
        )
        assert results == []

    def test_threshold_passes_high_score(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.8)]
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, threshold=0.5, tool_type="netphos"
        )
        assert len(results) == 1

    def test_yes_only_filters_no_answer(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.9, answer=".")]
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, yes_only=True, tool_type="netphos"
        )
        assert results == []

    def test_yes_only_passes_yes_answer(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.9, answer="YES")]
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, yes_only=True, tool_type="netphos"
        )
        assert len(results) == 1

    def test_wildtype_raises_not_implemented(self):
        with pytest.raises(NotImplementedError):
            parse_predictions_with_mutation_filtering(
                [], self.MAPPING, is_mutant=False, tool_type="netphos"
            )

    def test_netnglyc_threshold_filtering(self):
        preds = [_netnglyc_pred("GENE1-C330T", 110, "FXXX", potential=0.4)]
        mapping = {"C330T": "Y110F"}
        results = parse_predictions_with_mutation_filtering(
            preds, mapping, is_mutant=True, threshold=0.5, tool_type="netnglyc"
        )
        assert results == []

    def test_multiple_sequences_grouped(self):
        preds = [
            _netphos_pred("ZFP36-C330T", 110, "F", score=0.9),
            _netphos_pred("ZFP36-C330T", 50, "S", score=0.9),  # wrong pos
            _netphos_pred("BRCA1-A100G", 33, "V", score=0.9),
        ]
        mapping = {"C330T": "Y110F", "A100G": "V33A"}
        results = parse_predictions_with_mutation_filtering(
            preds, mapping, is_mutant=True, tool_type="netphos"
        )
        # Only ZFP36 pos=110 matches; BRCA1 target_aa='A' but pred_aa='V' → no match
        assert len(results) == 1

    def test_failure_map_filters_in_parse(self):
        preds = [_netphos_pred("ZFP36-C330T", 110, "F", score=0.9)]
        fm = {"ZFP36": {"C330T"}}
        results = parse_predictions_with_mutation_filtering(
            preds, self.MAPPING, is_mutant=True, tool_type="netphos", failure_map=fm
        )
        assert results == []

    def test_empty_predictions_returns_empty(self):
        results = parse_predictions_with_mutation_filtering(
            [], self.MAPPING, is_mutant=True, tool_type="netphos"
        )
        assert results == []
