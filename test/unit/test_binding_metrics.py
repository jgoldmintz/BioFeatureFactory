"""
Unit tests for alphafold3/bin/binding_metrics.py:
has_confident_binding, classify_binding_event, compute_delta_metrics

Run with: pytest test/unit/test_binding_metrics.py -v
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "alphafold3" / "bin"))
from binding_metrics import (
    BindingMetrics,
    BindingEventClass,
    RnaEditSpan,
    ThresholdConfig,
    aggregate_mutation_summary,
    has_confident_binding,
    classify_binding_event,
    compute_delta_metrics,
    format_sites_rows,
    qc_flag_for_deltas,
)



# helpers


def make_metrics(
    contacts=5,
    pae=5.0,
    plddt_rna=70.0,
    plddt_protein=70.0,
    has_binding=True,
):
    return BindingMetrics(
        rbp_name="RBP1",
        chain_pair_pae_min=pae,
        interface_contacts=contacts,
        interface_plddt_rna=plddt_rna,
        interface_plddt_protein=plddt_protein,
        has_binding=has_binding,
    )

CFG = ThresholdConfig()  # defaults: min_contacts=3, max_pae=10, min_plddt=50



# has_confident_binding


class TestHasConfidentBinding:

    def test_confident_binding(self):
        m = make_metrics(contacts=5, pae=5.0, plddt_rna=70.0)
        assert has_confident_binding(m, CFG) is True

    def test_too_few_contacts(self):
        m = make_metrics(contacts=2, pae=5.0, plddt_rna=70.0)
        assert has_confident_binding(m, CFG) is False

    def test_pae_too_high(self):
        m = make_metrics(contacts=5, pae=15.0, plddt_rna=70.0)
        assert has_confident_binding(m, CFG) is False

    def test_low_plddt_rna_but_high_protein_passes(self):
        # Either plddt_rna OR plddt_protein can satisfy the threshold
        m = make_metrics(contacts=5, pae=5.0, plddt_rna=20.0, plddt_protein=70.0)
        assert has_confident_binding(m, CFG) is True

    def test_both_plddt_below_threshold_fails(self):
        m = make_metrics(contacts=5, pae=5.0, plddt_rna=20.0, plddt_protein=20.0)
        assert has_confident_binding(m, CFG) is False

    def test_exactly_at_contact_threshold_passes(self):
        m = make_metrics(contacts=3, pae=5.0, plddt_rna=70.0)
        assert has_confident_binding(m, CFG) is True

    def test_exactly_at_pae_threshold_passes(self):
        m = make_metrics(contacts=5, pae=10.0, plddt_rna=70.0)
        assert has_confident_binding(m, CFG) is True



# classify_binding_event


class TestClassifyBindingEvent:

    BINDING = make_metrics(contacts=5, pae=5.0, plddt_rna=70.0)
    NO_BIND = make_metrics(contacts=0, pae=20.0, plddt_rna=10.0)

    def test_incomplete_when_wt_none(self):
        result = classify_binding_event(None, self.BINDING, CFG)
        assert result == BindingEventClass.INCOMPLETE

    def test_incomplete_when_mut_none(self):
        result = classify_binding_event(self.BINDING, None, CFG)
        assert result == BindingEventClass.INCOMPLETE

    def test_no_binding_when_both_absent(self):
        result = classify_binding_event(self.NO_BIND, self.NO_BIND, CFG)
        assert result == BindingEventClass.NO_BINDING

    def test_gained_when_only_mut_binds(self):
        result = classify_binding_event(self.NO_BIND, self.BINDING, CFG)
        assert result == BindingEventClass.GAINED

    def test_lost_when_only_wt_binds(self):
        result = classify_binding_event(self.BINDING, self.NO_BIND, CFG)
        assert result == BindingEventClass.LOST

    def test_strengthened_when_pae_drops_significantly(self):
        wt = make_metrics(pae=9.0)
        mut = make_metrics(pae=5.0)  # delta_pae = -4.0 < -2.0 threshold
        result = classify_binding_event(wt, mut, CFG)
        assert result == BindingEventClass.STRENGTHENED

    def test_weakened_when_pae_rises_significantly(self):
        wt = make_metrics(pae=5.0)
        mut = make_metrics(pae=9.0)  # delta_pae = +4.0 > +2.0 threshold
        result = classify_binding_event(wt, mut, CFG)
        assert result == BindingEventClass.WEAKENED

    def test_unchanged_when_no_significant_difference(self):
        wt = make_metrics(pae=5.0, contacts=5)
        mut = make_metrics(pae=5.5, contacts=5)  # delta_pae = 0.5, within threshold
        result = classify_binding_event(wt, mut, CFG)
        assert result == BindingEventClass.UNCHANGED

    def test_strengthened_by_contacts_when_pae_similar(self):
        wt = make_metrics(pae=5.0, contacts=3)
        mut = make_metrics(pae=5.5, contacts=6)  # delta_contacts = +3 >= 2
        result = classify_binding_event(wt, mut, CFG)
        assert result == BindingEventClass.STRENGTHENED



# compute_delta_metrics


class TestComputeDeltaMetrics:

    def test_both_present_computes_deltas(self):
        wt = make_metrics(pae=8.0, contacts=3)
        mut = make_metrics(pae=5.0, contacts=6)
        dm = compute_delta_metrics("RBP1", wt, mut)
        assert dm.delta_chain_pair_pae_min == pytest.approx(5.0 - 8.0)
        assert dm.delta_interface_contacts == 3

    def test_wt_none_uses_mut_values(self):
        mut = make_metrics(pae=5.0, contacts=4)
        dm = compute_delta_metrics("RBP1", None, mut)
        assert dm.delta_chain_pair_pae_min == pytest.approx(5.0)
        assert dm.delta_interface_contacts == 4

    def test_mut_none_negates_wt_values(self):
        wt = make_metrics(pae=5.0, contacts=4)
        dm = compute_delta_metrics("RBP1", wt, None)
        assert dm.delta_chain_pair_pae_min == pytest.approx(-5.0)
        assert dm.delta_interface_contacts == -4

    def test_both_none_returns_zeros(self):
        dm = compute_delta_metrics("RBP1", None, None)
        assert dm.delta_chain_pair_pae_min == pytest.approx(0.0)
        assert dm.delta_interface_contacts == 0

    def test_event_class_set(self):
        wt = make_metrics(pae=9.0)
        mut = make_metrics(pae=5.0)
        dm = compute_delta_metrics("RBP1", wt, mut)
        assert dm.event_class == BindingEventClass.STRENGTHENED

    def test_rbp_name_preserved(self):
        dm = compute_delta_metrics("MY_RBP", make_metrics(), make_metrics())
        assert dm.rbp_name == "MY_RBP"

    def test_low_confidence_contacts_are_not_counted_without_mutating_inputs(self):
        wt = make_metrics(contacts=8, pae=15.0, has_binding=True)
        mut = make_metrics(contacts=8, pae=5.0, has_binding=True)

        dm = compute_delta_metrics("RBP1", wt, mut, config=CFG)
        summary = aggregate_mutation_summary([dm])

        assert dm.wt_metrics is not wt
        assert dm.mut_metrics is not mut
        assert dm.wt_metrics.has_binding is False
        assert dm.mut_metrics.has_binding is True
        assert wt.has_binding is True
        assert mut.has_binding is True
        assert summary['n_rbps_binding_wt'] == 0
        assert summary['n_rbps_binding_mut'] == 1

    def test_has_binding_uses_the_supplied_threshold_config(self):
        config = ThresholdConfig(max_pae_binding=12.0)
        wt = make_metrics(pae=11.0, has_binding=False)
        mut = make_metrics(pae=11.0, has_binding=False)

        dm = compute_delta_metrics("RBP1", wt, mut, config=config)

        assert dm.wt_metrics.has_binding is True
        assert dm.mut_metrics.has_binding is True
        assert wt.has_binding is False
        assert mut.has_binding is False


class TestQcFlagForDeltas:

    def test_empty_is_no_rbps_tested(self):
        assert qc_flag_for_deltas([]) == 'no_rbps_tested'

    def test_all_complete_is_pass(self):
        deltas = [
            compute_delta_metrics("RBP1", make_metrics(), make_metrics()),
            compute_delta_metrics("RBP2", make_metrics(), make_metrics()),
        ]
        assert qc_flag_for_deltas(deltas) == 'PASS'

    def test_mixed_complete_and_failed_is_partial(self):
        deltas = [
            compute_delta_metrics("RBP1", make_metrics(), make_metrics()),
            compute_delta_metrics("RBP2", None, None),
        ]
        assert qc_flag_for_deltas(deltas) == 'PARTIAL'

    def test_one_sided_is_partial(self):
        deltas = [compute_delta_metrics("RBP1", make_metrics(), None)]
        assert qc_flag_for_deltas(deltas) == 'PARTIAL'

    def test_all_failed_is_all_failed(self):
        deltas = [
            compute_delta_metrics("RBP1", None, None),
            compute_delta_metrics("RBP2", None, None),
        ]
        assert qc_flag_for_deltas(deltas) == 'ALL_FAILED'


class TestFormatSitesRows:

    @staticmethod
    def site(res_id, chain="R"):
        return SimpleNamespace(
            chain=chain,
            res_id=res_id,
            res_name="A",
            plddt=90.0,
            is_contact=True,
            min_contact_distance=4.0,
        )

    def test_insertion_marks_new_mutant_bases_and_shifts_downstream_rows(self):
        span = RnaEditSpan(
            offset=2, ref_len=1, alt_len=3, wt_len=5, mut_len=7,
        )
        rows = format_sites_rows(
            "PAM-A3AGG", "RBP1", "MUT",
            [self.site(3), self.site(4), self.site(5), self.site(6)],
            edit_span=span,
        )

        assert [row["align_status"] for row in rows] == [
            "aligned", "inserted", "inserted", "aligned",
        ]
        assert [row["res_id_wt_frame"] for row in rows] == [3, "", "", 4]

    def test_deletion_marks_wild_type_bases_without_mutant_counterparts(self):
        span = RnaEditSpan(
            offset=1, ref_len=3, alt_len=1, wt_len=5, mut_len=3,
        )
        rows = format_sites_rows(
            "NPM1-ATC2A", "RBP1", "WT",
            [self.site(2), self.site(3), self.site(4), self.site(5)],
            edit_span=span,
        )

        assert [row["align_status"] for row in rows] == [
            "aligned", "deleted", "deleted", "aligned",
        ]
        assert [row["res_id_wt_frame"] for row in rows] == [2, 3, 4, 5]
