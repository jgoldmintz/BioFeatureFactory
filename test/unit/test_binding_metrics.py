"""
Unit tests for alphafold3/bin/binding_metrics.py:
has_confident_binding, classify_binding_event, compute_delta_metrics

Run with: pytest test/unit/test_binding_metrics.py -v
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "biofeaturefactory" / "alphafold3" / "bin"))
from binding_metrics import (
    BindingMetrics,
    BindingEventClass,
    ThresholdConfig,
    has_confident_binding,
    classify_binding_event,
    compute_delta_metrics,
)



# helpers


def make_metrics(contacts=5, pae=5.0, plddt_rna=70.0, plddt_protein=70.0):
    return BindingMetrics(
        rbp_name="RBP1",
        chain_pair_pae_min=pae,
        interface_contacts=contacts,
        interface_plddt_rna=plddt_rna,
        interface_plddt_protein=plddt_protein,
        has_binding=True,
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
