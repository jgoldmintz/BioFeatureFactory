# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Binding metrics computation and event classification.

Computes delta metrics between WT and MUT AF3 predictions
and classifies binding events.
"""

from dataclasses import dataclass
from typing import Optional, List
from enum import Enum


class BindingEventClass(Enum):
    """Classification of binding change events."""
    GAINED = "gained"           # No WT binding, MUT has binding
    LOST = "lost"               # WT has binding, MUT has none
    STRENGTHENED = "strengthened"  # Both have binding, MUT stronger
    WEAKENED = "weakened"       # Both have binding, MUT weaker
    UNCHANGED = "unchanged"     # No significant change
    NO_BINDING = "no_binding"   # Neither WT nor MUT has binding


@dataclass
class BindingMetrics:
    """Binding metrics for a single RBP at one position."""
    rbp_name: str
    chain_pair_pae_min: float
    interface_contacts: int
    interface_plddt_rna: float
    interface_plddt_protein: float
    has_binding: bool  # Based on confidence thresholds


@dataclass
class DeltaMetrics:
    """Delta metrics comparing WT and MUT binding."""
    rbp_name: str
    wt_metrics: Optional[BindingMetrics]
    mut_metrics: Optional[BindingMetrics]

    # Deltas (MUT - WT)
    delta_chain_pair_pae_min: float
    delta_interface_contacts: int
    delta_plddt_rna: float
    delta_plddt_protein: float

    # Classification
    event_class: BindingEventClass
    priority_score: float  # For ranking importance


@dataclass
class ThresholdConfig:
    """Thresholds for binding classification."""
    # Minimum interface contacts to consider binding
    min_contacts: int = 3
    # Maximum PAE to consider confident binding
    max_pae_binding: float = 10.0
    # Minimum pLDDT for confident interface
    min_plddt_interface: float = 50.0
    # Delta thresholds for classification
    delta_pae_significant: float = 2.0
    delta_contacts_significant: int = 2


def has_confident_binding(
    metrics: BindingMetrics,
    config: ThresholdConfig
) -> bool:
    """
    Determine if binding is confident based on metrics.

    Args:
        metrics: Binding metrics from AF3
        config: Threshold configuration

    Returns:
        True if binding appears confident
    """
    return (
        metrics.interface_contacts >= config.min_contacts and
        metrics.chain_pair_pae_min <= config.max_pae_binding and
        (metrics.interface_plddt_rna >= config.min_plddt_interface or
         metrics.interface_plddt_protein >= config.min_plddt_interface)
    )


def classify_binding_event(
    wt_metrics: Optional[BindingMetrics],
    mut_metrics: Optional[BindingMetrics],
    config: ThresholdConfig
) -> BindingEventClass:
    """
    Classify the binding change event.

    Args:
        wt_metrics: WT binding metrics (None if AF3 failed)
        mut_metrics: MUT binding metrics (None if AF3 failed)
        config: Threshold configuration

    Returns:
        Event classification
    """
    # Handle missing data
    if wt_metrics is None or mut_metrics is None:
        return BindingEventClass.UNCHANGED

    wt_binding = has_confident_binding(wt_metrics, config)
    mut_binding = has_confident_binding(mut_metrics, config)

    # No binding in either
    if not wt_binding and not mut_binding:
        return BindingEventClass.NO_BINDING

    # Gained binding
    if not wt_binding and mut_binding:
        return BindingEventClass.GAINED

    # Lost binding
    if wt_binding and not mut_binding:
        return BindingEventClass.LOST

    # Both have binding - compare strength
    delta_pae = mut_metrics.chain_pair_pae_min - wt_metrics.chain_pair_pae_min
    delta_contacts = mut_metrics.interface_contacts - wt_metrics.interface_contacts

    # Lower PAE = better binding, so positive delta = weaker
    if delta_pae < -config.delta_pae_significant:
        return BindingEventClass.STRENGTHENED
    elif delta_pae > config.delta_pae_significant:
        return BindingEventClass.WEAKENED

    # If PAE is similar, use contacts
    if delta_contacts >= config.delta_contacts_significant:
        return BindingEventClass.STRENGTHENED
    elif delta_contacts <= -config.delta_contacts_significant:
        return BindingEventClass.WEAKENED

    return BindingEventClass.UNCHANGED


def compute_delta_metrics(
    rbp_name: str,
    wt_metrics: Optional[BindingMetrics],
    mut_metrics: Optional[BindingMetrics],
    distance_to_mutation: int = 0,
    config: Optional[ThresholdConfig] = None
) -> DeltaMetrics:
    """
    Compute delta metrics between WT and MUT binding.

    Args:
        rbp_name: RBP identifier
        wt_metrics: WT binding metrics
        mut_metrics: MUT binding metrics
        distance_to_mutation: Distance from mutation to RBP binding site
        config: Threshold configuration

    Returns:
        DeltaMetrics with deltas and classification
    """
    if config is None:
        config = ThresholdConfig()

    # Compute deltas
    if wt_metrics and mut_metrics:
        delta_pae = mut_metrics.chain_pair_pae_min - wt_metrics.chain_pair_pae_min
        delta_contacts = mut_metrics.interface_contacts - wt_metrics.interface_contacts
        delta_plddt_rna = mut_metrics.interface_plddt_rna - wt_metrics.interface_plddt_rna
        delta_plddt_protein = mut_metrics.interface_plddt_protein - wt_metrics.interface_plddt_protein
    elif mut_metrics:
        delta_pae = mut_metrics.chain_pair_pae_min
        delta_contacts = mut_metrics.interface_contacts
        delta_plddt_rna = mut_metrics.interface_plddt_rna
        delta_plddt_protein = mut_metrics.interface_plddt_protein
    elif wt_metrics:
        delta_pae = -wt_metrics.chain_pair_pae_min
        delta_contacts = -wt_metrics.interface_contacts
        delta_plddt_rna = -wt_metrics.interface_plddt_rna
        delta_plddt_protein = -wt_metrics.interface_plddt_protein
    else:
        delta_pae = 0.0
        delta_contacts = 0
        delta_plddt_rna = 0.0
        delta_plddt_protein = 0.0

    # Classify event
    event_class = classify_binding_event(wt_metrics, mut_metrics, config)

    # Priority score: |delta| weighted by distance decay
    distance_weight = 1.0 / (1.0 + distance_to_mutation / 50.0)
    priority_score = abs(delta_pae) * distance_weight

    return DeltaMetrics(
        rbp_name=rbp_name,
        wt_metrics=wt_metrics,
        mut_metrics=mut_metrics,
        delta_chain_pair_pae_min=delta_pae,
        delta_interface_contacts=delta_contacts,
        delta_plddt_rna=delta_plddt_rna,
        delta_plddt_protein=delta_plddt_protein,
        event_class=event_class,
        priority_score=priority_score
    )


def aggregate_mutation_summary(
    delta_list: List[DeltaMetrics]
) -> dict:
    """
    Aggregate delta metrics across all RBPs for a single mutation.

    Args:
        delta_list: List of DeltaMetrics for all RBPs

    Returns:
        Summary dict for the mutation
    """
    n_tested = len(delta_list)
    n_binding_wt = sum(1 for d in delta_list if d.wt_metrics and d.wt_metrics.has_binding)
    n_binding_mut = sum(1 for d in delta_list if d.mut_metrics and d.mut_metrics.has_binding)

    # Count events
    event_counts = {e: 0 for e in BindingEventClass}
    for d in delta_list:
        event_counts[d.event_class] += 1

    # Find top event
    significant_events = [
        d for d in delta_list
        if d.event_class not in [BindingEventClass.UNCHANGED, BindingEventClass.NO_BINDING]
    ]

    if significant_events:
        top_event = max(significant_events, key=lambda d: d.priority_score)
        top_rbp = top_event.rbp_name
        top_class = top_event.event_class.value
        top_delta = top_event.delta_chain_pair_pae_min
    else:
        top_rbp = ""
        top_class = "none"
        top_delta = 0.0

    # Max absolute delta
    max_abs_delta = max(
        (abs(d.delta_chain_pair_pae_min) for d in delta_list),
        default=0.0
    )

    return {
        'n_rbps_tested': n_tested,
        'n_rbps_binding_wt': n_binding_wt,
        'n_rbps_binding_mut': n_binding_mut,
        'global_count_gained': event_counts[BindingEventClass.GAINED],
        'global_count_lost': event_counts[BindingEventClass.LOST],
        'global_count_strengthened': event_counts[BindingEventClass.STRENGTHENED],
        'global_count_weakened': event_counts[BindingEventClass.WEAKENED],
        'global_max_abs_delta_pae': round(max_abs_delta, 3),
        'top_event_rbp': top_rbp,
        'top_event_class': top_class,
        'top_event_delta_pae': round(top_delta, 3)
    }


def format_events_rows(
    pkey: str,
    delta_list: List[DeltaMetrics]
) -> List[dict]:
    """
    Format delta metrics as rows for events.tsv output.

    Args:
        pkey: Mutation pkey
        delta_list: List of DeltaMetrics

    Returns:
        List of row dicts for events.tsv
    """
    rows = []
    for d in delta_list:
        rows.append({
            'pkey': pkey,
            'rbp_name': d.rbp_name,
            'wt_chain_pair_pae_min': round(d.wt_metrics.chain_pair_pae_min, 3) if d.wt_metrics else '',
            'mut_chain_pair_pae_min': round(d.mut_metrics.chain_pair_pae_min, 3) if d.mut_metrics else '',
            'delta_chain_pair_pae_min': round(d.delta_chain_pair_pae_min, 3),
            'wt_interface_contacts': d.wt_metrics.interface_contacts if d.wt_metrics else 0,
            'mut_interface_contacts': d.mut_metrics.interface_contacts if d.mut_metrics else 0,
            'delta_interface_contacts': d.delta_interface_contacts,
            'cls': d.event_class.value,
            'priority': round(d.priority_score, 3)
        })
    return rows
