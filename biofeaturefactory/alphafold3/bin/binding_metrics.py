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
from typing import Dict, Optional, List
from enum import Enum

from biofeaturefactory.utils.utility import align_wt_to_mut


@dataclass(frozen=True)
class RnaEditSpan:
    """Where a variant sits inside one RNA window, in window-local 0-based coords.

    Both AF3 drivers (alphafold3_pipeline and burst) build this from the same
    three numbers -- the window-local offset of the first REF base and the two
    allele lengths -- so the alignment they report and the projection they join
    sites rows on cannot drift apart. It lives here, beside format_sites_rows,
    because that is the only function that consumes the projection.

    An SNV is offset/1/1 with wt_len == mut_len, for which the projection is the
    identity and `alignment` reports every base aligned. No SNV behaviour is
    conditional on this type.
    """
    offset: int      # window-local 0-based index of the first REF base
    ref_len: int
    alt_len: int
    wt_len: int      # WT window length
    mut_len: int     # MUT window length == wt_len + (alt_len - ref_len)

    @property
    def wt_to_mut(self) -> List[Optional[int]]:
        """WT index -> MUT index, None where the edit deleted that base."""
        return align_wt_to_mut(self.wt_len, self.offset, self.ref_len, self.alt_len)

    @property
    def length_delta(self) -> int:
        return self.alt_len - self.ref_len

    @property
    def alignment(self) -> str:
        """`aligned_a/u;deleted_d;inserted_i`, u = UNION of both alleles' bases.

        The denominator is (WT bases + inserted bases), never the WT length
        alone: every WT base keeps a counterpart under an insertion however
        large it is, so `aligned N/N` over WT positions would report a 20 nt
        insertion as fully aligned. The bases with no counterpart are all on
        the mutant side.
        """
        n_deleted = sum(1 for j in self.wt_to_mut if j is None)
        n_inserted = max(0, self.length_delta)
        return (f"aligned_{self.wt_len - n_deleted}/{self.wt_len + n_inserted};"
                f"deleted_{n_deleted};inserted_{n_inserted}")


class BindingEventClass(Enum):
    """Classification of binding change events."""
    GAINED = "gained"           # No WT binding, MUT has binding
    LOST = "lost"               # WT has binding, MUT has none
    STRENGTHENED = "strengthened"  # Both have binding, MUT stronger
    WEAKENED = "weakened"       # Both have binding, MUT weaker
    UNCHANGED = "unchanged"     # No significant change
    NO_BINDING = "no_binding"   # Neither WT nor MUT has binding
    INCOMPLETE = "incomplete"   # WT or MUT prediction missing — no comparison possible


@dataclass
class BindingMetrics:
    """Binding metrics for a single RBP at one position."""
    rbp_name: str
    chain_pair_pae_min: float
    interface_contacts: int
    interface_plddt_rna: float
    interface_plddt_protein: float
    has_binding: bool  # Based on confidence thresholds
    # Ensemble fields (None when single-model)
    n_samples: Optional[int] = None
    std_chain_pair_pae_min: Optional[float] = None
    std_interface_contacts: Optional[float] = None
    std_plddt_rna: Optional[float] = None
    std_plddt_protein: Optional[float] = None


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
    n_windows: Optional[int] = None  # Multi-window mode only


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
    # Handle missing data — one or both sides failed, no valid comparison
    if wt_metrics is None or mut_metrics is None:
        return BindingEventClass.INCOMPLETE

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

    # Find top event (exclude incomplete — no valid comparison happened)
    significant_events = [
        d for d in delta_list
        if d.event_class not in [BindingEventClass.UNCHANGED, BindingEventClass.NO_BINDING, BindingEventClass.INCOMPLETE]
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

    # Max absolute delta (exclude incomplete — one-sided deltas are not real comparisons)
    complete_deltas = [d for d in delta_list if d.event_class != BindingEventClass.INCOMPLETE]
    max_abs_delta = max(
        (abs(d.delta_chain_pair_pae_min) for d in complete_deltas),
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

    In multi-window mode the aggregation struct holds across-window statistics,
    so its count/std are emitted as n_windows_used_*/std_pae_across_windows_*
    and the per-sample n_samples_*/std_pae_* columns stay empty.
    """
    rows = []
    for d in delta_list:
        multi_window = bool(d.n_windows)
        n_wt = d.wt_metrics.n_samples if d.wt_metrics and d.wt_metrics.n_samples else ''
        n_mut = d.mut_metrics.n_samples if d.mut_metrics and d.mut_metrics.n_samples else ''
        std_wt = round(d.wt_metrics.std_chain_pair_pae_min, 3) if d.wt_metrics and d.wt_metrics.std_chain_pair_pae_min is not None else ''
        std_mut = round(d.mut_metrics.std_chain_pair_pae_min, 3) if d.mut_metrics and d.mut_metrics.std_chain_pair_pae_min is not None else ''
        row = {
            'pkey': pkey,
            'rbp_name': d.rbp_name,
            'wt_chain_pair_pae_min': round(d.wt_metrics.chain_pair_pae_min, 3) if d.wt_metrics else '',
            'mut_chain_pair_pae_min': round(d.mut_metrics.chain_pair_pae_min, 3) if d.mut_metrics else '',
            'delta_chain_pair_pae_min': round(d.delta_chain_pair_pae_min, 3),
            'wt_interface_contacts': d.wt_metrics.interface_contacts if d.wt_metrics else 0,
            'mut_interface_contacts': d.mut_metrics.interface_contacts if d.mut_metrics else 0,
            'delta_interface_contacts': d.delta_interface_contacts,
            'cls': d.event_class.value,
            'priority': round(d.priority_score, 3),
            'n_samples_wt': '' if multi_window else n_wt,
            'n_samples_mut': '' if multi_window else n_mut,
            'std_pae_wt': '' if multi_window else std_wt,
            'std_pae_mut': '' if multi_window else std_mut,
            'n_windows': d.n_windows if d.n_windows else '',
            'n_windows_used_wt': n_wt if multi_window else '',
            'n_windows_used_mut': n_mut if multi_window else '',
            'std_pae_across_windows_wt': std_wt if multi_window else '',
            'std_pae_across_windows_mut': std_mut if multi_window else '',
        }
        rows.append(row)
    return rows


def format_sites_rows(
    pkey: str,
    rbp_name: str,
    allele: str,
    sites: List,
    contact_frequency_rna: Optional[Dict[int, float]] = None,
    contact_frequency_protein: Optional[Dict[int, float]] = None,
    *,
    edit_span: RnaEditSpan,
    rna_chain: str = 'R'
) -> List[dict]:
    """
    Format interface sites as rows for sites.tsv.

    Args:
        pkey: Mutation pkey
        rbp_name: RBP identifier
        allele: 'WT' or 'MUT'
        sites: List of InterfaceSite objects
        contact_frequency_rna: Optional dict of res_id -> fraction for RNA chain
        contact_frequency_protein: Optional dict of res_id -> fraction for protein chain
        edit_span: REQUIRED RnaEditSpan for the window this allele was folded
            in. Supplies the WT->MUT projection and both window lengths.
        rna_chain: chain ID whose res_id lives in the RNA window frame. The
            protein chain is the same molecule in both alleles, so its res_id
            carries no edit and needs no projection.

    `res_id` is the ONE AF3 output column that encodes a residue correspondence;
    every other metric this module emits is a whole-complex scalar. For a
    residue 3' of an indel the MUT res_id is the WT res_id shifted by the length
    delta, so joining a WT row to a MUT row on res_id compares two different
    bases. edit_span is therefore REQUIRED rather than defaulted to the
    identity: a silent identity default is exactly the mis-join this argument
    exists to prevent, and the house convention for a parameter with no safe
    default is to fail loud (see utility.get_mutation_data_bioAccurate).
    For an SNV the projection IS the identity, so nothing changes.

    Two derived columns make the cross-allele join safe:
        res_id_wt_frame  WT-window coordinate of this residue -- the only key on
                         which a WT row and a MUT row may be joined. EMPTY for a
                         base the ALT inserted: it has no WT counterpart, and
                         reusing its own number would fabricate one.
        align_status     aligned | deleted (WT base the edit removed, no MUT
                         counterpart) | inserted (MUT base with no WT
                         counterpart) | res_id_outside_window (res_id not
                         addressable in its own window -- reported, not guessed)
    """
    wt_to_mut = edit_span.wt_to_mut
    mut_to_wt = {j: i for i, j in enumerate(wt_to_mut) if j is not None}
    rows = []
    for s in sites:
        freq = ''
        freq_dict = contact_frequency_rna if s.chain == rna_chain else contact_frequency_protein
        if freq_dict and s.res_id in freq_dict:
            freq = round(freq_dict[s.res_id], 3)

        if s.chain != rna_chain:
            # Protein chain: identical sequence in the WT and MUT jobs.
            res_id_wt_frame, align_status = s.res_id, 'aligned'
        elif allele == 'WT':
            i0 = s.res_id - 1
            if 0 <= i0 < len(wt_to_mut):
                res_id_wt_frame = s.res_id
                align_status = 'aligned' if wt_to_mut[i0] is not None else 'deleted'
            else:
                res_id_wt_frame, align_status = '', 'res_id_outside_window'
        else:
            j0 = s.res_id - 1
            if not (0 <= j0 < edit_span.mut_len):
                res_id_wt_frame, align_status = '', 'res_id_outside_window'
            else:
                i0 = mut_to_wt.get(j0)
                res_id_wt_frame = '' if i0 is None else i0 + 1
                align_status = 'inserted' if i0 is None else 'aligned'

        rows.append({
            'pkey': pkey,
            'rbp_name': rbp_name,
            'allele': allele,
            'chain': s.chain,
            'res_id': s.res_id,
            'res_name': s.res_name,
            'plddt': round(s.plddt, 1),
            'is_contact': 1 if s.is_contact else 0,
            'min_contact_distance': s.min_contact_distance,
            'contact_frequency': freq,
            'res_id_wt_frame': res_id_wt_frame,
            'align_status': align_status,
        })
    return rows
