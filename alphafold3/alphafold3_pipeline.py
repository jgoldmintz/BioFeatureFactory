#!/usr/bin/env python3
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
AlphaFold3 RNA-RBP Interaction Pipeline

Analyzes how synonymous mutations affect RNA-binding protein interactions
using AlphaFold3 structure predictions.

Outputs:
- summary.tsv: One row per mutation with aggregated RBP binding changes
- events.tsv: One row per mutation-RBP pair with detailed metrics
- sites.tsv: Per-position confidence scores at RNA-protein interface
"""

import argparse
import csv
import sys
from pathlib import Path
from concurrent.futures import Future
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field

# Local imports
from bin.rbp_database import POSTAR3Database, RBPBindingSite
from bin.rbp_sequence_mapper import RBPSequenceMapper
from bin.af3_runner import AF3Runner, AF3RunnerConfig, ExecutionMode, create_rna_protein_input
from bin.af3_parser import (
    AF3Parser, AF3Structure, analyze_binding, BindingAnalysis,
    parse_all_samples, aggregate_binding_analyses, AggregatedBindingAnalysis,
    extract_interface_sites
)
from bin.binding_metrics import (
    BindingMetrics, DeltaMetrics, ThresholdConfig,
    compute_delta_metrics, aggregate_mutation_summary,
    format_events_rows, format_sites_rows
)

from utils.utility import (
    read_fasta, trim_muts, get_mutation_data_bioAccurate,
    extract_gene_from_filename, subseq, load_mapping,
    _collect_failures_from_logs
)


def parse_vcf_chrom(vcf_path: str) -> Optional[str]:
    """Extract CHROM from a per-gene VCF (all data rows share same chromosome)."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                return fields[0]
    return None


@dataclass
class MutationContext:
    """Context for a single mutation."""
    pkey: str
    gene: str
    mutation: str
    nt_pos: int  # 1-based nucleotide position
    wt_nt: str
    mut_nt: str
    transcript_seq: str
    wt_rna_window: str
    mut_rna_window: str
    window_center: int  # Position of mutation in window (0-based)
    genomic_pos: Optional[int] = None  # Chromosomal position for RBP distance


@dataclass
class _ParsedResult:
    """Internal: parsed AF3 output with metrics, structures, and aggregation."""
    metrics: Optional[BindingMetrics]
    structures: List[AF3Structure]
    aggregation: Optional[AggregatedBindingAnalysis]


class AlphaFold3Pipeline:
    """
    Main pipeline for AF3 RNA-RBP interaction analysis.
    """

    def __init__(
        self,
        postar_db: str,
        rbp_mapping: str,
        output_dir: str,
        rbp_sequences: Optional[str] = None,
        msa_dir: Optional[str] = None,
        execution_mode: str = "local",
        af3_binary: str = "alphafold3",
        docker_image: str = "alphafold3",
        model_dir: Optional[str] = None,
        window_size: int = 101,
        rbp_window: int = 50,
        validation_log: Optional[str] = None,
        multi_window: bool = False,
        multi_window_offsets: Optional[List[float]] = None,
        max_gpus: Optional[int] = None
    ):
        """
        Initialize pipeline.

        Args:
            postar_db: Path to POSTAR3 database file
            rbp_mapping: Path to gene-UniProt mapping TSV
            output_dir: Output directory
            rbp_sequences: Path to protein sequences FASTA (optional if msa_dir provided)
            msa_dir: Directory containing A3M MSA files (preferred over rbp_sequences)
            execution_mode: 'local', 'batch', or 'cloud'
            af3_binary: Path to AF3 executable
            docker_image: Docker image name for AF3
            model_dir: Path to AF3 model weights directory
            window_size: RNA window size around mutation (odd number)
            rbp_window: Window to search for RBP binding sites (+/-bp)
            validation_log: Optional validation log for filtering mutations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.window_size = window_size
        self.rbp_window = rbp_window
        self.validation_log = validation_log
        self.multi_window = multi_window
        self.multi_window_offsets = multi_window_offsets or [0.3, 0.5, 0.7]

        # Initialize components
        print("Loading POSTAR3 database...", file=sys.stderr)
        self.rbp_db = POSTAR3Database(postar_db)

        print("Loading RBP sequence mapper...", file=sys.stderr)
        self.seq_mapper = RBPSequenceMapper(
            mapping_file=rbp_mapping,
            sequence_fasta=rbp_sequences,
            msa_dir=msa_dir
        )

        # AF3 runner config
        af3_config = AF3RunnerConfig(
            af3_binary=af3_binary,
            output_base_dir=str(self.output_dir / "af3_runs"),
            execution_mode=ExecutionMode(execution_mode),
            docker_image=docker_image,
            model_dir=model_dir,
            max_gpus=max_gpus
        )
        self.af3_runner = AF3Runner(af3_config)

        # Thresholds
        self.threshold_config = ThresholdConfig()

        # Results storage
        self.summary_rows: List[dict] = []
        self.events_rows: List[dict] = []
        self.sites_rows: List[dict] = []

    def process_gene(
        self,
        fasta_path: str,
        mutations_path: Optional[str] = None,
        gene_name: Optional[str] = None,
        chrom: Optional[str] = None,
        tx_start: Optional[int] = None,
        strand: str = "+",
        chrom_mapping: Optional[Dict[str, str]] = None
    ):
        """
        Process all mutations for a gene.

        Args:
            fasta_path: Path to ORF FASTA
            mutations_path: Path to mutations CSV (optional when chrom_mapping provided)
            gene_name: Gene symbol (extracted from filename if not provided)
            chrom: Chromosome for POSTAR3 lookup
            tx_start: Transcript start for coordinate conversion
            strand: Strand ('+' or '-')
            chrom_mapping: Dict mapping mutation -> chromosome entry (e.g., C123T -> C87504250T)
        """
        # Get gene name
        if gene_name is None:
            gene_name = extract_gene_from_filename(fasta_path)

        print(f"\nProcessing {gene_name}...", file=sys.stderr)

        # Load transcript
        fasta_data = read_fasta(fasta_path)
        if not fasta_data:
            print(f"  Error: Empty FASTA file", file=sys.stderr)
            return

        # Get transcript sequence (prefer 'ORF' or 'transcript' header)
        transcript_seq = None
        for key in ['ORF', 'orf', 'transcript', gene_name]:
            if key in fasta_data:
                transcript_seq = fasta_data[key]
                break
        if transcript_seq is None:
            transcript_seq = next(iter(fasta_data.values()))

        # Load mutations
        if mutations_path:
            mutations = trim_muts(mutations_path, self.validation_log, gene_name)
        elif chrom_mapping:
            mutations = list(chrom_mapping.keys())
            if self.validation_log:
                failures = _collect_failures_from_logs(self.validation_log)
                skip_set = failures.get(gene_name.upper(), set()) if gene_name else set()
                mutations = [m for m in mutations if m not in skip_set]
        else:
            print(f"  Error: No mutations source", file=sys.stderr)
            return

        print(f"  {len(mutations)} mutations to process", file=sys.stderr)

        # Process each mutation
        for mutation in mutations:
            try:
                self._process_mutation(
                    gene_name=gene_name,
                    mutation=mutation,
                    transcript_seq=transcript_seq,
                    chrom=chrom,
                    tx_start=tx_start,
                    strand=strand,
                    chrom_mapping=chrom_mapping
                )
            except Exception as e:
                print(f"  Error processing {mutation}: {e}", file=sys.stderr)
                self._add_failed_mutation(gene_name, mutation, str(e))

    def _process_mutation(
        self,
        gene_name: str,
        mutation: str,
        transcript_seq: str,
        chrom: Optional[str],
        tx_start: Optional[int],
        strand: str,
        chrom_mapping: Optional[Dict[str, str]] = None
    ):
        """Process a single mutation."""
        pkey = f"{gene_name}-{mutation}"

        # Parse mutation
        pos_data = get_mutation_data_bioAccurate(mutation)
        if pos_data[0] is None:
            raise ValueError(f"Invalid mutation format: {mutation}")

        nt_pos = pos_data[0]  # 1-based
        wt_nt, mut_nt = pos_data[1]

        # Validate mutation
        pos_0 = nt_pos - 1
        if pos_0 < 0 or pos_0 >= len(transcript_seq):
            raise ValueError(f"Position {nt_pos} out of range")

        if transcript_seq[pos_0].upper() != wt_nt.upper():
            raise ValueError(f"Reference mismatch at {nt_pos}: expected {wt_nt}, found {transcript_seq[pos_0]}")

        # Generate RNA windows
        windows = self._generate_windows(transcript_seq, pos_0, mut_nt)

        # Primary window (centered or first offset)
        wt_window, mut_window, window_center = windows[0]

        context = MutationContext(
            pkey=pkey,
            gene=gene_name,
            mutation=mutation,
            nt_pos=nt_pos,
            wt_nt=wt_nt,
            mut_nt=mut_nt,
            transcript_seq=transcript_seq,
            wt_rna_window=wt_window,
            mut_rna_window=mut_window,
            window_center=window_center
        )

        # Resolve chromosomal position
        genomic_pos = None
        if chrom_mapping and mutation in chrom_mapping:
            entry = chrom_mapping[mutation]
            try:
                genomic_pos = int(entry[1:-1])
            except (ValueError, IndexError):
                print(f"    Warning: Could not parse chromosome mapping for {mutation}: {entry}", file=sys.stderr)
        elif tx_start is not None:
            if strand == '+':
                genomic_pos = tx_start + nt_pos - 1
            else:
                genomic_pos = tx_start - nt_pos + 1

        context.genomic_pos = genomic_pos

        # Query RBPs near mutation
        rbp_sites = self._get_nearby_rbps(chrom, genomic_pos)

        if not rbp_sites:
            # No RBPs in region
            self._add_no_rbps_result(context)
            return

        # Group by RBP
        rbps_to_test = self.rbp_db.group_by_rbp(rbp_sites)
        print(f"    {pkey}: {len(rbps_to_test)} RBPs to test", file=sys.stderr)

        # Phase 1: Submit all RBP jobs (non-blocking)
        pending_list = []
        for rbp_name, sites in rbps_to_test.items():
            pending = self._submit_rbp_jobs(
                context, rbp_name, sites,
                windows=windows if len(windows) > 1 else None
            )
            if pending:
                pending_list.append(pending)

        # Phase 2: Collect all results (blocks on futures as they complete)
        delta_list = []
        for pending in pending_list:
            delta = self._collect_rbp_results(context, pending)
            if delta:
                delta_list.append(delta)

        # Aggregate results
        self._finalize_mutation_results(context, delta_list)

    def _get_nearby_rbps(
        self,
        chrom: Optional[str],
        genomic_pos: Optional[int]
    ) -> List[RBPBindingSite]:
        """Query POSTAR3 for RBPs near mutation position."""
        if chrom is None or genomic_pos is None:
            return []

        return self.rbp_db.query_position(chrom, genomic_pos, self.rbp_window)

    def _generate_windows(
        self,
        transcript_seq: str,
        pos_0: int,
        mut_nt: str
    ) -> List[Tuple[str, str, int]]:
        """
        Generate RNA windows around the mutation.

        Returns list of (wt_window, mut_window, window_center) tuples.
        Single window when multi_window is disabled.
        """
        mut_seq = transcript_seq[:pos_0] + mut_nt + transcript_seq[pos_0 + 1:]

        if not self.multi_window:
            wt_window = subseq(transcript_seq, pos_0, self.window_size)
            mut_window = subseq(mut_seq, pos_0, self.window_size)
            window_start = max(0, pos_0 - self.window_size // 2)
            return [(wt_window, mut_window, pos_0 - window_start)]

        seen = set()
        windows = []
        for frac in self.multi_window_offsets:
            target_center = int(frac * self.window_size)
            window_start = pos_0 - target_center
            window_start = max(0, min(window_start, len(transcript_seq) - self.window_size))
            window_end = window_start + self.window_size

            wt_win = transcript_seq[window_start:window_end]
            mut_win = mut_seq[window_start:window_end]
            center = pos_0 - window_start

            if wt_win not in seen:
                seen.add(wt_win)
                windows.append((wt_win, mut_win, center))

        return windows if windows else [(subseq(transcript_seq, pos_0, self.window_size),
                                          subseq(mut_seq, pos_0, self.window_size),
                                          pos_0 - max(0, pos_0 - self.window_size // 2))]

    @dataclass
    class _PendingRBPAnalysis:
        """Tracks submitted async jobs for one RBP."""
        rbp_name: str
        sites: List
        distance: int
        wt_future: Optional[Future] = None
        mut_future: Optional[Future] = None
        window_wt_futures: Optional[List[Future]] = None
        window_mut_futures: Optional[List[Future]] = None
        n_windows: int = 1

    def _submit_rbp_jobs(
        self,
        context: MutationContext,
        rbp_name: str,
        sites: List[RBPBindingSite],
        windows: Optional[List[Tuple[str, str, int]]] = None
    ) -> Optional['AlphaFold3Pipeline._PendingRBPAnalysis']:
        """Submit AF3 jobs for one RBP without blocking. Returns pending tracker."""
        rbp_data = self.seq_mapper.get_rbp_data(rbp_name)
        if not rbp_data:
            print(f"      {rbp_name}: sequence not found", file=sys.stderr)
            return None

        protein_seq = rbp_data.sequence
        protein_msa = rbp_data.msa_content

        total_tokens = len(context.wt_rna_window) + len(protein_seq)
        if total_tokens > 5000:
            print(f"      {rbp_name}: token limit exceeded ({total_tokens})", file=sys.stderr)
            return None

        distance = min(site.distance_to(context.genomic_pos) for site in sites) if context.genomic_pos else 0

        pending = self._PendingRBPAnalysis(
            rbp_name=rbp_name,
            sites=sites,
            distance=distance
        )

        if windows and len(windows) > 1:
            pending.n_windows = len(windows)
            pending.window_wt_futures = []
            pending.window_mut_futures = []
            for i, (wt_win, mut_win, _center) in enumerate(windows):
                wt_in = create_rna_protein_input(
                    job_name=f"{context.pkey}_{rbp_name}_WT_w{i}",
                    rna_seq=wt_win, protein_seq=protein_seq, protein_msa=protein_msa
                )
                mut_in = create_rna_protein_input(
                    job_name=f"{context.pkey}_{rbp_name}_MUT_w{i}",
                    rna_seq=mut_win, protein_seq=protein_seq, protein_msa=protein_msa
                )
                pending.window_wt_futures.append(
                    self.af3_runner.submit_job_async(wt_in, job_id=f"{context.pkey}_{rbp_name}_WT_w{i}")
                )
                pending.window_mut_futures.append(
                    self.af3_runner.submit_job_async(mut_in, job_id=f"{context.pkey}_{rbp_name}_MUT_w{i}")
                )
        else:
            wt_input = create_rna_protein_input(
                job_name=f"{context.pkey}_{rbp_name}_WT",
                rna_seq=context.wt_rna_window, protein_seq=protein_seq, protein_msa=protein_msa
            )
            mut_input = create_rna_protein_input(
                job_name=f"{context.pkey}_{rbp_name}_MUT",
                rna_seq=context.mut_rna_window, protein_seq=protein_seq, protein_msa=protein_msa
            )
            pending.wt_future = self.af3_runner.submit_job_async(
                wt_input, job_id=f"{context.pkey}_{rbp_name}_WT"
            )
            pending.mut_future = self.af3_runner.submit_job_async(
                mut_input, job_id=f"{context.pkey}_{rbp_name}_MUT"
            )

        return pending

    def _collect_rbp_results(
        self,
        context: MutationContext,
        pending: '_PendingRBPAnalysis'
    ) -> Optional[DeltaMetrics]:
        """Block on futures, parse results, compute delta metrics for one RBP."""
        rbp_name = pending.rbp_name
        distance = pending.distance

        if pending.n_windows > 1:
            wt_results = []
            mut_results = []
            for wt_f, mut_f in zip(pending.window_wt_futures, pending.window_mut_futures):
                wt_job = wt_f.result()
                mut_job = mut_f.result()
                if wt_job.status == "completed" and wt_job.result_path:
                    wt_results.append(self._parse_af3_output(wt_job.result_path, rbp_name))
                if mut_job.status == "completed" and mut_job.result_path:
                    mut_results.append(self._parse_af3_output(mut_job.result_path, rbp_name))

            wt_metrics = self._aggregate_parsed_results(wt_results, rbp_name)
            mut_metrics = self._aggregate_parsed_results(mut_results, rbp_name)

            for r in wt_results:
                if r and r.structures:
                    sites_data = extract_interface_sites(r.structures[0])
                    freq = r.aggregation.contact_frequency_rna if r.aggregation else None
                    self.sites_rows.extend(format_sites_rows(context.pkey, rbp_name, 'WT', sites_data, freq))
                    break
            for r in mut_results:
                if r and r.structures:
                    sites_data = extract_interface_sites(r.structures[0])
                    freq = r.aggregation.contact_frequency_rna if r.aggregation else None
                    self.sites_rows.extend(format_sites_rows(context.pkey, rbp_name, 'MUT', sites_data, freq))
                    break

            delta = compute_delta_metrics(
                rbp_name=rbp_name,
                wt_metrics=wt_metrics,
                mut_metrics=mut_metrics,
                distance_to_mutation=distance,
                config=self.threshold_config
            )
            delta.n_windows = pending.n_windows
            return delta

        # Single-window
        wt_job = pending.wt_future.result()
        mut_job = pending.mut_future.result()

        wt_parsed = self._parse_af3_output(wt_job.result_path, rbp_name) if wt_job.status == "completed" else None
        mut_parsed = self._parse_af3_output(mut_job.result_path, rbp_name) if mut_job.status == "completed" else None

        wt_metrics = wt_parsed.metrics if wt_parsed else None
        mut_metrics = mut_parsed.metrics if mut_parsed else None

        if wt_parsed and wt_parsed.structures:
            sites_data = extract_interface_sites(wt_parsed.structures[0])
            freq = wt_parsed.aggregation.contact_frequency_rna if wt_parsed.aggregation else None
            self.sites_rows.extend(format_sites_rows(context.pkey, rbp_name, 'WT', sites_data, freq))
        if mut_parsed and mut_parsed.structures:
            sites_data = extract_interface_sites(mut_parsed.structures[0])
            freq = mut_parsed.aggregation.contact_frequency_rna if mut_parsed.aggregation else None
            self.sites_rows.extend(format_sites_rows(context.pkey, rbp_name, 'MUT', sites_data, freq))

        return compute_delta_metrics(
            rbp_name=rbp_name,
            wt_metrics=wt_metrics,
            mut_metrics=mut_metrics,
            distance_to_mutation=distance,
            config=self.threshold_config
        )

    def _parse_af3_output(
        self,
        output_dir: Path,
        rbp_name: str
    ) -> Optional[_ParsedResult]:
        """Parse AF3 output with ensemble sample aggregation."""
        if output_dir is None:
            return None

        structures = parse_all_samples(str(output_dir))
        if not structures:
            return None

        analyses = [analyze_binding(s, rna_chain="R", protein_chain="P") for s in structures]

        if len(structures) == 1:
            binding = analyses[0]
            if not binding:
                return _ParsedResult(metrics=None, structures=structures, aggregation=None)
            metrics = BindingMetrics(
                rbp_name=rbp_name,
                chain_pair_pae_min=binding.chain_pair_pae_min,
                interface_contacts=binding.n_contacts,
                interface_plddt_rna=binding.interface_plddt_rna,
                interface_plddt_protein=binding.interface_plddt_protein,
                has_binding=binding.n_contacts >= self.threshold_config.min_contacts
            )
            return _ParsedResult(metrics=metrics, structures=structures, aggregation=None)

        # Multi-sample aggregation
        agg = aggregate_binding_analyses(analyses)
        if not agg:
            return _ParsedResult(metrics=None, structures=structures, aggregation=None)

        metrics = BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=agg.mean.chain_pair_pae_min,
            interface_contacts=agg.mean.n_contacts,
            interface_plddt_rna=agg.mean.interface_plddt_rna,
            interface_plddt_protein=agg.mean.interface_plddt_protein,
            has_binding=agg.mean.n_contacts >= self.threshold_config.min_contacts,
            n_samples=agg.n_samples,
            std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
            std_interface_contacts=agg.std_n_contacts,
            std_plddt_rna=agg.std_interface_plddt_rna,
            std_plddt_protein=agg.std_interface_plddt_protein
        )
        return _ParsedResult(metrics=metrics, structures=structures, aggregation=agg)

    def _aggregate_parsed_results(
        self,
        results: List[Optional[_ParsedResult]],
        rbp_name: str
    ) -> Optional[BindingMetrics]:
        """Aggregate metrics across multiple windows."""
        all_analyses = []
        for r in results:
            if r and r.metrics:
                # Build a BindingAnalysis-like from each result for aggregation
                all_analyses.append(BindingAnalysis(
                    rna_chain="R", protein_chain="P",
                    n_contacts=r.metrics.interface_contacts,
                    min_contact_distance=0.0,
                    mean_contact_distance=0.0,
                    interface_plddt_rna=r.metrics.interface_plddt_rna,
                    interface_plddt_protein=r.metrics.interface_plddt_protein,
                    chain_pair_pae_min=r.metrics.chain_pair_pae_min,
                    contact_residues_rna=[], contact_residues_protein=[]
                ))

        if not all_analyses:
            return None

        agg = aggregate_binding_analyses(all_analyses)
        if not agg:
            return None

        return BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=agg.mean.chain_pair_pae_min,
            interface_contacts=agg.mean.n_contacts,
            interface_plddt_rna=agg.mean.interface_plddt_rna,
            interface_plddt_protein=agg.mean.interface_plddt_protein,
            has_binding=agg.mean.n_contacts >= self.threshold_config.min_contacts,
            n_samples=agg.n_samples,
            std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
            std_interface_contacts=agg.std_n_contacts,
            std_plddt_rna=agg.std_interface_plddt_rna,
            std_plddt_protein=agg.std_interface_plddt_protein
        )

    def _finalize_mutation_results(
        self,
        context: MutationContext,
        delta_list: List[DeltaMetrics]
    ):
        """Aggregate and store results for a mutation."""
        # Summary row
        summary = aggregate_mutation_summary(delta_list)
        summary['pkey'] = context.pkey
        summary['Gene'] = context.gene

        # QC flags: check whether AF3 predictions actually succeeded
        has_complete = any(d.wt_metrics is not None and d.mut_metrics is not None for d in delta_list)
        has_partial = any(d.wt_metrics is not None or d.mut_metrics is not None for d in delta_list)
        if not delta_list:
            summary['qc_flags'] = 'no_rbps_tested'
        elif has_complete:
            summary['qc_flags'] = 'PASS'
        elif has_partial:
            summary['qc_flags'] = 'PARTIAL'
        else:
            summary['qc_flags'] = 'ALL_FAILED'

        self.summary_rows.append(summary)

        # Events rows
        events = format_events_rows(context.pkey, delta_list)
        self.events_rows.extend(events)

    def _add_no_rbps_result(self, context: MutationContext):
        """Add result for mutation with no nearby RBPs."""
        self.summary_rows.append({
            'pkey': context.pkey,
            'Gene': context.gene,
            'n_rbps_tested': 0,
            'n_rbps_binding_wt': 0,
            'n_rbps_binding_mut': 0,
            'global_count_gained': 0,
            'global_count_lost': 0,
            'global_count_strengthened': 0,
            'global_count_weakened': 0,
            'global_max_abs_delta_pae': 0,
            'top_event_rbp': '',
            'top_event_class': 'none',
            'top_event_delta_pae': 0,
            'qc_flags': 'no_rbps_in_region'
        })

    def _add_failed_mutation(self, gene: str, mutation: str, error: str):
        """Add result for failed mutation."""
        pkey = f"{gene}-{mutation}"
        self.summary_rows.append({
            'pkey': pkey,
            'Gene': gene,
            'n_rbps_tested': 0,
            'qc_flags': f'FAILED:{error[:50]}'
        })

    def write_outputs(self, prefix: str = "alphafold3"):
        """Write output TSV files."""
        # Summary
        summary_path = self.output_dir / f"{prefix}.summary.tsv"
        if self.summary_rows:
            fieldnames = list(self.summary_rows[0].keys())
            with open(summary_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(self.summary_rows)
            print(f"Wrote {len(self.summary_rows)} rows to {summary_path}", file=sys.stderr)

        # Events
        events_path = self.output_dir / f"{prefix}.events.tsv"
        if self.events_rows:
            fieldnames = list(self.events_rows[0].keys())
            with open(events_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(self.events_rows)
            print(f"Wrote {len(self.events_rows)} rows to {events_path}", file=sys.stderr)

        # Sites
        sites_path = self.output_dir / f"{prefix}.sites.tsv"
        if self.sites_rows:
            fieldnames = list(self.sites_rows[0].keys())
            with open(sites_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(self.sites_rows)
            print(f"Wrote {len(self.sites_rows)} rows to {sites_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='AlphaFold3 RNA-RBP Interaction Pipeline'
    )

    # Required inputs
    parser.add_argument('--postar-db', required=True,
                       help='POSTAR3 database file')
    parser.add_argument('--rbp-mapping', required=True,
                       help='Gene-UniProt mapping TSV')
    parser.add_argument('--rbp-sequences',
                       help='Protein sequences FASTA (optional if --msa-dir provided)')
    parser.add_argument('--msa-dir',
                       help='Directory with A3M MSA files (preferred over --rbp-sequences)')

    # Gene inputs (file or directory, auto-detected)
    parser.add_argument('--fasta',
                        help='Transcript FASTA file or directory of FASTA files')
    parser.add_argument('--mutations',
                        help='Mutations CSV file or directory (optional when --chromosome-mapping provided)')

    # VCF-based coordinate resolution (replaces --chrom/--tx-start/--strand)
    parser.add_argument('--vcf',
                        help='Per-gene VCF file or directory from vcf_converter.py (provides chromosome)')
    parser.add_argument('--chromosome-mapping',
                        help='Chromosome mapping CSV file or directory (provides mutations and chromosomal positions)')

    # Legacy genomic coordinates (still supported)
    parser.add_argument('--chrom', help='Chromosome (alternative to --vcf)')
    parser.add_argument('--tx-start', type=int, help='Transcript start position (alternative to --chromosome-mapping)')
    parser.add_argument('--strand', default='+', help='Strand (+/-)')

    # Output
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory')
    parser.add_argument('--prefix', default='alphafold3',
                       help='Output file prefix')

    # Execution
    parser.add_argument('--execution-mode', default='local',
                       choices=['local', 'batch', 'cloud'],
                       help='AF3 execution mode')
    parser.add_argument('--af3-binary', default='alphafold3',
                       help='Path to AF3 executable')
    parser.add_argument('--docker-image', default='alphafold3',
                       help='Docker image name for AF3')
    parser.add_argument('--model-dir',
                       help='Path to AF3 model weights directory')

    # Parameters
    parser.add_argument('--window-size', type=int, default=101,
                       help='RNA window size (odd number)')
    parser.add_argument('--rbp-window', type=int, default=50,
                       help='Window for RBP site lookup (+/-bp)')
    parser.add_argument('--validation-log',
                       help='Validation log for filtering mutations')

    # Multi-window mode (optional, multiplies AF3 runs)
    parser.add_argument('--multi-window', action='store_true', default=False,
                       help='Run multiple windows per mutation (multiplies AF3 runs)')
    parser.add_argument('--multi-window-offsets', type=str, default='0.3,0.5,0.7',
                       help='Mutation position as fraction of window (default: 0.3,0.5,0.7)')
    parser.add_argument('--max-gpus', type=int, default=None,
                       help='Max GPUs for parallel AF3 execution (default: auto-detect)')

    args = parser.parse_args()

    # Validate that we have either MSA dir or sequences
    if not args.msa_dir and not args.rbp_sequences:
        parser.error("Provide either --msa-dir or --rbp-sequences")

    # Model weights required for local execution
    if args.execution_mode == 'local' and not args.model_dir:
        parser.error("--model-dir is required for local execution mode")

    # Initialize pipeline
    pipeline = AlphaFold3Pipeline(
        postar_db=args.postar_db,
        rbp_mapping=args.rbp_mapping,
        output_dir=args.output,
        rbp_sequences=args.rbp_sequences,
        msa_dir=args.msa_dir,
        execution_mode=args.execution_mode,
        af3_binary=args.af3_binary,
        docker_image=args.docker_image,
        model_dir=args.model_dir,
        window_size=args.window_size,
        rbp_window=args.rbp_window,
        validation_log=args.validation_log,
        multi_window=args.multi_window,
        multi_window_offsets=[float(x) for x in args.multi_window_offsets.split(',')] if args.multi_window_offsets else None,
        max_gpus=args.max_gpus
    )

    # Resolve inputs
    fasta_input = Path(args.fasta) if args.fasta else None
    mutations_input = Path(args.mutations) if args.mutations else None
    vcf_input = Path(args.vcf) if args.vcf else None
    chrom_map_input = Path(args.chromosome_mapping) if args.chromosome_mapping else None

    if not fasta_input:
        parser.error("--fasta is required")
    if not mutations_input and not chrom_map_input:
        parser.error("Provide --mutations or --chromosome-mapping")

    # When using --mutations without --chromosome-mapping, require genomic coordinate flags
    if mutations_input and not chrom_map_input:
        if not args.chrom or args.tx_start is None:
            parser.error("--chrom and --tx-start are required when using --mutations without --chromosome-mapping")

    if fasta_input.is_dir():
        # --- Directory mode ---
        for fasta_file in sorted(fasta_input.glob('*.fasta')):
            gene_name = extract_gene_from_filename(str(fasta_file))

            # Find matching mutations file
            mut_path = None
            if mutations_input and mutations_input.is_dir():
                candidates = list(mutations_input.glob(f'*{gene_name}*.csv'))
                if candidates:
                    mut_path = str(candidates[0])

            # Resolve chromosome from VCF
            chrom = args.chrom
            if vcf_input:
                vcf_search = vcf_input if vcf_input.is_file() else None
                if vcf_input.is_dir():
                    vcf_candidates = list(vcf_input.glob(f'{gene_name}.vcf'))
                    if vcf_candidates:
                        vcf_search = vcf_candidates[0]
                if vcf_search:
                    chrom = parse_vcf_chrom(str(vcf_search))

            # Load chromosome mapping
            chrom_mapping = None
            if chrom_map_input:
                if chrom_map_input.is_file():
                    chrom_mapping = load_mapping(str(chrom_map_input), mapType="chromosome")
                elif chrom_map_input.is_dir():
                    map_candidates = list(chrom_map_input.glob(f'*{gene_name}*.csv'))
                    if map_candidates:
                        chrom_mapping = load_mapping(str(map_candidates[0]), mapType="chromosome")

            if not mut_path and not chrom_mapping:
                print(f"No mutations source for {gene_name}", file=sys.stderr)
                continue

            pipeline.process_gene(
                fasta_path=str(fasta_file),
                mutations_path=mut_path,
                chrom=chrom,
                tx_start=args.tx_start,
                strand=args.strand,
                chrom_mapping=chrom_mapping
            )

    else:
        # --- Single file mode ---
        chrom = args.chrom
        chrom_mapping = None

        if vcf_input:
            chrom = parse_vcf_chrom(str(vcf_input))

        if chrom_map_input:
            if chrom_map_input.is_file():
                chrom_mapping = load_mapping(str(chrom_map_input), mapType="chromosome")
            elif chrom_map_input.is_dir():
                gene_name = extract_gene_from_filename(str(fasta_input))
                map_candidates = list(chrom_map_input.glob(f'*{gene_name}*.csv'))
                if map_candidates:
                    chrom_mapping = load_mapping(str(map_candidates[0]), mapType="chromosome")

        mut_path = str(mutations_input) if mutations_input else None

        pipeline.process_gene(
            fasta_path=str(fasta_input),
            mutations_path=mut_path,
            chrom=chrom,
            tx_start=args.tx_start,
            strand=args.strand,
            chrom_mapping=chrom_mapping
        )

    # Write outputs
    pipeline.write_outputs(args.prefix)
    pipeline.af3_runner.shutdown()


if __name__ == '__main__':
    main()
