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
from typing import Dict, List, Optional, Set
from dataclasses import dataclass

# Local imports
from bin.rbp_database import POSTAR3Database, RBPBindingSite
from bin.rbp_sequence_mapper import RBPSequenceMapper
from bin.af3_runner import AF3Runner, AF3RunnerConfig, ExecutionMode, create_rna_protein_input
from bin.af3_parser import AF3Parser, analyze_binding, BindingAnalysis
from bin.binding_metrics import (
    BindingMetrics, DeltaMetrics, ThresholdConfig,
    compute_delta_metrics, aggregate_mutation_summary, format_events_rows
)

from utils.utility import (
    read_fasta, trim_muts, get_mutation_data_bioAccurate,
    extract_gene_from_filename, subseq
)


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
        window_size: int = 101,
        rbp_window: int = 50,
        validation_log: Optional[str] = None
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
            window_size: RNA window size around mutation (odd number)
            rbp_window: Window to search for RBP binding sites (+/-bp)
            validation_log: Optional validation log for filtering mutations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.window_size = window_size
        self.rbp_window = rbp_window
        self.validation_log = validation_log

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
            execution_mode=ExecutionMode(execution_mode)
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
        mutations_path: str,
        gene_name: Optional[str] = None,
        chrom: Optional[str] = None,
        tx_start: Optional[int] = None,
        strand: str = "+"
    ):
        """
        Process all mutations for a gene.

        Args:
            fasta_path: Path to ORF FASTA
            mutations_path: Path to mutations CSV
            gene_name: Gene symbol (extracted from filename if not provided)
            chrom: Chromosome for POSTAR3 lookup
            tx_start: Transcript start for coordinate conversion
            strand: Strand ('+' or '-')
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
        mutations = trim_muts(mutations_path, self.validation_log, gene_name)
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
                    strand=strand
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
        strand: str
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

        # Extract RNA windows
        wt_window = subseq(transcript_seq, pos_0, self.window_size)
        mut_seq = transcript_seq[:pos_0] + mut_nt + transcript_seq[pos_0+1:]
        mut_window = subseq(mut_seq, pos_0, self.window_size)

        # Find mutation position in window
        window_start = max(0, pos_0 - self.window_size // 2)
        window_center = pos_0 - window_start

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

        # Query RBPs near mutation
        rbp_sites = self._get_nearby_rbps(chrom, nt_pos, tx_start, strand)

        if not rbp_sites:
            # No RBPs in region
            self._add_no_rbps_result(context)
            return

        # Group by RBP
        rbps_to_test = self.rbp_db.group_by_rbp(rbp_sites)
        print(f"    {pkey}: {len(rbps_to_test)} RBPs to test", file=sys.stderr)

        # Run AF3 for each RBP
        delta_list = []
        for rbp_name, sites in rbps_to_test.items():
            delta = self._analyze_rbp_binding(context, rbp_name, sites)
            if delta:
                delta_list.append(delta)

        # Aggregate results
        self._finalize_mutation_results(context, delta_list)

    def _get_nearby_rbps(
        self,
        chrom: Optional[str],
        nt_pos: int,
        tx_start: Optional[int],
        strand: str
    ) -> List[RBPBindingSite]:
        """Query POSTAR3 for RBPs near mutation position."""
        if chrom is None:
            return []

        # Convert transcript position to genomic
        if tx_start is not None:
            if strand == '+':
                genomic_pos = tx_start + nt_pos - 1
            else:
                genomic_pos = tx_start - nt_pos + 1
        else:
            genomic_pos = nt_pos

        return self.rbp_db.query_position(chrom, genomic_pos, self.rbp_window)

    def _analyze_rbp_binding(
        self,
        context: MutationContext,
        rbp_name: str,
        sites: List[RBPBindingSite]
    ) -> Optional[DeltaMetrics]:
        """Run AF3 and analyze binding for one RBP."""
        # Get RBP data (sequence + MSA)
        rbp_data = self.seq_mapper.get_rbp_data(rbp_name)
        if not rbp_data:
            print(f"      {rbp_name}: sequence not found", file=sys.stderr)
            return None

        protein_seq = rbp_data.sequence
        protein_msa = rbp_data.msa_content  # May be None if no MSA available

        # Check token limit (rough estimate)
        total_tokens = len(context.wt_rna_window) + len(protein_seq)
        if total_tokens > 5000:
            print(f"      {rbp_name}: token limit exceeded ({total_tokens})", file=sys.stderr)
            return None

        # Create AF3 inputs (with MSA if available)
        wt_input = create_rna_protein_input(
            job_name=f"{context.pkey}_{rbp_name}_WT",
            rna_seq=context.wt_rna_window,
            protein_seq=protein_seq,
            protein_msa=protein_msa
        )

        mut_input = create_rna_protein_input(
            job_name=f"{context.pkey}_{rbp_name}_MUT",
            rna_seq=context.mut_rna_window,
            protein_seq=protein_seq,
            protein_msa=protein_msa
        )

        # Run AF3
        wt_job = self.af3_runner.submit_job(wt_input, job_id=f"{context.pkey}_{rbp_name}_WT")
        mut_job = self.af3_runner.submit_job(mut_input, job_id=f"{context.pkey}_{rbp_name}_MUT")

        # Parse results
        wt_metrics = self._parse_binding_metrics(wt_job.output_dir, rbp_name) if wt_job.status == "completed" else None
        mut_metrics = self._parse_binding_metrics(mut_job.output_dir, rbp_name) if mut_job.status == "completed" else None

        # Distance to nearest binding site
        distance = min(site.distance_to(context.nt_pos) for site in sites)

        return compute_delta_metrics(
            rbp_name=rbp_name,
            wt_metrics=wt_metrics,
            mut_metrics=mut_metrics,
            distance_to_mutation=distance,
            config=self.threshold_config
        )

    def _parse_binding_metrics(
        self,
        output_dir: Path,
        rbp_name: str
    ) -> Optional[BindingMetrics]:
        """Parse AF3 output to binding metrics."""
        parser = AF3Parser(str(output_dir))
        structure = parser.parse()

        if not structure:
            return None

        binding = analyze_binding(structure, rna_chain="R", protein_chain="P")
        if not binding:
            return None

        return BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=binding.chain_pair_pae_min,
            interface_contacts=binding.n_contacts,
            interface_plddt_rna=binding.interface_plddt_rna,
            interface_plddt_protein=binding.interface_plddt_protein,
            has_binding=binding.n_contacts >= self.threshold_config.min_contacts
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
        summary['qc_flags'] = 'PASS' if delta_list else 'no_rbps_tested'
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

    # Gene inputs
    parser.add_argument('--fasta', help='ORF FASTA file')
    parser.add_argument('--mutations', help='Mutations CSV file')
    parser.add_argument('--fasta-dir', help='Directory of FASTA files')
    parser.add_argument('--mutations-dir', help='Directory of mutation files')

    # Genomic coordinates (optional, for POSTAR3 lookup)
    parser.add_argument('--chrom', help='Chromosome')
    parser.add_argument('--tx-start', type=int, help='Transcript start position')
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

    # Parameters
    parser.add_argument('--window-size', type=int, default=101,
                       help='RNA window size (odd number)')
    parser.add_argument('--rbp-window', type=int, default=50,
                       help='Window for RBP site lookup (+/-bp)')
    parser.add_argument('--validation-log',
                       help='Validation log for filtering mutations')

    args = parser.parse_args()

    # Validate that we have either MSA dir or sequences
    if not args.msa_dir and not args.rbp_sequences:
        parser.error("Provide either --msa-dir or --rbp-sequences")

    # Initialize pipeline
    pipeline = AlphaFold3Pipeline(
        postar_db=args.postar_db,
        rbp_mapping=args.rbp_mapping,
        output_dir=args.output,
        rbp_sequences=args.rbp_sequences,
        msa_dir=args.msa_dir,
        execution_mode=args.execution_mode,
        af3_binary=args.af3_binary,
        window_size=args.window_size,
        rbp_window=args.rbp_window,
        validation_log=args.validation_log
    )

    # Process single gene or directory
    if args.fasta and args.mutations:
        pipeline.process_gene(
            fasta_path=args.fasta,
            mutations_path=args.mutations,
            chrom=args.chrom,
            tx_start=args.tx_start,
            strand=args.strand
        )
    elif args.fasta_dir and args.mutations_dir:
        # Process all genes in directories
        fasta_dir = Path(args.fasta_dir)
        mutations_dir = Path(args.mutations_dir)

        for fasta_file in sorted(fasta_dir.glob('*.fasta')):
            gene_name = extract_gene_from_filename(str(fasta_file))

            # Find matching mutations file
            mut_candidates = list(mutations_dir.glob(f'*{gene_name}*.csv'))
            if not mut_candidates:
                print(f"No mutations file for {gene_name}", file=sys.stderr)
                continue

            pipeline.process_gene(
                fasta_path=str(fasta_file),
                mutations_path=str(mut_candidates[0])
            )
    else:
        parser.error("Provide either --fasta/--mutations or --fasta-dir/--mutations-dir")

    # Write outputs
    pipeline.write_outputs(args.prefix)


if __name__ == '__main__':
    main()
