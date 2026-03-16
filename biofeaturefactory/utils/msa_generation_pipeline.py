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
BioFeatureFactory: MSA Generation Pipeline

Generates high-quality protein multiple sequence alignments suitable for
EVmutation/PLMC analysis using jackhmmer iterative search against UniRef90.

Output format: A2M (compatible with PLMC and EVmutation)
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from utility import (
    prepare_protein_query, extract_gene_from_filename,
    run_jackhmmer, parse_stockholm, stockholm_to_a2m,
    filter_msa_by_gaps, compute_sequence_weights, compute_neff,
    write_a2m,
)


def get_query_length(query_fasta):
    """Get length of query sequence from FASTA file."""
    record = next(SeqIO.parse(query_fasta, 'fasta'))
    return len(record.seq)


def get_focus_id(query_fasta):
    """Get sequence ID from query FASTA file."""
    record = next(SeqIO.parse(query_fasta, 'fasta'))
    return record.id


def generate_msa(query_fasta, database, output_dir, jackhmmer_binary,
                 iterations=5, evalue_inclusion=1e-3, bitscore_threshold=0.5,
                 min_neff_ratio=10, max_seq_gaps=0.4, max_col_gaps=0.6,
                 identity_threshold=0.8, threads=4, keep_intermediate=False):
    """
    Generate high-quality MSA using jackhmmer.

    Args:
        query_fasta: Path to query protein sequence
        database: Path to sequence database
        output_dir: Output base directory (writes {GENE}/MSA/{GENE}.msa.a2m and .stats.json)
        jackhmmer_binary: Path to jackhmmer executable
        iterations: Number of jackhmmer iterations
        evalue_inclusion: E-value for sequence inclusion
        bitscore_threshold: Minimum bits per residue (for quality gate)
        min_neff_ratio: Minimum N_eff / L ratio
        max_seq_gaps: Maximum fraction gaps per sequence
        max_col_gaps: Maximum fraction gaps per column
        identity_threshold: Clustering threshold for N_eff
        threads: Number of CPU threads
        keep_intermediate: Keep intermediate Stockholm file

    Returns:
        dict: MSA statistics
    """
    # Translate nucleotide query to protein if necessary
    protein_fasta, tmp_protein = prepare_protein_query(query_fasta)

    try:
        query_length = get_query_length(protein_fasta)
        focus_id = get_focus_id(protein_fasta)
        min_neff = min_neff_ratio * query_length

        print(f"Query: {focus_id}, Length: {query_length}")
        print(f"Target N_eff: >={min_neff:.0f} (ratio {min_neff_ratio}xL)")

        # Resolve nested output path
        gene = extract_gene_from_filename(query_fasta) or Path(query_fasta).stem
        nested_dir = Path(output_dir) / gene / "MSA"
        nested_dir.mkdir(parents=True, exist_ok=True)
        output_path = nested_dir / f"{gene}.msa.a2m"
        sto_path = str(nested_dir / f"{gene}.msa.sto")

        # Run jackhmmer
        run_jackhmmer(
            query_fasta=protein_fasta,
            database=database,
            output_sto=sto_path,
            jackhmmer_binary=jackhmmer_binary,
            iterations=iterations,
            evalue_inclusion=evalue_inclusion,
            threads=threads
        )
    finally:
        if tmp_protein and os.path.exists(tmp_protein):
            os.remove(tmp_protein)

    # Parse Stockholm
    print("Parsing Stockholm output...")
    msa, rf_annotation = parse_stockholm(sto_path)
    n_raw = len(msa)
    print(f"Raw sequences: {n_raw}")

    if n_raw == 0:
        raise RuntimeError("No sequences found in jackhmmer output")

    # Convert to A2M format
    print("Converting to A2M format...")
    a2m_msa = stockholm_to_a2m(msa, focus_id, rf_annotation)

    # Filter by gaps
    print(f"Filtering MSA (max_seq_gaps={max_seq_gaps}, max_col_gaps={max_col_gaps})...")
    filtered_msa = filter_msa_by_gaps(a2m_msa, max_seq_gaps, max_col_gaps, a2m_format=True, focus_id=focus_id)
    n_filtered = len(filtered_msa)
    print(f"Filtered sequences: {n_filtered}")

    if n_filtered == 0:
        raise RuntimeError("All sequences filtered out")

    # Rename focus sequence to gene name derived from query filename
    gene_name = extract_gene_from_filename(query_fasta) or focus_id
    if gene_name != focus_id and focus_id in filtered_msa:
        filtered_msa[gene_name] = filtered_msa.pop(focus_id)
        focus_id = gene_name

    # Compute N_eff
    print("Computing N_eff...")
    neff = compute_neff(filtered_msa, identity_threshold)
    neff_ratio = neff / query_length
    print(f"N_eff: {neff:.1f} (ratio: {neff_ratio:.1f}xL)")

    # Check quality threshold
    if neff < min_neff:
        print(f"Warning: N_eff ({neff:.1f}) below threshold ({min_neff:.0f})", file=sys.stderr)
        print("EVmutation predictions may have reduced accuracy", file=sys.stderr)

    # Compute coverage
    if filtered_msa:
        filtered_len = len(next(iter(filtered_msa.values())))
        coverage = filtered_len / query_length
    else:
        coverage = 0.0

    # Compute mean identity to focus
    focus_seq = filtered_msa.get(focus_id, '')
    identities = []
    for seq_id, seq in filtered_msa.items():
        if seq_id != focus_id and len(seq) == len(focus_seq):
            matches = sum(1 for a, b in zip(seq, focus_seq) if a == b and a not in '-.')
            aligned = sum(1 for a, b in zip(seq, focus_seq) if a not in '-.' and b not in '-.')
            if aligned > 0:
                identities.append(matches / aligned)
    mean_identity = np.mean(identities) if identities else 0.0

    # Write A2M output
    print(f"Writing A2M output to {output_path}...")
    write_a2m(filtered_msa, output_path, focus_id)

    # Clean up intermediate files
    if not keep_intermediate and os.path.exists(sto_path):
        os.remove(sto_path)

    # Compile statistics
    stats = {
        'query_id': focus_id,
        'query_length': query_length,
        'n_sequences_raw': n_raw,
        'n_sequences_filtered': n_filtered,
        'n_eff': float(round(neff, 1)),
        'n_eff_ratio': float(round(neff_ratio, 2)),
        'coverage': float(round(coverage, 3)),
        'mean_identity': float(round(mean_identity, 3)),
        'bitscore_threshold_used': bitscore_threshold,
        'iterations': iterations,
        'evalue_inclusion': evalue_inclusion,
        'max_seq_gaps': max_seq_gaps,
        'max_col_gaps': max_col_gaps,
        'identity_threshold': identity_threshold,
        'quality_pass': bool(neff >= min_neff)
    }

    # Write stats JSON
    stats_path = str(output_path).replace('.msa.a2m', '.msa.stats.json')
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Statistics written to {stats_path}")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: MSA Generation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python msa_generation_pipeline.py \\
    --query protein.fasta \\
    --database /path/to/uniref90.fasta \\
    --jackhmmer-binary jackhmmer \\
    --output results/
  # Writes: results/GENE/MSA/GENE.msa.a2m
  #         results/GENE/MSA/GENE.msa.stats.json

  # With custom parameters
  python msa_generation_pipeline.py \\
    --query protein.fasta \\
    --database /path/to/uniref90.fasta \\
    --jackhmmer-binary /usr/local/bin/jackhmmer \\
    --iterations 5 \\
    --min-neff-ratio 10 \\
    --threads 8 \\
    --output results/

Quality thresholds:
  N_eff >= 10L (L = protein length) is recommended for EVmutation.
  Lower N_eff may result in reduced prediction accuracy.
"""
    )

    # Required arguments
    parser.add_argument('--query', '-q', required=True,
                        help='Query protein sequence (FASTA format)')
    parser.add_argument('--database', '-d', required=True,
                        help='Sequence database (e.g., UniRef90)')
    parser.add_argument('--jackhmmer-binary', '-j', required=True,
                        help='Path to jackhmmer executable')
    parser.add_argument('--output', '-o', required=True,
                        help='Output base directory (writes {GENE}/MSA/{GENE}.msa.a2m and .stats.json)')

    # Optional parameters
    parser.add_argument('--iterations', '-N', type=int, default=5,
                        help='Number of jackhmmer iterations (default: 5)')
    parser.add_argument('--evalue-inclusion', '-E', type=float, default=1e-3,
                        help='E-value threshold for inclusion (default: 1e-3)')
    parser.add_argument('--bitscore-threshold', type=float, default=0.5,
                        help='Minimum bits per residue (default: 0.5)')
    parser.add_argument('--min-neff-ratio', type=float, default=10,
                        help='Minimum N_eff / L ratio (default: 10)')
    parser.add_argument('--max-seq-gaps', type=float, default=0.4,
                        help='Max fraction gaps per sequence (default: 0.4)')
    parser.add_argument('--max-col-gaps', type=float, default=0.6,
                        help='Max fraction gaps per column (default: 0.6)')
    parser.add_argument('--identity-threshold', type=float, default=0.8,
                        help='Clustering threshold for N_eff (default: 0.8)')
    parser.add_argument('--threads', '-t', type=int, default=4,
                        help='Number of CPU threads (default: 4)')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help='Keep intermediate Stockholm file')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.query):
        parser.error(f"Query file not found: {args.query}")
    if not os.path.exists(args.database):
        parser.error(f"Database not found: {args.database}")

    # Run pipeline
    stats = generate_msa(
        query_fasta=args.query,
        database=args.database,
        output_dir=args.output,
        jackhmmer_binary=args.jackhmmer_binary,
        iterations=args.iterations,
        evalue_inclusion=args.evalue_inclusion,
        bitscore_threshold=args.bitscore_threshold,
        min_neff_ratio=args.min_neff_ratio,
        max_seq_gaps=args.max_seq_gaps,
        max_col_gaps=args.max_col_gaps,
        identity_threshold=args.identity_threshold,
        threads=args.threads,
        keep_intermediate=args.keep_intermediate
    )

    # Print summary
    print("\n" + "="*50)
    print("MSA Generation Complete")
    print("="*50)
    print(f"Query: {stats['query_id']}")
    print(f"Query length: {stats['query_length']}")
    print(f"Sequences (raw): {stats['n_sequences_raw']}")
    print(f"Sequences (filtered): {stats['n_sequences_filtered']}")
    print(f"N_eff: {stats['n_eff']} ({stats['n_eff_ratio']}xL)")
    print(f"Coverage: {stats['coverage']:.1%}")
    print(f"Mean identity: {stats['mean_identity']:.1%}")
    print(f"Quality check: {'PASS' if stats['quality_pass'] else 'WARN (low N_eff)'}")
    print(f"Output dir: {args.output}")


if __name__ == "__main__":
    main()
