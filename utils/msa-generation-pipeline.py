#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023–2026  Jacob Goldmintz
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
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_jackhmmer(query_fasta, database, output_sto, jackhmmer_binary,
                  iterations=5, evalue_inclusion=1e-3, threads=4):
    """
    Run jackhmmer iterative search.

    Args:
        query_fasta: Path to query protein sequence (FASTA)
        database: Path to UniRef90 or similar database
        output_sto: Path for Stockholm output
        jackhmmer_binary: Path to jackhmmer executable
        iterations: Number of search iterations
        evalue_inclusion: E-value threshold for inclusion
        threads: Number of CPU threads

    Returns:
        Path to Stockholm output file
    """
    cmd = [
        jackhmmer_binary,
        '-N', str(iterations),
        '--incE', str(evalue_inclusion),
        '-A', output_sto,
        '--cpu', str(threads),
        '--noali',  # Don't output alignment to stdout
        query_fasta,
        database
    ]

    print(f"Running jackhmmer with {iterations} iterations...")
    print(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"jackhmmer stderr: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"jackhmmer failed with exit code {result.returncode}")

    print(f"jackhmmer complete. Output: {output_sto}")
    return Path(output_sto)


def parse_stockholm(stockholm_file):
    """
    Parse Stockholm format MSA file.

    Args:
        stockholm_file: Path to Stockholm file

    Returns:
        dict: {seq_id: sequence} mapping
    """
    msa = {}
    current_seqs = defaultdict(str)

    with open(stockholm_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Skip comments and empty lines
            if line.startswith('#') or line.startswith('//') or not line:
                continue

            # Parse sequence lines
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                current_seqs[seq_id] += seq

    # Convert to final dict
    for seq_id, seq in current_seqs.items():
        msa[seq_id] = seq

    return msa


def stockholm_to_a2m(msa, focus_seq_id):
    """
    Convert Stockholm MSA to A2M format.

    A2M format:
    - Uppercase: match states (aligned to query)
    - Lowercase: insertions relative to query
    - '-': deletions (gaps in sequence, not in query)
    - '.': gaps in query (insertions in other sequences)

    Args:
        msa: dict {seq_id: sequence}
        focus_seq_id: ID of the focus/query sequence

    Returns:
        dict: {seq_id: a2m_sequence}
    """
    if focus_seq_id not in msa:
        # Try partial match
        for seq_id in msa:
            if focus_seq_id in seq_id or seq_id in focus_seq_id:
                focus_seq_id = seq_id
                break
        else:
            raise ValueError(f"Focus sequence '{focus_seq_id}' not found in MSA")

    focus_seq = msa[focus_seq_id]

    # Identify match columns (non-gap in focus sequence)
    match_columns = [i for i, c in enumerate(focus_seq) if c not in '-.' ]

    a2m_msa = {}
    for seq_id, seq in msa.items():
        a2m_seq = []
        for i, c in enumerate(seq):
            if i in match_columns:
                # Match column: uppercase
                if c in '-.':
                    a2m_seq.append('-')
                else:
                    a2m_seq.append(c.upper())
            else:
                # Insert column: lowercase or '.'
                if c in '-.':
                    a2m_seq.append('.')
                else:
                    a2m_seq.append(c.lower())

        a2m_msa[seq_id] = ''.join(a2m_seq)

    return a2m_msa


def filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=0.6):
    """
    Filter MSA by removing gappy sequences and columns.

    Args:
        msa: dict {seq_id: sequence}
        max_seq_gaps: Maximum fraction of gaps allowed per sequence
        max_col_gaps: Maximum fraction of gaps allowed per column

    Returns:
        dict: Filtered MSA
    """
    if not msa:
        return {}

    # Step 1: Remove sequences with too many gaps
    filtered_seqs = {}
    for seq_id, seq in msa.items():
        gap_count = seq.count('-') + seq.count('.') + seq.count('!')
        gap_frac = gap_count / len(seq) if len(seq) > 0 else 1.0
        if gap_frac <= max_seq_gaps:
            filtered_seqs[seq_id] = seq

    if not filtered_seqs:
        print("Warning: All sequences filtered out by gap threshold", file=sys.stderr)
        return {}

    # Step 2: Identify columns with too many gaps
    seq_list = list(filtered_seqs.values())
    seq_len = len(seq_list[0])
    n_seqs = len(seq_list)

    cols_to_keep = []
    for i in range(seq_len):
        col = [s[i] for s in seq_list]
        gap_count = sum(1 for c in col if c in '-.')
        gap_frac = gap_count / n_seqs
        if gap_frac <= max_col_gaps:
            cols_to_keep.append(i)

    # Step 3: Remove gappy columns
    final_msa = {}
    for seq_id, seq in filtered_seqs.items():
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        final_msa[seq_id] = new_seq

    return final_msa


def compute_sequence_weights(msa, identity_threshold=0.8):
    """
    Compute sequence weights based on clustering at identity threshold.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold (0.8 = 80% identity)

    Returns:
        dict: {seq_id: weight}
    """
    seq_ids = list(msa.keys())
    n_seqs = len(seq_ids)

    if n_seqs == 0:
        return {}

    # Convert to numpy array for efficient comparison
    seqs = [msa[sid] for sid in seq_ids]
    seq_len = len(seqs[0])

    # Compute pairwise identities and count neighbors
    neighbor_counts = np.ones(n_seqs)  # Each sequence counts itself

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            # Compute identity
            matches = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a == b and a not in '-.')
            aligned = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a not in '-.' and b not in '-.')

            if aligned > 0:
                identity = matches / aligned
                if identity >= identity_threshold:
                    neighbor_counts[i] += 1
                    neighbor_counts[j] += 1

    # Weight = 1 / number of neighbors
    weights = {seq_ids[i]: 1.0 / neighbor_counts[i] for i in range(n_seqs)}
    return weights


def compute_neff(msa, identity_threshold=0.8):
    """
    Compute effective number of sequences (N_eff).

    N_eff = sum of sequence weights, where weight = 1/n_neighbors

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold

    Returns:
        float: N_eff value
    """
    weights = compute_sequence_weights(msa, identity_threshold)
    return sum(weights.values())


def get_query_length(query_fasta):
    """Get length of query sequence from FASTA file."""
    record = next(SeqIO.parse(query_fasta, 'fasta'))
    return len(record.seq)


def get_focus_id(query_fasta):
    """Get sequence ID from query FASTA file."""
    record = next(SeqIO.parse(query_fasta, 'fasta'))
    return record.id


def write_a2m(msa, output_path, focus_seq_id=None):
    """
    Write MSA to A2M format file.

    Args:
        msa: dict {seq_id: sequence}
        output_path: Output file path
        focus_seq_id: If provided, write focus sequence first
    """
    with open(output_path, 'w') as f:
        # Write focus sequence first if specified
        if focus_seq_id and focus_seq_id in msa:
            f.write(f">{focus_seq_id}\n{msa[focus_seq_id]}\n")

        # Write remaining sequences
        for seq_id, seq in msa.items():
            if seq_id != focus_seq_id:
                f.write(f">{seq_id}\n{seq}\n")


def generate_msa(query_fasta, database, output_path, jackhmmer_binary,
                 iterations=5, evalue_inclusion=1e-3, bitscore_threshold=0.5,
                 min_neff_ratio=10, max_seq_gaps=0.4, max_col_gaps=0.6,
                 identity_threshold=0.8, threads=4, keep_intermediate=False):
    """
    Generate high-quality MSA using jackhmmer.

    Args:
        query_fasta: Path to query protein sequence
        database: Path to sequence database
        output_path: Path for output A2M file
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
    query_length = get_query_length(query_fasta)
    focus_id = get_focus_id(query_fasta)
    min_neff = min_neff_ratio * query_length

    print(f"Query: {focus_id}, Length: {query_length}")
    print(f"Target N_eff: ≥{min_neff:.0f} (ratio {min_neff_ratio}×L)")

    # Create temporary file for Stockholm output
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    sto_path = str(output_path).replace('.a2m', '.sto')

    # Run jackhmmer
    run_jackhmmer(
        query_fasta=query_fasta,
        database=database,
        output_sto=sto_path,
        jackhmmer_binary=jackhmmer_binary,
        iterations=iterations,
        evalue_inclusion=evalue_inclusion,
        threads=threads
    )

    # Parse Stockholm
    print("Parsing Stockholm output...")
    msa = parse_stockholm(sto_path)
    n_raw = len(msa)
    print(f"Raw sequences: {n_raw}")

    if n_raw == 0:
        raise RuntimeError("No sequences found in jackhmmer output")

    # Convert to A2M format
    print("Converting to A2M format...")
    a2m_msa = stockholm_to_a2m(msa, focus_id)

    # Filter by gaps
    print(f"Filtering MSA (max_seq_gaps={max_seq_gaps}, max_col_gaps={max_col_gaps})...")
    filtered_msa = filter_msa_by_gaps(a2m_msa, max_seq_gaps, max_col_gaps)
    n_filtered = len(filtered_msa)
    print(f"Filtered sequences: {n_filtered}")

    if n_filtered == 0:
        raise RuntimeError("All sequences filtered out")

    # Compute N_eff
    print("Computing N_eff...")
    neff = compute_neff(filtered_msa, identity_threshold)
    neff_ratio = neff / query_length
    print(f"N_eff: {neff:.1f} (ratio: {neff_ratio:.1f}×L)")

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
        'n_eff': round(neff, 1),
        'n_eff_ratio': round(neff_ratio, 2),
        'coverage': round(coverage, 3),
        'mean_identity': round(mean_identity, 3),
        'bitscore_threshold_used': bitscore_threshold,
        'iterations': iterations,
        'evalue_inclusion': evalue_inclusion,
        'max_seq_gaps': max_seq_gaps,
        'max_col_gaps': max_col_gaps,
        'identity_threshold': identity_threshold,
        'quality_pass': neff >= min_neff
    }

    # Write stats JSON
    stats_path = str(output_path).replace('.a2m', '.stats.json')
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
  python msa-generation-pipeline.py \\
    --query protein.fasta \\
    --database /path/to/uniref90.fasta \\
    --jackhmmer-binary jackhmmer \\
    --output protein.msa.a2m

  # With custom parameters
  python msa-generation-pipeline.py \\
    --query protein.fasta \\
    --database /path/to/uniref90.fasta \\
    --jackhmmer-binary /usr/local/bin/jackhmmer \\
    --iterations 5 \\
    --min-neff-ratio 10 \\
    --threads 8 \\
    --output protein.msa.a2m

Quality thresholds:
  N_eff ≥ 10L (L = protein length) is recommended for EVmutation.
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
                        help='Output MSA file (A2M format)')

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
        output_path=args.output,
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
    print(f"N_eff: {stats['n_eff']} ({stats['n_eff_ratio']}×L)")
    print(f"Coverage: {stats['coverage']:.1%}")
    print(f"Mean identity: {stats['mean_identity']:.1%}")
    print(f"Quality check: {'PASS' if stats['quality_pass'] else 'WARN (low N_eff)'}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
