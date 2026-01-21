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
BioFeatureFactory: Rare Codon Enrichment Pipeline

Wrapper for cg_cotrans library (Copyright William M. Jacobs, GPL v3).
Detects regions enriched/depleted in rare codons using sliding window analysis.

Requires:
  - cg_cotrans library: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
  - Codon-aware MSA (FASTA)
  - Pre-computed codon usage file (.p.gz from calc_codon_usage.py)

Output format: TSV with BFF-standard pkey column
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.utility import (
    read_fasta,
    trim_muts,
    get_mutation_data_bioAccurate,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
)

# Import cg_cotrans library (GPL v3 licensed, Copyright 2017 William M. Jacobs)
# Download from: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
from cg_cotrans.calc_rare_enrichment import (
    calc_rare_enrichment,
    read_fasta as rc_read_fasta,
    clean_sequences,
    sorted_gis,
    get_codons,
    align_sequences,
    aa_identity,
    load_null_model,
    msa_rare_codon_analysis_wtalign_nseq,
)


FIELDNAMES = [
    'pkey',
    'Gene',
    'codon_position',
    'p_enriched',
    'p_depleted',
    'f_enriched_wt',
    'frac_seq_enriched',
    'frac_seq_depleted',
    'n_rare',
    'window_size',
    'qc_flags',
]


def run_rare_codon_analysis(gene, msa_path, usage_path, wt_gi, window_size=15,
                            rare_model='no_norm', rare_threshold=0.1,
                            null_model='genome', max_len_diff=0.2, min_aa_iden=0.5):
    """
    Run rare codon enrichment analysis on an MSA.

    Args:
        gene: Gene name
        msa_path: Path to codon-aware MSA FASTA
        usage_path: Path to pre-computed codon usage file (.p.gz)
        wt_gi: Identifier for the wildtype/focus sequence
        window_size: Sliding window width in codons
        rare_model: 'no_norm' or 'cmax_norm'
        rare_threshold: Threshold for defining rare codons
        null_model: 'genome', 'eq', or 'groups'
        max_len_diff: Max relative sequence length difference vs WT
        min_aa_iden: Min amino acid identity vs WT

    Returns:
        dict: Position -> {p_enriched, p_depleted, f_enriched_wt, etc.}
    """
    # Load and clean MSA
    seqs = rc_read_fasta(msa_path)
    seqs = clean_sequences(seqs)

    if wt_gi not in seqs:
        raise ValueError(f"WT sequence '{wt_gi}' not found in MSA")

    # Filter sequences by length
    wt_len = len(seqs[wt_gi]) - seqs[wt_gi].count('-')
    for gi in list(seqs.keys()):
        if gi == wt_gi:
            continue
        gi_len = len(seqs[gi]) - seqs[gi].count('-')
        wt_gi_overlap = sum(1 for i in range(len(seqs[wt_gi]))
                           if seqs[wt_gi][i] != '-' and seqs[gi][i] != '-')
        if abs(1. - gi_len / wt_len) > max_len_diff or \
           1. - wt_gi_overlap / wt_len > max_len_diff:
            del seqs[gi]

    # Filter by AA identity
    gis = sorted_gis(seqs, wt_gi)
    msa_codons = {gi: get_codons(seq) for gi, seq in seqs.items()}
    aa_perc_id = aa_identity(msa_codons, wt_gi)
    for gi in list(seqs.keys()):
        if gi == wt_gi:
            continue
        if aa_perc_id.get(gi, 0) < min_aa_iden:
            del seqs[gi]

    # Rebuild after filtering
    gis = sorted_gis(seqs, wt_gi)
    msa_codons = {gi: get_codons(seq) for gi, seq in seqs.items()}
    msa_index_dict = align_sequences(msa_codons, wt_gi)

    # Load null model
    rare_codons, gene_codon_usage, rare_codon_prob, _ = load_null_model(
        msa_codons, gis, usage_path, rare_model, False,
        rare_threshold, null_model, gene
    )

    # Run analysis
    rc_analysis = msa_rare_codon_analysis_wtalign_nseq(
        msa_codons, wt_gi, msa_index_dict, rare_codons, rare_codon_prob,
        L=window_size
    )

    # Extract per-position data
    results = {}
    for pos in rc_analysis['p_nseq_enriched'].keys():
        results[pos] = {
            'p_enriched': rc_analysis['p_nseq_enriched'][pos],
            'p_depleted': rc_analysis['p_nseq_depleted'][pos],
            'f_enriched_wt': rc_analysis['f_enriched'][wt_gi].get(pos),
            'frac_seq_enriched': rc_analysis['fseq_enriched'].get(pos),
            'frac_seq_depleted': rc_analysis['fseq_depleted'].get(pos),
            'n_rare': rc_analysis['n_rare'][wt_gi].get(pos, 0),
        }

    return results, len(seqs)


def process_mutations(mutations_list, gene, orf_sequence, rc_results, window_size, failure_map=None):
    """
    Map mutations to rare codon enrichment statistics.

    Args:
        mutations_list: List of mutation strings
        gene: Gene symbol
        orf_sequence: ORF nucleotide sequence
        rc_results: Dict from run_rare_codon_analysis
        window_size: Window size used in analysis
        failure_map: Optional validation failure map

    Returns:
        list: Result dictionaries
    """
    results = []
    failure_map = failure_map or {}

    for ntposnt in mutations_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        qc_flags = []
        pkey = f"{gene}-{ntposnt}"

        # Get mutation position
        pos_data = get_mutation_data_bioAccurate(ntposnt)
        if pos_data[0] is None:
            qc_flags.append('INVALID_MUTATION')
            results.append({
                'pkey': pkey,
                'Gene': gene,
                'codon_position': '',
                'p_enriched': '',
                'p_depleted': '',
                'f_enriched_wt': '',
                'frac_seq_enriched': '',
                'frac_seq_depleted': '',
                'n_rare': '',
                'window_size': window_size,
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
            })
            continue

        nt_pos = pos_data[0]
        codon_pos = (nt_pos - 1) // 3 + 1  # 1-based codon position

        # Look up enrichment data (center position of window)
        rc_data = rc_results.get(codon_pos)

        if rc_data is None:
            # Position may be near edges (not covered by window)
            qc_flags.append('POSITION_NOT_IN_WINDOW')
            results.append({
                'pkey': pkey,
                'Gene': gene,
                'codon_position': codon_pos,
                'p_enriched': '',
                'p_depleted': '',
                'f_enriched_wt': '',
                'frac_seq_enriched': '',
                'frac_seq_depleted': '',
                'n_rare': '',
                'window_size': window_size,
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
            })
            continue

        results.append({
            'pkey': pkey,
            'Gene': gene,
            'codon_position': codon_pos,
            'p_enriched': rc_data['p_enriched'],
            'p_depleted': rc_data['p_depleted'],
            'f_enriched_wt': rc_data['f_enriched_wt'],
            'frac_seq_enriched': rc_data['frac_seq_enriched'],
            'frac_seq_depleted': rc_data['frac_seq_depleted'],
            'n_rare': rc_data['n_rare'],
            'window_size': window_size,
            'qc_flags': 'PASS',
        })

    return results


def write_output(results, output_path):
    """Write results to TSV file."""
    if not results:
        print("No results to write")
        return

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Rare Codon Enrichment Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rare-codon-pipeline.py \\
    --msa /path/to/codon_msa.fasta \\
    --usage /path/to/codon_usage.p.gz \\
    --wt-gi "focus_sequence_id" \\
    --fasta /path/to/gene.fasta \\
    --mutations /path/to/mutations.csv \\
    --output output.tsv

Required preprocessing:
  1. Generate codon-aware MSA (use msa/codon-msa-pipeline.py)
  2. Generate codon usage file: python cg_cotrans/calc_codon_usage.py

Copyright notice:
  cg_cotrans library is Copyright (C) 2017 William M. Jacobs (GPL v3)
  Source: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
"""
    )

    # MSA and codon usage inputs
    parser.add_argument('--msa', required=True, help='Path to codon-aware MSA FASTA')
    parser.add_argument('--usage', required=True, help='Path to codon usage .p.gz file')
    parser.add_argument('--wt-gi', required=True, help='Identifier for WT/focus sequence in MSA')

    # Mutation inputs
    parser.add_argument('--fasta', required=True, help='FASTA file with ORF sequence')
    parser.add_argument('--mutations', required=True, help='Mutations CSV file')
    parser.add_argument('--validation-log', help='Validation log for filtering')

    # Analysis parameters
    parser.add_argument('--window-size', '-L', type=int, default=15,
                        help='Sliding window width in codons (default: 15)')
    parser.add_argument('--rare-model', choices=['no_norm', 'cmax_norm'], default='no_norm',
                        help='Rare codon definition model (default: no_norm)')
    parser.add_argument('--rare-threshold', type=float, default=0.1,
                        help='Threshold for codon rarity (default: 0.1)')
    parser.add_argument('--null-model', choices=['genome', 'eq', 'groups'], default='genome',
                        help='Null model for codon usage (default: genome)')
    parser.add_argument('--max-len-diff', type=float, default=0.2,
                        help='Max relative length difference vs WT (default: 0.2)')
    parser.add_argument('--min-aa-iden', type=float, default=0.5,
                        help='Min AA identity vs WT (default: 0.5)')

    # Output
    parser.add_argument('--output', '-o', required=True, help='Output TSV file path')

    args = parser.parse_args()

    # Load mutations
    gene = extract_gene_from_filename(args.fasta)
    fasta = read_fasta(args.fasta)

    if 'ORF' not in fasta:
        print(f"Error: No 'ORF' found in {args.fasta}", file=sys.stderr)
        sys.exit(1)

    orf_sequence = fasta['ORF']
    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    mut_list = trim_muts(args.mutations, log=args.validation_log, gene_name=gene)

    if not mut_list:
        print(f"Error: No mutations found in {args.mutations}", file=sys.stderr)
        sys.exit(1)

    # Run rare codon analysis
    print(f"Running rare codon enrichment analysis for {gene}...")
    print(f"  MSA: {args.msa}")
    print(f"  Usage: {args.usage}")
    print(f"  Window size: {args.window_size}")

    try:
        rc_results, n_seqs = run_rare_codon_analysis(
            gene=gene,
            msa_path=args.msa,
            usage_path=args.usage,
            wt_gi=args.wt_gi,
            window_size=args.window_size,
            rare_model=args.rare_model,
            rare_threshold=args.rare_threshold,
            null_model=args.null_model,
            max_len_diff=args.max_len_diff,
            min_aa_iden=args.min_aa_iden,
        )
        print(f"  MSA sequences used: {n_seqs}")
        print(f"  Positions analyzed: {len(rc_results)}")
    except Exception as e:
        print(f"Error in rare codon analysis: {e}", file=sys.stderr)
        sys.exit(1)

    # Map mutations to enrichment data
    results = process_mutations(
        mut_list, gene, orf_sequence, rc_results,
        args.window_size, failure_map
    )

    # Write output
    write_output(results, args.output)

    # Summary
    if results:
        n_enriched = sum(1 for r in results if r.get('p_enriched') and float(r['p_enriched']) < 0.05)
        n_depleted = sum(1 for r in results if r.get('p_depleted') and float(r['p_depleted']) < 0.05)
        print(f"\nSummary:")
        print(f"  Total mutations: {len(results)}")
        print(f"  In enriched regions (p<0.05): {n_enriched}")
        print(f"  In depleted regions (p<0.05): {n_depleted}")


if __name__ == "__main__":
    main()
