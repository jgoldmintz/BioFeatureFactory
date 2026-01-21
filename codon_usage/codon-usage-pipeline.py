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
    get_codon_counts,
    extract_codon_with_bicodons,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
    compute_cai,
    compute_tai,
    get_codon_tai,
    get_codon_cai_w,
)


FIELDNAMES = [
    'pkey',
    'Gene',
    'codon_number',
    'codon',
    'position_in_codon',
    'RSCU',
    'W',
    'CAI_W',
    'tAI',
    'CAI_gene',
    'tAI_gene',
    'bicodon_3prime',
    'RSCPU_3prime',
    'CPS_3prime',
    'noln_CPS_3prime',
    'W_CP_3prime',
    'bicodon_5prime',
    'RSCPU_5prime',
    'CPS_5prime',
    'noln_CPS_5prime',
    'W_CP_5prime',
    'bicodon_context',
    'qc_flags',
]


def process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene):
    """
    Process a single mutation and return codon usage statistics.

    Args:
        gene: Gene symbol
        ntposnt: Mutation string (e.g., "A123G")
        sequence: ORF nucleotide sequence
        codondata: Pre-computed codon statistics dict
        codonpairdata: Pre-computed codon pair statistics dict
        cai_gene: Pre-computed CAI for the gene
        tai_gene: Pre-computed tAI for the gene

    Returns:
        dict: Row data with codon usage metrics, or None if processing fails
    """
    result = extract_codon_with_bicodons(ntposnt, sequence)
    if result[0] is None:
        return None

    mutated_codon, forward_bicodon, reverse_bicodon, poc, pos, codon_number = result

    pkey = f"{gene}-{ntposnt}"
    qc_flags = []

    # Determine bicodon context
    if codon_number == 1 and forward_bicodon and not reverse_bicodon:
        bicodon_context = 'first_codon_3prime_only'
    elif not forward_bicodon and reverse_bicodon:
        bicodon_context = 'last_codon_5prime_only'
    elif forward_bicodon and reverse_bicodon:
        bicodon_context = 'middle_codon_both_directions'
    else:
        bicodon_context = 'insufficient_sequence'
        qc_flags.append('NO_BICODON')

    row_data = {
        'pkey': pkey,
        'Gene': gene,
        'codon_number': codon_number,
        'codon': mutated_codon,
        'position_in_codon': poc + 1,  # 1-based
        'RSCU': codondata['RSCU'].get(mutated_codon),
        'W': codondata['W'].get(mutated_codon),
        'CAI_W': get_codon_cai_w(mutated_codon),
        'tAI': get_codon_tai(mutated_codon),
        'CAI_gene': cai_gene,
        'tAI_gene': tai_gene,
        'bicodon_3prime': forward_bicodon if forward_bicodon else None,
        'RSCPU_3prime': codonpairdata['RSCPU'].get(forward_bicodon) if forward_bicodon else None,
        'CPS_3prime': codonpairdata['CPS'].get(forward_bicodon) if forward_bicodon else None,
        'noln_CPS_3prime': codonpairdata['noln CPS'].get(forward_bicodon) if forward_bicodon else None,
        'W_CP_3prime': codonpairdata['W_CP'].get(forward_bicodon) if forward_bicodon else None,
        'bicodon_5prime': reverse_bicodon if reverse_bicodon else None,
        'RSCPU_5prime': codonpairdata['RSCPU'].get(reverse_bicodon) if reverse_bicodon else None,
        'CPS_5prime': codonpairdata['CPS'].get(reverse_bicodon) if reverse_bicodon else None,
        'noln_CPS_5prime': codonpairdata['noln CPS'].get(reverse_bicodon) if reverse_bicodon else None,
        'W_CP_5prime': codonpairdata['W_CP'].get(reverse_bicodon) if reverse_bicodon else None,
        'bicodon_context': bicodon_context,
        'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
    }

    return row_data


def process_fasta_with_mutations(fasta_path, mutations_path, validation_log=None):
    """
    Process a FASTA file with mutations from a CSV file.

    Args:
        fasta_path: Path to FASTA file containing ORF sequence
        mutations_path: Path to CSV file containing mutations
        validation_log: Optional path to validation log for filtering

    Returns:
        list: List of row dictionaries
    """
    results = []
    fasta = read_fasta(fasta_path)

    if 'ORF' not in fasta:
        print(f"Warning: No 'ORF' found in {fasta_path}", file=sys.stderr)
        return results

    sequence = fasta['ORF']
    gene = extract_gene_from_filename(fasta_path)

    # Load validation failures if provided
    failure_map = load_validation_failures(validation_log) if validation_log else {}

    # Load mutations
    mut_list = trim_muts(mutations_path, log=validation_log, gene_name=gene)

    if not mut_list:
        print(f"Warning: No mutations found in {mutations_path}", file=sys.stderr)
        return results

    # Pre-compute codon statistics once per sequence
    codondata, codonpairdata = get_codon_counts(sequence)

    # Compute gene-level CAI and tAI
    cai_gene = compute_cai(sequence)
    tai_gene = compute_tai(sequence)

    for ntposnt in mut_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        row = process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene)
        if row:
            results.append(row)

    return results


def process_mutant_fasta(fasta_path):
    """
    Process a FASTA file where mutations are encoded in sequence names.

    Expected format: >GENE-MUTATION (e.g., >BRCA1-A123G)

    Args:
        fasta_path: Path to FASTA file

    Returns:
        list: List of row dictionaries
    """
    results = []
    fasta = read_fasta(fasta_path)

    for seq_name, sequence in fasta.items():
        try:
            gene = seq_name.rsplit('-', 1)[0]
            ntposnt = seq_name.rsplit('-', 1)[1]
        except IndexError:
            print(f"Warning: Could not parse mutation from '{seq_name}'", file=sys.stderr)
            continue

        # Compute codon statistics for this sequence
        codondata, codonpairdata = get_codon_counts(sequence)
        cai_gene = compute_cai(sequence)
        tai_gene = compute_tai(sequence)

        row = process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene)
        if row:
            results.append(row)

    return results


def process_directory(fasta_dir, mutations_dir=None, is_mutant=False, validation_log=None):
    """
    Process all FASTA files in a directory.

    Args:
        fasta_dir: Directory containing FASTA files
        mutations_dir: Directory containing mutation CSV files (for WT mode)
        is_mutant: If True, mutations are in sequence names; otherwise use CSV files
        validation_log: Optional path to validation log

    Returns:
        list: Combined results from all files
    """
    results = []
    fasta_extensions = ['*.fasta', '*.fa', '*.fna']

    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(Path(fasta_dir).glob(ext))

    print(f"Found {len(fasta_files)} FASTA file(s) to process")

    for fasta_file in sorted(fasta_files):
        print(f"Processing {fasta_file}...")

        if is_mutant:
            file_results = process_mutant_fasta(str(fasta_file))
        else:
            # Find matching mutations file
            gene = extract_gene_from_filename(str(fasta_file))
            mutations_file = None

            if mutations_dir:
                for csv_file in Path(mutations_dir).glob('*.csv'):
                    if gene.upper() in csv_file.stem.upper():
                        mutations_file = str(csv_file)
                        break

            if not mutations_file:
                print(f"Warning: No mutations file found for {gene}", file=sys.stderr)
                continue

            file_results = process_fasta_with_mutations(
                str(fasta_file),
                mutations_file,
                validation_log
            )

        results.extend(file_results)
        print(f"  Processed {len(file_results)} mutations")

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
        description="BioFeatureFactory: Codon Usage Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process WT FASTA files with mutation CSVs
  python codon-usage-pipeline.py --fasta-dir /path/to/fastas --mutations-dir /path/to/csvs --output output.tsv

  # Process mutant FASTA files (mutations in sequence names)
  python codon-usage-pipeline.py --fasta-dir /path/to/mutant_fastas --is-mutant --output output.tsv

  # Single file processing
  python codon-usage-pipeline.py --fasta /path/to/gene.fasta --mutations /path/to/mutations.csv --output output.tsv

Metrics:
  RSCU       - Relative Synonymous Codon Usage (gene-specific)
  W          - Relative adaptiveness within gene
  CAI_W      - Reference W value for CAI (human highly expressed genes)
  tAI        - tRNA Adaptation Index weight for the codon
  CAI_gene   - Codon Adaptation Index for entire gene
  tAI_gene   - tRNA Adaptation Index for entire gene
  RSCPU      - Relative Synonymous Codon Pair Usage
  CPS        - Codon Pair Score (ln of observed/expected)
"""
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--fasta-dir', help='Directory containing FASTA files')
    input_group.add_argument('--fasta', help='Single FASTA file')

    parser.add_argument('--mutations-dir', help='Directory containing mutation CSV files (for WT mode)')
    parser.add_argument('--mutations', help='Single mutations CSV file')
    parser.add_argument('--is-mutant', action='store_true',
                        help='Mutations are encoded in sequence names (GENE-MUTATION format)')
    parser.add_argument('--validation-log', help='Validation log for filtering failed mutations')

    # Output options
    parser.add_argument('--output', '-o', required=True, help='Output TSV file path')

    args = parser.parse_args()

    # Validate arguments
    if not args.is_mutant:
        if args.fasta_dir and not args.mutations_dir:
            parser.error("--mutations-dir required when using --fasta-dir without --is-mutant")
        if args.fasta and not args.mutations:
            parser.error("--mutations required when using --fasta without --is-mutant")

    # Process files
    if args.fasta_dir:
        results = process_directory(
            args.fasta_dir,
            args.mutations_dir,
            args.is_mutant,
            args.validation_log
        )
    else:
        if args.is_mutant:
            results = process_mutant_fasta(args.fasta)
        else:
            results = process_fasta_with_mutations(
                args.fasta,
                args.mutations,
                args.validation_log
            )

    # Write output
    write_output(results, args.output)

    # Summary
    if results:
        n_with_3prime = sum(1 for r in results if r.get('bicodon_3prime'))
        n_with_5prime = sum(1 for r in results if r.get('bicodon_5prime'))
        n_with_both = sum(1 for r in results if r.get('bicodon_3prime') and r.get('bicodon_5prime'))

        # Get gene-level stats from first result
        cai = results[0].get('CAI_gene')
        tai = results[0].get('tAI_gene')

        print(f"\nSummary:")
        print(f"  Total mutations processed: {len(results)}")
        print(f"  With 3' bicodon: {n_with_3prime}")
        print(f"  With 5' bicodon: {n_with_5prime}")
        print(f"  With both bicodons: {n_with_both}")
        if cai:
            print(f"  Gene CAI: {cai:.4f}")
        if tai:
            print(f"  Gene tAI: {tai:.4f}")


if __name__ == "__main__":
    main()
