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

import argparse
import csv
import os
import sys
from pathlib import Path

from biofeaturefactory.utils.utility import (
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
    write_tsv,
)


# write_tsv is called with this list and extrasaction='ignore', so a row key absent
# here is dropped silently — add the column here whenever the row dict gains one.
FIELDNAMES = [
    'pkey',
    'Gene',
    'codon_number',
    'position_in_codon',
    'codon_wt',
    'codon_mut',
    'RSCU_wt',
    'RSCU_mut',
    'delta_RSCU',
    'W_wt',
    'W_mut',
    'delta_W',
    'CAI_W_wt',
    'CAI_W_mut',
    'delta_CAI_W',
    'tAI_wt',
    'tAI_mut',
    'delta_tAI',
    'CAI_gene',
    'tAI_gene',
    'bicodon_3prime_wt',
    'bicodon_3prime_mut',
    'RSCPU_3prime_wt',
    'RSCPU_3prime_mut',
    'CPS_3prime_wt',
    'CPS_3prime_mut',
    'delta_CPS_3prime',
    'noln_CPS_3prime_wt',
    'noln_CPS_3prime_mut',
    'W_CP_3prime_wt',
    'W_CP_3prime_mut',
    'bicodon_5prime_wt',
    'bicodon_5prime_mut',
    'RSCPU_5prime_wt',
    'RSCPU_5prime_mut',
    'CPS_5prime_wt',
    'CPS_5prime_mut',
    'delta_CPS_5prime',
    'noln_CPS_5prime_wt',
    'noln_CPS_5prime_mut',
    'W_CP_5prime_wt',
    'W_CP_5prime_mut',
    'bicodon_context',
    'qc_flags',
]


def process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene,
                     wt_sequence=None):
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

    # WT counterpart. extract_codon_with_bicodons splices ALT into a copy, so the codon
    # it returns is the mutant; the wild-type triplet is the same slice of the unmutated
    # ORF. Safe because that function's REF guard already proved wt_sequence[pos] == REF.
    # Both alleles are looked up in the SAME codondata/codonpairdata, which are built once
    # from the WT ORF — a point mutation must not shift the gene's reference usage table.
    # wt_sequence is None when the caller only has an already-mutated sequence
    # (process_mutant_fasta): the WT is unrecoverable there, so the wt/delta columns are
    # left empty rather than filled with the mutant's own values, which would read as a
    # delta of 0.0 for every row.
    wt_codon = wt_bi3 = wt_bi5 = None
    if wt_sequence:
        start = (codon_number - 1) * 3
        wt_codon = wt_sequence[start:start + 3]
        if forward_bicodon:
            wt_bi3 = wt_codon + forward_bicodon[3:]
        if reverse_bicodon:
            wt_bi5 = reverse_bicodon[:3] + wt_codon
        if len(wt_codon) != 3:
            wt_codon = wt_bi3 = wt_bi5 = None
            qc_flags.append('NO_WT_CODON')

    def _delta(mut_val, wt_val):
        if mut_val is None or wt_val is None:
            return None
        return round(mut_val - wt_val, 6)

    rscu_mut = codondata['RSCU'].get(mutated_codon)
    rscu_wt = codondata['RSCU'].get(wt_codon) if wt_codon else None
    w_mut = codondata['W'].get(mutated_codon)
    w_wt = codondata['W'].get(wt_codon) if wt_codon else None
    caiw_mut = get_codon_cai_w(mutated_codon)
    caiw_wt = get_codon_cai_w(wt_codon) if wt_codon else None
    tai_mut = get_codon_tai(mutated_codon)
    tai_wt = get_codon_tai(wt_codon) if wt_codon else None
    cps3_mut = codonpairdata['CPS'].get(forward_bicodon) if forward_bicodon else None
    cps3_wt = codonpairdata['CPS'].get(wt_bi3) if wt_bi3 else None
    cps5_mut = codonpairdata['CPS'].get(reverse_bicodon) if reverse_bicodon else None
    cps5_wt = codonpairdata['CPS'].get(wt_bi5) if wt_bi5 else None

    row_data = {
        'pkey': pkey,
        'Gene': gene,
        'codon_number': codon_number,
        'position_in_codon': poc + 1,  # 1-based
        'codon_wt': wt_codon,
        'codon_mut': mutated_codon,
        'RSCU_wt': rscu_wt,
        'RSCU_mut': rscu_mut,
        'delta_RSCU': _delta(rscu_mut, rscu_wt),
        'W_wt': w_wt,
        'W_mut': w_mut,
        'delta_W': _delta(w_mut, w_wt),
        'CAI_W_wt': caiw_wt,
        'CAI_W_mut': caiw_mut,
        'delta_CAI_W': _delta(caiw_mut, caiw_wt),
        'tAI_wt': tai_wt,
        'tAI_mut': tai_mut,
        'delta_tAI': _delta(tai_mut, tai_wt),
        'CAI_gene': cai_gene,
        'tAI_gene': tai_gene,
        'bicodon_3prime_wt': wt_bi3,
        'bicodon_3prime_mut': forward_bicodon if forward_bicodon else None,
        'RSCPU_3prime_wt': codonpairdata['RSCPU'].get(wt_bi3) if wt_bi3 else None,
        'RSCPU_3prime_mut': codonpairdata['RSCPU'].get(forward_bicodon) if forward_bicodon else None,
        'CPS_3prime_wt': cps3_wt,
        'CPS_3prime_mut': cps3_mut,
        'delta_CPS_3prime': _delta(cps3_mut, cps3_wt),
        'noln_CPS_3prime_wt': codonpairdata['noln CPS'].get(wt_bi3) if wt_bi3 else None,
        'noln_CPS_3prime_mut': codonpairdata['noln CPS'].get(forward_bicodon) if forward_bicodon else None,
        'W_CP_3prime_wt': codonpairdata['W_CP'].get(wt_bi3) if wt_bi3 else None,
        'W_CP_3prime_mut': codonpairdata['W_CP'].get(forward_bicodon) if forward_bicodon else None,
        'bicodon_5prime_wt': wt_bi5,
        'bicodon_5prime_mut': reverse_bicodon if reverse_bicodon else None,
        'RSCPU_5prime_wt': codonpairdata['RSCPU'].get(wt_bi5) if wt_bi5 else None,
        'RSCPU_5prime_mut': codonpairdata['RSCPU'].get(reverse_bicodon) if reverse_bicodon else None,
        'CPS_5prime_wt': cps5_wt,
        'CPS_5prime_mut': cps5_mut,
        'delta_CPS_5prime': _delta(cps5_mut, cps5_wt),
        'noln_CPS_5prime_wt': codonpairdata['noln CPS'].get(wt_bi5) if wt_bi5 else None,
        'noln_CPS_5prime_mut': codonpairdata['noln CPS'].get(reverse_bicodon) if reverse_bicodon else None,
        'W_CP_5prime_wt': codonpairdata['W_CP'].get(wt_bi5) if wt_bi5 else None,
        'W_CP_5prime_mut': codonpairdata['W_CP'].get(reverse_bicodon) if reverse_bicodon else None,
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

        # F48: one malformed token must not discard the gene. get_mutation_data_bioAccurate
        # does int(ntposnt[1:-1]), and the is_nt guard only inspects the first/last chars,
        # so a token with nucleotide ends but a broken middle (c.76A>T, A76>T) still raises
        # ValueError. Unguarded, it propagated out of this loop and every mutation already
        # collected in `results` was lost with it.
        try:
            row = process_mutation(gene, ntposnt, sequence, codondata, codonpairdata,
                                   cai_gene, tai_gene, wt_sequence=sequence)
        except (ValueError, IndexError) as e:
            print(f"[codon_usage] {gene}-{ntposnt}: skipped malformed token ({e})", file=sys.stderr)
            continue
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

        try:
            row = process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene)
        except (ValueError, IndexError) as e:
            print(f"[codon_usage] {gene}-{ntposnt}: skipped malformed token ({e})", file=sys.stderr)
            continue
        if row:
            results.append(row)

    return results


def process_directory(fasta_dir, mutations_dir=None, is_mutant=False, validation_log=None, output_dir=None):
    """
    Process all FASTA files in a directory, writing per-gene output files.

    Args:
        fasta_dir: Directory containing FASTA files
        mutations_dir: Directory containing mutation CSV files (for WT mode)
        is_mutant: If True, mutations are in sequence names; otherwise use CSV files
        validation_log: Optional path to validation log
        output_dir: Base output directory for per-gene nested output

    Returns:
        list: Combined results from all files (for summary printing)
    """
    results = []
    fasta_extensions = ['*.fasta', '*.fa', '*.fna']

    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(Path(fasta_dir).glob(ext))

    print(f"Found {len(fasta_files)} FASTA file(s) to process")

    for fasta_file in sorted(fasta_files):
        print(f"Processing {fasta_file}...")

        gene = extract_gene_from_filename(str(fasta_file))

        if is_mutant:
            file_results = process_mutant_fasta(str(fasta_file))
        else:
            mutations_file = None

            if mutations_dir:
                gene_up = gene.upper()
                # F50: exact stem match over a SORTED glob. The old substring test
                # (`gene.upper() in stem.upper()`) over an unsorted glob with a
                # first-match break bound a prefix gene (F9) to a superset gene's
                # file (F9A), nondeterministically across machines — unlike the
                # sorted() fasta list above. Accepts `<GENE>.csv` and the
                # `<GENE>_mutations.csv` convention used in Bio_DBs/mappings/mutations.
                for csv_file in sorted(Path(mutations_dir).glob('*.csv')):
                    # Use the shared extractor rather than hardcoded stems: it resolves
                    # 59/59 production files AND the other layouts this repo uses
                    # (combined_<gene>.csv per main.nf:108-112, muts_<gene>.csv, <gene>.csv),
                    # and keeps hyphenated symbols intact (HLA-A, NKX2-1).
                    if (extract_gene_from_filename(csv_file.name) or '').upper() == gene_up:
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

        print(f"  Processed {len(file_results)} mutations")

        if output_dir and file_results:
            out_dir = Path(output_dir) / gene / "CodonUsage"
            out_dir.mkdir(parents=True, exist_ok=True)
            write_output(file_results, str(out_dir / f"{gene}.codon_usage.tsv"))

        results.extend(file_results)

    return results


def write_output(results, output_path):
    """Write results to TSV file."""
    if not results:
        print("No results to write")
        return

    write_tsv(results, output_path, FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(results)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Codon Usage Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file processing
  python codon_usage_pipeline.py --fasta /path/to/gene.fasta --mutations /path/to/mutations.csv --output results/

  # Directory processing
  python codon_usage_pipeline.py --fasta /path/to/fastas --mutations /path/to/mutations --output results/

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
    parser.add_argument('--fasta', required=True, help='FASTA file or directory of FASTA files')
    parser.add_argument('--mutations', help='Mutations CSV file or directory of CSV files')
    parser.add_argument('--mutant', action='store_true',
                        help='Input FASTAs contain already-mutated sequences, each named '
                             '>GENE-MUTATION (e.g. >BRCA1-A123G); no --mutations needed. '
                             'Requires --fasta to be a directory.')
    parser.add_argument('--validation-log', help='Validation log for filtering failed mutations')

    # Output options
    parser.add_argument('--output', '-o', required=True, help='Output base directory')

    args = parser.parse_args()

    # Validate arguments
    if args.mutant and not Path(args.fasta).is_dir():
        parser.error("--mutant requires --fasta to be a directory of mutated FASTA files")
    if not args.mutant:
        if not Path(args.fasta).is_dir() and not args.mutations:
            parser.error("--mutations required when using a single FASTA file")
        if Path(args.fasta).is_dir() and not args.mutations:
            parser.error("--mutations required when using a FASTA directory")

    # Process files
    if Path(args.fasta).is_dir():
        results = process_directory(
            args.fasta,
            args.mutations,
            args.mutant,
            args.validation_log,
            output_dir=args.output,
        )
    else:
        results = process_fasta_with_mutations(
            args.fasta,
            args.mutations,
            args.validation_log
        )

        # Write per-gene output for single-file mode
        gene = extract_gene_from_filename(args.fasta)
        out_dir = Path(args.output) / gene / "CodonUsage"
        out_dir.mkdir(parents=True, exist_ok=True)
        write_output(results, str(out_dir / f"{gene}.codon_usage.tsv"))

    # Summary
    if results:
        n_with_3prime = sum(1 for r in results if r.get('bicodon_3prime'))
        n_with_5prime = sum(1 for r in results if r.get('bicodon_5prime'))
        n_with_both = sum(1 for r in results if r.get('bicodon_3prime') and r.get('bicodon_5prime'))

        # F51: in directory mode `results` spans EVERY gene (results.extend per file),
        # so results[0] is just whichever gene sorted first — printing its CAI/tAI as
        # "Gene CAI/tAI" mislabels one gene's value as global. Report per gene instead.
        # (The per-gene TSVs were always correct; only this stdout line was wrong.)
        per_gene = {}
        for r in results:
            g = r.get('Gene')
            if g is not None and g not in per_gene:
                per_gene[g] = (r.get('CAI_gene'), r.get('tAI_gene'))

        print(f"\nSummary:")
        print(f"  Total mutations processed: {len(results)}")
        print(f"  With 3' bicodon: {n_with_3prime}")
        print(f"  With 5' bicodon: {n_with_5prime}")
        print(f"  With both bicodons: {n_with_both}")
        if len(per_gene) == 1:
            gene_only, (cai, tai) = next(iter(per_gene.items()))
            if cai:
                print(f"  Gene CAI ({gene_only}): {cai:.4f}")
            if tai:
                print(f"  Gene tAI ({gene_only}): {tai:.4f}")
        elif per_gene:
            print(f"  Genes processed: {len(per_gene)} "
                  f"(per-gene CAI/tAI are in each {{GENE}}.codon_usage.tsv)")


if __name__ == "__main__":
    main()
