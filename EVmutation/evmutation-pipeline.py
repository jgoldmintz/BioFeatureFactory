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
BioFeatureFactory: EVmutation/PLMC Pipeline

Computes evolutionary mutation effects using epistatic model from protein
multiple sequence alignments.

Requires:
  - EVmutation library: https://github.com/debbiemarkslab/EVmutation
  - plmc binary for model inference

Output format: TSV with BFF-standard pkey column
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# EVmutation library (MIT licensed, Copyright Thomas A. Hopf)
# Clone from: https://github.com/debbiemarkslab/EVmutation
from EVmutation.model import CouplingsModel
import EVmutation.tools as ev_tools
from utils.utility import (
    read_fasta,
    trim_muts,
    get_mutation_data_bioAccurate,
    get_mutant_aa,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
    codon_to_aa,
)


FIELDNAMES = [
    'pkey',
    'mutant',
    'pos',
    'wt',
    'subs',
    'prediction_epistatic',
    'prediction_independent',
    'frequency',
    'column_conservation',
    'qc_flags',
]


def run_plmc(fasta, focus, model_params, plmc_binary, ec_file=None,
             eij_lambda=16.2, hi_lambda=0.01, iterations=500, stepsize=0.2,
             alphabet="-ACDEFGHIKLMNPQRSTVWY", quiet=False):
    """
    Run plmc to generate model parameters.

    Args:
        fasta: Path to MSA FASTA file
        focus: Focus sequence identifier
        model_params: Output path for model parameters
        plmc_binary: Path to plmc binary
        ec_file: Output path for evolutionary couplings (optional)
        eij_lambda: Regularization strength for J_ij (default: 16.2)
        hi_lambda: Regularization strength for h_i (default: 0.01)
        iterations: Maximum iterations (default: 500)
        stepsize: Step size (default: 0.2)
        alphabet: Alphabet string (default: "-ACDEFGHIKLMNPQRSTVWY")
        quiet: Suppress output

    Returns:
        Path to generated model parameters file

    Raises:
        RuntimeError: If plmc fails
    """
    if not fasta:
        raise RuntimeError("Cannot run PLMC without a multiple sequence alignment.")
    if not focus:
        raise RuntimeError("Cannot run PLMC without a focus sequence.")
    if not plmc_binary:
        raise RuntimeError("Cannot run PLMC without specifying --plmc-binary path.")

    if not os.path.exists(plmc_binary):
        raise RuntimeError(f"PLMC binary not found at: {plmc_binary}")

    if ec_file is None:
        ec_file = os.path.join(os.path.dirname(model_params), focus + ".ec")

    if not quiet:
        print(f"Running PLMC on MSA: {fasta}")
        print(f"  Focus sequence: {focus}")
        print(f"  Output params: {model_params}")

    cmd = [
        plmc_binary,
        '-o', model_params,
        '-a', alphabet,
        '-c', ec_file,
        '-f', focus,
        '-le', str(eij_lambda),
        '-lh', str(hi_lambda),
        '-m', str(iterations),
        '-t', str(stepsize),
        '-g', fasta
    ]

    plmc_process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    stdout, stderr = plmc_process.communicate()

    if plmc_process.returncode != 0:
        print(f"PLMC failed with exit code {plmc_process.returncode}", file=sys.stderr)
        print(f"stderr: {stderr.decode()}", file=sys.stderr)
        raise RuntimeError("PLMC did not complete successfully.")

    if not quiet:
        print("PLMC complete.")
        if stdout:
            print(f"stdout: {stdout.decode()}")

    return model_params


def compute_column_conservation(model):
    """
    Compute column conservation from single-site frequencies.

    Conservation = max frequency at each position (excluding gaps).

    Args:
        model: CouplingsModel instance

    Returns:
        dict: {position: conservation_score}
    """
    conservation = {}
    gap_index = model.alphabet_map.get('-', model.alphabet_map.get('.', 0))

    for i, pos in enumerate(model.index_list):
        freqs = model.f_i[i].copy()
        # Exclude gap frequency from conservation calculation
        freqs[gap_index] = 0
        conservation[pos] = float(np.max(freqs))

    return conservation


def generate_single_mutant_predictions(model, gene=None):
    """
    Generate predictions for all possible single amino acid substitutions.

    Args:
        model: CouplingsModel instance
        gene: Gene name for pkey generation

    Returns:
        list: List of result dictionaries
    """
    results = []

    # Get epistatic predictions
    epistatic_df = ev_tools.single_mutant_matrix(model, output_column="prediction_epistatic")

    # Get independent model and predictions
    independent_model = model.to_independent_model()
    independent_df = ev_tools.single_mutant_matrix(independent_model, output_column="prediction_independent")

    # Compute column conservation
    conservation = compute_column_conservation(model)

    # Merge predictions
    for _, row in epistatic_df.iterrows():
        mutant = row['mutant']
        pos = int(row['pos'])
        wt = row['wt']
        subs = row['subs']

        # Get independent prediction for same mutation
        indep_row = independent_df[
            (independent_df['pos'] == pos) & (independent_df['subs'] == subs)
        ]
        indep_pred = indep_row['prediction_independent'].values[0] if len(indep_row) > 0 else np.nan

        pkey = f"{gene}-{mutant}" if gene else mutant

        results.append({
            'pkey': pkey,
            'mutant': mutant,
            'pos': pos,
            'wt': wt,
            'subs': subs,
            'prediction_epistatic': row['prediction_epistatic'],
            'prediction_independent': indep_pred,
            'frequency': row['frequency'],
            'column_conservation': conservation.get(pos, np.nan),
            'qc_flags': 'PASS',
        })

    return results


def map_nt_mutations_to_aa(mutations_list, gene, orf_sequence, evmut_results, failure_map=None):
    """
    Map nucleotide mutations to amino acid EVmutation predictions.

    Args:
        mutations_list: List of nucleotide mutation strings
        gene: Gene symbol
        orf_sequence: ORF nucleotide sequence
        evmut_results: List of EVmutation result dicts
        failure_map: Optional validation failure map

    Returns:
        list: Filtered results matching the mutations
    """
    results = []
    failure_map = failure_map or {}

    # Build lookup by AA mutation
    aa_lookup = {}
    for r in evmut_results:
        aa_lookup[r['mutant']] = r

    for ntposnt in mutations_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        qc_flags = []
        pkey = f"{gene}-{ntposnt}"

        # Get amino acid change from nucleotide mutation
        aa_result = get_mutant_aa((get_mutation_data_bioAccurate(ntposnt)), orf_sequence)

        if aa_result is None or aa_result[0] is None:
            qc_flags.append('INVALID_MUTATION')
            results.append({
                'pkey': pkey,
                'mutant': '',
                'pos': '',
                'wt': '',
                'subs': '',
                'prediction_epistatic': '',
                'prediction_independent': '',
                'frequency': '',
                'column_conservation': '',
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
            })
            continue

        aa_pos_info, _ = aa_result
        aa_pos = aa_pos_info[0]
        wt_aa, mut_aa = aa_pos_info[1]

        # Check if synonymous
        if wt_aa == mut_aa:
            qc_flags.append('SYNONYMOUS')
            # For synonymous, look up the position conservation
            pos_results = [r for r in evmut_results if r['pos'] == aa_pos]
            if pos_results:
                # Use first available for conservation info
                ref = pos_results[0]
                results.append({
                    'pkey': pkey,
                    'mutant': f"{wt_aa}{aa_pos}{mut_aa}",
                    'pos': aa_pos,
                    'wt': wt_aa,
                    'subs': mut_aa,
                    'prediction_epistatic': 0.0,  # No change for synonymous
                    'prediction_independent': 0.0,
                    'frequency': ref['frequency'] if wt_aa == ref['wt'] else '',
                    'column_conservation': ref['column_conservation'],
                    'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
                })
            else:
                results.append({
                    'pkey': pkey,
                    'mutant': f"{wt_aa}{aa_pos}{mut_aa}",
                    'pos': aa_pos,
                    'wt': wt_aa,
                    'subs': mut_aa,
                    'prediction_epistatic': 0.0,
                    'prediction_independent': 0.0,
                    'frequency': '',
                    'column_conservation': '',
                    'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
                })
            continue

        # Look up AA mutation
        aa_mutant = f"{wt_aa}{aa_pos}{mut_aa}"
        if aa_mutant in aa_lookup:
            r = aa_lookup[aa_mutant]
            results.append({
                'pkey': pkey,
                'mutant': aa_mutant,
                'pos': r['pos'],
                'wt': r['wt'],
                'subs': r['subs'],
                'prediction_epistatic': r['prediction_epistatic'],
                'prediction_independent': r['prediction_independent'],
                'frequency': r['frequency'],
                'column_conservation': r['column_conservation'],
                'qc_flags': 'PASS',
            })
        else:
            qc_flags.append('NOT_IN_MODEL')
            results.append({
                'pkey': pkey,
                'mutant': aa_mutant,
                'pos': aa_pos,
                'wt': wt_aa,
                'subs': mut_aa,
                'prediction_epistatic': '',
                'prediction_independent': '',
                'frequency': '',
                'column_conservation': '',
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
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


def write_matrix_output(results, output_path, alphabet="ACDEFGHIKLMNPQRSTVWY"):
    """
    Write results as position x substitution matrix.

    Args:
        results: List of result dictionaries
        output_path: Output file path
        alphabet: Amino acid alphabet (default: 20 standard AAs)
    """
    # Build matrix
    positions = sorted(set(r['pos'] for r in results if r['pos']))
    matrix = {pos: {aa: 0.0 for aa in alphabet} for pos in positions}

    for r in results:
        if r['pos'] and r['subs'] and r['prediction_epistatic']:
            try:
                matrix[r['pos']][r['subs']] = float(r['prediction_epistatic'])
            except (ValueError, KeyError):
                pass

    # Write as DataFrame
    df = pd.DataFrame(matrix).T
    df.index.name = 'position'
    df.to_csv(output_path, sep='\t')
    print(f"Wrote matrix to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: EVmutation/PLMC Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline with plmc
  python evmutation-pipeline.py \\
    --msa /path/to/protein.msa.fasta \\
    --focus GENE_HUMAN \\
    --plmc-binary /path/to/plmc \\
    --output /path/to/output.tsv

  # Use pre-computed model parameters
  python evmutation-pipeline.py \\
    --model-params /path/to/model.params \\
    --output /path/to/output.tsv

  # Map to specific nucleotide mutations
  python evmutation-pipeline.py \\
    --msa /path/to/protein.msa.fasta \\
    --focus GENE_HUMAN \\
    --plmc-binary /path/to/plmc \\
    --fasta /path/to/gene.fasta \\
    --mutations /path/to/mutations.csv \\
    --output /path/to/output.tsv

Notes:
  - MSA must be protein sequences (not nucleotides)
  - EVmutation predicts amino acid change effects
  - For synonymous mutations, prediction_epistatic = 0.0
  - column_conservation is derived from MSA single-site frequencies

References:
  - Hopf et al., Nature Biotechnology 2017
  - plmc: github.com/debbiemarkslab/plmc
"""
    )

    # MSA and model inputs
    msa_group = parser.add_mutually_exclusive_group(required=True)
    msa_group.add_argument('--msa', help='Path to protein MSA FASTA file')
    msa_group.add_argument('--model-params', help='Path to pre-computed plmc model parameters')

    parser.add_argument('--focus', help='Focus sequence identifier in MSA (required if --msa)')

    # plmc options
    parser.add_argument('--plmc-binary', help='Path to plmc binary executable')
    parser.add_argument('--alphabet', default='-ACDEFGHIKLMNPQRSTVWY',
                        help='Alphabet string (default: -ACDEFGHIKLMNPQRSTVWY)')
    parser.add_argument('--lambda-e', type=float, default=16.2,
                        help='J_ij regularization strength (default: 16.2)')
    parser.add_argument('--lambda-h', type=float, default=0.01,
                        help='h_i regularization strength (default: 0.01)')
    parser.add_argument('--skip-plmc', action='store_true',
                        help='Skip plmc, use existing model-params file')

    # Mutation mapping inputs (optional)
    parser.add_argument('--fasta', help='FASTA file with ORF sequence (for NT mutation mapping)')
    parser.add_argument('--mutations', help='Mutations CSV file (for NT mutation mapping)')
    parser.add_argument('--validation-log', help='Validation log for filtering')

    # Output options
    parser.add_argument('--output', '-o', required=True, help='Output TSV file path')
    parser.add_argument('--output-matrix', help='Output matrix file path (position x AA)')
    parser.add_argument('--quiet', '-q', action='store_true', help='Suppress verbose output')

    args = parser.parse_args()

    # Validate arguments
    if args.msa and not args.focus:
        parser.error("--focus required when using --msa")
    if args.msa and not args.skip_plmc and not args.plmc_binary:
        parser.error("--plmc-binary required when using --msa (or use --skip-plmc)")
    if args.fasta and not args.mutations:
        parser.error("--mutations required when using --fasta")

    # Determine model parameters path
    if args.model_params:
        model_params = args.model_params
    else:
        # Generate model params path from MSA
        msa_dir = os.path.dirname(args.msa)
        model_params = os.path.join(msa_dir, f"{args.focus}.model_params")

    # Run plmc if needed
    if args.msa and not args.skip_plmc:
        if not os.path.exists(model_params):
            print(f"Model params not found, running plmc...")
            run_plmc(
                fasta=args.msa,
                focus=args.focus,
                model_params=model_params,
                plmc_binary=args.plmc_binary,
                eij_lambda=args.lambda_e,
                hi_lambda=args.lambda_h,
                alphabet=args.alphabet,
                quiet=args.quiet
            )
        else:
            print(f"Using existing model params: {model_params}")

    # Load model
    if not os.path.exists(model_params):
        print(f"Error: Model params file not found: {model_params}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading CouplingsModel from {model_params}...")
    model = CouplingsModel(model_params)
    print(f"  Model length: {model.L}")
    print(f"  Alphabet: {''.join(model.alphabet)}")
    print(f"  N_eff: {model.N_eff}")

    # Determine gene name
    if args.fasta:
        gene = extract_gene_from_filename(args.fasta)
    elif args.focus:
        gene = args.focus.split('_')[0] if '_' in args.focus else args.focus
    else:
        gene = "GENE"

    # Generate all single mutant predictions
    print("Generating single mutant predictions...")
    all_results = generate_single_mutant_predictions(model, gene=gene)
    print(f"  Generated {len(all_results)} predictions")

    # Optionally filter to specific mutations
    if args.fasta and args.mutations:
        print(f"Mapping to nucleotide mutations from {args.mutations}...")
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

        results = map_nt_mutations_to_aa(mut_list, gene, orf_sequence, all_results, failure_map)
    else:
        results = all_results

    # Write output
    write_output(results, args.output)

    # Optionally write matrix output
    if args.output_matrix:
        write_matrix_output(all_results, args.output_matrix, alphabet=args.alphabet[1:])

    # Summary
    if results:
        n_with_pred = sum(1 for r in results if r.get('prediction_epistatic') not in ['', None])
        n_synonymous = sum(1 for r in results if 'SYNONYMOUS' in r.get('qc_flags', ''))
        n_not_in_model = sum(1 for r in results if 'NOT_IN_MODEL' in r.get('qc_flags', ''))

        print(f"\nSummary:")
        print(f"  Total mutations: {len(results)}")
        print(f"  With EVmutation prediction: {n_with_pred}")
        print(f"  Synonymous: {n_synonymous}")
        print(f"  Not in model: {n_not_in_model}")


if __name__ == "__main__":
    main()
