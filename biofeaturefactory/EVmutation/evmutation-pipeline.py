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

# Python 3.10+ compatibility: collections.Iterable moved to collections.abc
import collections
import collections.abc
if not hasattr(collections, 'Iterable'):
    collections.Iterable = collections.abc.Iterable

# EVmutation library (MIT licensed, Copyright Thomas A. Hopf)
# Clone from: https://github.com/debbiemarkslab/EVmutation
from EVmutation.model import CouplingsModel
import EVmutation.tools as ev_tools
from biofeaturefactory.utils.utility import (
    load_mapping,
    load_validation_failures,
    should_skip_mutation,
    extract_gene_from_filename,
)


FIELDNAMES = [
    'pkey',
    'mutant',
    'pos',
    'wt',
    'subs',
    'prediction_epistatic',
    'prediction_independent',
    'epistatic_contribution',
    'site_entropy',
    'mean_epistatic_at_pos',
    'std_epistatic_at_pos',
    'z_score_epistatic',
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


def compute_position_features(model):
    """
    Compute per-position features from single-site frequencies.

    Returns conservation (max freq excluding gaps) and Shannon entropy (bits).

    Args:
        model: CouplingsModel instance

    Returns:
        dict: {position: {'column_conservation': float, 'site_entropy': float}}
    """
    features = {}
    gap_index = model.alphabet_map.get('-', model.alphabet_map.get('.', 0))

    for i, pos in enumerate(model.index_list):
        freqs = model.f_i[i].copy()
        freqs[gap_index] = 0
        conservation = float(np.max(freqs))
        total = freqs.sum()
        if total > 0:
            probs = freqs / total
            entropy = -float(np.sum(probs[probs > 0] * np.log2(probs[probs > 0])))
        else:
            entropy = 0.0
        features[pos] = {'column_conservation': conservation, 'site_entropy': entropy}

    return features


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

    # Compute position features (conservation + entropy)
    pos_features = compute_position_features(model)

    # Build index for independent predictions: {(pos, subs): score}
    indep_index = {}
    for _, row in independent_df.iterrows():
        indep_index[(int(row['pos']), row['subs'])] = row['prediction_independent']

    # First pass: collect epistatic scores per position for summary stats
    pos_epistatic = {}
    for _, row in epistatic_df.iterrows():
        pos = int(row['pos'])
        pos_epistatic.setdefault(pos, []).append(float(row['prediction_epistatic']))

    pos_stats = {}
    for pos, scores in pos_epistatic.items():
        arr = np.array(scores)
        pos_stats[pos] = {'mean': float(np.mean(arr)), 'std': float(np.std(arr))}

    # Second pass: build result dicts with all fields
    for _, row in epistatic_df.iterrows():
        mutant = row['mutant']
        pos = int(row['pos'])
        wt = row['wt']
        subs = row['subs']

        indep_pred = indep_index.get((pos, subs), np.nan)
        epi_score = float(row['prediction_epistatic'])
        epistatic_contribution = epi_score - indep_pred if not np.isnan(indep_pred) else np.nan

        mean_epi = pos_stats[pos]['mean']
        std_epi = pos_stats[pos]['std']
        z_score = (epi_score - mean_epi) / std_epi if std_epi != 0 else np.nan

        pkey = f"{gene}-{mutant}" if gene else mutant

        pf = pos_features.get(pos, {})
        results.append({
            'pkey': pkey,
            'mutant': mutant,
            'pos': pos,
            'wt': wt,
            'subs': subs,
            'prediction_epistatic': epi_score,
            'prediction_independent': indep_pred,
            'epistatic_contribution': epistatic_contribution,
            'site_entropy': pf.get('site_entropy', np.nan),
            'mean_epistatic_at_pos': mean_epi,
            'std_epistatic_at_pos': std_epi,
            'z_score_epistatic': z_score,
            'frequency': row['frequency'],
            'column_conservation': pf.get('column_conservation', np.nan),
            'qc_flags': 'PASS',
        })

    return results



def map_mutations(mapping, gene, evmut_results, failure_map=None):
    """
    Rekey EVmutation results using a pre-computed nt→AA mapping file.

    Args:
        mapping: dict {ntposnt: aa_mutant_str} from load_mapping()
        gene: Gene symbol
        evmut_results: list of result dicts from generate_single_mutant_predictions
        failure_map: optional validation failure map

    Returns:
        list: result dicts with ntposnt-based pkeys
    """
    failure_map = failure_map or {}

    aa_lookup = {r['mutant']: r for r in evmut_results}

    pos_features = {}
    for r in evmut_results:
        if r['pos'] and r['pos'] not in pos_features:
            pos_features[r['pos']] = {
                'column_conservation': r['column_conservation'],
                'site_entropy': r['site_entropy'],
                'mean_epistatic_at_pos': r['mean_epistatic_at_pos'],
                'std_epistatic_at_pos': r['std_epistatic_at_pos'],
            }

    results = []
    aa_mutant_re = re.compile(r'^([A-Z*])(\d+)([A-Z*])$')

    for ntposnt, aa_mutant_str in mapping.items():
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        pkey = f"{gene}-{ntposnt}"
        qc_flags = []

        m = aa_mutant_re.match(aa_mutant_str)
        if not m:
            qc_flags.append('INVALID_MUTATION')
            results.append({
                'pkey': pkey, 'mutant': aa_mutant_str,
                'pos': '', 'wt': '', 'subs': '',
                'prediction_epistatic': '', 'prediction_independent': '',
                'epistatic_contribution': '', 'site_entropy': '',
                'mean_epistatic_at_pos': '', 'std_epistatic_at_pos': '',
                'z_score_epistatic': '', 'frequency': '', 'column_conservation': '',
                'qc_flags': ';'.join(qc_flags),
            })
            continue

        wt_aa, aa_pos, subs_aa = m.group(1), int(m.group(2)), m.group(3)

        if wt_aa == subs_aa:
            qc_flags.append('SYNONYMOUS')
            pf = pos_features.get(aa_pos, {})
            mean_epi = pf.get('mean_epistatic_at_pos', '')
            std_epi = pf.get('std_epistatic_at_pos', '')
            z_score = (0.0 - mean_epi) / std_epi if (std_epi not in ('', None) and std_epi != 0) else np.nan
            results.append({
                'pkey': pkey, 'mutant': aa_mutant_str,
                'pos': aa_pos, 'wt': wt_aa, 'subs': subs_aa,
                'prediction_epistatic': 0.0, 'prediction_independent': 0.0,
                'epistatic_contribution': 0.0,
                'site_entropy': pf.get('site_entropy', ''),
                'mean_epistatic_at_pos': mean_epi, 'std_epistatic_at_pos': std_epi,
                'z_score_epistatic': z_score, 'frequency': '',
                'column_conservation': pf.get('column_conservation', ''),
                'qc_flags': ';'.join(qc_flags),
            })
            continue

        if aa_mutant_str in aa_lookup:
            r = aa_lookup[aa_mutant_str]
            results.append({
                'pkey': pkey, 'mutant': aa_mutant_str,
                'pos': r['pos'], 'wt': r['wt'], 'subs': r['subs'],
                'prediction_epistatic': r['prediction_epistatic'],
                'prediction_independent': r['prediction_independent'],
                'epistatic_contribution': r['epistatic_contribution'],
                'site_entropy': r['site_entropy'],
                'mean_epistatic_at_pos': r['mean_epistatic_at_pos'],
                'std_epistatic_at_pos': r['std_epistatic_at_pos'],
                'z_score_epistatic': r['z_score_epistatic'],
                'frequency': r['frequency'],
                'column_conservation': r['column_conservation'],
                'qc_flags': 'PASS',
            })
        else:
            qc_flags.append('NOT_IN_MODEL')
            results.append({
                'pkey': pkey, 'mutant': aa_mutant_str,
                'pos': aa_pos, 'wt': wt_aa, 'subs': subs_aa,
                'prediction_epistatic': '', 'prediction_independent': '',
                'epistatic_contribution': '', 'site_entropy': '',
                'mean_epistatic_at_pos': '', 'std_epistatic_at_pos': '',
                'z_score_epistatic': '', 'frequency': '', 'column_conservation': '',
                'qc_flags': ';'.join(qc_flags),
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
    --gene SMN2 \\
    --plmc-binary /path/to/plmc \\
    --output /path/to/output.tsv

  # Use pre-computed model parameters, full AA matrix
  python evmutation-pipeline.py \\
    --model-params /path/to/model.params \\
    --gene SMN2 \\
    --output /path/to/output.tsv

  # Map to nt mutations using mapping file (ntposnt->aamutant)
  python evmutation-pipeline.py \\
    --model-params /path/to/model.params \\
    --gene SMN2 \\
    --map /path/to/combined_SMN2.csv \\
    --output /path/to/output.tsv

Notes:
  - MSA must be protein sequences (not nucleotides)
  - EVmutation predicts amino acid change effects
  - column_conservation is derived from MSA single-site frequencies
  - --map file must have 'mutant' (ntposnt) and 'aamutant' columns

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

    # Gene, mapping, and filtering
    parser.add_argument('--gene', help='Gene name for pkey generation (overrides focus-derived name)')
    parser.add_argument('--map', help='Mapping file (mutant->aamutant) for nt-based pkeys')
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
    if args.gene:
        gene = args.gene
    elif args.map:
        gene = extract_gene_from_filename(args.map)
    elif args.focus:
        gene = args.focus.split('_')[0] if '_' in args.focus else args.focus
    else:
        gene = "GENE"

    # Load validation failures if provided
    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}

    # Generate all single mutant predictions
    print("Generating single mutant predictions...")
    all_results = generate_single_mutant_predictions(model, gene=gene)
    print(f"  Generated {len(all_results)} predictions")

    # If mapping file provided, rekey to ntposnt-based pkeys
    if args.map:
        print(f"Loading mapping file: {args.map}")
        mapping = load_mapping(args.map, mapType='aamutant')
        if not mapping:
            print(f"Error: No entries loaded from mapping file {args.map}", file=sys.stderr)
            sys.exit(1)
        print(f"  Loaded {len(mapping)} mapping entries")
        results = map_mutations(mapping, gene, all_results, failure_map)
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
        if args.map:
            print(f"  Synonymous: {n_synonymous}")
            print(f"  Not in model: {n_not_in_model}")


if __name__ == "__main__":
    main()
