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
BioFeatureFactory: Codon-Aware EVmutation Pipeline

Single-gene mode:  --fasta GENE.fasta --mutations GENE.csv
Multi-gene mode:   --fasta /path/to/fastas/ --mutations /path/to/mutations/

Model params:
  - --model-params / --codon-model-params: direct file path, or directory
    containing {GENE}.model_params / {GENE}.codon_model_params.
  - --msa / --codon-msa: file or directory of per-gene MSA files; plmc is
    run to build params (skipped if params already exist).
  - Omit both to run without that model.

Outputs (per gene):
  <output>/<GENE>/EVmutation/<GENE>.protein.tsv   missense, stop, unclassified
  <output>/<GENE>/EVmutation/<GENE>.codon.tsv     synonymous mutations

Protein table: epistatic and independent scores are primary when
  N_eff >> 6,500 (q=20).

Codon table: prediction_codon_independent is the primary score.
  codon_epistatic_concordance (CONCORDANT / DISCORDANT / NEUTRAL) flags
  whether the regularisation-dominated epistatic score agrees.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

# Codon encoding utilities (bin/codon_encoding.py)
_BIN_DIR = Path(__file__).resolve().parent / "bin"
sys.path.insert(0, str(_BIN_DIR))
try:
    from codon_encoding import CHAR_TO_CODON, CODON_ALPHABET, encode_codon_msa
    _CODON_ENCODING_AVAILABLE = True
except ImportError:
    _CODON_ENCODING_AVAILABLE = False

import numpy as np

# Python 3.10+ compatibility
import collections
import collections.abc
if not hasattr(collections, "Iterable"):
    collections.Iterable = collections.abc.Iterable

from EVmutation.model import CouplingsModel
import EVmutation.tools as ev_tools
from biofeaturefactory.utils.utility import (
    codon_to_aa,
    discover_fasta_files,
    extract_gene_from_filename,
    load_validation_failures,
    read_fasta,
    should_skip_mutation,
    trim_muts,
)


# ---------------------------------------------------------------------------
# Output schemas
# ---------------------------------------------------------------------------

PROTEIN_FIELDNAMES = [
    "pkey",
    "nt_mutant",
    "codon_position",
    "wt_codon",
    "mut_codon",
    "mutant",
    "pos",
    "wt",
    "subs",
    "mutation_class",
    "prediction_epistatic",
    "prediction_independent",
    "epistatic_contribution",
    "site_entropy",
    "mean_epistatic_at_pos",
    "std_epistatic_at_pos",
    "z_score_epistatic",
    "frequency",
    "column_conservation",
    "qc_flags",
]

CODON_FIELDNAMES = [
    "pkey",
    "nt_mutant",
    "codon_position",
    "wt_codon",
    "mut_codon",
    "mutation_class",
    "prediction_codon_independent",
    "prediction_codon_epistatic",
    "codon_epistatic_contribution",
    "codon_epistatic_concordance",
    "codon_frequency",
    "qc_flags",
]


# ---------------------------------------------------------------------------
# plmc runner
# ---------------------------------------------------------------------------

def run_plmc(fasta, focus, model_params, plmc_binary, ec_file=None,
             eij_lambda=16.2, hi_lambda=0.01, iterations=500, stepsize=0.2,
             alphabet="-ACDEFGHIKLMNPQRSTVWY", quiet=False):
    """Run plmc to generate model parameters."""
    if not os.path.exists(plmc_binary):
        raise RuntimeError(f"PLMC binary not found at: {plmc_binary}")

    if ec_file is None:
        ec_file = os.path.join(os.path.dirname(model_params), focus + ".ec")

    if not quiet:
        print(f"  Running PLMC on MSA: {fasta}")
        print(f"    Focus: {focus}  ->  {model_params}")

    cmd = [
        plmc_binary,
        "-o", model_params,
        "-a", alphabet,
        "-c", ec_file,
        "-f", focus,
        "-le", str(eij_lambda),
        "-lh", str(hi_lambda),
        "-m", str(iterations),
        "-t", str(stepsize),
        "-g", fasta,
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        print(f"  PLMC failed (exit {proc.returncode})", file=sys.stderr)
        print(f"  stderr: {stderr.decode()}", file=sys.stderr)
        raise RuntimeError("PLMC did not complete successfully.")
    if not quiet:
        print("  PLMC complete.")
    return model_params


# ---------------------------------------------------------------------------
# Model building
# ---------------------------------------------------------------------------

def compute_position_features(model):
    """Compute per-position conservation and entropy from model frequencies."""
    features = {}
    gap_index = model.alphabet_map.get("-", model.alphabet_map.get(".", 0))
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
        features[pos] = {"column_conservation": conservation, "site_entropy": entropy}
    return features


def build_aa_prediction_lookup(model):
    """Build lookup keyed by amino-acid mutant string (e.g., M1V)."""
    epistatic_df = ev_tools.single_mutant_matrix(model, output_column="prediction_epistatic")
    independent_model = model.to_independent_model()
    independent_df = ev_tools.single_mutant_matrix(independent_model, output_column="prediction_independent")
    pos_features = compute_position_features(model)

    indep_index = {}
    for _, row in independent_df.iterrows():
        indep_index[(int(row["pos"]), row["subs"])] = row["prediction_independent"]

    pos_epistatic = {}
    for _, row in epistatic_df.iterrows():
        pos_epistatic.setdefault(int(row["pos"]), []).append(float(row["prediction_epistatic"]))

    pos_stats = {
        pos: {"mean": float(np.mean(s)), "std": float(np.std(s))}
        for pos, s in pos_epistatic.items()
    }

    lookup = {}
    for _, row in epistatic_df.iterrows():
        mutant = row["mutant"]
        pos = int(row["pos"])
        subs = row["subs"]
        indep_pred = indep_index.get((pos, subs), np.nan)
        epi_score = float(row["prediction_epistatic"])
        ec = epi_score - indep_pred if not np.isnan(indep_pred) else np.nan
        mean_epi = pos_stats[pos]["mean"]
        std_epi = pos_stats[pos]["std"]
        pf = pos_features.get(pos, {})
        lookup[mutant] = {
            "mutant": mutant, "pos": pos, "wt": row["wt"], "subs": subs,
            "prediction_epistatic": epi_score,
            "prediction_independent": indep_pred,
            "epistatic_contribution": ec,
            "site_entropy": pf.get("site_entropy", np.nan),
            "mean_epistatic_at_pos": mean_epi,
            "std_epistatic_at_pos": std_epi,
            "z_score_epistatic": (epi_score - mean_epi) / std_epi if std_epi != 0 else np.nan,
            "frequency": row["frequency"],
            "column_conservation": pf.get("column_conservation", np.nan),
        }
    return lookup


def build_codon_lookup(codon_model):
    """
    Build lookup keyed by (codon_pos, mut_codon) for synonymous scoring.

    codon_epistatic_concordance:
      CONCORDANT  - same sign and |contribution| > 0.5 * position std
      DISCORDANT  - opposite signs above the noise floor; flag for uncertainty
      NEUTRAL     - within noise floor or std==0; epistatic adds no information
    """
    epistatic_df = ev_tools.single_mutant_matrix(codon_model, output_column="prediction_codon_epistatic")
    independent_model = codon_model.to_independent_model()
    independent_df = ev_tools.single_mutant_matrix(independent_model, output_column="prediction_codon_independent")

    # plmc stores positions as alignment column indices (e.g. 143..857), not
    # sequential 1..L. Map each column index to its rank in index_list so that
    # the lookup keys match the sequential CDS positions (1..L) used by
    # score_nt_mutations: aa_pos = (codon_start // 3) + 1.
    col_to_seq = {int(col): (i + 1) for i, col in enumerate(codon_model.index_list)}

    indep_index = {}
    for _, row in independent_df.iterrows():
        indep_index[(int(row["pos"]), row["subs"])] = float(row["prediction_codon_independent"])

    raw = {}
    pos_contributions = {}
    for _, row in epistatic_df.iterrows():
        col_pos = int(row["pos"])
        seq_pos = col_to_seq[col_pos]   # sequential 1-based CDS codon position
        subs_char = row["subs"]
        mut_codon = CHAR_TO_CODON.get(subs_char, "???")
        epi = float(row["prediction_codon_epistatic"])
        indep = indep_index.get((col_pos, subs_char), np.nan)
        contrib = epi - indep if not np.isnan(indep) else np.nan
        raw[(seq_pos, mut_codon)] = {
            "prediction_codon_epistatic": epi,
            "prediction_codon_independent": indep,
            "codon_frequency": float(row["frequency"]),
            "codon_epistatic_contribution": contrib,
        }
        if not np.isnan(contrib):
            pos_contributions.setdefault(seq_pos, []).append(contrib)

    pos_contrib_std = {
        pos: float(np.std(c)) if len(c) > 1 else 0.0
        for pos, c in pos_contributions.items()
    }

    lookup = {}
    for (pos, mut_codon), data in raw.items():
        epi = data["prediction_codon_epistatic"]
        indep = data["prediction_codon_independent"]
        contrib = data["codon_epistatic_contribution"]
        std = pos_contrib_std.get(pos, 0.0)
        threshold = 0.5 * std
        if np.isnan(contrib) or std == 0.0 or abs(contrib) <= threshold:
            concordance = "NEUTRAL"
        elif (epi >= 0) == (indep >= 0):
            concordance = "CONCORDANT"
        else:
            concordance = "DISCORDANT"
        lookup[(pos, mut_codon)] = {
            "prediction_codon_epistatic": epi,
            "prediction_codon_independent": indep,
            "codon_epistatic_contribution": contrib,
            "codon_epistatic_concordance": concordance,
            "codon_frequency": data["codon_frequency"],
        }
    return lookup


# ---------------------------------------------------------------------------
# Mutation scoring
# ---------------------------------------------------------------------------

def read_orf_sequence(fasta_path):
    """
    Load ORF sequence from FASTA.

    Prefers the 'ORF' key (exon_aware_mapping.py convention).
    Falls back to the first record if 'ORF' is absent.
    """
    seqs = read_fasta(fasta_path)
    if not seqs:
        raise RuntimeError(f"No sequences found in FASTA: {fasta_path}")
    for sid, seq in seqs.items():
        if sid.upper() == "ORF":
            return sid, seq.upper()
    first_id, first_seq = next(iter(seqs.items()))
    return first_id, first_seq.upper()


def aa_symbol(aa):
    return "*" if aa == "Stop" else aa


def mutation_class(wt_aa, mut_aa):
    if wt_aa == "X" or mut_aa == "X":
        return "UNKNOWN"
    if wt_aa == mut_aa:
        return "SYNONYMOUS"
    if mut_aa == "Stop" and wt_aa != "Stop":
        return "STOP_GAIN"
    if wt_aa == "Stop" and mut_aa != "Stop":
        return "STOP_LOSS"
    return "MISSENSE"


def score_nt_mutations(nt_mutations, gene, orf_seq, aa_lookup, failure_map=None, codon_lookup=None):
    """
    Map and score NT mutations.

    Returns:
        tuple: (protein_rows, codon_rows)
          protein_rows  missense, stop, unknown, invalid
          codon_rows    synonymous
    """
    failure_map = failure_map or {}
    nt_re = re.compile(r"^([ACGT])(\d+)([ACGT])$")
    protein_rows = []
    codon_rows = []

    for nt_mut in nt_mutations:
        if should_skip_mutation(gene, nt_mut, failure_map):
            continue

        pkey = f"{gene}-{nt_mut}"
        qc_flags = []

        m = nt_re.match(nt_mut)
        if not m:
            prow = {f: "" for f in PROTEIN_FIELDNAMES}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "INVALID_MUTATION"})
            protein_rows.append(prow)
            continue

        ref_nt, pos_str, alt_nt = m.groups()
        nt_pos = int(pos_str)
        idx = nt_pos - 1
        if idx < 0 or idx >= len(orf_seq):
            prow = {f: "" for f in PROTEIN_FIELDNAMES}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "OUT_OF_RANGE"})
            protein_rows.append(prow)
            continue

        if orf_seq[idx] != ref_nt:
            qc_flags.append("REF_MISMATCH")

        codon_start = (idx // 3) * 3
        if codon_start + 3 > len(orf_seq):
            prow = {f: "" for f in PROTEIN_FIELDNAMES}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "PARTIAL_CODON"})
            protein_rows.append(prow)
            continue

        wt_codon = orf_seq[codon_start:codon_start + 3]
        mut_codon_list = list(wt_codon)
        mut_codon_list[idx % 3] = alt_nt
        mut_codon = "".join(mut_codon_list)

        wt_aa  = codon_to_aa.get(wt_codon, "X")
        mut_aa = codon_to_aa.get(mut_codon, "X")
        aa_pos = (codon_start // 3) + 1
        aa_mutant = f"{aa_symbol(wt_aa)}{aa_pos}{aa_symbol(mut_aa)}"
        mclass = mutation_class(wt_aa, mut_aa)

        shared = {
            "codon_position": aa_pos, "wt_codon": wt_codon, "mut_codon": mut_codon,
            "mutant": aa_mutant, "pos": aa_pos, "wt": aa_symbol(wt_aa),
            "mutation_class": mclass,
        }

        if mclass == "SYNONYMOUS":
            crow = {f: "" for f in CODON_FIELDNAMES}
            crow.update({"pkey": pkey, "nt_mutant": nt_mut})
            crow.update(shared)
            if codon_lookup is not None:
                scored = codon_lookup.get((aa_pos, mut_codon))
                if scored is None:
                    qc_flags.append("SYNONYMOUS_NOT_IN_CODON_MODEL")
                else:
                    crow.update({
                        "prediction_codon_epistatic":   scored["prediction_codon_epistatic"],
                        "prediction_codon_independent": scored["prediction_codon_independent"],
                        "codon_epistatic_contribution": scored["codon_epistatic_contribution"],
                        "codon_epistatic_concordance":  scored["codon_epistatic_concordance"],
                        "codon_frequency":              scored["codon_frequency"],
                    })
                    qc_flags.append("SYNONYMOUS_SCORED")
            else:
                qc_flags.append("SYNONYMOUS_UNSCORED")
            crow["qc_flags"] = ";".join(qc_flags)
            codon_rows.append(crow)

        elif mclass in {"STOP_GAIN", "STOP_LOSS"}:
            crow = {f: "" for f in CODON_FIELDNAMES}
            crow.update({"pkey": pkey, "nt_mutant": nt_mut})
            crow.update(shared)
            qc_flags.append(mclass)
            crow["qc_flags"] = ";".join(qc_flags)
            codon_rows.append(crow)

        else:
            prow = {f: "" for f in PROTEIN_FIELDNAMES}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut})
            prow.update(shared)
            prow["subs"] = aa_symbol(mut_aa)
            if mclass == "UNKNOWN":
                qc_flags.append("UNKNOWN_CODON")
            else:
                scored = aa_lookup.get(aa_mutant)
                if scored is None:
                    qc_flags.append("NOT_IN_MODEL")
                else:
                    prow.update({
                        "prediction_epistatic":      scored["prediction_epistatic"],
                        "prediction_independent":    scored["prediction_independent"],
                        "epistatic_contribution":    scored["epistatic_contribution"],
                        "site_entropy":              scored["site_entropy"],
                        "mean_epistatic_at_pos":     scored["mean_epistatic_at_pos"],
                        "std_epistatic_at_pos":      scored["std_epistatic_at_pos"],
                        "z_score_epistatic":         scored["z_score_epistatic"],
                        "frequency":                 scored["frequency"],
                        "column_conservation":       scored["column_conservation"],
                    })
                    qc_flags.append("PASS")
            prow["qc_flags"] = ";".join(qc_flags) if qc_flags else "PASS"
            protein_rows.append(prow)

    return protein_rows, codon_rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_protein_output(results, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PROTEIN_FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"  Wrote {len(results)} rows -> {output_path}")


def write_codon_output(results, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CODON_FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"  Wrote {len(results)} rows -> {output_path}")


# ---------------------------------------------------------------------------
# Discovery / resolution helpers
# ---------------------------------------------------------------------------

def _find_file_for_gene(gene, path_arg, glob_patterns):
    """
    Resolve a per-gene file from a direct path or directory.

    Direct file: returned as-is.
    Directory: returns the first file matching any glob pattern whose filename
               extracts to the given gene name.
    """
    if not path_arg:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        for pattern in glob_patterns:
            for f in sorted(p.glob(pattern)):
                if extract_gene_from_filename(f.name).upper() == gene.upper():
                    return str(f)
    return None


def _resolve_model_params(gene, path_arg, suffix):
    """
    Resolve a model params file from a direct path or directory.

    File: use directly.
    Directory: looks for {gene}{suffix} inside it.
    """
    if not path_arg:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        candidate = p / f"{gene}{suffix}"
        if candidate.exists():
            return str(candidate)
    return None


def _build_model_params(gene, msa_file, focus, output_path, args, codon=False):
    """
    Encode (codon only) and run plmc to build model params.

    --skip-plmc skips inference only; encoding always runs for codon MSAs.
    Returns the output_path (may not exist if --skip-plmc was set and no
    pre-existing params are found).
    """
    if codon:
        if not _CODON_ENCODING_AVAILABLE:
            raise RuntimeError("codon_encoding module not found in bin/")
        encoded = msa_file + ".encoded.fasta"
        if not args.quiet:
            print(f"  Encoding codon MSA: {msa_file} -> {encoded}")
        n_seqs = encode_codon_msa(msa_file, encoded)
        if not args.quiet:
            print(f"    Encoded {n_seqs} sequences")
        effective_msa = encoded
        alphabet = CODON_ALPHABET
    else:
        effective_msa = msa_file
        alphabet = args.alphabet

    if os.path.exists(output_path):
        if not args.quiet:
            label = "codon " if codon else ""
            print(f"  Using existing {label}model params: {output_path}")
        return output_path

    if args.skip_plmc:
        return output_path  # caller checks existence

    if not args.plmc_binary:
        raise RuntimeError("--plmc-binary required to run plmc")

    run_plmc(
        fasta=effective_msa,
        focus=focus,
        model_params=output_path,
        plmc_binary=args.plmc_binary,
        alphabet=alphabet,
        eij_lambda=args.lambda_e,
        hi_lambda=args.lambda_h,
        quiet=args.quiet,
    )
    return output_path


def _resolve_per_gene_models(gene, args):
    """
    Resolve (and optionally build) model params for one gene.

    Resolution order per model type:
      1. --msa / --codon-msa provided:
           - If --model-params / --codon-model-params is also given, use that path
             for the params file (so params and MSAs can live in different directories).
           - Otherwise derive params path from the MSA's directory.
           - Build params via plmc if the resolved path does not yet exist.
      2. --model-params / --codon-model-params alone -> direct file or {gene} in directory.
      3. Neither -> None (that model is skipped).

    Returns (model_params_path, codon_model_params_path).
    Both may be None.
    """
    # -- Protein model --
    model_params = None
    if args.msa:
        msa_file = _find_file_for_gene(
            gene, args.msa, ["*.a2m", "*.msa.fasta", "*.msa.fa", "*.fasta", "*.fa", "*.fas"]
        )
        if msa_file:
            focus = args.focus or gene
            if args.model_params:
                target = _resolve_model_params(gene, args.model_params, ".model_params") \
                         or str(Path(args.model_params) / f"{gene}.model_params")
            else:
                target = str(Path(msa_file).parent / f"{gene}.model_params")
            _build_model_params(gene, msa_file, focus, target, args, codon=False)
            model_params = target if os.path.exists(target) else None
        elif not args.quiet:
            print(f"  Warning: no protein MSA found for {gene} in {args.msa}")
    elif args.model_params:
        model_params = _resolve_model_params(gene, args.model_params, ".model_params")

    # -- Codon model --
    codon_model_params = None
    if args.codon_msa:
        codon_msa_file = _find_file_for_gene(
            gene, args.codon_msa,
            ["*.codon.msa.fasta", "*.codon.fasta", "*.codon.fa",
             "*.codon.msa.fa", "*.fasta", "*.fa", "*.fas"]
        )
        if codon_msa_file:
            codon_focus = args.codon_focus or "ORF"
            if args.codon_model_params:
                target_codon = _resolve_model_params(
                    gene, args.codon_model_params, ".codon_model_params"
                ) or str(Path(args.codon_model_params) / f"{gene}.codon_model_params")
            else:
                target_codon = str(Path(codon_msa_file).parent / f"{gene}.codon_model_params")
            _build_model_params(gene, codon_msa_file, codon_focus, target_codon, args, codon=True)
            codon_model_params = target_codon if os.path.exists(target_codon) else None
        elif not args.quiet:
            print(f"  Warning: no codon MSA found for {gene} in {args.codon_msa}")
    elif args.codon_model_params:
        codon_model_params = _resolve_model_params(
            gene, args.codon_model_params, ".codon_model_params"
        )

    return model_params, codon_model_params


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def _process_gene(gene, fasta_file, mutations_file, model_params_path,
                   codon_model_params_path, output_dir, args):
    """
    Score all mutations for one gene and write the two output tables.

    Returns (n_protein_rows, n_codon_rows, n_pass, n_syn_scored).
    """
    _, orf_seq = read_orf_sequence(fasta_file)

    nt_mutations = trim_muts(mutations_file, log=args.validation_log, gene_name=gene)
    if not nt_mutations:
        print("  (no mutations) -> skipping")
        return 0, 0, 0, 0
    print(f"  Loaded {len(nt_mutations)} NT mutations")

    # Protein model
    aa_lookup = {}
    if model_params_path:
        if os.path.exists(model_params_path):
            if not args.quiet:
                print(f"  Loading protein model: {model_params_path}")
            model = CouplingsModel(model_params_path)
            if not args.quiet:
                print(f"    N_eff: {model.N_eff:.1f}  L: {model.L}")
            aa_lookup = build_aa_prediction_lookup(model)
        else:
            print(f"  Warning: protein model params not found: {model_params_path}",
                  file=sys.stderr)

    # Codon model
    codon_lookup = None
    if codon_model_params_path:
        if not _CODON_ENCODING_AVAILABLE:
            print("  Warning: codon_encoding module not available; synonymous unscored",
                  file=sys.stderr)
        elif os.path.exists(codon_model_params_path):
            if not args.quiet:
                print(f"  Loading codon model: {codon_model_params_path}")
            codon_model = CouplingsModel(codon_model_params_path)
            if not args.quiet:
                print(f"    N_eff: {codon_model.N_eff:.1f}  L: {codon_model.L}")
            codon_lookup = build_codon_lookup(codon_model)
        else:
            print(f"  Warning: codon model params not found: {codon_model_params_path}",
                  file=sys.stderr)

    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    protein_rows, codon_rows = score_nt_mutations(
        nt_mutations, gene, orf_seq, aa_lookup,
        failure_map=failure_map, codon_lookup=codon_lookup,
    )

    gene_dir = Path(output_dir) / gene / "EVmutation"
    gene_dir.mkdir(parents=True, exist_ok=True)
    write_protein_output(protein_rows, str(gene_dir / f"{gene}.protein.tsv"))
    write_codon_output(codon_rows,   str(gene_dir / f"{gene}.codon.tsv"))

    n_pass       = sum(1 for r in protein_rows if "PASS" in (r.get("qc_flags") or ""))
    n_syn_scored = sum(1 for r in codon_rows   if "SYNONYMOUS_SCORED" in (r.get("qc_flags") or ""))

    if codon_lookup is not None:
        n_concordant = sum(1 for r in codon_rows if r.get("codon_epistatic_concordance") == "CONCORDANT")
        n_discordant = sum(1 for r in codon_rows if r.get("codon_epistatic_concordance") == "DISCORDANT")
        print(f"  Protein: {len(protein_rows)} rows ({n_pass} scored missense) | "
              f"Codon: {len(codon_rows)} rows ({n_syn_scored} scored, "
              f"{n_concordant} concordant, {n_discordant} discordant)")
    else:
        print(f"  Protein: {len(protein_rows)} rows ({n_pass} scored missense) | "
              f"Codon: {len(codon_rows)} rows (unscored)")

    return len(protein_rows), len(codon_rows), n_pass, n_syn_scored


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Codon-Aware EVmutation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output structure (per gene):
  <output>/<GENE>/EVmutation/<GENE>.protein.tsv
  <output>/<GENE>/EVmutation/<GENE>.codon.tsv

Model params resolution (--model-params / --codon-model-params):
  - File: used directly.
  - Directory: looks for {GENE}.model_params / {GENE}.codon_model_params inside it.
  - Omit: uses --msa / --codon-msa to build, or skips that model.

MSA resolution (--msa / --codon-msa):
  - File: used directly (single-gene mode).
  - Directory: finds {GENE}.* files by gene name (multi-gene mode).
  - --focus defaults to gene name; --codon-focus defaults to ORF.

Examples:
  # Single gene, pre-built params
  python evmutation-pipeline.py \\
      --fasta SMN2.fasta --mutations SMN2.csv \\
      --model-params SMN2.model_params \\
      --codon-model-params SMN2.codon_model_params \\
      --output results/

  # Single gene, build from MSAs
  python evmutation-pipeline.py \\
      --fasta SMN2.fasta --mutations SMN2.csv \\
      --msa SMN2.msa.a2m \\
      --codon-msa SMN2.codon.msa.fasta \\
      --plmc-binary /usr/local/bin/plmc \\
      --output results/

  # Multi-gene directory
  python evmutation-pipeline.py \\
      --fasta /data/fastas/ --mutations /data/mutations/ \\
      --msa /data/msas/ --codon-msa /data/codon_msas/ \\
      --plmc-binary /usr/local/bin/plmc \\
      --output results/
""",
    )

    parser.add_argument("--fasta", required=True,
                        help="ORF FASTA file (single gene) or directory (multi-gene)")
    parser.add_argument("--mutations", required=True,
                        help="Mutations CSV file or directory of per-gene CSVs")

    # Model params (file or directory; optional)
    parser.add_argument("--model-params",
                        help="Protein model params file or directory "
                             "containing {GENE}.model_params")
    parser.add_argument("--codon-model-params",
                        help="Codon model params file or directory "
                             "containing {GENE}.codon_model_params")

    # MSA / plmc args (file or directory)
    parser.add_argument("--msa",
                        help="Protein MSA file or directory of per-gene MSA files; "
                             "triggers plmc if model params absent")
    parser.add_argument("--focus",
                        help="Focus sequence ID in protein MSA "
                             "(default: gene name)")
    parser.add_argument("--codon-msa",
                        help="Codon MSA file or directory of per-gene codon MSA files; "
                             "encoded then run through plmc if params absent")
    parser.add_argument("--codon-focus",
                        help="Focus sequence ID in codon MSA (default: ORF)")

    parser.add_argument("--gene",
                        help="Gene name override (single-gene mode only)")
    parser.add_argument("--plmc-binary",
                        help="Path to plmc binary (required when running plmc)")
    parser.add_argument("--alphabet", default="-ACDEFGHIKLMNPQRSTVWY",
                        help="Protein alphabet for plmc (default: standard 20 AA)")
    parser.add_argument("--lambda-e", type=float, default=16.2,
                        help="J_ij regularisation strength (default: 16.2)")
    parser.add_argument("--lambda-h", type=float, default=0.01,
                        help="h_i regularisation strength (default: 0.01)")
    parser.add_argument("--skip-plmc", action="store_true",
                        help="Skip plmc inference (encoding still runs for codon MSAs)")
    parser.add_argument("--validation-log",
                        help="Validation log from exon_aware_mapping for mutation filtering")
    parser.add_argument("--output", "-o", default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress verbose output")

    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Directory / multi-gene mode
    # ------------------------------------------------------------------
    if fasta_path.is_dir():
        if args.gene:
            parser.error("--gene is only supported in single-gene mode")

        fasta_map = discover_fasta_files(str(fasta_path))
        if not fasta_map:
            print(f"No FASTA files found in {fasta_path}", file=sys.stderr)
            sys.exit(1)

        metrics = {
            "total": 0, "ok": 0, "fail": 0,
            "protein_rows": 0, "codon_rows": 0, "pass": 0, "syn_scored": 0,
        }
        failures = []

        for gene, fasta_file in sorted(fasta_map.items()):
            metrics["total"] += 1
            print(f"\n== {gene} ==")
            try:
                mutations_file = _find_file_for_gene(gene, args.mutations, ["*.csv"])
                if not mutations_file:
                    print("  No mutations file found -> skipping")
                    metrics["fail"] += 1
                    failures.append((gene, "No mutations file found"))
                    continue

                model_params, codon_model_params = _resolve_per_gene_models(gene, args)

                if not model_params and not codon_model_params:
                    print("  Warning: no model params resolved; output will contain "
                          "codon annotations only")

                n_p, n_c, n_pass, n_syn = _process_gene(
                    gene, fasta_file, mutations_file,
                    model_params, codon_model_params, output_dir, args,
                )
                metrics["ok"]          += 1
                metrics["protein_rows"] += n_p
                metrics["codon_rows"]   += n_c
                metrics["pass"]         += n_pass
                metrics["syn_scored"]   += n_syn

            except Exception as exc:
                print(f"  ERROR: {exc}", file=sys.stderr)
                metrics["fail"] += 1
                failures.append((gene, str(exc)))

        print(f"\nDone. Success: {metrics['ok']}/{metrics['total']}")
        if failures:
            print("Failures:")
            for name, msg in failures:
                print(f"  - {name}: {msg}")
        if metrics["total"] > 0:
            print("\nGrand totals:")
            print(f"  Genes:        {metrics['total']}  "
                  f"(passed: {metrics['ok']}, failed: {metrics['fail']})")
            print(f"  Protein rows: {metrics['protein_rows']}  "
                  f"(scored missense: {metrics['pass']})")
            print(f"  Codon rows:   {metrics['codon_rows']}  "
                  f"(synonymous scored: {metrics['syn_scored']})")
        return

    # ------------------------------------------------------------------
    # Single-gene mode
    # ------------------------------------------------------------------
    if not fasta_path.is_file():
        print(f"Error: --fasta path not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    gene = args.gene or extract_gene_from_filename(str(fasta_path)) or "GENE"
    print(f"\n== {gene} ==")

    mutations_file = _find_file_for_gene(gene, args.mutations, ["*.csv"])
    if not mutations_file:
        print(f"Error: mutations file not found (gene={gene}, arg={args.mutations})",
              file=sys.stderr)
        sys.exit(1)

    try:
        model_params, codon_model_params = _resolve_per_gene_models(gene, args)
    except RuntimeError as exc:
        parser.error(str(exc))

    if not model_params and not codon_model_params:
        print("Warning: no model params resolved; output will contain codon "
              "annotations only (no EVmutation scores)")

    _process_gene(gene, str(fasta_path), mutations_file,
                   model_params, codon_model_params, output_dir, args)


if __name__ == "__main__":
    main()
