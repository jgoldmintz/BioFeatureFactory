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
    from codon_encoding import (CHAR_TO_CODON, CODON_ALPHABET, CODON_TO_CHAR,
                                encode_codon_msa)
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
    format_aa_token,
    load_validation_failures,
    mutation_class,
    parse_variant,
    protein_consequence,
    read_fasta,
    should_skip_mutation,
    splice_seq,
    trim_muts,
    write_tsv,
    split_intronic_tokens,
    warn_intronic_unsupported,
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
    """Run plmc to generate model parameters.

    -g (--gapignore) excludes the gap (first alphabet symbol) from the model:
    the standard EVmutation/plmc residue-level workflow. Applied to BOTH the AA
    and codon Potts models (both are residue-level alphabets where a symbol is
    an amino acid / a codon); it is NOT the nucleotide workflow.
    """
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
        # -g/--gapignore excludes the gap (first alphabet symbol) from the
        # model; no-argument flag, so fasta stays plmc's positional last arg.
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
    gap_index = model.alphabet_map.get("-", model.alphabet_map.get("."))
    for i, pos in enumerate(model.index_list):
        freqs = model.f_i[i].copy()
        if gap_index is not None:
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
    """Build lookup keyed by amino-acid mutant string (e.g., M1V).

    Returns (lookup, independent_model). The independent model is handed back
    rather than discarded because the non-SNV path scores multi-site changes by
    evaluating both Hamiltonians directly, and refitting it costs one
    fmin_bfgs per position.
    """
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
    return lookup, independent_model


def build_codon_lookup(codon_model):
    """
    Build lookup keyed by (codon_pos, mut_codon) for synonymous scoring.

    Returns (lookup, independent_model), for the same reason
    build_aa_prediction_lookup does: a synonymous MNV spanning two codons is a
    multi-site change, and a single-mutant lookup cannot answer it.

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
    return lookup, independent_model


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


# ---------------------------------------------------------------------------
# Non-SNV support
#
# A Potts model has a FIXED number of sites L and a FIXED symbol alphabet, so
# what a variant class costs is decided by whether the change can be written as
# "these sites now carry these symbols":
#
#   MNV / delins that shortens or preserves the protein  -> writable. It is a
#       MULTI-SITE substitution, which is exactly what delta_hamiltonian was
#       built for; a deleted residue is written as the gap symbol.
#   insertion (mut peptide longer than wt)               -> NOT writable. There
#       is no site to put the extra residue on and L cannot change.
#   frameshift                                           -> NOT writable. Every
#       downstream site changes identity and the product's length changes.
#
# Neither refusal is a defect to fix; both are named in qc_flags with a reason
# code so the row says WHY rather than disappearing. Whether the gap symbol is
# available is a property of the LOADED MODEL, not a user preference: plmc's
# -g/--gapreduce (run_plmc below always passes it) writes the params file
# without the gap (plmc/src/plm.c:638-640, plmc/src/inference.c:82-83), while
# params built without -g keep it. So it is read off model.alphabet_map at
# runtime rather than assumed, and there is no CLI flag for any of this.
# ---------------------------------------------------------------------------

# plmc writes '-' ; some alignments use '.'. compute_position_features above
# already probes both when it blanks the gap column, so match that.
_GAP_SYMBOLS = ("-", ".")

# protein_consequence's class -> the mutation_class column. The four classes the
# original vocabulary cannot express are carried through under their own names
# rather than being flattened into MISSENSE, which means "one residue swapped"
# and is what downstream filters select on. 'mnv' IS a substitution (no residue
# added or removed), so it stays MISSENSE; the multiplicity is visible in the
# comma-joined mutant/pos columns and named in qc_flags as AA:mnv.
_AA_CONSEQUENCE_TO_CLASS = {
    "synonymous": "SYNONYMOUS",
    "stop_gained": "STOP_GAIN",
    "stop_lost": "STOP_LOSS",
    "snv": "MISSENSE",
    "mnv": "MISSENSE",
    "inframe_del": "INFRAME_DEL",
    "inframe_ins": "INFRAME_INS",
    "inframe_delins": "INFRAME_DELINS",
    "frameshift": "FRAMESHIFT",
}


def _model_gap_symbol(model):
    """Return the gap symbol this model actually carries, or None."""
    for sym in _GAP_SYMBOLS:
        if sym in model.alphabet_map:
            return sym
    return None


def _delta_hamiltonian_multi(model, plan):
    """Delta statistical energy for SIMULTANEOUS substitutions at several sites.

    plan is [(model_pos, wt_symbol, mut_symbol), ...] in the model's own
    numbering. Returns np.array([full, couplings, fields]).

    This duplicates EVmutation.model._delta_hamiltonian on purpose. The vendored
    CouplingsModel.delta_hamiltonian allocates its index vectors with
    `np.int` (model.py:608-609), removed in numpy>=1.24, so on this stack it
    raises AttributeError before it reaches any model logic -- verified by
    execution, and EVmutation/ is a vendored read-only tree so it cannot be
    patched. The single-mutant path the SNV code uses is unaffected because it
    goes through the numba kernel (model.smm), which is why nothing noticed.

    Reproduces model.smm exactly for one substitution and model.dmm exactly for
    two (max |diff| 0.000e+00 over 200 random draws each), which is what pins
    this to the vendored arithmetic rather than merely resembling it.

    Raises KeyError for a position outside the model or a symbol outside its
    alphabet, ValueError when the model's own residue disagrees with the ORF's.
    """
    seq = model.target_seq_mapped
    J, h = model.J_ij, model.h_i
    L = h.shape[0]
    all_j = np.arange(L)

    idx, sym = [], []
    for pos, wt_symbol, mut_symbol in plan:
        i = int(model.index_map[pos])
        a = int(model.alphabet_map[mut_symbol])
        if str(model.target_seq[i]) != str(wt_symbol):
            raise ValueError(
                f"model position {pos} carries {model.target_seq[i]!s}, "
                f"ORF declares {wt_symbol!s}"
            )
        idx.append(i)
        sym.append(a)

    delta_h = 0.0
    delta_J = 0.0
    for m, (i, a) in enumerate(zip(idx, sym)):
        delta_h += float(h[i, a] - h[i, seq[i]])
        # couplings to every other site, in the wild-type background
        new = J[i, all_j, a, seq]
        old = J[i, all_j, seq[i], seq]
        delta_J += float(new.sum() - old.sum() - (new[i] - old[i]))  # drop j == i
        # correct the pairs where BOTH sites moved: remove the double count and
        # re-add the coupling in the new background (vendored comment at
        # model.py:159-176 explains the cancellation).
        for n in range(m + 1, len(idx)):
            j, b = idx[n], sym[n]
            delta_J -= float(J[i, j, a, seq[j]])
            delta_J -= float(J[i, j, seq[i], b])
            delta_J += float(J[i, j, seq[i], seq[j]])
            delta_J += float(J[i, j, a, b])

    return np.array([delta_J + delta_h, delta_J, delta_h])


def _position_stats(lookup):
    """Per-position MSA statistics, read back off an already-built lookup.

    site_entropy / column_conservation / mean & std of the position's single
    mutants are properties of the ALIGNMENT COLUMN, not of any one substitution,
    so they are defined for every variant that lands on a modelled position --
    including the ones whose score is refused. Taken from the lookup rather than
    recomputed so the two can never disagree.
    """
    stats = {}
    for rec in lookup.values():
        pos = rec["pos"]
        if pos not in stats:
            stats[pos] = {
                "site_entropy": rec["site_entropy"],
                "mean_epistatic_at_pos": rec["mean_epistatic_at_pos"],
                "std_epistatic_at_pos": rec["std_epistatic_at_pos"],
                "column_conservation": rec["column_conservation"],
            }
    return stats


def _join(values):
    """Comma-join per-site values.

    EVmutation's own multi-mutant convention: tools.extract_mutations splits on
    ',' and split_mutants documents "if higher-order mutations, symbols/numbers
    are comma-separated". A single-site change produces no comma, so an SNV row
    is unchanged.
    """
    return ",".join(str(v) for v in values)


def _substitution_plan(consequence, aa_pos, wt_aa, mut_aa):
    """Write a protein change as fixed-L site substitutions, or refuse by name.

    Returns (plan, refusal). Exactly one is meaningful: a non-empty plan with
    refusal None, or an empty plan with a reason code. That holds on a
    PRECONDITION the sole caller enforces -- wt_aa and mut_aa are never BOTH
    empty here, because protein_consequence classifies that case as 'synonymous'
    (utility.py:738) and _non_snv_rows routes synonymous to the codon table
    before reaching this function. Both empty would produce an empty plan with no
    refusal, and the caller would then score a no-op plan as a real 0.0. Anything
    that reorders those branches has to re-check this.

    The mutant residues are laid down from the 5' end of the changed span and
    any residues left over on the wild-type side become gaps. That is the same
    left-anchored reading _trim_aa_span already produced the span in, so a
    2-residue deletion written 'MK'->'' and one written 'MKL'->'M' land on the
    same sites.

    A deleted residue is written '-' regardless of what the loaded model turns
    out to carry: the plan describes the CHANGE, and whether a given model can
    score it is a separate question answered by the caller. That keeps one
    change with one mutant string whether or not it ends up scorable.
    """
    if consequence == "frameshift":
        # wt_aa/mut_aa are deliberately empty here (protein_consequence declines
        # to pair residues across a frame change), and every downstream site
        # changes identity, so there is nothing to write even in principle.
        return [], "FRAMESHIFT_NOT_REPRESENTABLE_FIXED_L"
    if len(mut_aa) > len(wt_aa):
        return [], "INSERTION_NOT_REPRESENTABLE_FIXED_L"

    plan = [(aa_pos + i, wt_aa[i], mut_aa[i]) for i in range(len(mut_aa))]
    plan += [(aa_pos + i, wt_aa[i], "-") for i in range(len(mut_aa), len(wt_aa))]
    return plan, None


def _wt_residue(orf_seq, aa_pos):
    """The wild-type residue at a 1-based codon position, '' if there is no codon."""
    codon = orf_seq[(aa_pos - 1) * 3:aa_pos * 3]
    if len(codon) != 3:
        return ""
    return aa_symbol(codon_to_aa.get(codon, "X"))


def _anchor_aa_alleles(aa_pos, wt_aa, mut_aa, orf_seq):
    """Re-anchor a one-sided aa change on a neighbouring residue.

    protein_consequence returns the MINIMAL span, so a pure insertion has an
    empty wild-type allele. An empty allele does not render as a token, and a
    consumer handed one drops the row -- the silent drop this path exists to
    prevent. Anchor 5' where there is a preceding residue, 3' at residue 1; the
    convention infer_aavariant_from_nt already uses for the mapping CSVs.
    """
    if wt_aa and mut_aa:
        return aa_pos, wt_aa, mut_aa
    if aa_pos > 1:
        anchor = _wt_residue(orf_seq, aa_pos - 1)
        if anchor:
            return aa_pos - 1, anchor + wt_aa, anchor + mut_aa
    anchor = _wt_residue(orf_seq, aa_pos + len(wt_aa))
    if anchor:
        return aa_pos, wt_aa + anchor, mut_aa + anchor
    return aa_pos, wt_aa, mut_aa


def _stats_columns(positions, pos_stats):
    """The four alignment-column statistics, comma-joined over the sites moved.

    A key with no value at ANY site stays a single empty cell rather than a run
    of bare commas, so "no model loaded" reads as empty instead of as data.
    """
    stats = [pos_stats.get(p, {}) for p in positions]
    columns = {}
    for key in ("site_entropy", "mean_epistatic_at_pos",
                "std_epistatic_at_pos", "column_conservation"):
        values = [s.get(key, "") for s in stats]
        columns[key] = _join(values) if any(v != "" for v in values) else ""
    return columns


def _score_substitution_plan(models, plan, pos_stats):
    """Score a substitution plan against the epistatic and independent models.

    Returns (columns, reason). `columns` is empty when the plan cannot be
    scored, and `reason` is the qc code in that case.
    """
    model, independent_model = models

    gap_symbol = _model_gap_symbol(model)
    if any(s == "-" for _, _, s in plan) and gap_symbol is None:
        # The headline refusal. plmc's -g strips the gap out of the params file
        # it writes, and run_plmc always passes -g, so a deletion is unscorable
        # against params this pipeline built -- but scorable against params
        # built without it. Read off the model, never assumed.
        return {}, "DELETION_NO_GAP_SYMBOL_IN_MODEL"
    if gap_symbol is not None and gap_symbol != "-":
        plan = [(p, w, gap_symbol if s == "-" else s) for p, w, s in plan]

    missing = [s for _, _, s in plan if s not in model.alphabet_map]
    if missing:
        return {}, f"SYMBOL_NOT_IN_MODEL:{','.join(sorted(set(missing)))}"
    if any(p not in model.index_map for p, _, _ in plan):
        return {}, "NOT_IN_MODEL"

    try:
        epistatic = float(_delta_hamiltonian_multi(model, plan)[0])
        independent = float(_delta_hamiltonian_multi(independent_model, plan)[0])
    except ValueError:
        # The model's focus sequence and this ORF disagree about the wild-type
        # residue. Naming it beats scoring a substitution the model never saw.
        return {}, "MODEL_WT_MISMATCH"

    columns = {
        "prediction_epistatic": epistatic,
        "prediction_independent": independent,
        "epistatic_contribution": epistatic - independent,
        "frequency": _join(float(model.fi(p, s)) for p, _, s in plan),
    }
    columns.update(_stats_columns([p for p, _, _ in plan], pos_stats))
    if len(plan) == 1:
        # z is against THIS column's single-substitution distribution. With more
        # than one site moving there is no such distribution to stand it
        # against, and inventing one would be the fabricated number the empty
        # cell exists to avoid.
        site = pos_stats.get(plan[0][0], {})
        mean, std = site.get("mean_epistatic_at_pos"), site.get("std_epistatic_at_pos")
        if mean is not None and std:
            columns["z_score_epistatic"] = (epistatic - mean) / std
    return columns, None


def _blank_row(fieldnames, pkey, nt_mut, qc_flags):
    row = {f: "" for f in fieldnames}
    row.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": ";".join(qc_flags)})
    return row


def _score_codon_plan(codon_lookup, codon_models, codon_plan):
    """Score whole-codon substitutions against the codon Potts model.

    codon_plan is [(codon_pos, wt_codon, mut_codon), ...]. Only length-preserving
    changes reach here: a codon-level row exists for synonymous variants, and the
    only non-SNV class that can be protein-synonymous is an MNV, which adds and
    removes nothing. So this never needs the gap symbol.

    One codon changed is the SNV case and goes through the precomputed lookup, so
    an MNV confined to a single codon scores identically to the SNV that produces
    the same codon. More than one codon changed is a multi-site evaluation, which
    the lookup (single mutants only) cannot answer.
    """
    if len(codon_plan) == 1:
        pos, _, mut_codon = codon_plan[0]
        if codon_lookup is None:
            return {}, "SYNONYMOUS_UNSCORED"
        scored = codon_lookup.get((pos, mut_codon))
        if scored is None:
            return {}, "SYNONYMOUS_NOT_IN_CODON_MODEL"
        return {
            "prediction_codon_epistatic": scored["prediction_codon_epistatic"],
            "prediction_codon_independent": scored["prediction_codon_independent"],
            "codon_epistatic_contribution": scored["codon_epistatic_contribution"],
            "codon_epistatic_concordance": scored["codon_epistatic_concordance"],
            "codon_frequency": scored["codon_frequency"],
        }, None

    if codon_models is None or not _CODON_ENCODING_AVAILABLE:
        return {}, "SYNONYMOUS_UNSCORED"

    model, independent_model = codon_models
    # codon_plan positions are SEQUENTIAL CDS codon positions -- the space
    # build_codon_lookup's KEYS live in, because it maps every alignment column to
    # its RANK in index_list before keying. The model's own APIs (index_map, fi,
    # target_seq) are keyed by the COLUMN, not the rank, and plmc numbers columns
    # from wherever the focus starts in the alignment (143..857, not 1..L). So the
    # position has to be carried back through that same map before the model sees
    # it: handing a rank straight to index_map silently evaluates a DIFFERENT site
    # whenever the two spaces disagree, and the wild-type check only catches it
    # when the two sites happen to carry different codons.
    if any(p < 1 or p > model.L for p, _, _ in codon_plan):
        return {}, "SYNONYMOUS_NOT_IN_CODON_MODEL"
    plan = [(int(model.index_list[pos - 1]),
             CODON_TO_CHAR.get(wt, ""), CODON_TO_CHAR.get(mut, ""))
            for pos, wt, mut in codon_plan]
    missing = [s for _, _, s in plan if s not in model.alphabet_map]
    if missing:
        return {}, "SYNONYMOUS_NOT_IN_CODON_MODEL"
    try:
        epistatic = float(_delta_hamiltonian_multi(model, plan)[0])
        independent = float(_delta_hamiltonian_multi(independent_model, plan)[0])
    except ValueError:
        return {}, "MODEL_WT_MISMATCH"
    return {
        "prediction_codon_epistatic": epistatic,
        "prediction_codon_independent": independent,
        "codon_epistatic_contribution": epistatic - independent,
        # concordance is a comparison against ONE position's contribution noise
        # floor (build_codon_lookup). Several positions moved, so that floor does
        # not exist; left empty and named rather than picked from one of them.
        "codon_epistatic_concordance": "",
        "codon_frequency": _join(float(model.fi(p, s)) for p, _, s in plan),
    }, "CONCORDANCE_UNDEFINED_MULTICODON"


def _non_snv_rows(pkey, nt_mut, variant, orf_seq, aa_pos_stats, aa_models,
                  codon_lookup, codon_models, skip_codon):
    """Build the row for one non-SNV token. Returns (protein_row, codon_row).

    Exactly one of the two is not None -- the routing mirrors the SNV path:
    synonymous and stop variants belong to the codon table unless --skip-codon
    sends them to the protein table.
    """
    qc = [f"NON_SNV:{variant.kind}"]

    # Contract: bound the END of the REF span, not its start. A multi-base REF
    # can begin inside the ORF and run off the end, which seq[pos] never catches.
    if variant.pos0 + len(variant.ref) > len(orf_seq):
        return _blank_row(PROTEIN_FIELDNAMES, pkey, nt_mut, qc + ["OUT_OF_RANGE"]), None

    first_codon = variant.pos0 // 3
    last_codon = (variant.pos0 + len(variant.ref) - 1) // 3
    # Bound the codon the span ENDS in, not just the one it starts in. A REF that
    # begins in a whole codon can finish inside a ragged final one, and checking
    # only first_codon lets that through as a 4- or 5-base "codon" in wt_codon.
    # last_codon >= first_codon, so this subsumes the start check it replaces.
    if last_codon * 3 + 3 > len(orf_seq):
        return _blank_row(PROTEIN_FIELDNAMES, pkey, nt_mut, qc + ["PARTIAL_CODON"]), None

    # Guard the WHOLE REF span, not just its first base: a multi-base REF whose
    # first base happens to match would otherwise pass. Flag and keep going,
    # which is what the SNV path does -- the ORF's own bases still describe a
    # real codon, and dropping the row would hide the disagreement.
    observed = orf_seq[variant.pos0:variant.pos0 + len(variant.ref)].upper()
    if observed != variant.ref.upper():
        qc.append("REF_MISMATCH")

    # No None check: protein_consequence returns None only for a span outside the
    # ORF, which the bound above has already rejected.
    cons = protein_consequence(variant, orf_seq)
    consequence = cons["aa_consequence"]
    aa_pos, wt_aa, mut_aa = cons["aa_pos"], cons["wt_aa"], cons["mut_aa"]
    qc.append(f"AA:{consequence}")
    if consequence == "frameshift" and cons["new_stop_aa_pos"] is not None:
        # The only informative number a frameshift row can carry: where
        # translation now terminates.
        qc.append(f"NEW_STOP_AA:{cons['new_stop_aa_pos']}")

    mut_orf = splice_seq(orf_seq, variant.pos0, variant.ref, variant.alt,
                         validate=False)
    wt_span = orf_seq[first_codon * 3:(last_codon + 1) * 3]
    delta = variant.length_delta
    # A frame-preserving edit shortens or lengthens the mutant span by exactly
    # its length delta. A frameshift does not: the codons after it are re-read
    # rather than removed, so the mutant span is the SAME width as the wild-type
    # one and reports what the ribosome now reads over those bases.
    mut_width = len(wt_span) + (delta if delta % 3 == 0 else 0)
    mut_span = mut_orf[first_codon * 3:first_codon * 3 + mut_width]

    shared = {
        "codon_position": first_codon + 1,
        "wt_codon": wt_span,
        "mut_codon": mut_span,
        "mutation_class": _AA_CONSEQUENCE_TO_CLASS.get(consequence, "UNKNOWN"),
    }

    # ---- codon-table classes: synonymous and stop, exactly as for an SNV ----
    if consequence in ("synonymous", "stop_gained", "stop_lost"):
        if consequence == "synonymous" and skip_codon:
            prow = _blank_row(PROTEIN_FIELDNAMES, pkey, nt_mut, qc)
            prow.update(shared)
            # h_i(M) - h_i(M) = 0 for an unchanged protein: a measured zero, not
            # a coalesced one, which is why it is written rather than blanked.
            prow.update({"prediction_epistatic": 0.0,
                         "prediction_independent": 0.0,
                         "epistatic_contribution": 0.0})
            qc.append("SYNONYMOUS_PROTEIN_LEVEL")
            prow["qc_flags"] = ";".join(qc)
            return prow, None
        if consequence != "synonymous" and skip_codon:
            prow = _blank_row(PROTEIN_FIELDNAMES, pkey, nt_mut, qc)
            prow.update(shared)
            qc.append(_AA_CONSEQUENCE_TO_CLASS[consequence])
            prow["qc_flags"] = ";".join(qc)
            return prow, None

        crow = _blank_row(CODON_FIELDNAMES, pkey, nt_mut, qc)
        crow.update({k: v for k, v in shared.items() if k in CODON_FIELDNAMES})
        if consequence != "synonymous":
            qc.append(_AA_CONSEQUENCE_TO_CLASS[consequence])
            crow["qc_flags"] = ";".join(qc)
            return None, crow

        codon_plan = [
            (c + 1, orf_seq[c * 3:c * 3 + 3], mut_orf[c * 3:c * 3 + 3])
            for c in range(first_codon, last_codon + 1)
            if orf_seq[c * 3:c * 3 + 3] != mut_orf[c * 3:c * 3 + 3]
        ]
        if not codon_plan:
            # Length-preserving and every codon identical: the token names no
            # change at all. Say so instead of emitting a row of zeros.
            qc.append("NO_CODON_CHANGED")
            crow["qc_flags"] = ";".join(qc)
            return None, crow
        columns, reason = _score_codon_plan(codon_lookup, codon_models, codon_plan)
        crow.update(columns)
        qc.append("SYNONYMOUS_SCORED" if columns else reason)
        if columns and reason:
            qc.append(reason)
        crow["qc_flags"] = ";".join(qc)
        return None, crow

    # ---- protein-table classes ----
    prow = _blank_row(PROTEIN_FIELDNAMES, pkey, nt_mut, qc)
    prow.update(shared)

    plan, refusal = _substitution_plan(consequence, aa_pos, wt_aa, mut_aa)
    if plan:
        # The model form: one 'M30V' per site moved, comma-joined, a deleted
        # residue written 'K61-'. Round-trips through the vendored
        # tools.extract_mutations, which is what reads mutant strings back.
        positions = [p for p, _, _ in plan]
        prow.update({"pos": _join(positions),
                     "wt": _join(w for _, w, _ in plan),
                     "subs": _join(s for _, _, s in plan),
                     "mutant": _join(f"{w}{p}{s}" for p, w, s in plan)})
    else:
        # No model form exists -- that IS the refusal. Fall back to the aa-level
        # token so the row still names the change. A frameshift has no residue
        # pair at all and takes the shared layer's 'K57fs' spelling; a pure
        # insertion has an empty wild-type allele and is anchored on the
        # preceding residue, the convention infer_aavariant_from_nt uses.
        if consequence == "frameshift":
            wt_res = _wt_residue(orf_seq, aa_pos)
            positions = [aa_pos]
            prow.update({"pos": aa_pos, "wt": wt_res, "subs": "",
                         "mutant": format_aa_token(aa_pos, wt_res, "", "frameshift")})
        else:
            token_pos, token_wt, token_mut = _anchor_aa_alleles(
                aa_pos, wt_aa, mut_aa, orf_seq)
            # The statistics below are read off the position this row NAMES, so
            # they have to follow the anchor. A pure insertion renders on the
            # PRECEDING residue, and leaving positions at aa_pos published the
            # next column's conservation and entropy against a pos that says
            # otherwise -- a measured number filed under the wrong site.
            positions = [token_pos]
            prow.update({"pos": token_pos, "wt": token_wt, "subs": token_mut,
                         "mutant": f"{token_wt}{token_pos}{token_mut}"})

    # The alignment column's statistics are properties of the POSITION, not of
    # the substitution, so they stay true even where no score exists. They are
    # what keeps a refused row informative instead of blank.
    prow.update(_stats_columns(positions, aa_pos_stats))

    if refusal:
        qc.append(refusal)
        prow["qc_flags"] = ";".join(qc)
        return prow, None
    if aa_models is None:
        qc.append("NOT_IN_MODEL")
        prow["qc_flags"] = ";".join(qc)
        return prow, None

    columns, reason = _score_substitution_plan(aa_models, plan, aa_pos_stats)
    if reason:
        qc.append(reason)
    else:
        prow.update(columns)
        if len(plan) > 1:
            qc.append("Z_SCORE_UNDEFINED_MULTISITE")
        qc.append("PASS")
    prow["qc_flags"] = ";".join(qc)
    return prow, None


def score_nt_mutations(nt_mutations, gene, orf_seq, aa_lookup, failure_map=None,
                       codon_lookup=None, skip_codon=False,
                       aa_models=None, codon_models=None):
    """
    Map and score NT mutations.

    Returns:
        tuple: (protein_rows, codon_rows)
          protein_rows  missense, stop, unknown, invalid
                        (also synonymous + stop when skip_codon=True)
          codon_rows    synonymous (empty when skip_codon=True)

    When skip_codon=True:
      synonymous variants land in the protein TSV with all scores set to 0
        (h_i(M) - h_i(M) = 0 by construction; flagged SYNONYMOUS_PROTEIN_LEVEL).
      stop_gain / stop_loss land in the protein TSV with empty scores
        (model lacks '*' symbol).

    aa_models / codon_models are (epistatic_model, independent_model) pairs, the
    two objects build_aa_prediction_lookup / build_codon_lookup already build.
    They are needed only by the non-SNV path, which evaluates the Hamiltonian
    directly: the precomputed lookups hold SINGLE mutants and skip the gap
    symbol outright (EVmutation/tools.py:149-151), so neither a multi-residue
    change nor a deletion can be answered from them.

    Non-SNV tokens are handled by DEFAULT -- there is no flag. Which columns
    survive is decided per column from the variant's own shape (a frameshift has
    no fixed-L representation; an MNV has one) and from what the loaded model
    carries (the gap symbol), both facts of the record rather than user
    preferences. A per-run boolean could only permit everything or refuse
    everything, and neither is what the data needs.
    """
    failure_map = failure_map or {}
    nt_re = re.compile(r"^([ACGT])(\d+)([ACGT])$")
    protein_rows = []
    codon_rows = []
    aa_pos_stats = _position_stats(aa_lookup)

    for nt_mut in nt_mutations:
        if should_skip_mutation(gene, nt_mut, failure_map):
            continue

        pkey = f"{gene}-{nt_mut}"
        qc_flags = []

        m = nt_re.match(nt_mut)
        if not m:
            # nt_re is uppercase-ACGT only, so a perfectly valid SNV written
            # lowercase or with U (RNA notation) missed it and fell through to
            # INVALID_MUTATION -- parse_variant accepts all three spellings and
            # reports is_snv=True for each, so the token was well formed and the
            # label was false. Canonicalise for the GATE only; pkey and
            # nt_mutant keep the user's original spelling, and an already-
            # uppercase token cannot reach here, so SNV output is byte-identical.
            canon = nt_mut.upper().replace("U", "T")
            m = nt_re.match(canon)
        if not m:
            # nt_re stays the SNV gate and stays FIRST, so every token it
            # accepts takes byte-identical the path it always did. Only tokens
            # it rejects can reach here, and of those only the ones that parse
            # as a genuine multi-base variant are treated differently: an
            # SNV-shaped token the SNV gate declines (lowercase, U instead of T)
            # still falls through to INVALID_MUTATION exactly as before.
            #
            # parse_variant is what makes "indel" distinguishable from
            # "garbage"; try/except around the legacy parser cannot, because
            # both raise the same ValueError.
            variant = parse_variant(nt_mut, is_nt=True)
            if variant is not None and not variant.is_snv:
                prow, crow = _non_snv_rows(
                    pkey, nt_mut, variant, orf_seq, aa_pos_stats, aa_models,
                    codon_lookup, codon_models, skip_codon)
                if prow is not None:
                    protein_rows.append(prow)
                else:
                    codon_rows.append(crow)
                continue
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
            if skip_codon:
                prow = {f: "" for f in PROTEIN_FIELDNAMES}
                prow.update({"pkey": pkey, "nt_mutant": nt_mut})
                prow.update(shared)
                prow["subs"] = aa_symbol(mut_aa)
                prow.update({
                    "prediction_epistatic":   0.0,
                    "prediction_independent": 0.0,
                    "epistatic_contribution": 0.0,
                })
                qc_flags.append("SYNONYMOUS_PROTEIN_LEVEL")
                prow["qc_flags"] = ";".join(qc_flags) if qc_flags else "SYNONYMOUS_PROTEIN_LEVEL"
                protein_rows.append(prow)
            else:
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
            if skip_codon:
                prow = {f: "" for f in PROTEIN_FIELDNAMES}
                prow.update({"pkey": pkey, "nt_mutant": nt_mut})
                prow.update(shared)
                prow["subs"] = aa_symbol(mut_aa)
                qc_flags.append(mclass)
                prow["qc_flags"] = ";".join(qc_flags)
                protein_rows.append(prow)
            else:
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
    write_tsv(results, output_path, PROTEIN_FIELDNAMES, extrasaction='ignore')
    print(f"  Wrote {len(results)} rows -> {output_path}")


def write_codon_output(results, output_path):
    write_tsv(results, output_path, CODON_FIELDNAMES, extrasaction='ignore')
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

    def _write(protein_rows, codon_rows):
        """Write both TSVs and return the counts tuple.

        ALWAYS called, including on the empty paths. Returning before the writes
        left the gene with NO FILE, which downstream is indistinguishable from a
        gene the pipeline never reached -- and leaves any stale file from a prior
        run standing in its place. write_tsv emits a header-only file for an
        empty list, which is what every other BFF pipeline produces.
        """
        gene_dir = Path(output_dir) / gene / "EVmutation"
        gene_dir.mkdir(parents=True, exist_ok=True)
        write_protein_output(protein_rows, str(gene_dir / f"{gene}.protein.tsv"))
        if not getattr(args, "skip_codon", False):
            write_codon_output(codon_rows, str(gene_dir / f"{gene}.codon.tsv"))
        n_pass = sum(1 for r in protein_rows if "PASS" in (r.get("qc_flags") or ""))
        n_syn = sum(1 for r in codon_rows if "SYNONYMOUS_SCORED" in (r.get("qc_flags") or ""))
        return len(protein_rows), len(codon_rows), n_pass, n_syn

    if not nt_mutations:
        print("  (no mutations)")
        return _write([], [])

    # Intronic gate, before any model is loaded or any row is built. A Potts
    # model has a fixed number of SITES -- residues in the protein model, codons
    # in the codon model -- and an intron occupies neither, so there is no index
    # to evaluate the Hamiltonian at.
    #
    # Unguarded these tokens do not crash: nt_re declines them, parse_variant
    # returns None, and they land in the protein TSV flagged INVALID_MUTATION.
    # That label is false -- 'gd.T5000C' is well formed in a coordinate space
    # this model has no sites for -- and the row's presence implies it was
    # evaluated and found unscorable rather than never being scoreable at all.
    nt_mutations, intronic = split_intronic_tokens(nt_mutations)
    warn_intronic_unsupported(
        'evmutation', gene, intronic,
        "A Potts model is indexed by residue or codon site; an intron occupies "
        "neither. Score these with RNAfold, miranda, genesplicer or AlphaFold3.")
    # A named row per excluded token, the same contract netNglyc, netphos,
    # netMHC and NetSurfP3 keep. The token is unscorable, not absent: dropping it
    # entirely leaves a hole at EVmutation for anyone joining the five pipelines
    # on pkey, and is indistinguishable from a mutation that was never submitted.
    # Every metric column stays EMPTY -- there is no site index to evaluate.
    intronic_rows = [
        _blank_row(PROTEIN_FIELDNAMES, f"{gene}-{tok}", tok,
                   ["NON_ORF_TOKEN:no_residue_or_codon_site_in_potts_model"])
        for tok in intronic
    ]

    if not nt_mutations:
        print("  (every mutation was intronic)")
        return _write(intronic_rows, [])

    print(f"  Loaded {len(nt_mutations)} NT mutations")

    # Protein model
    aa_lookup = {}
    aa_models = None
    if model_params_path:
        if os.path.exists(model_params_path):
            if not args.quiet:
                print(f"  Loading protein model: {model_params_path}")
            model = CouplingsModel(model_params_path)
            if not args.quiet:
                print(f"    N_eff: {model.N_eff:.1f}  L: {model.L}")
            aa_lookup, aa_independent = build_aa_prediction_lookup(model)
            aa_models = (model, aa_independent)
        else:
            print(f"  Warning: protein model params not found: {model_params_path}",
                  file=sys.stderr)

    # Codon model
    codon_lookup = None
    codon_models = None
    skip_codon = getattr(args, "skip_codon", False)
    if codon_model_params_path and not skip_codon:
        if not _CODON_ENCODING_AVAILABLE:
            print("  Warning: codon_encoding module not available; synonymous unscored",
                  file=sys.stderr)
        elif os.path.exists(codon_model_params_path):
            if not args.quiet:
                print(f"  Loading codon model: {codon_model_params_path}")
            codon_model = CouplingsModel(codon_model_params_path)
            if not args.quiet:
                print(f"    N_eff: {codon_model.N_eff:.1f}  L: {codon_model.L}")
            codon_lookup, codon_independent = build_codon_lookup(codon_model)
            codon_models = (codon_model, codon_independent)
        else:
            print(f"  Warning: codon model params not found: {codon_model_params_path}",
                  file=sys.stderr)

    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    protein_rows, codon_rows = score_nt_mutations(
        nt_mutations, gene, orf_seq, aa_lookup,
        failure_map=failure_map, codon_lookup=codon_lookup,
        skip_codon=getattr(args, "skip_codon", False),
        aa_models=aa_models, codon_models=codon_models,
    )

    protein_rows.extend(intronic_rows)
    n_protein, n_codon, n_pass, n_syn_scored = _write(protein_rows, codon_rows)

    if codon_lookup is not None:
        n_concordant = sum(1 for r in codon_rows if r.get("codon_epistatic_concordance") == "CONCORDANT")
        n_discordant = sum(1 for r in codon_rows if r.get("codon_epistatic_concordance") == "DISCORDANT")
        print(f"  Protein: {len(protein_rows)} rows ({n_pass} scored missense) | "
              f"Codon: {len(codon_rows)} rows ({n_syn_scored} scored, "
              f"{n_concordant} concordant, {n_discordant} discordant)")
    else:
        print(f"  Protein: {len(protein_rows)} rows ({n_pass} scored missense) | "
              f"Codon: {len(codon_rows)} rows (unscored)")

    return n_protein, n_codon, n_pass, n_syn_scored


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
  python evmutation_pipeline.py \\
      --fasta SMN2.fasta --mutations SMN2.csv \\
      --model-params SMN2.model_params \\
      --codon-model-params SMN2.codon_model_params \\
      --output results/

  # Single gene, build from MSAs
  python evmutation_pipeline.py \\
      --fasta SMN2.fasta --mutations SMN2.csv \\
      --msa SMN2.msa.a2m \\
      --codon-msa SMN2.codon.msa.fasta \\
      --plmc-binary /usr/local/bin/plmc \\
      --output results/

  # Multi-gene directory
  python evmutation_pipeline.py \\
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
    parser.add_argument("--skip-codon", action="store_true",
                        help="Skip codon-level scoring; route synonymous + stop variants to the protein TSV "
                             "(synonymous score is 0 by construction; stop has no AA-level score)")
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
