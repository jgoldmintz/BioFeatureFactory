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
BioFeatureFactory: Rare Codon Enrichment Pipeline

Wrapper for cg_cotrans library (Copyright William M. Jacobs, GPL v3).
Detects regions enriched/depleted in rare codons using sliding window analysis.

Requires:
  - cg_cotrans library: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
  - Codon-aware MSA (FASTA), from core/codon_msa_pipeline.py

Optional, but the null model depends on it:
  - --reference-codon-usage: genome-wide codon usage TSV, built by
    scripts/build_db.sh at Bio_DBs/cocoputs/human_GRCh38_codon_usage.tsv.
    Without it "rare" is derived from the single gene under analysis, which makes
    the reference distribution and the test sequence the same data.
  - --usage: codon usage .p.gz; auto-built when absent.

Output format: TSV with BFF-standard pkey column
"""

import argparse
import csv
import gzip
import os
import pickle
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
CG_DIR = SCRIPT_DIR / "cg_cotrans"

from biofeaturefactory.lib.utility import (
    derive_mutations_root,
    discover_msa_files,
    read_fasta,
    mint_pkey,
    trim_muts,
    parse_variant,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
    write_tsv,
    split_intronic_tokens,
    warn_intronic_unsupported,
)

# Import cg_cotrans library (GPL v3 licensed, Copyright 2017 William M. Jacobs)
# Download from: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis

try:
    from calc_rare_enrichment import (
        read_fasta as rc_read_fasta,
        clean_sequences,
        sorted_gis,
        get_codons,
        align_sequences,
        aa_identity,
        load_null_model,
        msa_rare_codon_analysis_wtalign_nseq,
    )
    from codons import codon_to_aa
except ImportError:
    if not CG_DIR.exists():
        raise FileNotFoundError(f"cg_cotrans not found: {CG_DIR} you must download cg_cotrans from https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis")
    py_files = sorted(CG_DIR.glob("*.py"))
    if not py_files:
        raise FileNotFoundError(f"No *.py files found in {CG_DIR}")
      # copy into pipeline script directory (parent of cg_cotrans)
    import subprocess
    subprocess.run(
        ["cp", "-f", *[str(p) for p in py_files], str(SCRIPT_DIR)],
        check=True,
    )
    import importlib
    importlib.invalidate_caches()

    from calc_rare_enrichment import (
        read_fasta as rc_read_fasta,
        clean_sequences,
        sorted_gis,
        get_codons,
        align_sequences,
        aa_identity,
        load_null_model,
        msa_rare_codon_analysis_wtalign_nseq,
    )
    from codons import codon_to_aa


# ---------------------------------------------------------------------------
# cg_cotrans prob_ntuple: exponential -> O(N^2)
# ---------------------------------------------------------------------------
# calc_rare_enrichment.prob_ntuple(p, n) returns P(X >= n) for a Poisson-binomial
# X, by enumerating EVERY subset of `p` (itertools.combinations, two
# functools.reduce per subset). Cost is sum_m C(N, m) subsets at O(N) apiece.
#
# msa_rare_codon_analysis_wtalign_nseq calls it two ways. The per-window call
# (calc_rare_enrichment.py:210, :218) passes one entry per codon -- N ~ 15,
# 2^15 subsets, instant, and the only shape cg_cotrans was exercised at. The
# per-window-across-sequences call (:222, :224) passes one entry per MSA
# SEQUENCE, and that is the one that does not terminate:
#
#   N=35   worst case 2^34 = 17,179,869,184 subsets = 17.8 h for ONE call
#   N=709  worst case 2^708                         = not a runnable quantity
#
# A run against the 35-sequence SMN2 codon MSA (280 windows, 560 such calls)
# was killed after 19 h 21 m at 100% CPU having completed zero windows.
#
# BOTH branches of the original return the same quantity, P(X >= n); the
# `if n < len(p)/2` split only picks whichever half is cheaper to enumerate.
# So the standard Poisson-binomial recurrence below is a drop-in with identical
# semantics, not an approximation -- verified against the original over
# N in {1,2,3,5,8,12,16,18}, every n in 0..N+1, plus degenerate p=0 / p=1 /
# mixed: 267 comparisons, max absolute difference 4.8e-15 (float64 round-off).
#
#   N=35   51.9 us/call    (vs 17.8 h)
#   N=556   1.07 ms/call
#   N=709   1.41 ms/call
#
# Rebound onto the module rather than edited in place: calc_rare_enrichment.py
# is vendored cg_cotrans (GPL v3, Copyright 2017 William M. Jacobs) and is not
# ours to modify. msa_rare_codon_analysis_wtalign_nseq resolves `prob_ntuple`
# as a module global at CALL time, so assigning it here takes effect inside the
# vendored caller with the vendored file untouched.
import math as _math
import numpy as _np
import calc_rare_enrichment as _cre


def _prob_ntuple_poisson_binomial(p, n):
    """P(X >= n) where X ~ Poisson-binomial with per-trial probabilities `p`.

    dp[k] holds P(exactly k successes) over the trials consumed so far. Adding
    trial i moves probability mass up one slot with weight p_i and leaves it in
    place with weight 1 - p_i:

        dp'[k] = dp[k-1] * p_i + dp[k] * (1 - p_i)

    The slice assignment evaluates its whole right-hand side before storing, so
    every k >= 1 is updated from the PRE-update dp and no reverse iteration is
    needed. dp[0] is scaled afterwards because dp[:-1] above still had to read
    its old value.
    """
    N = len(p)
    if n <= 0:
        return 1.0
    if n > N:
        return 0.0
    dp = _np.zeros(N + 1)
    dp[0] = 1.0
    for pi in p:
        dp[1:] = dp[1:] * (1.0 - pi) + dp[:-1] * pi
        dp[0] *= (1.0 - pi)
    return float(dp[n:].sum())


_cre.prob_ntuple = _prob_ntuple_poisson_binomial


# ---------------------------------------------------------------------------
# cg_cotrans msa_rare_codon_analysis_wtalign_nseq: guard the empty window
# ---------------------------------------------------------------------------
# VERBATIM COPY of calc_rare_enrichment.py:158-249 with ONE behavioural change,
# marked <<GUARD>> below. Keep it a verbatim copy: if cg_cotrans is ever updated,
# re-copy the upstream body and re-apply the guard rather than merging by hand.
#
# Upstream divides by the number of codons a sequence has in the current window:
#
#     f_enriched[gi][center] = n_rare[gi][center] / len(all_indices)      (:207)
#
# len(all_indices) is 0 whenever a sequence is entirely gapped across the window.
# Measured: 105 of F9's 531 sequences have a gap run >= L, and one of them starts
# at window 0, so the run dies immediately and writes nothing. SMN2 has 2 of 35
# (both 85% coverage with the same 32-codon deletion -- exon-skip isoforms).
#
# Upstream INTENDED to handle this: :215 already guards `if len(all_indices) > 0`
# for nseq_possible_depleted_center, eight lines below the line that divides. So
# the fix is per-WINDOW exclusion, which is what that guard implies -- a sequence
# sits out the windows where it has no data and counts in the ones where it does.
#
# The alternative, dropping such sequences from the alignment entirely, costs far
# more: 19.8% of F9's (sequence, window) observations versus the 1.49% that are
# genuinely empty.
#
# This is a copy rather than a rebind of the inner arithmetic because the division
# is an inline statement, not a call that could be swapped -- unlike prob_ntuple.
# calc_rare_enrichment.py is vendored cg_cotrans (GPL v3, (C) 2017 William M.
# Jacobs) and is not ours to edit.

def _msa_rare_codon_analysis_guarded(msa_codons, wtgi, msa_index_dict,
                                     rare_codons, rare_codon_prob, L=10, zsig=1.,
                                     verbose=True):
    sorted_gis_ = _cre.sorted_gis
    prob_ntuple = _cre.prob_ntuple          # the O(N^2) recurrence bound above
    gis = sorted_gis_(msa_codons, wtgi)
    wt_ncodons = sum(1 for i in range(len(msa_codons[wtgi]))
                     if msa_codons[wtgi][i] != '---' and codon_to_aa[msa_codons[wtgi][i]] != 'Stop')
    f_enriched_avg = {gi: (sum(1 for i in range(wt_ncodons) for j in msa_index_dict[i]
                               if msa_codons[gi][j] != '---'
                               and msa_codons[gi][j] in rare_codons[gi]) /
                           sum(1 for i in range(wt_ncodons) for j in msa_index_dict[i]
                               if msa_codons[gi][j] != '---')) for gi in gis}
    f_gi_avg = _np.mean(list(f_enriched_avg.values()))

    def prob_poisson(l, n):
        return l ** n * _math.exp(-l) / _math.factorial(n)

    def min_n_poisson_cum(l, z):
        pz = _math.erf(z)
        s = 0
        n = -1
        while s < pz:
            n += 1
            s += prob_poisson(l, n)
        return n

    def max_n_poisson_cum(l, z):
        pz = _math.erfc(z)
        s = 0
        n = -1
        while s < pz:
            n += 1
            s += prob_poisson(l, n)
        return n

    p_enriched = {gi: {} for gi in gis}
    p_depleted = {gi: {} for gi in gis}
    f_enriched = {gi: {} for gi in gis}
    n_rare = {gi: {} for gi in gis}
    nseq_enriched = {}
    nseq_depleted = {}
    fseq_enriched = {}
    fseq_depleted = {}
    p_nseq_enriched = {}
    p_nseq_depleted = {}
    for i in range(wt_ncodons - L + 1):
        center = i + L // 2
        nseq_enriched[center] = nseq_depleted[center] = 0
        nseq_possible_enriched_center = nseq_possible_depleted_center = 0
        present = []                                            # <<GUARD>>
        for gi in gis:
            all_indices = sorted(k for j in range(i, i + L) for k in msa_index_dict[j]
                                 if msa_codons[gi][k] != '---')
            # <<GUARD>> This sequence has no codons at all in this window, so every
            # quantity below is undefined for it -- not zero. Leave n_rare and
            # f_enriched UNSET for (gi, center) so a consumer reading them with
            # .get() sees None and renders an empty cell, and leave gi out of the
            # cross-sequence tallies. Upstream instead divided by 0 here, and
            # counted the sequence as BOTH enriched and depleted when it didn't:
            # with len(all_indices) == 0 the thresholds collapse to
            # min_n_poisson_cum(0, z) == max_n_poisson_cum(0, z) == 0, so
            # `n_rare >= 0` and `n_rare <= 0` are both true.
            if not all_indices:
                continue
            present.append(gi)
            indices = [j for j in all_indices
                       if rare_codon_prob[gi][codon_to_aa[msa_codons[gi][j]]] > 0]
            n_rare[gi][center] = sum(1 for j in indices if msa_codons[gi][j] in rare_codons[gi])
            f_enriched[gi][center] = n_rare[gi][center] / len(all_indices)
            p_rc = _np.array([rare_codon_prob[gi][codon_to_aa[msa_codons[gi][j]]] for j in indices])
            nmin_enriched = min_n_poisson_cum(f_gi_avg * len(all_indices), zsig)
            p_enriched[gi][center] = prob_ntuple(p_rc, nmin_enriched)
            if n_rare[gi][center] >= nmin_enriched:
                nseq_enriched[center] += 1
            if len(p_rc) >= nmin_enriched:
                nseq_possible_enriched_center += 1
            nseq_possible_depleted_center += 1      # <<GUARD>> was `if len(all_indices) > 0`
            nmax_depleted = max_n_poisson_cum(f_gi_avg * len(all_indices), zsig)
            p_depleted[gi][center] = prob_ntuple(1. - p_rc, len(p_rc) - nmax_depleted)
            if n_rare[gi][center] <= nmax_depleted:
                nseq_depleted[center] += 1
        # <<GUARD>> aggregate over the sequences that had data, not over all gis:
        # an absent sequence has no p_enriched[gi][center] to read, and including
        # a placeholder would let "no data" vote in the cross-sequence test.
        p_enriched_center = _np.array([p_enriched[gi][center] for gi in present])
        p_depleted_center = _np.array([p_depleted[gi][center] for gi in present])
        if len(present):
            p_nseq_enriched[center] = prob_ntuple(p_enriched_center, nseq_enriched[center])
            p_nseq_depleted[center] = prob_ntuple(p_depleted_center, nseq_depleted[center])
        else:
            p_nseq_enriched[center] = p_nseq_depleted[center] = _np.nan
        if nseq_possible_enriched_center > 0:
            fseq_enriched[center] = nseq_enriched[center] / nseq_possible_enriched_center
        else:
            fseq_enriched[center] = _np.nan
        if nseq_possible_depleted_center > 0:
            fseq_depleted[center] = nseq_depleted[center] / nseq_possible_depleted_center
        else:
            fseq_depleted[center] = _np.nan
        if verbose:
            print("> %4d %2d %2d %6.3f %6.3f" %
                  (center, nseq_enriched[center], nseq_depleted[center],
                   _np.mean(list(f_enriched[g][center] for g in present)) if present else _np.nan,
                   f_gi_avg))
            print(' p =', ' '.join("%5.3f" % x for x in p_enriched_center))
    return {'nmin_enriched': min_n_poisson_cum(f_gi_avg * L, zsig),
            'nmax_depleted': max_n_poisson_cum(f_gi_avg * L, zsig),
            'p_enriched': p_enriched,
            'f_enriched': f_enriched,
            'f_enriched_avg': f_enriched_avg,
            'n_rare': n_rare,
            'nseq_enriched': nseq_enriched,
            'fseq_enriched': fseq_enriched,
            'p_nseq_enriched': p_nseq_enriched,
            'nseq_depleted': nseq_depleted,
            'fseq_depleted': fseq_depleted,
            'p_nseq_depleted': p_nseq_depleted}


# Shadow the name this module imported at the top, so run_rare_codon_analysis's
# call site picks up the guarded version without any change there.
msa_rare_codon_analysis_wtalign_nseq = _msa_rare_codon_analysis_guarded


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
                            null_model='genome', max_len_diff=0.2, min_aa_iden=0.5,
                            reference_usage=None):
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
        reference_usage: Path to a genome-wide codon usage TSV (see
            _reference_usage_pgz). When given, it defines "rare" for every sequence
            instead of the per-gene table auto-built from this MSA. None keeps the
            historical self-referential behaviour.

    Returns:
        dict: Position -> {p_enriched, p_depleted, f_enriched_wt, etc.}
    """
    # Load and clean MSA; sanitize ambiguous codons before clean_sequences
    seqs = rc_read_fasta(msa_path)
    # F54: codon_to_aa is uppercase-only (utility.py), so a soft-masked lowercase
    # codon ('atg') failed the membership test and was silently blanked to '---',
    # deleting real codons from the analysis; a fully lowercase WT collapsed to
    # wt_len=0. Upper-case before the test and store the upper-cased codon.
    seqs = {gi: ''.join(c.upper() if (c == '---' or c.upper() in codon_to_aa) else '---'
                        for c in (seq[i:i+3] for i in range(0, len(seq), 3)))
            for gi, seq in seqs.items()}
    seqs = clean_sequences(seqs)

    if wt_gi not in seqs:
        # F53: do NOT silently fall back to next(iter(seqs)). wt_gi anchors wt_len,
        # the overlap/identity filters, sorted_gis and aa_identity, so an arbitrary
        # homolog re-bases the entire enrichment analysis with no diagnostic.
        # The original raise was removed in e97c9e37 ("add directory input mode")
        # because one --wt-gi cannot match every gene and the raise aborted the run.
        # That is no longer a problem: _run_single_gene wraps this call in
        # try/except Exception and returns, so raising skips ONE gene loudly and
        # directory mode continues with the rest.
        raise ValueError(
            f"WT sequence '{wt_gi}' not found in MSA {msa_path} "
            f"({len(seqs)} records; first: {next(iter(seqs), 'none')}). "
            f"Pass --wt-gi matching the MSA header for this gene."
        )
    # Filter sequences by length
    wt_len = len(seqs[wt_gi]) - seqs[wt_gi].count('-')
    if wt_len == 0:
        # F56: every 3-mer was blanked to '---' above, i.e. this is not a codon MSA
        # (a protein MSA does exactly this). Previously surfaced as an opaque
        # ZeroDivisionError at the `gi_len / wt_len` length filter below.
        raise ValueError(
            f"{gene}: WT '{wt_gi}' has 0 valid codons after sanitization -- "
            f"{msa_path} does not look like a codon-aware nucleotide MSA "
            f"(a protein MSA produces this)."
        )
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
        val = aa_perc_id.get(gi, 0)
        # np.nan < threshold is False; use "not >=" to treat nan as failing
        if not (val >= min_aa_iden):
            del seqs[gi]

    # Rebuild after filtering
    gis = sorted_gis(seqs, wt_gi)
    msa_codons = {gi: get_codons(seq) for gi, seq in seqs.items()}

    # Remove sequences with no valid sense codons after position 30 (nstart used
    # by gene_avg_codon_probabilities) to prevent division by zero there.
    _NSTART = 30
    msa_codons = {gi: codons for gi, codons in msa_codons.items()
                  if gi == wt_gi or
                  sum(1 for c in codons[_NSTART:]
                      if c != '---' and c in codon_to_aa and codon_to_aa[c] != 'Stop') > 0}
    gis = sorted_gis(msa_codons, wt_gi)

    msa_index_dict = align_sequences(msa_codons, wt_gi)

    # Load null model.
    #
    # With a genome-wide reference the rare set is defined ONCE, from the reference,
    # and applied to every gi -- which is what cg_cotrans' use_wt_rare_codons flag
    # does. Without one we keep the historical per-gi behaviour so existing runs are
    # byte-identical.
    #
    # cl_wt_gi is set on the MODULE, not passed: calc_rare_enrichment.py:262 reads it
    # as a free name inside the use_wt_rare_codons branch, but it is not a parameter of
    # load_null_model and not a module-level global (verified: hasattr is False), so
    # that branch raises NameError on any input as shipped. It has never fired because
    # this call hardcoded the flag to False. Assigning the attribute here makes the
    # name resolve without editing the vendored file -- the same technique used for
    # the prob_ntuple rebind above.
    use_wt_rare_codons = False
    if reference_usage:
        usage_path = _reference_usage_pgz(
            reference_usage, gene, gis,
            Path(usage_path).parent / f"{gene}.reference_codon_usage.p.gz")
        use_wt_rare_codons = True
        _cre.cl_wt_gi = wt_gi

    rare_codons, gene_codon_usage, rare_codon_prob, _ = load_null_model(
        msa_codons, gis, usage_path, rare_model, use_wt_rare_codons,
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
            # No default. This is the direct structural sibling of
            # f_enriched_wt above -- same rc_analysis, same [wt_gi] indexing,
            # same key space -- and that one already returns None for an absent
            # position, which _rare_codon_row renders as empty. Defaulting this
            # one to 0 instead made the SAME missing key read as "measured, no
            # rare codons in this window", which is indistinguishable downstream
            # from a real zero. If cg_cotrans genuinely omits zero-count windows
            # then f_enriched_wt already loses them the same way, so the fix
            # belongs there and in both columns at once, not in a lone default
            # that hides the asymmetry.
            'n_rare': rc_analysis['n_rare'][wt_gi].get(pos),
        }

    return results, len(seqs)


def _rare_codon_row(gene, ntposnt, codon_position, rc_data, window_size, qc_flags):
    """Build the full 11-column row. Single builder for every outcome.

    The three exit paths of process_mutations used to construct three separate
    literals with the same keys, which is how a column comes to exist on one
    path and not another. They now share this one.

    rc_data is None when no enrichment window covers the codon. The six
    window-derived columns are then EMPTY, never 0: a 0 in p_enriched reads as
    "measured, maximally enriched" and a 0 in n_rare as "measured, no rare
    codons", and neither is distinguishable afterwards from a real null.
    qc_flags always names why.

    A per-key None inside an otherwise-present rc_data is emptied for the same
    reason AND named in qc_flags. Emptying it alone was still a silent loss: the
    row said 'PASS' while carrying blank cells, so a reader could not tell a
    column cg_cotrans never produced from one this pipeline failed to carry.
    """
    absent = []

    def cell(key):
        if rc_data is None:
            return ''
        value = rc_data.get(key)
        if value is None:
            absent.append(key)
            return ''
        return value

    row = {
        # {GENE}-{sha}. ntposnt comes straight from trim_muts, which strips only
        # '*' and whitespace, so it is the VERBATIM token variant_mapping hashed.
        'pkey': mint_pkey(gene, ntposnt),
        'Gene': gene,
        'codon_position': '' if codon_position is None else codon_position,
        'p_enriched': cell('p_enriched'),
        'p_depleted': cell('p_depleted'),
        'f_enriched_wt': cell('f_enriched_wt'),
        'frac_seq_enriched': cell('frac_seq_enriched'),
        'frac_seq_depleted': cell('frac_seq_depleted'),
        'n_rare': cell('n_rare'),
        'window_size': window_size,
        'qc_flags': '',      # filled below, once `absent` is complete
    }
    # Built after the dict literal because `absent` is only complete once every
    # cell() call has run.
    flags = list(qc_flags or [])
    if absent:
        flags.append(f"VALUE_ABSENT:{','.join(absent)}")
    row['qc_flags'] = ';'.join(flags) if flags else 'PASS'
    return row


def _variant_flags(variant):
    """qc_flags describing a non-SNV token, or [] for an SNV.

    Says what the row is and what it cannot say:
      NON_SNV:<kind>            substitution class from the shared Variant record
      length_delta:<+/-N>nt     omitted when 0 (an MNV is a non-SNV of equal length)
      FRAMESHIFT:...            same reason code codon_usage uses, because the
                                scope statement is the same: every downstream
                                codon is re-read, and one window row cannot
                                express that
      span_codons_<a>-<b>;centred_codon_<c>
                                only when the REF span crosses a codon boundary,
                                so the reader can see which codon the single
                                reported window actually describes
    """
    if variant.is_snv:
        return []
    flags = [f"NON_SNV:{variant.kind}"]
    if variant.length_delta:
        flags.append(f"length_delta:{variant.length_delta:+d}nt")
    if variant.length_delta % 3 != 0:
        flags.append('FRAMESHIFT:downstream_codons_also_change')
    first_codon = (variant.pos - 1) // 3 + 1
    last_codon = (variant.pos + len(variant.ref) - 2) // 3 + 1
    if last_codon != first_codon:
        flags.append(f"span_codons_{first_codon}-{last_codon}")
        flags.append(f"centred_codon_{_centre_codon(variant)}")
    return flags


def _centre_codon(variant):
    """1-based codon holding the MIDPOINT of the REF span.

    Anchoring on the first REF base instead would push the reported window
    steadily 5' of the region the variant actually disturbs as len(REF) grows:
    a 30 nt deletion would be described by the window centred 5 codons before
    its middle. For an SNV len(REF) is 1, the midpoint IS the first base, and
    this returns exactly what (nt_pos - 1) // 3 + 1 returned before.
    """
    centre_nt = variant.pos + (len(variant.ref) - 1) // 2
    return (centre_nt - 1) // 3 + 1


def process_mutations(mutations_list, gene, orf_sequence, rc_results, window_size, failure_map=None):
    """
    Annotate mutations with the rare codon enrichment of the window they fall in.

    SCOPE: the reported enrichment is a windowed WILD-TYPE property, derived from
    the WT sequence and its ortholog MSA. It is NOT allele-specific: only the
    codon position of the mutation is used, the mutant allele is not. Every
    mutation mapping to the same codon window therefore receives identical
    p_enriched/p_depleted/f_enriched_wt/n_rare values, and no WT-vs-MUT delta is
    computed or implied by these columns.

    Because the value is allele-independent, every column is defined for every
    variant class: a deletion, an insertion and a frameshift all sit at a codon
    position and that codon's window has an enrichment. Non-SNV tokens are
    therefore processed by default and fully populated; what changes is the
    qc_flags cell, which names the class and, for a multi-codon REF span, which
    codon the single reported window describes.

    Args:
        mutations_list: List of mutation strings
        gene: Gene symbol
        orf_sequence: ORF nucleotide sequence. Deliberately NOT used as a REF
            guard: the codon frame here is the MSA's WT record (wt_gi), and
            nothing in this pipeline establishes that the ORF FASTA and that MSA
            record are the same sequence rather than two isoforms. Guarding
            against the ORF would therefore flag correct rows whenever they
            differ. Out-of-range positions are already caught, in the right
            frame, by the rc_results lookup below.
        rc_results: Dict from run_rare_codon_analysis
        window_size: Window size used in analysis
        failure_map: Optional validation failure map

    Returns:
        list: Result dictionaries
    """
    results = []
    failure_map = failure_map or {}
    warned_missing = False

    for ntposnt in mutations_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        # parse_variant NEVER raises and is length-aware. The former
        # get_mutation_data_bioAccurate call did `int(ntposnt[1:-1])`, which on
        # any multi-base token (ACAA112A -> int("CAA112")) raised ValueError.
        # That exception escaped this function, escaped _run_single_gene, and
        # escaped main's gene loop: one indel token anywhere aborted the entire
        # run and lost every gene after it. Its F45 alphabet guard could not
        # catch it either, because it inspects only ntposnt[0] and ntposnt[-1],
        # which on an indel are both legal bases.
        variant = parse_variant(ntposnt, is_nt=True)
        if variant is None:
            results.append(_rare_codon_row(gene, ntposnt, None, None, window_size,
                                           ['INVALID_MUTATION']))
            continue

        qc_flags = _variant_flags(variant)
        codon_pos = _centre_codon(variant)

        # Look up enrichment data (center position of window)
        rc_data = rc_results.get(codon_pos)

        if rc_data is None:
            # Position may be near edges (not covered by window)
            if not warned_missing:
                keys = sorted(rc_results)
                observed = f"{keys[0]}..{keys[-1]}" if keys else "none"
                print(f"{gene}: codon_position {codon_pos} absent from rare codon "
                      f"results (observed keys {observed}, n={len(rc_results)}). "
                      f"Expected keys are 1-based window-centre codon positions; "
                      f"if most mutations miss, the cg_cotrans key convention "
                      f"(0- vs 1-based, centre vs start) does not match.",
                      file=sys.stderr)
                warned_missing = True
            qc_flags.append('POSITION_NOT_IN_WINDOW')

        results.append(_rare_codon_row(gene, ntposnt, codon_pos, rc_data,
                                       window_size, qc_flags))

    return results


def write_output(results, output_path):
    """Write results to TSV file."""
    if not results:
        print("No results to write")
        return

    write_tsv(results, output_path, FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(results)} rows to {output_path}")


def _resolve_per_gene(path_arg, gene, extensions):
    """Return a matching file for *gene* from a file or directory path.

    If path_arg is a file, return it.  If it is a directory, glob for
    ``{gene}.*`` files whose suffix (case-insensitive) is in *extensions*
    and return the first match, or None.
    """
    if path_arg is None:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        # Per-gene layout first. variant_mapping writes
        # <root>/<GENE>/mappings/mutations/<GENE>_mutations.csv, so a flat scan of
        # the root sees only directories and reports "Skipping <GENE>: missing
        # mutations" for a gene whose mutations file is right there. Search inside
        # <root>/<GENE>/ when that directory exists; the flat scan below still
        # handles a plain directory of per-gene files.
        gene_dir = next((c for c in sorted(p.iterdir())
                         if c.is_dir() and c.name.upper() == gene.upper()), None)
        if gene_dir is not None:
            exact = [f for f in sorted(gene_dir.rglob("*"))
                     if f.is_file() and f.suffix.lower() in extensions
                     and f.stem.upper() == gene.upper()]
            if exact:
                return str(exact[0])
            loose = [f for f in sorted(gene_dir.rglob("*"))
                     if f.is_file() and f.suffix.lower() in extensions
                     and gene.upper() in f.stem.upper()]
            if loose:
                return str(loose[0])

        for f in sorted(p.iterdir()):
            if f.stem.upper() == gene.upper() and f.suffix.lower() in extensions:
                return str(f)
        # broader search: gene anywhere in stem
        for f in sorted(p.iterdir()):
            if gene.upper() in f.stem.upper() and f.suffix.lower() in extensions:
                return str(f)
    return None


def _collect_msa_inputs(msa_arg):
    """Resolve one or more codon MSA files from --msa."""
    p = Path(msa_arg)
    if p.is_file():
        return [str(p.resolve())]
    if not p.is_dir():
        return []

    # Per-gene layout first. codon_msa_pipeline.py:1019 writes
    # <out>/<GENE>/CodonMSA/<GENE>.codon.msa.fasta, so a flat glob of the output
    # root finds nothing -- that root holds only gene directories. Measured
    # before this: `--msa out/` raised
    #     FileNotFoundError: Could not resolve MSA inputs from --msa: out/
    # The flat glob below still handles a plain directory of MSAs.
    disc = discover_msa_files(str(p), "codon")
    if disc:
        files = [Path(v) for _, v in sorted(disc.items())]
    else:
        patterns = ("*.msa.fasta", "*.fasta", "*.fa")
        files = []
        for pat in patterns:
            files.extend(sorted(p.glob(pat)))
    # stable de-dup preserving order
    seen = set()
    out = []
    for f in files:
        rf = str(f.resolve())
        if rf not in seen:
            seen.add(rf)
            out.append(rf)
    return out


def _ensure_codon_usage(args):
    """Return a usable codon usage .p.gz path, auto-building if missing."""
    if args.usage and Path(args.usage).is_file():
        return args.usage

    if args.usage and not Path(args.usage).exists():
        print(f"Usage file not found: {args.usage}; auto-building codon usage .p.gz", file=sys.stderr)
    elif not args.usage:
        print("--usage not provided; auto-building codon usage .p.gz")

    msa_inputs = _collect_msa_inputs(args.msa)
    if not msa_inputs:
        raise FileNotFoundError(f"Could not resolve MSA inputs from --msa: {args.msa}")

    usage_workdir = Path(args.output).resolve() / "_rare_codon_usage_cache"
    usage_workdir.mkdir(parents=True, exist_ok=True)

    abund_path = usage_workdir / "auto_abundances.tsv"
    # Empty abundances is valid; calc_codon_usage.py falls back to unit weights.
    abund_path.write_text("")

    calc_script = SCRIPT_DIR / "calc_codon_usage.py"
    if not calc_script.exists():
        raise FileNotFoundError(f"Required script missing: {calc_script}")

    cmd = [sys.executable, str(calc_script), *msa_inputs, str(abund_path)]
    # cg_cotrans implementation may require two passes before codon_usage.p.gz exists.
    try:
        for i in range(2):
            print(f"Building codon usage pass {i + 1}/2...")
            subprocess.run(
                cmd,
                cwd=str(usage_workdir),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        err_tail = (e.stderr or "").strip().splitlines()
        err_msg = err_tail[-1] if err_tail else str(e)
        print(
            "Auto-build via calc_codon_usage.py failed; "
            f"falling back to internal usage builder for sparse/synthetic MSAs. ({err_msg})",
            file=sys.stderr,
        )
        usage_path = usage_workdir / "codon_usage.p.gz"
        _build_usage_pgz_from_msa(msa_inputs, usage_path)
        print(f"Using codon usage: {usage_path}")
        return str(usage_path)

    usage_path = usage_workdir / "codon_usage.p.gz"
    if not usage_path.exists():
        # Defensive fallback in case calc script exits 0 but does not emit codon_usage.p.gz.
        _build_usage_pgz_from_msa(msa_inputs, usage_path)
    print(f"Using codon usage: {usage_path}")
    return str(usage_path)


def _build_usage_pgz_from_msa(msa_inputs, usage_path, pseudocount=1.0):
    """Build a minimal codon_usage.p.gz directly from codon MSA(s).

    This fallback is robust for very small/synthetic alignments where the original
    cg_cotrans builder can hit divide-by-zero for absent amino-acid classes.
    """
    aa_codons = defaultdict(list)
    sense_codons = []
    for codon, aa in codon_to_aa.items():
        if aa != "Stop":
            sense_codons.append(codon)
            aa_codons[aa].append(codon)

    # gi -> codon count table
    gi_counts = {}
    gene_groups = {}

    def _init_counts():
        d = defaultdict(float)
        for c in sense_codons:
            d[c] = float(pseudocount)
        return d

    def _sanitize_seqs(seqs):
        """Replace ambiguous codons (not in codon_to_aa) with '---' so that
        clean_sequences does not KeyError on IUPAC ambiguity codes (e.g. CCN)."""
        out = {}
        for gi, seq in seqs.items():
            sanitized = []
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3].upper()
                sanitized.append(codon if (codon == '---' or codon in codon_to_aa) else '---')
            out[gi] = ''.join(sanitized)
        return out

    for msa_path in msa_inputs:
        gene = extract_gene_from_filename(msa_path)
        gene_groups[gene] = 0
        seqs = clean_sequences(_sanitize_seqs(rc_read_fasta(msa_path)))
        for gi, seq in seqs.items():
            if gi not in gi_counts:
                gi_counts[gi] = _init_counts()
            if len(seq) % 3 != 0:
                raise ValueError(f"MSA sequence length must be multiple of 3 for {gi} in {msa_path}")
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                if codon == "---":
                    continue
                if codon in codon_to_aa and codon_to_aa[codon] != "Stop":
                    gi_counts[gi][codon] += 1.0

    if not gi_counts:
        raise ValueError("No usable GI sequences found while building fallback codon usage")

    relative_usage = {}
    for gi, counts in gi_counts.items():
        rel = {}
        for aa, codons in aa_codons.items():
            denom = sum(counts[c] for c in codons)
            if denom <= 0:
                # Should not happen due to pseudocount, but keep stable.
                denom = float(len(codons))
            for c in codons:
                rel[c] = counts[c] / denom
        relative_usage[gi] = rel

    usage_obj = {
        "groups": ["all"],
        "gene_groups": gene_groups,
        "overall_codon_usage": relative_usage,
        "unweighted_codon_usage": relative_usage,
        "gene_group_codon_usage": {gi: {0: relative_usage[gi]} for gi in relative_usage},
    }

    with gzip.open(usage_path, "wb") as f:
        pickle.dump(usage_obj, f)


def _reference_usage_pgz(ref_tsv, gene, gis, usage_path):
    """Build a codon_usage.p.gz from a genome-wide reference table.

    cg_cotrans' README specifies GENOME-WIDE codon statistics, built by running
    calc_codon_usage.py across many genes ("input MSA fasta file(s)", plural) before
    running calc_rare_enrichment.py once per gene. BFF never does that: the auto-build
    fails on every real input tried (KeyError 'gi|556503834|ref|NC_000913.3|' -- their
    default E. coli WT GI -- on SMN2, KeyError 'NNN' on F9), so _build_usage_pgz_from_msa
    takes over and derives "rare codon" from the ONE gene under analysis. That makes the
    reference distribution and the test sequence the same data, and it is wrong by a lot:
    against CoCoPUTs human genome-wide usage it produced 7 false rare codons and missed 2
    real ones out of a true set of 4, including calling GCC rare at 9.5% when it is 39.9%
    of all alanine codons in the human transcriptome.

    ref_tsv is the derived table at <Bio_DBs>/cocoputs/human_GRCh38_codon_usage.tsv,
    built by scripts/build_cocoputs_cut.py and installed by bootstrap step 9b alongside
    the other prepared databases (19,322 genes, 10,908,112 codon observations; see the
    README written beside it for derivation and validation). Its relative_usage_within_aa column uses the same within-amino-acid-family
    normalisation as calc_codon_usage.py:172, so the values drop straight into the slots
    load_null_model reads.

    The SAME table is assigned to every gi. cg_cotrans gives each organism its own rare
    set because a codon rare in E. coli is not rare in another bacterium; BFF scores human
    disease variants, so the question is whether a position is rare IN HUMANS, and the
    orthologs are evidence about that position rather than organisms to be scored on their
    own terms. The caller pairs this with use_wt_rare_codons=True so the rare set is
    likewise defined once.
    """
    rel = {}
    with open(ref_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rel[row["codon"]] = float(row["relative_usage_within_aa"])

    missing = [c for c, aa in codon_to_aa.items() if aa != "Stop" and c not in rel]
    if missing:
        raise ValueError(
            f"reference codon usage {ref_tsv} is missing {len(missing)} sense codons "
            f"({','.join(sorted(missing)[:6])}...); load_null_model indexes every sense "
            f"codon and would KeyError"
        )

    # One shared table per gi. load_null_model only reads these, never mutates them.
    relative_usage = {gi: rel for gi in gis}
    usage_obj = {
        "groups": ["all"],
        "gene_groups": {gene: 0},
        "overall_codon_usage": relative_usage,
        "unweighted_codon_usage": relative_usage,
        "gene_group_codon_usage": {gi: {0: rel} for gi in gis},
    }
    Path(usage_path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(usage_path, "wb") as f:
        pickle.dump(usage_obj, f)
    return str(usage_path)


def _resolve_wt_gi(seqs, wt_gi, gene):
    """Which MSA record is the WT/focus sequence.

    Precedence when --wt-gi is None: 'ORF', then the gene name, then the FIRST
    record. main passes args.wt_gi RAW (not `args.wt_gi or gene`) so that None
    still means "not given" here -- pre-substituting the gene name made an
    unrequested value appear in the error text and left this function unable to
    tell a default from a choice. codon_msa_pipeline names the focus record 'ORF', so defaulting
    to the gene name alone failed on every MSA it produces -- measured: SMN2
    MSA, 35 records, first 'ORF', lookup for 'SMN2' missed. The first-record
    fallback is last because an MSA whose focus row is neither named is still
    conventionally focus-first.

    An EXPLICIT --wt-gi is never overridden: if the caller named a record, a
    silent substitution would score a different sequence than they asked for.
    """
    if wt_gi:
        # EXPLICIT --wt-gi. Honour it or fail; never substitute another record.
        if wt_gi in seqs:
            return wt_gi, None
        return None, (f"WT record '{wt_gi}' not in MSA ({len(seqs)} records; "
                      f"first: {next(iter(seqs), 'none')}).")
    for cand in ('ORF', gene):
        if cand in seqs:
            return cand, None
    if seqs:
        return next(iter(seqs)), None
    return None, "MSA is empty."


def _orf_from_msa(msa_file, wt_gi, gene):
    """The WT ORF, degapped from the MSA's own WT record.

    Read RAW, never through _sanitize_seqs/clean_sequences. Those two rewrite the
    WT row -- _sanitize_seqs blanks any codon outside codon_to_aa to '---', and
    clean_sequences truncates from the first stop codon to the end -- so the
    sanitized row is a TRUNCATED, hole-punched ORF. Validating a mutation token
    against it would silently reject every position past the first stop.

    Returns (orf, resolved_wt_gi, error_message); error is None on success.
    """
    seqs = rc_read_fasta(msa_file)
    resolved, err = _resolve_wt_gi(seqs, wt_gi, gene)
    if err:
        return None, None, f"{err} In {msa_file}. Pass --wt-gi matching the header for {gene}."
    orf = seqs[resolved].replace('-', '').upper()
    if not orf:
        return None, None, f"WT record '{resolved}' in {msa_file} is all gaps; no ORF to score."
    return orf, resolved, None


def _run_single_gene(gene, msa_file, mut_file, args, wt_gi, output_dir):
    """Run rare codon analysis for a single gene and write per-gene output."""
    # wt_gi is REBOUND to whatever _orf_from_msa actually used: the ORF and the
    # enrichment analysis below must index the same MSA row, or the codon
    # positions the tokens are validated against are not the ones scored.
    orf_sequence, wt_gi, orf_err = _orf_from_msa(msa_file, wt_gi, gene)
    if orf_err:
        print(f"Error: {orf_err}", file=sys.stderr)
        return
    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    mut_list = trim_muts(mut_file, log=args.validation_log, gene_name=gene)

    if not mut_list:
        print(f"Error: No mutations found in {mut_file}", file=sys.stderr)
        return

    # Intronic gate. Every column here is indexed by CODON position: _centre_codon
    # converts an nt position to a codon, and rc_results is keyed by window-centre
    # codon of the MSA's WT record. An intron has no codon, so there is no key to
    # look up and no defensible value to report.
    #
    # Unguarded these tokens do not crash -- parse_variant returns None -- but
    # they produce a ROW flagged 'INVALID_MUTATION' with every window column
    # empty. That label is false: 'gd.T5000C' is a well-formed token in a
    # coordinate space this pipeline cannot index, and the row's presence in the
    # TSV implies it was analysed and found null. Exclude and say so instead.
    mut_list, intronic = split_intronic_tokens(mut_list)
    warn_intronic_unsupported(
        'rare_codon', gene, intronic,
        "Rare codon enrichment is indexed by codon position; an intron has none. "
        "Score these with RNAfold or miranda instead.")

    if not mut_list:
        print(f"Error: {gene}: every mutation was intronic; nothing to analyse",
              file=sys.stderr)
        return

    print(f"Running rare codon enrichment analysis for {gene}...")
    print(f"  MSA: {msa_file}")
    print(f"  Usage: {args.usage}")
    print(f"  Window size: {args.window_size}")

    try:
        rc_results, n_seqs = run_rare_codon_analysis(
            gene=gene,
            msa_path=msa_file,
            usage_path=args.usage,
            wt_gi=wt_gi,
            window_size=args.window_size,
            rare_model=args.rare_model,
            rare_threshold=args.rare_threshold,
            null_model=args.null_model,
            max_len_diff=args.max_len_diff,
            min_aa_iden=args.min_aa_iden,
            reference_usage=args.reference_codon_usage,
        )
        print(f"  MSA sequences used: {n_seqs}")
        print(f"  Positions analyzed: {len(rc_results)}")
    except Exception as e:
        print(f"Error in rare codon analysis for {gene}: {e}", file=sys.stderr)
        return

    results = process_mutations(
        mut_list, gene, orf_sequence, rc_results,
        args.window_size, failure_map
    )

    out_dir = Path(output_dir) / gene / "RareCodon"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_output(results, str(out_dir / f"{gene}.rare_codon.tsv"))

    if results:
        n_enriched = sum(1 for r in results
                         if r.get('p_enriched') not in (None, '')
                         and float(r['p_enriched']) < 0.05)
        n_depleted = sum(1 for r in results
                         if r.get('p_depleted') not in (None, '')
                         and float(r['p_depleted']) < 0.05)
        print(f"\nSummary for {gene}:")
        print(f"  Total mutations: {len(results)}")
        print(f"  In enriched regions (p<0.05): {n_enriched}")
        print(f"  In depleted regions (p<0.05): {n_depleted}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Rare Codon Enrichment Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single gene
  python rare_codon_pipeline.py \\
    --msa /path/to/codon_msa.fasta \\
    --usage /path/to/codon_usage.p.gz \\
    --wt-gi "focus_sequence_id" \\
    --fasta /path/to/gene.fasta \\
    --mutations /path/to/mutations.csv \\
    --output results/

  # Directory mode
  python rare_codon_pipeline.py \\
    --fasta fastas/ --msa msas/ --usage codon_usage.p.gz \\
    --mutations mutations/ --output results/

Required preprocessing:
  1. Generate codon-aware MSA (use msa/codon_msa_pipeline.py)
  2. Optional: provide --usage codon_usage.p.gz (auto-generated if omitted)

Copyright notice:
  cg_cotrans library is Copyright (C) 2017 William M. Jacobs (GPL v3)
  Source: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
"""
    )

    # MSA and codon usage inputs
    parser.add_argument('-a', '--msa', required=True,
                        help='Path to codon-aware MSA FASTA file or directory')
    parser.add_argument('-u', '--usage', help='Path to codon usage .p.gz file (optional; auto-built if missing)')
    parser.add_argument('-wg', '--wt-gi', help='Identifier for WT/focus sequence in MSA (single-gene mode)')

    # Mutation inputs
    # --fasta was removed: its only informational job was to supply the WT ORF,
    # and the MSA already carries it as the wt_gi record. Everything else the flag
    # did (enumerating genes in directory mode) the MSA resolver does too.
    parser.add_argument('-m', '--mutations', required=False,
                        help='Mutations CSV file or directory of CSV files')
    parser.add_argument('-vl', '--validation-log', help='Validation log for filtering')

    # Analysis parameters
    parser.add_argument('--window-size', '-L', type=int, default=15,
                        help='Sliding window width in codons (default: 15)')
    parser.add_argument('-rm', '--rare-model', choices=['no_norm', 'cmax_norm'], default='no_norm',
                        help='Rare codon definition model (default: no_norm)')
    parser.add_argument('-rt', '--rare-threshold', type=float, default=0.1,
                        help='Threshold for codon rarity (default: 0.1)')
    parser.add_argument('-nm', '--null-model', choices=['genome', 'eq', 'groups'], default='genome',
                        help='Null model for codon usage (default: genome)')
    parser.add_argument('-rcu', '--reference-codon-usage',
                        help='Genome-wide codon usage TSV defining which codons are rare '
                             '(columns: codon, amino_acid, count, relative_usage_within_aa). '
                             'Without it, "rare" is derived from the single gene under '
                             'analysis, which makes the null model self-referential -- see '
                             'the README beside the table in <Bio_DBs>/cocoputs/ for the '
                             'measured error.')
    parser.add_argument('--max-len-diff', type=float, default=0.2,
                        help='Max relative length difference vs WT (default: 0.2)')
    parser.add_argument('--min-aa-iden', type=float, default=0.5,
                        help='Min AA identity vs WT (default: 0.5)')

    # Output
    parser.add_argument('--output', '-o', required=True, help='Output base directory')

    args = parser.parse_args()

    # Directory mode: <root>/<GENE>/mappings/mutations/ sits beside the input,

    # so the root supplies both. Explicit --mutations always wins; FILE MODE and

    # any layout outside the tree are unaffected.

    args.mutations = derive_mutations_root(args.mutations, args.msa, "rare_codon")

    if not args.mutations:

        parser.error("--mutations is required (no <GENE>/mappings/mutations/ under "

                     f"--msa {args.msa})")

    args.usage = _ensure_codon_usage(args)

    msa_path = Path(args.msa)
    if msa_path.is_dir():
        # The MSA is now the gene enumerator. _collect_msa_inputs already globs
        # *.msa.fasta / *.fasta / *.fa and de-dups, which is exactly the listing
        # the FASTA glob used to provide.
        msa_files = _collect_msa_inputs(args.msa)
        if not msa_files:
            print(f"Error: No MSA files found in {args.msa}", file=sys.stderr)
            sys.exit(1)
        failed_genes = []
        for msa_file in msa_files:
            gene = extract_gene_from_filename(str(msa_file))
            mut_file = _resolve_per_gene(args.mutations, gene, ('.csv',))
            if not mut_file:
                print(f"  Skipping {gene}: missing mutations")
                continue
            # Per-gene guard. _run_single_gene already catches failures of the
            # enrichment analysis itself, but nothing outside that inner try was
            # guarded: read_fasta, trim_muts, process_mutations and write_output
            # all ran bare, so a single bad input killed the run and every gene
            # after it. A gene that fails is now named and counted, and the rest
            # of the directory still runs.
            try:
                _run_single_gene(gene, str(msa_file), mut_file,
                                 args, args.wt_gi, args.output)
            except Exception as exc:
                failed_genes.append((gene, f"{type(exc).__name__}: {exc}"))
                print(f"Error processing {gene}: {type(exc).__name__}: {exc}",
                      file=sys.stderr)
        if failed_genes:
            print(f"\n{len(failed_genes)}/{len(msa_files)} genes failed:", file=sys.stderr)
            for gene, reason in failed_genes:
                print(f"  {gene}: {reason}", file=sys.stderr)
    else:
        # Single-gene mode is deliberately NOT guarded: there is no later gene to
        # protect, and a swallowed exception here would exit 0 with no output.
        gene = extract_gene_from_filename(args.msa)
        _run_single_gene(gene, args.msa, args.mutations,
                         args, args.wt_gi, args.output)


if __name__ == "__main__":
    main()
