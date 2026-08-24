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
BioFeatureFactory: Codon Metrics

Codon-level translational-efficiency metrics: codon counts, CAI, tAI and
bicodon context. Consumers: codon_usage and rare_codon.

Split out of utility.py, which had grown to 92 symbols. utility.py re-exports
every name here lazily, so existing `from ...utility import X` callers are
unaffected.
"""

import csv
import os
import re
import sys
import math
import shutil
import tempfile
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Dict, Optional, Tuple

from Bio.Seq import Seq

from biofeaturefactory.lib.primitives import (
    HUMAN_REFERENCE_W,
    HUMAN_TAI_WEIGHTS,
    codon_table,
    codon_to_aa,
    get_mutation_data_bioAccurate,
)


def get_codon_counts(seq):
    """
    Compute codon and codon-pair statistics from a nucleotide sequence.

    Returns:
        tuple: (codondata, codonpairdata) dictionaries containing:
            - codondata['counts']: Raw codon counts
            - codondata['RSCU']: Relative Synonymous Codon Usage
            - codondata['W']: Relative adaptiveness (codon frequency / max synonymous frequency)
            - codonpairdata['counts']: Raw bicodon counts
            - codonpairdata['RSCPU']: Relative Synonymous Codon Pair Usage
            - codonpairdata['CPS']: Codon Pair Score (ln of observed/expected)
            - codonpairdata['noln CPS']: CPS without natural log
            - codonpairdata['W_CP']: Relative adaptiveness for codon pairs
    """
    import numpy as np

    codondata = {
        "counts": {codon: 0 for codon in codon_to_aa.keys()},
        "RSCU": {codon: 0 for codon in codon_to_aa.keys()}
    }
    codonpairdata = {
        "counts": {codon1 + codon2: 0 for codon1 in codon_to_aa.keys() for codon2 in codon_to_aa.keys()},
        "RSCPU": {codon1 + codon2: 0 for codon1 in codon_to_aa.keys() for codon2 in codon_to_aa.keys()},
        "noln CPS": {},
        "CPS": {}
    }

    # Ensure sequence length is a multiple of 3
    seq_len_multiple_of_3 = (len(seq) // 3) * 3

    for i in range(0, seq_len_multiple_of_3, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and codon in codondata["counts"]:
            codondata["counts"][codon] += 1

        if i + 6 <= seq_len_multiple_of_3:
            bicodon = seq[i:i + 6]
            if bicodon in codonpairdata["counts"]:
                codonpairdata["counts"][bicodon] += 1

    # Calculate RSCU for each codon
    for codon1 in codondata["counts"].keys():
        aa = codon_to_aa.get(codon1)
        if not aa or aa == 'Stop' or aa == '-':
            codondata["RSCU"][codon1] = np.nan
            continue
        syn_codons = codon_table.get(aa, [])
        numsyn = sum(codondata["counts"].get(c, 0) for c in syn_codons)
        try:
            codondata["RSCU"][codon1] = codondata["counts"][codon1] / (numsyn / len(syn_codons))
        except ZeroDivisionError:
            codondata["RSCU"][codon1] = np.nan

    # Calculate W (relative adaptiveness) for codons
    codondata["W"] = {}
    for codon1 in codondata["RSCU"].keys():
        aa = codon_to_aa.get(codon1)
        if not aa or aa == 'Stop' or aa == '-' or np.isnan(codondata["RSCU"].get(codon1, np.nan)):
            codondata["W"][codon1] = np.nan
            continue
        syn_codons = codon_table.get(aa, [])
        max_rscu = max([codondata["RSCU"].get(c, 0) for c in syn_codons] or [1])
        codondata["W"][codon1] = codondata["RSCU"][codon1] / max_rscu if max_rscu > 0 else np.nan

    # Calculate RSCPU for codon pairs
    for cp1 in codonpairdata["counts"].keys():
        codon_a, codon_b = cp1[:3], cp1[3:]
        aa1 = codon_to_aa.get(codon_a)
        aa2 = codon_to_aa.get(codon_b)

        if not aa1 or not aa2 or aa1 in ('Stop', '-') or aa2 in ('Stop', '-'):
            codonpairdata["RSCPU"][cp1] = np.nan
            continue

        syn_cp_codons1 = codon_table.get(aa1, [])
        syn_cp_codons2 = codon_table.get(aa2, [])

        numsyn = sum(codonpairdata["counts"].get(cpa + cpb, 0)
                     for cpa in syn_cp_codons1 for cpb in syn_cp_codons2)
        syn_cp_count = len(syn_cp_codons1) * len(syn_cp_codons2)

        try:
            codonpairdata["RSCPU"][cp1] = codonpairdata["counts"][cp1] / (numsyn / syn_cp_count)
        except ZeroDivisionError:
            codonpairdata["RSCPU"][cp1] = np.nan

        # Calculate CPS
        try:
            expected_count = codondata["counts"].get(codon_a, 0) * codondata["counts"].get(codon_b, 0)
            if expected_count == 0:
                raise ZeroDivisionError
            noln_cps = codonpairdata["counts"][cp1] / expected_count
            codonpairdata["noln CPS"][cp1] = noln_cps
            codonpairdata["CPS"][cp1] = math.log(noln_cps)
        except (ZeroDivisionError, ValueError):
            codonpairdata["noln CPS"][cp1] = np.nan
            codonpairdata["CPS"][cp1] = np.nan

    # Calculate W_CP for codon pairs
    codonpairdata["W_CP"] = {}
    for cp1 in codonpairdata["RSCPU"].keys():
        if np.isnan(codonpairdata["RSCPU"].get(cp1, np.nan)):
            codonpairdata["W_CP"][cp1] = np.nan
            continue
        codon_a, codon_b = cp1[:3], cp1[3:]
        aa1 = codon_to_aa.get(codon_a)
        aa2 = codon_to_aa.get(codon_b)
        if not aa1 or not aa2:
            codonpairdata["W_CP"][cp1] = np.nan
            continue
        syn_cp_codons1 = codon_table.get(aa1, [])
        syn_cp_codons2 = codon_table.get(aa2, [])
        max_rscpu = max([codonpairdata["RSCPU"].get(cpa + cpb, 0)
                         for cpa in syn_cp_codons1 for cpb in syn_cp_codons2] or [1])
        codonpairdata["W_CP"][cp1] = codonpairdata["RSCPU"][cp1] / max_rscpu if max_rscpu > 0 else np.nan

    return codondata, codonpairdata


def compute_cai(seq, w_values=None):
    """
    Compute Codon Adaptation Index (CAI) for a sequence.

    CAI = geometric mean of W values across all codons
    CAI = exp((1/L) * sum(ln(W_i)))

    Args:
        seq: Nucleotide sequence (in-frame ORF)
        w_values: Dict of codon -> W values. If None, uses HUMAN_REFERENCE_W

    Returns:
        float: CAI value (0-1), or None if cannot compute
    """
    import numpy as np

    if w_values is None:
        w_values = HUMAN_REFERENCE_W

    seq = seq.upper().replace('U', 'T')
    seq_len = (len(seq) // 3) * 3

    log_w_sum = 0.0
    codon_count = 0

    for i in range(0, seq_len, 3):
        codon = seq[i:i+3]
        if codon in ('TAA', 'TAG', 'TGA', '---'):  # Skip stops and gaps
            continue
        w = w_values.get(codon)
        if w is None or w <= 0:
            continue
        log_w_sum += np.log(w)
        codon_count += 1

    if codon_count == 0:
        return None

    return np.exp(log_w_sum / codon_count)


def compute_tai(seq, tai_weights=None):
    """
    Compute tRNA Adaptation Index (tAI) for a sequence.

    tAI = geometric mean of tRNA adaptation weights across all codons.

    Args:
        seq: Nucleotide sequence (in-frame ORF)
        tai_weights: Dict of codon -> tAI weights. If None, uses HUMAN_TAI_WEIGHTS

    Returns:
        float: tAI value (0-1), or None if cannot compute
    """
    import numpy as np

    if tai_weights is None:
        tai_weights = HUMAN_TAI_WEIGHTS

    seq = seq.upper().replace('U', 'T')
    seq_len = (len(seq) // 3) * 3

    log_tai_sum = 0.0
    codon_count = 0

    for i in range(0, seq_len, 3):
        codon = seq[i:i+3]
        if codon in ('TAA', 'TAG', 'TGA', '---'):  # Skip stops and gaps
            continue
        w = tai_weights.get(codon)
        if w is None or w <= 0:
            continue
        log_tai_sum += np.log(w)
        codon_count += 1

    if codon_count == 0:
        return None

    return np.exp(log_tai_sum / codon_count)


def get_codon_tai(codon, tai_weights=None):
    """Get tAI weight for a single codon."""
    if tai_weights is None:
        tai_weights = HUMAN_TAI_WEIGHTS
    return tai_weights.get(codon.upper().replace('U', 'T'), None)


def get_codon_cai_w(codon, w_values=None):
    """Get CAI W value (reference adaptiveness) for a single codon."""
    if w_values is None:
        w_values = HUMAN_REFERENCE_W
    return w_values.get(codon.upper().replace('U', 'T'), None)


def extract_codon_with_bicodons(ntposnt, seq):
    """
    Extract codon and bicodons for a given SNP, respecting biological constraints.

    Biology rules:
    - First codon: only forward bicodon possible (codon1 + codon2)
    - Last codon: only reverse bicodon possible (codon_n-1 + codon_n)
    - Middle codons: both forward and reverse bicodons possible

    Args:
        ntposnt (str): SNP notation (e.g., "A123G") - 1-based position
        seq (str): DNA sequence

    Returns:
        tuple: (original_codon, forward_bicodon, reverse_bicodon, pos_in_codon, pos, codon_number)
               where bicodons may be empty strings if not biologically possible
    """
    pos, mut = get_mutation_data_bioAccurate(ntposnt, is_nt=True)
    if pos is None:
        return None, "", "", 0, 0, 0

    # Convert to 0-based indexing
    pos_0_indexed = pos - 1

    # F46: apply the ALT base before slicing. Previously `mut` was bound and never
    # used, so every codon/bicodon (and every RSCU/W/CAI/tAI derived from them) was
    # sliced from the unmodified WT ORF while the caller labeled it 'mutated_codon'
    # — the mutation was invisible on the only CLI-reachable path. Only the mutated
    # codon changes; neighbouring codons in the bicodons stay WT, which is correct.
    if pos_0_indexed < 0 or pos_0_indexed >= len(seq):
        return None, "", "", 0, 0, 0
    ref_nt, alt_nt = mut
    # F46 REF GUARD: splicing ALT onto a base that is not REF invents a codon present
    # in NEITHER the WT nor the true mutant sequence. Measured on the production corpus
    # (59 genes / 35,089 tokens vs Bio_DBs/fastas ORF records): 1,537 tokens (4.38%)
    # have a REF that disagrees with the ORF base — concentrated in MECP2 295/423,
    # NOD2 271/387, FGFR2 215/290, MPZ 195/270 — i.e. an ORF/isoform mismatch, not noise.
    # Reject rather than fabricate; the caller turns None into a skipped row.
    if ref_nt and seq[pos_0_indexed].upper() != ref_nt.upper():
        return None, "", "", 0, 0, 0
    seq = seq[:pos_0_indexed] + alt_nt.upper() + seq[pos_0_indexed + 1:]

    pos_in_codon = pos_0_indexed % 3

    # Find the codon start position and codon number
    codon_start_pos = (pos_0_indexed // 3) * 3
    codon_number = (pos_0_indexed // 3) + 1  # 1-based codon numbering

    # Calculate total number of complete codons in sequence
    total_codons = len(seq) // 3

    # Extract the original codon containing the mutation
    original_codon = seq[codon_start_pos:codon_start_pos + 3]

    # Initialize bicodons
    forward_bicodon = ""
    reverse_bicodon = ""

    # Determine which bicodons are biologically possible
    is_first_codon = (codon_number == 1)
    is_last_codon = (codon_number == total_codons)

    # Forward bicodon (current codon + following codon)
    # Possible for first and middle codons, but not last codon
    if not is_last_codon and codon_start_pos + 6 <= len(seq):
        following_codon = seq[codon_start_pos + 3:codon_start_pos + 6]
        if len(following_codon) == 3:
            forward_bicodon = original_codon + following_codon

    # Reverse bicodon (preceding codon + current codon)
    # Possible for middle and last codons, but not first codon
    if not is_first_codon and codon_start_pos >= 3:
        preceding_codon = seq[codon_start_pos - 3:codon_start_pos]
        if len(preceding_codon) == 3:
            reverse_bicodon = preceding_codon + original_codon

    return original_codon, forward_bicodon, reverse_bicodon, pos_in_codon, pos, codon_number

