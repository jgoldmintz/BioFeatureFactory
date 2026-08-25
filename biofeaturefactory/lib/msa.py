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
BioFeatureFactory: Msa

Multiple-sequence-alignment construction: jackhmmer search, Stockholm/A2M
conversion, gap filtering and sequence-weight statistics. Consumers: EVmutation
and msa_generation_pipeline.

Split out of utility.py, which had grown to 92 symbols. utility.py re-exports
every name here with a plain eager import, so existing
`from ...utility import X` callers are unaffected.
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
    detect_alphabet,
)


def prepare_protein_query(query_fasta):
    """Return a protein FASTA path suitable for jackhmmer.

    If the query sequence is nucleotide, translates to protein (to stop codon)
    and writes the result to a temporary file.

    Returns:
        (protein_fasta_path, tmp_path_or_None)
        Caller must delete tmp_path when finished if it is not None.
    """
    import tempfile
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    record = next(SeqIO.parse(query_fasta, 'fasta'))
    alphabet = detect_alphabet(str(record.seq))

    if alphabet == 'protein':
        return query_fasta, None

    print(f"Nucleotide query detected for '{record.id}' -- translating to protein.")
    aa_seq = record.seq.translate(to_stop=True)
    if not aa_seq:
        raise ValueError(f"Translation of '{record.id}' produced an empty protein sequence.")

    aa_record = SeqRecord(aa_seq, id=record.id, description='translated')
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    SeqIO.write([aa_record], tmp, 'fasta')
    tmp.close()
    print(f"Translated protein length: {len(aa_seq)} aa")
    return tmp.name, tmp.name


def run_jackhmmer(query_fasta, database, output_sto, jackhmmer_binary,
                  iterations=5, evalue_inclusion=1e-3, threads=4, max_seqs=10000):
    """
    Run jackhmmer iterative search against a sequence database.

    Args:
        query_fasta: Path to query protein sequence (FASTA)
        database: Path to UniRef90 or similar database
        output_sto: Path for Stockholm output
        jackhmmer_binary: Path to jackhmmer executable
        iterations: Number of search iterations
        evalue_inclusion: E-value threshold for inclusion
        threads: Number of CPU threads
        max_seqs: Maximum number of sequences to include in the alignment (default: 10000).
                  Prevents profile drift for genes in large superfamilies.

    Returns:
        Path to Stockholm output file
    """
    cmd = [
        jackhmmer_binary,
        '-N', str(iterations),
        '--incE', str(evalue_inclusion),
        '-A', output_sto,
        '--cpu', str(threads),
        '--noali',
        query_fasta,
        database
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"jackhmmer failed: {result.stderr}")

    return Path(output_sto)


def parse_stockholm(stockholm_file):
    """
    Parse Stockholm format MSA file.

    Args:
        stockholm_file: Path to Stockholm file

    Returns:
        tuple: (dict {seq_id: sequence}, rf_annotation str or None)
            rf_annotation is the concatenated #=GC RF string where 'x' marks
            match columns and '.' marks insert columns. None if absent.
    """
    from collections import defaultdict
    current_seqs = defaultdict(str)
    rf_parts = []

    with open(stockholm_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Capture RF column annotation before skipping other # lines
            if line.startswith('#=GC RF'):
                parts = line.split()
                if len(parts) >= 3:
                    rf_parts.append(parts[2])
                continue

            # Skip remaining comments and empty lines
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                current_seqs[seq_id] += seq

    rf_annotation = ''.join(rf_parts) if rf_parts else None
    return dict(current_seqs), rf_annotation


def stockholm_to_a2m(msa, focus_seq_id, rf_annotation=None):
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
        rf_annotation: #=GC RF string from parse_stockholm where 'x' = match
            column and '.' = insert column. When provided, this is used as the
            authoritative match/insert column definition. Falls back to focus
            sequence non-gap positions when None.

    Returns:
        dict: {seq_id: a2m_sequence}
    """
    if focus_seq_id not in msa:
        for seq_id in msa:
            if focus_seq_id in seq_id or seq_id in focus_seq_id:
                focus_seq_id = seq_id
                break
        else:
            raise ValueError(f"Focus sequence '{focus_seq_id}' not found in MSA")

    # Identify match columns using RF annotation when available
    if rf_annotation is not None:
        match_columns = {i for i, c in enumerate(rf_annotation) if c == 'x'}
    else:
        focus_seq = msa[focus_seq_id]
        match_columns = {i for i, c in enumerate(focus_seq) if c not in '-.'}

    a2m_msa = {}
    for seq_id, seq in msa.items():
        a2m_seq = []
        for i, c in enumerate(seq):
            if i in match_columns:
                if c in '-.':
                    a2m_seq.append('-')
                else:
                    a2m_seq.append(c.upper())
            else:
                if c in '-.':
                    a2m_seq.append('.')
                else:
                    a2m_seq.append(c.lower())
        a2m_msa[seq_id] = ''.join(a2m_seq)

    return a2m_msa


def filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=0.6, a2m_format=False, focus_id=None):
    """
    Filter MSA by removing gappy sequences and columns.

    Args:
        msa: dict {seq_id: sequence}
        max_seq_gaps: Maximum fraction of gaps allowed per sequence
        max_col_gaps: Maximum fraction of gaps allowed per column
        a2m_format: When True, gap fractions are computed over match-state
            columns only (uppercase + '-'). Insert columns ('.' and lowercase)
            are structural in A2M and excluded from gap accounting. Column
            filtering also operates on match-state columns only.
        focus_id: When provided, columns where the focus sequence has a residue
            (non-gap match state) are always retained regardless of gap fraction.

    Returns:
        dict: Filtered MSA
    """
    if not msa:
        return {}

    if a2m_format:
        # Identify match-state column indices (uppercase or '-' in any sequence)
        sample = next(iter(msa.values()))
        match_cols = [i for i, c in enumerate(sample) if c == '-' or c.isupper()]

        # Remove sequences with too many deletions in match-state columns
        filtered_seqs = {}
        for seq_id, seq in msa.items():
            match_states = [seq[i] for i in match_cols]
            gap_count = sum(1 for c in match_states if c == '-')
            gap_frac = gap_count / len(match_states) if match_states else 1.0
            if gap_frac <= max_seq_gaps:
                filtered_seqs[seq_id] = seq

        if not filtered_seqs:
            return {}

        focus_seq = filtered_seqs.get(focus_id) if focus_id else None

        # Remove match-state columns with too many gaps; always retain focus residue columns
        seq_list = list(filtered_seqs.values())
        n_seqs = len(seq_list)
        cols_to_keep = []
        for i in match_cols:
            if focus_seq is not None and focus_seq[i] != '-':
                cols_to_keep.append(i)
                continue
            col = [s[i] for s in seq_list]
            gap_frac = sum(1 for c in col if c == '-') / n_seqs
            if gap_frac <= max_col_gaps:
                cols_to_keep.append(i)

        # Reconstruct full sequences keeping only retained match columns;
        # insert columns are dropped since downstream tools use match states only
        final_msa = {}
        for seq_id, seq in filtered_seqs.items():
            final_msa[seq_id] = ''.join(seq[i] for i in cols_to_keep)

        return final_msa

    # Stockholm / non-A2M path
    filtered_seqs = {}
    for seq_id, seq in msa.items():
        gap_count = seq.count('-') + seq.count('.') + seq.count('!')
        gap_frac = gap_count / len(seq) if len(seq) > 0 else 1.0
        if gap_frac <= max_seq_gaps:
            filtered_seqs[seq_id] = seq

    if not filtered_seqs:
        return {}

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

    final_msa = {}
    for seq_id, seq in filtered_seqs.items():
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        final_msa[seq_id] = new_seq

    return final_msa


def compute_sequence_weights(msa, identity_threshold=0.8, codon_mode=False):
    """
    Compute sequence weights based on clustering at identity threshold.

    Used for N_eff calculation. Weight = 1 / number of neighbors.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold (0.8 = 80% identity)
        codon_mode: If True, compare codon triplets instead of single characters.
                    Use for codon MSAs where each unit is a 3-nt codon.

    Returns:
        dict: {seq_id: weight}
    """
    seq_ids = list(msa.keys())
    n_seqs = len(seq_ids)

    if n_seqs == 0:
        return {}

    if codon_mode:
        seqs = [_chunk_codons(msa[sid]) for sid in seq_ids]
        gap_tokens = {'---', '...'}
    else:
        seqs = [msa[sid] for sid in seq_ids]
        gap_tokens = {'-', '.'}

    neighbor_counts = [1] * n_seqs

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            matches = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a == b and a not in gap_tokens)
            aligned = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a not in gap_tokens and b not in gap_tokens)
            if aligned > 0:
                identity = matches / aligned
                if identity >= identity_threshold:
                    neighbor_counts[i] += 1
                    neighbor_counts[j] += 1

    weights = {seq_ids[i]: 1.0 / neighbor_counts[i] for i in range(n_seqs)}
    return weights


def compute_neff(msa, identity_threshold=0.8, codon_mode=False):
    """
    Compute effective number of sequences (N_eff).

    N_eff = sum of sequence weights, where weight = 1/n_neighbors.
    Higher N_eff indicates more diverse MSA with better evolutionary signal.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold
        codon_mode: If True, compare codon triplets instead of single characters.

    Returns:
        float: N_eff value
    """
    weights = compute_sequence_weights(msa, identity_threshold, codon_mode=codon_mode)
    return sum(weights.values())


def _chunk_codons(seq):
    """Split a nucleotide sequence into codon triplets.

    Any triplet containing a gap character ('-' or '.') is normalized to '---'
    so the gap_tokens filter in compute_sequence_weights treats partial-gap
    codons as gaps rather than as informative residues. Without this, a column
    drop in filter_msa_by_gaps that breaks codon alignment leaks partial-gap
    codons (e.g. 'A-G') into the identity calculation, inflating both the
    match and aligned counts and skewing N_eff.

    Trailing 1-2 nt that don't complete a codon are dropped (codon MSAs are
    expected to have lengths that are multiples of 3).
    """
    codons = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        c = seq[i:i+3]
        codons.append('---' if any(ch in '-.' for ch in c) else c)
    return codons


def write_a2m(msa, output_path, focus_seq_id=None):
    """
    Write MSA to A2M format file.

    Args:
        msa: dict {seq_id: sequence}
        output_path: Output file path
        focus_seq_id: If provided, write focus sequence first
    """
    with open(output_path, 'w') as f:
        if focus_seq_id and focus_seq_id in msa:
            f.write(f">{focus_seq_id}\n{msa[focus_seq_id]}\n")
        for seq_id, seq in msa.items():
            if seq_id != focus_seq_id:
                f.write(f">{seq_id}\n{seq}\n")

