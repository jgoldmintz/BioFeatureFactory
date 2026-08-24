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
NetSurfP-3.0 Pipeline for Protein Structure Prediction

Predicts surface accessibility, secondary structure, and disorder for WT and mutant sequences.
Uses nsp3 Python library with trained model for direct prediction.
Generates ensemble TSV outputs with WT vs MUT structural comparisons.

Key features:
- Surface accessibility prediction (RSA)
- Secondary structure prediction (Q8/Q3)
- Disorder region prediction (disorder_pf, disorder_pt)
- Backbone torsion angles (phi, psi)
- WT vs mutant comparison with delta scores
- Three-tier output: summary, per-residue, local context
- Integration with nucleotide mutation mapping
- Sequence chunking for long proteins
- Non-SNV variants (indel / MNV / delins / frameshift) processed BY DEFAULT

Non-SNV handling, in brief:
  * There is no flag. Whether a column survives is decided per column from the
    record, not from a per-run preference.
  * WT and mutant residues are joined through align_wt_to_mut, never by integer
    position equality -- after an indel, residue i of the two alleles are
    different sites. residues.tsv and local.tsv carry that projection in
    wt_pos/mut_pos plus an align_status of aligned/deleted/inserted.
  * A value that cannot be computed is EMPTY with a named reason in qc_flags,
    never 0.0, which would read as "measured, no change".
  * A token that yields no mutant sequence still gets a summary row naming why
    (REF_MISMATCH:orf_has_X, REF_SPANS_PAST_ORF, UNPARSEABLE_TOKEN, ...).
"""

import os
import csv
import tempfile
from pathlib import Path
import sys
from typing import Optional, Dict, List, Tuple
# Import utility functions
from biofeaturefactory.lib.utility import (
    read_fasta,
    mint_pkey,
    token_from_name,
    get_mutation_data_bioAccurate,
    load_validation_failures,
    should_skip_mutation,
    get_mutant_aa,
    extract_gene_from_filename,
    discover_fasta_files,
    translate_orf_sequence,
    build_mutant_sequences_for_gene,
    detect_alphabet,
    write_tsv,
    # Length-aware layer. The SNV path still routes through
    # get_mutation_data_bioAccurate + get_mutant_aa above and is NOT re-plumbed:
    # infer_aavariant_from_nt returns the TRIMMED span, so a synonymous SNV comes
    # back anchored on the PRECEDING residue (utility.py:2106-2111) and would move
    # mutation_pos by one. Non-SNV tokens take the new path, SNVs keep the old one.
    parse_variant,
    splice_seq,
    align_wt_to_mut,
    infer_edit_span,
    infer_aavariant_from_nt,
    is_intronic_token,
)

# Import for NSP3 prediction
import pandas as pd
import numpy as np

# Add local nsp3 git clone to sys.path so it's importable from any working directory
_nsp3_package_dir = Path(__file__).resolve().parent / "nsp3" / "nsp3"
if str(_nsp3_package_dir) not in sys.path:
    sys.path.insert(0, str(_nsp3_package_dir))

from nsp3 import main as nsp3_main
from nsp3.cli import load_config

# Patch nsp3's BasePredict to use strict=False when loading state_dict.
# The ESM-1b checkpoint may contain extra keys (e.g., emb_layer_norm_before)
# not present in the model constructed from config.
# Wrapped in try/except so future nsp3 updates that fix this upstream won't break the pipeline.
try:
    import nsp3.base.base_predict as _base_predict

    def _patched_base_init(self, model, model_data, *args, **kwargs):
        import torch as _torch
        super(_base_predict.BasePredict, self).__init__()
        self.model = model
        # Load the checkpoint onto the available device: CUDA if present, else
        # CPU. (Inference device itself is governed by nsp3's setup_device /
        # SecondaryFeatures.inference; this just avoids a needless CPU->GPU copy
        # on GPU hosts and keeps deserialization device-aware.)
        _device = _torch.device('cuda' if _torch.cuda.is_available() else 'cpu')
        print(f"Loading model onto {_device.type}... \n")
        data = _torch.load(model_data, map_location=_device)
        self.model.load_state_dict(data['state_dict'], strict=False)
        self.model.eval()

    _base_predict.BasePredict.__init__ = _patched_base_init
except Exception:
    pass


# Q8/Q3 class labels
Q8_LABELS = "GHIBESTC"
Q3_LABELS = "HEC"

# Max ASA values per residue (Tien et al. 2013) for computing ASA from RSA
MAX_ASA = {
    'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167,
    'E': 223, 'Q': 225, 'G': 104, 'H': 224, 'I': 197,
    'L': 201, 'K': 236, 'M': 224, 'F': 240, 'P': 159,
    'S': 155, 'T': 172, 'W': 285, 'Y': 263, 'V': 174
}


def extract_residue_predictions(predictions_batch, seq_idx, pos_idx, residue):
    """
    Extract per-residue predictions from model output tensors.

    Model outputs 6 tensors after postprocessing:
        predictions_batch[0] = ss8  (batch, seq, 8) - Q8 secondary structure probs
        predictions_batch[1] = ss3  (batch, seq, 3) - Q3 secondary structure probs
        predictions_batch[2] = dis  (batch, seq, 2) - disorder scores (flDPnn, flDPthr)
        predictions_batch[3] = rsa  (batch, seq, 1) - relative solvent accessibility
        predictions_batch[4] = phi  (batch, seq, 1) - phi dihedral angle (post arctan)
        predictions_batch[5] = psi  (batch, seq, 1) - psi dihedral angle (post arctan)

    Args:
        predictions_batch: List of 6 numpy arrays from model
        seq_idx: Sequence index in batch
        pos_idx: Position index in sequence
        residue: Amino acid at this position (for ASA calculation)

    Returns:
        dict with all prediction values for this residue

    Raises:
        IndexError, ValueError: If tensor indexing fails (indicates model output mismatch)
    """
    # Extract Q8 probabilities (8 classes: G,H,I,B,E,S,T,C)
    ss8_tensor = predictions_batch[0]
    q8_probs = ss8_tensor[seq_idx, pos_idx, :].flatten().tolist()
    if len(q8_probs) != 8:
        raise ValueError(f"Expected 8 Q8 probabilities, got {len(q8_probs)} at seq={seq_idx}, pos={pos_idx}")

    # Extract Q3 probabilities (3 classes: H,E,C)
    ss3_tensor = predictions_batch[1]
    q3_probs = ss3_tensor[seq_idx, pos_idx, :].flatten().tolist()
    if len(q3_probs) != 3:
        raise ValueError(f"Expected 3 Q3 probabilities, got {len(q3_probs)} at seq={seq_idx}, pos={pos_idx}")

    # Extract disorder scores (2 values: flDPnn, flDPthr)
    dis_tensor = predictions_batch[2]
    dis_values = dis_tensor[seq_idx, pos_idx, :].flatten()
    if len(dis_values) != 2:
        raise ValueError(f"Expected 2 disorder values, got {len(dis_values)} at seq={seq_idx}, pos={pos_idx}")
    disorder_pf = float(dis_values[0])
    disorder_pt = float(dis_values[1])

    # Extract RSA (single value, sigmoid output 0-1)
    rsa_tensor = predictions_batch[3]
    rsa = float(rsa_tensor[seq_idx, pos_idx, 0])

    # Extract phi angle (single value after arctan conversion)
    phi_tensor = predictions_batch[4]
    phi = float(phi_tensor[seq_idx, pos_idx, 0])

    # Extract psi angle (single value after arctan conversion)
    psi_tensor = predictions_batch[5]
    psi = float(psi_tensor[seq_idx, pos_idx, 0])

    # Determine predicted classes
    q8_class = Q8_LABELS[np.argmax(q8_probs)]
    q3_class = Q3_LABELS[np.argmax(q3_probs)]

    # Compute ASA from RSA
    max_asa = MAX_ASA.get(residue.upper(), 200)
    asa = rsa * max_asa

    return {
        'residue': residue,
        'q8_g': q8_probs[0],
        'q8_h': q8_probs[1],
        'q8_i': q8_probs[2],
        'q8_b': q8_probs[3],
        'q8_e': q8_probs[4],
        'q8_s': q8_probs[5],
        'q8_t': q8_probs[6],
        'q8_c': q8_probs[7],
        'q3_h': q3_probs[0],
        'q3_e': q3_probs[1],
        'q3_c': q3_probs[2],
        'disorder_pf': disorder_pf,
        'disorder_pt': disorder_pt,
        'rsa': rsa,
        'asa': asa,
        'phi': phi,
        'psi': psi,
        'q8': q8_class,
        'q3': q3_class,
    }


def discover_mutation_files(mutation_path_str):
    """
    Discover mutation files for each gene.

    Handles both:
    - Single mutation file: extracts gene name, returns {gene: file}
    - Directory: scans for *.csv files matching patterns

    Looks for files matching:
    - <GENE>_mutations.csv
    - combined_<GENE>.csv
    - <GENE>.csv

    Args:
        mutation_path_str: Path to mutation file or directory

    Returns:
        dict: {gene_name: file_path}
    """
    if not mutation_path_str or not Path(mutation_path_str).exists():
        return {}

    mutation_path = Path(mutation_path_str)
    mutation_files = {}

    # Handle single file
    if mutation_path.is_file():
        gene_name = extract_gene_from_filename(mutation_path.name)
        mutation_files[gene_name] = str(mutation_path)
        return mutation_files

    # Handle directory
    for csv_file in mutation_path.glob("*.csv"):
        gene_name = extract_gene_from_filename(csv_file.name)
        mutation_files[gene_name] = str(csv_file)

    return mutation_files


def run_nsp3_prediction(fasta_file, model_path, config_path, batch_size=100, verbose=False, max_seq_length=1500):
    """
    Run NetSurfP-3.0 prediction using the nsp3 Python library.

    Processes sequences in batches and chunks long sequences to avoid LSTM dimension errors.

    Args:
        fasta_file: Input FASTA file with AA sequences
        model_path: Path to trained NSP3 model checkpoint
        config_path: Path to NSP3 config YAML file
        batch_size: Number of sequences to process per batch (default: 100)
        verbose: Print batch progress
        max_seq_length: Maximum sequence length before chunking (default: 1500)

    Returns:
        dict: {sequence_id: {pos: {residue, q8_probs, q3_probs, disorder_pf, disorder_pt, rsa, asa, phi, psi, q8_class, q3_class}}}
    """
    # Read all sequences from FASTA
    sequences_dict = read_fasta(fasta_file)

    # Check for long sequences and chunk if necessary
    processed_sequences = []
    chunk_mapping = {}  # Maps chunk_id -> (original_id, start_pos)

    for seq_id, seq in sequences_dict.items():
        seq_len = len(seq)

        if seq_len <= max_seq_length:
            # Short sequence - process as is
            processed_sequences.append((seq_id, seq))
            chunk_mapping[seq_id] = (seq_id, 0, seq_len)
        else:
            # Long sequence - chunk with overlap
            if verbose:
                print(f"      Chunking {seq_id} ({seq_len} AA) into segments of {max_seq_length} AA")

            overlap = 50  # Overlap between chunks
            chunk_size = max_seq_length - overlap

            for i, start in enumerate(range(0, seq_len, chunk_size)):
                end = min(start + max_seq_length, seq_len)
                chunk_id = f"{seq_id}_chunk{i}"
                chunk_seq = seq[start:end]
                processed_sequences.append((chunk_id, chunk_seq))
                chunk_mapping[chunk_id] = (seq_id, start, end)

                if end >= seq_len:
                    break

    total_chunks = len(processed_sequences)
    all_predictions_raw = {}
    config = load_config(config_path)

    # Process chunks in batches
    num_batches = (total_chunks + batch_size - 1) // batch_size

    for batch_num, batch_start in enumerate(range(0, total_chunks, batch_size), 1):
        batch_end = min(batch_start + batch_size, total_chunks)
        batch_sequences = processed_sequences[batch_start:batch_end]

        if verbose:
            print(
                f"    Batch {batch_num}/{num_batches}: processing chunks {batch_start + 1}-{batch_end} of {total_chunks}")

        # Create temporary FASTA for this batch
        batch_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        try:
            for chunk_id, chunk_seq in batch_sequences:
                batch_fasta.write(f">{chunk_id}\n{chunk_seq}\n")
            batch_fasta.close()

            # Run prediction on batch with error handling
            try:
                result = nsp3_main.predict(config, "SecondaryFeatures", model_path, batch_fasta.name)

                # Process results
                # NSP3 returns (identifiers_list, sequences_list, predictions_list)
                # SecondaryFeatures.__call__ internally chunks by 25, returning:
                #   identifiers = [[chunk0_ids], [chunk1_ids], ...]
                #   sequences   = [[chunk0_seqs], [chunk1_seqs], ...]
                #   predictions = [[chunk0_tensors], [chunk1_tensors], ...]
                # Each chunk's tensors have their own seq_idx space (0..len(chunk)-1)
                identifiers_all = result[0]
                sequences_all = result[1]
                predictions_all = result[2]

                # Normalize to list-of-chunks if not already nested
                if identifiers_all and not isinstance(identifiers_all[0], list):
                    identifiers_all = [identifiers_all]
                    sequences_all = [sequences_all]
                    predictions_all = [predictions_all]

                # Process each NSP3 internal chunk
                for chunk_idx in range(len(identifiers_all)):
                    identifiers_chunk = identifiers_all[chunk_idx]
                    sequences_chunk = sequences_all[chunk_idx]
                    predictions_chunk = predictions_all[chunk_idx]

                    for seq_idx in range(len(identifiers_chunk)):
                        chunk_id = identifiers_chunk[seq_idx]
                        sequence = sequences_chunk[seq_idx]
                        seq_len = len(sequence)

                        per_residue_predictions = {}

                        for pos_idx in range(seq_len):
                            residue = sequence[pos_idx]
                            per_residue_predictions[pos_idx] = extract_residue_predictions(
                                predictions_chunk, seq_idx, pos_idx, residue
                            )

                        all_predictions_raw[chunk_id] = per_residue_predictions

            except RuntimeError as e:
                if "exceeds dimension size" in str(e) or "start" in str(e):
                    print(
                        f"    Warning: Batch {batch_num} contains sequences too long for model, reducing batch size...")
                    # Try processing sequences one by one
                    for chunk_id, chunk_seq in batch_sequences:
                        single_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
                        try:
                            single_fasta.write(f">{chunk_id}\n{chunk_seq}\n")
                            single_fasta.close()

                            try:
                                result = nsp3_main.predict(config, "SecondaryFeatures", model_path, single_fasta.name)
                                # Process single result (same chunk iteration as above)
                                identifiers_all = result[0]
                                sequences_all = result[1]
                                predictions_all = result[2]

                                if identifiers_all and not isinstance(identifiers_all[0], list):
                                    identifiers_all = [identifiers_all]
                                    sequences_all = [sequences_all]
                                    predictions_all = [predictions_all]

                                for ci in range(len(identifiers_all)):
                                    ids_c = identifiers_all[ci]
                                    seqs_c = sequences_all[ci]
                                    preds_c = predictions_all[ci]

                                    for seq_idx in range(len(ids_c)):
                                        chunk_id = ids_c[seq_idx]
                                        sequence = seqs_c[seq_idx]
                                        seq_len = len(sequence)

                                        per_residue_predictions = {}

                                        for pos_idx in range(seq_len):
                                            residue = sequence[pos_idx]
                                            per_residue_predictions[pos_idx] = extract_residue_predictions(
                                                preds_c, seq_idx, pos_idx, residue
                                            )

                                        all_predictions_raw[chunk_id] = per_residue_predictions

                            except Exception as e2:
                                print(f"      Skipping {chunk_id}: {e2}")

                        finally:
                            if os.path.exists(single_fasta.name):
                                os.unlink(single_fasta.name)
                else:
                    raise

        finally:
            # Clean up temporary batch file
            if os.path.exists(batch_fasta.name):
                os.unlink(batch_fasta.name)

    # Reassemble chunked predictions
    all_predictions = {}

    for chunk_id, chunk_predictions in all_predictions_raw.items():
        original_id, start_pos, end_pos = chunk_mapping[chunk_id]

        if original_id not in all_predictions:
            all_predictions[original_id] = {}

        # Map chunk positions back to original positions
        for chunk_pos, pred_data in chunk_predictions.items():
            original_pos = start_pos + chunk_pos + 1  # 1-indexed

            # For overlapping regions, prefer the prediction from the middle of a chunk
            if original_pos not in all_predictions[original_id]:
                all_predictions[original_id][original_pos] = pred_data
            else:
                # Already have a prediction for this position from another chunk
                # Keep the one that's further from chunk boundaries
                pass

    return all_predictions


def classify_burial_change(wt_rsa, mut_rsa):
    """
    Classify RSA change as buried/intermediate/exposed transition.

    RSA categories (encoded as levels 0, 1, 2):
    - buried (0): RSA < 0.25
    - intermediate (1): 0.25 <= RSA < 0.50
    - exposed (2): RSA >= 0.50

    Returns:
        int: Transition score from -2 to +2
             Positive = toward exposed, Negative = toward buried
             +2 = buried->exposed, -2 = exposed->buried, 0 = no change
    """

    def rsa_level(rsa):
        if rsa < 0.25:
            return 0  # buried
        elif rsa < 0.50:
            return 1  # intermediate
        else:
            return 2  # exposed

    wt_level = rsa_level(wt_rsa)
    mut_level = rsa_level(mut_rsa)

    return mut_level - wt_level


def classify_disorder_change(wt_disorder, mut_disorder):
    """
    Classify disorder change as ordered/intermediate/disordered transition.

    Disorder categories (encoded as levels 0, 1, 2):
    - ordered (0): disorder < 0.3
    - intermediate (1): 0.3 <= disorder < 0.7
    - disordered (2): disorder >= 0.7

    Returns:
        int: Transition score from -2 to +2
             Positive = toward disordered, Negative = toward ordered
             +2 = ordered->disordered, -2 = disordered->ordered, 0 = no change
    """

    def disorder_level(d):
        if d < 0.3:
            return 0  # ordered
        elif d < 0.7:
            return 1  # intermediate
        else:
            return 2  # disordered

    wt_level = disorder_level(wt_disorder)
    mut_level = disorder_level(mut_disorder)

    return mut_level - wt_level


def compute_qc_flags(wt_aa_expected, wt_aa_actual, mut_aa_expected, mut_aa_actual, phi, psi):
    """
    Compute QC flags for structural predictions.

    Flags:
    - PASS: Normal prediction
    - ALIGNMENT_MISMATCH: WT AA doesn't match expected from sequence
    - STRUCTURAL_ANOMALY: Extreme phi/psi values
    - MUTANT_AA_MISMATCH: Mutant AA doesn't match expected

    Returns:
        str: Comma-separated list of flags or "PASS"
    """
    flags = []

    if wt_aa_expected and wt_aa_actual and wt_aa_expected.upper() != wt_aa_actual.upper():
        flags.append("ALIGNMENT_MISMATCH")

    if mut_aa_expected and mut_aa_actual and mut_aa_expected.upper() != mut_aa_actual.upper():
        flags.append("MUTANT_AA_MISMATCH")

    # Check for extreme phi/psi angles (should be in range -180 to 180)
    if abs(phi) > 180 or abs(psi) > 180:
        flags.append("STRUCTURAL_ANOMALY")

    return ",".join(flags) if flags else "PASS"


# ---------------------------------------------------------------------------
# Non-SNV support.
#
# There is no CLI flag. Non-SNV tokens are processed by default: the token
# grammar is uniquely decodable, so no parser mode has to be selected, and
# whether a given column survives is decided PER COLUMN from the record itself
# (does this WT residue have a mutant counterpart?) rather than from a per-run
# preference that could only permit everything or refuse everything.
#
# The one thing an indel breaks here is the assumption that residue i of the WT
# protein and residue i of the mutant protein are the same site. Every
# cross-allele read below therefore goes through align_wt_to_mut, never through
# integer position equality.
# ---------------------------------------------------------------------------

# Same key list build_mutant_sequences_for_gene uses (utility.py:2249) so token
# enumeration here sees exactly the tokens the mutant builder saw.
_MUTANT_COLUMN_KEYS = ('mutant', 'mutation', 'nt_mutation', 'ntmutant')


def _iter_mutation_tokens(mapping_file):
    """Yield the mutation tokens in a mapping file, in order.

    Mirrors the format detection in build_mutant_sequences_for_gene
    (utility.py:2242-2256): a first line containing a comma AND one of the
    mutant/mutation/aamutant keywords means CSV-with-header, anything else is
    treated as a single column whose literal 'mutant' header line is dropped.

    trim_muts is deliberately NOT used for this: it drops line 0 unconditionally
    (utility.py:206-207), which silently eats the first mutation of a header-less
    single-column file, and this list is the denominator for "was every token
    accounted for".
    """
    tokens = []
    try:
        with open(mapping_file, 'r') as handle:
            lines = handle.readlines()
    except OSError:
        return tokens

    is_single_column = True
    if lines and ',' in lines[0]:
        first_line_lower = lines[0].lower()
        if any(k in first_line_lower for k in ('mutant', 'mutation', 'aamutant')):
            is_single_column = False

    if is_single_column:
        for line in lines:
            token = line.strip()
            if not token or token.lower() == 'mutant':
                continue
            tokens.append(token.replace(" ", ""))
    else:
        reader = csv.DictReader(lines)
        for row in reader:
            for key in _MUTANT_COLUMN_KEYS:
                if row.get(key):
                    tokens.append(row[key].strip().replace(" ", ""))
                    break
    return tokens


def _build_aa_nonsnv_mutants(gene_name, wt_aa_seq, tokens, failure_map, pkey_map=None):
    """Mutant proteins for non-SNV AMINO-ACID tokens. Returns {header: sequence}.

    BLOCKER WORKAROUND. build_mutant_sequences_for_gene cannot do this: its
    non-SNV branch is _non_snv_mutant_protein, which returns None when there is
    no nucleotide sequence (utility.py:2168) because its whole method is "splice
    the ORF and retranslate". In --input-type aa there is no ORF, so an aa-level
    indel token (KE100K) falls through to get_mutation_data_bioAccurate, whose
    int(token[1:-1]) raises, and the token is dropped.

    The correct fix is a shared aa-level builder in utility.py, which is
    read-only for this change, so the splice is done here from exported
    primitives only. splice_seq is called with validate=False and the REF span is
    checked explicitly instead: splice_seq's own validator normalises U to T
    (utility.py:541-542), which is right for nucleotides and wrong for a protein,
    where U is selenocysteine and T is threonine.
    """
    built = {}
    for token in tokens:
        variant = parse_variant(token, is_nt=False)
        if variant is None or variant.is_snv:
            continue  # SNVs are the shared builder's job; unparseable is reported elsewhere
        if should_skip_mutation(gene_name, token, failure_map):
            continue
        # Bound on the END of the REF span, not its start: a multi-residue REF can
        # begin inside the protein and run off the end.
        if variant.pos0 + len(variant.ref) > len(wt_aa_seq):
            continue
        # Compare the WHOLE REF span. Checking wt_aa_seq[pos0] alone passes on any
        # multi-residue REF whose first residue happens to match.
        if wt_aa_seq[variant.pos0:variant.pos0 + len(variant.ref)].upper() != variant.ref.upper():
            continue
        key = mint_pkey(gene_name, token)
        built[key] = splice_seq(
            wt_aa_seq, variant.pos0, variant.ref, variant.alt, validate=False)
        if pkey_map is not None:
            pkey_map[key] = token
    return built


def _aa_edit_record(token, wt_nt_seq, wt_aa_seq, mut_aa_seq):
    """Resolve a token to its amino-acid-level edit. Returns a dict or None.

    Keys:
        mutation_pos  1-based residue where the declared change starts
        wt_aa/mut_aa  declared alleles, MULTI-CHARACTER strings ('' for a
                      frameshift, where no bounded residue pair exists)
        site_wt_aa/site_mut_aa  what to QC against the single residue predicted
                      at the site
        consequence   snv | mnv | inframe_del | inframe_ins | inframe_delins |
                      frameshift | stop_gained | stop_lost | synonymous
        is_snv        True for a single-residue substitution
        offset0/ref_len/alt_len  inputs to align_wt_to_mut

    The SNV branch is the ORIGINAL code path, untouched. It is NOT re-plumbed
    through infer_aavariant_from_nt, because that returns protein_consequence's
    TRIMMED span: a synonymous SNV trims to ''/'' and is then re-anchored on the
    PRECEDING residue (utility.py:2106-2111), which would move mutation_pos by one
    for every synonymous row this pipeline has ever emitted.

    The projection (offset0/ref_len/alt_len) is always taken from the two actual
    protein sequences via infer_edit_span rather than from the token. That is the
    only source that knows about truncation: _non_snv_mutant_protein cuts the
    mutant at its first stop (utility.py:2187), so a stop-gaining delins loses its
    whole C-terminus, and no token-derived length can express that.
    """
    is_nt = wt_nt_seq is not None
    variant = parse_variant(token, is_nt=is_nt)

    if variant is not None and not variant.is_snv:
        if is_nt:
            record = infer_aavariant_from_nt(token, wt_nt_seq)
            if record is None:
                return None
            aa_pos, wt_aa, mut_aa, consequence = record
        else:
            # aa-level token: the residues ARE the alleles. There is no reading
            # frame at this level, so no aa-level token can be a frameshift; the
            # class comes from the observed edit below.
            aa_pos, wt_aa, mut_aa, consequence = variant.pos, variant.ref, variant.alt, None
        is_snv = False
    else:
        # ---- SNV path: byte-for-byte the original main() logic ----
        try:
            nt_info = get_mutation_data_bioAccurate(token, is_nt=is_nt)
        except (ValueError, IndexError):
            return None
        if nt_info[0] is None:
            return None
        if not is_nt:
            # aa mode: the token IS the amino-acid change, so
            # get_mutation_data_bioAccurate already yields (aa_pos, (wt_aa, mut_aa)).
            # get_mutant_aa is an nt->aa helper and cannot run without an nt sequence.
            aa_info = ((nt_info[0], nt_info[1]), None)
        else:
            aa_info = get_mutant_aa(nt_info, wt_nt_seq)
        if not aa_info:
            return None
        (aa_pos, (wt_aa, mut_aa)), _ = aa_info
        consequence = 'snv'
        is_snv = True

    offset0, ref_len, alt_len = infer_edit_span(wt_aa_seq, mut_aa_seq)
    if consequence == 'frameshift':
        # Nothing downstream of a frameshift aligns, and that cannot be inferred
        # from the sequences: prefix/suffix trimming returns the MINIMAL edit, so
        # MKEWLTCD -> MKNG reads as a 6->2 replacement and would pair E with N.
        offset0, ref_len, alt_len = infer_edit_span(wt_aa_seq, mut_aa_seq, frameshift=True)
    elif consequence is None:
        consequence = _classify_observed_edit(ref_len, alt_len, mut_aa_seq, offset0)

    return {
        'mutation_pos': int(aa_pos),
        'wt_aa': wt_aa,
        'mut_aa': mut_aa,
        # QC compares the declared allele against the ONE residue predicted at the
        # site, so a multi-residue span is compared on its first residue. SNV rows
        # pass the full value so multi-character 'Stop' keeps its current handling.
        'site_wt_aa': wt_aa if is_snv else wt_aa[:1],
        'site_mut_aa': mut_aa if is_snv else mut_aa[:1],
        'consequence': consequence,
        'is_snv': is_snv,
        'offset0': offset0,
        'ref_len': ref_len,
        'alt_len': alt_len,
    }


def _classify_observed_edit(ref_len, alt_len, mut_aa_seq, offset0):
    """Consequence class for an aa-level token, read off the observed edit.

    Uses the same vocabulary protein_consequence produces for nucleotide tokens
    so the column has one domain regardless of input type.
    """
    if ref_len == 0 and alt_len == 0:
        return 'synonymous'
    if '*' in mut_aa_seq[offset0:offset0 + alt_len]:
        return 'stop_gained'
    if alt_len == 0:
        return 'inframe_del'
    if ref_len == 0:
        return 'inframe_ins'
    if ref_len == alt_len:
        return 'snv' if ref_len == 1 else 'mnv'
    return 'inframe_delins'


def _rows_for_mutant(gene_name, mut_id, wt_pred, mut_pred, edit, local_window=5):
    """Build (summary_row, residue_rows, local_rows) for one mutant.

    Every WT<->MUT read goes through `mut_of_wt`, the projection produced by
    align_wt_to_mut. Indexing mut_pred by a WT position instead would compare a
    residue to whatever the indel slid into that slot.
    """
    # max(), not len(): the projection must span every position that HAS a
    # prediction, and a dropped chunk leaves the dict sparse. The accounting
    # below uses len(), which is the count of residues actually predicted.
    wt_span = max(wt_pred) if wt_pred else 0
    projection = align_wt_to_mut(wt_span, edit['offset0'], edit['ref_len'], edit['alt_len'])

    # WT 1-based position -> MUT 1-based position, restricted to positions both
    # alleles actually have a prediction for.
    mut_of_wt = {}
    for i, j in enumerate(projection):
        if j is not None and (i + 1) in wt_pred and (j + 1) in mut_pred:
            mut_of_wt[i + 1] = j + 1
    wt_of_mut = {v: k for k, v in mut_of_wt.items()}

    n_aligned = len(mut_of_wt)
    n_deleted = len(wt_pred) - n_aligned
    n_inserted = len(mut_pred) - n_aligned

    flags = []
    site_wt = edit['mutation_pos']
    site_mut = mut_of_wt.get(site_wt)

    # ---- mutation site ----
    if site_wt not in wt_pred:
        site = None
        flags.append('SITE_OUTSIDE_WT_PREDICTION')
    elif site_mut is None:
        # The edit removed this residue (or a frameshift re-read it), so there is
        # nothing to subtract. EMPTY, never 0.0 -- a zero here reads downstream as
        # "measured, no change" and is indistinguishable from a real null.
        site = None
        flags.append('SITE_HAS_NO_MUT_COUNTERPART')
    else:
        wt_res, mut_res = wt_pred[site_wt], mut_pred[site_mut]
        site = {
            'delta_rsa': mut_res['rsa'] - wt_res['rsa'],
            'delta_disorder_pf': mut_res['disorder_pf'] - wt_res['disorder_pf'],
            'delta_disorder_pt': mut_res['disorder_pt'] - wt_res['disorder_pt'],
            'ss3_change': 0 if wt_res['q3'] == mut_res['q3'] else 1,
            'ss8_change': 0 if wt_res['q8'] == mut_res['q8'] else 1,
            'burial_classification': classify_burial_change(wt_res['rsa'], mut_res['rsa']),
            'disorder_classification': classify_disorder_change(
                wt_res['disorder_pf'], mut_res['disorder_pf']),
        }
        base_qc = compute_qc_flags(edit['site_wt_aa'], wt_res['residue'],
                                   edit['site_mut_aa'], mut_res['residue'],
                                   mut_res['phi'], mut_res['psi'])
        if base_qc != 'PASS':
            flags.extend(base_qc.split(','))

    # ---- local context ----
    # Centre on the MIDPOINT OF THE REF SPAN, not on its first residue. Anchoring
    # on the first residue gives asymmetric flanks that degrade as the span grows.
    # A frameshift has no bounded span (wt_aa is deliberately empty), so its
    # centre is the shift point itself. len(wt_aa)//2 == 0 for an SNV, so this is
    # identical to the previous behaviour there.
    span_len = len(edit['wt_aa']) if edit['consequence'] != 'frameshift' else 0
    centre = edit['mutation_pos'] + span_len // 2
    site_span = range(edit['mutation_pos'], edit['mutation_pos'] + max(1, span_len))

    local_rows = []
    local_impact = 0.0
    n_local_aligned = 0
    claimed_mut = set()
    for offset in range(-local_window, local_window + 1):
        pos = centre + offset
        if pos not in wt_pred:
            continue  # outside the protein entirely -- not a row, same as before
        j = mut_of_wt.get(pos)
        if j is not None:
            claimed_mut.add(j)
        wt_local = wt_pred[pos]
        row = {
            'pkey': mut_id, 'gene': gene_name,
            'relative_pos': offset, 'absolute_pos': pos,
            'mut_pos': j if j is not None else '',
            'align_status': 'aligned' if j is not None else 'deleted',
            'wt_ss3': wt_local['q3'], 'wt_ss8': wt_local['q8'],
            'is_mutation_site': (pos in site_span),
        }
        if j is None:
            # WT side is known and is reported; the deltas and the mutant's
            # secondary structure have no referent, so they stay empty.
            row.update({'delta_rsa': '', 'delta_disorder_pf': '', 'delta_disorder_pt': '',
                        'delta_phi': '', 'delta_psi': '', 'mut_ss3': '', 'mut_ss8': ''})
        else:
            mut_local = mut_pred[j]
            delta_local_rsa = mut_local['rsa'] - wt_local['rsa']
            row.update({
                'delta_rsa': delta_local_rsa,
                'delta_disorder_pf': mut_local['disorder_pf'] - wt_local['disorder_pf'],
                'delta_disorder_pt': mut_local['disorder_pt'] - wt_local['disorder_pt'],
                'delta_phi': mut_local['phi'] - wt_local['phi'],
                'delta_psi': mut_local['psi'] - wt_local['psi'],
                'mut_ss3': mut_local['q3'], 'mut_ss8': mut_local['q8'],
            })
            local_impact += abs(delta_local_rsa)
            n_local_aligned += 1
        local_rows.append(row)

    # Residues the edit INSERTED have no WT coordinate, so none of the rows above
    # can carry them and an insertion is invisible in this table however large it
    # is. RNAfold emits exactly these rows for the same reason
    # (run_viennaRNA_pipeline.py:138-147): keyed on the mutant position, with the
    # WT-frame and delta columns empty because there is no counterpart to
    # subtract. The mutant span of this window is the interval between the
    # projected endpoints of its WT positions; an unclaimed mutant position inside
    # it was inserted there. Empty for an SNV, an MNV and a pure deletion, so the
    # substitution output is unchanged.
    for j in range(min(claimed_mut), max(claimed_mut) + 1) if claimed_mut else ():
        if j in claimed_mut or j not in mut_pred:
            continue
        mut_local = mut_pred[j]
        local_rows.append({
            'pkey': mut_id, 'gene': gene_name,
            'relative_pos': '', 'absolute_pos': '',
            'mut_pos': j, 'align_status': 'inserted',
            'wt_ss3': '', 'wt_ss8': '',
            'mut_ss3': mut_local['q3'], 'mut_ss8': mut_local['q8'],
            'delta_rsa': '', 'delta_disorder_pf': '', 'delta_disorder_pt': '',
            'delta_phi': '', 'delta_psi': '',
            # An inserted residue is part of the edit by definition.
            'is_mutation_site': True,
        })

    # ---- global metrics, over the aligned pairs only ----
    aligned = sorted(mut_of_wt.items())
    if aligned:
        rsa_deltas = [mut_pred[j]['rsa'] - wt_pred[i]['rsa'] for i, j in aligned]
        total_delta_rsa = sum(abs(d) for d in rsa_deltas)
        global_block = {
            'global_mean_delta_rsa': total_delta_rsa / len(aligned),
            'global_ss_changes': sum(1 for i, j in aligned
                                     if wt_pred[i]['q3'] != mut_pred[j]['q3']),
            'max_abs_delta_rsa': max(abs(d) for d in rsa_deltas),
            'max_abs_delta_disorder': max(
                abs(mut_pred[j]['disorder_pf'] - wt_pred[i]['disorder_pf'])
                for i, j in aligned),
        }
    else:
        # No residue of the WT survives with a counterpart. There is nothing to
        # average, so these are empty and the reason is named.
        global_block = {'global_mean_delta_rsa': '', 'global_ss_changes': '',
                        'max_abs_delta_rsa': '', 'max_abs_delta_disorder': ''}
        flags.append('NO_ALIGNED_RESIDUES')

    if edit['consequence'] == 'frameshift':
        flags.append('FRAMESHIFT:downstream_residues_also_change')
    elif not edit['is_snv'] and edit['consequence'] != 'snv':
        # A nucleotide MNV confined to one codon has protein consequence 'snv';
        # flagging that row 'SNV' says nothing, and the aa_consequence column
        # already carries the class for every row. Only classes that change the
        # residue COUNT, or that the reader must not treat as a substitution, are
        # worth a flag.
        flags.append(edit['consequence'].upper())
    if n_deleted or n_inserted:
        # Count over the UNION of both alleles' residues. aligned_N/N over WT
        # positions alone reports a 20-residue insertion as fully aligned, because
        # every WT residue does keep a counterpart -- the residues with no
        # counterpart are all on the mutant side.
        flags.append(f"length_changed:{len(mut_pred) - len(wt_pred):+d}res;"
                     f"aligned_{n_aligned}/{len(wt_pred) + n_inserted};"
                     f"deleted_{n_deleted};inserted_{n_inserted}")

    summary_row = {
        'pkey': mut_id,
        'gene': gene_name,
        'mutation_pos': edit['mutation_pos'],
        'wt_aa': edit['wt_aa'],
        'mut_aa': edit['mut_aa'],
        'aa_consequence': edit['consequence'],
        'delta_rsa': site['delta_rsa'] if site else '',
        'delta_disorder_pf': site['delta_disorder_pf'] if site else '',
        'delta_disorder_pt': site['delta_disorder_pt'] if site else '',
        'ss3_change': site['ss3_change'] if site else '',
        'ss8_change': site['ss8_change'] if site else '',
        'burial_classification': site['burial_classification'] if site else '',
        'disorder_classification': site['disorder_classification'] if site else '',
        'local_structural_impact': local_impact if n_local_aligned else '',
        'n_aa_wt': len(edit['wt_aa']),
        'n_aa_mut': len(edit['mut_aa']),
        'wt_len': len(wt_pred),
        'mut_len': len(mut_pred),
        'qc_flags': ','.join(flags) if flags else 'PASS',
    }
    summary_row.update(global_block)

    residues_rows = []
    for allele, pred, counterpart in (('wt', wt_pred, mut_of_wt),
                                      ('mut', mut_pred, wt_of_mut)):
        for pos, res_data in pred.items():
            other = counterpart.get(pos)
            residues_rows.append({
                'pkey': mut_id,
                'gene': gene_name,
                'allele': allele,
                'pos': pos,
                # The join key across alleles. (pkey, pos) alone mis-joins after an
                # indel: WT residue 100 and MUT residue 100 are different sites.
                'wt_pos': pos if allele == 'wt' else (other if other is not None else ''),
                'mut_pos': (other if other is not None else '') if allele == 'wt' else pos,
                'align_status': ('aligned' if other is not None
                                 else ('deleted' if allele == 'wt' else 'inserted')),
                'residue': res_data['residue'],
                'rsa': res_data['rsa'],
                'asa': res_data['asa'],
                'disorder_pf': res_data['disorder_pf'],
                'disorder_pt': res_data['disorder_pt'],
                'phi': res_data['phi'],
                'psi': res_data['psi'],
                'q8_class': res_data['q8'],
                'q3_class': res_data['q3'],
                'q8_g': res_data['q8_g'],
                'q8_h': res_data['q8_h'],
                'q8_i': res_data['q8_i'],
                'q8_b': res_data['q8_b'],
                'q8_e': res_data['q8_e'],
                'q8_s': res_data['q8_s'],
                'q8_t': res_data['q8_t'],
                'q8_c': res_data['q8_c'],
                'q3_h': res_data['q3_h'],
                'q3_e': res_data['q3_e'],
                'q3_c': res_data['q3_c'],
            })

    return summary_row, residues_rows, local_rows


def _rejected_summary_row(gene_name, token, reason, edit=None, is_nt=True):
    """A NAMED row for a token that produced no prediction pair.

    Returning None (or `continue`) instead drops the mutation from the output
    entirely, which downstream is indistinguishable from "this gene never had
    that mutation". Every column that IS known is filled; the rest are empty and
    qc_flags carries the reason.
    """
    row = {
        'pkey': mint_pkey(gene_name, token),
        'gene': gene_name,
        'qc_flags': reason,
    }
    if edit is not None:
        row.update({
            'mutation_pos': edit['mutation_pos'],
            'wt_aa': edit['wt_aa'],
            'mut_aa': edit['mut_aa'],
            'aa_consequence': edit['consequence'],
            'n_aa_wt': len(edit['wt_aa']),
            'n_aa_mut': len(edit['mut_aa']),
        })
        return row

    # No resolved effect, but the position the token DECLARES is still a fact of
    # the record and is reported -- codon_usage fills the same coordinate columns
    # on its REF_MISMATCH / REF_SPANS_PAST_ORF rows
    # (codon_usage_pipeline.py:292-300). A nucleotide token names a codon.
    variant = parse_variant(token, is_nt=is_nt)
    if variant is not None:
        row['mutation_pos'] = (variant.pos0 // 3) + 1 if is_nt else variant.pos
    return row


def _rejection_reason(token, wt_nt_seq, wt_aa_seq):
    """Name why a token never produced a mutant sequence.

    Re-runs the same guards the builders apply, in the same order, so the reason
    is the true first cause rather than a generic 'build failed'.
    """
    # A gd./ch. token is WELL FORMED and out of scope for a protein tool, not
    # malformed. UNPARSEABLE_TOKEN gives it the same reason as garbage;
    # build_mutant_sequences_for_gene names it correctly on stderr, and this
    # string is the only record that reaches the TSV.
    if is_intronic_token(token):
        return 'NON_ORF_TOKEN:no_reading_frame_at_protein_level'
    is_nt = wt_nt_seq is not None
    reference = wt_nt_seq if is_nt else wt_aa_seq
    variant = parse_variant(token, is_nt=is_nt)
    if variant is None:
        return 'UNPARSEABLE_TOKEN'
    if variant.pos0 + len(variant.ref) > len(reference):
        return 'REF_SPANS_PAST_ORF' if is_nt else 'REF_SPANS_PAST_PROTEIN'
    observed = reference[variant.pos0:variant.pos0 + len(variant.ref)].upper()
    if observed != variant.ref.upper():
        return f'REF_MISMATCH:{"orf" if is_nt else "protein"}_has_{observed}'
    return 'NO_MUTANT_SEQUENCE_BUILT'


def write_summary_tsv(summary_rows, output_file):
    """
    Write summary TSV with per-mutation delta metrics.

    wt_aa/mut_aa are MULTI-CHARACTER strings for a non-SNV (e.g. 'KE'/'K' for a
    one-residue in-frame deletion) and are empty for a frameshift, where no
    bounded residue pair exists and n_aa_wt/n_aa_mut are 0. Any delta column may
    be EMPTY -- never 0.0 -- when the site has no mutant counterpart; qc_flags
    then names why.

    write_tsv passes this list to csv.DictWriter, so a row key that is absent
    here is dropped silently -- add the column here whenever a row dict gains one.
    """
    fieldnames = [
        'pkey', 'gene', 'mutation_pos', 'wt_aa', 'mut_aa', 'aa_consequence',
        'delta_rsa', 'delta_disorder_pf', 'delta_disorder_pt',
        'ss3_change', 'ss8_change',
        'burial_classification', 'disorder_classification',
        'local_structural_impact', 'global_mean_delta_rsa', 'global_ss_changes',
        'max_abs_delta_rsa', 'max_abs_delta_disorder',
        'n_aa_wt', 'n_aa_mut', 'wt_len', 'mut_len', 'qc_flags'
    ]

    write_tsv(summary_rows, output_file, fieldnames)
    print(f"Wrote {len(summary_rows)} summary entries to {output_file}")


def write_residues_tsv(residues_rows, output_file):
    """
    Write per-residue TSV with all residue predictions for WT and mutants.

    `pos` is the position in that row's OWN allele. Joining WT to MUT on
    (pkey, pos) is only correct while both alleles have the same length: after an
    indel, WT residue 100 and MUT residue 100 are different sites. wt_pos/mut_pos
    carry the projection, so the correct cross-allele join is on
    (pkey, wt_pos) or (pkey, mut_pos), and align_status says which side a row has
    no counterpart on:
        aligned  -- this residue has a counterpart in the other allele
        deleted  -- WT residue the edit removed (mut_pos empty)
        inserted -- MUT residue the edit added (wt_pos empty)
    """
    fieldnames = [
        'pkey', 'gene', 'allele', 'pos', 'wt_pos', 'mut_pos', 'align_status',
        'residue',
        'rsa', 'asa', 'disorder_pf', 'disorder_pt', 'phi', 'psi',
        'q8_class', 'q3_class',
        'q8_g', 'q8_h', 'q8_i', 'q8_b', 'q8_e', 'q8_s', 'q8_t', 'q8_c',
        'q3_h', 'q3_e', 'q3_c'
    ]

    write_tsv(residues_rows, output_file, fieldnames)
    print(f"Wrote {len(residues_rows)} residue entries to {output_file}")


def write_local_tsv(local_rows, output_file):
    """
    Write local changes TSV with +/-5 residue window per mutation.

    Rows are emitted in the WT frame: absolute_pos is a WT residue number and
    relative_pos is its offset from the centre of the REF span. mut_pos is the
    projected mutant residue, empty when the edit deleted this one -- such a row
    keeps its wt_ss3/wt_ss8 and leaves every delta EMPTY rather than 0.0, which
    would read as "measured, unchanged".

    align_status says which frame a row belongs to:
        aligned  -- both alleles have this residue; every delta is populated
        deleted  -- WT residue the edit removed; mut_pos and the deltas are empty
        inserted -- MUT residue the edit added, so it has no WT coordinate:
                    relative_pos/absolute_pos/wt_ss3/wt_ss8 and the deltas are
                    empty and mut_pos carries the position. Without these rows an
                    insertion is invisible in this table.

    local_structural_impact in the summary is the sum of |delta_rsa| over the
    ALIGNED rows only, so it is a sum over a variable number of residues whenever
    the edit removed part of the window; count the aligned rows here before
    comparing it across variant classes.
    """
    fieldnames = [
        'pkey', 'gene', 'relative_pos', 'absolute_pos', 'mut_pos', 'align_status',
        'delta_rsa', 'delta_disorder_pf', 'delta_disorder_pt',
        'delta_phi', 'delta_psi',
        'wt_ss3', 'mut_ss3', 'wt_ss8', 'mut_ss8',
        'is_mutation_site'
    ]

    write_tsv(local_rows, output_file, fieldnames)
    print(f"Wrote {len(local_rows)} local change entries to {output_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="NetSurfP-3.0 pipeline for protein structure prediction with WT/mutant comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('input', nargs='?',
                        help='Input: WT FASTA file or directory of FASTA files (nucleotide or amino acid)')
    parser.add_argument('output', nargs='?',
                        help='Output base directory')

    # Input type
    parser.add_argument('--input-type', choices=['nt', 'aa'], default=None,
                        help='Input sequence type: "nt" for nucleotide (will translate), "aa" for amino acid. '
                             'Optional — when omitted, auto-detected from the WT sequence via detect_alphabet.')

    # Processing options
    parser.add_argument('-m', '--mutation-dir', required=True,
                        help='Mutation file or directory. For --input-type=nt: NT mutations '
                             '(A1002T, and non-SNV forms such as ACAA1002A, T28TGGT). '
                             'For --input-type=aa: AA mutations (M334V, KE100K). '
                             'Non-SNV tokens are processed by default (REQUIRED)')
    parser.add_argument('-M', '--model', required=True,
                        help='Path to trained NSP3 model checkpoint (REQUIRED)')
    parser.add_argument('-c', '--config', required=True,
                        help='Path to NSP3 config YAML file (REQUIRED)')
    parser.add_argument('--log',
                        help='Validation log file to skip failed mutations')
    parser.add_argument('--batch-size', type=int, default=100,
                        help='Number of sequences to process per NSP3 batch (default: 100)')
    parser.add_argument('--max-seq-length', type=int, default=1500,
                        help='Maximum sequence length before chunking (default: 1500)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()

    # Validate arguments
    if not args.input or not args.output:
        parser.error("input and output arguments are required")

    # Implementation
    if args.verbose:
        print("NetSurfP-3.0 Pipeline - Delta-based structural analysis")
        print(f"Input: {args.input}")
        print(f"Input type: {args.input_type or 'auto-detect'}")
        print(f"Output: {args.output}")
        print(f"Model: {args.model}")
        print(f"Config: {args.config}")
        print(f"Max sequence length: {args.max_seq_length}")

    # Load validation failures if provided
    failure_map = load_validation_failures(args.log) if args.log else None

    # Discover WT FASTA files
    input_path = Path(args.input)
    if input_path.is_file():
        fasta_files = {extract_gene_from_filename(input_path.name): str(input_path)}
    elif input_path.is_dir():
        fasta_files = discover_fasta_files(str(input_path))
    else:
        print(f"Error: Input path not found: {args.input}")
        return 1

    if args.verbose:
        print(f"Found {len(fasta_files)} FASTA files to process")

    # Discover mutation files
    mutation_files = discover_mutation_files(args.mutation_dir) if args.mutation_dir else {}

    # Collect all results
    summary_rows = []
    residues_rows = []
    local_rows = []

    # Process each gene
    for gene_name, fasta_path in fasta_files.items():
        if args.verbose:
            print(f"\nProcessing gene: {gene_name}")

        # Load WT sequence
        wt_sequences = read_fasta(fasta_path)
        if not wt_sequences:
            print(f"Warning: No sequences in {fasta_path}, skipping")
            continue

        # Handle multiple sequences in FASTA - prefer 'ORF' header like netMHC
        if len(wt_sequences) > 1:
            if 'ORF' in wt_sequences:
                wt_header, wt_seq = 'ORF', wt_sequences['ORF']
            else:
                wt_header, wt_seq = next(iter(wt_sequences.items()))
        else:
            wt_header, wt_seq = next(iter(wt_sequences.items()))

        # Resolve input alphabet: an explicit --input-type overrides; otherwise
        # auto-detect from the WT sequence. NetSurfP3 is a protein tool, so only
        # nt (translated) and aa are valid; codon-encoded input is skipped.
        if args.input_type is not None:
            seq_type = args.input_type
        else:
            try:
                detected = detect_alphabet(wt_seq)
            except ValueError:
                print(f"Warning: {gene_name} WT sequence is empty, skipping")
                continue
            if detected == 'nucleotide':
                seq_type = 'nt'
            elif detected == 'protein':
                seq_type = 'aa'
            else:  # codon-encoded
                print(f"Warning: {gene_name} WT looks codon-encoded; NetSurfP3 needs nt or aa, skipping")
                continue

        # Handle based on resolved input type
        if seq_type == 'nt':
            # Nucleotide input - translate to amino acids
            wt_nt_seq = wt_seq
            wt_aa_seq = translate_orf_sequence(wt_nt_seq)
            if not wt_aa_seq:
                print(f"Warning: Could not translate {gene_name}, skipping")
                continue
        else:
            # Amino acid input - use directly
            wt_nt_seq = None
            wt_aa_seq = wt_seq
            # Basic validation - check for non-AA characters. U/O are real residues
            # and X/B/Z are IUPAC ambiguity codes; ESM-1b has a learned token for each,
            # so NSP3 scores them normally and one ambiguous residue no longer costs
            # the gene every one of its mutations.
            valid_aa = set('ACDEFGHIKLMNPQRSTVWYUOXBZ*')
            invalid_chars = sorted({c.upper() for c in wt_aa_seq if c.upper() not in valid_aa})
            if invalid_chars:
                print(f"Warning: {gene_name} contains non-amino acid characters ({''.join(invalid_chars)}), skipping")
                continue

        if args.verbose:
            print(f"  WT sequence length: {len(wt_aa_seq)} AA")

        # Build mutant sequences
        mutation_file = mutation_files.get(gene_name)
        if not mutation_file:
            print(f"Warning: No mutation file for {gene_name}, skipping")
            continue

        # non_snp=True is not a user choice and gets no flag. With it False the
        # shared builder routes every multi-base token into
        # get_mutation_data_bioAccurate, whose int(token[1:-1]) raises, and the
        # token dies inside utility before this pipeline ever sees it. An indel's
        # mutant protein is the translation of the edited ORF, which is what
        # non_snp=True builds; update_str on the aa sequence cannot express one.
        # Sequence names are {GENE}-{sha} now; pkey_map carries name -> token so
        # token_from_name recovers the spelling without parsing the name.
        pkey_map = {}
        mutant_seqs = build_mutant_sequences_for_gene(
            gene_name, wt_nt_seq, wt_aa_seq, mutation_file, args.log, failure_map,
            input_type=seq_type, non_snp=True, pkey_map=pkey_map
        )

        # The token list is the denominator for "was every mutation accounted
        # for". Tokens that never became a sequence get a NAMED row below rather
        # than vanishing.
        tokens = _iter_mutation_tokens(mutation_file)

        if wt_nt_seq is None:
            # BLOCKER WORKAROUND (see _build_aa_nonsnv_mutants): the shared builder
            # cannot splice a non-SNV amino-acid token, because its non-SNV branch
            # needs an ORF to re-translate and aa mode has none.
            for header, seq in _build_aa_nonsnv_mutants(
                    gene_name, wt_aa_seq, tokens, failure_map, pkey_map).items():
                mutant_seqs.setdefault(header, seq)

        if args.verbose:
            print(f"  Generated {len(mutant_seqs)} mutant sequences")

        # Account for every token that never became a mutant sequence, BEFORE the
        # early exit below. A gene whose tokens all fail still owes the caller one
        # named row per token; skipping the gene silently would make those
        # mutations indistinguishable from mutations that were never submitted.
        # Tokens the validation log deliberately EXCLUDED are not failures, and
        # should_skip_mutation is exactly that set: the mutations the log names for
        # this gene (utility.py:2262).
        #
        # The shared builder's own `allowed_mutations` set (utility.py:2221-2230)
        # is deliberately NOT mirrored here. It is trim_muts' output, and trim_muts
        # drops line 0 unconditionally (utility.py:206-207) and strips '*'
        # (utility.py:208). Filtering on it therefore suppressed the row for the
        # FIRST mutation of a header-less single-column file, and for any token
        # carrying '*', on every --log run -- the token was neither built nor
        # reported, which is the silent drop this accounting exists to remove.
        # Verified: 3 tokens in, 2 rows out, T11G gone without a trace.
        built_tokens = {token_from_name(mut_id, gene_name, pkey_map) for mut_id in mutant_seqs}
        for token in tokens:
            if token in built_tokens:
                continue
            if should_skip_mutation(gene_name, token, failure_map):
                continue
            summary_rows.append(_rejected_summary_row(
                gene_name, token, _rejection_reason(token, wt_nt_seq, wt_aa_seq),
                is_nt=wt_nt_seq is not None))

        if not mutant_seqs:
            print(f"Warning: No valid mutants for {gene_name}, skipping")
            continue

        # Create combined FASTA with WT + all mutants
        combined_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        try:
            # Write WT sequence
            combined_fasta.write(f">{gene_name}-WT\n{wt_aa_seq}\n")
            # Write mutant sequences
            for mut_id, mut_seq in mutant_seqs.items():
                combined_fasta.write(f">{mut_id}\n{mut_seq}\n")
            combined_fasta.close()

            # Run NSP3 prediction
            if args.verbose:
                print(f"  Running NSP3 prediction on {len(mutant_seqs) + 1} sequences...")

            all_predictions = run_nsp3_prediction(
                combined_fasta.name,
                args.model,
                args.config,
                batch_size=args.batch_size,
                verbose=args.verbose,
                max_seq_length=args.max_seq_length
            )

            if args.verbose:
                print(f"  Prediction complete, processing results...")

            # Extract WT predictions
            wt_key = f"{gene_name}-WT"
            if wt_key not in all_predictions:
                print(f"Warning: WT predictions not found for {gene_name}, skipping")
                # Contract D is per token, not per gene. Without the WT allele
                # nothing can be compared, but every mutant that WAS built is a
                # mutation this pipeline was asked about; dropping them all behind
                # one gene-level warning makes them indistinguishable downstream
                # from mutations that were never submitted. The per-mutant sibling
                # of this branch (NO_PREDICTION, below) already emits a named row.
                # The aa effect needs no prediction to resolve, so it is reported
                # here too rather than left blank.
                for mut_id, mut_seq in mutant_seqs.items():
                    token = token_from_name(mut_id, gene_name, pkey_map)
                    summary_rows.append(_rejected_summary_row(
                        gene_name, token, 'NO_WT_PREDICTION',
                        _aa_edit_record(token, wt_nt_seq, wt_aa_seq, mut_seq),
                        is_nt=wt_nt_seq is not None))
                continue

            wt_pred = all_predictions[wt_key]

            # Process each mutant
            is_nt_input = wt_nt_seq is not None
            for mut_id, mut_seq in mutant_seqs.items():
                mutation_str = token_from_name(mut_id, gene_name, pkey_map)

                # Resolved before the prediction check so a mutant that failed
                # prediction still reports the aa effect it was going to be scored
                # for, instead of an anonymous failure row.
                edit = _aa_edit_record(mutation_str, wt_nt_seq, wt_aa_seq, mut_seq)

                if mut_id not in all_predictions:
                    print(f"Warning: Predictions not found for {mut_id}, skipping")
                    summary_rows.append(_rejected_summary_row(
                        gene_name, mutation_str, 'NO_PREDICTION', edit, is_nt_input))
                    continue

                if edit is None:
                    # Named, not dropped. A bare `continue` here is
                    # indistinguishable downstream from "this gene never had that
                    # mutation".
                    summary_rows.append(_rejected_summary_row(
                        gene_name, mutation_str, 'AA_EFFECT_UNRESOLVED',
                        is_nt=is_nt_input))
                    continue

                summary_row, residue_block, local_block = _rows_for_mutant(
                    gene_name, mut_id, wt_pred, all_predictions[mut_id], edit)
                summary_rows.append(summary_row)
                residues_rows.extend(residue_block)
                local_rows.extend(local_block)

        finally:
            # Clean up temp file
            if os.path.exists(combined_fasta.name):
                os.unlink(combined_fasta.name)

    # Write output TSVs
    input_path = Path(args.input)
    if input_path.is_file():
        gene = extract_gene_from_filename(args.input) or input_path.stem
    else:
        gene = input_path.stem
    out_dir = Path(args.output) / gene / "NetSurfP3"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path  = out_dir / f"{gene}.netsurfp3.summary.tsv"
    residues_path = out_dir / f"{gene}.netsurfp3.residues.tsv"
    local_path    = out_dir / f"{gene}.netsurfp3.local.tsv"

    write_summary_tsv(summary_rows, str(summary_path))
    write_residues_tsv(residues_rows, str(residues_path))
    write_local_tsv(local_rows, str(local_path))

    if args.verbose:
        print(f"\nPipeline complete!")
        print(f"  Summary: {summary_path}")
        print(f"  Residues: {residues_path}")
        print(f"  Local: {local_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())