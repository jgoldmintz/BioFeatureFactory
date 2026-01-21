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
"""

import os
import csv
import tempfile
from pathlib import Path
import sys
from typing import Optional, Dict, List, Tuple
from Bio.Seq import Seq
import re

# Import utility functions
sys.path.append(os.path.join(os.path.dirname(__file__), '../utils'))

from utility import (
    read_fasta,
    get_mutation_data_bioAccurate,
    load_validation_failures,
    should_skip_mutation,
    trim_muts,
    get_mutant_aa,
    update_str,
    extract_gene_from_filename,
)

# Import for NSP3 prediction
import pandas as pd
import numpy as np
from nsp3 import main as nsp3_main
from nsp3.cli import load_config


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


def discover_fasta_files(fasta_dir):
    """
    Discover nucleotide FASTA files for each gene.

    Args:
        fasta_dir: Directory containing FASTA files

    Returns:
        dict: {gene_name: file_path}
    """
    if not fasta_dir or not Path(fasta_dir).exists():
        return {}

    fasta_files = {}
    fasta_path = Path(fasta_dir)

    for ext in ['*.fasta', '*.fa', '*.fna', '*.faa']:
        for fasta_file in fasta_path.glob(ext):
            gene_name = extract_gene_from_filename(fasta_file.stem)
            fasta_files[gene_name] = str(fasta_file)

    return fasta_files


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


def translate_orf_sequence(nt_sequence: str) -> str:
    """Translate a nucleotide ORF into an amino acid sequence, trimming trailing stops."""
    if not nt_sequence:
        return ""
    cleaned = nt_sequence.strip().upper().replace("U", "T")
    if not cleaned:
        return ""
    aa_seq = str(Seq(cleaned).translate(to_stop=False))
    return aa_seq.rstrip('*').strip()


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
                # Each element may be nested - extract properly
                identifiers_batch = result[0]
                sequences_batch = result[1]
                predictions_batch = result[2]

                # Handle nested structure - results might be [[items]] instead of [items]
                if identifiers_batch and isinstance(identifiers_batch[0], list):
                    identifiers_batch = identifiers_batch[0]
                if sequences_batch and isinstance(sequences_batch[0], list):
                    sequences_batch = sequences_batch[0]
                # Predictions are nested as [[tensor1, tensor2, ...]] - flatten to [tensor1, tensor2, ...]
                if predictions_batch and isinstance(predictions_batch[0], list):
                    predictions_batch = predictions_batch[0]

                # Process each sequence in the batch
                for seq_idx in range(len(identifiers_batch)):
                    chunk_id = identifiers_batch[seq_idx]
                    sequence = sequences_batch[seq_idx]
                    seq_len = len(sequence)

                    per_residue_predictions = {}

                    for pos_idx in range(seq_len):
                        residue = sequence[pos_idx]
                        per_residue_predictions[pos_idx] = extract_residue_predictions(
                            predictions_batch, seq_idx, pos_idx, residue
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
                                # Process single result (same as above)
                                identifiers_batch = result[0]
                                sequences_batch = result[1]
                                predictions_batch = result[2]

                                # Handle nested structure
                                if identifiers_batch and isinstance(identifiers_batch[0], list):
                                    identifiers_batch = identifiers_batch[0]
                                if sequences_batch and isinstance(sequences_batch[0], list):
                                    sequences_batch = sequences_batch[0]
                                if predictions_batch and isinstance(predictions_batch[0], list):
                                    predictions_batch = predictions_batch[0]

                                for seq_idx in range(len(identifiers_batch)):
                                    chunk_id = identifiers_batch[seq_idx]
                                    sequence = sequences_batch[seq_idx]
                                    seq_len = len(sequence)

                                    per_residue_predictions = {}

                                    for pos_idx in range(seq_len):
                                        residue = sequence[pos_idx]
                                        per_residue_predictions[pos_idx] = extract_residue_predictions(
                                            predictions_batch, seq_idx, pos_idx, residue
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

    RSA categories:
    - buried: RSA < 0.25
    - intermediate: 0.25 <= RSA < 0.50
    - exposed: RSA >= 0.50

    Returns:
        str: Transition string like "buried→exposed" or "stable"
    """

    def rsa_category(rsa):
        if rsa < 0.25:
            return "buried"
        elif rsa < 0.50:
            return "intermediate"
        else:
            return "exposed"

    wt_cat = rsa_category(wt_rsa)
    mut_cat = rsa_category(mut_rsa)

    if wt_cat == mut_cat:
        return "stable"
    return f"{wt_cat}→{mut_cat}"


def classify_disorder_change(wt_disorder, mut_disorder):
    """
    Classify disorder change as ordered/intermediate/disordered transition.

    Disorder categories:
    - ordered: disorder < 0.3
    - intermediate: 0.3 <= disorder < 0.7
    - disordered: disorder >= 0.7

    Returns:
        str: Transition string like "ordered→disordered" or "stable"
    """

    def disorder_category(d):
        if d < 0.3:
            return "ordered"
        elif d < 0.7:
            return "intermediate"
        else:
            return "disordered"

    wt_cat = disorder_category(wt_disorder)
    mut_cat = disorder_category(mut_disorder)

    if wt_cat == mut_cat:
        return "stable"
    return f"{wt_cat}→{mut_cat}"


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


def compare_wt_mut_structures(wt_predictions, mut_predictions, mutation_pos):
    """
    Compare WT and mutant structural predictions to identify changes.

    Args:
        wt_predictions: Parsed WT predictions
        mut_predictions: Parsed mutant predictions
        mutation_pos: Position of the mutation (1-indexed)

    Returns:
        dict: Comparison results with delta scores and classifications
    """
    comparison = {
        'mutation_site': {},
        'global_changes': {},
        'local_changes': {}
    }

    wt_residues = {r['pos']: r for r in wt_predictions.get('residues', [])}
    mut_residues = {r['pos']: r for r in mut_predictions.get('residues', [])}

    # Analyze mutation site
    if mutation_pos in wt_residues and mutation_pos in mut_residues:
        wt_res = wt_residues[mutation_pos]
        mut_res = mut_residues[mutation_pos]

        comparison['mutation_site'] = {
            'pos': mutation_pos,
            'wt_aa': wt_res.get('aa', ''),
            'mut_aa': mut_res.get('aa', ''),
            'delta_rsa': mut_res.get('rsa', 0) - wt_res.get('rsa', 0),
            'delta_asa': mut_res.get('asa', 0) - wt_res.get('asa', 0),
            'delta_disorder': mut_res.get('disorder', 0) - wt_res.get('disorder', 0),
            'delta_binding': mut_res.get('binding', 0) - wt_res.get('binding', 0),
            'ss3_change': 0 if wt_res.get('ss3', '') == mut_res.get('ss3', '') else 1,
            'ss8_change': 0 if wt_res.get('ss8', '') == mut_res.get('ss8', '') else 1,
        }

    # Analyze local context (±5 residues)
    local_window = 5
    local_changes = []
    for pos in range(max(1, mutation_pos - local_window),
                     min(len(wt_residues), mutation_pos + local_window) + 1):
        if pos in wt_residues and pos in mut_residues:
            wt_res = wt_residues[pos]
            mut_res = mut_residues[pos]

            delta_rsa = abs(mut_res.get('rsa', 0) - wt_res.get('rsa', 0))
            delta_disorder = abs(mut_res.get('disorder', 0) - wt_res.get('disorder', 0))

            if delta_rsa > 0.1 or delta_disorder > 0.1:  # Significant changes
                local_changes.append({
                    'pos': pos,
                    'delta_rsa': mut_res.get('rsa', 0) - wt_res.get('rsa', 0),
                    'delta_disorder': mut_res.get('disorder', 0) - wt_res.get('disorder', 0),
                })

    comparison['local_changes'] = local_changes

    # Global metrics
    all_positions = set(wt_residues.keys()) & set(mut_residues.keys())
    if all_positions:
        total_delta_rsa = sum(
            abs(mut_residues[pos].get('rsa', 0) - wt_residues[pos].get('rsa', 0))
            for pos in all_positions
        )

        ss_changes = sum(
            1 for pos in all_positions
            if wt_residues[pos].get('ss3', '') != mut_residues[pos].get('ss3', '')
        )

        comparison['global_changes'] = {
            'total_delta_rsa': total_delta_rsa,
            'mean_delta_rsa': total_delta_rsa / len(all_positions),
            'ss3_changes': ss_changes,
            'ss3_change_fraction': ss_changes / len(all_positions),
        }

    return comparison


def write_summary_tsv(summary_rows, output_file):
    """
    Write summary TSV with per-mutation delta metrics.

    Columns: pkey, gene, mutation_pos, wt_aa, mut_aa, delta_rsa, delta_disorder_pf,
             delta_disorder_pt, ss3_change, ss8_change, burial_classification,
             disorder_classification, local_structural_impact, global_mean_delta_rsa,
             global_ss_changes, max_abs_delta_rsa, max_abs_delta_disorder, qc_flags
    """
    fieldnames = [
        'pkey', 'gene', 'mutation_pos', 'wt_aa', 'mut_aa',
        'delta_rsa', 'delta_disorder_pf', 'delta_disorder_pt',
        'ss3_change', 'ss8_change',
        'burial_classification', 'disorder_classification',
        'local_structural_impact', 'global_mean_delta_rsa', 'global_ss_changes',
        'max_abs_delta_rsa', 'max_abs_delta_disorder', 'qc_flags'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {len(summary_rows)} summary entries to {output_file}")


def write_residues_tsv(residues_rows, output_file):
    """
    Write per-residue TSV with all residue predictions for WT and mutants.

    Columns: pkey, gene, allele, pos, residue, rsa, asa, disorder_pf, disorder_pt, phi, psi,
             q8_class, q3_class, q8_g, q8_h, q8_i, q8_b, q8_e, q8_s, q8_t, q8_c,
             q3_h, q3_e, q3_c
    """
    fieldnames = [
        'pkey', 'gene', 'allele', 'pos', 'residue',
        'rsa', 'asa', 'disorder_pf', 'disorder_pt', 'phi', 'psi',
        'q8_class', 'q3_class',
        'q8_g', 'q8_h', 'q8_i', 'q8_b', 'q8_e', 'q8_s', 'q8_t', 'q8_c',
        'q3_h', 'q3_e', 'q3_c'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(residues_rows)

    print(f"Wrote {len(residues_rows)} residue entries to {output_file}")


def write_local_tsv(local_rows, output_file):
    """
    Write local changes TSV with ±5 residue window per mutation.

    Columns: pkey, gene, relative_pos, absolute_pos, delta_rsa, delta_disorder_pf,
             delta_disorder_pt, delta_phi, delta_psi, wt_ss3, mut_ss3, wt_ss8, mut_ss8,
             is_mutation_site
    """
    fieldnames = [
        'pkey', 'gene', 'relative_pos', 'absolute_pos',
        'delta_rsa', 'delta_disorder_pf', 'delta_disorder_pt',
        'delta_phi', 'delta_psi',
        'wt_ss3', 'mut_ss3', 'wt_ss8', 'mut_ss8',
        'is_mutation_site'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(local_rows)

    print(f"Wrote {len(local_rows)} local change entries to {output_file}")


def build_mutant_sequences_for_gene(
        gene_name: str,
        nt_sequence: Optional[str],
        aa_sequence: str,
        mutation_file: Optional[str],
        log_path: Optional[str],
        failure_map: Optional[dict],
        input_type: str = 'nt',
):
    """
    Build mutant amino acid sequences from mutations.

    Handles both CSV format (with headers) and single-column format:
    - Single column: mutant\n{wt}{pos}{mut}\n{wt}{pos}{mut}...
    - CSV format: mutant,aamutant (with optional headers)

    Args:
        gene_name: Gene identifier
        nt_sequence: Wild-type nucleotide sequence (None if input_type='aa')
        aa_sequence: Wild-type amino acid sequence
        mutation_file: Path to mutation file (single-column or CSV)
        log_path: Optional validation log path
        failure_map: Optional map of failed mutations to skip
        input_type: 'nt' for nucleotide mutations (e.g., A1002T), 'aa' for amino acid mutations (e.g., M334V)

    Returns:
        dict: {mutation_id: mutant_aa_sequence}
    """
    if not mutation_file or not os.path.exists(mutation_file):
        return {}

    allowed_mutations = None
    if log_path:
        try:
            allowed_mutations = {
                entry.split(',')[0].strip().upper()
                for entry in trim_muts(mutation_file, log=log_path, gene_name=gene_name)
                if entry
            }
        except Exception:
            allowed_mutations = None

    mutant_sequences = {}
    try:
        with open(mutation_file, 'r') as handle:
            lines = handle.readlines()

        # Detect format: single-column or CSV
        is_single_column = True
        if lines and ',' in lines[0]:
            # Check if first line looks like CSV header
            first_line_lower = lines[0].lower()
            if any(keyword in first_line_lower for keyword in ['mutant', 'mutation', 'aamutant']):
                is_single_column = False

        if is_single_column:
            # Single-column format: each line is a mutation
            for line in lines:
                mutant_id = line.strip()
                if not mutant_id or mutant_id.lower() == 'mutant':  # Skip empty lines and header if present
                    continue

                mutant_clean = mutant_id.replace(" ", "")
                if allowed_mutations and mutant_clean.upper() not in allowed_mutations:
                    continue
                if should_skip_mutation(gene_name, mutant_clean, failure_map):
                    continue

                pos = None
                wt_aa = mut_aa = None

                if input_type == 'aa':
                    # Parse AA mutation directly (e.g., M334V)
                    aa_info = get_mutation_data_bioAccurate(mutant_clean)
                    if aa_info[0] is not None and aa_info[1]:
                        pos = aa_info[0]
                        wt_aa, mut_aa = aa_info[1]
                else:
                    # Infer AA mutation from NT mutation
                    nt_info = get_mutation_data_bioAccurate(mutant_clean)
                    if nt_info[0] is not None:
                        aa_info = get_mutant_aa(nt_info, nt_sequence)
                        if aa_info:
                            (pos, (wt_aa, mut_aa)), _ = aa_info

                if pos is None or not wt_aa or not mut_aa:
                    continue

                idx = int(pos) - 1
                if idx < 0 or idx >= len(aa_sequence):
                    continue
                if wt_aa and aa_sequence[idx].upper() != wt_aa.upper():
                    continue

                header = f"{gene_name}-{mutant_clean}"
                mutant_sequences[header] = update_str(aa_sequence, mut_aa, idx)

        else:
            # CSV format with headers
            with open(mutation_file, 'r') as handle:
                reader = csv.DictReader(handle)
                mutant_keys = ['mutant', 'mutation', 'nt_mutation', 'ntmutant']
                aa_keys = ['aamutant', 'aa_mutation', 'amino_acid_mutation', 'protein_mutation']

                for row in reader:
                    mutant_id = ""
                    for key in mutant_keys:
                        if key in row and row[key]:
                            mutant_id = row[key].strip()
                            break
                    if not mutant_id:
                        continue

                    mutant_clean = mutant_id.replace(" ", "")
                    if allowed_mutations and mutant_clean.upper() not in allowed_mutations:
                        continue
                    if should_skip_mutation(gene_name, mutant_clean, failure_map):
                        continue

                    aa_string = ""
                    for key in aa_keys:
                        if key in row and row[key]:
                            aa_string = row[key].strip()
                            break

                    pos = None
                    wt_aa = mut_aa = None
                    if aa_string:
                        pos, nts = get_mutation_data_bioAccurate(aa_string)
                        if pos is not None and nts:
                            wt_aa, mut_aa = nts

                    if pos is None or not wt_aa or not mut_aa:
                        if input_type == 'aa':
                            # For AA input, treat mutant_id as AA mutation directly
                            aa_info = get_mutation_data_bioAccurate(mutant_clean)
                            if aa_info[0] is not None and aa_info[1]:
                                pos = aa_info[0]
                                wt_aa, mut_aa = aa_info[1]
                        else:
                            # Try to infer from nucleotide mutation
                            nt_info = get_mutation_data_bioAccurate(mutant_clean)
                            if nt_info[0] is not None and nt_sequence:
                                aa_info = get_mutant_aa(nt_info, nt_sequence)
                                if aa_info:
                                    (pos, (wt_aa, mut_aa)), _ = aa_info

                    if pos is None or not wt_aa or not mut_aa:
                        continue

                    idx = int(pos) - 1
                    if idx < 0 or idx >= len(aa_sequence):
                        continue
                    if wt_aa and aa_sequence[idx].upper() != wt_aa.upper():
                        continue

                    header = f"{gene_name}-{mutant_clean}"
                    mutant_sequences[header] = update_str(aa_sequence, mut_aa, idx)

    except Exception as exc:
        print(f"Warning: Failed to synthesize mutants for {gene_name} ({mutation_file}): {exc}")
        return {}

    return mutant_sequences


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="NetSurfP-3.0 pipeline for protein structure prediction with WT/mutant comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('input', nargs='?',
                        help='Input: WT FASTA file or directory of FASTA files (nucleotide or amino acid)')
    parser.add_argument('output', nargs='?',
                        help='Output path: base filename (e.g., "results" → results.summary.tsv) or directory (uses "nsp3_output" as base name)')

    # Input type
    parser.add_argument('--input-type', choices=['nt', 'aa'], default='nt',
                        help='Input sequence type: "nt" for nucleotide (will translate), "aa" for amino acid (default: nt)')

    # Processing options
    parser.add_argument('--mutation-dir', required=True,
                        help='Mutation file or directory. For --input-type=nt: NT mutations (e.g., A1002T). For --input-type=aa: AA mutations (e.g., M334V) (REQUIRED)')
    parser.add_argument('--model', required=True,
                        help='Path to trained NSP3 model checkpoint (REQUIRED)')
    parser.add_argument('--config', required=True,
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
        print(f"Input type: {args.input_type} ({'nucleotide - will translate' if args.input_type == 'nt' else 'amino acid - no translation'})")
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

        # Handle based on input type
        if args.input_type == 'nt':
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
            # Basic validation - check for non-AA characters
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
            if not all(c.upper() in valid_aa for c in wt_aa_seq):
                print(f"Warning: {gene_name} contains non-amino acid characters, skipping")
                continue

        if args.verbose:
            print(f"  WT sequence length: {len(wt_aa_seq)} AA")

        # Build mutant sequences
        mutation_file = mutation_files.get(gene_name)
        if not mutation_file:
            print(f"Warning: No mutation file for {gene_name}, skipping")
            continue

        mutant_seqs = build_mutant_sequences_for_gene(
            gene_name, wt_nt_seq, wt_aa_seq, mutation_file, args.log, failure_map,
            input_type=args.input_type
        )

        if args.verbose:
            print(f"  Generated {len(mutant_seqs)} mutant sequences")

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
                continue

            wt_pred = all_predictions[wt_key]

            # Process each mutant
            for mut_id, mut_seq in mutant_seqs.items():
                if mut_id not in all_predictions:
                    print(f"Warning: Predictions not found for {mut_id}, skipping")
                    continue

                mut_pred = all_predictions[mut_id]

                # Extract mutation info from ID (e.g., "ABCB1-A1002T")
                mutation_str = mut_id.split('-', 1)[1] if '-' in mut_id else mut_id
                nt_info = get_mutation_data_bioAccurate(mutation_str)
                if nt_info[0] is None:
                    continue

                aa_info = get_mutant_aa(nt_info, wt_nt_seq)
                if not aa_info:
                    continue

                (mutation_pos, (wt_aa, mut_aa)), _ = aa_info
                mutation_pos = int(mutation_pos)

                # Get predictions at mutation site
                if mutation_pos not in wt_pred or mutation_pos not in mut_pred:
                    continue

                wt_res = wt_pred[mutation_pos]
                mut_res = mut_pred[mutation_pos]

                # Compute deltas
                delta_rsa = mut_res['rsa'] - wt_res['rsa']
                delta_disorder_pf = mut_res['disorder_pf'] - wt_res['disorder_pf']
                delta_disorder_pt = mut_res['disorder_pt'] - wt_res['disorder_pt']
                # Binary: 0 if same, 1 if different
                ss3_change = 0 if wt_res['q3'] == mut_res['q3'] else 1
                ss8_change = 0 if wt_res['q8'] == mut_res['q8'] else 1

                # Classifications
                burial_class = classify_burial_change(wt_res['rsa'], mut_res['rsa'])
                disorder_class = classify_disorder_change(wt_res['disorder_pf'], mut_res['disorder_pf'])

                # QC flags
                qc_flags = compute_qc_flags(
                    wt_aa, wt_res['residue'],
                    mut_aa, mut_res['residue'],
                    mut_res['phi'], mut_res['psi']
                )

                # Local structural impact (±5 residues)
                local_window = 5
                local_impact = 0.0
                local_changes_list = []

                for offset in range(-local_window, local_window + 1):
                    pos = mutation_pos + offset
                    if pos in wt_pred and pos in mut_pred:
                        wt_local = wt_pred[pos]
                        mut_local = mut_pred[pos]

                        delta_local_rsa = mut_local['rsa'] - wt_local['rsa']
                        delta_local_disorder_pf = mut_local['disorder_pf'] - wt_local['disorder_pf']
                        delta_local_disorder_pt = mut_local['disorder_pt'] - wt_local['disorder_pt']
                        delta_local_phi = mut_local['phi'] - wt_local['phi']
                        delta_local_psi = mut_local['psi'] - wt_local['psi']

                        local_impact += abs(delta_local_rsa)

                        local_changes_list.append({
                            'pkey': mut_id,
                            'gene': gene_name,
                            'relative_pos': offset,
                            'absolute_pos': pos,
                            'delta_rsa': delta_local_rsa,
                            'delta_disorder_pf': delta_local_disorder_pf,
                            'delta_disorder_pt': delta_local_disorder_pt,
                            'delta_phi': delta_local_phi,
                            'delta_psi': delta_local_psi,
                            'wt_ss3': wt_local['q3'],
                            'mut_ss3': mut_local['q3'],
                            'wt_ss8': wt_local['q8'],
                            'mut_ss8': mut_local['q8'],
                            'is_mutation_site': (offset == 0)
                        })

                local_rows.extend(local_changes_list)

                # Global metrics
                all_positions = set(wt_pred.keys()) & set(mut_pred.keys())
                total_delta_rsa = sum(
                    abs(mut_pred[pos]['rsa'] - wt_pred[pos]['rsa'])
                    for pos in all_positions
                )
                mean_delta_rsa = total_delta_rsa / len(all_positions) if all_positions else 0.0

                ss3_changes = sum(
                    1 for pos in all_positions
                    if wt_pred[pos]['q3'] != mut_pred[pos]['q3']
                )

                max_abs_delta_rsa = max(
                    abs(mut_pred[pos]['rsa'] - wt_pred[pos]['rsa'])
                    for pos in all_positions
                ) if all_positions else 0.0

                max_abs_delta_disorder = max(
                    abs(mut_pred[pos]['disorder_pf'] - wt_pred[pos]['disorder_pf'])
                    for pos in all_positions
                ) if all_positions else 0.0

                # Summary row
                summary_rows.append({
                    'pkey': mut_id,
                    'gene': gene_name,
                    'mutation_pos': mutation_pos,
                    'wt_aa': wt_aa,
                    'mut_aa': mut_aa,
                    'delta_rsa': delta_rsa,
                    'delta_disorder_pf': delta_disorder_pf,
                    'delta_disorder_pt': delta_disorder_pt,
                    'ss3_change': ss3_change,
                    'ss8_change': ss8_change,
                    'burial_classification': burial_class,
                    'disorder_classification': disorder_class,
                    'local_structural_impact': local_impact,
                    'global_mean_delta_rsa': mean_delta_rsa,
                    'global_ss_changes': ss3_changes,
                    'max_abs_delta_rsa': max_abs_delta_rsa,
                    'max_abs_delta_disorder': max_abs_delta_disorder,
                    'qc_flags': qc_flags
                })

                # Per-residue rows for WT
                for pos, res_data in wt_pred.items():
                    residues_rows.append({
                        'pkey': mut_id,
                        'gene': gene_name,
                        'allele': 'wt',
                        'pos': pos,
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
                        'q3_c': res_data['q3_c']
                    })

                # Per-residue rows for mutant
                for pos, res_data in mut_pred.items():
                    residues_rows.append({
                        'pkey': mut_id,
                        'gene': gene_name,
                        'allele': 'mut',
                        'pos': pos,
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
                        'q3_c': res_data['q3_c']
                    })

        finally:
            # Clean up temp file
            if os.path.exists(combined_fasta.name):
                os.unlink(combined_fasta.name)

    # Write output TSVs
    output_path = Path(args.output)

    # Handle directory vs file base name
    if output_path.is_dir() or str(args.output) == '.':
        # Extract gene name from input path using utility function
        base_name = extract_gene_from_filename(args.input)

        # Fallback to nsp3_out if extraction fails (empty or looks like a path component)
        if not base_name or base_name in ('orf', 'fasta', 'sequences', 'input', 'data'):
            base_name = 'nsp3_out'

        output_path = output_path / base_name

        # Avoid overwriting: increment suffix if files exist
        candidate = output_path
        n = 1
        while (candidate.with_suffix('.summary.tsv').exists() or
               candidate.with_suffix('.residues.tsv').exists() or
               candidate.with_suffix('.local.tsv').exists()):
            candidate = output_path.parent / f"{base_name}_{n}"
            n += 1
        output_path = candidate

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary_path = output_path.with_suffix('.summary.tsv')
    residues_path = output_path.with_suffix('.residues.tsv')
    local_path = output_path.with_suffix('.local.tsv')

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