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
BioFeatureFactory: Dtu Outputs

Parsing and recombination of DTU tool output (NetNGlyc, NetPhos), including
per-mutation prediction filtering. Consumers: netNglyc and netphos.

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
    mint_pkey,
    ExtractGeneFromFASTA,
    extract_mutation_from_sequence_name,
    get_mutation_data_bioAccurate,
    should_skip_mutation,
)


def _combine_glycosylation_outputs(batch_output_files, final_output_file, original_fasta_file=None):
    """Combine glycosylation prediction batch outputs (NetNGlyc format)"""
    try:
        # Count total sequences for header
        total_sequences = 0
        all_sequence_sections = []
        all_prediction_lines = []

        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()


                import os
                seperator = ['-','_']
                # Use original FASTA file if provided, otherwise fallback to counting Name: lines
                if original_fasta_file and os.path.exists(original_fasta_file):
                    try:
                        #fasta_sequences = read_fasta(original_fasta_file)
                        gene_name, total_sequences = ExtractGeneFromFASTA(original_fasta_file, count=True)
                        # Determine if this is wildtype or mutant based on sequence names
                            # For the first batch, use total sequences from original file
                        if i == 0:
                            print(f"Using original FASTA file for sequence count: {total_sequences} sequences from {gene_name}")
                            # For subsequent batches, the total was previously counted

                    except Exception as e:
                        print(f"Warning: Error reading original FASTA file {original_fasta_file}: {e}")
                        # Fallback to counting Name: lines in NetNGlyc output
                        lines = content.split('\n')
                        name_count = sum(1 for line in lines if line.startswith('Name:'))
                        total_sequences += name_count if name_count > 0 else 1

                else:
                    # Fallback: Count sequences directly from NetNGlyc output
                    lines = content.split('\n')
                    name_count = sum(1 for line in lines if line.startswith('Name:'))

                    if name_count > 0:
                        total_sequences += name_count
                    else:
                        # Parse header line for sequence count: ">debug-GENE-aa-netnglyc\t5 amino acids"
                        for line in lines:
                            if 'amino acids' in line:
                                try:
                                    parts = line.split()
                                    for j, part in enumerate(parts):
                                        if part.isdigit() and j+1 < len(parts) and 'amino' in parts[j+1]:
                                            total_sequences += int(part)
                                            break
                                    break
                                except:
                                    total_sequences += 1  # Ultimate fallback
                            else:
                                total_sequences += 1  # Fallback if no header found

                # Collect sequence display sections and prediction lines
                in_sequence_section = False
                in_prediction_section = False

                lines = content.split('\n')
                for line in lines:
                    if line.strip().startswith('>') and 'amino acids' in line:
                        continue  # Skip individual headers

                    if 'SeqName' in line and 'Position' in line:
                        in_prediction_section = True
                        continue

                    if line.startswith('    ') and len(line.strip()) > 10:
                        in_sequence_section = True
                        all_sequence_sections.append(line)
                    elif in_prediction_section and line.strip() and not line.startswith('#'):
                        all_prediction_lines.append(line)

            except Exception as e:
                print(f"Error reading batch file {batch_file}: {e}")
                continue

        # Write combined output
        with open(final_output_file, 'w') as f:
            # Write header with total sequence count
            f.write(f">{os.path.basename(final_output_file).replace('.out', '')}\t{total_sequences} amino acids\n\n")

            # Write prediction header
            f.write("SeqName                 Position  Potential  N-Glyc result  Comment\n")
            f.write("=" * 70 + "\n")

            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")

            # Write sequence sections
            if all_sequence_sections:
                f.write("\n")
                for line in all_sequence_sections:
                    f.write(line + "\n")

        print(f"Combined {len(batch_output_files)} batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")

        return True

    except Exception as e:
        print(f"Error combining glycosylation outputs: {e}")
        return False


def _combine_phosphorylation_outputs(batch_output_files, final_output_file):
    """Combine phosphorylation prediction batch outputs (NetPhos format)"""
    try:
        total_sequences = 0
        all_prediction_lines = []

        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()

                # Extract sequences count from header
                lines = content.split('\n')
                for line in lines:
                    if 'amino acids' in line and line.startswith('>'):
                        try:
                            seq_count = int(line.split()[1])
                            total_sequences += seq_count
                        except:
                            total_sequences += 1  # Fallback
                        break

                # Collect prediction lines
                for line in lines:
                    if line.startswith('# ') and len(line.split()) >= 7:
                        # This is a prediction line
                        all_prediction_lines.append(line)

            except Exception as e:
                print(f"Error reading phosphorylation batch file {batch_file}: {e}")
                continue

        # Write combined phosphorylation output
        with open(final_output_file, 'w') as f:
            # Write header
            gene_name = os.path.basename(final_output_file).replace('.out', '').replace('-netphos', '')
            f.write(f">{gene_name}\t{total_sequences} amino acids\n")
            f.write("#\n")
            f.write("#  prediction results\n")
            f.write("#\n")
            f.write("# Sequence\t\t   # x   Context     Score   Kinase    Answer\n")
            f.write("# " + "-" * 67 + "\n")

            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")

        print(f"Combined {len(batch_output_files)} phosphorylation batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")

        return True

    except Exception as e:
        print(f"Error combining phosphorylation outputs: {e}")
        return False


def process_single_mutation_for_sequence(seq_name, predictions, mapping_dict, is_mutant=True, tool_type='netphos', failure_map=None):
    """Process predictions for one sequence against its specific mutation

    Args:
        seq_name: Sequence name (e.g., 'ZFP36-C330T')
        predictions: List of predictions for this sequence
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant sequences (should be True for single-mutation processing)
        tool_type: 'netphos' or 'netnglyc' for tool-specific field handling

    Returns:
        list: Filtered predictions with pkeys for the specific mutation
    """
    if not is_mutant:
        raise ValueError("process_single_mutation_for_sequence should only be used for mutant sequences")

    # Extract mutation ID from sequence name
    gene, mutation_id = extract_mutation_from_sequence_name(seq_name)
    if mutation_id is None:
        return []

    # Sequence names are {GENE}-{sha} now, so what comes back is the DIGEST, not
    # the token -- and mapping_dict is keyed by the verbatim token, so the lookup
    # below missed on every sequence and this function returned [] for all of
    # them. No new parameter is needed: mapping_dict already holds the tokens, so
    # mint each and match the digest. A {GENE}-{token} name still hits the direct
    # lookup first and never reaches this.
    if mutation_id not in mapping_dict:
        for _tok in mapping_dict:
            if mint_pkey(gene, _tok).rsplit('-', 1)[-1] == mutation_id:
                mutation_id = _tok
                break

    # Look up this mutation in the mapping
    if mutation_id not in mapping_dict:
        return []

    if should_skip_mutation(gene, mutation_id, failure_map):
        return []

    aaposaa = mapping_dict[mutation_id]  # e.g., "Y110F"

    # Parse amino acid position and mutation info
    position_data = get_mutation_data_bioAccurate(aaposaa, is_nt=False)
    if position_data[0] is None:
        return []

    aa_pos = position_data[0]  # e.g., 110
    aa_tuple = position_data[1]  # e.g., ('Y', 'F')
    target_aa = aa_tuple[1]  # F for mutant

    # Filter predictions for this specific mutation
    results = []
    for pred in predictions:
        # Get position from prediction
        if tool_type == 'netphos':
            pred_pos = pred['pos']
            pred_aa = pred['amino_acid']
        elif tool_type == 'netnglyc':
            pred_pos = pred['position']
            pred_aa = pred['sequon'][0] if pred['sequon'] else None
        else:
            raise ValueError(f"Unsupported tool_type: {tool_type}")

        # NetNGlyc is a DELTA tool and keeps every mutant site; netphos keeps only
        # the site AT the mutated residue.
        #
        # The shared condition below asks "is there a predicted site exactly at the
        # mutated position, on the mutant residue". For netphos that is the
        # question. For netnglyc it is unanswerable by construction: pred_aa is
        # sequon[0], which is ALWAYS 'N', so `pred_aa == target_aa` can only pass
        # for a mutation TO asparagine. MEASURED on PAM: D563G (aa_pos 563,
        # target_aa 'G') and S539W (539, 'W') against predictions at 411 and 762 --
        # no position match, and 'N' == 'G' / 'N' == 'W' never holds. Every mutant
        # returned [], surfacing as missing_mut + no_comparable_site, and the
        # count_gained / count_lost / delta columns could never populate.
        #
        # What netnglyc actually needs is ALL mutant sites, so the WT<->MUT
        # alignment downstream can pair them and compute deltas -- a variant that
        # weakens a site 200 residues away is precisely what the delta columns are
        # for. netphos and netMHC keep the original behaviour.
        keep = (pred_pos == aa_pos and pred_aa == target_aa)
        if tool_type == 'netnglyc':
            keep = True
        if keep:
            # Create pkey for this match. mutation_id is the VERBATIM token by
            # this point (resolved above), so minting reproduces exactly the pkey
            # variant_mapping wrote.
            pkey = mint_pkey(gene, mutation_id)

            # Add pkey and fix Gene field to prediction
            result_pred = pred.copy()
            result_pred['pkey'] = pkey
            # Carry the token explicitly. The filter below used to recover it
            # with pkey.split('-')[-1], which is right only while the pkey ends
            # in the token; once the pkey is a digest that returns the digest and
            # should_skip_mutation silently stops matching the validation log, so
            # nothing is ever skipped. It is in scope here -- no parsing needed.
            result_pred['mutation_id'] = mutation_id
            # Fix Gene field to just gene name (not gene-mutation)
            result_pred['Gene'] = gene

            # Map field names to match CSV writer expectations
            if tool_type == 'netnglyc':
                # Map NetNGlyc field names to CSV format
                if 'position' in result_pred:
                    result_pred['pos'] = result_pred.pop('position')
                if 'sequon' in result_pred:
                    result_pred['Sequon'] = result_pred.pop('sequon')
                # Remove seq_name as it's not needed in CSV
                if 'seq_name' in result_pred:
                    result_pred.pop('seq_name')

            results.append(result_pred)

    return results


def parse_predictions_with_mutation_filtering(predictions, mapping_dict, is_mutant, threshold=0.0, yes_only=False, tool_type='netphos', failure_map=None):
    """Universal prediction filtering logic for both NetPhos and NetNGlyc

    Args:
        predictions: List of all predictions
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant (True) or wildtype (False) sequences
        threshold: Score threshold for filtering
        yes_only: Only include predictions with 'YES' answer (NetPhos only)
        tool_type: 'netphos' or 'netnglyc' for tool-specific handling

    Returns:
        list: Filtered predictions with pkeys
    """
    results = []

    failure_map = failure_map or {}

    if is_mutant:
        # Group predictions by sequence name for single-mutation processing
        seq_predictions = {}
        for pred in predictions:
            seq_name = pred.get('Gene', '') if tool_type == 'netphos' else pred.get('seq_name', '')
            if seq_name not in seq_predictions:
                seq_predictions[seq_name] = []
            seq_predictions[seq_name].append(pred)

        # Process each sequence separately with its specific mutation
        for seq_name, seq_preds in seq_predictions.items():
            seq_results = process_single_mutation_for_sequence(
                seq_name, seq_preds, mapping_dict, is_mutant=True, tool_type=tool_type,
                failure_map=failure_map
            )

            # Apply additional filters
            for result in seq_results:
                # Apply threshold filter
                score_field = 'score' if tool_type == 'netphos' else 'potential'
                if result[score_field] < threshold:
                    continue

                if should_skip_mutation(result['Gene'], result.get('mutation_id'), failure_map):
                    continue

                # Apply yes_only filter (NetPhos only)
                if yes_only and tool_type == 'netphos' and result.get('answer') != 'YES':
                    continue

                results.append(result)

    else:
        # Wildtype processing - use existing bulk logic (not implemented here)
        # This should use the existing wildtype processing logic from each pipeline
        raise NotImplementedError("Wildtype processing should use existing pipeline-specific logic")

    return results

