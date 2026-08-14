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
SpliceAI VCF Parser
Parses SpliceAI-annotated VCF files and extracts splice predictions with pkey mapping.
Uses transcript-level mutation mappings for nucleotide-level analysis.
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict
from pathlib import Path

from biofeaturefactory.utils.utility import _collect_failures_from_logs, load_mapping, write_tsv

LABELS = ['A', 'B', 'C', 'D']


def parse_spliceai_entries(info_field):
    """Return unique SpliceAI transcript blocks with block labels."""
    if not info_field.startswith('SpliceAI='):
        return []

    entries = info_field[len('SpliceAI='):].split(',')
    grouped: "OrderedDict[tuple, list[dict]]" = OrderedDict()

    for raw in entries:
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split('|')
        if len(parts) != 10:
            continue
        try:
            entry = {
                'allele': parts[0],
                'symbol': parts[1],
                'ds_ag': float(parts[2]),
                'ds_al': float(parts[3]),
                'ds_dg': float(parts[4]),
                'ds_dl': float(parts[5]),
                'dp_ag': int(parts[6]),
                'dp_al': int(parts[7]),
                'dp_dg': int(parts[8]),
                'dp_dl': int(parts[9]),
            }
        except (ValueError, IndexError):
            continue

        key = (
            entry['allele'],
            entry['symbol'],
            entry['ds_ag'],
            entry['ds_al'],
            entry['ds_dg'],
            entry['ds_dl'],
            entry['dp_ag'],
            entry['dp_al'],
            entry['dp_dg'],
            entry['dp_dl'],
        )
        grouped.setdefault(key, []).append(entry)

    results = []
    label_idx = 0
    for key, items in grouped.items():
        entry = items[0]
        if len(items) > 1 or label_idx >= len(LABELS):
            entry['block_label'] = 'dup'
        else:
            entry['block_label'] = LABELS[label_idx]
            label_idx += 1
        results.append(entry)

    return results


def generate_pkey_with_mapping(pos, ref, alt, gene_context, chromosome_mapping, transcript_mapping, skip_mutations=None):
    """
    Generate pkey using dual mapping approach.

    Args:
        pos (str): Position from VCF file
        ref (str): Reference allele from VCF file
        alt (str): Alternate allele from VCF file
        gene_context (str): Gene name from VCF header
        chromosome_mapping (dict): Transcript mutation ID to chromosome notation mapping
        transcript_mapping (dict): Transcript mutation to ORF mutation mapping
        skip_mutations (set): Optional set of ORF mutations to skip based on logs

    Returns:
        str: Generated pkey
    """
    # Step 1: Construct chromosome notation (REF + POS + ALT)
    chromosome_notation = f"{ref}{pos}{alt}"
    skip_upper = {m.upper() for m in skip_mutations} if skip_mutations else None

    # Step 2: Find the transcript mutation ID that maps to this chromosome notation
    candidates = [
        mutant for mutant, cn in chromosome_mapping.items()
        if cn == chromosome_notation and not (skip_upper and mutant.upper() in skip_upper)
    ]
    if not candidates:
        return None
    # F43: deterministic pick. A genomic REF+POS+ALT maps to multiple transcript keys
    # for a multi-isoform gene; the old dict-order `break` made the chosen ORF label
    # arbitrary. sorted()[0] is stable; the single-match case is byte-identical.
    transcript_mutation_id = sorted(candidates)[0]
    if len(candidates) > 1:
        print(f"[WARN] {chromosome_notation}: {len(candidates)} transcript keys share this "
              f"genomic notation; ORF label from '{transcript_mutation_id}'", file=sys.stderr)

    # Step 3: Look up ORF mutation ID
    orf_mutation_id = transcript_mapping.get(transcript_mutation_id)

    # Generate pkey
    if orf_mutation_id and gene_context:
        return f"{gene_context}-{orf_mutation_id}"

    return None

def parse_vcf_header(vcf_file):
    """
    Extract gene context and other metadata from VCF header.

    Args:
        vcf_file (str): Path to VCF file

    Returns:
        dict: Header metadata
    """
    metadata = {}

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                if '=' in line:
                    key_value = line[2:].strip().split('=', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        metadata[key] = value
            elif line.startswith('#CHROM'):
                break

    return metadata

def parse_spliceai_vcf(
    vcf_file,
    output_file,
    chromosome_mapping_file=None,
    transcript_mapping_file=None,
    threshold=0.0,
    failure_log=None,
):
    """
    Parse SpliceAI VCF file and extract predictions above threshold.

    Args:
        vcf_file (str): Path to input VCF file
        output_file (str): Path to output TSV file
        chromosome_mapping_file (str): Path to chromosome mapping file
        transcript_mapping_file (str): Path to transcript mapping file
        threshold (float): Minimum delta score threshold
        failure_log (str): Optional validation log used to filter failed mutations

    Returns:
        tuple: (success_bool, processed_count, predictions_count)
    """
    try:
        # Parse header for metadata
        metadata = parse_vcf_header(vcf_file)
        gene_context = metadata.get('gene_context', '').strip()
        if not gene_context:
            print(f"[WARN] {vcf_file}: ##gene_context header missing; every pkey will be None -> 0 rows", file=sys.stderr)

        # Load mapping files
        chromosome_mapping = {}
        transcript_mapping = {}
        #inverted_cmap = {}

        failures = _collect_failures_from_logs(failure_log) if failure_log else {}

        if chromosome_mapping_file and os.path.exists(chromosome_mapping_file):
            chromosome_mapping = load_mapping(chromosome_mapping_file, mapType='chromosome')
            print(f"Loaded {len(chromosome_mapping)} chromosome mappings")

        if transcript_mapping_file and os.path.exists(transcript_mapping_file):
            transcript_mapping = load_mapping(transcript_mapping_file)
            print(f"Loaded {len(transcript_mapping)} transcript mappings")
        else:
            print("Warning: No transcript mapping file found, using dual mapping approach")

        predictions = []
        processed_count = 0
        dropped_no_pkey = 0
        # F42: split the drop counter. A drop because the mutation is on the
        # validation skip-list is EXPECTED; a drop because the genomic notation has
        # no chromosome_mapping entry is a contract break. Reporting them together
        # made a total mapping failure look like routine filtering.
        dropped_skipped = 0
        dropped_unmapped = 0

        skip_set = failures.get(gene_context.upper(), set()) if failures and gene_context else set()
        # O(1) membership test for "did this genomic notation map to anything at all"
        mapped_notations = set(chromosome_mapping.values())

        with open(vcf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip header lines
                if line.startswith('#'):
                    continue

                if not line:
                    continue

                # Parse VCF line
                fields = line.split('\t')
                if len(fields) < 8:
                    print(f"Warning: Malformed VCF line {line_num}: insufficient fields", file=sys.stderr)
                    continue

                chrom, pos, variant_id, ref, alt, qual, filter_field, info = fields[:8]
                processed_count += 1

                # Parse SpliceAI INFO field
                spliceai_calls = parse_spliceai_entries(info)
                if not spliceai_calls:
                    continue

                for call in spliceai_calls:
                    max_delta_score = max(
                        call['ds_ag'],
                        call['ds_al'],
                        call['ds_dg'],
                        call['ds_dl']
                    )

                    if max_delta_score < threshold:
                        continue

                    pkey = generate_pkey_with_mapping(
                        pos,
                        ref,
                        alt,
                        gene_context,
                        chromosome_mapping,
                        transcript_mapping,
                        skip_set,
                    )

                    if not pkey:
                        dropped_no_pkey += 1
                        if f"{ref}{pos}{alt}" in mapped_notations:
                            dropped_skipped += 1      # matched a mapping key; skip-listed or no ORF id
                        else:
                            dropped_unmapped += 1     # genomic notation absent from chromosome_mapping
                        continue

                    prediction = {
                        'pkey': pkey,
                        # F44 REVERTED 2026-08-12: `gene` MUST stay the per-block SpliceAI
                        # SYMBOL. SpliceAI emits one block per annotation record spanning the
                        # position (spliceai/utils.py get_delta_scores -> genes[i]), and 2,349
                        # of 35,078 production variants (6.70%) lie inside >=2 distinct #NAME
                        # spans (CDSN 615/615, CFTR 466/1495, PKD1 277/2915). Overwriting this
                        # with gene_context made every block of such a locus claim the target
                        # gene and left block_label as the only difference. `symbol` is also
                        # part of the row-identity key in parse_spliceai_entries.
                        'gene': call['symbol'],
                        'gene_context': gene_context,   # pipeline namespace (matches the pkey prefix)
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'allele': call['allele'],
                        'block_label': call['block_label'],
                        'ds_ag': call['ds_ag'],
                        'ds_al': call['ds_al'],
                        'ds_dg': call['ds_dg'],
                        'ds_dl': call['ds_dl'],
                        'dp_ag': call['dp_ag'],
                        'dp_al': call['dp_al'],
                        'dp_dg': call['dp_dg'],
                        'dp_dl': call['dp_dl'],
                        'max_delta_score': max_delta_score,
                    }
                    predictions.append(prediction)

        # Write TSV output
        fieldnames = [
            'pkey', 'gene', 'gene_context', 'chrom', 'pos', 'ref', 'alt', 'allele', 'block_label',
            'ds_ag', 'ds_al', 'ds_dg', 'ds_dl',
            'dp_ag', 'dp_al', 'dp_dg', 'dp_dl',
            'max_delta_score'
        ]
        write_tsv(predictions, output_file, fieldnames)

        if dropped_unmapped:
            print(f"[WARN] {vcf_file}: {dropped_unmapped} above-threshold predictions dropped because the "
                  f"genomic notation is absent from chromosome_mapping (contract break)", file=sys.stderr)
        if dropped_skipped:
            print(f"[INFO] {vcf_file}: {dropped_skipped} above-threshold predictions dropped as expected "
                  f"(skip-listed mutation or no ORF id)", file=sys.stderr)
        if processed_count and not predictions and dropped_unmapped:
            print(f"[ERROR] {vcf_file}: {processed_count} variants processed but ZERO rows written and "
                  f"{dropped_unmapped} unmapped — the mapping file almost certainly does not match this VCF.",
                  file=sys.stderr)

        print(f"Processed {processed_count} variants, found {len(predictions)} predictions above threshold {threshold}")
        return True, processed_count, len(predictions)

    except FileNotFoundError:
        error_msg = f"Error: VCF file not found at '{vcf_file}'"
        print(error_msg, file=sys.stderr)
        return False, 0, 0
    except Exception as e:
        error_msg = f"An unexpected error occurred processing {vcf_file}: {e}"
        print(error_msg, file=sys.stderr)
        return False, 0, 0

def main():
    """
    Main function to drive the script from command line.
    """
    parser = argparse.ArgumentParser(
        description="Parse SpliceAI-annotated VCF files and extract splice predictions with transcript-level pkey mapping."
    )
    parser.add_argument("-i", "--input", required=True, help="Input SpliceAI VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("--chromosome-mapping", help="Path to chromosome mapping file (mutations/combined/chromosome/combined_GENE.csv)")
    parser.add_argument("--transcript-mapping", help="Path to transcript mapping file (mutations/combined/transcript/combined_GENE.csv)")
    parser.add_argument("--log", help="Path to validation log (or directory of logs) used to filter failed mutations")
    parser.add_argument("-t", "--threshold", type=float, default=0.0,
                        help="Minimum delta score threshold (default: 0.0)")

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parse VCF file
    success, processed_count, predictions_count = parse_spliceai_vcf(
        args.input,
        args.output,
        args.chromosome_mapping,
        args.transcript_mapping,
        args.threshold,
        args.log,
    )

    if success:
        print(f"[OK] Successfully parsed {args.input}")
        print(f"[OK] Output saved to {args.output}")
        print(f"[OK] {processed_count} variants processed, {predictions_count} predictions extracted")
    else:
        print(f"[X] Failed to parse {args.input}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
