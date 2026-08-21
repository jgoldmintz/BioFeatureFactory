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

from biofeaturefactory.utils.utility import (
    _collect_failures_from_logs,
    Variant,
    canonical_token,
    load_mapping,
    parse_variant,
    write_tsv,
)

LABELS = ['A', 'B', 'C', 'D']

# Characters an allele may contain and still be a real sequence. Anything else
# (<DEL>, *, .) is symbolic and is left to exact string matching rather than
# normalized, because a symbolic allele has no parsimonious form.
_NT_CHARS = frozenset("ACGTUacgtu")


def spliceai_info_value(info_field):
    """Return the SpliceAI INFO value, or None when the record carries none.

    INFO is a ';'-delimited list of key=value fields, so the SpliceAI field is
    not necessarily the first one. Matching only `info_field.startswith(...)`
    discarded every record of any VCF whose INFO carries anything else --
    'AC=1;AN=2;SpliceAI=...' is what the SpliceAI CLI emits when it annotates a
    VCF that already had INFO keys, and main.nf accepts such a VCF directly via
    --input_vcf_file/--input_vcf_dir. Splitting first is byte-identical for the
    BFF-generated VCFs (vcf_converter writes INFO='.', so the annotated field is
    the only one) and additionally correct for externally annotated input.

    The SpliceAI value itself is comma-delimited and contains no ';', so a
    single split on ';' cannot cut a block in half.
    """
    for field in info_field.split(';'):
        if field.startswith('SpliceAI='):
            return field[len('SpliceAI='):]
    return None


def parse_spliceai_entries(info_field, stats=None):
    """Return unique SpliceAI transcript blocks with block labels.

    stats is an optional dict the caller owns. Blocks this function refuses are
    counted into it under 'malformed_blocks' instead of vanishing: a block with
    the wrong arity or an unparseable score used to be `continue`d silently, so
    a whole VCF of subtly corrupt INFO fields reported zero predictions and zero
    problems. Passing nothing keeps the old single-argument call shape.
    """
    value = spliceai_info_value(info_field)
    if value is None:
        return []

    entries = value.split(',')
    grouped: "OrderedDict[tuple, list[dict]]" = OrderedDict()

    for raw in entries:
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split('|')
        if len(parts) != 10:
            if stats is not None:
                stats['malformed_blocks'] = stats.get('malformed_blocks', 0) + 1
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
            if stats is not None:
                stats['malformed_blocks'] = stats.get('malformed_blocks', 0) + 1
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


def as_variant(pos, ref, alt):
    """Build a length-aware Variant from VCF fields, or None if not one.

    Returns None for a symbolic or absent allele (<DEL>, *, .) and for a POS that
    is not a positive integer, rather than guessing: those records still join by
    exact string equality, they just cannot be normalized or classified.
    """
    if not ref or not alt:
        return None
    if not (set(ref) <= _NT_CHARS and set(alt) <= _NT_CHARS):
        return None
    try:
        pos_i = int(pos)
    except (TypeError, ValueError):
        return None
    if pos_i < 1:
        return None
    return Variant(pos=pos_i, ref=ref, alt=alt)


def build_notation_index(chromosome_mapping):
    """canonical genomic notation -> list of transcript mutation IDs.

    Built ONCE per VCF by parse_spliceai_vcf, because pkey resolution runs once
    per (variant x SpliceAI block) and rebuilding this per call is O(mapping) on
    every one of them.

    The index exists because the join at the heart of this parser is exact string
    equality on {REF}{POS}{ALT}. An SNV has exactly one spelling, so exact
    equality is total for it. A multi-base allele does not: ACG100ACT and
    ATGC100ATTC are both the substitution G102T carrying redundant shared flank,
    and acaa112217430a is ACAA112217430A. canonical_token uppercases and
    parsimoniously trims, collapsing all of those onto one key, so a mapping and
    a VCF that disagree only on how much flank they carry still join.
    (Verified by execution: canonical_token maps ACG100ACT and ATGC100ATTC to
    G102T, and leaves anchor-form pure indels such as ACAA112217430A untouched.)

    LIMIT, stated because it is invisible otherwise: canonical_token(v) without a
    sequence normalizes REPRESENTATION, not PLACEMENT. Two left-alignment
    spellings of the same indel inside a tandem repeat remain distinct keys, and
    this parser is handed no reference FASTA with which to left-align them.
    Repeat-adjacent indels therefore still depend on producer and consumer
    agreeing on placement.
    """
    index = {}
    for mutant, cn in chromosome_mapping.items():
        v = parse_variant(cn, is_nt=True)
        if v is None:
            continue
        index.setdefault(canonical_token(v), []).append(mutant)
    return index


def resolve_pkey(pos, ref, alt, gene_context, chromosome_mapping, transcript_mapping,
                 skip_mutations=None, notation_index=None):
    """Resolve a VCF record to a BFF pkey AND name the outcome.

    Returns (pkey_or_None, mode). mode is one of:
        'exact'                 matched a chromosome_mapping value verbatim
        'canonical'             matched only after representation normalization
        'no_chromosome_match'   genomic notation absent from chromosome_mapping
        'skip_listed'           matched, but every candidate is on the skip-list
        'no_orf_id'             matched, but transcript_mapping has no ORF label
        'no_gene_context'       resolved, but the VCF header named no gene

    The mode is the point: a dropped record used to be indistinguishable from a
    record that was never there. Every caller-visible drop now carries a reason.
    """
    # Step 1: Construct chromosome notation (REF + POS + ALT)
    chromosome_notation = f"{ref}{pos}{alt}"
    skip_upper = {m.upper() for m in skip_mutations} if skip_mutations else None

    # Step 2: Find the transcript mutation ID that maps to this chromosome notation.
    # Exact equality is tried FIRST and alone decides the outcome whenever it
    # matches, so every record that resolved before resolves identically now.
    candidates = [
        mutant for mutant, cn in chromosome_mapping.items()
        if cn == chromosome_notation
    ]
    mode = 'exact'
    if not candidates:
        variant = as_variant(pos, ref, alt)
        if variant is not None:
            if notation_index is None:
                notation_index = build_notation_index(chromosome_mapping)
            candidates = list(notation_index.get(canonical_token(variant), ()))
            mode = 'canonical'
    if not candidates:
        return None, 'no_chromosome_match'

    candidates = [m for m in candidates
                  if not (skip_upper and m.upper() in skip_upper)]
    if not candidates:
        return None, 'skip_listed'

    # F43: deterministic pick. A genomic REF+POS+ALT maps to multiple transcript keys
    # for a multi-isoform gene; the old dict-order `break` made the chosen ORF label
    # arbitrary. sorted()[0] is stable; the single-match case is byte-identical.
    transcript_mutation_id = sorted(candidates)[0]
    if len(candidates) > 1:
        print(f"[WARN] {chromosome_notation}: {len(candidates)} transcript keys share this "
              f"genomic notation; ORF label from '{transcript_mutation_id}'", file=sys.stderr)

    # Step 3: Look up ORF mutation ID
    orf_mutation_id = transcript_mapping.get(transcript_mutation_id)
    if not orf_mutation_id:
        return None, 'no_orf_id'
    if not gene_context:
        return None, 'no_gene_context'

    return f"{gene_context}-{orf_mutation_id}", mode


def generate_pkey_with_mapping(pos, ref, alt, gene_context, chromosome_mapping, transcript_mapping, skip_mutations=None):
    """
    Generate pkey using dual mapping approach.

    Args:
        pos (str): Position from VCF file
        ref (str): Reference allele from VCF file (may be multi-base)
        alt (str): Alternate allele from VCF file (may be multi-base)
        gene_context (str): Gene name from VCF header
        chromosome_mapping (dict): Transcript mutation ID to chromosome notation mapping
        transcript_mapping (dict): Transcript mutation to ORF mutation mapping
        skip_mutations (set): Optional set of ORF mutations to skip based on logs

    Returns:
        str: Generated pkey, or None if the record does not resolve.

    Use resolve_pkey directly when the REASON for a None matters; this wrapper
    keeps the pkey-only contract for callers that only need the key.
    """
    return resolve_pkey(pos, ref, alt, gene_context, chromosome_mapping,
                        transcript_mapping, skip_mutations)[0]


def variant_qc_flags(record_ref, record_alt, allele, match_mode):
    """qc_flags cell for one emitted row.

    Names, per row, every fact about the record that a consumer of DS_*/DP_*
    would otherwise have to rediscover:

      * whether the row is an SNV, and if not, which class and by how many nt --
        SpliceAI scores an indel by collapsing the mutant into the reference
        frame (deletions zero-padded, insertions max-collapsed), so DP_* offsets
        on a non-SNV row are reference-frame offsets from POS, not mutant-frame
        ones. A consumer that does not know the row is an indel will read them
        as if it were.
      * whether this row is one allele of a multi-allelic record, and which.
      * whether the pkey was reached by exact or normalized matching.

    'PASS' means plain bi-allelic SNV, exact key match -- the only case that
    needed no qualification before this column existed.
    """
    flags = []
    variant = as_variant('1', record_ref, allele)   # pos is irrelevant to kind
    if variant is None:
        flags.append('SYMBOLIC_ALLELE')
    elif not variant.is_snv:
        flags.append(f"NON_SNV:{variant.kind}")
        # Only when the length actually changed. An MNV is a non-SNV with
        # length_delta 0, and printing 'length_delta:+0nt' on it reads as a
        # measured null rather than as "this class does not change length".
        if variant.length_delta:
            flags.append(f"length_delta:{variant.length_delta:+d}nt")
    if allele != record_alt:
        flags.append(f"MULTIALLELIC:allele_{allele}_of_{record_alt}")
    if match_mode == 'canonical':
        flags.append('PKEY_MATCH:canonical')
    return ';'.join(flags) if flags else 'PASS'

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
        # F42: split the drop counter. A drop because the mutation is on the
        # validation skip-list is EXPECTED; a drop because the genomic notation has
        # no chromosome_mapping entry is a contract break. Reporting them together
        # made a total mapping failure look like routine filtering.
        #
        # The split is now taken from resolve_pkey's own verdict instead of being
        # reconstructed afterwards from a set-membership test, so 'matched but had
        # no ORF label' and 'matched but was skip-listed' stop being reported as
        # one bucket. dropped_non_snv is tracked in parallel because the failure
        # this work exists to catch -- a mapping that carries only SNVs -- shows up
        # as an all-non-SNV drop bucket and is invisible in the total.
        dropped = {}
        dropped_non_snv = {}
        dropped_examples = {}
        block_stats = {}
        malformed_lines = 0

        skip_set = failures.get(gene_context.upper(), set()) if failures and gene_context else set()
        # Canonical index of the mapping's genomic notations, built once per VCF
        # and consulted only when exact matching fails (see build_notation_index).
        notation_index = build_notation_index(chromosome_mapping)

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
                    # Counted as well as printed: a per-line warning scrolls past
                    # and leaves no trace in the totals, so a VCF that is wholly
                    # the wrong shape produced 0 rows and 0 processed with nothing
                    # in the summary saying why.
                    malformed_lines += 1
                    continue

                chrom, pos, variant_id, ref, alt, qual, filter_field, info = fields[:8]
                processed_count += 1

                # Parse SpliceAI INFO field
                spliceai_calls = parse_spliceai_entries(info, block_stats)
                if not spliceai_calls:
                    # NOT a non-event. SpliceAI declines to score whole classes of
                    # record (anything outside an annotated gene span, and indels
                    # it cannot collapse into the reference frame), and it says so
                    # by emitting no INFO field rather than a zero-score block. A
                    # bare `continue` here therefore dropped exactly the records
                    # this non-SNV work exists to carry, and the summary line
                    # ("Processed N variants, found M predictions") gave no way to
                    # tell an unannotated record from one that was never in the
                    # VCF. Reason codes are record-level; the pkey reasons below
                    # are per (record x block).
                    reason = ('spliceai_blocks_all_refused'
                              if spliceai_info_value(info) is not None
                              else 'no_spliceai_annotation')
                    dropped[reason] = dropped.get(reason, 0) + 1
                    # ALT may list several alleles; the record is non-SNV if any
                    # one of them is, because a single ALT column is scored as a
                    # set and the whole record is what was dropped.
                    if any((v is not None and not v.is_snv)
                           for v in (as_variant(pos, ref, a) for a in alt.split(','))):
                        dropped_non_snv[reason] = dropped_non_snv.get(reason, 0) + 1
                    examples = dropped_examples.setdefault(reason, [])
                    if len(examples) < 5:
                        examples.append(f"{ref}{pos}{alt}")
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

                    # Mint the key from THIS BLOCK'S allele, not from the record's
                    # raw ALT column. SpliceAI emits one block per (alt, gene)
                    # pair, so on a multi-allelic record the ALT column is
                    # 'A,AT' while each block scores exactly one of them; keying
                    # on the column produced the notation 'ATG100A,AT', which
                    # matches nothing and dropped every allele of every
                    # multi-allelic record. For a bi-allelic record
                    # call['allele'] IS the ALT column, so nothing changes there.
                    allele = call['allele']
                    pkey, mode = resolve_pkey(
                        pos,
                        ref,
                        allele,
                        gene_context,
                        chromosome_mapping,
                        transcript_mapping,
                        skip_set,
                        notation_index,
                    )

                    if not pkey:
                        # Counted skip with a NAMED reason, per drop class, and
                        # separately for non-SNV records -- never a silent drop.
                        dropped[mode] = dropped.get(mode, 0) + 1
                        dropped_variant = as_variant(pos, ref, allele)
                        if dropped_variant is not None and not dropped_variant.is_snv:
                            dropped_non_snv[mode] = dropped_non_snv.get(mode, 0) + 1
                        examples = dropped_examples.setdefault(mode, [])
                        if len(examples) < 5:
                            examples.append(f"{ref}{pos}{allele}")
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
                        # Trailing column. Every column above is defined and
                        # populated for every variant class -- they come straight
                        # out of the SpliceAI block, which the tool computes for
                        # indels natively. What a non-SNV row needs is not a
                        # blanked value but a statement of what it is; that is
                        # this column.
                        'qc_flags': variant_qc_flags(ref, alt, allele, mode),
                    }
                    predictions.append(prediction)

        # Write TSV output
        fieldnames = [
            'pkey', 'gene', 'gene_context', 'chrom', 'pos', 'ref', 'alt', 'allele', 'block_label',
            'ds_ag', 'ds_al', 'ds_dg', 'ds_dl',
            'dp_ag', 'dp_al', 'dp_dg', 'dp_dl',
            'max_delta_score', 'qc_flags'
        ]
        write_tsv(predictions, output_file, fieldnames)

        dropped_unmapped = dropped.get('no_chromosome_match', 0)
        if dropped_unmapped:
            print(f"[WARN] {vcf_file}: {dropped_unmapped} above-threshold predictions dropped because the "
                  f"genomic notation is absent from chromosome_mapping (contract break)", file=sys.stderr)
        dropped_skipped = dropped.get('skip_listed', 0) + dropped.get('no_orf_id', 0)
        if dropped_skipped:
            print(f"[INFO] {vcf_file}: {dropped_skipped} above-threshold predictions dropped as expected "
                  f"(skip-listed mutation or no ORF id)", file=sys.stderr)
        if dropped.get('no_gene_context'):
            print(f"[WARN] {vcf_file}: {dropped['no_gene_context']} above-threshold predictions resolved to "
                  f"an ORF mutation but were dropped for want of a ##gene_context header", file=sys.stderr)
        # Itemize every drop class. dropped_examples always holds >=1 notation for
        # any reason present in `dropped` (both are written together), so this is
        # the full picture of what was refused and a sample of what it looked like.
        for reason in sorted(dropped):
            print(f"[INFO] {vcf_file}: dropped[{reason}]={dropped[reason]} "
                  f"(non-SNV {dropped_non_snv.get(reason, 0)}); "
                  f"e.g. {', '.join(dropped_examples[reason])}", file=sys.stderr)
        # Only meaningful when SOME row survived: if nothing resolved at all the
        # mapping simply does not match this VCF, which the [ERROR] below says.
        #
        # Scoped to the MAPPING drop reasons. dropped_non_snv also carries
        # record-level reasons now ('no SpliceAI INFO at all'), and those say
        # nothing about the mapping file — blaming it for them would send the
        # reader to the wrong file.
        mapping_non_snv = sum(dropped_non_snv.get(r, 0) for r in
                              ('no_chromosome_match', 'skip_listed', 'no_orf_id'))
        if predictions and mapping_non_snv and not any('NON_SNV' in p['qc_flags'] for p in predictions):
            print(f"[WARN] {vcf_file}: every non-SNV record was dropped "
                  f"({mapping_non_snv} of them) while SNV records resolved — the mapping "
                  f"file carries no indel entries.", file=sys.stderr)
        if block_stats.get('malformed_blocks'):
            print(f"[WARN] {vcf_file}: {block_stats['malformed_blocks']} SpliceAI INFO blocks refused "
                  f"(wrong arity or unparseable score) and are absent from the output", file=sys.stderr)
        if malformed_lines:
            print(f"[WARN] {vcf_file}: {malformed_lines} VCF lines refused for having fewer than 8 "
                  f"columns; they are absent from the {processed_count} processed", file=sys.stderr)
        unannotated = dropped.get('no_spliceai_annotation', 0)
        if unannotated:
            print(f"[WARN] {vcf_file}: {unannotated} of {processed_count} records carry no SpliceAI "
                  f"INFO field at all ({dropped_non_snv.get('no_spliceai_annotation', 0)} of them "
                  f"non-SNV) — SpliceAI scored no transcript for them", file=sys.stderr)
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
