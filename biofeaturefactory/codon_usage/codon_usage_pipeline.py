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

import argparse
import csv
import os
import sys
from pathlib import Path

from biofeaturefactory.utils.utility import (
    read_fasta,
    trim_muts,
    get_codon_counts,
    extract_codon_with_bicodons,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
    compute_cai,
    compute_tai,
    get_codon_tai,
    get_codon_cai_w,
    parse_variant,
    splice_seq,
    write_tsv,
    split_intronic_tokens,
    warn_intronic_unsupported,
)


# write_tsv is called with this list and extrasaction='ignore', so a row key absent
# here is dropped silently — add the column here whenever the row dict gains one.
FIELDNAMES = [
    'pkey',
    'Gene',
    'codon_number',
    'position_in_codon',
    'codon_wt',
    'codon_mut',
    'RSCU_wt',
    'RSCU_mut',
    'delta_RSCU',
    'W_wt',
    'W_mut',
    'delta_W',
    'CAI_W_wt',
    'CAI_W_mut',
    'delta_CAI_W',
    'tAI_wt',
    'tAI_mut',
    'delta_tAI',
    'CAI_gene',
    'tAI_gene',
    'bicodon_3prime_wt',
    'bicodon_3prime_mut',
    'RSCPU_3prime_wt',
    'RSCPU_3prime_mut',
    'CPS_3prime_wt',
    'CPS_3prime_mut',
    'delta_CPS_3prime',
    'noln_CPS_3prime_wt',
    'noln_CPS_3prime_mut',
    'W_CP_3prime_wt',
    'W_CP_3prime_mut',
    'bicodon_5prime_wt',
    'bicodon_5prime_mut',
    'RSCPU_5prime_wt',
    'RSCPU_5prime_mut',
    'CPS_5prime_wt',
    'CPS_5prime_mut',
    'delta_CPS_5prime',
    'noln_CPS_5prime_wt',
    'noln_CPS_5prime_mut',
    'W_CP_5prime_wt',
    'W_CP_5prime_mut',
    'bicodon_context',
    'qc_flags',
]


def _codon_context(seq, codon_idx0):
    """Return (codon, bicodon_3prime, bicodon_5prime) at a 0-based codon index.

    Each is None when the sequence does not supply a whole triplet there, which
    is the same condition the SNV path expresses through forward/reverse_bicodon.
    """
    if not seq or codon_idx0 < 0:
        return None, None, None
    start = codon_idx0 * 3
    codon = seq[start:start + 3]
    if len(codon) != 3:
        return None, None, None
    nxt = seq[start + 3:start + 6]
    prv = seq[start - 3:start] if start >= 3 else ''
    return (codon,
            codon + nxt if len(nxt) == 3 else None,
            prv + codon if len(prv) == 3 else None)


def _row_from_codons(gene, ntposnt, codon_number, position_in_codon,
                     wt_codon, wt_bi3, wt_bi5,
                     mut_codon, mut_bi3, mut_bi5,
                     codondata, codonpairdata, cai_gene, tai_gene, qc_flags):
    """Build the full 44-column row from codon/bicodon STRINGS.

    Every metric here is a table lookup keyed on a triplet or a hexamer. None of
    them needs a WT<->MUT positional correspondence, so once both alleles' codon
    context is named, every column is computable for ANY variant class --
    substitution, in-frame indel, or frameshift. Shared by the SNV and non-SNV
    paths so the two cannot drift apart.
    """
    def _delta(mut_val, wt_val):
        if mut_val is None or wt_val is None:
            return None
        return round(mut_val - wt_val, 6)

    def _cp(table, key):
        return codonpairdata[table].get(key) if key else None

    rscu_wt = codondata['RSCU'].get(wt_codon) if wt_codon else None
    rscu_mut = codondata['RSCU'].get(mut_codon) if mut_codon else None
    w_wt = codondata['W'].get(wt_codon) if wt_codon else None
    w_mut = codondata['W'].get(mut_codon) if mut_codon else None
    caiw_wt = get_codon_cai_w(wt_codon) if wt_codon else None
    caiw_mut = get_codon_cai_w(mut_codon) if mut_codon else None
    tai_wt = get_codon_tai(wt_codon) if wt_codon else None
    tai_mut = get_codon_tai(mut_codon) if mut_codon else None
    cps3_wt, cps3_mut = _cp('CPS', wt_bi3), _cp('CPS', mut_bi3)
    cps5_wt, cps5_mut = _cp('CPS', wt_bi5), _cp('CPS', mut_bi5)

    if mut_bi3 and not mut_bi5:
        bicodon_context = 'first_codon_3prime_only'
    elif mut_bi5 and not mut_bi3:
        bicodon_context = 'last_codon_5prime_only'
    elif mut_bi3 and mut_bi5:
        bicodon_context = 'middle_codon_both_directions'
    else:
        bicodon_context = 'insufficient_sequence'
        qc_flags = (qc_flags + ';NO_BICODON') if qc_flags else 'NO_BICODON'

    return {
        'pkey': f"{gene}-{ntposnt}",
        'Gene': gene,
        'codon_number': codon_number,
        'position_in_codon': position_in_codon,
        'codon_wt': wt_codon,
        'codon_mut': mut_codon,
        'RSCU_wt': rscu_wt,
        'RSCU_mut': rscu_mut,
        'delta_RSCU': _delta(rscu_mut, rscu_wt),
        'W_wt': w_wt,
        'W_mut': w_mut,
        'delta_W': _delta(w_mut, w_wt),
        'CAI_W_wt': caiw_wt,
        'CAI_W_mut': caiw_mut,
        'delta_CAI_W': _delta(caiw_mut, caiw_wt),
        'tAI_wt': tai_wt,
        'tAI_mut': tai_mut,
        'delta_tAI': _delta(tai_mut, tai_wt),
        'CAI_gene': cai_gene,
        'tAI_gene': tai_gene,
        'bicodon_3prime_wt': wt_bi3,
        'bicodon_3prime_mut': mut_bi3,
        'RSCPU_3prime_wt': _cp('RSCPU', wt_bi3),
        'RSCPU_3prime_mut': _cp('RSCPU', mut_bi3),
        'CPS_3prime_wt': cps3_wt,
        'CPS_3prime_mut': cps3_mut,
        'delta_CPS_3prime': _delta(cps3_mut, cps3_wt),
        'noln_CPS_3prime_wt': _cp('noln CPS', wt_bi3),
        'noln_CPS_3prime_mut': _cp('noln CPS', mut_bi3),
        'W_CP_3prime_wt': _cp('W_CP', wt_bi3),
        'W_CP_3prime_mut': _cp('W_CP', mut_bi3),
        'bicodon_5prime_wt': wt_bi5,
        'bicodon_5prime_mut': mut_bi5,
        'RSCPU_5prime_wt': _cp('RSCPU', wt_bi5),
        'RSCPU_5prime_mut': _cp('RSCPU', mut_bi5),
        'CPS_5prime_wt': cps5_wt,
        'CPS_5prime_mut': cps5_mut,
        'delta_CPS_5prime': _delta(cps5_mut, cps5_wt),
        'noln_CPS_5prime_wt': _cp('noln CPS', wt_bi5),
        'noln_CPS_5prime_mut': _cp('noln CPS', mut_bi5),
        'W_CP_5prime_wt': _cp('W_CP', wt_bi5),
        'W_CP_5prime_mut': _cp('W_CP', mut_bi5),
        'bicodon_context': bicodon_context,
        'qc_flags': qc_flags or 'PASS',
    }


def _change_start0(variant):
    """0-based index of the FIRST BASE THE EDIT ACTUALLY CHANGES.

    NOT variant.pos0. Tokens are written in anchored (VCF) form, where a pure
    indel carries a retained anchor base and the change begins AFTER it: deleting
    codon 50 of TP53 is written TATT147T, whose pos0 is 146 -- the last base of
    codon 49. Every codon-frame question this module asks is about the changed
    bases, so asking it of the anchor answers about the wrong codon.

    Measured before this existed: `pos0 % 3` classified a clean whole-codon
    deletion as INDEL_NOT_CODON_ALIGNED and a boundary-straddling one as aligned
    -- inverted in all 9 anchor offsets sampled -- and codon_number named the
    codon BEFORE the deletion, so the row reported codon_wt == codon_mut and
    delta 0.0 for a codon that had been removed outright.

    The common prefix is the general rule: 0 for an MNV or delins (which change
    from their first base), 1 for the usual single-anchor indel, and correct for
    any longer shared prefix a producer happens to emit.
    """
    ref, alt = variant.ref, variant.alt
    k = 0
    while k < min(len(ref), len(alt)) and ref[k].upper() == alt[k].upper():
        k += 1
    return variant.pos0 + k


def _non_snv_row(gene, ntposnt, variant, wt_sequence, codondata, codonpairdata,
                 cai_gene, tai_gene, consequence):
    """Full row for a non-SNV token -- every column computed, none blanked.

    The mutant codon at a given position is whatever occupies that position AFTER
    the edit. That is well defined for a deletion (the following codon slides in),
    an insertion (the inserted codon), and a frameshift (the first re-framed
    codon). So codon and bicodon lookups work for all of them, and delta_* is a
    genuine "what sits here before vs after" comparison rather than an
    approximation.

    `consequence` records the SCOPE of the change, which is what the row cannot
    express: for a frameshift this row describes the first affected codon while
    every downstream codon has also changed identity.
    """
    change0 = _change_start0(variant)
    codon_idx0 = change0 // 3
    # The mutant side is unavailable when the token is being reported precisely
    # BECAUSE it does not fit the ORF (REF mismatch, span past the end). Emit the
    # WT context anyway rather than dropping the row -- the WT codon at that
    # position is a real, correct value regardless of what the token claims.
    try:
        mut_sequence = splice_seq(wt_sequence, variant.pos0, variant.ref,
                                  variant.alt, validate=False)
    except ValueError:
        mut_sequence = None
    wt_codon, wt_bi3, wt_bi5 = _codon_context(wt_sequence, codon_idx0)
    mut_codon, mut_bi3, mut_bi5 = _codon_context(mut_sequence, codon_idx0)

    # CAI_gene / tAI_gene are whole-ORF indices. Copying the WT value onto a
    # frameshift row is wrong: a frameshift re-reads every downstream codon, so
    # the mutant ORF's gene-level indices genuinely differ. Recompute them on the
    # mutant ORF and report the shift, rather than silently reusing the WT number.
    gene_cai, gene_tai = cai_gene, tai_gene
    if mut_sequence and variant.length_delta % 3 != 0:
        try:
            # Same call shape the pipeline uses for the WT gene indices.
            gene_cai = compute_cai(mut_sequence)
            gene_tai = compute_tai(mut_sequence)
            consequence = f"{consequence};gene_indices_recomputed_on_mutant_orf"
        except Exception as exc:
            # Name the failure. A bare 'failed' marker hides which call broke and
            # leaves the WT indices in place looking like real mutant values.
            consequence = (f"{consequence};gene_index_recompute_failed:"
                           f"{type(exc).__name__}")

    return _row_from_codons(
        gene, ntposnt, codon_idx0 + 1, (change0 % 3) + 1,
        wt_codon, wt_bi3, wt_bi5, mut_codon, mut_bi3, mut_bi5,
        codondata, codonpairdata, gene_cai, gene_tai, consequence)


def process_mutation(gene, ntposnt, sequence, codondata, codonpairdata, cai_gene, tai_gene,
                     wt_sequence):
    """
    Process a single mutation and return codon usage statistics.

    Args:
        gene: Gene symbol
        ntposnt: Mutation string (e.g., "A123G")
        sequence: WILD-TYPE ORF nucleotide sequence
        codondata: Pre-computed codon statistics dict
        codonpairdata: Pre-computed codon pair statistics dict
        cai_gene: Pre-computed CAI for the gene
        tai_gene: Pre-computed tAI for the gene
        wt_sequence: The wild-type ORF. REQUIRED -- both alleles are derived from
            it, so there is no sensible default. It used to default to None for
            the removed --mutant mode, which handed in an already-mutated
            sequence; that mode silently dropped every SNV, because
            extract_codon_with_bicodons REF-guards the token against whatever
            sequence it is given and a mutated one carries the ALT.

    Returns:
        dict: Row data with codon usage metrics, or None if processing fails
    """
    # Non-SNV gate. extract_codon_with_bicodons parses the token itself and
    # raises on a multi-base one, so the frame question has to be settled first.
    #
    # A codon-usage row is one codon. That survives an edit only when whole
    # codons are removed or added at a codon boundary:
    #   frameshift (len_delta % 3 != 0) -- every downstream codon changes identity;
    #                                      one row cannot express it.
    #   not codon-aligned (pos0 % 3 != 0) -- a 3 bp deletion straddling a boundary
    #                                      fuses two codons into one, so
    #                                      `mutated_codon` has no referent.
    # Both are refused BY NAME rather than approximated.
    variant = parse_variant(ntposnt, is_nt=True)
    if variant is not None and not variant.is_snv:
        # A bad token gets a NAMED row, not a silent None. Returning None here
        # drops the mutation from the output entirely, which is indistinguishable
        # from "this gene had no such mutation" -- the exact failure mode the
        # reason codes exist to prevent.
        ref_seq = wt_sequence
        if variant.pos0 + len(variant.ref) > len(ref_seq):
            return _non_snv_row(gene, ntposnt, variant, ref_seq, codondata,
                                codonpairdata, cai_gene, tai_gene,
                                'REF_SPANS_PAST_ORF')
        observed = ref_seq[variant.pos0:variant.pos0 + len(variant.ref)].upper()
        if observed != variant.ref.upper():
            return _non_snv_row(gene, ntposnt, variant, ref_seq, codondata,
                                codonpairdata, cai_gene, tai_gene,
                                f'REF_MISMATCH:orf_has_{observed}')
        if variant.length_delta % 3 != 0:
            consequence = 'FRAMESHIFT:downstream_codons_also_change'
        elif variant.length_delta == 0:
            consequence = 'MNV'
        elif _change_start0(variant) % 3 != 0:
            consequence = ('INDEL_NOT_CODON_ALIGNED:'
                           + ('deletion' if variant.length_delta < 0 else 'insertion'))
        else:
            consequence = 'CODON_DELETED' if variant.length_delta < 0 else 'CODON_INSERTED'
        return _non_snv_row(gene, ntposnt, variant, ref_seq, codondata,
                            codonpairdata, cai_gene, tai_gene, consequence)

    result = extract_codon_with_bicodons(ntposnt, sequence)
    if result[0] is None:
        return None

    mutated_codon, forward_bicodon, reverse_bicodon, poc, pos, codon_number = result

    pkey = f"{gene}-{ntposnt}"
    qc_flags = []

    # Determine bicodon context
    if codon_number == 1 and forward_bicodon and not reverse_bicodon:
        bicodon_context = 'first_codon_3prime_only'
    elif not forward_bicodon and reverse_bicodon:
        bicodon_context = 'last_codon_5prime_only'
    elif forward_bicodon and reverse_bicodon:
        bicodon_context = 'middle_codon_both_directions'
    else:
        bicodon_context = 'insufficient_sequence'
        qc_flags.append('NO_BICODON')

    # WT counterpart. extract_codon_with_bicodons splices ALT into a copy, so the codon
    # it returns is the mutant; the wild-type triplet is the same slice of the unmutated
    # ORF. Safe because that function's REF guard already proved wt_sequence[pos] == REF.
    # Both alleles are looked up in the SAME codondata/codonpairdata, which are built once
    # from the WT ORF — a point mutation must not shift the gene's reference usage table.
    # wt_sequence is None when the caller only has an already-mutated sequence
    # (process_mutant_fasta): the WT is unrecoverable there, so the wt/delta columns are
    # left empty rather than filled with the mutant's own values, which would read as a
    # delta of 0.0 for every row.
    wt_codon = wt_bi3 = wt_bi5 = None
    if wt_sequence:
        start = (codon_number - 1) * 3
        wt_codon = wt_sequence[start:start + 3]
        if forward_bicodon:
            wt_bi3 = wt_codon + forward_bicodon[3:]
        if reverse_bicodon:
            wt_bi5 = reverse_bicodon[:3] + wt_codon
        if len(wt_codon) != 3:
            wt_codon = wt_bi3 = wt_bi5 = None
            qc_flags.append('NO_WT_CODON')

    def _delta(mut_val, wt_val):
        if mut_val is None or wt_val is None:
            return None
        return round(mut_val - wt_val, 6)

    # SNV path routed through the same builder the non-SNV path uses, so the two
    # cannot drift. Every column below is a codon/bicodon table lookup.
    return _row_from_codons(
        gene, ntposnt, codon_number, poc + 1,
        wt_codon, wt_bi3, wt_bi5,
        mutated_codon, forward_bicodon or None, reverse_bicodon or None,
        codondata, codonpairdata, cai_gene, tai_gene,
        ';'.join(qc_flags) if qc_flags else '')



def process_fasta_with_mutations(fasta_path, mutations_path, validation_log=None):
    """
    Process a FASTA file with mutations from a CSV file.

    Args:
        fasta_path: Path to FASTA file containing ORF sequence
        mutations_path: Path to CSV file containing mutations
        validation_log: Optional path to validation log for filtering

    Returns:
        list: List of row dictionaries
    """
    results = []
    fasta = read_fasta(fasta_path)

    if 'ORF' not in fasta:
        print(f"Warning: No 'ORF' found in {fasta_path}", file=sys.stderr)
        return results

    sequence = fasta['ORF']
    gene = extract_gene_from_filename(fasta_path)

    # Load validation failures if provided
    failure_map = load_validation_failures(validation_log) if validation_log else {}

    # Load mutations
    mut_list = trim_muts(mutations_path, log=validation_log, gene_name=gene)

    if not mut_list:
        print(f"Warning: No mutations found in {mutations_path}", file=sys.stderr)
        return results

    # Intronic gate. Codon usage is frame-dependent: every column in FIELDNAMES
    # is a codon or bicodon table lookup, and an intron has no reading frame, so
    # none of them has a defined value for an intronic variant.
    #
    # Without this gate the tokens do NOT reach a codon lookup anyway -- they
    # raise inside extract_codon_with_bicodons and are caught by the F48 handler
    # below -- but that handler reports them as "skipped malformed token", which
    # is false. 'gd.T5000C' is a well-formed token in a coordinate space this
    # pipeline cannot use, and calling it malformed sends the reader looking for
    # a corrupt input file.
    mut_list, intronic = split_intronic_tokens(mut_list)
    warn_intronic_unsupported(
        'codon_usage', gene, intronic,
        "Codon usage requires a reading frame; an intron has none. "
        "Score these with RNAfold or miranda instead.")

    if not mut_list:
        print(f"Warning: {gene}: every mutation was intronic; no codon rows to write",
              file=sys.stderr)
        return results

    # Pre-compute codon statistics once per sequence
    codondata, codonpairdata = get_codon_counts(sequence)

    # Compute gene-level CAI and tAI
    cai_gene = compute_cai(sequence)
    tai_gene = compute_tai(sequence)

    for ntposnt in mut_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        # F48: one malformed token must not discard the gene. get_mutation_data_bioAccurate
        # does int(ntposnt[1:-1]), and the is_nt guard only inspects the first/last chars,
        # so a token with nucleotide ends but a broken middle (c.76A>T, A76>T) still raises
        # ValueError. Unguarded, it propagated out of this loop and every mutation already
        # collected in `results` was lost with it.
        try:
            row = process_mutation(gene, ntposnt, sequence, codondata, codonpairdata,
                                   cai_gene, tai_gene, wt_sequence=sequence)
        except (ValueError, IndexError) as e:
            print(f"[codon_usage] {gene}-{ntposnt}: skipped malformed token ({e})", file=sys.stderr)
            continue
        if row:
            results.append(row)

    return results


def process_directory(fasta_dir, mutations_dir=None, validation_log=None, output_dir=None):
    """
    Process all FASTA files in a directory, writing per-gene output files.

    Args:
        fasta_dir: Directory containing FASTA files
        mutations_dir: Directory containing mutation CSV files (for WT mode)
        validation_log: Optional path to validation log
        output_dir: Base output directory for per-gene nested output

    Returns:
        list: Combined results from all files (for summary printing)
    """
    results = []
    fasta_extensions = ['*.fasta', '*.fa', '*.fna']

    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(Path(fasta_dir).glob(ext))

    print(f"Found {len(fasta_files)} FASTA file(s) to process")

    for fasta_file in sorted(fasta_files):
        print(f"Processing {fasta_file}...")

        gene = extract_gene_from_filename(str(fasta_file))

        mutations_file = None

        if mutations_dir:
            gene_up = gene.upper()
            # F50: exact stem match over a SORTED glob. The old substring test
            # (`gene.upper() in stem.upper()`) over an unsorted glob with a
            # first-match break bound a prefix gene (F9) to a superset gene's
            # file (F9A), nondeterministically across machines — unlike the
            # sorted() fasta list above. Accepts `<GENE>.csv` and the
            # `<GENE>_mutations.csv` convention used in Bio_DBs/mappings/mutations.
            for csv_file in sorted(Path(mutations_dir).glob('*.csv')):
                # Use the shared extractor rather than hardcoded stems: it resolves
                # 59/59 production files AND the other layouts this repo uses
                # (combined_<gene>.csv per main.nf:108-112, muts_<gene>.csv, <gene>.csv),
                # and keeps hyphenated symbols intact (HLA-A, NKX2-1).
                if (extract_gene_from_filename(csv_file.name) or '').upper() == gene_up:
                    mutations_file = str(csv_file)
                    break

        if not mutations_file:
            print(f"Warning: No mutations file found for {gene}", file=sys.stderr)
            continue

        file_results = process_fasta_with_mutations(
            str(fasta_file),
            mutations_file,
            validation_log
        )

        print(f"  Processed {len(file_results)} mutations")

        if output_dir and file_results:
            out_dir = Path(output_dir) / gene / "CodonUsage"
            out_dir.mkdir(parents=True, exist_ok=True)
            write_output(file_results, str(out_dir / f"{gene}.codon_usage.tsv"))

        results.extend(file_results)

    return results


def write_output(results, output_path):
    """Write results to TSV file."""
    if not results:
        print("No results to write")
        return

    write_tsv(results, output_path, FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(results)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Codon Usage Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file processing
  python codon_usage_pipeline.py --fasta /path/to/gene.fasta --mutations /path/to/mutations.csv --output results/

  # Directory processing
  python codon_usage_pipeline.py --fasta /path/to/fastas --mutations /path/to/mutations --output results/

Metrics:
  RSCU       - Relative Synonymous Codon Usage (gene-specific)
  W          - Relative adaptiveness within gene
  CAI_W      - Reference W value for CAI (human highly expressed genes)
  tAI        - tRNA Adaptation Index weight for the codon
  CAI_gene   - Codon Adaptation Index for entire gene
  tAI_gene   - tRNA Adaptation Index for entire gene
  RSCPU      - Relative Synonymous Codon Pair Usage
  CPS        - Codon Pair Score (ln of observed/expected)
"""
    )

    # Input options
    parser.add_argument('--fasta', required=True, help='FASTA file or directory of FASTA files')
    parser.add_argument('--mutations', help='Mutations CSV file or directory of CSV files')
    parser.add_argument('--validation-log', help='Validation log for filtering failed mutations')

    # Output options
    parser.add_argument('--output', '-o', required=True, help='Output base directory')

    args = parser.parse_args()

    # Validate arguments
    if not Path(args.fasta).is_dir() and not args.mutations:
        parser.error("--mutations required when using a single FASTA file")
    if Path(args.fasta).is_dir() and not args.mutations:
        parser.error("--mutations required when using a FASTA directory")

    # Process files
    if Path(args.fasta).is_dir():
        results = process_directory(
            args.fasta,
            args.mutations,
            args.validation_log,
            output_dir=args.output,
        )
    else:
        results = process_fasta_with_mutations(
            args.fasta,
            args.mutations,
            args.validation_log
        )

        # Write per-gene output for single-file mode
        gene = extract_gene_from_filename(args.fasta)
        out_dir = Path(args.output) / gene / "CodonUsage"
        out_dir.mkdir(parents=True, exist_ok=True)
        write_output(results, str(out_dir / f"{gene}.codon_usage.tsv"))

    # Summary
    if results:
        n_with_3prime = sum(1 for r in results if r.get('bicodon_3prime'))
        n_with_5prime = sum(1 for r in results if r.get('bicodon_5prime'))
        n_with_both = sum(1 for r in results if r.get('bicodon_3prime') and r.get('bicodon_5prime'))

        # F51: in directory mode `results` spans EVERY gene (results.extend per file),
        # so results[0] is just whichever gene sorted first — printing its CAI/tAI as
        # "Gene CAI/tAI" mislabels one gene's value as global. Report per gene instead.
        # (The per-gene TSVs were always correct; only this stdout line was wrong.)
        per_gene = {}
        for r in results:
            g = r.get('Gene')
            if g is not None and g not in per_gene:
                per_gene[g] = (r.get('CAI_gene'), r.get('tAI_gene'))

        print(f"\nSummary:")
        print(f"  Total mutations processed: {len(results)}")
        print(f"  With 3' bicodon: {n_with_3prime}")
        print(f"  With 5' bicodon: {n_with_5prime}")
        print(f"  With both bicodons: {n_with_both}")
        if len(per_gene) == 1:
            gene_only, (cai, tai) = next(iter(per_gene.items()))
            if cai:
                print(f"  Gene CAI ({gene_only}): {cai:.4f}")
            if tai:
                print(f"  Gene tAI ({gene_only}): {tai:.4f}")
        elif per_gene:
            print(f"  Genes processed: {len(per_gene)} "
                  f"(per-gene CAI/tAI are in each {{GENE}}.codon_usage.tsv)")


if __name__ == "__main__":
    main()
