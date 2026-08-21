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
Build exon- and intron-aware mappings + sequences for genes:
- FASTA per gene:
      ORF         CDS only, mRNA orientation
      transcript  spliced, mRNA orientation
      genomic     tx_start..tx_end, GENOMIC orientation (strand not applied)
      pre_mRNA    the same span in TRANSCRIPT orientation (strand applied)
      intron<N>|<g_start>-<g_end>   one per mutation-carrying intron,
                  transcript orientation
  genomic and pre_mRNA are the same bases and differ only by strand; both are
  kept because genesplicer takes genomic-forward sequence with a separate strand
  flag, while RNAfold, miranda and AlphaFold3 are all strand-blind and must be
  handed the molecule's own orientation.
- Mapping CSVs:
    * chromosome: mutant -> genomic-orientation "REF<abs_pos>ALT"
    * genomic (optional): mutant -> genomic-slice-relative "REF<gDNA_pos>ALT"
    * transcript (optional): mutant -> transcript-orientation "REF<tx_pos>ALT"
    * amino-acid (optional): mutant -> AA mutation "REF<aa_pos>ALT"
    * intron (optional): mutant -> "intron<N>:REF<intron_pos>ALT"
    * pre_mRNA (optional): mutant -> "REF<premrna_pos>ALT", transcript
      orientation, covering BOTH exonic and intronic variants

Input tokens are ORF-relative by default. An intronic (or otherwise
gDNA-relative) token carries the prefix 'gd.' -- e.g. 'gd.T5000C' -- because a
bare token is ambiguous between the two spaces whenever a gene's 5'UTR is
shorter than its ORF (SMN2: ORF 1-879, gDNA slice starting at 18). Intronic
tokens appear in the chromosome, genomic, intron and pre_mRNA CSVs; they get no
transcript row and no amino-acid row, because neither coordinate exists for
them and a zero-delta row would read as a measured 'no effect'.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path
import pysam
import warnings

from biofeaturefactory.utils.utility import (
    trim_muts,
    extract_gene_from_filename,
    get_genome_loc,
    read_fasta,
    write_fasta,
    parse_variant,
    revcomp_seq,
    infer_aavariant_from_nt,
    format_aa_token,
    INTRONIC_PREFIX,
    split_intronic_tokens,
    mint_pkey,
)


# Verbose logging helper

def emit_verbose(message: str, verbose: bool, collector: list[str] | None) -> None:
    if not verbose:
        return
    if collector is not None:
        collector.append(message)
    else:
        print(message)



# ORF helpers

def load_supplied_orfs(path: str) -> dict[str, str]:
    """Load ORF sequences from a directory or FASTA file.

    Returns dict mapping UPPERCASE gene symbol -> uppercase ORF sequence.
    """
    src = Path(path)
    if not src.exists():
        raise FileNotFoundError(f"ORF path '{path}' does not exist")

    fasta_exts = {".fa", ".fasta", ".fna", ".fas", ".fn"}
    orfs: dict[str, str] = {}

    def _record_orf(gene_key: str, seq: str, origin: str) -> None:
        key = gene_key.upper()
        if key in orfs:
            warnings.warn(
                f"Duplicate ORF sequence for {gene_key} from '{origin}', keeping the first one.",
                RuntimeWarning,
            )
            return
        orfs[key] = seq.upper()

    if src.is_dir():
        for fasta_file in sorted(src.glob("*")):
            if not fasta_file.is_file() or fasta_file.suffix.lower() not in fasta_exts:
                continue
            records = read_fasta(str(fasta_file))
            if not records:
                continue
            first_name, first_seq = next(iter(records.items()))
            gene_name = extract_gene_from_filename(first_name)
            if gene_name and gene_name.upper() in {"ORF", "TRANSCRIPT", "GENOMIC"}:
                gene_name = None
            if not gene_name:
                gene_name = extract_gene_from_filename(fasta_file.name)
            if not gene_name:
                warnings.warn(
                    f"Could not infer gene name for ORF file '{fasta_file.name}', skipping.",
                    RuntimeWarning,
                )
                continue
            _record_orf(gene_name, first_seq, fasta_file.name)
    else:
        records = read_fasta(str(src))
        if not records:
            raise ValueError(f"No sequences found in ORF FASTA '{path}'")
        for name, seq in records.items():
            gene_name = extract_gene_from_filename(name)
            if gene_name and gene_name.upper() in {"ORF", "TRANSCRIPT", "GENOMIC"}:
                gene_name = None
            if not gene_name:
                gene_name = extract_gene_from_filename(src.name)
            if not gene_name:
                warnings.warn(
                    f"Could not infer gene name for ORF record '{name}' in '{src.name}', skipping.",
                    RuntimeWarning,
                )
                continue
            _record_orf(gene_name, seq, src.name)

    return orfs


def parse_force_cds(value: str) -> tuple[str | None, dict[str, str]]:
    """Parse --force-cds argument.

    Returns:
        (global_transcript, per_gene_map)
        - If single accession: (accession, {})
        - If CSV file: (None, {GENE: accession, ...})
    """
    # Check if value matches RefSeq accession pattern (e.g., NM_022162.3, XM_123456.1)
    accession_pattern = re.compile(r'^[NX][MR]_\d+\.\d+$')
    if accession_pattern.match(value):
        return (value, {})

    # Otherwise treat as file path
    path = Path(value)
    if not path.exists():
        raise FileNotFoundError(
            f"--force-cds value '{value}' is neither a valid accession (e.g., NM_022162.3) "
            f"nor an existing file."
        )

    transcript_map: dict[str, str] = {}
    with open(path, 'r') as f:
        import csv
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"CSV file '{value}' has no header row.")
        # Normalize fieldnames to lowercase for flexible matching
        fieldnames_lower = [fn.lower() for fn in reader.fieldnames]
        gene_col = None
        tid_col = None
        for fn, fn_lower in zip(reader.fieldnames, fieldnames_lower):
            if fn_lower in ('gene', 'gene_name', 'gene_symbol'):
                gene_col = fn
            if fn_lower in ('transcript_id', 'transcript', 'accession'):
                tid_col = fn
        if gene_col is None or tid_col is None:
            raise ValueError(
                f"CSV file '{value}' must have columns for gene (gene/gene_name/gene_symbol) "
                f"and transcript ID (transcript_id/transcript/accession). Found: {reader.fieldnames}"
            )
        for row in reader:
            gene = row.get(gene_col, '').strip().upper()
            tid = row.get(tid_col, '').strip()
            if gene and tid:
                transcript_map[gene] = tid

    return (None, transcript_map)


STOP_CODONS = {"TAA", "TAG", "TGA"}


def collect_orf_mutation_positions(orf_mutations: list[str]) -> tuple[list[int], dict[int, list[str]]]:
    """Collect the ORF position each mutation must be covered by.

    The reported position is the LAST base of the REF span (pos + len(ref) - 1),
    not the first, because ORF-coverage is only satisfied when the whole REF
    span fits inside the ORF. For SNVs len(ref) == 1 so this is identical to the
    previous first-base behaviour.

    parse_variant never raises and returns None for anything it cannot decode
    (stop-codon tokens such as 'A4Stop', empty alleles such as 'ACG30'). The
    previous get_mutation_data_bioAccurate call raised ValueError on every
    non-SNV token, and that exception escaped all the way to the per-gene
    handler in main() -- a single indel discarded every other mutation in the
    gene and produced no mapping CSVs at all.
    """
    positions: list[int] = []
    pos_to_mut: dict[int, list[str]] = {}
    for mut in orf_mutations:
        variant = parse_variant(mut, is_nt=True)
        if variant is None:
            continue
        span_end = variant.pos + len(variant.ref) - 1
        positions.append(span_end)
        pos_to_mut.setdefault(span_end, []).append(mut)
    return positions, pos_to_mut


def locate_orf_in_transcript(orf_seq: str, transcript_seq: str, gene_name: str) -> int:
    idx = transcript_seq.upper().find(orf_seq.upper())
    if idx == -1:
        raise ValueError(
            f"Supplied ORF for {gene_name} does not align with the spliced transcript."
        )
    return idx


def validate_mutations_against_orf(
    gene_name: str,
    orf_seq: str,
    orf_mutations: list[str],
) -> list[str]:
    if not orf_mutations:
        return []
    seq = orf_seq.upper()
    # (token, message) rather than a bare message. main() counts failures by
    # summing these lists against nt_dropped, and a token that is BOTH a REF
    # mismatch and a mapping drop was counted twice: measured on a 1-token file,
    # "Mutations: 1 (passed: 0, failed: 2)". Carrying the token lets the caller
    # take a set union instead of a sum.
    issues: list[tuple[str, str]] = []
    mismatches: list[tuple[str, str]] = []
    for mut in orf_mutations:
        variant = parse_variant(mut, is_nt=True)
        if variant is None:
            continue
        ref = variant.ref.upper()
        idx = variant.pos - 1
        end = idx + len(ref)
        # The whole REF span must fit, not just its first base. A 4 nt deletion
        # anchored 2 nt before the ORF end is out of range even though its start
        # position is in range.
        if idx < 0 or end > len(seq):
            issues.append((
                mut,
                f"{gene_name}: mutation {mut} spans ORF positions {variant.pos}-{variant.pos + len(ref) - 1},"
                f" but ORF length is {len(seq)}.",
            ))
            continue
        observed = seq[idx:end]
        if observed != ref:
            mismatches.append((
                mut,
                f"{gene_name}: mutation {mut} expects '{ref}' at ORF position(s)"
                f" {variant.pos}-{variant.pos + len(ref) - 1}, but ORF sequence has '{observed}'.",
            ))
    return issues, mismatches


def derive_orf_from_transcript(
    gene_name: str,
    transcript_seq: str,
    mutation_positions: list[int],
    pos_to_mut: dict[int, list[str]],
    verbose: bool,
    log_messages: list[str] | None,
) -> tuple[str, int]:
    seq = transcript_seq.upper()
    best_candidate: tuple[bool, int, int] | None = None  # (covers_mutations, length, -start)
    best_range: tuple[int, int] | None = None

    for start in range(len(seq) - 2):
        if seq[start:start+3] != "ATG":
            continue
        stop_idx = None
        for idx in range(start, len(seq) - 2, 3):
            codon = seq[idx:idx+3]
            if codon in STOP_CODONS:
                stop_idx = idx
                break
        if stop_idx is None:
            continue
        orf_len = (stop_idx - start) + 3
        covers = all(pos <= orf_len for pos in mutation_positions) if mutation_positions else True
        priority = (covers, orf_len, -start)
        if best_candidate is None or priority > best_candidate:
            best_candidate = priority
            best_range = (start, stop_idx + 3)

    if best_range is None:
        raise ValueError(
            f"Could not derive an ORF for {gene_name}: no in-frame start/stop codon found."
        )

    start, end = best_range
    orf_seq = seq[start:end]
    if mutation_positions and not best_candidate[0]:
        uncovered = sorted({pos for pos in mutation_positions if pos > len(orf_seq)})
        if uncovered:
            emit_verbose(
                f"{gene_name}: derived ORF length {len(orf_seq)} cannot cover mutation position(s) {uncovered}.",
                verbose,
                log_messages,
            )
            for pos in uncovered:
                muts = pos_to_mut.get(pos, [])
                for mut in muts:
                    emit_verbose(
                        f"  - {mut}: ORF length {len(orf_seq)} < position {pos}",
                        verbose,
                        log_messages,
                    )
        # let validation step raise the final error after reporting specifics

    return orf_seq, start


def resolve_orf_from_sources(
    gene_name: str,
    transcript_seq: str,
    _annotation_info: dict,
    orf_mutations: list[str],
    supplied_orf: str | None,
    verbose: bool = False,
    log_messages: list[str] | None = None,
) -> tuple[int, int, str, list[str], list[str]]:
    mutation_positions, pos_to_mut = collect_orf_mutation_positions(orf_mutations)

    if supplied_orf:
        orf_seq = supplied_orf.upper()
        start_idx = locate_orf_in_transcript(orf_seq, transcript_seq, gene_name)
    else:
        try:
            orf_seq, start_idx = derive_orf_from_transcript(
                gene_name,
                transcript_seq,
                mutation_positions,
                pos_to_mut,
                verbose,
                log_messages,
            )
        except ValueError as err:
            raise ValueError(
                f"{gene_name}: {err} Provide --orf to supply an explicit ORF sequence."
            ) from err

    # validate_mutations_against_orf returns a BARE [] for an empty token list and
    # a 2-tuple otherwise -- a return type that depends on the input. That was
    # unreachable while main() skipped any gene with no mutations, but a gene whose
    # tokens are ALL intronic now arrives here with orf_mutations == [] after
    # split_intronic_tokens, and the unpack below raised
    # "not enough values to unpack (expected 2, got 0)" on the first real run.
    # Guarded at the call site rather than changing the function, whose bare-[]
    # contract is pinned by test_exon_mapping.py::test_empty_mutations_returns_empty.
    if orf_mutations:
        out_of_range, mismatches = validate_mutations_against_orf(gene_name, orf_seq, orf_mutations)
    else:
        out_of_range, mismatches = [], []
    if out_of_range or mismatches:
        if verbose:
            for _mut, msg in out_of_range:
                emit_verbose(msg, True, log_messages)
            for _mut, msg in mismatches:
                emit_verbose(msg, True, log_messages)
            total = len({m for m, _ in out_of_range} | {m for m, _ in mismatches})
            matched = max(len(orf_mutations) - total, 0)
            emit_verbose(
                f"{gene_name} validation summary:\n"
                f"  - Total mutations: {len(orf_mutations)}\n"
                f"  - Passed: {matched}\n"
                f"  - Failed: {total} (length: {len(out_of_range)}, base mismatch: {len(mismatches)})",
                True,
                log_messages,
            )
        if supplied_orf:
            raise ValueError(
                f"{gene_name}: ORF/mutation validation failed for {len(out_of_range) + len(mismatches)} mutation(s)."
            )

    cds_tx_start = start_idx + 1  # convert to 1-based
    cds_tx_end = cds_tx_start + len(orf_seq) - 1
    return cds_tx_start, cds_tx_end, orf_seq.upper(), out_of_range, mismatches


# Helpers

def rc_base(b: str) -> str:
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(b.upper(), b)

def _canonical_chrom_token(name: str | None) -> str | None:
    if not name:
        return None
    token = name.strip().lower()
    if not token:
        return None

    if token.startswith("chr"):
        token = token[3:]

    base = token.replace("_", "")
    base = base.split(".", 1)[0]

    if base.startswith("nc") and len(base) > 2:
        match = re.match(r"nc0*(\d+)$", base)
        if match:
            number = str(int(match.group(1)))
            if number == "23":
                return "x"
            if number == "24":
                return "y"
            if number == "12920":
                return "mt"
            return number
    if base in {"m", "mt"}:
        return "mt"
    if base in {"x", "y"}:
        return base
    if base.isdigit():
        number = str(int(base))
        if number == "23":
            return "x"
        if number == "24":
            return "y"
        return number

    return token


def pick_chr_name(fasta: pysam.FastaFile, chrom: str) -> str:
    # Try simple, UCSC, already prefixed
    cand = [chrom, f"chr{chrom}", chrom.replace("chr", "")]
    chrom_core = chrom.replace("chr", "").strip()
    chrom_upper = chrom_core.upper()
    if chrom_upper == "X":
        cand.extend(["23", "chr23"])
    elif chrom_upper == "Y":
        cand.extend(["24", "chr24"])
    elif chrom_upper in {"MT", "M"}:
        cand.extend(["MT", "chrM", "chrMT", "M", "12920", "chr12920"])
    canon = _canonical_chrom_token(chrom)
    if canon:
        cand.append(canon)
        cand.append(f"chr{canon}")
    candidates = [c for c in {c for c in cand if c}]
    for ref in fasta.references:
        if ref in candidates:
            return ref
    # fallback: compare canonical tokens derived from FASTA reference names
    ref_aliases = {}
    for ref in fasta.references:
        alias = _canonical_chrom_token(ref)
        if alias:
            ref_aliases.setdefault(alias, ref)
    for candidate in candidates:
        alias = _canonical_chrom_token(candidate)
        if alias and alias in ref_aliases:
            return ref_aliases[alias]
    raise ValueError(f"Chromosome '{chrom}' not found in FASTA. Available: {list(fasta.references)[:6]} ...")

def write_dict_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    """Write dict rows under an explicit header. Missing keys become empty cells.

    Separate from write_mapping_csv because the decomposed mappings are no longer
    two columns: a variant spanning a splice site has a genomic address AND a
    per-feature breakdown, and the `mutant` column now carries the user's
    verbatim token so the {GENE}-{token} pkey still resolves after decomposition.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def write_mapping_csv(path: Path, header_label: str, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mutant", header_label])
        w.writerows(rows)


# Core builders

def build_transcript_seq_and_map(
    gene_name: str,
    annotation_file: str,
    reference_fasta: str,
    orf_mutations: list[str],
    supplied_orf: str | None = None,
    transcript_id: str | None = None,
    verbose: bool = False,
    log_messages: list[str] | None = None,
):
    """
    Returns dict with:
      chrom, strand, tx_start, tx_end, exons (list[[start,end],...]),
      transcript_seq (mRNA orientation),
      tx_to_genome (list[int], 1-based genomic coord per transcript index),
      cds_tx_start (1-based), cds_tx_end (1-based, inclusive)
      orf_seq (CDS only, mRNA orientation)

    Args:
        transcript_id: Optional specific transcript ID to force selection of
            (e.g., NM_022162.3). Overrides auto-selection if provided.
    """
    info = get_genome_loc(gene_name, annotation_file, transcript_id=transcript_id)
    if not info:
        raise ValueError(f"Gene '{gene_name}' not found in annotation '{annotation_file}'.")

    chrom      = info["chrom"]
    strand     = info["strand"]
    tx_start   = int(info["tx_start"])
    tx_end     = int(info["tx_end"])
    exons      = [(int(s), int(e)) for s, e in info["exons"]]

    fasta = pysam.FastaFile(reference_fasta)
    chr_name = pick_chr_name(fasta, chrom)

    # Build spliced transcript and tx->genome mapping (mRNA orientation)
    exon_blocks = sorted(exons, key=lambda x: x[0], reverse=(strand == "-"))
    tx_seq_chunks = []
    tx_to_genome = []

    if strand == "+":
        for s, e in exon_blocks:
            seg = fasta.fetch(chr_name, s-1, e).upper()
            tx_seq_chunks.append(seg)
            tx_to_genome.extend(range(s, e+1))
    else:
        for s, e in exon_blocks:
            seg = fasta.fetch(chr_name, s-1, e).upper()
            seg_rc = "".join(rc_base(b) for b in seg[::-1])
            tx_seq_chunks.append(seg_rc)
            tx_to_genome.extend(range(e, s-1, -1))

    transcript_seq = "".join(tx_seq_chunks)

    # Locate CDS boundaries in transcript space
    # cds_anchor approach: find exact genomic coords within tx_to_genome
    cds_tx_start, cds_tx_end, orf_seq, len_issues, mismatches = resolve_orf_from_sources(
        gene_name,
        transcript_seq,
        info,
        orf_mutations,
        supplied_orf,
        verbose=verbose,
        log_messages=log_messages,
    )

    # Genomic contiguous slice (genome orientation)
    gdna_seq = fasta.fetch(chr_name, tx_start-1, tx_end).upper()
    fasta.close()

    return {
        "chrom": chrom,
        "chr_name": chr_name,
        "strand": strand,
        "tx_start": tx_start,
        "tx_end": tx_end,
        "exons": exons,
        "transcript_seq": transcript_seq,
        "tx_to_genome": tx_to_genome,      # len == len(transcript_seq)
        "cds_tx_start": cds_tx_start,      # 1-based
        "cds_tx_end": cds_tx_end,          # 1-based
        "orf_seq": orf_seq,
        "gdna_seq": gdna_seq,
        "validation_length_issues": len_issues,
        "validation_mismatches": mismatches,
    }

def map_orf_mutations_to_transcript_and_genome(orf_mutations: list[str], tx_map: dict):
    """
    Map ORF mutations (relative to CDS) to:
      - transcript coordinates (mRNA orientation)
      - absolute genomic coordinates (genomic orientation; REF/ALT complemented on '-' strand)

    Returns (tx_rows, chrom_rows, gdna_rows, aa_rows, dropped):
      tx_rows:     [(input_mut, "REF<tx_pos>ALT"), ...]
      chrom_rows: [(input_mut, "REF<abs_pos>ALT"), ...]
      gdna_rows:  [(input_mut, "REF<gdna_pos>ALT"), ...]
      premrna_rows: [(input_mut, "REF<premrna_pos>ALT"), ...]
      aa_rows:    [(input_mut, "AA_REF<aa_pos>AA_ALT"), ...]
      dropped:    [(input_mut, stage, reason), ...]  -- never silently discarded.
                  stage == "nt" means no mapping row of any kind was written;
                  stage == "coordinate_space" means the token was supplied in ORF
                  space but describes a genomically discontiguous event -- almost
                  always a gd./chromosomal variant written in the wrong space. Its
                  transcript and amino-acid rows are written; the genomic ones are
                  not, because no correct single-token form exists;
                  stage == "aa" means the nt rows were written but the
                  amino-acid row was not.

    Token semantics are BFF replace-span, matching utility.splice_seq: the
    token "REF<pos>ALT" means "replace the len(REF) bases starting at pos with
    ALT". That is self-consistent for every variant class in both orientations
    and round-trips through splice_seq. It is deliberately NOT VCF anchoring --
    converting to VCF anchor-first form is vcf_converter's job, not this file's.
    The grammar cannot express an empty allele in any case, so every token
    carries its own anchor base already.
    """
    transcript_seq = tx_map["transcript_seq"]
    tx_to_genome   = tx_map["tx_to_genome"]
    cds_tx_start   = tx_map["cds_tx_start"]
    strand         = tx_map["strand"]
    orf_seq        = tx_map["orf_seq"]

    tx_rows, chrom_rows, gdna_rows, premrna_rows, aa_rows = [], [], [], [], []
    dropped: list[tuple[str, str, str]] = []

    for mut in orf_mutations:
        variant = parse_variant(mut, is_nt=True)
        if variant is None:
            dropped.append((mut, "nt", "token could not be parsed as a nucleotide variant"))
            continue
        ref = variant.ref.upper()
        alt = variant.alt.upper()

        # Transcript span (1-based inclusive, mRNA orientation)
        tx_pos = cds_tx_start + variant.pos - 1
        tx_end = tx_pos + len(ref) - 1
        if tx_pos < 1 or tx_end > len(transcript_seq):
            dropped.append((
                mut, "nt",
                f"REF span tx {tx_pos}-{tx_end} falls outside transcript (length {len(transcript_seq)})",
            ))
            continue

        observed = transcript_seq[tx_pos - 1:tx_end]
        if observed != ref:
            dropped.append((
                mut, "nt",
                f"REF '{ref}' does not match transcript '{observed}' at tx {tx_pos}-{tx_end}",
            ))
            continue

        # Transcript-oriented mapping string
        tx_rows.append((mut, f"{ref}{tx_pos}{alt}"))

        # Amino-acid mapping. Computed HERE, before any genomic work, because it
        # depends only on the token and the ORF -- infer_aavariant_from_nt takes
        # nothing else. It used to sit after the genomic block, below the
        # exon-junction `continue`, so a deletion straddling a junction lost its
        # aa row as a side effect of not being representable as one contiguous
        # GENOMIC token. Those are unrelated facts: the protein consequence of
        # such a deletion is perfectly well defined. Measured on the real run:
        # 8 transcript rows but only 7 aa rows for both SMN2 and TP53.
        #
        # get_mutant_aa was SNV-only and wrapped in a bare except, so netNglyc /
        # netphos / netMHC -- which read this 'aamutant' column -- received no row
        # at all for any non-SNV. infer_aavariant_from_nt handles snv, mnv,
        # inframe_ins, inframe_del(ins) and frameshift.
        try:
            aa_pos, wt_aa, mut_aa, consequence = infer_aavariant_from_nt(mut, orf_seq)
            if aa_pos is None:
                dropped.append((mut, "aa",
                                "amino-acid consequence could not be determined; "
                                "nt mappings still written"))
            else:
                aa_rows.append((mut, format_aa_token(aa_pos, wt_aa, mut_aa, consequence)))
        except Exception as exc:
            dropped.append((mut, "aa",
                            f"amino-acid inference failed ({exc}); nt mappings still written"))

        # Every base of the REF span, not just the first. On the minus strand
        # the transcript-first base is the genomic-LAST base, so a single
        # tx_to_genome lookup anchored the token at the wrong end for any
        # multi-base REF.
        gpos_span = tx_to_genome[tx_pos - 1:tx_end]

        # tx_to_genome jumps at every EXON-EXON junction (it is built per exon
        # block), so a REF span crossing one is contiguous in transcript space and
        # disjoint in genomic space.
        #
        # This is a COORDINATE-SPACE ERROR ON INPUT, not a representation gap. An
        # ORF-space token spanning a junction asserts that N bases were removed
        # from the MATURE transcript while the intron between them stayed intact
        # -- two independent edits at the two exon ends, not one DNA variant. A
        # real genomic deletion across that junction removes the intron as well
        # and IS contiguous in genomic space, so it maps without complaint when
        # supplied as a gd. token. Measured on SMN2 CAGAG79C: 3 bases in exon 1,
        # 2 in exon 2, ZERO intronic bases, leaping 13,648 nt of intron 1.
        #
        # So the drop names the space the token should have been written in, and
        # computes the contiguous genomic span the author most likely meant.
        # Transcript and amino-acid rows are still written: they are true
        # statements about the token as given, and are what an RNA-level edit
        # would legitimately produce.
        step = 1 if strand == "+" else -1
        contiguous = all(
            gpos_span[i + 1] - gpos_span[i] == step for i in range(len(gpos_span) - 1)
        )
        if not contiguous:
            g_lo, g_hi = min(gpos_span), max(gpos_span)
            n_skipped = (g_hi - g_lo + 1) - len(gpos_span)
            gd_lo = g_lo - tx_map["tx_start"] + 1
            dropped.append((
                mut, "coordinate_space",
                f"ORF-space token spans an exon-exon junction: its {len(ref)} REF bases "
                f"are contiguous in the transcript but occupy disjoint genomic ranges "
                f"({g_lo}..{g_hi}, skipping {n_skipped} nt of intron). No DNA variant "
                f"removes transcript bases while leaving the intron between them "
                f"intact. If a genomic deletion across this junction was meant, it "
                f"spans genomic {g_lo}-{g_hi} ({g_hi - g_lo + 1} nt) and is contiguous "
                f"there -- supply it as gd.<REF>{gd_lo}<ALT> or as a chromosomal token. "
                f"Transcript and amino-acid rows were still written",
            ))
            continue

        # Genomic-lowest coordinate is the anchor in genomic orientation. On the
        # plus strand that is the span's first base, on the minus strand its
        # last; min() is correct for both.
        g_start = min(gpos_span)

        # Genomic-orientation alleles. rc_base is a single-base helper --
        # rc_base('ACG') returns 'ACG' unchanged, which silently emitted an
        # uncomplemented, unreversed allele for every minus-strand indel.
        if strand == "+":
            g_ref, g_alt = ref, alt
        else:
            g_ref, g_alt = revcomp_seq(ref), revcomp_seq(alt)

        chrom_rows.append((mut, f"{g_ref}{g_start}{g_alt}"))

        # Relative coordinate within the gene's contiguous genomic slice (1-based)
        gdna_pos = g_start - tx_map["tx_start"] + 1
        gdna_rows.append((mut, f"{g_ref}{gdna_pos}{g_alt}"))

        # pre-mRNA is transcript-oriented, so the alleles are the ORF-space ones
        # (ref/alt), NOT the genomic-forward ones (g_ref/g_alt). Only the
        # coordinate changes.
        premrna_rows.append(
            (mut, f"{ref}{premrna_pos(g_start, max(gpos_span), tx_map)}{alt}"))

    return tx_rows, chrom_rows, gdna_rows, premrna_rows, aa_rows, dropped


# Intron-aware layer

GD_PREFIX = INTRONIC_PREFIX   # single source of truth lives in utility


def parse_gd_token(token: str):
    """Recognise a gDNA-space token, e.g. 'gd.T5000C'. Returns a Variant or None.

    The Variant's pos is 1-based within gdna_seq (the tx_start..tx_end slice),
    NOT within the ORF.

    The prefix is not decoration. An ORF token and a gDNA token are otherwise
    indistinguishable: SMN2's ORF spans positions 1-879 and its gDNA slice
    starts at 18, so 'T500C' is a valid coordinate in BOTH spaces and nothing in
    the bare grammar says which was meant. (APC happens not to collide -- its
    5'UTR offset is 17006, past its 5880 nt ORF -- but that is a per-gene
    accident, not a rule.)

    parse_variant is deliberately NOT taught this prefix. It returns None for
    'gd.T5000C' today, which means every one of its other callers -- codon_usage,
    RNAfold, the DTU pipelines -- FAILS CLOSED on an intronic token rather than
    silently reading 5000 as an ORF position. Teaching the shared parser would
    turn a dropped row into a wrong number.
    """
    if not token.startswith(GD_PREFIX):
        return None
    return parse_variant(token[len(GD_PREFIX):], is_nt=True)


CH_PREFIX = "ch."


def parse_ch_token(token: str):
    """Recognise an absolute-chromosomal token, e.g. 'ch.AACTTG564988GTT'.

    Returns a Variant whose pos is an ABSOLUTE chromosome coordinate. Separate
    from parse_gd_token because the two spaces differ by tx_start: the same
    variant is 'gd.T599A' and 'ch.T70050267A' on SMN2, and nothing in a bare
    token says which was meant.
    """
    if not token.startswith(CH_PREFIX):
        return None
    return parse_variant(token[len(CH_PREFIX):], is_nt=True)


def token_space(token: str) -> str:
    """Which coordinate space a user token is written in: orf | gdna | chrom."""
    if token.startswith(CH_PREFIX):
        return "chrom"
    if token.startswith(GD_PREFIX):
        return "gdna"
    return "orf"


def classify_genomic_span(g_start: int, g_end: int, tx_map: dict,
                          intron_index: list[dict]) -> list[dict]:
    """Every annotated feature a genomic span touches, in TRANSCRIPT order.

    Returns [{feature, n, g_start, g_end, length}, ...] where feature is 'exon'
    or 'intron'. The list has one entry for a variant inside a single feature and
    several for one that spans a boundary.

    Both the ordinals and the list order are transcript-relative, matching
    build_intron_index. Exons were previously numbered and sorted by ascending
    GENOMIC coordinate while introns were numbered in transcript order, so on a
    minus-strand gene the two conventions disagreed inside a single CSV row:
    measured on TP53, a span across intron 1 reported ['exon10','intron1',
    'exon11'] when the features either side of intron 1 are exons 1 and 2.
    The list order matters just as much -- callers append per-feature pieces in
    this order, and a bracketed multi-piece cell read left-to-right must
    reconstruct the edit in the direction the molecule is read.

    This is the detection step the whole intron layer turns on. A variant is not
    'exonic' or 'intronic' as a property of how it was written -- it is a
    property of WHERE ITS REF SPAN LANDS, which is only knowable after resolving
    the token to genomic coordinates. Three outcomes, all legitimate:

        wholly exonic    has an ORF/transcript address; scoreable everywhere
        wholly intronic  no ORF address; scoreable on the intron and pre-mRNA
        spans features   a splice-site variant. Genomic addresses are exact; the
                         MATURE transcript is a splicing prediction, not
                         arithmetic, so no transcript or aa row is emitted.

    The third case is why a 'wholly inside one intron' rule was wrong: a deletion
    removing the last 3 bases of an exon plus the first 6 of an intron is the
    canonical splice-donor variant, and it was being dropped as unmappable.
    """
    minus = tx_map["strand"] == "-"
    exons = sorted(tx_map["exons"])
    n_exons = len(exons)
    out = []
    for k, (s, e) in enumerate(exons, 1):
        if g_start <= e and g_end >= s:
            lo, hi = max(g_start, s), min(g_end, e)
            # transcript ordinal: on '-' the genomic-last exon is exon 1
            n = (n_exons - k + 1) if minus else k
            out.append({"feature": "exon", "n": n, "g_start": lo, "g_end": hi,
                        "length": hi - lo + 1})
    for i in intron_index:
        if g_start <= i["g_end"] and g_end >= i["g_start"]:
            lo, hi = max(g_start, i["g_start"]), min(g_end, i["g_end"])
            out.append({"feature": "intron", "n": i["n"], "g_start": lo,
                        "g_end": hi, "length": hi - lo + 1})
    return sorted(out, key=lambda f: f["g_start"], reverse=minus)


def span_class(features: list[dict]) -> str:
    """'exonic' | 'intronic' | 'spanning' for a classified span."""
    kinds = {f["feature"] for f in features}
    if not features:
        return "outside"
    if kinds == {"exon"}:
        return "exonic" if len(features) == 1 else "spanning"
    if kinds == {"intron"}:
        return "intronic" if len(features) == 1 else "spanning"
    return "spanning"


def build_intron_index(tx_map: dict) -> list[dict]:
    """WT introns as INDIVIDUAL molecules, numbered in transcript order.

    The difference between the genomic and transcript records is the SUM of the
    introns, not one sequence -- for SMN2 that is 26,445 nt across 8 introns
    ranging from 159 nt to 13,648 nt, with intron 1 alone accounting for 51.6%.
    Any per-sequence statistic computed on the concatenation would be an
    intron-1-weighted average of molecules that never coexist.

    Each entry carries:
        n            1-based ordinal in TRANSCRIPT order. On the minus strand
                     intron 1 is at the HIGHEST genomic coordinate.
        g_start/end  1-based inclusive absolute genomic bounds (always ascending)
        slice_start/end  1-based inclusive bounds within gdna_seq
        seq          TRANSCRIPT orientation -- reverse-complemented on '-'.

    The orientation is not cosmetic. Neither RNAfold nor miranda is strand-aware;
    both consume whatever orientation they are handed. Measured on this machine:
    miranda finds a perfect let-7a site on the forward sequence (score 200.00,
    dG -39.30) and ZERO hits on its reverse complement, and RNAfold returns
    -27.60 vs -27.30 kcal/mol for a sequence and its reverse complement (the
    gap is the G-U wobble asymmetry, since G-U reverse-complements to A-C, which
    does not pair). Handing either tool a genomic-forward minus-strand intron
    yields a plausible-looking wrong structure, or zero miRNA hits that read as
    'no interaction' rather than 'wrong strand'.
    """
    exons = sorted(tx_map["exons"], key=lambda x: x[0])   # genomic ascending
    strand = tx_map["strand"]
    tx_start = tx_map["tx_start"]
    gdna = tx_map["gdna_seq"]

    gaps: list[tuple[int, int]] = []
    for i in range(len(exons) - 1):
        g_start = exons[i][1] + 1
        g_end = exons[i + 1][0] - 1
        if g_end < g_start:
            # Abutting or overlapping exon records leave no intron between them.
            continue
        gaps.append((g_start, g_end))

    if strand == "-":
        gaps = gaps[::-1]        # transcript order runs high -> low genomic

    index: list[dict] = []
    for n, (g_start, g_end) in enumerate(gaps, 1):
        # gdna_seq[0] is genomic position tx_start (fetch is tx_start-1 .. tx_end)
        s0 = g_start - tx_start
        e0 = g_end - tx_start + 1
        seq = gdna[s0:e0]
        if strand == "-":
            seq = revcomp_seq(seq)
        index.append({
            "n": n,
            "g_start": g_start,
            "g_end": g_end,
            "slice_start": s0 + 1,
            "slice_end": e0,
            "length": g_end - g_start + 1,
            "seq": seq,
        })
    return index


def build_pre_mrna(tx_map: dict) -> tuple[str, list[tuple[int, int]]]:
    """The pre-mRNA: the full tx_start..tx_end span in TRANSCRIPT orientation.

    Returns (seq, exon_spans) where exon_spans are 0-based [start, end) offsets
    into seq, in transcript order.

    This is NOT a new fetch. It is the existing `genomic` record with the strand
    applied -- identical on '+', reverse-complemented on '-'. The two records are
    kept side by side rather than one replacing the other because they answer
    different questions and different tools want different ones: genesplicer
    scores genomic-forward sequence with the strand supplied separately, while
    AlphaFold3, RNAfold and miranda are all strand-blind and must be handed the
    orientation the molecule actually has (measured: miranda finds a perfect
    let-7a site on the forward sequence and ZERO on its reverse complement).

    The pre-mRNA is the correct substrate for RNA-binding-protein work. Intronic
    splicing regulatory elements -- branch point, polypyrimidine tract, ISE/ISS --
    exist only here, and the 5' splice site spans the exon|intron junction, so a
    bare intron record truncates the exonic half of every U1 snRNP site. The
    pre-mRNA carries both halves.
    """
    strand = tx_map["strand"]
    tx_start = tx_map["tx_start"]
    tx_end = tx_map["tx_end"]
    gdna = tx_map["gdna_seq"]

    exons = sorted(tx_map["exons"], key=lambda x: x[0])
    if strand == "+":
        seq = gdna
        spans = [(s - tx_start, e - tx_start + 1) for s, e in exons]
    else:
        seq = revcomp_seq(gdna)
        # A genomic position g sits at pre-mRNA offset (tx_end - g), so an exon
        # [s, e] runs from (tx_end - e) to (tx_end - s) inclusive.
        spans = [(tx_end - e, tx_end - s + 1) for s, e in exons]
        spans.sort()          # ascending pre-mRNA offset == transcript order
    return seq, spans


def premrna_pos(g_start: int, g_end: int, tx_map: dict) -> int:
    """1-based pre-mRNA position of a genomic span's first base in TRANSCRIPT order.

    On '+' that is the genomic-lowest base and pre_mRNA[0] is tx_start; on '-' it
    is the genomic-highest base and pre_mRNA[0] is tx_end.

    This lives here, not in the consumer pipelines, because it is the only place
    that holds the strand. RNAfold reads a FASTA and a mapping CSV and has no
    strand field; making it re-derive this would duplicate the one piece of
    arithmetic whose sign flips per gene.
    """
    if tx_map["strand"] == "+":
        return g_start - tx_map["tx_start"] + 1
    return tx_map["tx_end"] - g_end + 1


def check_pre_mrna_splices_to_transcript(tx_map: dict) -> tuple[bool, str]:
    """Splice the pre-mRNA at its exon spans and compare against transcript_seq.

    Returns (ok, detail). This is a genuine end-to-end assertion on the whole
    coordinate chain -- exon bounds from the annotation, strand handling, and the
    tx_start offset all have to be right simultaneously for it to pass, and any
    one of them being wrong produces a mismatch rather than a plausible-looking
    sequence.
    """
    seq, spans = build_pre_mrna(tx_map)
    spliced = "".join(seq[s:e] for s, e in spans)
    expected = tx_map["transcript_seq"]
    if spliced == expected:
        return True, f"pre-mRNA {len(seq)} nt splices to transcript {len(expected)} nt"
    return False, (
        f"pre-mRNA splice check FAILED: spliced {len(spliced)} nt != transcript "
        f"{len(expected)} nt"
        + ("" if len(spliced) != len(expected) else "; same length, different bases")
    )


PIECE_RE = re.compile(r"^(?:(intron\d+):)?([ACGTU]+)(\d+)(?:([ACGTU]+)|del)$", re.I)


def parse_piece_token(token: str):
    """Read one decomposed piece back: 'CAG79A', 'GTGAGGTCG1del', 'intron1:T501A'.

    Returns {record, ref, pos, alt, deleted} or None. `alt` is '' when the piece
    is wholly deleted, and `deleted` says so explicitly rather than making the
    caller test for an empty string.

    This is deliberately NOT parse_variant. A piece is a FRAGMENT of one edit,
    not a variant: a spanning deletion retains a single anchor base that belongs
    to exactly one piece, so the others have no ALT at all. Variant refuses an
    empty allele at construction (utility.py:342) precisely so that a pure
    deletion cannot be represented without the VCF anchor convention -- which is
    right for whole variants and impossible for pieces.

    Keeping them separate preserves the fail-closed property: parse_variant
    returns None for 'GTGAGGTCG1del', so a pipeline that does not understand
    decomposition drops the piece instead of scoring a fragment as if it were
    the whole edit.
    """
    if not isinstance(token, str):
        return None
    m = PIECE_RE.match(token.strip())
    if not m:
        return None
    record, ref, pos, alt = m.group(1), m.group(2).upper(), int(m.group(3)), m.group(4)
    if pos < 1:
        return None
    return {"record": record, "ref": ref, "pos": pos,
            "alt": (alt or "").upper(), "deleted": alt is None}


def piece_fields(token: str):
    """(pos1, ref, alt, kind, length_delta) for a piece token, or None.

    The Variant-shaped view of a piece, so a consumer can treat a decomposed
    piece and a whole variant through one code path. `kind` uses exactly the
    vocabulary of Variant.kind (utility.py) -- snv, mnv, insertion, deletion,
    delins -- with one extension: Variant cannot hold an empty allele, so its
    rule for 'deletion' is len(alt) == 1, whereas a wholly deleted piece has
    len(alt) == 0. Both are deletions.
    """
    d = parse_piece_token(token)
    if d is None:
        return None
    ref, alt = d["ref"], d["alt"]
    if len(ref) == 1 and len(alt) == 1:
        kind = "snv"
    elif len(ref) == len(alt):
        kind = "mnv"
    elif len(ref) == 1:
        kind = "insertion"
    elif len(alt) <= 1:
        kind = "deletion"
    else:
        kind = "delins"
    return d["pos"], ref, alt, kind, len(alt) - len(ref)


def split_piece_cell(cell: str) -> list[str]:
    """Split one mapping cell into its pieces: 'a' -> ['a'], '[a,b]' -> ['a','b'].

    A variant spanning several features occupies several addresses in the orf,
    intron and pre-mRNA columns, written as a bracketed list. Reading a cell is
    therefore two steps -- split, then parse_piece_token each piece -- and the
    split half was reimplemented per consumer (miranda inline, RNAfold in
    _load_intron_premrna) while only parse_piece_token was shared. Two readers
    of one format drift; this is the half that was missing.

    Order is TRANSCRIPT order, matching classify_genomic_span, so the first
    piece is the one carrying the edit's retained anchor base.
    """
    if not isinstance(cell, str):
        return []
    cell = cell.strip()
    if not cell:
        return []
    if cell.startswith("[") and cell.endswith("]"):
        return [p.strip() for p in cell[1:-1].split(",") if p.strip()]
    return [cell]


def genomic_to_orf(g_pos: int, tx_map: dict):
    """ORF position (1-based) for an absolute genomic coordinate, or None.

    None means the base is not in the CDS -- intronic, UTR, or outside the
    transcript. tx_to_genome is the authority; it already encodes the strand.
    """
    try:
        tx_idx = tx_map["tx_to_genome"].index(g_pos)
    except ValueError:
        return None
    orf_pos = (tx_idx + 1) - tx_map["cds_tx_start"] + 1
    return orf_pos if 1 <= orf_pos <= len(tx_map["orf_seq"]) else None


def _piece_tokens(feat, ref_full, alt_full, g_lo, g_hi, tx_map, intron_index):
    """Per-feature tokens for one slice of a variant's REF span.

    Returns (orf_token, intron_token, premrna_token). Any of the three is None
    where that space has no address for this piece.

    The ALT is carried on the FIRST piece only. A multi-feature variant is ONE
    edit: a 12 nt deletion retains exactly one anchor base, so only one piece can
    hold it. Repeating it on every piece would say the anchor survived twice.

    A piece with no ALT is written 'REF<pos>del'. The bare grammar cannot express
    it -- REF<pos>ALT requires both sides non-empty, so 'GTGAGGTCG1' returns None
    from parse_variant -- and the two alternatives are both worse: repeating the
    anchor states something untrue, and leaving the cell blank makes the orf
    column understate how much exonic sequence a splice variant removes.
    """
    def _tok(r, pos, a):
        return f"{r}{pos}{a}" if a else f"{r}{pos}del"
    strand = tx_map["strand"]
    seg_ref = ref_full[feat["g_start"] - g_lo: feat["g_end"] - g_lo + 1]
    # TRANSCRIPT-first, not genomic-first. Every piece token below is written in
    # transcript orientation, so the retained anchor base belongs to the piece
    # the molecule is read into first -- the genomic-LOWEST on '+', the
    # genomic-HIGHEST on '-'. Comparing against g_lo on both strands put the
    # anchor on the wrong end of every minus-strand multi-feature edit: measured
    # on TP53 intron 5, the ALT landed on ORF 560 (transcript-last) while ORF 557
    # (transcript-first) was written 'del'.
    is_first = (feat["g_end"] == g_hi) if strand == "-" else (feat["g_start"] == g_lo)
    seg_alt = alt_full if is_first else ""

    # transcript orientation for everything pre-mRNA-facing
    t_ref = seg_ref if strand == "+" else revcomp_seq(seg_ref)
    t_alt = seg_alt if strand == "+" else revcomp_seq(seg_alt)
    p_pos = premrna_pos(feat["g_start"], feat["g_end"], tx_map)
    premrna_token = _tok(t_ref, p_pos, t_alt)

    orf_token = None
    intron_token = None
    if feat["feature"] == "exon":
        anchor = feat["g_start"] if strand == "+" else feat["g_end"]
        o_pos = genomic_to_orf(anchor, tx_map)
        if o_pos is not None:
            orf_token = _tok(t_ref, o_pos, t_alt)
    else:
        host = next(i for i in intron_index if i["n"] == feat["n"])
        if strand == "+":
            i_pos = feat["g_start"] - host["g_start"] + 1
        else:
            i_pos = host["g_end"] - feat["g_end"] + 1
        body = _tok(t_ref, i_pos, t_alt)
        intron_token = f"intron{feat['n']}:{body}"
    return orf_token, intron_token, premrna_token


def map_genomic_variants(tokens: list[str], tx_map: dict, intron_index: list[dict]):
    """Map gd./ch. tokens by WHERE THEIR REF SPAN LANDS, not by their prefix.

    Returns (chrom_rows, gdna_rows, ipm_rows, dropped, used_introns) where every
    row is a dict, because a variant spanning several features has several
    addresses in some spaces and one in others.

        chrom_rows / gdna_rows  {mutant, orf, value}
        ipm_rows                {mutant, orf, intron, premrna}

    `mutant` is always the verbatim user token, so the {GENE}-{token} pkey every
    downstream pipeline builds still resolves. `orf` is the derived ORF-space
    token, which exists only for exonic pieces.

    Routing follows classify_genomic_span:

        exonic     full mapping INCLUDING the ORF/transcript address. Previously
                   dropped as "not wholly inside one intron" -- a genomic SNV in
                   an exon has an exact ORF position and was being discarded.
        intronic   intron + pre-mRNA rows; no ORF address exists.
        spanning   genomic addresses are exact and are emitted; the per-feature
                   pieces are carried as lists. No transcript or aa row: the
                   variant removes a splice site, so the mature transcript is a
                   PREDICTION (genesplicer / spliceai), not arithmetic. Emitting
                   a naively-spliced consequence would be a fabricated value.
    """
    gdna = tx_map["gdna_seq"]
    tx_start = tx_map["tx_start"]

    chrom_rows, gdna_rows, ipm_rows = [], [], []
    dropped: list[tuple[str, str, str]] = []
    used: dict[int, dict] = {}

    for tok in tokens:
        space = token_space(tok)
        variant = parse_ch_token(tok) if space == "chrom" else parse_gd_token(tok)
        if variant is None:
            dropped.append((tok, "nt", "token could not be parsed after its coordinate prefix"))
            continue
        ref, alt = variant.ref.upper(), variant.alt.upper()

        # normalise to absolute genomic, whichever space it arrived in
        g_lo = variant.pos if space == "chrom" else tx_start + variant.pos - 1
        g_hi = g_lo + len(ref) - 1
        s0 = g_lo - tx_start

        if s0 < 0 or s0 + len(ref) > len(gdna):
            dropped.append((tok, "nt",
                            f"REF span genomic {g_lo}-{g_hi} falls outside the "
                            f"transcript bounds {tx_start}-{tx_map['tx_end']}"))
            continue
        observed = gdna[s0:s0 + len(ref)]
        if observed != ref:
            dropped.append((tok, "nt",
                            f"REF '{ref}' does not match the genome '{observed}' "
                            f"at genomic {g_lo}-{g_hi}"))
            continue

        feats = classify_genomic_span(g_lo, g_hi, tx_map, intron_index)
        kind = span_class(feats)
        if kind == "outside":
            dropped.append((tok, "nt",
                            f"REF span genomic {g_lo}-{g_hi} touches no annotated "
                            f"exon or intron of this transcript"))
            continue

        orf_toks, intron_toks, pre_toks = [], [], []
        for f in feats:
            o, i, pre = _piece_tokens(f, ref, alt, g_lo, g_hi, tx_map, intron_index)
            if o: orf_toks.append(o)
            if i:
                intron_toks.append(i)
                used[f["n"]] = next(x for x in intron_index if x["n"] == f["n"])
            pre_toks.append(pre)

        def cell(vals):
            """One value stays scalar; several become a bracketed list."""
            if not vals:
                return ""
            return vals[0] if len(vals) == 1 else "[" + ",".join(vals) + "]"

        orf_cell = cell(orf_toks)
        gdna_pos = g_lo - tx_start + 1
        # chromosome / gDNA are genomic-forward and always contiguous here --
        # the span is one block in genomic space whatever features it crosses.
        chrom_rows.append({"mutant": tok, "orf": orf_cell,
                           "value": f"{ref}{g_lo}{alt}"})
        gdna_rows.append({"mutant": tok, "orf": orf_cell,
                          "value": f"{ref}{gdna_pos}{alt}"})

        # The combined file carries one row per SPACE the variant occupies:
        # exonic pieces on the orf row, intronic pieces on the intron row.
        #
        # Gated on which FEATURES the span touches, not on which produced a
        # token. An exonic piece has no ORF address whenever it lands in a UTR
        # or a non-coding exon -- genomic_to_orf correctly returns None -- and
        # gating on orf_toks then discarded the pre-mRNA token that had already
        # been computed for it, with no drop reason recorded. Measured: every
        # 5'UTR and 3'UTR token on SMN2 and TP53 produced chromosome and gDNA
        # rows and NO pre-mRNA row, and a canonical donor-site deletion at TP53
        # intron 1 lost its exonic half entirely, because TP53 exon 1 is
        # non-coding. pre-mRNA is the one space in which every piece of every
        # class has a coordinate, so its column is populated whenever the
        # feature is touched.
        exon_pre = [p for f, p in zip(feats, pre_toks) if f["feature"] == "exon"]
        intron_pre = [p for f, p in zip(feats, pre_toks) if f["feature"] == "intron"]
        if exon_pre:
            ipm_rows.append({"mutant": tok, "orf": orf_cell, "intron": "",
                             "premrna": cell(exon_pre)})
        if intron_pre:
            ipm_rows.append({"mutant": tok, "orf": "", "intron": cell(intron_toks),
                             "premrna": cell(intron_pre)})

        if kind == "spanning":
            dropped.append((tok, "splice_site",
                            f"REF span covers "
                            + "+".join(f"{f['feature']}{f['n']}:{f['length']}nt" for f in feats)
                            + ". Genomic, pre-mRNA and intron rows written; NO transcript or "
                              "amino-acid row -- the variant removes a splice site, so the "
                              "mature transcript is a splicing prediction, not arithmetic"))

    used_introns = [used[n] for n in sorted(used)]
    return chrom_rows, gdna_rows, ipm_rows, dropped, used_introns


# CLI

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Generate exon-aware mapping CSVs and FASTAs (ORF/transcript/genomic) per gene."
    )
    p.add_argument("-m", "--mutations", required=True,
                   help="Path to a mutations CSV or a directory of CSVs (expects '<GENE>_mutations.csv').")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA (indexed).")
    p.add_argument("-a", "--annotation", required=True, help="Gene annotation file consumed by utility.get_genome_loc.")
    p.add_argument("-o", "--output", help="Output directory; defaults to current working directory.")
    p.add_argument("-Ec", "--exclude-chromosome", action="store_true", help="Exclude the chromosome mapping CSV.")
    p.add_argument("-Eg", "--exclude-genomic", action="store_true", help="Exclude the gDNA mapping CSVs (relative to gene genomic slice).")
    p.add_argument("-Et", "--exclude-transcript", action="store_true", help="Exclude the transcript mapping CSVs.")
    p.add_argument("-EA", "--exclude-aa", action="store_true", help="Exclude the amino-acid mapping CSVs.")
    p.add_argument("-Ei", "--exclude-intron", action="store_true",
                   help="Exclude the intron mapping CSVs and the per-intron FASTA records.")
    p.add_argument("-Ep", "--exclude-premrna", action="store_true",
                   help="Exclude the pre-mRNA mapping CSVs.")
    p.add_argument("--orf", help="Optional ORF FASTA (file or directory). If omitted, ORF is inferred from transcript.")
    p.add_argument("--force-cds",
        help="Force specific transcript: single accession (e.g., NM_022162.3) for all genes, "
             "or CSV file mapping genes to transcript IDs (gene,transcript_id columns).")
    p.add_argument("--verbose", action="store_true", help="Print detailed ORF/mutation validation messages.")
    args = p.parse_args(argv)

    mut_path = Path(args.mutations)
    output = Path(args.output) if args.output else Path.cwd()

    do_chrom = not args.exclude_chromosome
    do_genomic = not args.exclude_genomic
    do_transcript = not args.exclude_transcript
    do_aa = not args.exclude_aa
    do_intron = not args.exclude_intron
    do_premrna = not args.exclude_premrna

    capture_verbose = args.verbose
    verbose_log: list[str] | None = [] if capture_verbose else None
    gene_metrics = {
        "total_genes": 0,
        "passed_genes": 0,
        "failed_genes": 0,
        "total_mutations": 0,
        "passed_mutations": 0,
        "failed_mutations": 0,
        "failed_length": 0,
        "failed_mismatch": 0,
        "dropped_mapping": 0,
        "dropped_aa_only": 0,
        "dropped_coordinate_space": 0,
        "splice_site_variants": 0,
    }
    orf_lookup = {}
    if args.orf:
        try:
            orf_lookup = load_supplied_orfs(args.orf)
            print(f"Loaded {len(orf_lookup)} ORF sequence(s) from {args.orf}")
        except Exception as exc:
            print(f"Error loading supplied ORFs: {exc}", file=sys.stderr)
            sys.exit(1)

    # Load transcript ID overrides from --force-cds
    global_transcript, transcript_lookup = None, {}
    if args.force_cds:
        try:
            global_transcript, transcript_lookup = parse_force_cds(args.force_cds)
            if global_transcript:
                print(f"Forcing transcript {global_transcript} for all genes")
            else:
                print(f"Loaded {len(transcript_lookup)} transcript override(s) from {args.force_cds}")
        except Exception as exc:
            print(f"Error loading --force-cds: {exc}", file=sys.stderr)
            sys.exit(1)

    if mut_path.is_file():
        files = [mut_path]
    else:
        files = sorted(mut_path.glob("*.csv"))

    ok, fail = 0, []

    for f in files:
        gene: str | None = None
        muts: list[str] = []
        try:
            gene = extract_gene_from_filename(str(f))
            print(f"\n== {gene} ==")

            # Load mutations, then separate coordinate spaces BEFORE any ORF work.
            muts = trim_muts(str(f))
            if not muts:
                print("  (no mutations) -> skipping")
                continue
            # This split MUST precede resolve_orf_from_sources. A gDNA position
            # is typically far larger than the ORF, and derive_orf_from_transcript
            # scores ORF candidates on whether they cover every mutation position
            # -- one gDNA token would force selection of the longest ORF in the
            # transcript and corrupt the coordinates of every real ORF mutation in
            # the gene. It would also fail validate_mutations_against_orf, which
            # marks the WHOLE GENE failed.
            orf_muts, gd_muts = split_intronic_tokens(muts)
            if gd_muts:
                print(f"  Coordinate spaces: {len(orf_muts)} ORF, {len(gd_muts)} gDNA/intronic")

            # Mint the primary key ONCE, here, for every token of this gene.
            # This is the only point in the system holding both halves before
            # anything downstream runs -- gene comes from the filename, token
            # from the row -- so every pipeline can read a key it never has to
            # derive. Minting per-pipeline is what produced four independent
            # '{GENE}-{token}' mints that all had to agree.
            pkey_of = {tok: mint_pkey(gene, tok) for tok in muts}

            gene_key = gene.upper()
            gene_metrics["total_genes"] += 1
            supplied_orf_seq = None
            if orf_lookup:
                supplied_orf_seq = orf_lookup.get(gene_key)
                if supplied_orf_seq is None:
                    warnings.warn(
                        f"No supplied ORF found for {gene_key} in {args.orf}; falling back to transcript-derived ORF.",
                        RuntimeWarning,
                    )

            # Resolve transcript override: global takes precedence, then per-gene lookup
            transcript_override = global_transcript or transcript_lookup.get(gene_key)

            # Build sequences + maps. Only ORF-space tokens participate -- see
            # split_intronic_tokens for why this separation cannot happen later.
            tx_map = build_transcript_seq_and_map(
                gene,
                args.annotation,
                args.reference,
                orf_muts,
                supplied_orf=supplied_orf_seq,
                transcript_id=transcript_override,
                verbose=args.verbose,
                log_messages=verbose_log,
            )

            # Fail the gene BEFORE writing anything. ORF validation is already
            # decided by this point -- build_transcript_seq_and_map stores its
            # result on tx_map -- but the check used to sit after every write, so
            # a gene that failed still emitted a full set of header-only CSVs
            # (measured: 4 files, 1 line each). A consumer reading the mappings
            # directory later, without the process exit status, cannot tell those
            # from a gene that genuinely had no mutations. Writing nothing is
            # unambiguous.
            len_issues = tx_map.get("validation_length_issues") or []
            mismatches = tx_map.get("validation_mismatches") or []
            if len_issues or mismatches:
                failed_tokens = {m for m, _ in len_issues} | {m for m, _ in mismatches}
                err_msg = (
                    f"{gene}: ORF/mutation validation failed for "
                    f"{len(failed_tokens)} mutation(s); no files written."
                )
                print(f"  ERROR: {err_msg}", file=sys.stderr)
                for _mut, msg in list(len_issues) + list(mismatches):
                    print(f"    - {msg}", file=sys.stderr)
                fail.append((f.name, err_msg))
                gene_metrics["failed_genes"] += 1
                gene_metrics["total_mutations"] += len(muts)
                gene_metrics["failed_mutations"] += len(failed_tokens)
                gene_metrics["passed_mutations"] += max(len(muts) - len(failed_tokens), 0)
                gene_metrics["failed_length"] += len(len_issues)
                gene_metrics["failed_mismatch"] += len(mismatches)
                continue

            # Build mappings
            tx_rows, chrom_rows, gdna_rows, premrna_rows, aa_rows, dropped_rows = \
                map_orf_mutations_to_transcript_and_genome(orf_muts, tx_map)

            # ORF-space rows carry the user token as both mutant and orf: the
            # token IS its own ORF address, so no decomposition applies.
            chrom_dicts = [{"mutant": m, "orf": m, "value": v} for m, v in chrom_rows]
            gdna_dicts = [{"mutant": m, "orf": m, "value": v} for m, v in gdna_rows]
            ipm_dicts = [{"mutant": m, "orf": m, "intron": "", "premrna": v}
                         for m, v in premrna_rows]

            used_introns: list[dict] = []
            if gd_muts:
                intron_index = build_intron_index(tx_map)
                g_chrom, g_gdna, g_ipm, g_dropped, used_introns = \
                    map_genomic_variants(gd_muts, tx_map, intron_index)
                chrom_dicts += g_chrom
                gdna_dicts += g_gdna
                ipm_dicts += g_ipm
                dropped_rows = dropped_rows + g_dropped
                print(f"  Introns: {len(intron_index)} in transcript,"
                      f" {len(used_introns)} carrying a mutation")

            # Write FASTA: ORF / transcript / genomic, plus one record per
            # mutation-carrying intron.
            #
            # The intron header has NO SPACE. read_fasta keys on the first
            # whitespace-delimited word, so '>intron 5|...' collapses every
            # intron onto the key 'intron' and silently keeps only the last --
            # measured: 6 records written, 4 parsed back.
            fasta_path = output / "fastas" / f"{gene}.fasta"
            fasta_path.parent.mkdir(parents=True, exist_ok=True)

            # pre_mRNA is the genomic span with the strand applied. Both records
            # are emitted: 'genomic' stays genomic-forward for tools that take a
            # strand flag, 'pre_mRNA' is the molecule's own orientation for the
            # strand-blind ones.
            pre_mrna_seq, _ = build_pre_mrna(tx_map)
            pre_ok, pre_detail = check_pre_mrna_splices_to_transcript(tx_map)
            if not pre_ok:
                print(f"  WARNING: {gene}: {pre_detail}", file=sys.stderr)
                if verbose_log is not None:
                    verbose_log.append(f"{gene}: {pre_detail}")

            records = {
                "ORF": tx_map["orf_seq"],
                "transcript": tx_map["transcript_seq"],
                "genomic": tx_map["gdna_seq"],
                "pre_mRNA": pre_mrna_seq,
            }
            if do_intron:
                for intron in used_introns:
                    key = f"intron{intron['n']}|{intron['g_start']}-{intron['g_end']}"
                    records[key] = intron["seq"]
            write_fasta(fasta_path, records)
            print(f"  FASTA: {fasta_path}  ({len(records)} records)")

            nt_dropped = {mut for mut, stage, _ in dropped_rows if stage == "nt"}
            aa_dropped = {mut for mut, stage, _ in dropped_rows if stage == "aa"}
            # Genomic-only: transcript and amino-acid rows exist, but the variant
            # has no single contiguous genomic token. NOT a failure -- it is a
            # partial mapping, and counting it as one would understate a gene
            # whose rows are mostly present.
            coord_dropped = {mut for mut, stage, _ in dropped_rows if stage == "coordinate_space"}
            # Splice-site variants are NOT failures: their genomic, pre-mRNA and
            # intron rows are all written. Only the transcript/aa rows are
            # withheld, because the mature transcript is a splicing prediction.
            splice_noted = {mut for mut, stage, _ in dropped_rows if stage == "splice_site"}
            if dropped_rows:
                print(
                    f"  Dropped during mapping: {len(nt_dropped)} nt, "
                    f"{len(coord_dropped)} wrong-coordinate-space, {len(aa_dropped)} aa-only"
                    f"; {len(splice_noted)} splice-site (rows written)",
                    file=sys.stderr,
                )
                for mut, stage, reason in dropped_rows:
                    print(f"    - [{stage}] {mut}: {reason}", file=sys.stderr)
                    if verbose_log is not None:
                        verbose_log.append(
                            f"{gene}: mutation {mut} dropped during mapping [{stage}] ({reason})."
                        )

            # mutant = the user's verbatim token (pkey stays resolvable);
            # orf = the derived ORF-space token, present only where an ORF
            # address exists. A variant spanning a splice site carries several,
            # bracketed.
            if do_chrom:
                chrom_csv = output / "mappings" / "chromosome" / f"chr_mapping_{gene}.csv"
                write_dict_csv(chrom_csv, ["pkey", "mutant", "orf", "chromosome"],
                               [{**r, "pkey": pkey_of.get(r["mutant"], ""),
                                 "chromosome": r["value"]} for r in chrom_dicts])
                print(f"  Chromosome mapping: {chrom_csv}  ({len(chrom_dicts)} rows)")

            if do_genomic:
                gdna_csv = output / "mappings" / "gDNA" / f"genomic_mapping_{gene}.csv"
                write_dict_csv(gdna_csv, ["pkey", "mutant", "orf", "genomic"],
                               [{**r, "pkey": pkey_of.get(r["mutant"], ""),
                                 "genomic": r["value"]} for r in gdna_dicts])
                print(f"  Genomic mapping: {gdna_csv}  ({len(gdna_dicts)} rows)")

            if do_transcript:
                tx_csv = output / "mappings" / "transcript" / f"transcript_mapping_{gene}.csv"
                write_dict_csv(tx_csv, ["pkey", "mutant", "transcript"],
                               [{"pkey": pkey_of.get(m, ""), "mutant": m,
                                 "transcript": v} for m, v in tx_rows])
                print(f"  Transcript mapping: {tx_csv}  ({len(tx_rows)} rows)")

            if do_aa:
                aa_csv = output / "mappings" / "aa" / f"{gene}_aa_mapping.csv"
                write_dict_csv(aa_csv, ["pkey", "mutant", "aamutant"],
                               [{"pkey": pkey_of.get(m, ""), "mutant": m,
                                 "aamutant": v} for m, v in aa_rows])
                print(f"  Amino-acid mapping: {aa_csv}  ({len(aa_rows)} rows)")

            # One combined file replaces the separate intron and pre-mRNA CSVs.
            # `orf` and `intron` are mutually exclusive per row, so a variant with
            # both exonic and intronic pieces contributes TWO rows under one
            # `mutant`. pre-mRNA is the only nucleotide space in which every
            # piece of every class has a coordinate, so that column is always
            # populated.
            if do_premrna and ipm_dicts:
                ipm_csv = (output / "mappings" / "intron_premRNA"
                           / f"intron_premRNA_mapping_{gene}.csv")
                write_dict_csv(ipm_csv,
                               ["pkey", "mutant", "orf", "intron", "pre-mRNA-Transcript"],
                               [{**r, "pkey": pkey_of.get(r["mutant"], ""),
                                 "pre-mRNA-Transcript": r["premrna"]} for r in ipm_dicts])
                print(f"  Intron/pre-mRNA mapping: {ipm_csv}  ({len(ipm_dicts)} rows)")

            # The mutant -> pkey inversion table. The four pipelines that get
            # only a sequence name back from their tool (netNglyc, netphos,
            # netMHC via FASTA headers; miranda via output filenames) need to go
            # from that name to the token, and a hashed key cannot be split back
            # apart. Kept in its own subdirectory because discover_mapping_files
            # rglobs a tree and keys purely on gene name, so a sixth mapping file
            # sharing a directory with the coordinate mappings would overwrite
            # one of them at random.
            pkey_csv = output / "mappings" / "pkey" / f"pkey_mapping_{gene}.csv"
            write_dict_csv(pkey_csv, ["pkey", "mutant"],
                           [{"pkey": pkey_of[t], "mutant": t} for t in muts])
            print(f"  Pkey mapping: {pkey_csv}  ({len(muts)} rows)")

            mut_dest = output / "mappings" / "mutations" / f"{gene}_mutations.csv"
            mut_dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, mut_dest)
            print(f"  Mutations copy: {mut_dest}")
            total_mut = len(muts)
            # A mutation dropped in mapping produced no nt row at all, so it is
            # a failure regardless of whether ORF validation flagged it. Counting
            # only ORF-validation failures reported such a mutation as PASSED
            # while its row was absent from every CSV.
            #
            # Union, not sum. One token can be BOTH an ORF-validation failure and
            # a mapping drop, and adding the list lengths counted it twice --
            # measured on a 1-token file: "Mutations: 1 (passed: 0, failed: 2)".
            # len_issues/mismatches carry their token for exactly this reason;
            # genes reaching this point have neither, but the union keeps the
            # arithmetic correct if that ever changes.
            failed_tokens = ({m for m, _ in len_issues}
                             | {m for m, _ in mismatches}
                             | set(nt_dropped))
            failed_mut = len(failed_tokens)
            passed_mut = max(total_mut - failed_mut, 0)
            gene_metrics["total_mutations"] += total_mut
            gene_metrics["passed_mutations"] += passed_mut
            gene_metrics["failed_mutations"] += failed_mut
            gene_metrics["failed_length"] += len(len_issues)
            gene_metrics["failed_mismatch"] += len(mismatches)
            gene_metrics["dropped_mapping"] += len(nt_dropped)
            gene_metrics["dropped_aa_only"] += len(aa_dropped)
            gene_metrics["dropped_coordinate_space"] += len(coord_dropped)
            gene_metrics["splice_site_variants"] += len(splice_noted)

            ok += 1
            gene_metrics["passed_genes"] += 1
        except Exception as e:
            fail_msg = str(e)
            fail.append((f.name, fail_msg))
            print(f"  ERROR: {fail_msg}", file=sys.stderr)
            gene_metrics["failed_genes"] += 1
            if 'muts' in locals() and isinstance(muts, list):
                gene_metrics["total_mutations"] += len(muts)
                gene_metrics["failed_mutations"] += len(muts)
                if verbose_log is not None and muts:
                    label = gene or extract_gene_from_filename(str(f)) or f.name
                    verbose_log.append(
                        f"{label}: processing failed ('{fail_msg}'); all {len(muts)} mutation(s) marked as failed."
                    )
                    for mut in muts:
                        verbose_log.append(
                            f"{label}: mutation {mut} skipped due to gene-level failure."
                        )
            elif verbose_log is not None:
                label = gene or extract_gene_from_filename(str(f)) or f.name
                verbose_log.append(f"{label}: processing failed ('{fail_msg}').")

    print(f"\nDone. Success: {ok}/{len(files)}")
    if fail:
        print("Failures:")
        for name, msg in fail:
            print(f"  - {name}: {msg}")

    if verbose_log:
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = output / f"validation_{timestamp}.log"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "w") as log_file:
                log_file.write("\n".join(verbose_log) + "\n")
            print(f"Verbose validation log written to {log_path}")
        except Exception as exc:
            print(f"Warning: Unable to write verbose log ({exc})", file=sys.stderr)

    if gene_metrics["total_genes"] > 0:
        print("\nGrand totals:")
        print(f"  Genes processed: {gene_metrics['total_genes']} (passed: {gene_metrics['passed_genes']}, failed: {gene_metrics['failed_genes']})")
        print(
            "  Mutations: {total} (passed: {passed}, failed: {failed}; length issues: {length}, "
            "mismatches: {mismatch}, dropped in mapping: {dropped})".format(
                total=gene_metrics["total_mutations"],
                passed=gene_metrics["passed_mutations"],
                failed=gene_metrics["failed_mutations"],
                length=gene_metrics["failed_length"],
                mismatch=gene_metrics["failed_mismatch"],
                dropped=gene_metrics["dropped_mapping"],
            )
        )
        if gene_metrics["splice_site_variants"]:
            print(
                f"  {gene_metrics['splice_site_variants']} splice-site variant(s):"
                f" REF span crosses an exon/intron boundary. Genomic, pre-mRNA and"
                f" intron rows written; transcript and amino-acid rows withheld"
                f" because the mature transcript is a splicing prediction."
            )
        if gene_metrics["dropped_coordinate_space"]:
            print(
                f"  WARNING: {gene_metrics['dropped_coordinate_space']} mutation(s) were"
                f" supplied in ORF space but span an exon-exon junction, which no DNA"
                f" variant does. They are almost certainly gDNA or chromosomal variants"
                f" written in the wrong coordinate space -- re-supply them as gd. tokens."
                f" Their transcript and amino-acid rows were still written."
            )
        if gene_metrics["dropped_aa_only"]:
            print(
                f"  Amino-acid rows missing for {gene_metrics['dropped_aa_only']} mutation(s)"
                f" that did receive nt mappings."
            )

    # Non-zero when any gene failed. The failures were printed but never
    # propagated: main() returned None and was called bare, so a run in which
    # EVERY gene failed exited 0 and any orchestrating shell script or Nextflow
    # process read it as success and carried on with absent or partial CSVs.
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
