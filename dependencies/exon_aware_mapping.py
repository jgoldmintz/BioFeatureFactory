#!/usr/bin/env python3
"""
Build exon-aware mappings + sequences for genes:
- FASTA per gene: ORF, transcript (spliced, mRNA orientation), genomic (tx_start..tx_end, genomic orientation)
- Mapping CSVs:
    * chromosome: mutant -> genomic-orientation "REF<abs_pos>ALT"
    * genomic (optional): mutant -> genomic-slice-relative "REF<gDNA_pos>ALT"
    * transcript (optional): mutant -> transcript-orientation "REF<tx_pos>ALT"
    * amino-acid (optional): mutant -> AA mutation "REF<aa_pos>ALT"
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path
import pysam
import warnings

# This script is placed alongside utility.py
sys.path.append(str(Path(__file__).parent))
from utility import (
    get_mutation_data,
    get_mutation_data_bioAccurate,
    get_mutant_aa,
    trim_muts,
    extract_gene_from_filename,
    get_genome_loc,
    read_fasta,
    write_fasta,
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


STOP_CODONS = {"TAA", "TAG", "TGA"}


def collect_orf_mutation_positions(orf_mutations: list[str]) -> tuple[list[int], dict[int, list[str]]]:
    positions: list[int] = []
    pos_to_mut: dict[int, list[str]] = {}
    for mut in orf_mutations:
        rel_pos, nts = get_mutation_data_bioAccurate(mut)
        if rel_pos is None or not nts:
            continue
        try:
            pos = int(rel_pos)
        except ValueError:
            continue
        positions.append(pos)
        pos_to_mut.setdefault(pos, []).append(mut)
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
    issues: list[str] = []
    mismatches: list[str] = []
    for mut in orf_mutations:
        rel_pos, nts = get_mutation_data_bioAccurate(mut)
        if rel_pos is None or not nts:
            continue
        wt_nt, _ = nts
        idx = int(rel_pos) - 1
        if idx < 0 or idx >= len(seq):
            issues.append(
                f"{gene_name}: mutation {mut} references position {rel_pos}, but ORF length is {len(seq)}."
            )
            continue
        if seq[idx] != wt_nt.upper():
            mismatches.append(
                f"{gene_name}: mutation {mut} expects '{wt_nt.upper()}' at ORF position {rel_pos},"
                f" but ORF sequence has '{seq[idx]}'."
            )
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

    out_of_range, mismatches = validate_mutations_against_orf(gene_name, orf_seq, orf_mutations)
    if out_of_range or mismatches:
        if verbose:
            for msg in out_of_range:
                emit_verbose(msg, True, log_messages)
            for msg in mismatches:
                emit_verbose(msg, True, log_messages)
            total = len(out_of_range) + len(mismatches)
            matched = max(len(orf_mutations) - total, 0)
            emit_verbose(
                f"{gene_name} validation summary:\n"
                f"  • Total mutations: {len(orf_mutations)}\n"
                f"  • Passed: {matched}\n"
                f"  • Failed: {total} (length: {len(out_of_range)}, base mismatch: {len(mismatches)})",
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
    """
    info = get_genome_loc(gene_name, annotation_file)
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

    Returns (tx_rows, chrom_rows, gdna_rows, aa_rows):
      tx_rows:     [(input_mut, "REF<tx_pos>ALT"), ...]
      chrom_rows: [(input_mut, "REF<abs_pos>ALT"), ...]
      gdna_rows:  [(input_mut, "REF<gdna_pos>ALT"), ...]
      aa_rows:    [(input_mut, "AA_REF<aa_pos>AA_ALT"), ...]
    """
    transcript_seq = tx_map["transcript_seq"]
    tx_to_genome   = tx_map["tx_to_genome"]
    cds_tx_start   = tx_map["cds_tx_start"]
    strand         = tx_map["strand"]
    orf_seq        = tx_map["orf_seq"]

    tx_rows, chrom_rows, gdna_rows, aa_rows = [], [], [], []

    for mut in orf_mutations:
        rel_pos, nts = get_mutation_data_bioAccurate(mut)
        if rel_pos is None or not nts:
            continue
        wt_nt, alt_nt = nts

        # Transcript position (1-based, mRNA orientation)
        tx_pos = cds_tx_start + int(rel_pos) - 1
        if tx_pos < 1 or tx_pos > len(transcript_seq):
            # skip out-of-bounds cleanly
            continue

        # Transcript-oriented mapping string
        tx_rows.append((mut, f"{wt_nt}{tx_pos}{alt_nt}"))

        # Absolute genomic coordinate
        gpos = tx_to_genome[tx_pos - 1]

        # Genomic-orientation nucleotides
        if strand == "+":
            g_wt, g_alt = wt_nt, alt_nt
        else:
            g_wt, g_alt = rc_base(wt_nt), rc_base(alt_nt)

        chrom_rows.append((mut, f"{g_wt}{gpos}{g_alt}"))

        # Relative coordinate within the gene's contiguous genomic slice (1-based)
        gdna_pos = gpos - tx_map["tx_start"] + 1
        gdna_rows.append((mut, f"{g_wt}{gdna_pos}{g_alt}"))

        # Amino-acid mapping via shared helper (handles validation/warnings)
        try:
            nt_mut = get_mutation_data(mut)
            aa_info = get_mutant_aa(nt_mut, orf_seq, aaseq=None, index=0)
            if aa_info:
                (aa_pos, (wt_aa, mut_aa)), _ = aa_info
                aa_rows.append((mut, f"{wt_aa}{aa_pos}{mut_aa}"))
        except Exception:
            continue

    return tx_rows, chrom_rows, gdna_rows, aa_rows


# CLI

def main():
    p = argparse.ArgumentParser(
        description="Generate exon-aware mapping CSVs and FASTAs (ORF/transcript/genomic) per gene."
    )
    p.add_argument("-m", "--mutations", required=True,
                   help="Path to a mutations CSV or a directory of CSVs (expects '<GENE>_mutations.csv').")
    p.add_argument("-r", "--reference", required=True, help="Reference genome FASTA (indexed).")
    p.add_argument("-a", "--annotation", required=True, help="Gene annotation file consumed by utility.get_genome_loc.")
    p.add_argument("--out-fasta", required=True, help="Output directory for gene FASTAs.")
    p.add_argument("--out-chromosome-mapping", required=True, help="Output directory for chromosome mapping CSVs (combined_<GENE>.csv with absolute genomic coordinates).")
    p.add_argument("--out-genomic-mapping", help="Optional output directory for gDNA mapping CSVs (relative to gene genomic slice).")
    p.add_argument("--out-transcript-mapping", help="Optional output dir for transcript mapping CSVs.")
    p.add_argument("--out-aa-mapping", help="Optional output dir for amino-acid mapping CSVs.")
    p.add_argument("--orf", help="Optional ORF FASTA (file or directory). If omitted, ORF is inferred from transcript.")
    p.add_argument("--verbose", action="store_true", help="Print detailed ORF/mutation validation messages.")
    args = p.parse_args()

    mut_path = Path(args.mutations)
    fasta_out = Path(args.out_fasta)
    chrom_out = Path(args.out_chromosome_mapping)
    gdna_out = Path(args.out_genomic_mapping) if args.out_genomic_mapping else None
    tmap_out = Path(args.out_transcript_mapping) if args.out_transcript_mapping else None
    aa_out = Path(args.out_aa_mapping) if args.out_aa_mapping else None
    aa_out = Path(args.out_aa_mapping) if args.out_aa_mapping else None

    capture_verbose = args.verbose #and mut_path.is_dir()
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
    }
    orf_lookup = {}
    if args.orf:
        try:
            orf_lookup = load_supplied_orfs(args.orf)
            print(f"Loaded {len(orf_lookup)} ORF sequence(s) from {args.orf}")
        except Exception as exc:
            print(f"Error loading supplied ORFs: {exc}", file=sys.stderr)
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

            # Load ORF mutations
            muts = trim_muts(str(f))
            if not muts:
                print("  (no mutations) -> skipping")
                continue

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

            # Build sequences + maps
            tx_map = build_transcript_seq_and_map(
                gene,
                args.annotation,
                args.reference,
                muts,
                supplied_orf=supplied_orf_seq,
                verbose=args.verbose,
                log_messages=verbose_log,
            )

            # Write FASTA with ORF / transcript / genomic
            write_fasta(
                fasta_out / f"{gene}.fasta",
                {
                    "ORF": tx_map["orf_seq"],
                    "transcript": tx_map["transcript_seq"],
                    "genomic": tx_map["gdna_seq"],
                },
            )
            print(f"  FASTA: {fasta_out / (gene + '.fasta')}")

            # Build mappings
            tx_rows, chrom_rows, gdna_rows, aa_rows = map_orf_mutations_to_transcript_and_genome(muts, tx_map)

            # Chromosome mapping CSV (required)
            write_mapping_csv(chrom_out / f"chr_mapping_{gene}.csv", "chromosome", chrom_rows)
            print(f"  Chromosome mapping: {chrom_out / ('chr_mapping_' + gene + '.csv')}  ({len(chrom_rows)} rows)")

            # Genomic (relative gDNA) mapping CSV (optional)
            if gdna_out:
                write_mapping_csv(gdna_out / f"genomic_mapping_{gene}.csv", "genomic", gdna_rows)
                print(f"  Genomic mapping: {gdna_out / ('genomic_mapping_' + gene + '.csv')}  ({len(gdna_rows)} rows)")

            # Transcript mapping CSV (optional)
            if tmap_out:
                write_mapping_csv(tmap_out / f"transcript_mapping_{gene}.csv", "transcript", tx_rows)
                print(f"  Transcript mapping: {tmap_out / ('transcript_mapping_' + gene + '.csv')}  ({len(tx_rows)} rows)")

            # Amino-acid mapping CSV (optional)
            if aa_out:
                write_mapping_csv(aa_out / f"{gene}_nt_to_aa_mapping.csv", "aamutant", aa_rows)
                print(f"  Amino-acid mapping: {aa_out / (gene + '_nt_to_aa_mapping.csv')}  ({len(aa_rows)} rows)")
            len_issues = tx_map.get("validation_length_issues") or []
            mismatches = tx_map.get("validation_mismatches") or []
            total_mut = len(muts)
            failed_mut = len(len_issues) + len(mismatches)
            passed_mut = max(total_mut - failed_mut, 0)
            gene_metrics["total_mutations"] += total_mut
            gene_metrics["passed_mutations"] += passed_mut
            gene_metrics["failed_mutations"] += failed_mut
            gene_metrics["failed_length"] += len(len_issues)
            gene_metrics["failed_mismatch"] += len(mismatches)
            if len_issues or mismatches:
                err_msg = (
                    f"{gene}: ORF/mutation validation failed for "
                    f"{len(len_issues) + len(mismatches)} mutation(s)."
                )
                print(f"  ERROR: {err_msg}", file=sys.stderr)
                fail.append((f.name, err_msg))
                gene_metrics["failed_genes"] += 1
                continue

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
        log_path = fasta_out / f"validation_{timestamp}.log"
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
            "  Mutations: {total} (passed: {passed}, failed: {failed}; length issues: {length}, mismatches: {mismatch})".format(
                total=gene_metrics["total_mutations"],
                passed=gene_metrics["passed_mutations"],
                failed=gene_metrics["failed_mutations"],
                length=gene_metrics["failed_length"],
                mismatch=gene_metrics["failed_mismatch"],
            )
        )

if __name__ == "__main__":
    main()
