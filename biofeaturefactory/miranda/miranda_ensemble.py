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
Unified MirandA pipeline (WT + MUT in one run) with comparative outputs.

Key changes
- WT: run MirandA once per gene -> {GENE}-wt-miranda.out (no per-mutation duplication).
- MUT: per-gene parallel execution; logs every 100 processed with true done/total seeded from disk.
- Parser: joins WT and MUT by (pkey, mirna_id, locus_id); WT rows are materialized per mutation from the single per-gene WT output.
- Intermediate .out files are deleted after parsing to reclaim disk space.
- Redundancy reduction: leverage utility.load_mapping, extract_mutation_from_sequence_name, extract_gene_from_filename, parse_variant, splice_seq, align_wt_to_mut, load_wt_sequences, load_validation_failures, should_skip_mutation.

Variant classes
- Every token the nt grammar accepts is processed BY DEFAULT: SNV, MNV, insertion,
  deletion and delins. There is no flag. The grammar is uniquely decodable, so no
  parser mode has to be selected, and whether a given column survives is decided
  per record from len(ref) != len(alt) -- a fact of the mutation, not a user
  preference.
- The WT/MUT join is made on a PROJECTED coordinate (`join_pos`), not on raw
  position equality. Under a length change, position i of the WT and position i of
  the MUT are different bases, so the old pooled proximity chain split one physical
  site into a phantom lost+gained pair. `align_wt_to_mut` maps the WT anchor into
  the mutant frame; for len(ref) == len(alt) that map is the identity, so the SNV
  path is unchanged by construction rather than by a branch.
- A token that cannot be turned into a mutant sequence is COUNTED AND NAMED in
  miranda.run_summary.json. It is never dropped silently.
"""

import os
import sys
import csv
import argparse
import math
import json
import re
import concurrent.futures
import multiprocessing
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterable

import pandas as pd

# -------------------------------------------------------------------------
# utility imports
# -------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
from biofeaturefactory.lib.utility import (
    parse_piece_token,
    mint_pkey,
    piece_fields,
    split_piece_cell,
    load_validation_failures,
    should_skip_mutation,
    extract_gene_from_filename,
    extract_mutation_from_sequence_name,
    load_wt_sequences,
    load_mapping,
    parse_variant,
    splice_seq,
    align_wt_to_mut,
    read_fasta,
)

# -------------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------------
VISIBILITY_THRESHOLD = 140.0
HIGH_CUTOFF = 150.0
MERGE_WINDOW_NT = 15
SEGMENT_WINDOW_NT = 25
SHIFT_NT = 4
REPORT_RADIUS = 40
DISTANCE_K = 25.0
MAX_PARALLEL_CAP = 16
GENE_PROGRESS_STEP = 100
PROGRESS_SAVE_EVERY = 100  # persist progress every N completions

# -------------------------------------------------------------------------
# UTILS
# -------------------------------------------------------------------------
# A run is keyed on (gene, substrate). The transcript substrate keeps the BARE
# gene name, so every existing filename, bucket key and pkey is byte-identical to
# before -- only the new substrates carry a suffix.
#
# '~' is the separator because it appears in neither a HGNC symbol (which admits
# letters, digits and '-') nor a variant token (ACGTU, digits, and the 'gd.'
# prefix). '-' could not be used: _collect_gene_outputs recovers the gene by
# rsplit('-', 1), which already has to tolerate hyphenated symbols like HLA-A.
SUBSTRATE_SEP = "~"


# The key must be UPPER-CASE-IDEMPOTENT: run_wt_phase and run_mut_phase both do
# gene.upper() before looking the sequence up, so a key carrying a lower-case
# substrate ('TESTG~intron1') became 'TESTG~INTRON1' and missed every lookup --
# observed as "Missing WT sequence" for all three substrates while the transcript
# ran fine. The display name is recovered on the way out instead.
_SUBSTRATE_DISPLAY = {"PRE_MRNA": "pre_mRNA"}


def _gene_key(gene: str, substrate: str) -> str:
    g = gene.upper()
    return g if substrate == "transcript" else f"{g}{SUBSTRATE_SEP}{substrate.upper()}"


def _split_gene_key(key: str) -> Tuple[str, str]:
    if SUBSTRATE_SEP in key:
        gene, substrate = key.split(SUBSTRATE_SEP, 1)
        up = substrate.upper()
        return gene, _SUBSTRATE_DISPLAY.get(up, up.lower())
    return key, "transcript"


def _substrate_of(key: str) -> str:
    return _split_gene_key(key)[1]


def safe_div(num: float, den: float) -> float:
    return num / den if den not in (0, None, 0.0) else 0.0

def exp_weight(d: Optional[float], k: float) -> float:
    if d is None:
        return 1.0
    return math.exp(-(d / k)) if k else 1.0

def distance_to_variant(site_pos: Optional[int],
                        var_pos1: Optional[int],
                        span_len: Optional[int]) -> Optional[int]:
    """1-based distance from a site anchor to the variant's ALLELE SPAN.

    Zero inside the span, otherwise the distance to the nearer edge. The single
    -base form this replaces (`abs(site_pos - var_pos1)`) measured to the FIRST
    base of the span, so for a multi-base allele the reported distance grew with
    the allele's own length -- a 50 nt deletion put every site inside itself 0..49
    nt "away". For span_len == 1 the two forms are identical, which is why the SNV
    path is unaffected.

    span_len is len(REF) for a WT-frame site and len(ALT) for a MUT-frame one:
    bases 5' of the edit share a coordinate in both frames, so the span starts at
    var_pos1 in either, but its length is whichever allele that frame carries.
    """
    if site_pos is None or var_pos1 is None:
        return None
    span_end = var_pos1 + (span_len if span_len else 1) - 1
    if site_pos < var_pos1:
        return var_pos1 - site_pos
    if site_pos > span_end:
        return site_pos - span_end
    return 0

def _existing_mut_outputs(outdir: str, gene_upper: str) -> set[str]:
    """Return seq_ids like GENE-<mut>-mut that already exist on disk."""
    existing = set()
    prefix = f"{gene_upper}-"
    suffix = "-mut-miranda.out"
    try:
        for name in os.listdir(outdir):
            if not (name.startswith(prefix) and name.endswith(suffix)):
                continue
            # name: GENE-<mut>-mut-miranda.out
            core = name[:-len(suffix)]         # GENE-<mut>-mut
            mut = core[len(prefix):-len("-mut")]
            existing.add(f"{gene_upper}-{mut}-mut")
    except FileNotFoundError:
        pass
    return existing

def _write_gene_progress(outdir: str, gene: str, processed: int, total: int, completed_ids: Iterable[str]):
    p = os.path.join(outdir, f"{gene}.mut.progress.json")
    data = {
        "gene": gene,
        "processed": processed,
        "total": total,
        "completed_ids": sorted(set(completed_ids))[:5000],
    }
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f)
    os.replace(tmp, p)

# -------------------------------------------------------------------------
# MAPPING LOAD
# -------------------------------------------------------------------------
def load_transcript_mappings(mapping_dir: str, map_type: str = "transcript") -> Dict[str, pd.DataFrame]:
    path = Path(mapping_dir)
    if not path.exists():
        raise FileNotFoundError(f"Mapping path not found: {mapping_dir}")

    out: Dict[str, pd.DataFrame] = {}
    if path.is_file():
        files = [path] if path.suffix.lower() in (".csv", ".tsv", ".txt") else []
    else:
        files = sorted(f for f in path.rglob("*") if f.is_file() and f.suffix.lower() in (".csv", ".tsv", ".txt"))
    for f in files:
        # F36: rglob("*") + suffix filter admits .csv/.tsv/.txt (was rglob("*.csv"),
        # which silently dropped .tsv/.txt mappings). F37: no local pre-validation —
        # load_mapping now self-guards and returns {} on an invalid/single-column file,
        # which the `if not mp` skip already handles. (The old try/except was a no-op:
        # validate_mapping_content returns False rather than raising, so it caught
        # nothing and the ignored False later reached load_mapping's unpack → TypeError.)
        mp = load_mapping(str(f), mapType=map_type)
        if not mp:
            continue

        df = pd.DataFrame(
            [(str(k).strip(), str(v).strip()) for k, v in mp.items()],
            columns=["mutant", map_type],
        ).dropna()

        if map_type != "transcript":
            df["transcript"] = df[map_type]
        cols = ["mutant", "transcript"] + ([map_type] if map_type != "transcript" else [])
        df = df[cols].drop_duplicates(subset=["mutant", "transcript"]).reset_index(drop=True)

        gene = (extract_gene_from_filename(str(f)) or f.stem).upper()
        out[gene] = (pd.concat([out[gene], df], ignore_index=True)
                     if gene in out else df)
    return out

# -------------------------------------------------------------------------
# MIRANDA EXECUTION
# -------------------------------------------------------------------------
def _run_single_miranda_task(args: Tuple[str, str, str, str, str, bool]) -> Tuple[str, str]:
    """
    args = (seq_id, seq_string, miranda_dir, outdir, mirna_db_abs, strict)

    miranda_dir may be empty string/None to use miranda from PATH.

    `strict` maps to miranda's own -strict ("Demand strict 5' seed pairing",
    default off). It is intended for intronic substrates: an intron is long and
    unfiltered relative to a 3'UTR, so the default permissive seed model returns
    a large low-confidence hit set there.
    """
    import shutil as _shutil
    import tempfile as _tempfile
    from subprocess import run, TimeoutExpired
    seq_id, seq_str, miranda_dir, outdir, mirna_db_abs, strict = args

    # Resolve binary and working directory
    if miranda_dir:
        if not os.path.isdir(miranda_dir):
            return (seq_id, "error:miranda_dir_missing")
        bin_path = os.path.join(miranda_dir, "miranda")
        if not (os.path.isfile(bin_path) and os.access(bin_path, os.X_OK)):
            return (seq_id, "error:miranda_binary_not_exec")
        cmd_binary = "./miranda"
        cwd = miranda_dir
        tmp_base = miranda_dir
    else:
        if not _shutil.which("miranda"):
            return (seq_id, "error:miranda_not_on_PATH")
        cmd_binary = "miranda"
        cwd = None
        tmp_base = _tempfile.gettempdir()

    if not os.path.isfile(mirna_db_abs):
        return (seq_id, "error:mirna_db_missing")

    # Skip if already produced
    out_file = os.path.join(outdir, f"{seq_id}-miranda.out")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return (seq_id, "already")

    # F40: atomically-unique temp path. The prior tmp_{pid}_{hash&0xFFFFF} scheme
    # collided across ThreadPoolExecutor threads (constant pid) on any 20-bit hash
    # clash, racing open/run/unlink between two concurrent mutations.
    fd, tmp_file = _tempfile.mkstemp(prefix="tmp_", suffix=".fasta", dir=tmp_base)
    os.close(fd)
    try:
        with open(tmp_file, "w") as f:
            f.write(f">{seq_id}\n{seq_str}")

        # When using cwd, miranda expects basename; when on PATH, use absolute path
        fasta_arg = os.path.basename(tmp_file) if miranda_dir else tmp_file
        cmd = [cmd_binary, mirna_db_abs, fasta_arg]
        if strict:
            cmd.append("-strict")
        result = run(cmd, cwd=cwd, capture_output=True, text=True, timeout=1200)

        # Always write stdout for inspection
        try:
            with open(out_file, "w") as f:
                f.write(result.stdout or "")
        except Exception as e:
            return (seq_id, f"error:write_out:{e}")

        # Miranda (v1.9) crashes in its exit cleanup handler (SIGBUS) on macOS
        # but produces fully valid scan output before that. Treat any run
        # where stdout was produced as successful.
        if result.stdout:
            status = "ok"
        else:
            status = f"rc={result.returncode}"
    except TimeoutExpired:
        status = "timeout"
    except Exception as e:
        status = f"error:{e}"
    finally:
        try:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        except Exception:
            pass

    return (seq_id, status)

def run_wt_phase(wt_sequences: Dict[str, str],
                 outdir: str,
                 miranda_dir: str,
                 mirna_db: str,
                 strict_substrates: Optional[set] = None):
    """
    One WT file per gene: {GENE}-wt-miranda.out
    miranda_dir may be empty string/None to use miranda from PATH.
    """
    import tempfile as _tempfile
    os.makedirs(outdir, exist_ok=True)
    print("[WT] Starting WT MirandA phase...")

    genes = sorted(wt_sequences.keys())
    total_genes = len(genes)
    print(f"[WT] Loaded {total_genes} WT transcripts")

    from subprocess import run

    db_abs = os.path.abspath(mirna_db)

    if miranda_dir:
        cmd_binary = "./miranda"
        cwd = miranda_dir
        tmp_base = miranda_dir
    else:
        cmd_binary = "miranda"
        cwd = None
        tmp_base = _tempfile.gettempdir()

    for idx, gene in enumerate(genes, 1):
        gene_upper = gene.upper()
        out_file = os.path.join(outdir, f"{gene_upper}-wt-miranda.out")

        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            print(f"[WT] [{idx}/{total_genes}] {gene_upper}: exists on disk, skip")
            continue

        print(f"[WT] [{idx}/{total_genes}] Running WT MirandA: {gene_upper}")

        wt_seq = wt_sequences.get(gene_upper)
        if not wt_seq:
            print(f"[WT]   Missing WT sequence for {gene_upper}, skip")
            continue

        wt_tmp_file = os.path.join(tmp_base, f"tmp_{gene_upper}_WT.fasta")
        with open(wt_tmp_file, "w") as f:
            f.write(f">transcript\n{wt_seq}")

        fasta_arg = os.path.basename(wt_tmp_file) if miranda_dir else wt_tmp_file
        cmd = [cmd_binary, db_abs, fasta_arg]
        # WT and MUT must be scanned under IDENTICAL settings, or every delta is
        # partly an artefact of the two runs using different seed models.
        if _substrate_of(gene_upper) in (strict_substrates or set()):
            cmd.append("-strict")
        result = run(cmd, cwd=cwd, capture_output=True, text=True, timeout=3000)
        miranda_output = result.stdout

        try:
            os.remove(wt_tmp_file)
        except Exception:
            pass

        with open(out_file, "w") as f:
            f.write(miranda_output)

        print(f"[WT]   Wrote {out_file} (len={len(miranda_output)})")

    print("[WT] Done.")

def _per_gene_mut_tasks(gene_upper: str,
                        wt_seq: str,
                        mapping_df: pd.DataFrame,
                        outdir: str,
                        miranda_dir: str,
                        mirna_db: str,
                        failure_map: Optional[Dict],
                        already_on_disk: set,
                        strict: bool = False) -> Tuple[List[Tuple[str, str, str, str, str, bool]], List[Dict]]:
    """
    Build per-gene list of MirandA tasks for MUT allele, skipping disk.

    Returns (tasks, rejected). `rejected` names every token that produced no task
    for a reason other than "the output already exists": the previous version
    `continue`d on all of them, including every non-SNV token, so an indel was
    indistinguishable from a mutation the gene never had.

    Non-SNV tokens are built here by default. `parse_variant` is what makes that
    possible: `get_mutation_data` raised the SAME ValueError for a valid indel and
    for a corrupt token, so the try/except it was wrapped in could not tell them
    apart and had to discard both.
    """
    tasks: List[Tuple[str, str, str, str, str, bool]] = []
    rejected: List[Dict] = []
    db_abs = os.path.abspath(mirna_db)

    def _reject(mutant_id, token, reason):
        rejected.append({"gene": gene_upper, "mutant": mutant_id,
                         "transcript_token": token, "reason": reason})

    # Prebuild lookup once per gene (no per-file DataFrame filtering)
    mp = dict(zip(mapping_df["mutant"].astype(str), mapping_df["transcript"].astype(str)))

    for mutant_id, transcript_token in mp.items():
        mutant_id = str(mutant_id).strip()
        transcript_token = str(transcript_token).strip()
        if not mutant_id or not transcript_token:
            _reject(mutant_id, transcript_token, "empty_mapping_row")
            continue

        # The seq_id becomes the OUTPUT FILENAME ({seq_id}-miranda.out), so the
        # token cannot go in it verbatim: a knockout token is unbounded and
        # overran the 224-char filename cap as [Errno 63] File name too long.
        # mint_pkey is 18 chars regardless of allele size. The header written to
        # the temp FASTA is this same id, and no parser reads it back
        # (parse_miranda_text ignores 'Read Sequence:'), so the id's only job is
        # to name the file -- and it is now the pkey the rest of the repo joins on.
        seq_id = f"{mint_pkey(gene_upper, mutant_id)}-mut"
        if seq_id in already_on_disk:
            continue

        if failure_map and should_skip_mutation(gene_upper, mutant_id, failure_map):
            _reject(mutant_id, transcript_token, "validation_log_filtered")
            continue

        # Which parser depends on the SUBSTRATE, because the two mappings speak
        # different grammars. The transcript mapping carries whole variants
        # (REF<pos>ALT, both alleles non-empty). The intron/pre-mRNA mapping
        # carries decomposed PIECES of one edit, and a piece that is wholly
        # deleted is written 'REF<pos>del' -- a form parse_variant refuses by
        # design, because Variant cannot hold an empty allele. Sending those
        # through parse_variant rejected them as "unparseable_token": the
        # deletion that removes a splice site, which is precisely the variant
        # class the intron substrate exists to score, was the one class never
        # scored on it. splice_seq accepts an empty ALT and produces the deletion
        # directly (verified: 'GTATTG5del' on a 14 nt sequence -> 8 nt).
        if _substrate_of(gene_upper) == "transcript":
            variant = parse_variant(transcript_token, is_nt=True)
            if variant is None:
                _reject(mutant_id, transcript_token, "unparseable_token")
                continue
            v_pos0, v_ref, v_alt = variant.pos0, variant.ref, variant.alt
        else:
            piece = parse_piece_token(transcript_token)
            if piece is None:
                _reject(mutant_id, transcript_token, "unparseable_piece_token")
                continue
            v_pos0, v_ref, v_alt = piece["pos"] - 1, piece["ref"], piece["alt"]

        # Bound the END of the REF span, not its start: a multi-base REF can begin
        # inside the transcript and run off the 3' end.
        span_end = v_pos0 + len(v_ref)
        if not (0 <= v_pos0 and span_end <= len(wt_seq)):
            _reject(mutant_id, transcript_token,
                    f"position_out_of_range:span_{v_pos0}..{span_end}_of_{len(wt_seq)}")
            continue

        # Compare the WHOLE REF span. `wt_seq[pos] != ref[0]` passed any multi-base
        # REF whose first base happened to match, splicing a mutant that does not
        # correspond to the token.
        observed = wt_seq[v_pos0:span_end].upper()
        if observed != v_ref.upper():
            _reject(mutant_id, transcript_token, f"REF_MISMATCH:transcript_has_{observed}")
            continue

        # splice_seq honours len(ref); update_str hardcodes a one-base stride and
        # cannot express an indel. REF was just verified, so skip the re-check.
        mut_seq = splice_seq(wt_seq, v_pos0, v_ref, v_alt, validate=False)
        out_file = os.path.join(outdir, f"{seq_id}-miranda.out")
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            continue

        tasks.append((seq_id, mut_seq, miranda_dir, outdir, db_abs, strict))
    return tasks, rejected

def run_mut_phase(wt_sequences: Dict[str, str],
                  mapping_dict: Dict[str, pd.DataFrame],
                  outdir: str,
                  miranda_dir: str,
                  mirna_db: str,
                  use_parallel: bool = True,
                  max_workers: Optional[int] = None,
                  failure_map: Optional[Dict] = None,
                  strict_substrates: Optional[set] = None) -> Dict:
    """
    MUT: per-gene parallel; progress prints seeded with prior completions.

    Returns the run-level accounting: every token that yielded no MirandA task,
    named with its reason, plus the run/failure tallies. Written to
    miranda.run_summary.json by main().
    """
    print("[MUT] Starting MUT MirandA phase...")
    os.makedirs(outdir, exist_ok=True)

    mut_summary: Dict = {
        "genes_total": 0, "genes_skipped": 0,
        "tokens_total": 0, "tokens_built": 0, "tokens_rejected": 0,
        "runs_ok": 0, "runs_failed": 0,
        "rejected": [], "run_failures": [], "skipped_genes": [],
    }

    # A gene with a WT sequence but no usable mapping was previously dropped by
    # this comprehension without a word -- no print, no counter, no output file --
    # so "gene absent from the results" meant either "no mapping" or "no hits" and
    # nothing distinguished them. Same set, same order, now named.
    genes = []
    for g in sorted(wt_sequences.keys()):
        mdf = mapping_dict.get(g)
        if mdf is None:
            mut_summary["skipped_genes"].append({"gene": g, "reason": "no_mapping_file"})
        elif mdf.empty:
            mut_summary["skipped_genes"].append({"gene": g, "reason": "empty_mapping"})
        else:
            genes.append(g)
    total_genes = len(genes)
    mut_summary["genes_total"] = len(wt_sequences)
    mut_summary["genes_skipped"] = len(mut_summary["skipped_genes"])
    if mut_summary["skipped_genes"]:
        print(f"[MUT] {mut_summary['genes_skipped']} gene(s) skipped for lack of a mapping: "
              + ", ".join(s["gene"] for s in mut_summary["skipped_genes"][:10]))

    total_rows = sum(len(mapping_dict[g]) for g in genes)
    mut_summary["tokens_total"] = total_rows
    print(f"[MUT] Prepared estimates using mapping rows: {total_rows}")

    if max_workers is None:
        max_workers = min(max(multiprocessing.cpu_count() - 1, 1), MAX_PARALLEL_CAP)
    print(f"[MUT] Using {max_workers} workers (per batch)")

    for idx, gene in enumerate(genes, 1):
        gene_upper = gene.upper()
        print(f"[MUT] [{idx}/{total_genes}] Starting {gene_upper} processing ...")

        # Both of the guards below are unreachable from main(): load_wt_sequences
        # upper-cases every key (utility.py:1980) so gene == gene_upper and the
        # lookup cannot miss, and `genes` above was already filtered on a non-empty
        # mapping. They are pre-existing and left in place for direct callers; the
        # gene-level skip accounting lives at the filter, where it can actually fire.
        wt_seq = wt_sequences.get(gene_upper)
        if not wt_seq:
            print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: no WT sequence, skip")
            continue

        mapping_df = mapping_dict.get(gene_upper)
        if mapping_df is None or mapping_df.empty:
            print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: no mapping, skip")
            continue

        # Prior completions
        all_mut_ids = [str(m).strip() for m in mapping_df["mutant"].dropna().astype(str)]
        total_all = len(all_mut_ids)
        already_on_disk = _existing_mut_outputs(outdir, gene_upper)
        done_prior = len(already_on_disk)

        tasks, rejected = _per_gene_mut_tasks(
            gene_upper=gene_upper,
            wt_seq=wt_seq,
            mapping_df=mapping_df,
            outdir=outdir,
            miranda_dir=miranda_dir,
            mirna_db=mirna_db,
            failure_map=failure_map,
            already_on_disk=already_on_disk,
            strict=_substrate_of(gene_upper) in (strict_substrates or set()),
        )
        remaining = len(tasks)
        mut_summary["tokens_built"] += remaining
        mut_summary["tokens_rejected"] += len(rejected)
        mut_summary["rejected"].extend(rejected)
        if rejected:
            by_reason: Dict[str, int] = {}
            for r in rejected:
                # Group on the reason CLASS; REF_MISMATCH and position_out_of_range
                # carry the offending detail after ':' and would otherwise print one
                # line per mutation.
                by_reason[r["reason"].split(":", 1)[0]] = by_reason.get(r["reason"].split(":", 1)[0], 0) + 1
            print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: rejected {len(rejected)} token(s): "
                  + ", ".join(f"{k}={v}" for k, v in sorted(by_reason.items())))

        if remaining == 0:
            print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: already complete [{done_prior}/{total_all}]")
            continue

        print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: remaining {remaining} of {total_all} (done={done_prior})")

        processed = done_prior
        completed_ids: List[str] = []

        def _record(seq_id: str, status: str):
            """Tally one completed MirandA run and emit the periodic progress lines.

            A non-ok status used to be printed and then forgotten; it is now also
            counted and named, so a gene whose runs all failed is distinguishable
            from a gene that had nothing to run. Shared by the serial and parallel
            branches so the two cannot drift.
            """
            nonlocal processed
            if status not in ("ok", "already"):
                print(f"[MUT]   {seq_id}: {status}")
                mut_summary["runs_failed"] += 1
                mut_summary["run_failures"].append(
                    {"gene": gene_upper, "seq_id": seq_id, "status": status})
            else:
                completed_ids.append(seq_id)
                mut_summary["runs_ok"] += 1
            processed += 1
            if processed % GENE_PROGRESS_STEP == 0 or processed == total_all:
                print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: Processed mutations - [{processed}/{total_all}]")
            if processed % PROGRESS_SAVE_EVERY == 0 or processed == total_all:
                _write_gene_progress(outdir, gene_upper, processed, total_all, completed_ids)

        if not use_parallel:
            for t in tasks:
                seq_id, status = _run_single_miranda_task(t)
                _record(seq_id, status)
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                # map with chunksize=1 for steady progress feedback
                for seq_id, status in ex.map(_run_single_miranda_task, tasks):
                    _record(seq_id, status)

    print("[MUT] Done.")
    print(f"[MUT] Built {mut_summary['tokens_built']} mutant sequence(s); "
          f"rejected {mut_summary['tokens_rejected']} token(s); "
          f"{mut_summary['runs_failed']} run failure(s)")
    return mut_summary

# -------------------------------------------------------------------------
# PARSING
# -------------------------------------------------------------------------
def parse_miranda_text(miranda_text: str) -> List[Dict]:
    """
    Parse miRanda stdout -> per-hit site rows.

    F35: emit one row per INDIVIDUAL hit, read from the single-'>' "Scores for
    this hit" lines, which carry that site's own score/energy/position. The '>>'
    summary line is a transcript-wide aggregate (Tot Score = SUM over all hits,
    Max Score = best single hit); the previous code exploded it across the
    position field, stamping the aggregate onto every site (e.g. tot=473
    replicated on the three SMN2/miR-608 sites at 589/673/736 instead of each
    site's own score).

    '>' layout after whitespace->tab:
      [0]>mirna [1]target [2]score [3]energy [4]_ [5]q_start [6]q_end
      [7]ref_start [8]ref_end [9]aln_len [10]%id [11]%sim
    """
    sites: List[Dict] = []
    query_seq = ""
    ref_seq = ""
    # strand / len_mirna / len_target are transcript-wide descriptors that miRanda
    # only prints on the '>>' summary (Strand=parts[6], Len1=parts[7], Len2=parts[8]).
    # They are collected per miRNA here and stamped onto that miRNA's hits after the
    # scan. An earlier revision inferred lengths from the 'Read Sequence:' lines by
    # testing `"transcript" in line`, which holds only for the WT FASTA (run_wt_phase
    # writes '>transcript'); MUT files are headed '>{GENE}-{MUT}-mut', so every MUT
    # row got len_mirna = transcript length and len_target = None.
    summary_by_mirna: Dict[str, Dict] = {}

    for raw in miranda_text.splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("Read Sequence:"):
            continue

        if "Query:" in line:
            # Alignment display line — extract sequence for query_seq only.
            seg = line.split("Query:", 1)[1].strip()
            if "'" in seg:
                parts = seg.split("'")
                token = parts[1].strip() if len(parts) >= 2 else seg.split()[0]
            else:
                token = seg.split()[0]
            query_seq = token.replace("5", "").replace("3", "").replace(" ", "")
            continue

        if "Ref:" in line:
            seg = line.split("Ref:", 1)[1].strip()
            if "'" in seg:
                parts = seg.split("'")
                token = parts[1].strip() if len(parts) >= 2 else seg
            else:
                token = seg
            ref_seq = token.replace("5", "").replace("3", "").replace(" ", "")
            continue

        # Transcript-wide summary. Per-site scores/positions are NOT taken from here
        # (Tot Score is a SUM over hits), but Strand/Len1/Len2 are only available here.
        if line.startswith(">>"):
            sparts = re.sub(r"\s+", "\t", line).split("\t")
            if len(sparts) >= 9:
                sid = sparts[0][2:]
                if sid:
                    def _i(v):
                        try: return int(v)
                        except Exception: return None
                    summary_by_mirna[sid] = {
                        "strand": sparts[6],
                        "len_mirna": _i(sparts[7]),
                        "len_target": _i(sparts[8]),
                    }
            continue

        if line.startswith(">"):
            norm = re.sub(r"\s+", "\t", line)
            parts = norm.split("\t")
            if len(parts) < 9:
                continue

            mirna_id = parts[0][1:]  # strip the single '>'
            if not mirna_id:
                continue

            def f(i):
                try: return float(parts[i])
                except Exception: return None

            def i_(i):
                try: return int(parts[i])
                except Exception: return None

            score     = f(2)   # THIS hit's score (was the '>>' transcript sum)
            energy    = f(3)
            ref_start = i_(7)  # R: start = 1-based site position
            if ref_start is None:
                continue

            sites.append({
                "mirna_id": mirna_id,
                "site_pos": ref_start,
                "tot_score": score,
                "tot_energy": energy,
                "max_score": score,     # single hit → its own score is its max
                "max_energy": energy,
                "strand": "",
                "len_mirna": None,
                "len_target": None,
                "query_seq": query_seq,
                "ref_seq": ref_seq,
                "parser_confidence": 1.0 if (mirna_id and score is not None) else 0.5,
            })

    # Stamp the '>>' descriptors onto each hit (the summary follows its hits, so this
    # is done after the scan rather than inline).
    for s in sites:
        meta = summary_by_mirna.get(s["mirna_id"])
        if meta:
            s["strand"] = meta["strand"]
            s["len_mirna"] = meta["len_mirna"]
            s["len_target"] = meta["len_target"]
    return sites

def _collect_gene_outputs(outdir: str) -> Tuple[Dict[str, Dict[str, List[Path]]], List[Dict]]:
    """
    ({ GENE: { 'wt': [Path], 'mut': [List[Path]] } }, unrecognised)

    WT:  GENE-wt-miranda.out
    MUT: GENE-<mut>-mut-miranda.out

    `unrecognised` names every *-miranda.out this function declines to bucket.
    A file that reaches neither branch used to vanish here without a print, a
    counter or a row -- and the caller's own "unrecognised_filename" rejection
    could not compensate, because it re-tested the identical predicate on a file
    this filter had already discarded, so it was unreachable. The drop is named
    where it actually happens.
    """
    buckets: Dict[str, Dict[str, List[Path]]] = {}
    unrecognised: List[Dict] = []
    for f in os.scandir(outdir):
        if not f.is_file() or not f.name.endswith("-miranda.out"):
            continue
        base = f.name[:-len("-miranda.out")]

        if base.endswith("-wt"):
            gene = base[:-3].upper()
            buckets.setdefault(gene, {}).setdefault("wt", []).append(Path(f.path))
        elif base.endswith("-mut"):
            core = base[:-4]  # GENE-<mut>
            if "-" not in core:
                unrecognised.append({"file": f.name, "reason": "unrecognised_filename"})
                continue
            # rsplit on the LAST hyphen (F547): the gene may itself contain hyphens
            # (HLA-A, NKX2-1). split("-",1) truncated "HLA-A-<mut>" to "HLA", bucketing
            # MUT files under a different key than the WT branch's suffix-strip, silently
            # losing all MUT/delta for hyphenated genes. Mutant nt tokens have no hyphen.
            gene = core.rsplit("-", 1)[0].upper()
            buckets.setdefault(gene, {}).setdefault("mut", []).append(Path(f.path))
        else:
            unrecognised.append({"file": f.name, "reason": "unrecognised_filename"})
    return buckets, unrecognised

def build_sites_table_from_outputs(outdir: str,
                                   mapping_dict: Dict[str, pd.DataFrame],
                                   failure_map: Dict) -> Tuple[pd.DataFrame, List[Dict]]:
    """
    Build sites by joining single per-gene WT output with each mutation-specific MUT output.
    WT rows are materialized per pkey (allele='WT'); MUT rows parsed from files (allele='MUT').
    Intermediate .out files are deleted after parsing to reclaim disk space.

    Returns (sites_df, rejected). Every MUT output that produced no rows is named
    in `rejected` with a reason, instead of only bumping an anonymous counter.

    Each row carries `join_pos`, the coordinate the WT/MUT join is made on. For a
    WT row that is the WT anchor PROJECTED INTO THE MUTANT FRAME by
    align_wt_to_mut; for a MUT row it is the site position itself. Joining on raw
    `site_pos` is only correct when the two alleles share a coordinate system,
    which a length-changing variant breaks: every site 3' of the edit is renumbered
    by length_delta, so the proximity chain split one physical site into a WT-only
    and a MUT-only locus once |length_delta| exceeded MERGE_WINDOW_NT. For
    len(ref) == len(alt) the projection is the identity, so SNV and MNV rows get
    join_pos == site_pos and nothing about their handling changes.
    """
    print("[PARSE] Building sites table from MirandA outputs...")
    records: List[Dict] = []
    rejected: List[Dict] = []

    buckets, unrecognised = _collect_gene_outputs(outdir)
    rejected.extend(unrecognised)
    genes = sorted(buckets.keys())
    print(f"[PARSE] Found outputs for {len(genes)} genes")

    miss_map = miss_hits = 0
    parsed_files = 0

    for gi, gene_key in enumerate(genes, 1):
        # buckets are keyed on (gene, substrate); pkey and the mapping lookup use
        # the REAL gene so an intronic row joins the same mutation the transcript
        # rows do.
        gene, substrate = _split_gene_key(gene_key)
        m = buckets[gene_key]
        wt_files = m.get("wt", [])
        mut_files = m.get("mut", [])
        if not wt_files:
            print(f"[PARSE]   {gene_key}: missing WT output; skipping {len(mut_files)} MUT files")
            for mf in mut_files:
                rejected.append({"gene": gene, "file": mf.name, "reason": "missing_wt_output"})
            continue
        wt_file = wt_files[0]

        wt_text = wt_file.read_text(encoding="utf-8", errors="ignore")
        wt_hits = parse_miranda_text(wt_text)

        mapping_df = mapping_dict.get(gene_key)
        if mapping_df is None or mapping_df.empty:
            print(f"[PARSE]   {gene_key}: no mapping, skip")
            for mf in mut_files:
                rejected.append({"gene": gene, "file": mf.name, "reason": "no_mapping_for_gene"})
            continue

        # Precompute mutant->transcript map once
        mp = dict(zip(mapping_df["mutant"].astype(str), mapping_df["transcript"].astype(str)))
        # {sha -> verbatim token} for this gene. _per_gene_mut_tasks minted the
        # filename from gene_key + the same token, so this inverts it exactly.
        sha_to_token = {}
        for _t in mp:
            sha_to_token[mint_pkey(gene_key, _t).rsplit("-", 1)[-1]] = _t

        # The projection is built over WT coordinates, so it must span the furthest
        # WT site. Computed once per gene; the variant-dependent part is per pkey.
        wt_span = max((h["site_pos"] for h in wt_hits if h["site_pos"] is not None),
                      default=0)

        print(f"[PARSE]   {gene_key}: WT hits={len(wt_hits)}, MUT files={len(mut_files)}")
        for mf in mut_files:
            parsed_files += 1

            # GENE-<mut>-mut. `core` is guaranteed to contain a hyphen: that is
            # the predicate _collect_gene_outputs bucketed this file on, and it
            # names the files it rejects itself. The local re-test that used to
            # stand here computed the identical `core` from the identical name,
            # so it could never fire; it has been removed rather than left
            # looking like live coverage.
            base = mf.name[:-len("-miranda.out")]
            core = base[:-4]  # GENE-<mut>
            # Use utility extractor (rsplit-based, hyphen-safe). It cannot raise
            # for a hyphen-containing str (guaranteed by the bucketing filter), so
            # the former try/except fallback (`core.split('-',1)[1]`, first-hyphen
            # unsafe) was dead code and has been removed.
            # `core` is {gene_key}-{sha}. The sha is not the token, so resolve it
            # through the per-gene reverse map built from this gene's own mapping;
            # mp, should_skip_mutation and the `rejected` rows all need the
            # VERBATIM spelling.
            gene_tok, mut_tok = extract_mutation_from_sequence_name(core)
            mutation_token = sha_to_token.get(mut_tok, mut_tok)

            transcript_nt = mp.get(mutation_token, None)
            if not transcript_nt:
                miss_map += 1
                rejected.append({"gene": gene, "mutant": mutation_token,
                                 "file": mf.name, "reason": "no_mapping_for_mutant"})
                continue

            if should_skip_mutation(gene, mutation_token, failure_map):
                rejected.append({"gene": gene, "mutant": mutation_token,
                                 "file": mf.name, "reason": "validation_log_filtered"})
                continue

            # parse_variant supersedes get_mutation_data here for the same reason
            # it does in _per_gene_mut_tasks: get_mutation_data raised ValueError
            # for a valid indel and for a corrupt token alike, and the except arm
            # set tx_pos_1 = None, which silently voided distance_to_snv for EVERY
            # row of that pkey rather than only for genuinely bad tokens.
            # Same substrate split as _per_gene_mut_tasks, and for the same
            # reason: the intron/pre-mRNA mapping carries decomposed PIECES, and
            # a wholly deleted piece is 'REF<pos>del', which parse_variant
            # refuses because Variant cannot hold an empty allele. Patching only
            # the build path left these tokens building, running through miranda,
            # and then being discarded HERE -- measured: 4 outputs rejected as
            # unparseable_token (GTATTG1del, GTATTG33931del, GTAAGG1del,
            # GTAAGG17467del), i.e. the splice-site deletions had their miranda
            # run computed and thrown away.
            if _substrate_of(gene_key) == "transcript":
                variant = parse_variant(str(transcript_nt), is_nt=True)
                fields = (None if variant is None else
                          (variant.pos, variant.ref, variant.alt,
                           variant.kind, variant.length_delta))
            else:
                fields = piece_fields(str(transcript_nt))
            if fields is None:
                tx_pos_1 = None
                ref_len = alt_len = 1
                variant_kind = ""
                length_delta = None
                wt_to_mut = None
                rejected.append({"gene": gene, "mutant": mutation_token,
                                 "file": mf.name, "reason": "unparseable_token",
                                 "transcript_token": str(transcript_nt)})
            else:
                tx_pos_1, v_ref, v_alt, variant_kind, length_delta = fields
                ref_len, alt_len = len(v_ref), len(v_alt)
                # The WT->MUT projection. Identity when ref_len == alt_len.
                wt_to_mut = align_wt_to_mut(wt_span, tx_pos_1 - 1, ref_len, alt_len)

            pkey = mint_pkey(gene, mutation_token)

            # WT rows for this pkey
            for h in wt_hits:
                site_pos = h["site_pos"]
                # Distance is measured to the REF SPAN in the WT frame.
                dist = distance_to_variant(site_pos, tx_pos_1, ref_len)
                join_pos, align_status = site_pos, "aligned"
                if wt_to_mut is None:
                    # No variant record -> no projection exists. Say so rather than
                    # letting the raw position masquerade as a projected one.
                    align_status = "unprojected"
                elif site_pos is not None and 0 <= site_pos - 1 < len(wt_to_mut):
                    j = wt_to_mut[site_pos - 1]
                    if j is None:
                        # The edit deleted this anchor: no mutant counterpart. The
                        # deleted span collapses to the edit start, so that is the
                        # only defensible join coordinate -- and the row is marked
                        # so downstream can refuse to read a shift out of it.
                        join_pos, align_status = tx_pos_1, "deleted"
                    else:
                        join_pos = j + 1
                records.append({
                    "pkey": pkey,
                    "allele": "WT",
                    "mirna_id": h["mirna_id"],
                    "site_pos": site_pos,
                    "tot_score": h["tot_score"],
                    "tot_energy": h["tot_energy"],
                    "max_score": h["max_score"],
                    "max_energy": h["max_energy"],
                    "strand": h["strand"],
                    "len_mirna": h["len_mirna"],
                    "len_target": h["len_target"],
                    "visibility_flag": 1 if (h["tot_score"] is not None and h["tot_score"] >= VISIBILITY_THRESHOLD) else 0,
                    "distance_to_snv": dist,
                    "locus_id": None,
                    "segment_id": None,
                    "parser_confidence": h["parser_confidence"],
                    "run_meta": f"{gene_key}-wt",
                    "substrate": substrate,
                    "join_pos": join_pos,
                    "align_status": align_status,
                    "variant_kind": variant_kind,
                    "length_delta": length_delta,
                })

            # MUT rows
            mut_text = mf.read_text(encoding="utf-8", errors="ignore")
            mut_hits = parse_miranda_text(mut_text)

            # Delete MUT .out file after parsing
            try:
                mf.unlink(missing_ok=True)
            except Exception:
                pass

            if not mut_hits:
                miss_hits += 1
                rejected.append({"gene": gene, "mutant": mutation_token,
                                 "file": mf.name, "reason": "no_hits_parsed"})
                continue

            # Mutant-frame indices with no WT counterpart, i.e. bases the ALT
            # inserted. WT index offset+k maps to MUT index offset+k only while
            # k < min(ref_len, alt_len); everything the ALT adds beyond that is new
            # sequence. Empty whenever alt_len <= ref_len.
            # None, not 0: with no variant record there is no inserted region at
            # all, and 0 would read as a real coordinate at the 5' end.
            ins_lo = ins_hi = None
            if variant is not None:
                ins_lo = variant.pos0 + min(ref_len, alt_len)
                ins_hi = variant.pos0 + alt_len

            for h in mut_hits:
                site_pos = h["site_pos"]
                # MUT sites are in the mutant frame, where the variant's span is
                # the ALT allele: the two frames agree 5' of the edit, so the span
                # still starts at tx_pos_1, but its length is len(ALT) here.
                dist = distance_to_variant(site_pos, tx_pos_1, alt_len)
                if variant is None:
                    align_status = "unprojected"
                elif site_pos is not None and ins_lo <= site_pos - 1 < ins_hi:
                    align_status = "inserted"
                else:
                    align_status = "aligned"
                records.append({
                    "pkey": pkey,
                    "allele": "MUT",
                    "mirna_id": h["mirna_id"],
                    "site_pos": site_pos,
                    "tot_score": h["tot_score"],
                    "tot_energy": h["tot_energy"],
                    "max_score": h["max_score"],
                    "max_energy": h["max_energy"],
                    "strand": h["strand"],
                    "len_mirna": h["len_mirna"],
                    "len_target": h["len_target"],
                    "visibility_flag": 1 if (h["tot_score"] is not None and h["tot_score"] >= VISIBILITY_THRESHOLD) else 0,
                    "distance_to_snv": dist,
                    "locus_id": None,
                    "segment_id": None,
                    "parser_confidence": h["parser_confidence"],
                    "run_meta": f"{gene_key}-mut",
                    "substrate": substrate,
                    # A MUT row is already in the frame the join is made in.
                    "join_pos": site_pos,
                    "align_status": align_status,
                    "variant_kind": variant_kind,
                    "length_delta": length_delta,
                })

        # Delete WT .out file after all MUT files for this gene are parsed
        try:
            wt_file.unlink(missing_ok=True)
        except Exception:
            pass

        if gi % 25 == 0:
            print(f"[PARSE]   Progress: {gi}/{len(genes)} genes processed")

    if miss_map or miss_hits:
        print(f"[PARSE] Mapping misses: {miss_map}")
        print(f"[PARSE] Files with no parsed hits: {miss_hits}")

    sites_df = pd.DataFrame.from_records(records)
    if sites_df.empty:
        print("[PARSE] Sites table built with 0 rows")
    else:
        print(f"[PARSE] Sites table built with {len(sites_df)} rows")
    if rejected:
        print(f"[PARSE] Rejected {len(rejected)} MUT output(s); reasons in miranda.run_summary.json")
    return sites_df, rejected

# Both chaining passes below run on `join_pos`, never on `site_pos`. The two are
# equal for every length-preserving variant, so SNV and MNV behaviour is
# unchanged; under an indel `join_pos` puts both alleles of one physical site in
# the same coordinate system, which is what stops a renumbered site from being
# chained into two separate loci and reported as a lost/gained pair.
def assign_loci(sites_df: pd.DataFrame) -> pd.DataFrame:
    print("[PARSE] Assigning per-miRNA loci...")
    if sites_df is None or sites_df.empty or "pkey" not in sites_df.columns:
        print("[PARSE] No sites to assign loci")
        return sites_df
    updates = 0
    # Grouped by (pkey, substrate), never pkey alone. One variant now has rows on
    # several molecules -- the excised intron and the pre-mRNA -- and chaining
    # across them would merge sites from two different sequences into one locus.
    gcols = ["pkey", "substrate"] if "substrate" in sites_df.columns else ["pkey"]
    for _gk, sub in sites_df.groupby(gcols, dropna=False):
        for mirna_id, g in sub.groupby("mirna_id", dropna=False):
            g = g.sort_values("join_pos")
            current = 1
            last_pos = None
            last_unpaired = False
            for idx, row in g.iterrows():
                pos = row["join_pos"]
                # A WT anchor the edit DELETED has no mutant counterpart at all --
                # align_wt_to_mut returns None for that index precisely to say so.
                # Its join_pos was collapsed to the edit start only so the row
                # still sorts; chaining on that coordinate would marry it to
                # whatever the mutant happens to carry at the junction, which is a
                # different physical site, and would merge every anchor the same
                # deletion swallowed into one locus -- where head(1) then discards
                # all but the best-scoring one. So it breaks the chain on BOTH
                # sides and occupies a locus of its own.
                # No row is ever 'deleted' under a length-preserving variant, so
                # both predicates are constant-False on the SNV/MNV path.
                unpaired = (row.get("align_status") == "deleted")
                if (last_pos is None or unpaired or last_unpaired
                        or (pos - last_pos) > MERGE_WINDOW_NT):
                    locus_id = f"m{current}"
                    current += 1
                else:
                    locus_id = f"m{current-1}"
                last_pos = pos
                last_unpaired = unpaired
                sites_df.at[idx, "locus_id"] = locus_id
                updates += 1
    print(f"[PARSE] Loci assigned for {updates} rows")
    return sites_df

def assign_segments(sites_df: pd.DataFrame) -> pd.DataFrame:
    print("[PARSE] Assigning cross-miRNA segments...")
    if sites_df is None or sites_df.empty or "pkey" not in sites_df.columns:
        print("[PARSE] No sites to assign segments")
        return sites_df
    updates = 0
    gcols = ["pkey", "substrate"] if "substrate" in sites_df.columns else ["pkey"]
    for _gk, sub in sites_df.groupby(gcols, dropna=False):
        g = sub.sort_values("join_pos")
        current = 1
        last_pos = None
        for idx, row in g.iterrows():
            pos = row["join_pos"]
            if last_pos is None or (pos - last_pos) > SEGMENT_WINDOW_NT:
                seg_id = f"s{current}"
                current += 1
            else:
                seg_id = f"s{current-1}"
            last_pos = pos
            sites_df.at[idx, "segment_id"] = seg_id
            updates += 1
    print(f"[PARSE] Segments assigned for {updates} rows")
    return sites_df

def build_events(sites_df: pd.DataFrame) -> pd.DataFrame:
    """Join the two alleles per (pkey, mirna_id, locus_id) and score the change.

    `dpos` is the site's displacement AFTER the renumbering the variant forces on
    every downstream coordinate has been removed: it is measured against the WT
    anchor's projection into the mutant frame (`join_pos_wt`), not against the raw
    WT position. Without that, a 20 nt insertion made every site 3' of it report
    dpos = 20 and classify as `shifted` although its sequence was never touched.
    For a length-preserving variant join_pos == site_pos, so dpos is unchanged.

    Three new columns: `variant_kind`, `length_delta` and `qc_flags`. `qc_flags`
    names anything the row could not measure, so an absent value is never
    indistinguishable from a measured one.
    """
    print("[PARSE] Building events...")
    cols = [
        "pkey","mirna_id","locus_id","segment_id",
        "wt_pos","mut_pos","dpos",
        "wt_tot_score","mut_tot_score","delta_tot_score","pct_delta",
        "wt_tot_energy","mut_tot_energy","delta_energy",
        "distance_to_snv",
        "rank_wt","rank_mut",
        "conf_wt","conf_mut","conf_weighted_delta",
        "cls","is_high_impact","priority","in_radius",
        "variant_kind","length_delta","qc_flags","substrate"
    ]
    # join_pos/align_status are required, not optional: without them dpos would
    # fall back to comparing two frames that an indel has pulled apart, and the
    # error would be silent numbers rather than a missing table.
    required = {"pkey", "mirna_id", "locus_id", "allele", "join_pos", "align_status"}
    if sites_df is None or sites_df.empty or not required.issubset(sites_df.columns):
        print("[PARSE] No sites to build events")
        return pd.DataFrame(columns=cols)

    # substrate joins the key for the same reason it joins the locus grouping: a
    # WT site on the intron record must never be paired with a MUT site on the
    # pre-mRNA record.
    key = (["pkey","substrate","mirna_id","locus_id"] if "substrate" in sites_df.columns
           else ["pkey","mirna_id","locus_id"])
    wt_top = (sites_df[sites_df["allele"]=="WT"]
              .sort_values("tot_score", ascending=False)
              .groupby(key, dropna=False).head(1))
    mut_top = (sites_df[sites_df["allele"]=="MUT"]
               .sort_values("tot_score", ascending=False)
               .groupby(key, dropna=False).head(1))
    merged = pd.merge(wt_top, mut_top, on=key, how="outer", suffixes=("_wt","_mut"))

    out = []
    for _, r in merged.iterrows():
        pkey = r["pkey"]; mirna = r["mirna_id"]; locus = r["locus_id"]
        substrate = r.get("substrate", "transcript")
        seg = r.get("segment_id_mut") if pd.notnull(r.get("segment_id_mut")) else r.get("segment_id_wt")
        flags: List[str] = []

        wt_pos  = int(r["site_pos_wt"]) if pd.notnull(r.get("site_pos_wt")) else None
        mut_pos = int(r["site_pos_mut"]) if pd.notnull(r.get("site_pos_mut")) else None

        # The WT anchor expressed in the mutant frame. This is what mut_pos is
        # comparable to; wt_pos is not, once the variant changes length.
        wt_join = int(r["join_pos_wt"]) if pd.notnull(r.get("join_pos_wt")) else None
        wt_align = r.get("align_status_wt")
        mut_align = r.get("align_status_mut")

        if wt_align == "deleted":
            # The variant removed this anchor outright, so it has no mutant
            # coordinate and no displacement. A number here would be a distance
            # from the point the deletion collapsed to, not a shift of the site.
            dpos = None
            flags.append("WT_ANCHOR_DELETED")
        elif wt_join is not None and mut_pos is not None:
            dpos = mut_pos - wt_join
        else:
            dpos = None
        if mut_align == "inserted":
            flags.append("MUT_SITE_IN_INSERTED_SEQUENCE")
        if wt_align == "unprojected" or mut_align == "unprojected":
            flags.append("NO_VARIANT_RECORD:join_on_raw_position")

        # An absent hit and an unreadable score are NOT the same thing.
        #   absent  -- the allele was scanned and MirandA reported nothing here,
        #              i.e. the site did not clear its reporting threshold. That is
        #              an observation, and the pipeline's existing model for it is
        #              0.0; it is kept, but it is now NAMED so a consumer can tell
        #              it from a measured zero.
        #   unread  -- the hit row exists but parse_miranda_text could not read its
        #              score field. Nothing was measured, so the score and every
        #              quantity derived from it stay EMPTY. Coalescing this to 0.0
        #              fabricated a full-magnitude event out of a parser failure.
        def _score(side):
            present = pd.notnull(r.get(f"allele_{side}"))
            raw = r.get(f"tot_score_{side}")
            if not present:
                flags.append(f"{side.upper()}_HIT_ABSENT")
                return 0.0
            if pd.isnull(raw):
                flags.append(f"{side.upper()}_SCORE_UNPARSED")
                return None
            return float(raw)

        wt_score = _score("wt")
        mut_score = _score("mut")
        measurable = wt_score is not None and mut_score is not None

        if measurable:
            delta = mut_score - wt_score
            denom = max(wt_score, mut_score) if (wt_score or mut_score) else 0.0
            pct_delta = (delta / denom) if denom else 0.0
        else:
            delta = None
            pct_delta = None

        wt_energy = float(r["tot_energy_wt"]) if pd.notnull(r.get("tot_energy_wt")) else None
        mut_energy = float(r["tot_energy_mut"]) if pd.notnull(r.get("tot_energy_mut")) else None
        delta_energy = (mut_energy - wt_energy) if (mut_energy is not None and wt_energy is not None) else None

        dist = (float(r["distance_to_snv_mut"]) if pd.notnull(r.get("distance_to_snv_mut"))
                else (float(r["distance_to_snv_wt"]) if pd.notnull(r.get("distance_to_snv_wt")) else None))

        wt_vis = (wt_score is not None and wt_score >= VISIBILITY_THRESHOLD)
        mut_vis = (mut_score is not None and mut_score >= VISIBILITY_THRESHOLD)

        if not measurable:
            # Classifying against a score that was never read would assert a
            # direction the data does not support.
            cls = "undetermined"
        elif not wt_vis and mut_vis:
            cls = "gained"
        elif wt_vis and not mut_vis:
            cls = "lost"
        elif wt_vis and mut_vis and dpos is not None and abs(dpos) >= SHIFT_NT:
            cls = "shifted"
        elif wt_vis and mut_vis and abs(delta) >= 1.0 and (dpos is None or abs(dpos) < SHIFT_NT):
            cls = "strengthened" if delta > 0 else "weakened"
        else:
            cls = "none"

        in_radius = 1 if (dist is not None and dist <= REPORT_RADIUS) else 0

        class_bonus = {"gained":3.0,"lost":3.0,"shifted":1.0,"strengthened":0.5,
                       "weakened":0.5,"none":0.0,"undetermined":0.0}[cls]
        if measurable:
            is_high = 1 if (abs(delta) >= HIGH_CUTOFF) or cls in ("gained","lost") else 0
            priority = abs(delta) * exp_weight(dist, DISTANCE_K) + class_bonus
            conf_weighted = (1.0 if mut_vis else 0.5)*mut_score - (1.0 if wt_vis else 0.5)*wt_score
        else:
            is_high = None
            priority = None
            conf_weighted = None

        out.append({
            "pkey": pkey,
            "mirna_id": mirna,
            "locus_id": locus,
            "segment_id": seg,
            "wt_pos": wt_pos,
            "mut_pos": mut_pos,
            "dpos": dpos,
            "wt_tot_score": wt_score,
            "mut_tot_score": mut_score,
            "delta_tot_score": delta,
            "pct_delta": pct_delta,
            "wt_tot_energy": wt_energy,
            "mut_tot_energy": mut_energy,
            "delta_energy": delta_energy,
            "distance_to_snv": dist,
            "rank_wt": None,
            "rank_mut": None,
            # 0.5 is this pipeline's code for "measured, below the visibility
            # threshold". On a score that was never read that asserts an
            # observation which did not happen -- the same fabrication as the 0.0
            # removed from the score itself. Empty, like every other quantity
            # derived from an unread score.
            "conf_wt": (1.0 if wt_vis else 0.5) if wt_score is not None else None,
            "conf_mut": (1.0 if mut_vis else 0.5) if mut_score is not None else None,
            "conf_weighted_delta": conf_weighted,
            "cls": cls,
            "is_high_impact": is_high,
            "priority": priority,
            "in_radius": in_radius,
            "variant_kind": (r.get("variant_kind_wt") if pd.notnull(r.get("variant_kind_wt"))
                             else r.get("variant_kind_mut")),
            "length_delta": (r.get("length_delta_wt") if pd.notnull(r.get("length_delta_wt"))
                             else r.get("length_delta_mut")),
            "qc_flags": ";".join(flags),
            "substrate": substrate,
        })

    events_df = pd.DataFrame(out, columns=cols)
    if events_df.empty:
        print("[PARSE] Events table built with 0 rows")
        return events_df

    # Force the numeric columns to a numeric dtype so an EMPTY value is NaN rather
    # than a Python None. A pkey whose every score was unreadable leaves an
    # all-None column, which pandas infers as object dtype -- and on object dtype
    # `.abs()` raises TypeError and `>= HIGH_CUTOFF` raises on the None instead of
    # evaluating False. Both crash build_summary. Where the values are already
    # numeric (every SNV row) this is a no-op and the dtype is unchanged.
    # rank_wt/rank_mut are deliberately excluded: they are None here and are filled
    # with ints below, so coercing them now would blank them.
    for c in ("wt_pos", "mut_pos", "dpos",
              "wt_tot_score", "mut_tot_score", "delta_tot_score", "pct_delta",
              "wt_tot_energy", "mut_tot_energy", "delta_energy",
              "distance_to_snv", "conf_wt", "conf_mut", "conf_weighted_delta",
              "is_high_impact", "priority", "length_delta"):
        events_df[c] = pd.to_numeric(events_df[c], errors="coerce")

    # ranks within pkey
    for pkey, sub in events_df.groupby("pkey", dropna=False):
        ms = sub.sort_values("mut_tot_score", ascending=False)
        for i, (idx, _) in enumerate(ms.iterrows(), 1):
            events_df.at[idx, "rank_mut"] = i
        ws = sub.sort_values("wt_tot_score", ascending=False)
        for i, (idx, _) in enumerate(ws.iterrows(), 1):
            events_df.at[idx, "rank_wt"] = i

    print(f"[PARSE] Events table built with {len(events_df)} rows")
    return events_df

def build_summary(events_df: pd.DataFrame, sites_df: pd.DataFrame) -> pd.DataFrame:
    print("[PARSE] Building summary...")
    if events_df is None or events_df.empty:
        print("[PARSE] No events; writing empty summary")
        return pd.DataFrame(columns=[
            "pkey", "n_hits_wt", "n_hits_mut", "n_mirna", "n_loci", "n_segments",
            "n_competitive_segments", "n_new_competitors",
            "global_count_gained_high", "global_count_lost_high", "global_count_shifted",
            "global_max_abs_delta_score", "global_sum_weighted_abs_delta",
            "nearest_event_bp_any",
            "local_count_gained_high", "local_count_lost_high", "local_count_shifted",
            "local_max_abs_delta_score", "nearest_event_bp_local", "frac_effect_in_radius",
            "max_segment_abs_delta_best", "frac_effect_in_competitive_segments",
            "top_event_mirna", "top_event_class", "top_event_delta_score", "top_event_pos", "qc_flags",
            "variant_kind", "length_delta", "n_sites_aligned", "n_sites_deleted",
            "n_sites_inserted", "variant_qc", "substrate"
        ])

    summary_records: List[Dict] = []

    mut_vis_counts = {}
    wt_vis_counts = {}
    if sites_df is not None and not sites_df.empty:
        mut_vis = sites_df[(sites_df["allele"]=="MUT") & (sites_df["visibility_flag"]==1)].groupby("pkey", dropna=False).size()
        wt_vis  = sites_df[(sites_df["allele"]=="WT") & (sites_df["visibility_flag"]==1)].groupby("pkey", dropna=False).size()
        mut_vis_counts = mut_vis.to_dict()
        wt_vis_counts = wt_vis.to_dict()

    has_sub = "substrate" in events_df.columns
    for _gk, ev in events_df.groupby(["pkey", "substrate"] if has_sub else ["pkey"],
                                     dropna=False):
        pkey, substrate = (_gk if has_sub else (_gk, "transcript"))
        if sites_df is not None and not sites_df.empty:
            sites_sub = sites_df[sites_df["pkey"] == pkey]
            if has_sub and "substrate" in sites_df.columns:
                sites_sub = sites_sub[sites_sub["substrate"] == substrate]
        else:
            sites_sub = pd.DataFrame()

        n_hits_wt = int((sites_sub["allele"].eq("WT") & sites_sub["visibility_flag"].eq(1)).sum()) if not sites_sub.empty else 0
        n_hits_mut = int((sites_sub["allele"].eq("MUT") & sites_sub["visibility_flag"].eq(1)).sum()) if not sites_sub.empty else 0
        n_mirna = ev["mirna_id"].nunique()
        n_loci = len(ev)
        segments = [] if sites_sub.empty else sites_sub["segment_id"].dropna().unique().tolist()
        n_segments = len(segments)

        n_competitive_segments = 0
        n_new_competitors = 0
        if not sites_sub.empty:
            for seg in segments:
                seg_sites = sites_sub[(sites_sub["segment_id"] == seg) & (sites_sub["visibility_flag"] == 1)]
                comp_mut = seg_sites[seg_sites["allele"]=="MUT"]["mirna_id"].nunique()
                comp_wt  = seg_sites[seg_sites["allele"]=="WT"]["mirna_id"].nunique()
                if (comp_mut + comp_wt) >= 2:
                    n_competitive_segments += 1
                if comp_mut > comp_wt:
                    n_new_competitors += 1

        # `ev` comes from groupby, which never yields an empty group, so the
        # former `if len(ev): ... else: <zeros>` arm was unreachable and has been
        # removed rather than left looking like a live fallback.
        global_max_abs_delta_score = ev["delta_tot_score"].abs().max()
        weighted_abs = 0.0
        total_abs = 0.0
        n_undefined_delta = 0
        for _, r in ev.iterrows():
            d_raw = r["delta_tot_score"]
            if pd.isnull(d_raw):
                # An event whose delta could not be measured contributes nothing.
                # Treating it as 0.0 would quietly deflate both aggregates below
                # AND the frac_effect_in_radius denominator.
                n_undefined_delta += 1
                continue
            absd = abs(d_raw)
            total_abs += absd
            weighted_abs += absd * exp_weight(r["distance_to_snv"], DISTANCE_K)
        global_sum_weighted_abs_delta = weighted_abs

        nearest_event_bp_any = ev["distance_to_snv"].min() if ev["distance_to_snv"].notnull().any() else None

        ev_local = ev[(ev["distance_to_snv"].notnull()) & (ev["distance_to_snv"] <= REPORT_RADIUS)]
        local_count_gained_high = len(ev_local[(ev_local["cls"] == "gained") & (ev_local["mut_tot_score"] >= HIGH_CUTOFF)])
        local_count_lost_high = len(ev_local[(ev_local["cls"] == "lost") & (ev_local["wt_tot_score"] >= HIGH_CUTOFF)])
        local_count_shifted = len(ev_local[(ev_local["cls"] == "shifted")])

        if len(ev_local):
            local_max_abs_delta_score = ev_local["delta_tot_score"].abs().max()
            nearest_event_bp_local = ev_local["distance_to_snv"].min()
            local_abs = ev_local["delta_tot_score"].abs().sum()
        else:
            local_max_abs_delta_score = 0.0
            nearest_event_bp_local = None
            local_abs = 0.0

        frac_effect_in_radius = safe_div(local_abs, total_abs)

        max_segment_abs_delta_best = 0.0
        for seg in segments:
            ev_seg = ev[ev["segment_id"] == seg]
            if len(ev_seg):
                seg_delta = ev_seg["delta_tot_score"].abs().max()
                if seg_delta > max_segment_abs_delta_best:
                    max_segment_abs_delta_best = seg_delta

        # Same reasoning as above: groupby cannot hand back an empty `ev`, so the
        # zero-filled else arm this replaces was unreachable.
        top = ev.sort_values("priority", ascending=False).iloc[0]
        top_event_mirna = top["mirna_id"]
        top_event_class = top["cls"]
        top_event_delta_score = top["delta_tot_score"]
        top_event_pos = top["mut_pos"] if pd.notnull(top.get("mut_pos")) else top.get("wt_pos")

        qc_flags = []
        if (n_hits_wt + n_hits_mut) == 0:
            qc_flags.append("no_hits")
        if nearest_event_bp_any is not None and nearest_event_bp_any > 2000:
            qc_flags.append("far_event>2kb")
        if global_max_abs_delta_score == 0.0:
            qc_flags.append("low_signal_only")

        # Variant-class provenance. Kept OUT of qc_flags: that column's existing
        # vocabulary is consumed downstream, and codes like NO_LOCAL_EVENTS fire on
        # SNV rows too, so appending them there would silently rewrite output for
        # mutations this work is not supposed to touch.
        variant_kind = ""
        length_delta = None
        if "variant_kind" in ev.columns and ev["variant_kind"].notnull().any():
            variant_kind = ev["variant_kind"].dropna().iloc[0]
        if "length_delta" in ev.columns and ev["length_delta"].notnull().any():
            length_delta = ev["length_delta"].dropna().iloc[0]

        # Site-level alignment accounting over the UNION of both alleles. Counting
        # WT sites alone would report a 20 nt insertion as fully aligned, because
        # every WT site does keep a counterpart -- the sites with no counterpart
        # are all on the mutant side.
        n_sites_deleted = n_sites_inserted = 0
        n_sites_wt = 0
        n_sites_unprojected = 0
        if not sites_sub.empty and "align_status" in sites_sub.columns:
            wt_sub = sites_sub[sites_sub["allele"] == "WT"]
            n_sites_wt = len(wt_sub)
            n_sites_deleted = int((wt_sub["align_status"] == "deleted").sum())
            n_sites_inserted = int((sites_sub[sites_sub["allele"] == "MUT"]["align_status"]
                                    == "inserted").sum())
            n_sites_unprojected = int((sites_sub["align_status"] == "unprojected").sum())
        n_sites_aligned = n_sites_wt - n_sites_deleted

        variant_qc = []
        # Emitted FIRST because it is the reason variant_kind, length_delta and
        # every distance on this row are empty. The events table and
        # run_summary.json already named it; the per-mutation summary is the table
        # a consumer actually reads, and it showed four empty columns and no cause.
        if n_sites_unprojected:
            variant_qc.append("NO_VARIANT_RECORD:token_unparsed")
        if length_delta:
            variant_qc.append(
                f"length_changed:{int(length_delta):+d}nt;"
                f"aligned_{n_sites_aligned}/{n_sites_wt + n_sites_inserted};"
                f"deleted_{n_sites_deleted};inserted_{n_sites_inserted}")
        if n_undefined_delta:
            variant_qc.append(f"delta_undefined_events:{n_undefined_delta}")
        if not len(ev_local):
            # Documents that local_max_abs_delta_score below is a max over an
            # EMPTY set, not a measured zero.
            variant_qc.append("NO_LOCAL_EVENTS")
        if n_sites_deleted:
            variant_qc.append(f"wt_sites_without_mutant_coordinate:{n_sites_deleted}")
        if n_sites_inserted:
            variant_qc.append(f"mut_sites_in_inserted_sequence:{n_sites_inserted}")

        summary_records.append({
            "pkey": pkey,
            "n_hits_wt": n_hits_wt,
            "n_hits_mut": n_hits_mut,
            "n_mirna": n_mirna,
            "n_loci": n_loci,
            "n_segments": n_segments,
            "n_competitive_segments": n_competitive_segments,
            "n_new_competitors": n_new_competitors,
            "global_count_gained_high": len(ev[(ev["cls"] == "gained") & (ev["mut_tot_score"] >= HIGH_CUTOFF)]),
            "global_count_lost_high": len(ev[(ev["cls"] == "lost") & (ev["wt_tot_score"] >= HIGH_CUTOFF)]),
            "global_count_shifted": len(ev[(ev["cls"] == "shifted")]),
            "global_max_abs_delta_score": global_max_abs_delta_score,
            "global_sum_weighted_abs_delta": global_sum_weighted_abs_delta,
            "nearest_event_bp_any": nearest_event_bp_any,
            "local_count_gained_high": local_count_gained_high,
            "local_count_lost_high": local_count_lost_high,
            "local_count_shifted": local_count_shifted,
            "local_max_abs_delta_score": local_max_abs_delta_score,
            "nearest_event_bp_local": nearest_event_bp_local,
            "frac_effect_in_radius": frac_effect_in_radius,
            "max_segment_abs_delta_best": max_segment_abs_delta_best,
            "frac_effect_in_competitive_segments": safe_div(n_competitive_segments, max(n_segments, 1)),
            "top_event_mirna": top_event_mirna,
            "top_event_class": top_event_class,
            "top_event_delta_score": top_event_delta_score,
            "top_event_pos": top_event_pos,
            "qc_flags": ";".join(qc_flags) if qc_flags else "",
            "variant_kind": variant_kind,
            "length_delta": length_delta,
            "n_sites_aligned": n_sites_aligned,
            "n_sites_deleted": n_sites_deleted,
            "n_sites_inserted": n_sites_inserted,
            "variant_qc": ";".join(variant_qc),
            "substrate": substrate,
        })

    summary_df = pd.DataFrame.from_records(summary_records)
    print(f"[PARSE] Summary table built with {len(summary_df)} rows")
    return summary_df

# -------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------
def _load_substrate_sequences(input_path: str) -> Dict[str, str]:
    """{gene_key: sequence} for every NON-transcript substrate in the input FASTAs.

    Reads the pre_mRNA record and each intron<N>|... record that
    variant_mapping emits. Both are already in TRANSCRIPT orientation, which
    is the only orientation miranda can use: it scans the strand it is handed and
    does not try the reverse complement. Measured on this machine -- a perfect
    let-7a site scores 200.00 on the forward sequence and returns ZERO hits on its
    reverse complement, so a genomic-forward minus-strand record would report "no
    miRNA interaction" for a gene that has plenty.
    """
    path = Path(input_path)
    files = ([path] if path.is_file()
             else sorted(f for f in path.iterdir()
                         if f.suffix.lower() in (".fasta", ".fa", ".fna", ".faa")))
    out: Dict[str, str] = {}
    for f in files:
        gene = (extract_gene_from_filename(str(f)) or f.stem).upper()
        try:
            recs = read_fasta(str(f))
        except Exception:
            continue
        for name, seq in recs.items():
            head = name.split("|", 1)[0]
            if head == "pre_mRNA" or head.startswith("intron"):
                out[_gene_key(gene, head)] = seq
    return out


def main():
    parser = argparse.ArgumentParser(description="Unified MirandA pipeline (WT+MUT) with comparative outputs")
    # IO
    parser.add_argument("-i", "--input", required=True, help="WT FASTA file, or directory of FASTA files")
    parser.add_argument("-o", "--output", required=True, help="output directory of miranda")
    # MirandA
    parser.add_argument("-m", "--miranda_dir", default=None,
                        help="path to directory containing the miranda executable; "
                             "omit to use miranda from PATH (e.g. conda install -c bioconda miranda)")
    parser.add_argument("-d", "--mirna_db", required=True, help="path to mirna database")
    # Mapping
    parser.add_argument("--mapping-dir", required=True, help="transcript mapping CSV/TSV file, or directory of CSV/TSV files")
    parser.add_argument("--intron-premrna-mapping",
                        help="intron_premRNA_mapping_<GENE>.csv (file or directory). Intronic "
                             "variants are scanned against BOTH the intron record and the "
                             "pre_mRNA record.")
    parser.add_argument("--strict-introns", action="store_true",
                        help="Pass miranda -strict (strict 5' seed pairing) for the intron and "
                             "pre_mRNA substrates only. WT and MUT always use the same setting.")
    # Performance
    parser.add_argument("--no-parallel", action='store_true', help="disable parallel processing")
    parser.add_argument("--max-workers", type=int, help="maximum number of parallel workers")
    parser.add_argument("--log", help="Validation log file or directory to skip failed mutations")
    # WT header
    parser.add_argument("--wt-header", default="transcript", help="Preferred WT FASTA header")
    args = parser.parse_args()

    # Validate miranda availability
    if args.miranda_dir:
        bin_path = os.path.join(args.miranda_dir, "miranda")
        if not (os.path.isfile(bin_path) and os.access(bin_path, os.X_OK)):
            parser.error(f"miranda executable not found or not executable at {bin_path}")
    else:
        import shutil
        if not shutil.which("miranda"):
            parser.error("miranda not found on PATH; install via: conda install -c bioconda miranda")

    failure_map = load_validation_failures(args.log) if args.log else {}

    mapping_dict = load_transcript_mappings(args.mapping_dir)
    if not mapping_dict:
        raise ValueError(f"No valid mapping files found in {args.mapping_dir}")

    wt_sequences = load_wt_sequences(args.input, wt_header=args.wt_header)
    if not wt_sequences:
        raise ValueError(f"No WT sequences found in {args.input}")

    # ---- intronic substrates -------------------------------------------
    # An intronic variant has no transcript coordinate, so it can never appear in
    # the transcript mapping -- it is absent from this pipeline entirely rather
    # than rejected. It arrives through its own mappings and is scanned against
    # the intron record and the pre-mRNA record, each keyed GENE~<substrate>.
    strict_substrates: set = set()
    substrate_seqs = _load_substrate_sequences(args.input)
    if args.intron_premrna_mapping:
        # One file per gene, columns: mutant, orf, intron, pre-mRNA-Transcript.
        # `orf` and `intron` are mutually exclusive per row, so a variant spanning
        # a splice site occupies TWO rows under one mutant. Exonic rows are
        # skipped here -- they are already covered on the transcript substrate,
        # and scanning them again would double every exonic result.
        ipm_root = Path(args.intron_premrna_mapping)
        files = ([ipm_root] if ipm_root.is_file()
                 else sorted(f for f in ipm_root.rglob("*")
                             if f.is_file() and f.suffix.lower() in (".csv", ".tsv", ".txt")))
        per_gene: Dict[str, Dict[str, List[Tuple[str, str]]]] = {}
        for f in files:
            gene = (extract_gene_from_filename(str(f)) or f.stem).upper()
            try:
                with open(f, newline="") as fh:
                    for row in csv.DictReader(fh):
                        mut = (row.get("mutant") or "").strip()
                        orf = (row.get("orf") or "").strip()
                        intron = (row.get("intron") or "").strip()
                        pre = (row.get("pre-mRNA-Transcript") or "").strip()
                        if not mut:
                            continue
                        # split_piece_cell is the shared reader; this used to be a
                        # local reimplementation, so the cell format had two
                        # independent parsers and only one of them knew 'del'.
                        for piece in split_piece_cell(intron):
                            if ":" not in piece:
                                continue
                            label, tok = piece.split(":", 1)
                            per_gene.setdefault(gene, {}).setdefault(label, []).append((mut, tok))
                        # "Already covered on the transcript substrate" means HAS
                        # AN ORF ADDRESS, not "is exonic". The transcript mapping
                        # is built from ORF-space tokens only, so a gd./ch. variant
                        # in a UTR or non-coding exon has no transcript row and no
                        # aa row -- pre-mRNA is its ONLY nucleotide coordinate.
                        # Gating on a non-empty `intron` cell skipped exactly those,
                        # and since they still appear in the chromosome and gDNA
                        # CSVs they did not look dropped. An exonic row that DOES
                        # carry an orf token is still excluded, so nothing is
                        # double-counted.
                        if intron or not orf:
                            for piece in split_piece_cell(pre):
                                per_gene.setdefault(gene, {}).setdefault("pre_mRNA", []).append((mut, piece))
            except Exception as exc:
                print(f"[MUT] {f.name}: unreadable ({exc})", file=sys.stderr)

        for gene, by_sub in per_gene.items():
            for label, rows in by_sub.items():
                key = _gene_key(gene, label)
                if key not in substrate_seqs:
                    print(f"[MUT] {key}: no FASTA record, skipping {len(rows)} token(s)",
                          file=sys.stderr)
                    continue
                mapping_dict[key] = pd.DataFrame(rows, columns=["mutant", "transcript"])
                wt_sequences[key] = substrate_seqs[key]
                if args.strict_introns:
                    strict_substrates.add(label)

    os.makedirs(args.output, exist_ok=True)

    # 1) WT
    run_wt_phase(
        wt_sequences=wt_sequences,
        outdir=args.output,
        miranda_dir=args.miranda_dir,
        mirna_db=args.mirna_db,
        strict_substrates=strict_substrates,
    )

    # 2) MUT
    mut_summary = run_mut_phase(
        wt_sequences=wt_sequences,
        mapping_dict=mapping_dict,
        outdir=args.output,
        miranda_dir=args.miranda_dir,
        mirna_db=args.mirna_db,
        use_parallel=not args.no_parallel,
        max_workers=args.max_workers,
        failure_map=failure_map,
        strict_substrates=strict_substrates,
    )

    # 3) PARSE
    print("\n[PARSE] Starting ...")
    sites_df, parse_rejected = build_sites_table_from_outputs(args.output, mapping_dict, failure_map)
    sites_df = assign_loci(sites_df)
    sites_df = assign_segments(sites_df)

    # join_pos/align_status/variant_kind/length_delta are appended, never inserted:
    # a consumer that reads these TSVs positionally keeps working.
    sites_cols = [
        "pkey", "allele", "mirna_id", "site_pos",
        "tot_score", "tot_energy", "max_score", "max_energy",
        "strand", "len_mirna", "len_target",
        "visibility_flag", "distance_to_snv",
        "locus_id", "segment_id",
        "parser_confidence", "run_meta",
        "join_pos", "align_status", "variant_kind", "length_delta", "substrate"
    ]

    events_df = build_events(sites_df)
    events_cols = [
        "pkey","mirna_id","locus_id","segment_id",
        "wt_pos","mut_pos","dpos",
        "wt_tot_score","mut_tot_score","delta_tot_score","pct_delta",
        "wt_tot_energy","mut_tot_energy","delta_energy",
        "distance_to_snv",
        "rank_wt","rank_mut",
        "conf_wt","conf_mut","conf_weighted_delta",
        "cls","is_high_impact","priority","in_radius",
        "variant_kind","length_delta","qc_flags","substrate"
    ]

    summary_df = build_summary(events_df, sites_df)
    summary_cols = [
        "pkey", "n_hits_wt", "n_hits_mut", "n_mirna", "n_loci", "n_segments",
        "n_competitive_segments", "n_new_competitors",
        "global_count_gained_high", "global_count_lost_high", "global_count_shifted",
        "global_max_abs_delta_score", "global_sum_weighted_abs_delta",
        "nearest_event_bp_any",
        "local_count_gained_high", "local_count_lost_high", "local_count_shifted",
        "local_max_abs_delta_score", "nearest_event_bp_local", "frac_effect_in_radius",
        "max_segment_abs_delta_best", "frac_effect_in_competitive_segments",
        "top_event_mirna", "top_event_class", "top_event_delta_score", "top_event_pos", "qc_flags",
        "variant_kind", "length_delta", "n_sites_aligned", "n_sites_deleted",
        "n_sites_inserted", "variant_qc", "substrate"
    ]

    # Write per-gene output files
    # One output directory per REAL gene: an intronic row belongs to the same
    # gene's Miranda folder as its exonic rows, distinguished by the substrate
    # column rather than by a separate directory.
    gene_names = sorted({_split_gene_key(k)[0] for k in wt_sequences})
    for gene_upper in gene_names:
        out_dir = Path(args.output) / gene_upper / "Miranda"
        out_dir.mkdir(parents=True, exist_ok=True)
        summary_path = str(out_dir / f"{gene_upper}.tsv")
        events_path  = str(out_dir / f"{gene_upper}.events.tsv")
        sites_path   = str(out_dir / f"{gene_upper}.sites.tsv")

        gene_prefix = f"{gene_upper}-"

        g_sites = (sites_df[sites_df['pkey'].str.startswith(gene_prefix)]
                   if sites_df is not None and not sites_df.empty
                   else pd.DataFrame(columns=sites_cols))
        g_events = (events_df[events_df['pkey'].str.startswith(gene_prefix)]
                    if events_df is not None and not events_df.empty
                    else pd.DataFrame(columns=events_cols))
        g_summary = (summary_df[summary_df['pkey'].str.startswith(gene_prefix)]
                     if summary_df is not None and not summary_df.empty
                     else pd.DataFrame(columns=summary_cols))

        for c in sites_cols:
            if c not in g_sites.columns:
                g_sites[c] = None
        g_sites[sites_cols].to_csv(sites_path, sep="\t", index=False)

        for c in events_cols:
            if c not in g_events.columns:
                g_events[c] = None
        g_events[events_cols].to_csv(events_path, sep="\t", index=False)

        for c in summary_cols:
            if c not in g_summary.columns:
                g_summary[c] = None
        g_summary[summary_cols].to_csv(summary_path, sep="\t", index=False)

        print(f"[OUT] Wrote {gene_upper}: {summary_path}")

    # Run-level accounting. Every token that produced no row is named here with
    # its reason; the pipeline previously `continue`d past all of them, so an
    # indel, a REF mismatch and a mutation the gene never had were
    # indistinguishable from each other and from success.
    run_summary = {
        "genes_total": mut_summary["genes_total"],
        "genes_skipped": mut_summary["genes_skipped"],
        "tokens_total": mut_summary["tokens_total"],
        "tokens_built": mut_summary["tokens_built"],
        "tokens_rejected": mut_summary["tokens_rejected"],
        "runs_ok": mut_summary["runs_ok"],
        "runs_failed": mut_summary["runs_failed"],
        "parse_rejected": len(parse_rejected),
        "pkeys_with_rows": int(summary_df["pkey"].nunique()) if (summary_df is not None and not summary_df.empty) else 0,
        "skipped_genes": mut_summary["skipped_genes"],
        "rejected_tokens": mut_summary["rejected"],
        "run_failures": mut_summary["run_failures"],
        "rejected_outputs": parse_rejected,
    }
    with open(Path(args.output) / "miranda.run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)
    print(f"[SUMMARY] genes {run_summary['genes_total'] - run_summary['genes_skipped']}"
          f"/{run_summary['genes_total']} with a mapping, tokens {run_summary['tokens_built']} built / "
          f"{run_summary['tokens_rejected']} rejected of {run_summary['tokens_total']}, "
          f"{run_summary['runs_failed']} run failure(s), {run_summary['parse_rejected']} "
          f"output(s) rejected at parse -> miranda.run_summary.json")

    print("\n[DONE] Unified MirandA pipeline completed.")

if __name__ == "__main__":
    main()
