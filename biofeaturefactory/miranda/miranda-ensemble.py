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
- Redundancy reduction: leverage utility.load_mapping, extract_mutation_from_sequence_name, extract_gene_from_filename, get_mutation_data, update_str, load_wt_sequences, load_validation_failures, should_skip_mutation, validate_mapping_content.
"""

import os
import sys
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
from biofeaturefactory.utils.utility import (
    read_fasta,
    get_mutation_data,
    load_validation_failures,
    should_skip_mutation,
    extract_gene_from_filename,
    extract_mutation_from_sequence_name,
    update_str,
    validate_mapping_content,
    load_wt_sequences,
    load_mapping,
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
def safe_div(num: float, den: float) -> float:
    return num / den if den not in (0, None, 0.0) else 0.0

def exp_weight(d: Optional[float], k: float) -> float:
    if d is None:
        return 1.0
    return math.exp(-(d / k)) if k else 1.0

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
        files = sorted(f for f in path.rglob("*.csv") if f.is_file() and f.suffix.lower() in (".csv", ".tsv", ".txt"))
    for f in files:
        try:
            validate_mapping_content(str(f))
        except Exception:
            continue

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
def _run_single_miranda_task(args: Tuple[str, str, str, str, str]) -> Tuple[str, str]:
    """
    args = (seq_id, seq_string, miranda_dir, outdir, mirna_db_abs)
    miranda_dir may be empty string/None to use miranda from PATH.
    """
    import shutil as _shutil
    import tempfile as _tempfile
    from subprocess import run, TimeoutExpired
    seq_id, seq_str, miranda_dir, outdir, mirna_db_abs = args

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

    tmp_file = os.path.join(tmp_base, f"tmp_{os.getpid()}_{hash(seq_id) & 0xFFFFF}.fasta")
    try:
        with open(tmp_file, "w") as f:
            f.write(f">{seq_id}\n{seq_str}")

        # When using cwd, miranda expects basename; when on PATH, use absolute path
        fasta_arg = os.path.basename(tmp_file) if miranda_dir else tmp_file
        cmd = [cmd_binary, mirna_db_abs, fasta_arg]
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
                 mirna_db: str):
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
                        already_on_disk: set) -> List[Tuple[str, str, str, str, str]]:
    """
    Build per-gene list of MirandA tasks for MUT allele, skipping disk.
    """
    tasks: List[Tuple[str, str, str, str, str]] = []
    db_abs = os.path.abspath(mirna_db)

    # Prebuild lookup once per gene (no per-file DataFrame filtering)
    mp = dict(zip(mapping_df["mutant"].astype(str), mapping_df["transcript"].astype(str)))

    for mutant_id, transcript_token in mp.items():
        mutant_id = str(mutant_id).strip()
        transcript_token = str(transcript_token).strip()
        if not mutant_id or not transcript_token:
            continue

        seq_id = f"{gene_upper}-{mutant_id}-mut"
        if seq_id in already_on_disk:
            continue

        if failure_map and should_skip_mutation(gene_upper, mutant_id, failure_map):
            continue

        try:
            tx_pos, (ref_nt, alt_nt) = get_mutation_data(transcript_token)
        except Exception:
            continue

        if not (0 <= tx_pos < len(wt_seq)):
            continue

        if ref_nt and wt_seq[tx_pos].upper() != ref_nt.upper():
            continue

        mut_seq = update_str(wt_seq, alt_nt, tx_pos)
        out_file = os.path.join(outdir, f"{seq_id}-miranda.out")
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            continue

        tasks.append((seq_id, mut_seq, miranda_dir, outdir, db_abs))
    return tasks

def run_mut_phase(wt_sequences: Dict[str, str],
                  mapping_dict: Dict[str, pd.DataFrame],
                  outdir: str,
                  miranda_dir: str,
                  mirna_db: str,
                  use_parallel: bool = True,
                  max_workers: Optional[int] = None,
                  failure_map: Optional[Dict] = None):
    """
    MUT: per-gene parallel; progress prints seeded with prior completions.
    """
    print("[MUT] Starting MUT MirandA phase...")
    os.makedirs(outdir, exist_ok=True)

    genes = sorted([g for g in wt_sequences.keys() if g in mapping_dict and not mapping_dict[g].empty])
    total_genes = len(genes)

    total_rows = sum(len(mapping_dict[g]) for g in genes)
    print(f"[MUT] Prepared estimates using mapping rows: {total_rows}")

    if max_workers is None:
        max_workers = min(max(multiprocessing.cpu_count() - 1, 1), MAX_PARALLEL_CAP)
    print(f"[MUT] Using {max_workers} workers (per batch)")

    for idx, gene in enumerate(genes, 1):
        gene_upper = gene.upper()
        print(f"[MUT] [{idx}/{total_genes}] Starting {gene_upper} processing ...")

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

        tasks = _per_gene_mut_tasks(
            gene_upper=gene_upper,
            wt_seq=wt_seq,
            mapping_df=mapping_df,
            outdir=outdir,
            miranda_dir=miranda_dir,
            mirna_db=mirna_db,
            failure_map=failure_map,
            already_on_disk=already_on_disk,
        )
        remaining = len(tasks)

        if remaining == 0:
            print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: already complete [{done_prior}/{total_all}]")
            continue

        print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: remaining {remaining} of {total_all} (done={done_prior})")

        processed = done_prior
        completed_ids: List[str] = []

        if not use_parallel:
            for t in tasks:
                seq_id, status = _run_single_miranda_task(t)
                if status not in ("ok", "already"):
                    print(f"[MUT]   {seq_id}: {status}")
                else:
                    completed_ids.append(seq_id)
                processed += 1
                if processed % GENE_PROGRESS_STEP == 0 or processed == total_all:
                    print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: Processed mutations - [{processed}/{total_all}]")
                if processed % PROGRESS_SAVE_EVERY == 0 or processed == total_all:
                    _write_gene_progress(outdir, gene_upper, processed, total_all, completed_ids)
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
                # map with chunksize=1 for steady progress feedback
                for seq_id, status in ex.map(_run_single_miranda_task, tasks):
                    if status not in ("ok", "already"):
                        print(f"[MUT]   {seq_id}: {status}")
                    else:
                        completed_ids.append(seq_id)
                    processed += 1
                    if processed % GENE_PROGRESS_STEP == 0 or processed == total_all:
                        print(f"[MUT] [{idx}/{total_genes}] {gene_upper}: Processed mutations - [{processed}/{total_all}]")
                    if processed % PROGRESS_SAVE_EVERY == 0 or processed == total_all:
                        _write_gene_progress(outdir, gene_upper, processed, total_all, completed_ids)

    print("[MUT] Done.")

# -------------------------------------------------------------------------
# PARSING
# -------------------------------------------------------------------------
def parse_miranda_text(miranda_text: str) -> List[Dict]:
    """
    Parse MirandA stdout -> site rows.

    miRNA ID is taken from the '>>' summary line (parts[0][2:] strips '>>'),
    NOT from the 'Query: 3'..5'' alignment line which contains the alignment
    sequence (not the ID) and produces a spurious trailing ' 5' artifact.
    """
    sites: List[Dict] = []
    query_seq = ""
    ref_seq = ""

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

        if line.startswith(">>"):
            norm = re.sub(r"\s+", "\t", line)
            parts = norm.split("\t")
            if len(parts) < 10:
                continue

            # miRNA ID is the first field with '>>' stripped.
            mirna_id = parts[0][2:] if parts[0].startswith(">>") else parts[0]
            if not mirna_id:
                continue

            def f(i):
                try: return float(parts[i])
                except Exception: return None

            def i_(i):
                try: return int(parts[i])
                except Exception: return None

            tot_score  = f(2)
            tot_energy = f(3)
            max_score  = f(4)
            max_energy = f(5)
            strand     = parts[6] if len(parts) > 6 else ""
            len1       = i_(7)
            len2       = i_(8)
            # re.sub converted spaces inside the positional field to tabs,
            # so rejoin parts[9:] to recover all space-separated positions.
            pos_field  = " ".join(parts[9:])
            pos_tokens = [p for p in pos_field.strip().split() if p.isdigit()]
            if not pos_tokens:
                continue

            for p in pos_tokens:
                site_pos = int(p)
                sites.append({
                    "mirna_id": mirna_id,
                    "site_pos": site_pos,
                    "tot_score": tot_score,
                    "tot_energy": tot_energy,
                    "max_score": max_score,
                    "max_energy": max_energy,
                    "strand": strand,
                    "len_mirna": len1,
                    "len_target": len2,
                    "query_seq": query_seq,
                    "ref_seq": ref_seq,
                    "parser_confidence": 1.0 if (mirna_id and tot_score is not None) else 0.5,
                })
    return sites

def _collect_gene_outputs(outdir: str) -> Dict[str, Dict[str, List[Path]]]:
    """
    { GENE: { 'wt': [Path], 'mut': [List[Path]] } }
    WT:  GENE-wt-miranda.out
    MUT: GENE-<mut>-mut-miranda.out
    """
    buckets: Dict[str, Dict[str, List[Path]]] = {}
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
                continue
            gene = core.split("-", 1)[0].upper()
            buckets.setdefault(gene, {}).setdefault("mut", []).append(Path(f.path))
    return buckets

def build_sites_table_from_outputs(outdir: str,
                                   mapping_dict: Dict[str, pd.DataFrame],
                                   failure_map: Dict) -> pd.DataFrame:
    """
    Build sites by joining single per-gene WT output with each mutation-specific MUT output.
    WT rows are materialized per pkey (allele='WT'); MUT rows parsed from files (allele='MUT').
    Intermediate .out files are deleted after parsing to reclaim disk space.
    """
    print("[PARSE] Building sites table from MirandA outputs...")
    records: List[Dict] = []

    buckets = _collect_gene_outputs(outdir)
    genes = sorted(buckets.keys())
    print(f"[PARSE] Found outputs for {len(genes)} genes")

    miss_map = miss_hits = 0
    parsed_files = 0

    for gi, gene in enumerate(genes, 1):
        m = buckets[gene]
        wt_files = m.get("wt", [])
        mut_files = m.get("mut", [])
        if not wt_files:
            print(f"[PARSE]   {gene}: missing WT output; skipping {len(mut_files)} MUT files")
            continue
        wt_file = wt_files[0]

        wt_text = wt_file.read_text(encoding="utf-8", errors="ignore")
        wt_hits = parse_miranda_text(wt_text)

        mapping_df = mapping_dict.get(gene)
        if mapping_df is None or mapping_df.empty:
            print(f"[PARSE]   {gene}: no mapping, skip")
            continue

        # Precompute mutant->transcript map once
        mp = dict(zip(mapping_df["mutant"].astype(str), mapping_df["transcript"].astype(str)))

        print(f"[PARSE]   {gene}: WT hits={len(wt_hits)}, MUT files={len(mut_files)}")
        for mf in mut_files:
            parsed_files += 1

            # GENE-<mut>-mut
            base = mf.name[:-len("-miranda.out")]
            core = base[:-4]  # GENE-<mut>
            if "-" not in core:
                miss_map += 1
                continue
            # Use utility extractor when possible
            try:
                gene_tok, mut_tok = extract_mutation_from_sequence_name(core)
                mutation_token = mut_tok
            except Exception:
                # Fallback split
                mutation_token = core.split("-", 1)[1]

            transcript_nt = mp.get(mutation_token, None)
            if not transcript_nt:
                miss_map += 1
                continue

            if should_skip_mutation(gene, mutation_token, failure_map):
                continue

            try:
                tx_pos_0, _ = get_mutation_data(str(transcript_nt))
                tx_pos_1 = tx_pos_0 + 1
            except Exception:
                tx_pos_1 = None

            # WT rows for this pkey
            for h in wt_hits:
                dist = abs(h["site_pos"] - tx_pos_1) if (tx_pos_1 is not None and h["site_pos"] is not None) else None
                records.append({
                    "pkey": f"{gene}-{mutation_token}",
                    "allele": "WT",
                    "mirna_id": h["mirna_id"],
                    "site_pos": h["site_pos"],
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
                    "run_meta": f"{gene}-wt",
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
                continue

            for h in mut_hits:
                dist = abs(h["site_pos"] - tx_pos_1) if (tx_pos_1 is not None and h["site_pos"] is not None) else None
                records.append({
                    "pkey": f"{gene}-{mutation_token}",
                    "allele": "MUT",
                    "mirna_id": h["mirna_id"],
                    "site_pos": h["site_pos"],
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
                    "run_meta": f"{gene}-mut",
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
    return sites_df

def assign_loci(sites_df: pd.DataFrame) -> pd.DataFrame:
    print("[PARSE] Assigning per-miRNA loci...")
    if sites_df is None or sites_df.empty or "pkey" not in sites_df.columns:
        print("[PARSE] No sites to assign loci")
        return sites_df
    updates = 0
    for pkey, sub in sites_df.groupby("pkey", dropna=False):
        for mirna_id, g in sub.groupby("mirna_id", dropna=False):
            g = g.sort_values("site_pos")
            current = 1
            last_pos = None
            for idx, row in g.iterrows():
                pos = row["site_pos"]
                if last_pos is None or (pos - last_pos) > MERGE_WINDOW_NT:
                    locus_id = f"m{current}"
                    current += 1
                else:
                    locus_id = f"m{current-1}"
                last_pos = pos
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
    for pkey, sub in sites_df.groupby("pkey", dropna=False):
        g = sub.sort_values("site_pos")
        current = 1
        last_pos = None
        for idx, row in g.iterrows():
            pos = row["site_pos"]
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
    print("[PARSE] Building events...")
    cols = [
        "pkey","mirna_id","locus_id","segment_id",
        "wt_pos","mut_pos","dpos",
        "wt_tot_score","mut_tot_score","delta_tot_score","pct_delta",
        "wt_tot_energy","mut_tot_energy","delta_energy",
        "distance_to_snv",
        "rank_wt","rank_mut",
        "conf_wt","conf_mut","conf_weighted_delta",
        "cls","is_high_impact","priority","in_radius"
    ]
    if sites_df is None or sites_df.empty or not {"pkey","mirna_id","locus_id","allele"}.issubset(sites_df.columns):
        print("[PARSE] No sites to build events")
        return pd.DataFrame(columns=cols)

    key = ["pkey","mirna_id","locus_id"]
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
        seg = r.get("segment_id_mut") if pd.notnull(r.get("segment_id_mut")) else r.get("segment_id_wt")

        wt_pos  = int(r["site_pos_wt"]) if pd.notnull(r.get("site_pos_wt")) else None
        mut_pos = int(r["site_pos_mut"]) if pd.notnull(r.get("site_pos_mut")) else None
        dpos = (mut_pos - wt_pos) if (wt_pos is not None and mut_pos is not None) else None

        wt_score = float(r["tot_score_wt"]) if pd.notnull(r.get("tot_score_wt")) else 0.0
        mut_score = float(r["tot_score_mut"]) if pd.notnull(r.get("tot_score_mut")) else 0.0
        delta = mut_score - wt_score
        denom = max(wt_score, mut_score) if (wt_score or mut_score) else 0.0
        pct_delta = (delta / denom) if denom else 0.0

        wt_energy = float(r["tot_energy_wt"]) if pd.notnull(r.get("tot_energy_wt")) else None
        mut_energy = float(r["tot_energy_mut"]) if pd.notnull(r.get("tot_energy_mut")) else None
        delta_energy = (mut_energy - wt_energy) if (mut_energy is not None and wt_energy is not None) else None

        dist = (float(r["distance_to_snv_mut"]) if pd.notnull(r.get("distance_to_snv_mut"))
                else (float(r["distance_to_snv_wt"]) if pd.notnull(r.get("distance_to_snv_wt")) else None))

        wt_vis = (wt_score >= VISIBILITY_THRESHOLD)
        mut_vis = (mut_score >= VISIBILITY_THRESHOLD)

        if not wt_vis and mut_vis:
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
        is_high = 1 if (abs(delta) >= HIGH_CUTOFF) or cls in ("gained","lost") else 0

        w = exp_weight(dist, DISTANCE_K)
        class_bonus = {"gained":3.0,"lost":3.0,"shifted":1.0,"strengthened":0.5,"weakened":0.5,"none":0.0}[cls]
        priority = abs(delta) * w + class_bonus

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
            "conf_wt": 1.0 if wt_vis else 0.5,
            "conf_mut": 1.0 if mut_vis else 0.5,
            "conf_weighted_delta": (1.0 if mut_vis else 0.5)*mut_score - (1.0 if wt_vis else 0.5)*wt_score,
            "cls": cls,
            "is_high_impact": is_high,
            "priority": priority,
            "in_radius": in_radius,
        })

    events_df = pd.DataFrame(out, columns=cols)
    if events_df.empty:
        print("[PARSE] Events table built with 0 rows")
        return events_df

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
            "top_event_mirna", "top_event_class", "top_event_delta_score", "top_event_pos", "qc_flags"
        ])

    summary_records: List[Dict] = []

    mut_vis_counts = {}
    wt_vis_counts = {}
    if sites_df is not None and not sites_df.empty:
        mut_vis = sites_df[(sites_df["allele"]=="MUT") & (sites_df["visibility_flag"]==1)].groupby("pkey", dropna=False).size()
        wt_vis  = sites_df[(sites_df["allele"]=="WT") & (sites_df["visibility_flag"]==1)].groupby("pkey", dropna=False).size()
        mut_vis_counts = mut_vis.to_dict()
        wt_vis_counts = wt_vis.to_dict()

    for pkey, ev in events_df.groupby("pkey", dropna=False):
        sites_sub = sites_df[sites_df["pkey"] == pkey] if (sites_df is not None and not sites_df.empty) else pd.DataFrame()

        n_hits_wt = wt_vis_counts.get(pkey, 0)
        n_hits_mut = mut_vis_counts.get(pkey, 0)
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

        if len(ev):
            global_max_abs_delta_score = ev["delta_tot_score"].abs().max()
            weighted_abs = 0.0
            total_abs = 0.0
            for _, r in ev.iterrows():
                absd = abs(r["delta_tot_score"])
                total_abs += absd
                d = r["distance_to_snv"]
                weighted_abs += absd * exp_weight(d, DISTANCE_K)
            global_sum_weighted_abs_delta = weighted_abs
        else:
            global_max_abs_delta_score = 0.0
            global_sum_weighted_abs_delta = 0.0
            total_abs = 0.0

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

        if len(ev):
            top = ev.sort_values("priority", ascending=False).iloc[0]
            top_event_mirna = top["mirna_id"]
            top_event_class = top["cls"]
            top_event_delta_score = top["delta_tot_score"]
            top_event_pos = top["mut_pos"] if pd.notnull(top.get("mut_pos")) else top.get("wt_pos")
        else:
            top_event_mirna = None
            top_event_class = "none"
            top_event_delta_score = 0.0
            top_event_pos = None

        qc_flags = []
        if (n_hits_wt + n_hits_mut) == 0:
            qc_flags.append("no_hits")
        if nearest_event_bp_any is not None and nearest_event_bp_any > 2000:
            qc_flags.append("far_event>2kb")
        if global_max_abs_delta_score == 0.0:
            qc_flags.append("low_signal_only")

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
        })

    summary_df = pd.DataFrame.from_records(summary_records)
    print(f"[PARSE] Summary table built with {len(summary_df)} rows")
    return summary_df

# -------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------
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

    os.makedirs(args.output, exist_ok=True)

    # 1) WT
    run_wt_phase(
        wt_sequences=wt_sequences,
        outdir=args.output,
        miranda_dir=args.miranda_dir,
        mirna_db=args.mirna_db,
    )

    # 2) MUT
    run_mut_phase(
        wt_sequences=wt_sequences,
        mapping_dict=mapping_dict,
        outdir=args.output,
        miranda_dir=args.miranda_dir,
        mirna_db=args.mirna_db,
        use_parallel=not args.no_parallel,
        max_workers=args.max_workers,
        failure_map=failure_map,
    )

    # 3) PARSE
    print("\n[PARSE] Starting ...")
    sites_df = build_sites_table_from_outputs(args.output, mapping_dict, failure_map)
    sites_df = assign_loci(sites_df)
    sites_df = assign_segments(sites_df)

    sites_cols = [
        "pkey", "allele", "mirna_id", "site_pos",
        "tot_score", "tot_energy", "max_score", "max_energy",
        "strand", "len_mirna", "len_target",
        "visibility_flag", "distance_to_snv",
        "locus_id", "segment_id",
        "parser_confidence", "run_meta"
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
        "cls","is_high_impact","priority","in_radius"
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
        "top_event_mirna", "top_event_class", "top_event_delta_score", "top_event_pos", "qc_flags"
    ]

    # Write per-gene output files
    gene_names = list(wt_sequences.keys())
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

    print("\n[DONE] Unified MirandA pipeline completed.")

if __name__ == "__main__":
    main()
