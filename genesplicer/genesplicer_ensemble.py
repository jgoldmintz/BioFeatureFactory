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
GeneSplicer WT↔ALT ensemble delta caller
- Runs GeneSplicer on full genomic sequences (or windowed/custom) for each gene.
"""

import os
import sys
import argparse
import tempfile
import subprocess
from pathlib import Path
import pandas as pd
import math
import multiprocessing
import concurrent.futures

# ---------------------------------------------------------------------------
# imports from ../utils/utility.py
# ---------------------------------------------------------------------------
DEP_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "utils"))
if DEP_DIR not in sys.path:
    sys.path.insert(0, DEP_DIR)

from utility import (
    read_fasta,
    update_str,
    get_mutation_data,
    load_validation_failures,
    should_skip_mutation,
    extract_gene_from_filename,
)

# ---------------------------------------------------------------------------
# defaults
# ---------------------------------------------------------------------------
DEFAULT_WINDOW = 151
DEFAULT_REPORT_RADIUS = None  # fallback to window
DEFAULT_DISTANCE_K = 75
DEFAULT_VISIBILITY_THRESHOLD = 1.0
DEFAULT_HIGH_CUTOFF = 5.0
DEFAULT_SHIFT_BP = 3
DEFAULT_MAX_WORKERS = 8

# ---------------------------------------------------------------------------
# helpers: run genesplicer
# ---------------------------------------------------------------------------

def _run_genesplicer_on_seq(seq_name: str, seq: str, genesplicer_dir: str, model_rel: str = "../human") -> pd.DataFrame:
    """
    Run GeneSplicer on a single in-memory sequence.
    Returns a DataFrame with columns:
        End5, End3, Score, confidence, splice_site_type
    """
    os.makedirs(genesplicer_dir, exist_ok=True)
    # use temp file in genesplicer_dir to avoid path issues
    with tempfile.NamedTemporaryFile(mode="w", dir=genesplicer_dir, delete=False, suffix=".fasta") as tmpf:
        tmp_name = os.path.basename(tmpf.name)
        tmpf.write(f">{seq_name}\n{seq}\n")
        tmp_path = tmpf.name

    out = None
    try:
        cmd = f"cd {genesplicer_dir} && ./genesplicer {tmp_name} {model_rel}"
        # capture stdout
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        out = result.stdout.strip()
    finally:
        # cleanup temp file
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    if not out:
        return pd.DataFrame(columns=["End5", "End3", "Score", "confidence", "splice_site_type"])

    # GeneSplicer output is space-separated: End5 End3 Score confidence splice_site_type
    rows = []
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        toks = line.split()
        if len(toks) < 5:
            # malformed line, skip
            continue
        End5, End3, Score, conf, sstype = toks[0], toks[1], toks[2], toks[3], toks[4]
        try:
            rows.append({
                "End5": int(End5),
                "End3": int(End3),
                "Score": float(Score),
                "confidence": conf,
                "splice_site_type": sstype,
            })
        except ValueError:
            # ignore bad lines
            continue

    if not rows:
        return pd.DataFrame(columns=["End5", "End3", "Score", "confidence", "splice_site_type"])

    return pd.DataFrame(rows)


def _confidence_to_weight(conf: str) -> float:
    if conf is None:
        return 0.75
    c = conf.lower()
    if c.startswith("h"):
        return 1.0
    if c.startswith("m"):
        return 0.75
    if c.startswith("l"):
        return 0.5
    return 0.75


# ---------------------------------------------------------------------------
# mapping loader
# ---------------------------------------------------------------------------

def load_mapping_dir(mapping_dir: str):
    """
    Load all mapping files in mapping_dir into a dict: {gene_name: DataFrame}
    Expected columns include at least: 'mutant' and 'genomic' as in your example.
    """
    mapping_dir = Path(mapping_dir)
    if not mapping_dir.exists():
        raise FileNotFoundError(f"Mapping dir not found: {mapping_dir}")
    maps = {}
    for f in mapping_dir.iterdir():
        if not f.is_file():
            continue
        if f.suffix not in (".tsv", ".csv", ".txt"):
            continue
        df = pd.read_csv(f, sep="\t" if f.suffix == ".tsv" else ",")
        gene = extract_gene_from_filename(str(f))
        if not gene:
            # fall back to stem
            gene = f.stem
        maps[gene] = df
    return maps


# ---------------------------------------------------------------------------
# site builders
# ---------------------------------------------------------------------------

def _build_sites_for_allele(pkey: str,
                            allele: str,
                            snv_pos: int,
                            df: pd.DataFrame,
                            visibility_threshold: float,
                            report_radius: int,
                            scan_mode: str) -> pd.DataFrame:
    """
    Convert GeneSplicer raw DF into the canonical per-site schema.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=[
            "pkey", "allele", "type", "site_pos", "score", "confidence",
            "rank", "distance_to_snv", "visible_flag", "cluster_id", "in_radius"
        ])

    rows = []
    for _, r in df.iterrows():
        sstype = r["splice_site_type"]
        score = float(r["Score"])
        conf = r["confidence"]
        conf_w = _confidence_to_weight(conf)
        # GeneSplicer: donor keyed by End5; acceptor keyed by End3
        if sstype == "donor":
            site_pos = int(r["End5"])
        else:
            site_pos = int(r["End3"])

        dist = abs(site_pos - snv_pos)
        in_radius = 0
        if report_radius is not None:
            in_radius = 1 if dist <= report_radius else 0

        rows.append({
            "pkey": pkey,
            "allele": allele,
            "type": sstype,
            "site_pos": site_pos,
            "score": score,
            "confidence": conf_w,
            "distance_to_snv": dist,
            "visible_flag": 1 if score >= visibility_threshold else 0,
            "cluster_id": None,  # filled later
            "in_radius": in_radius,
        })

    # rank within allele×type by score
    df_sites = pd.DataFrame(rows)
    if not df_sites.empty:
        for t in ("donor", "acceptor"):
            mask = df_sites["type"] == t
            df_sites.loc[mask, "rank"] = (
                df_sites[mask]
                .sort_values("score", ascending=False)
                .reset_index(drop=True)
                .reset_index()
                .set_index(df_sites[mask].index)["index"] + 1
            )
    else:
        df_sites["rank"] = pd.NA

    return df_sites


# ---------------------------------------------------------------------------
# clustering and event detection
# ---------------------------------------------------------------------------

def _cluster_sites(sites_df: pd.DataFrame, cluster_radius: int) -> pd.DataFrame:
    """
    Assign cluster_id per (pkey, allele, type) using single-linkage by position.
    cluster_id is the same across alleles later when we merge WT and MUT.
    """
    if sites_df is None or sites_df.empty:
        return sites_df

    sites_df = sites_df.sort_values(["pkey", "type", "site_pos", "allele"]).reset_index(drop=True)

    cluster_ids = []
    current_key = None
    current_cluster = 0
    last_pos = None

    for idx, row in sites_df.iterrows():
        key = (row["pkey"], row["type"])
        pos = int(row["site_pos"])
        if key != current_key:
            current_key = key
            current_cluster = 1
            last_pos = pos
            cluster_ids.append(f"{row['type'][0]}{current_cluster}")
            continue

        # same key
        if abs(pos - last_pos) <= cluster_radius:
            # same cluster
            cluster_ids.append(f"{row['type'][0]}{current_cluster}")
        else:
            current_cluster += 1
            cluster_ids.append(f"{row['type'][0]}{current_cluster}")
        last_pos = pos

    sites_df["cluster_id"] = cluster_ids
    return sites_df


def _pair_clusters_to_events(sites_df: pd.DataFrame,
                             visibility_threshold: float,
                             shift_bp: int) -> pd.DataFrame:
    """
    From per-site rows, produce per-cluster per-pkey events.
    For each (pkey, type, cluster_id):
        - take top WT (by score) and top MUT
        - compute deltas
        - classify
    """
    if sites_df is None or sites_df.empty:
        return pd.DataFrame(columns=[
            "pkey", "type", "cluster_id",
            "wt_pos", "mut_pos", "dpos",
            "wt_score", "mut_score", "dscore", "pct_delta",
            "distance_to_snv",
            "rank_wt", "rank_mut",
            "conf_wt", "conf_mut", "conf_weighted_delta",
            "cls", "is_high_impact", "priority", "in_radius",
        ])

    events = []
    grouped = sites_df.groupby(["pkey", "type", "cluster_id"], dropna=False)

    for (pkey, sstype, cluster_id), grp in grouped:
        # split by allele
        wt_grp = grp[grp["allele"] == "WT"].sort_values("score", ascending=False)
        mut_grp = grp[grp["allele"] == "MUT"].sort_values("score", ascending=False)

        wt_site = wt_grp.iloc[0] if not wt_grp.empty else None
        mut_site = mut_grp.iloc[0] if not mut_grp.empty else None

        wt_visible = wt_site is not None and float(wt_site["score"]) >= visibility_threshold
        mut_visible = mut_site is not None and float(mut_site["score"]) >= visibility_threshold

        wt_pos = int(wt_site["site_pos"]) if wt_site is not None else None
        mut_pos = int(mut_site["site_pos"]) if mut_site is not None else None
        wt_score = float(wt_site["score"]) if wt_site is not None else None
        mut_score = float(mut_site["score"]) if mut_site is not None else None
        conf_wt = float(wt_site["confidence"]) if wt_site is not None else 0.0
        conf_mut = float(mut_site["confidence"]) if mut_site is not None else 0.0
        dist = None
        if wt_site is not None and mut_site is not None:
            dist = min(int(wt_site["distance_to_snv"]), int(mut_site["distance_to_snv"]))
        elif wt_site is not None:
            dist = int(wt_site["distance_to_snv"])
        elif mut_site is not None:
            dist = int(mut_site["distance_to_snv"])

        dpos = None
        if wt_pos is not None and mut_pos is not None:
            dpos = mut_pos - wt_pos

        dscore = None
        if wt_score is not None and mut_score is not None:
            dscore = mut_score - wt_score
        elif wt_score is None and mut_score is not None:
            dscore = mut_score
        elif wt_score is not None and mut_score is None:
            dscore = -wt_score

        # pct delta
        if wt_score is not None and abs(wt_score) > 1e-6 and dscore is not None:
            pct_delta = dscore / abs(wt_score)
        else:
            pct_delta = None

        # ranks
        rank_wt = int(wt_site["rank"]) if wt_site is not None and not pd.isna(wt_site["rank"]) else None
        rank_mut = int(mut_site["rank"]) if mut_site is not None and not pd.isna(mut_site["rank"]) else None

        # class
        if (not wt_visible) and mut_visible:
            cls = "gained"
        elif wt_visible and (not mut_visible):
            cls = "lost"
        elif wt_visible and mut_visible and dpos is not None and abs(dpos) >= shift_bp:
            cls = "shifted"
        elif wt_visible and mut_visible and dscore is not None and abs(dscore) >= 1.0:
            cls = "strengthened" if dscore > 0 else "weakened"
        else:
            cls = "none"

        # in_radius
        in_radius = int(any(grp["in_radius"] == 1)) if "in_radius" in grp.columns else 0

        # conf-weighted delta
        conf_weighted_delta = None
        if mut_score is not None or wt_score is not None:
            conf_weighted_delta = (conf_mut * (mut_score or 0.0)) - (conf_wt * (wt_score or 0.0))

        events.append({
            "pkey": pkey,
            "type": sstype,
            "cluster_id": cluster_id,
            "wt_pos": wt_pos,
            "mut_pos": mut_pos,
            "dpos": dpos,
            "wt_score": wt_score,
            "mut_score": mut_score,
            "dscore": dscore,
            "pct_delta": pct_delta,
            "distance_to_snv": dist,
            "rank_wt": rank_wt,
            "rank_mut": rank_mut,
            "conf_wt": conf_wt,
            "conf_mut": conf_mut,
            "conf_weighted_delta": conf_weighted_delta,
            "cls": cls,
            "is_high_impact": 0,  # fill later
            "priority": 0.0,      # fill later
            "in_radius": in_radius,
        })

    return pd.DataFrame(events)


# ---------------------------------------------------------------------------
# summarization with HARDENING
# ---------------------------------------------------------------------------

def _summarize_variant(events_df: pd.DataFrame,
                       sites_df: pd.DataFrame,
                       pkey: str,
                       report_radius: int,
                       policy: dict) -> dict:
    """
    Build the summary row for a single pkey.
    Hardened to return a valid row even if there are *no* sites/events.
    """
    visibility_threshold = policy["visibility_threshold"]
    high_cutoff = policy["high_cutoff"]
    shift_bp = policy["shift_bp"]
    distance_k = policy["distance_k"]

    if events_df is None or events_df.empty or "pkey" not in events_df.columns:
        return {
            "pkey": pkey,
            "n_sites_wt": 0,
            "n_sites_mut": 0,
            "n_clusters": 0,
            "global_count_gained_high": 0,
            "global_count_lost_high": 0,
            "global_count_shifted": 0,
            "global_max_abs_Δscore": 0.0,
            "global_sum_weighted_abs_Δ": 0.0,
            "nearest_event_bp_any": None,
            "local_count_gained_high": 0,
            "local_count_lost_high": 0,
            "local_count_shifted": 0,
            "local_max_abs_Δscore": 0.0,
            "nearest_event_bp_local": None,
            "frac_effect_in_radius": 0.0,
            "top_event_type": "none",
            "top_event_Δscore": 0.0,
            "top_event_pos": None,
            "dominant_boundary": None,
            "qc_flags": "no_sites",
        }

    sub_events = events_df[events_df["pkey"] == pkey]

    if sites_df is not None and not sites_df.empty and "pkey" in sites_df.columns:
        sub_sites = sites_df[sites_df["pkey"] == pkey]
    else:
        sub_sites = pd.DataFrame()

    if sub_events.empty and sub_sites.empty:
        return {
            "pkey": pkey,
            "n_sites_wt": 0,
            "n_sites_mut": 0,
            "n_clusters": 0,
            "global_count_gained_high": 0,
            "global_count_lost_high": 0,
            "global_count_shifted": 0,
            "global_max_abs_Δscore": 0.0,
            "global_sum_weighted_abs_Δ": 0.0,
            "nearest_event_bp_any": None,
            "local_count_gained_high": 0,
            "local_count_lost_high": 0,
            "local_count_shifted": 0,
            "local_max_abs_Δscore": 0.0,
            "nearest_event_bp_local": None,
            "frac_effect_in_radius": 0.0,
            "top_event_type": "none",
            "top_event_Δscore": 0.0,
            "top_event_pos": None,
            "dominant_boundary": None,
            "qc_flags": "no_sites",
        }

    # counts
    if not sub_sites.empty:
        n_sites_wt = len(sub_sites[(sub_sites["allele"] == "WT") & (sub_sites["visible_flag"] == 1)])
        n_sites_mut = len(sub_sites[(sub_sites["allele"] == "MUT") & (sub_sites["visible_flag"] == 1)])
    else:
        n_sites_wt = 0
        n_sites_mut = 0

    n_clusters = len(sub_events)

    # global metrics
    global_count_gained_high = len(sub_events[(sub_events["cls"] == "gained") &
                                              (sub_events[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    global_count_lost_high = len(sub_events[(sub_events["cls"] == "lost") &
                                            (sub_events[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    global_count_shifted = len(sub_events[sub_events["cls"] == "shifted"])

    # max abs delta
    if not sub_events["dscore"].isna().all():
        global_max_abs_Δscore = float(sub_events["dscore"].abs().max(skipna=True))
    else:
        global_max_abs_Δscore = 0.0

    # sum weighted abs delta
    swa = 0.0
    for _, ev in sub_events.iterrows():
        d = ev["distance_to_snv"]
        if pd.isna(d) or ev["dscore"] is None or pd.isna(ev["dscore"]):
            continue
        w = math.exp(- float(d) / float(distance_k)) if distance_k > 0 else 1.0
        swa += w * abs(float(ev["dscore"]))
    global_sum_weighted_abs_Δ = swa

    # nearest event
    if not sub_events["distance_to_snv"].isna().all():
        nearest_event_bp_any = int(sub_events["distance_to_snv"].min(skipna=True))
    else:
        nearest_event_bp_any = None

    # local metrics
    if report_radius is None:
        report_radius = 151  # fallback

    in_local = sub_events[(sub_events["distance_to_snv"] <= report_radius)]
    local_count_gained_high = len(in_local[(in_local["cls"] == "gained") &
                                           (in_local[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    local_count_lost_high = len(in_local[(in_local["cls"] == "lost") &
                                         (in_local[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    local_count_shifted = len(in_local[in_local["cls"] == "shifted"])
    if not in_local.empty and not in_local["dscore"].isna().all():
        local_max_abs_Δscore = float(in_local["dscore"].abs().max(skipna=True))
        nearest_event_bp_local = int(in_local["distance_to_snv"].min(skipna=True))
    else:
        local_max_abs_Δscore = 0.0
        nearest_event_bp_local = None

    # frac_effect_in_radius
    num = 0.0
    den = 0.0
    for _, ev in sub_events.iterrows():
        ds = ev["dscore"]
        if ds is None or pd.isna(ds):
            continue
        den += abs(float(ds))
        if ev["distance_to_snv"] <= report_radius:
            num += abs(float(ds))
    frac_effect_in_radius = (num / den) if den > 0 else 0.0

    # top event (by priority)
    if "priority" in sub_events.columns and not sub_events["priority"].isna().all():
        top_ev = sub_events.sort_values("priority", ascending=False).iloc[0]
        top_event_type = top_ev["cls"]
        top_event_Δscore = float(top_ev["dscore"]) if not pd.isna(top_ev["dscore"]) else 0.0
        top_event_pos = top_ev["mut_pos"] if pd.notna(top_ev["mut_pos"]) else top_ev["wt_pos"]
    else:
        top_event_type = "none"
        top_event_Δscore = 0.0
        top_event_pos = None

    # dominant boundary
    dom = None
    if not sub_events.empty:
        # sum |dscore| per type
        agg = (
            sub_events
            .dropna(subset=["dscore"])
            .assign(absd=lambda x: x["dscore"].abs())
            .groupby("type")["absd"].sum()
        )
        if not agg.empty:
            dom = agg.idxmax()

    # qc_flags
    flags = []
    if n_clusters == 0:
        flags.append("no_sites")
    if nearest_event_bp_any is not None and nearest_event_bp_any > 2000:
        flags.append("far_event>2kb")
    if global_max_abs_Δscore < 1.0:
        flags.append("low_signal_only")
    qc_flags = ";".join(flags) if flags else ""

    return {
        "pkey": pkey,
        "n_sites_wt": n_sites_wt,
        "n_sites_mut": n_sites_mut,
        "n_clusters": n_clusters,
        "global_count_gained_high": global_count_gained_high,
        "global_count_lost_high": global_count_lost_high,
        "global_count_shifted": global_count_shifted,
        "global_max_abs_Δscore": global_max_abs_Δscore,
        "global_sum_weighted_abs_Δ": global_sum_weighted_abs_Δ,
        "nearest_event_bp_any": nearest_event_bp_any,
        "local_count_gained_high": local_count_gained_high,
        "local_count_lost_high": local_count_lost_high,
        "local_count_shifted": local_count_shifted,
        "local_max_abs_Δscore": local_max_abs_Δscore,
        "nearest_event_bp_local": nearest_event_bp_local,
        "frac_effect_in_radius": frac_effect_in_radius,
        "top_event_type": top_event_type,
        "top_event_Δscore": top_event_Δscore,
        "top_event_pos": top_event_pos,
        "dominant_boundary": dom,
        "qc_flags": qc_flags,
    }


# ---------------------------------------------------------------------------
# main per-gene worker
# ---------------------------------------------------------------------------

def _process_gene(fasta_path: Path,
                  gene_mapping_df: pd.DataFrame,
                  genesplicer_dir: str,
                  pipeline: str,
                  window: int,
                  report_radius: int,
                  visibility_threshold: float,
                  high_cutoff: float,
                  shift_bp: int,
                  distance_k: int,
                  failure_map: dict):
    """
    Per-gene worker:
    - read WT genomic
    - run genesplicer on WT once
    - loop over each mutation row:
        - skip if in failure_map
        - synthesize ALT
        - run genesplicer
        - build sites
    - return events and sites
    """
    gene_name = extract_gene_from_filename(str(fasta_path)) or fasta_path.stem
    stats = {
        "gene": gene_name,
        "total_rows": int(gene_mapping_df.shape[0]) if gene_mapping_df is not None else 0,
        "processed": 0,
        "skipped_validation": 0,
        "skipped_invalid": 0,
        "skipped_out_of_range": 0,
        "skipped_runtime": 0,
    }
    fa = read_fasta(str(fasta_path))
    if not fa:
        return [], [], [], stats

    # adopt first sequence or key "genomic"
    if "genomic" in fa:
        wt_seq = fa["genomic"]
    else:
        wt_seq = next(iter(fa.values()))
    wt_len = len(wt_seq)

    # run WT once (full)
    wt_df = _run_genesplicer_on_seq(f"{gene_name}_WT", wt_seq, genesplicer_dir)

    # if mapping df is empty or missing, nothing to do
    if gene_mapping_df is None or gene_mapping_df.empty:
        return [], [], [], stats

    # harmonize column names
    cols = {c.lower(): c for c in gene_mapping_df.columns}
    # expect 'mutant' and 'genomic'
    mutant_col = cols.get("mutant") or cols.get("mutation") or cols.get("mut")
    genomic_col = cols.get("genomic") or cols.get("genomic_nt") or cols.get("genomic_mut") or cols.get("genomicmutation")

    if not mutant_col or not genomic_col:
        # skip gene
        return [], [], [], stats

    events_all = []
    sites_all = []
    pkeys = []

    # allow fallback for report radius
    if report_radius is None:
        report_radius = window

    for _, row in gene_mapping_df.iterrows():
        mutant_tok = str(row[mutant_col]).strip()
        genomic_tok = str(row[genomic_col]).strip()

        pkey = f"{gene_name}-{mutant_tok}"

        # skip by validation logs if needed
        if should_skip_mutation(gene_name, mutant_tok, failure_map):
            stats["skipped_validation"] += 1
            continue

        # parse genomic token e.g. T57261C
        try:
            snv_pos, (ref_nt, alt_nt) = get_mutation_data(genomic_tok)
        except Exception:
            # malformed, skip
            stats["skipped_invalid"] += 1
            continue

        # build ALT sequence in memory
        if snv_pos < 1 or snv_pos > wt_len:
            # out of range, skip
            stats["skipped_out_of_range"] += 1
            continue
        alt_seq = update_str(wt_seq, alt_nt, snv_pos)

        # run GeneSplicer on ALT
        try:
            mut_df = _run_genesplicer_on_seq(f"{gene_name}_{mutant_tok}", alt_seq, genesplicer_dir)
        except Exception:
            stats["skipped_runtime"] += 1
            continue

        stats["processed"] += 1
        pkeys.append(pkey)

        # convert WT ALT to canonical sites
        wt_sites = _build_sites_for_allele(
            pkey=pkey,
            allele="WT",
            snv_pos=snv_pos,
            df=wt_df,
            visibility_threshold=visibility_threshold,
            report_radius=report_radius,
            scan_mode=pipeline,
        )

        mut_sites = _build_sites_for_allele(
            pkey=pkey,
            allele="MUT",
            snv_pos=snv_pos,
            df=mut_df,
            visibility_threshold=visibility_threshold,
            report_radius=report_radius,
            scan_mode=pipeline,
        )

        # merge sites
        sites_concat = pd.concat([wt_sites, mut_sites], ignore_index=True)
        # cluster
        sites_clustered = _cluster_sites(sites_concat, cluster_radius=3)
        # pair
        events_df = _pair_clusters_to_events(sites_clustered,
                                             visibility_threshold=visibility_threshold,
                                             shift_bp=shift_bp)

        events_all.append(events_df)
        sites_all.append(sites_clustered)

    return events_all, sites_all, pkeys, stats


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="GeneSplicer WT↔ALT ensemble delta caller")
    parser.add_argument("-i", "--input", required=True, help="Directory of genomic FASTA files")
    parser.add_argument("-m", "--mapping-dir", required=True, help="Directory of genomic mutation/mapping files")
    parser.add_argument("-g", "--genesplicer-dir", required=True, help="Directory containing GeneSplicer binary")
    parser.add_argument("-o", "--output", required=True, help="Output summary TSV basename (e.g. genesplicer.summary.tsv)")
    parser.add_argument("--pipeline", choices=["full", "window", "custom"], default="full",
                        help="Pipeline mode. 'full' = run on whole genomic sequence once and mutate in-memory.")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW, help="Window/reporting size (default 151)")
    parser.add_argument("--report-radius", type=int, default=DEFAULT_REPORT_RADIUS,
                        help="Radius (bp) to count as local; default = --window")
    parser.add_argument("--distance-k", type=int, default=DEFAULT_DISTANCE_K,
                        help="Distance-decay kernel k (default 75)")
    parser.add_argument("--visibility-threshold", type=float, default=DEFAULT_VISIBILITY_THRESHOLD,
                        help="Minimum score to treat a site as visible (default 1.0)")
    parser.add_argument("--high-cutoff", type=float, default=DEFAULT_HIGH_CUTOFF,
                        help="High-confidence score cutoff (default 5.0)")
    parser.add_argument("--shift-bp", type=int, default=DEFAULT_SHIFT_BP,
                        help="Minimum bp difference to classify as a shifted site (default 3)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Max parallel workers (default: half cores, capped at 8)")
    parser.add_argument("--log", help="Validation log file/dir to skip failed mutations")
    args = parser.parse_args()

    input_dir = Path(args.input)
    mapping_dir = args.mapping_dir
    genesplicer_dir = args.genesplicer_dir
    output_base = args.output
    pipeline = args.pipeline
    window = args.window
    report_radius = args.report_radius
    distance_k = args.distance_k
    visibility_threshold = args.visibility_threshold
    high_cutoff = args.high_cutoff
    shift_bp = args.shift_bp

    if report_radius is None:
        report_radius = window

    failure_map = load_validation_failures(args.log) if args.log else {}

    # load mappings
    mappings = load_mapping_dir(mapping_dir)

    # discover FASTA files
    fasta_paths = [p for p in input_dir.iterdir() if p.is_file() and p.suffix in (".fa", ".fasta", ".fna")]
    if not fasta_paths:
        print(f"No FASTA files found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    # policy bundle
    policy = {
        "visibility_threshold": visibility_threshold,
        "high_cutoff": high_cutoff,
        "shift_bp": shift_bp,
        "distance_k": distance_k,
    }

    # parallel per-gene
    if args.workers is None:
        n_cpu = os.cpu_count() or multiprocessing.cpu_count() or 1
        max_workers = max(1, min(n_cpu // 2 if n_cpu > 1 else 1, DEFAULT_MAX_WORKERS))
    else:
        max_workers = args.workers

    all_events = []
    all_sites = []
    all_pkeys = []

    total_genes = len(fasta_paths)
    genes_completed = 0
    total_variants_seen = 0
    total_variants_processed = 0
    total_skipped_validation = 0
    total_skipped_invalid = 0
    total_skipped_out_of_range = 0
    total_skipped_runtime = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        futs = []
        fut_meta = {}
        for fasta_path in fasta_paths:
            gene_name = extract_gene_from_filename(str(fasta_path)) or fasta_path.stem
            gene_map_df = mappings.get(gene_name, pd.DataFrame())
            fut = ex.submit(
                _process_gene,
                fasta_path,
                gene_map_df,
                genesplicer_dir,
                pipeline,
                window,
                report_radius,
                visibility_threshold,
                high_cutoff,
                shift_bp,
                distance_k,
                failure_map,
            )
            futs.append(fut)
            fut_meta[fut] = gene_name

        for fut in concurrent.futures.as_completed(futs):
            gene_events, gene_sites, gene_pkeys, gene_stats = fut.result()
            all_events.extend(gene_events)
            all_sites.extend(gene_sites)
            all_pkeys.extend(gene_pkeys)

            genes_completed += 1
            total_variants_seen += gene_stats.get("total_rows", 0)
            total_variants_processed += gene_stats.get("processed", 0)
            total_skipped_validation += gene_stats.get("skipped_validation", 0)
            total_skipped_invalid += gene_stats.get("skipped_invalid", 0)
            total_skipped_out_of_range += gene_stats.get("skipped_out_of_range", 0)
            total_skipped_runtime += gene_stats.get("skipped_runtime", 0)

            gene_name = gene_stats.get("gene") or fut_meta.get(fut, "unknown")
            skips = []
            if gene_stats.get("skipped_validation", 0):
                skips.append(f"validation={gene_stats['skipped_validation']}")
            if gene_stats.get("skipped_invalid", 0):
                skips.append(f"invalid={gene_stats['skipped_invalid']}")
            if gene_stats.get("skipped_out_of_range", 0):
                skips.append(f"out_of_range={gene_stats['skipped_out_of_range']}")
            if gene_stats.get("skipped_runtime", 0):
                skips.append(f"runtime={gene_stats['skipped_runtime']}")
            skip_msg = ", ".join(skips) if skips else "none"
            print(
                f"[{genes_completed}/{total_genes}] {gene_name}: "
                f"processed {gene_stats.get('processed', 0)}/{gene_stats.get('total_rows', 0)} variants "
                f"(skipped: {skip_msg}). "
                f"Total processed so far: {total_variants_processed}"
            )
            sys.stdout.flush()

    # concat events/sites with hardening
    if all_events:
        events_df = pd.concat(all_events, ignore_index=True)
    else:
        events_df = pd.DataFrame()

    if all_sites:
        sites_df = pd.concat(all_sites, ignore_index=True)
    else:
        sites_df = pd.DataFrame()

    # ensure columns exist
    events_required = [
        "pkey", "type", "cluster_id",
        "wt_pos", "mut_pos", "dpos",
        "wt_score", "mut_score", "dscore", "pct_delta",
        "distance_to_snv",
        "rank_wt", "rank_mut",
        "conf_wt", "conf_mut", "conf_weighted_delta",
        "cls", "is_high_impact", "priority", "in_radius",
    ]
    if events_df.empty:
        events_df = pd.DataFrame(columns=events_required)
    else:
        for c in events_required:
            if c not in events_df.columns:
                events_df[c] = pd.NA

    sites_required = [
        "pkey", "allele", "type", "site_pos", "score", "confidence",
        "rank", "distance_to_snv", "visible_flag", "cluster_id", "in_radius",
    ]
    if sites_df.empty:
        sites_df = pd.DataFrame(columns=sites_required)
    else:
        for c in sites_required:
            if c not in sites_df.columns:
                sites_df[c] = pd.NA

    # set event-level priority now that everything exists
    def _calc_priority(row):
        if pd.isna(row["dscore"]) or row["dscore"] is None:
            base = 0.0
        else:
            base = abs(float(row["dscore"]))
        d = row["distance_to_snv"]
        if pd.isna(d) or d is None:
            w = 1.0
        else:
            w = math.exp(- float(d) / float(distance_k)) if distance_k > 0 else 1.0
        bonus = 0.0
        if row["cls"] in ("gained", "lost"):
            # check high cutoff
            m = max(
                float(row["wt_score"]) if not pd.isna(row["wt_score"]) else 0.0,
                float(row["mut_score"]) if not pd.isna(row["mut_score"]) else 0.0,
            )
            if m >= high_cutoff:
                bonus += 2.0
        if row["cls"] == "shifted" and row["dpos"] is not None and not pd.isna(row["dpos"]) and abs(int(row["dpos"])) >= shift_bp:
            bonus += 1.0
        return base * w + bonus

    if not events_df.empty:
        events_df["priority"] = events_df.apply(_calc_priority, axis=1)
        # set is_high_impact
        def _is_hi(row):
            if pd.isna(row["dscore"]):
                return 0
            if abs(float(row["dscore"])) >= 5.0:
                return 1
            if row["cls"] in ("gained", "lost"):
                m = max(
                    float(row["wt_score"]) if not pd.isna(row["wt_score"]) else 0.0,
                    float(row["mut_score"]) if not pd.isna(row["mut_score"]) else 0.0,
                )
                if m >= high_cutoff:
                    return 1
            return 0
        events_df["is_high_impact"] = events_df.apply(_is_hi, axis=1)

    # summarize per pkey
    summary_rows = []
    # unique pkeys from run; if empty, derive from events_df
    if not all_pkeys:
        all_pkeys = list(events_df["pkey"].unique()) if "pkey" in events_df.columns else []
    for pkey in all_pkeys:
        summary_rows.append(
            _summarize_variant(
                events_df=events_df,
                sites_df=sites_df,
                pkey=pkey,
                report_radius=report_radius,
                policy=policy,
            )
        )
    summary_df = pd.DataFrame(summary_rows)

    # write outputs
    summary_path = Path(output_base)
    events_path = summary_path.with_suffix(".events.tsv")
    sites_path = summary_path.with_suffix(".sites.tsv")

    summary_df.to_csv(summary_path, sep="\t", index=False)
    events_df.to_csv(events_path, sep="\t", index=False)
    sites_df.to_csv(sites_path, sep="\t", index=False)

    print(f"wrote summary to {summary_path}")
    print(f"wrote events to  {events_path}")
    print(f"wrote sites to   {sites_path}")
    print(
        "Run summary: "
        f"{genes_completed}/{total_genes} genes processed, "
        f"{total_variants_processed}/{total_variants_seen} variants analysed; "
        f"skipped validation={total_skipped_validation}, "
        f"invalid={total_skipped_invalid}, "
        f"out_of_range={total_skipped_out_of_range}, "
        f"runtime={total_skipped_runtime}."
    )


if __name__ == "__main__":
    main()
