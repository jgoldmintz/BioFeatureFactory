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
GeneSplicer WT<->ALT ensemble delta caller
- Runs GeneSplicer on full genomic sequences for each gene.

GeneSplicer is run once per allele and reports each site in the coordinate frame
of the sequence it was handed. For an SNV those two frames coincide, so the WT
and MUT site lists can be joined on raw position. For anything that changes
length they do not: every MUT site 3' of the edit sits at wt_pos + length_delta.
Joining on raw position there reports a phantom displacement for every
downstream site -- large enough to be called `shifted`, or to split into a
`lost`/`gained` pair that the redirect pass then links as one relocated site.

So every cross-allele comparison in this module is done in the WT frame: each
site carries `ref_frame_pos` (its WT-frame coordinate) alongside `site_pos` (its
real coordinate in its own allele), and clustering, dpos and the redirect pass
all key on the former. `frame_status` names the three cases -- aligned, a WT site
the edit deleted, and a MUT site sitting on inserted sequence that has no WT
coordinate at all.
"""

import os
import sys
import argparse
import json
import shutil
import tempfile
import subprocess
from pathlib import Path
import pandas as pd
import math
import multiprocessing
import concurrent.futures

# ---------------------------------------------------------------------------
# imports from biofeaturefactory.utils.utility
# ---------------------------------------------------------------------------
from biofeaturefactory.utils.utility import (
    read_fasta,
    parse_variant,
    splice_seq,
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

def _resolve_model_dir(model_dir_arg: str | None, genesplicer_dir: str | None) -> str:
    """Return an absolute GeneSplicer model directory (one containing config_file).

    Precedence:
      1. explicit --model-dir;
      2. source-build layout: {genesplicer_dir}/../human (binary in .../sources);
      3. PATH/conda layout: probed relative to the resolved `genesplicer` binary —
         {bindir}/human (conda share layout, model beside the binary) or
         {bindir}/../human (source layout).
    Exits with a clear message if none contains config_file, rather than letting
    a missing model surface as silent empty GeneSplicer output.
    """
    candidates = []
    if model_dir_arg:
        candidates.append(model_dir_arg)
    elif genesplicer_dir:
        candidates.append(os.path.join(genesplicer_dir, "..", "human"))
    else:
        gs = shutil.which("genesplicer")
        if gs:
            bindir = os.path.dirname(os.path.realpath(gs))
            candidates.append(os.path.join(bindir, "human"))         # conda share layout
            candidates.append(os.path.join(bindir, "..", "human"))   # source-build layout
    tried = []
    for c in candidates:
        c_abs = os.path.abspath(c)
        tried.append(c_abs)
        if os.path.isfile(os.path.join(c_abs, "config_file")):
            return c_abs
    sys.exit(f"GeneSplicer human model not found (need a directory containing config_file). "
             f"Tried: {tried or 'none'}. Pass --model-dir with an absolute path.")


def _run_genesplicer_on_seq(seq_name: str, seq: str, genesplicer_dir: str | None, model_dir: str) -> pd.DataFrame:
    """
    Run GeneSplicer on a single in-memory sequence.
    Returns a DataFrame with columns:
        End5, End3, Score, confidence, splice_site_type

    genesplicer_dir may be None to use genesplicer from PATH.
    """
    if genesplicer_dir:
        os.makedirs(genesplicer_dir, exist_ok=True)
        tmp_dir = genesplicer_dir
    else:
        tmp_dir = tempfile.gettempdir()

    with tempfile.NamedTemporaryFile(mode="w", dir=tmp_dir, delete=False, suffix=".fasta") as tmpf:
        tmp_name = os.path.basename(tmpf.name)
        # NO trailing newline, deliberately. GeneSplicer's reader
        # (sources/genesplicer.cpp:376-385) strips the newline, decrements the
        # length, and then decrements AGAIN -- `length-=1` plus
        # strncat(Data,line,lcopy-1) -- so it silently drops the last base of
        # every newline-terminated line. Measured on a 2096 nt sequence:
        # with a trailing newline it reports "Done 2095bp", without one
        # "Done 2096bp". Writing the sequence as ONE line already keeps the loss
        # to a single base instead of one per wrapped line (a 60-column FASTA
        # loses ~1.7% of its bases at scattered positions and shifts every
        # downstream coordinate); dropping the final newline removes it entirely.
        tmpf.write(f">{seq_name}\n{seq}")
        tmp_path = tmpf.name

    out = None
    rc = None
    err = ""
    cmd = ""
    try:
        if genesplicer_dir:
            cmd = f"cd {genesplicer_dir} && ./genesplicer {tmp_name} {model_dir}"
        else:
            cmd = f"genesplicer {tmp_path} {model_dir}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        out = result.stdout.strip()
        rc = result.returncode
        err = (result.stderr or "").strip()
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    # F28: a nonzero exit returns empty stdout that is otherwise indistinguishable
    # from "no splice sites". Raise so the caller records a failure instead of
    # silently classifying every counterpart allele's sites as gained/lost.
    if rc != 0:
        raise RuntimeError(f"GeneSplicer failed for {seq_name} (rc={rc}): {err or 'no stderr'}")

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
    Accepts a single CSV file or a directory of CSV files.
    Expected columns include at least: 'mutant' and 'genomic'.
    """
    mapping_path = Path(mapping_dir)
    if not mapping_path.exists():
        raise FileNotFoundError(f"Mapping path not found: {mapping_path}")
    maps = {}
    files = [mapping_path] if mapping_path.is_file() else mapping_path.rglob("*.csv")
    for f in files:
        if not f.is_file() or f.suffix not in (".tsv", ".csv", ".txt"):
            continue
        df = pd.read_csv(f, sep="\t" if f.suffix == ".tsv" else ",")
        gene = extract_gene_from_filename(str(f))
        if not gene:
            gene = f.stem
        maps[gene] = df
    return maps


# ---------------------------------------------------------------------------
# site builders
# ---------------------------------------------------------------------------

SITE_COLUMNS = [
    "pkey", "allele", "type", "site_pos", "ref_frame_pos", "frame_status",
    "score", "confidence", "rank", "distance_to_snv", "visible_flag",
    "cluster_id", "in_radius",
]


def _wt_frame_position(site_pos: int, allele: str, variant_pos1: int,
                       ref_len: int, alt_len: int) -> tuple[int, str]:
    """Project a 1-based site coordinate into the WT frame.

    Returns (ref_frame_pos, frame_status). This is the exact inverse of
    utility.align_wt_to_mut, region for region, so the two can never disagree
    about which MUT base corresponds to which WT base:

        before the edit   identity
        inside the REF    the first min(ref_len, alt_len) bases pair up
        after the edit    shifted by (alt_len - ref_len)

    frame_status:
        aligned    the site has a counterpart coordinate in the other allele
        deleted    (WT only) the edit removed the base this site is keyed on --
                   align_wt_to_mut maps it to None
        inserted   (MUT only) the site sits on bases the ALT inserted, which have
                   no WT coordinate at all

    An inserted site is reported at the edit anchor, because under the SpliceAI
    convention align_wt_to_mut implements an insertion collapses to the edit
    boundary -- that is its honest WT-frame position, not a fabricated one. It is
    still tagged `inserted` so clustering keeps it out of the aligned namespace.

    That tag covers the ALT bases BEYOND min(ref_len, alt_len) and no others,
    because that is the extent align_wt_to_mut leaves without a WT index. Inside
    the first min(ref_len, alt_len) bases the two alleles are paired by
    coordinate whatever their content -- deliberately, since that is the
    measurement a delta caller exists to make ("did the site at this locus
    survive the edit?") and it is what already happens for the substituted base
    of an SNV. For a delins those paired ALT bases are novel sequence, so a site
    called on them WILL be joined to a WT site at the same coordinate and
    reported as one site weakened/strengthened rather than as a loss plus a gain.
    Measured: WT donor 760 (6.22) vs a donor built entirely from a 20 nt ALT,
    called at MUT 759, comes out `weakened` dscore -2.36. The pairing is correct
    under the convention and the scores are real; only the one-site reading of it
    is an inference the caller has to make for itself.

    Note that the anchor locates the site but does NOT measure displacement to it:
    every base of the insertion shares that one coordinate, so a difference taken
    against it says nothing about how far into the insertion the site sits. See
    _reconcile_redirects, which is why it refuses the register-shift test there.

    For an SNV (ref_len == alt_len == 1) every branch is the identity and every
    status is `aligned`, so nothing about the SNV path changes.
    """
    edit_off = variant_pos1 - 1     # 0-based WT index of the first REF base
    i = site_pos - 1                # 0-based index in this allele's own sequence

    if allele == "WT":
        # A WT site is already in the WT frame; only its survival is in question.
        if edit_off + alt_len <= i < edit_off + ref_len:
            return site_pos, "deleted"
        return site_pos, "aligned"

    if i < edit_off:
        return i + 1, "aligned"
    k = i - edit_off
    if k < alt_len:
        if k < ref_len:
            return edit_off + k + 1, "aligned"
        return variant_pos1, "inserted"
    return i - (alt_len - ref_len) + 1, "aligned"


def _build_sites_for_allele(pkey: str,
                            allele: str,
                            df: pd.DataFrame,
                            visibility_threshold: float,
                            report_radius: int,
                            variant_pos1: int,
                            ref_len: int,
                            alt_len: int) -> pd.DataFrame:
    """
    Convert GeneSplicer raw DF into the canonical per-site schema.

    variant_pos1 is 1-based at the first REF base (Variant.pos). ref_len/alt_len
    are the allele lengths; the span this allele actually carries is ref_len bases
    for WT and alt_len for MUT.

    `scan_mode` used to be accepted here and never read -- the only value ever
    passed was the literal "full" -- so it is gone along with the `pipeline`
    parameter that existed solely to supply it.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=SITE_COLUMNS)

    # The variant's own span, in THIS allele's coordinates. WT carries the REF,
    # MUT carries the ALT, and they differ in length whenever the edit does.
    span_start = variant_pos1
    span_end = variant_pos1 + (ref_len if allele == "WT" else alt_len) - 1

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

        # Distance to the variant's SPAN, measured in this allele's own frame:
        # zero inside the span, else to the nearer edge. A multi-base REF has no
        # single base to measure to, and measuring to its first base alone makes
        # the distance grow with len(REF) for sites lying 3' of it. For an SNV
        # both edges are the same base, so this is exactly the old
        # abs(site_pos - (snv_pos + 1)).
        if span_start <= site_pos <= span_end:
            dist = 0
        else:
            dist = min(abs(site_pos - span_start), abs(site_pos - span_end))
        in_radius = 0
        if report_radius is not None:
            in_radius = 1 if dist <= report_radius else 0

        ref_frame_pos, frame_status = _wt_frame_position(
            site_pos, allele, variant_pos1, ref_len, alt_len)

        rows.append({
            "pkey": pkey,
            "allele": allele,
            "type": sstype,
            "site_pos": site_pos,
            "ref_frame_pos": ref_frame_pos,
            "frame_status": frame_status,
            "score": score,
            "confidence": conf_w,
            "distance_to_snv": dist,
            "visible_flag": 1 if score >= visibility_threshold else 0,
            "cluster_id": None,  # filled later
            "in_radius": in_radius,
        })

    # rank within allele×type by score (F33). The old sort/reset/set_index chain
    # reattached score-ordered rank values onto position-ordered index labels, so
    # rank tracked row position, not score. pd.Series.rank does it directly;
    # method="first" breaks ties by order to keep ranks distinct 1..n.
    df_sites = pd.DataFrame(rows)
    if not df_sites.empty:
        for t in ("donor", "acceptor"):
            mask = df_sites["type"] == t
            df_sites.loc[mask, "rank"] = (
                df_sites.loc[mask, "score"]
                .rank(ascending=False, method="first")
                .astype(int)
            )
    else:
        df_sites["rank"] = pd.NA

    return df_sites


# ---------------------------------------------------------------------------
# clustering and event detection
# ---------------------------------------------------------------------------

def _cluster_sites(sites_df: pd.DataFrame, cluster_radius: int, max_cluster_span: int = 0) -> pd.DataFrame:
    """
    Assign cluster_id per (pkey, type, inserted-or-not) using single-linkage by
    position. Both alleles are clustered together and NOT split by allele -- that
    shared cluster_id is exactly what joins a WT site to its MUT counterpart in
    _pair_clusters_to_events.

    Linkage is on `ref_frame_pos`, NOT on `site_pos`. The cluster is the unit the
    WT and MUT alleles are joined on, so it has to be built in a frame both
    alleles share. Under an indel a MUT site 3' of the edit has
    site_pos == wt_site_pos + length_delta, so linking on site_pos compares a WT
    coordinate to a MUT coordinate and reads the frame offset as a displacement
    of the site. With the defaults that is not a subtle error: cluster_radius and
    --shift-bp are both 3, so a 3 nt indel puts every downstream site just close
    enough to co-cluster and just far enough to be classified `shifted`, and a
    larger one splits it into a lost/gained pair that _reconcile_redirects then
    links as a relocated site. In the WT frame an unaffected downstream site has
    dpos == 0 and no event at all, which is the truth.

    Sites on inserted sequence have no WT coordinate and are clustered in a
    separate namespace (cluster ids `di1`, `ai1`, ... rather than `d1`, `a1`), so
    they surface as genuinely gained sites instead of being welded onto whichever
    WT site the edit anchor happens to sit next to.

    max_cluster_span is accepted but INERT (0 = unbounded), and must stay that way
    until a non-fabricating bound exists. Single linkage is unbounded — 101 sites
    3 bp apart at radius 3 weld into one 300 bp cluster — but a left-anchored fixed
    window is not the fix: it cuts at an absolute boundary rather than at a gap, so
    a genuine 1 bp WT/MUT register pair straddling the boundary is split into
    separate clusters. The orphaned WT site then pairs with a bystander present in
    both alleles and a high-magnitude 'shifted'/'gained' event is fabricated
    (measured: bystanders at 100/103/106/109, WT 112 -> MUT 113, span 12).
    A correct bound must cut at the largest gap and never separate a cross-allele
    pair. See PIPELINE_ERROR_AUDIT.md finding 30(c).
    """
    if sites_df is None or sites_df.empty:
        return sites_df

    sites_df = sites_df.copy()
    # Inserted sites cluster among themselves. Carrying that as a sort/group key
    # keeps this a single pass, so the ordering of the aligned rows -- and every
    # cluster id built from it -- is exactly what it was before.
    sites_df["_ins"] = (sites_df["frame_status"] == "inserted").astype(int)
    sites_df = sites_df.sort_values(
        ["pkey", "type", "_ins", "ref_frame_pos", "allele"]).reset_index(drop=True)

    cluster_ids = []
    current_key = None
    current_cluster = 0
    last_pos = None
    cluster_start = None

    for _, row in sites_df.iterrows():
        key = (row["pkey"], row["type"], int(row["_ins"]))
        pos = int(row["ref_frame_pos"])
        prefix = f"{row['type'][0]}i" if row["_ins"] else row["type"][0]
        if key != current_key:
            current_key = key
            current_cluster = 1
            last_pos = pos
            cluster_start = pos
            cluster_ids.append(f"{prefix}{current_cluster}")
            continue

        # same key
        if abs(pos - last_pos) <= cluster_radius and (
                not max_cluster_span or abs(pos - cluster_start) <= max_cluster_span):
            # same cluster
            cluster_ids.append(f"{prefix}{current_cluster}")
        else:
            current_cluster += 1
            cluster_start = pos
            cluster_ids.append(f"{prefix}{current_cluster}")
        last_pos = pos

    sites_df["cluster_id"] = cluster_ids
    return sites_df.drop(columns=["_ins"])


EVENT_BASE_COLUMNS = [
    "pkey", "type", "cluster_id",
    "wt_pos", "mut_pos", "dpos", "dpos_raw",
    "wt_frame_pos", "mut_frame_pos", "wt_frame_status", "mut_frame_status",
    "wt_score", "mut_score", "dscore", "pct_delta",
    "distance_to_snv",
    "rank_wt", "rank_mut",
    "conf_wt", "conf_mut", "conf_weighted_delta",
    "cls", "is_high_impact", "priority", "in_radius",
]


def _pair_clusters_to_events(sites_df: pd.DataFrame,
                             visibility_threshold: float,
                             shift_bp: int) -> tuple[pd.DataFrame, int]:
    """
    From per-site rows, produce per-cluster per-pkey events.
    For each (pkey, type, cluster_id):
        - take top WT (by score) and top MUT
        - compute deltas
        - classify
    Returns the events frame and the number of non-top cluster members dropped
    by the top-of-cluster selection.

    `wt_pos`/`mut_pos` stay in their own allele's coordinates -- that is where the
    site really is, and it is what a consumer needs to look the site up. `dpos` is
    computed from the WT-frame positions instead, because it answers "did this
    site move?", and under an indel the raw difference answers "how long was the
    indel?" for every site 3' of it. `dpos_raw` keeps the raw difference so
    nothing is lost; the two are equal for every length-preserving variant.
    """
    if sites_df is None or sites_df.empty:
        return pd.DataFrame(columns=EVENT_BASE_COLUMNS), 0

    events = []
    n_discarded = 0
    grouped = sites_df.groupby(["pkey", "type", "cluster_id"], dropna=False)

    for (pkey, sstype, cluster_id), grp in grouped:
        # split by allele
        wt_grp = grp[grp["allele"] == "WT"].sort_values("score", ascending=False)
        mut_grp = grp[grp["allele"] == "MUT"].sort_values("score", ascending=False)

        wt_site = wt_grp.iloc[0] if not wt_grp.empty else None
        mut_site = mut_grp.iloc[0] if not mut_grp.empty else None
        n_discarded += max(len(wt_grp) - 1, 0) + max(len(mut_grp) - 1, 0)

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

        # WT-frame positions: the only frame in which the two alleles' coordinates
        # are comparable. Equal to site_pos for every site of a length-preserving
        # variant, so dpos is unchanged for SNVs and MNVs.
        wt_frame_pos = int(wt_site["ref_frame_pos"]) if wt_site is not None else None
        mut_frame_pos = int(mut_site["ref_frame_pos"]) if mut_site is not None else None
        wt_frame_status = str(wt_site["frame_status"]) if wt_site is not None else None
        mut_frame_status = str(mut_site["frame_status"]) if mut_site is not None else None

        dpos = None
        if wt_frame_pos is not None and mut_frame_pos is not None:
            dpos = mut_frame_pos - wt_frame_pos
        dpos_raw = None
        if wt_pos is not None and mut_pos is not None:
            dpos_raw = mut_pos - wt_pos

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
            "dpos_raw": dpos_raw,
            "wt_frame_pos": wt_frame_pos,
            "mut_frame_pos": mut_frame_pos,
            "wt_frame_status": wt_frame_status,
            "mut_frame_status": mut_frame_status,
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

    return pd.DataFrame(events), n_discarded


# ---------------------------------------------------------------------------
# summarization with HARDENING
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# redirect reconciliation
# ---------------------------------------------------------------------------

DEFAULT_REDIRECT_MAX_BP = 100   # link a lost/gained pair within this distance
DEFAULT_REGISTER_MAX_BP = 12    # <= this is a register shift, not a new site
DEFAULT_REDIRECT_SNV_BP = 80    # variant must be this close to the lost site
DEFAULT_CLUSTER_RADIUS = 3
DEFAULT_MAX_CLUSTER_SPAN = 0    # 0 = unbounded; a fixed window splits real pairs (see _cluster_sites)


def _reconcile_redirects(events_df: pd.DataFrame,
                         max_bp: int = DEFAULT_REDIRECT_MAX_BP,
                         register_bp: int = DEFAULT_REGISTER_MAX_BP,
                         snv_max_bp: int = DEFAULT_REDIRECT_SNV_BP) -> pd.DataFrame:
    """Link `lost` + `gained` event pairs that represent one redirected site.

    A single physical redirection whose displacement exceeds the cluster radius
    is emitted by _pair_clusters_to_events as two independent events, which
    double-counts it in every summary metric. This pass identifies those pairs
    and tags them with a shared redirect_id so the summariser can count them
    once.

    Rows are NEVER merged and dscore is NEVER rewritten. `cls` is preserved.
    The relabel lives entirely in the added columns, so it is reversible and
    auditable, and a downstream consumer can re-threshold on redirect_dpos or
    score_transfer_ratio without re-running GeneSplicer.

    Matching is type-matched, one-to-one, nearest-first (ties broken by the
    stronger gained site, then by position for determinism).

    Displacement is measured in the WT frame (mut_frame_pos - wt_frame_pos), for
    the same reason clustering is: on raw coordinates an indel makes every
    downstream lost/gained pair look like a relocation by exactly length_delta,
    so a 20 nt deletion would report every site 3' of it as a redirected site.
    A gained site sitting on inserted sequence is still linked -- it is one
    physical redirection and must be counted once -- but it is always classed
    `redirected` and carries no redirect_dpos; see the comment at the assignment
    for why a displacement measured across an insertion is not a register shift.

    Adds columns:
        redirect_id           shared id for the two rows of a linked pair
        redirect_dpos         mut_frame_pos - wt_frame_pos, signed (WT frame);
                              EMPTY when the gained site is on inserted sequence
        score_transfer_ratio  mut_score(gain) / wt_score(loss); diagnostic only
        redirect_class        'shifted' (<= register_bp) or 'redirected'
        cls_reconciled        e.g. 'lost_redirected', 'gained_shifted'
    """
    if events_df is None or events_df.empty:
        return events_df

    ev = events_df.copy()
    for c in ("redirect_id", "redirect_dpos", "score_transfer_ratio", "redirect_class"):
        if c not in ev.columns:
            ev[c] = pd.NA
    ev["cls_reconciled"] = ev["cls"]

    for (_pkey, sstype), grp in ev.groupby(["pkey", "type"], dropna=False):
        losses = grp[(grp["cls"] == "lost") & grp["wt_frame_pos"].notna()]
        gains = grp[(grp["cls"] == "gained") & grp["mut_frame_pos"].notna()]
        if losses.empty or gains.empty:
            continue

        cands = []
        for li, lrow in losses.iterrows():
            d_snv = lrow["distance_to_snv"]
            if pd.notna(d_snv) and float(d_snv) > snv_max_bp:
                continue  # variant too far away to have plausibly caused the loss
            wt_pos = int(lrow["wt_pos"])
            wt_frame = int(lrow["wt_frame_pos"])
            wt_s = float(lrow["wt_score"]) if pd.notna(lrow["wt_score"]) else 0.0
            for gi, grow in gains.iterrows():
                mut_pos = int(grow["mut_pos"])
                dpos = int(grow["mut_frame_pos"]) - wt_frame
                if abs(dpos) > max_bp:
                    continue
                mut_s = float(grow["mut_score"]) if pd.notna(grow["mut_score"]) else 0.0
                ratio = (mut_s / wt_s) if wt_s > 1e-6 else float("nan")
                gained_inserted = str(grow["mut_frame_status"]) == "inserted"
                # sort key: nearest, then strongest gain, then position (determinism)
                cands.append((abs(dpos), -mut_s, wt_pos, mut_pos, li, gi, dpos, ratio,
                              gained_inserted))

        cands.sort(key=lambda t: t[:4])
        used_l, used_g = set(), set()
        for _absd, _negs, wt_pos, mut_pos, li, gi, dpos, ratio, gained_inserted in cands:
            if li in used_l or gi in used_g:
                continue
            used_l.add(li)
            used_g.add(gi)
            # deterministic id: stable between the global run and the per-gene run
            tag = f"{sstype[0]}:{wt_pos}>{mut_pos}"
            if gained_inserted:
                # The gained site is on inserted sequence, which has no WT extent:
                # every base of the insertion collapses to the edit boundary, so
                # the WT-frame displacement is the distance from the lost site to
                # the edit and says nothing about how far into the insertion the
                # new site actually sits. Measured on it, a site 14 bp inside a
                # 90 bp insertion comes out as a 3 bp "register shift" -- the one
                # thing it certainly is not. So the register test is skipped (a
                # site in novel sequence is a redirection by construction) and
                # redirect_dpos is left EMPTY rather than reporting a number that
                # re-thresholds, per this function's own docstring, straight back
                # into the wrong class. `mut_frame_status` on the row names the
                # reason.
                rclass = "redirected"
                dpos_out = pd.NA
            else:
                rclass = "shifted" if abs(dpos) <= register_bp else "redirected"
                dpos_out = dpos
            for ix in (li, gi):
                ev.at[ix, "redirect_id"] = tag
                ev.at[ix, "redirect_dpos"] = dpos_out
                ev.at[ix, "score_transfer_ratio"] = ratio
                ev.at[ix, "redirect_class"] = rclass
            ev.at[li, "cls_reconciled"] = f"lost_{rclass}"
            ev.at[gi, "cls_reconciled"] = f"gained_{rclass}"

    return ev


def _empty_summary_row(pkey: str, variant_meta: dict) -> dict:
    """Summary row for a variant GeneSplicer found nothing to compare.

    Written twice verbatim before this existed, which meant any column added to
    the schema had to be added to both copies or the summary frame's columns
    would depend on which branch happened to fire.

    The zeros here are real measurements -- GeneSplicer ran on both alleles and
    reported no sites -- and `qc_flags` says so by name. The variant's own
    identity is still known and is still reported: which class it was and how the
    length changed do not depend on any site being found.
    """
    delta = variant_meta["length_delta"]
    return {
        "pkey": pkey,
        "n_sites_wt": 0,
        "n_sites_mut": 0,
        "n_clusters": 0,
        "global_count_gained_high": 0,
        "global_count_lost_high": 0,
        "global_count_shifted": 0,
        "global_count_redirected": 0,
        "global_max_abs_deltascore": 0.0,
        "global_sum_weighted_abs_delta": 0.0,
        "nearest_event_bp_any": None,
        "local_count_gained_high": 0,
        "local_count_lost_high": 0,
        "local_count_shifted": 0,
        "local_count_redirected": 0,
        "local_max_abs_deltascore": 0.0,
        "nearest_event_bp_local": None,
        "frac_effect_in_radius": 0.0,
        "top_event_type": "none",
        "top_event_deltascore": 0.0,
        "top_event_pos": None,
        "top_event_ref_frame_pos": None,
        "dominant_boundary": None,
        "qc_flags": "no_sites",
        "variant_class": variant_meta["variant_class"],
        "length_delta": delta,
        "align_qc": f"length_changed:{delta:+d}nt;no_sites" if delta else "",
    }


def _summarize_variant(events_df: pd.DataFrame,
                       sites_df: pd.DataFrame,
                       pkey: str,
                       report_radius: int,
                       policy: dict,
                       variant_meta: dict) -> dict:
    """
    Build the summary row for a single pkey.
    Hardened to return a valid row even if there are *no* sites/events.

    variant_meta carries the variant's own identity (`variant_class`,
    `length_delta`), which no amount of inspecting the events frame can recover:
    a deletion that changed nothing produces the same events as an SNV that
    changed nothing, and the reader of the summary needs to be able to tell them
    apart.
    """
    # Only these two are used here. `visibility_threshold` and `shift_bp` were
    # unpacked from policy and never read -- visibility is already resolved into
    # visible_flag by _build_sites_for_allele and into cls by
    # _pair_clusters_to_events, both of which run before this function sees the
    # frames.
    high_cutoff = policy["high_cutoff"]
    distance_k = policy["distance_k"]

    if events_df is None or events_df.empty or "pkey" not in events_df.columns:
        return _empty_summary_row(pkey, variant_meta)

    sub_events = events_df[events_df["pkey"] == pkey]

    if sites_df is not None and not sites_df.empty and "pkey" in sites_df.columns:
        sub_sites = sites_df[sites_df["pkey"] == pkey]
    else:
        sub_sites = pd.DataFrame()

    if sub_events.empty and sub_sites.empty:
        return _empty_summary_row(pkey, variant_meta)

    # counts
    if not sub_sites.empty:
        n_sites_wt = len(sub_sites[(sub_sites["allele"] == "WT") & (sub_sites["visible_flag"] == 1)])
        n_sites_mut = len(sub_sites[(sub_sites["allele"] == "MUT") & (sub_sites["visible_flag"] == 1)])
    else:
        n_sites_wt = 0
        n_sites_mut = 0

    n_clusters = len(sub_events)

    # global metrics
    # rows that belong to a linked redirect pair are counted via redirect_id,
    # not twice as an independent gain and an independent loss
    if "redirect_id" in sub_events.columns:
        linked = sub_events["redirect_id"].notna()
    else:
        linked = pd.Series(False, index=sub_events.index)
    unlinked = sub_events[~linked]
    linked_ev = sub_events[linked]

    global_count_gained_high = len(unlinked[(unlinked["cls"] == "gained") &
                                            (unlinked[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    global_count_lost_high = len(unlinked[(unlinked["cls"] == "lost") &
                                          (unlinked[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    global_count_shifted = (
        len(unlinked[unlinked["cls"] == "shifted"])
        + linked_ev.loc[linked_ev["redirect_class"] == "shifted", "redirect_id"].nunique()
    )
    global_count_redirected = (
        linked_ev.loc[linked_ev["redirect_class"] == "redirected", "redirect_id"].nunique()
    )

    # max abs delta
    if not sub_events["dscore"].isna().all():
        global_max_abs_deltascore = float(sub_events["dscore"].abs().max(skipna=True))
    else:
        global_max_abs_deltascore = 0.0

    # sum weighted abs delta
    swa = 0.0
    for _, ev in sub_events.iterrows():
        d = ev["distance_to_snv"]
        if pd.isna(d) or ev["dscore"] is None or pd.isna(ev["dscore"]):
            continue
        w = math.exp(- float(d) / float(distance_k)) if distance_k > 0 else 1.0
        swa += w * abs(float(ev["dscore"]))
    global_sum_weighted_abs_delta = swa

    # nearest event
    if not sub_events["distance_to_snv"].isna().all():
        nearest_event_bp_any = int(sub_events["distance_to_snv"].min(skipna=True))
    else:
        nearest_event_bp_any = None

    # local metrics
    if report_radius is None:
        report_radius = 151  # fallback

    in_local = sub_events[(sub_events["distance_to_snv"] <= report_radius)]
    if "redirect_id" in in_local.columns:
        l_linked = in_local["redirect_id"].notna()
    else:
        l_linked = pd.Series(False, index=in_local.index)
    l_unlinked = in_local[~l_linked]
    l_linked_ev = in_local[l_linked]

    local_count_gained_high = len(l_unlinked[(l_unlinked["cls"] == "gained") &
                                             (l_unlinked[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    local_count_lost_high = len(l_unlinked[(l_unlinked["cls"] == "lost") &
                                           (l_unlinked[["wt_score", "mut_score"]].max(axis=1) >= high_cutoff)])
    local_count_shifted = (
        len(l_unlinked[l_unlinked["cls"] == "shifted"])
        + l_linked_ev.loc[l_linked_ev["redirect_class"] == "shifted", "redirect_id"].nunique()
    )
    local_count_redirected = (
        l_linked_ev.loc[l_linked_ev["redirect_class"] == "redirected", "redirect_id"].nunique()
    )
    if not in_local.empty and not in_local["dscore"].isna().all():
        local_max_abs_deltascore = float(in_local["dscore"].abs().max(skipna=True))
        nearest_event_bp_local = int(in_local["distance_to_snv"].min(skipna=True))
    else:
        local_max_abs_deltascore = 0.0
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
        top_event_deltascore = float(top_ev["dscore"]) if not pd.isna(top_ev["dscore"]) else 0.0
        top_event_pos = top_ev["mut_pos"] if pd.notna(top_ev["mut_pos"]) else top_ev["wt_pos"]
        # top_event_pos is in whichever allele's coordinates the top event happened
        # to have -- MUT for anything with a mutant site, WT only for a pure loss --
        # and the column does not say which. Under an indel those two frames differ
        # by length_delta, so the summary's only coordinate silently changes frame
        # with the event class, and a consumer joining it back to the genomic
        # coordinates the mapping file is written in lands on the wrong base. The
        # WT-frame coordinate of the SAME site is reported alongside it; the two
        # are equal for every length-preserving variant. `inserted` sites carry the
        # edit anchor, which is their honest WT-frame position (see
        # _wt_frame_position) -- top_event_pos still gives the real mutant offset.
        if pd.notna(top_ev["mut_pos"]):
            trf = top_ev["mut_frame_pos"]
        else:
            trf = top_ev["wt_frame_pos"]
        top_event_ref_frame_pos = trf if pd.notna(trf) else None
    else:
        top_event_type = "none"
        top_event_deltascore = 0.0
        top_event_pos = None
        top_event_ref_frame_pos = None

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

    # Cross-allele pairing accounting, over the UNION of both alleles' clusters.
    # Counting WT clusters alone would report a 20 nt insertion as fully paired,
    # because every WT cluster does keep a counterpart -- the clusters with no
    # counterpart are all on the mutant side. Reported only when the length
    # actually changed; for an SNV or an MNV there is nothing to reconcile.
    length_delta = variant_meta["length_delta"]
    align_qc = ""
    if length_delta:
        both = sub_events["wt_pos"].notna() & sub_events["mut_pos"].notna()
        wt_only = int((sub_events["wt_pos"].notna() & sub_events["mut_pos"].isna()).sum())
        mut_only = int((sub_events["mut_pos"].notna() & sub_events["wt_pos"].isna()).sum())
        paired = int(both.sum())
        union = paired + wt_only + mut_only
        align_qc = (f"length_changed:{length_delta:+d}nt;"
                    f"clusters_paired_{paired}/{union};"
                    f"wt_only_{wt_only};mut_only_{mut_only}")
        if not sub_sites.empty and "frame_status" in sub_sites.columns:
            n_del = int((sub_sites["frame_status"] == "deleted").sum())
            n_ins = int((sub_sites["frame_status"] == "inserted").sum())
            align_qc += f";sites_in_deleted_{n_del};sites_in_inserted_{n_ins}"

    # qc_flags
    flags = []
    if n_clusters == 0:
        flags.append("no_sites")
    if nearest_event_bp_any is not None and nearest_event_bp_any > 2000:
        flags.append("far_event>2kb")
    if global_max_abs_deltascore < 1.0:
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
        "global_count_redirected": global_count_redirected,
        "global_max_abs_deltascore": global_max_abs_deltascore,
        "global_sum_weighted_abs_delta": global_sum_weighted_abs_delta,
        "nearest_event_bp_any": nearest_event_bp_any,
        "local_count_gained_high": local_count_gained_high,
        "local_count_lost_high": local_count_lost_high,
        "local_count_shifted": local_count_shifted,
        "local_count_redirected": local_count_redirected,
        "local_max_abs_deltascore": local_max_abs_deltascore,
        "nearest_event_bp_local": nearest_event_bp_local,
        "frac_effect_in_radius": frac_effect_in_radius,
        "top_event_type": top_event_type,
        "top_event_deltascore": top_event_deltascore,
        "top_event_pos": top_event_pos,
        "top_event_ref_frame_pos": top_event_ref_frame_pos,
        "dominant_boundary": dom,
        "qc_flags": qc_flags,
        "variant_class": variant_meta["variant_class"],
        "length_delta": length_delta,
        "align_qc": align_qc,
    }


# ---------------------------------------------------------------------------
# main per-gene worker
# ---------------------------------------------------------------------------

def _process_gene(fasta_path: Path,
                  gene_mapping_df: pd.DataFrame,
                  genesplicer_dir: str,
                  model_dir: str,
                  window: int,
                  report_radius: int,
                  visibility_threshold: float,
                  high_cutoff: float,
                  shift_bp: int,
                  distance_k: int,
                  cluster_radius: int,
                  max_cluster_span: int,
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
    - return events, sites, per-variant metadata, stats

    Non-SNV tokens are processed BY DEFAULT and there is no flag for it. The nt
    token grammar is uniquely decodable (the base and digit classes are disjoint),
    so no parser mode has to be selected, and whether a given variant changes the
    frame is decided per variant from len(ref) != len(alt) -- a fact of the
    record, not a user preference.

    Every rejected token is counted under its OWN reason and its token is kept in
    `skipped_detail`, which is written to genesplicer.run_summary.json. Previously
    an indel, an off-alphabet string and a token with a broken middle all landed
    in one `invalid` counter, so a run that silently dropped every indel in the
    mapping was indistinguishable from one whose mapping was malformed.

    A gene can also be abandoned WHOLESALE, before any token is looked at: an
    unreadable FASTA, a failed WT GeneSplicer run, an empty mapping, or a mapping
    with no recognisable mutant/genomic columns. Those rows were counted into
    `total_rows` and then dropped against no reason at all, so the run summary
    showed a variant count that did not reconcile with anything. `gene_skipped`
    names which of the four it was; the per-token counters cannot express it,
    because the tokens were never parsed.
    """
    gene_name = extract_gene_from_filename(str(fasta_path)) or fasta_path.stem
    stats = {
        "gene": gene_name,
        "total_rows": int(gene_mapping_df.shape[0]) if gene_mapping_df is not None else 0,
        "processed": 0,
        "skipped_validation": 0,
        "skipped_invalid": 0,
        "skipped_out_of_range": 0,
        "skipped_ref_mismatch": 0,
        "skipped_runtime": 0,
        "cluster_members_discarded": 0,
        "skipped_detail": [],
        "gene_skipped": None,
    }

    def _abandon_gene(reason):
        """Every row of this gene is being dropped before it is parsed."""
        stats["gene_skipped"] = reason
        return [], [], [], stats

    fa = read_fasta(str(fasta_path))
    if not fa:
        return _abandon_gene("fasta_unreadable")

    # adopt first sequence or key "genomic"
    if "genomic" in fa:
        wt_seq = fa["genomic"]
    else:
        wt_seq = next(iter(fa.values()))
    wt_len = len(wt_seq)

    # run WT once (full)
    try:
        wt_df = _run_genesplicer_on_seq(f"{gene_name}_WT", wt_seq, genesplicer_dir, model_dir)
    except Exception as e:
        # F28: a failed WT run must not fabricate "all gained" for the gene; skip it.
        print(f"[SKIP] {gene_name}: WT GeneSplicer run failed ({e}); skipping gene", file=sys.stderr)
        return _abandon_gene("wt_run_failed")

    # if mapping df is empty or missing, nothing to do
    if gene_mapping_df is None or gene_mapping_df.empty:
        return _abandon_gene("no_mapping_rows")

    # harmonize column names
    cols = {c.lower(): c for c in gene_mapping_df.columns}
    # expect 'mutant' and 'genomic'
    mutant_col = cols.get("mutant") or cols.get("mutation") or cols.get("mut")
    genomic_col = cols.get("genomic") or cols.get("genomic_nt") or cols.get("genomic_mut") or cols.get("genomicmutation")

    if not mutant_col or not genomic_col:
        return _abandon_gene("mapping_missing_mutant_or_genomic_column")

    events_all = []
    sites_all = []
    variants_meta = []

    def _reject(reason, mutant_tok, genomic_tok):
        """Record a dropped token by name. A bare counter bump loses which token
        it was, and returning without recording anything -- the behaviour this
        replaces for indels -- is indistinguishable from the mapping never having
        listed the variant."""
        stats[f"skipped_{reason}"] += 1
        stats["skipped_detail"].append(
            {"gene": gene_name, "mutant": mutant_tok, "genomic": genomic_tok,
             "reason": reason})

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

        # Parse the genomic token, e.g. T57261C or ACAA112217430A. parse_variant
        # never raises and returns None for a token it cannot decode, which is
        # what makes an indel distinguishable from garbage here: the old
        # try/except around get_mutation_data could not tell them apart, because
        # a valid indel and a malformed token raised the SAME ValueError, and
        # both were counted as `invalid`.
        variant = parse_variant(genomic_tok, is_nt=True)
        if variant is None:
            _reject("invalid", mutant_tok, genomic_tok)
            continue
        pos0, ref_nt, alt_nt = variant.pos0, variant.ref, variant.alt
        ref_len, alt_len = len(ref_nt), len(alt_nt)

        # Bound the END of the REF span, not its start: a multi-base REF can begin
        # inside the sequence and run off the end, and splicing it would then
        # truncate silently. (pos0 >= 0 needs no check -- Variant rejects pos < 1
        # at construction, so a guard for it here could never fire.)
        if pos0 + ref_len > wt_len:
            _reject("out_of_range", mutant_tok, genomic_tok)
            continue

        # REF guard on the WHOLE span. This pipeline had no REF guard at all --
        # ref_nt was unpacked and never read -- so a mapping whose coordinate
        # disagreed with the FASTA produced a complete, plausible, entirely
        # fictional result set. Checking wt_seq[pos0] alone would not be enough
        # either: it passes on any multi-base REF whose first base happens to
        # match.
        observed = wt_seq[pos0:pos0 + ref_len].upper()
        if observed != ref_nt.upper():
            _reject("ref_mismatch", mutant_tok, genomic_tok)
            continue

        # splice_seq honours len(ref); update_str hardcodes a one-base stride and
        # cannot express an indel. REF was just verified, so skip the re-check.
        alt_seq = splice_seq(wt_seq, pos0, ref_nt, alt_nt, validate=False)

        # run GeneSplicer on ALT
        try:
            mut_df = _run_genesplicer_on_seq(f"{gene_name}_{mutant_tok}", alt_seq, genesplicer_dir, model_dir)
        except Exception:
            _reject("runtime", mutant_tok, genomic_tok)
            continue

        stats["processed"] += 1
        variants_meta.append({
            "pkey": pkey,
            "variant_class": variant.kind,
            "length_delta": variant.length_delta,
        })

        # convert WT ALT to canonical sites
        wt_sites = _build_sites_for_allele(
            pkey=pkey,
            allele="WT",
            df=wt_df,
            visibility_threshold=visibility_threshold,
            report_radius=report_radius,
            variant_pos1=variant.pos,
            ref_len=ref_len,
            alt_len=alt_len,
        )

        mut_sites = _build_sites_for_allele(
            pkey=pkey,
            allele="MUT",
            df=mut_df,
            visibility_threshold=visibility_threshold,
            report_radius=report_radius,
            variant_pos1=variant.pos,
            ref_len=ref_len,
            alt_len=alt_len,
        )

        # merge sites
        sites_concat = pd.concat([wt_sites, mut_sites], ignore_index=True)
        # cluster
        sites_clustered = _cluster_sites(sites_concat,
                                         cluster_radius=cluster_radius,
                                         max_cluster_span=max_cluster_span)
        # pair
        events_df, n_discarded = _pair_clusters_to_events(sites_clustered,
                                                          visibility_threshold=visibility_threshold,
                                                          shift_bp=shift_bp)
        stats["cluster_members_discarded"] += n_discarded

        events_all.append(events_df)
        sites_all.append(sites_clustered)

    return events_all, sites_all, variants_meta, stats


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="GeneSplicer WT<->ALT ensemble delta caller")
    parser.add_argument("-i", "--input", required=True, help="Directory of genomic FASTA files, or a single FASTA file")
    parser.add_argument("-m", "--mapping-dir", required=True, help="Genomic mutation mapping CSV file, or directory of CSV files")
    parser.add_argument("-g", "--genesplicer-dir", default=None,
                        help="Directory containing GeneSplicer binary; "
                             "omit to use genesplicer from PATH (e.g. after conda install)")
    parser.add_argument("--model-dir", default=None,
                        help="Absolute path to a GeneSplicer model directory (one containing config_file). "
                             "If omitted, resolved from --genesplicer-dir or the genesplicer binary location "
                             "(conda ships it beside the binary as share/genesplicer-*/human).")
    parser.add_argument("-o", "--output", required=True, help="Output base directory (writes {GENE}/GeneSplicer/{GENE}.tsv, .events.tsv, .sites.tsv)")
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
    parser.add_argument("--cluster-radius", type=int, default=DEFAULT_CLUSTER_RADIUS,
                        help="Single-linkage radius (bp) for grouping sites into clusters (default 3)")
    parser.add_argument("--max-cluster-span", type=int, default=DEFAULT_MAX_CLUSTER_SPAN,
                        help="Max total span (bp) of one cluster; a site that would push the "
                             "cluster past this starts a new one. Default 0 = UNBOUNDED and "
                             "inert: a left-anchored fixed window splits genuine cross-allele "
                             "pairs that straddle its boundary, so it is off until a bound that "
                             "cuts at the largest gap exists (see _cluster_sites)")
    parser.add_argument("--redirect-max-bp", type=int, default=DEFAULT_REDIRECT_MAX_BP,
                        help="Max bp between a lost and a gained site to link them as one "
                             "redirected site (default 100)")
    parser.add_argument("--register-max-bp", type=int, default=DEFAULT_REGISTER_MAX_BP,
                        help="Linked pairs within this distance are labelled 'shifted' "
                             "(register shift) rather than 'redirected' (default 12)")
    parser.add_argument("--redirect-snv-max-bp", type=int, default=DEFAULT_REDIRECT_SNV_BP,
                        help="Max distance from the variant span (zero inside it) to the lost "
                             "site for a redirect link to be considered variant-caused "
                             "(default 80, = GeneSplicer's coding/non-coding context window)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Max parallel workers (default: half cores, capped at 8)")
    parser.add_argument("--log", help="Validation log file/dir to skip failed mutations")
    args = parser.parse_args()

    if args.genesplicer_dir:
        bin_path = os.path.join(args.genesplicer_dir, "genesplicer")
        if not os.access(bin_path, os.X_OK):
            parser.error(f"genesplicer executable not found or not executable at {bin_path}")
    elif not shutil.which("genesplicer"):
        parser.error("genesplicer not found on PATH. Install via conda (conda install -c bioconda genesplicer) "
                      "or provide --genesplicer-dir")

    # F29: resolve the human model to an absolute path (fail fast if absent) so the
    # PATH/conda default (no cd, model beside the binary) is not silently broken.
    model_dir = _resolve_model_dir(args.model_dir, args.genesplicer_dir)

    input_path = Path(args.input)
    mapping_dir = args.mapping_dir
    genesplicer_dir = args.genesplicer_dir
    output_base = args.output
    window = args.window
    report_radius = args.report_radius
    distance_k = args.distance_k
    visibility_threshold = args.visibility_threshold
    high_cutoff = args.high_cutoff
    shift_bp = args.shift_bp
    cluster_radius = args.cluster_radius
    max_cluster_span = args.max_cluster_span

    if report_radius is None:
        report_radius = window

    failure_map = load_validation_failures(args.log) if args.log else {}

    # load mappings
    mappings = load_mapping_dir(mapping_dir)

    # discover FASTA files
    if input_path.is_file():
        fasta_paths = [input_path]
    else:
        fasta_paths = [p for p in input_path.iterdir() if p.is_file() and p.suffix in (".fa", ".fasta", ".fna")]
    if not fasta_paths:
        print(f"No FASTA files found in {args.input}", file=sys.stderr)
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

    # (The former all_events/all_sites/all_pkeys accumulators are gone: they were
    # the inputs to the global pass removed in F34 and had been extended-but-never-
    # read ever since.)
    events_by_gene = {}
    sites_by_gene = {}
    variants_by_gene = {}

    total_genes = len(fasta_paths)
    genes_completed = 0
    total_variants_seen = 0
    total_variants_processed = 0
    total_skipped_validation = 0
    total_skipped_invalid = 0
    total_skipped_out_of_range = 0
    total_skipped_ref_mismatch = 0
    total_skipped_runtime = 0
    total_cluster_members_discarded = 0
    skipped_detail = []
    genes_skipped = []

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
                model_dir,
                window,
                report_radius,
                visibility_threshold,
                high_cutoff,
                shift_bp,
                distance_k,
                cluster_radius,
                max_cluster_span,
                failure_map,
            )
            futs.append(fut)
            fut_meta[fut] = gene_name

        for fut in concurrent.futures.as_completed(futs):
            gene_events, gene_sites, gene_variants, gene_stats = fut.result()

            genes_completed += 1
            total_variants_seen += gene_stats.get("total_rows", 0)
            total_variants_processed += gene_stats.get("processed", 0)
            total_skipped_validation += gene_stats.get("skipped_validation", 0)
            total_skipped_invalid += gene_stats.get("skipped_invalid", 0)
            total_skipped_out_of_range += gene_stats.get("skipped_out_of_range", 0)
            total_skipped_ref_mismatch += gene_stats.get("skipped_ref_mismatch", 0)
            total_skipped_runtime += gene_stats.get("skipped_runtime", 0)
            total_cluster_members_discarded += gene_stats.get("cluster_members_discarded", 0)
            skipped_detail.extend(gene_stats.get("skipped_detail", []))

            gene_name = gene_stats.get("gene") or fut_meta.get(fut, "unknown")
            events_by_gene[gene_name] = gene_events
            sites_by_gene[gene_name] = gene_sites
            variants_by_gene[gene_name] = gene_variants
            # A gene abandoned before its tokens were parsed contributes rows to
            # variants_total that no per-token counter can account for. Without
            # this the run summary simply does not add up and the reader cannot
            # tell an empty mapping from an unreadable FASTA.
            if gene_stats.get("gene_skipped"):
                genes_skipped.append({
                    "gene": gene_name,
                    "reason": gene_stats["gene_skipped"],
                    "rows_dropped": gene_stats.get("total_rows", 0),
                })
            skips = []
            if gene_stats.get("skipped_validation", 0):
                skips.append(f"validation={gene_stats['skipped_validation']}")
            if gene_stats.get("skipped_invalid", 0):
                skips.append(f"invalid={gene_stats['skipped_invalid']}")
            if gene_stats.get("skipped_out_of_range", 0):
                skips.append(f"out_of_range={gene_stats['skipped_out_of_range']}")
            if gene_stats.get("skipped_ref_mismatch", 0):
                skips.append(f"ref_mismatch={gene_stats['skipped_ref_mismatch']}")
            if gene_stats.get("skipped_runtime", 0):
                skips.append(f"runtime={gene_stats['skipped_runtime']}")
            if gene_stats.get("gene_skipped"):
                skips.append(f"GENE ABANDONED: {gene_stats['gene_skipped']}")
            skip_msg = ", ".join(skips) if skips else "none"
            print(
                f"[{genes_completed}/{total_genes}] {gene_name}: "
                f"processed {gene_stats.get('processed', 0)}/{gene_stats.get('total_rows', 0)} variants "
                f"(skipped: {skip_msg}). "
                f"Total processed so far: {total_variants_processed}"
            )
            sys.stdout.flush()

    # F34: the former global events_df/sites_df/summary_df built here (concat +
    # reconcile + priority + is_high_impact + per-pkey summarize across all genes)
    # were never written — only the per-gene recomputation below ships — so that
    # entire pass is removed. What the per-gene loop still needs is kept: the two
    # column contracts and the two scorer closures (lifted out of the old
    # `if not events_df.empty` block so they are always defined, even when a run
    # produces zero events).
    events_required = EVENT_BASE_COLUMNS + [
        "redirect_id", "redirect_dpos", "score_transfer_ratio",
        "redirect_class", "cls_reconciled",
    ]
    sites_required = list(SITE_COLUMNS)

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
        rc = row.get("redirect_class")
        if (pd.notna(rc) and rc == "shifted") or (
            row["cls"] == "shifted" and row["dpos"] is not None
            and not pd.isna(row["dpos"]) and abs(int(row["dpos"])) >= shift_bp
        ):
            bonus += 1.0
        return base * w + bonus

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

    # write per-gene outputs
    for gname in events_by_gene:
        g_events = events_by_gene[gname]
        g_sites = sites_by_gene[gname]
        g_variants = variants_by_gene[gname]

        g_events_df = pd.concat(g_events, ignore_index=True) if g_events else pd.DataFrame(columns=events_required)
        g_sites_df = pd.concat(g_sites, ignore_index=True) if g_sites else pd.DataFrame(columns=sites_required)
        for c in events_required:
            if c not in g_events_df.columns:
                g_events_df[c] = pd.NA
        for c in sites_required:
            if c not in g_sites_df.columns:
                g_sites_df[c] = pd.NA

        # F30: same reconciliation on the per-gene frame; redirect_id is
        # position-derived so it matches the global run's ids exactly.
        g_events_df = _reconcile_redirects(
            g_events_df,
            max_bp=args.redirect_max_bp,
            register_bp=args.register_max_bp,
            snv_max_bp=args.redirect_snv_max_bp,
        )

        if not g_events_df.empty:
            g_events_df["priority"] = g_events_df.apply(_calc_priority, axis=1)
            g_events_df["is_high_impact"] = g_events_df.apply(_is_hi, axis=1)

        g_summary_rows = []
        for vmeta in g_variants:
            g_summary_rows.append(
                _summarize_variant(
                    events_df=g_events_df,
                    sites_df=g_sites_df,
                    pkey=vmeta["pkey"],
                    report_radius=report_radius,
                    policy=policy,
                    variant_meta=vmeta,
                )
            )
        g_summary_df = pd.DataFrame(g_summary_rows)

        gene_out_dir = Path(output_base) / gname / "GeneSplicer"
        gene_out_dir.mkdir(parents=True, exist_ok=True)
        g_summary_df.to_csv(str(gene_out_dir / f"{gname}.tsv"), sep="\t", index=False)
        g_events_df.to_csv(str(gene_out_dir / f"{gname}.events.tsv"), sep="\t", index=False)
        g_sites_df.to_csv(str(gene_out_dir / f"{gname}.sites.tsv"), sep="\t", index=False)
        print(f"wrote {gname} outputs to {gene_out_dir}")

    # Durable skip accounting. The printed line below is ephemeral, and a token
    # that was dropped is exactly the thing a reader of the outputs cannot
    # otherwise discover -- there is no row for it anywhere.
    run_summary = {
        "genes_total": total_genes,
        "genes_processed": genes_completed,
        "variants_total": total_variants_seen,
        "variants_processed": total_variants_processed,
        "skipped_validation": total_skipped_validation,
        "skipped_invalid": total_skipped_invalid,
        "skipped_out_of_range": total_skipped_out_of_range,
        "skipped_ref_mismatch": total_skipped_ref_mismatch,
        "skipped_runtime": total_skipped_runtime,
        "cluster_members_discarded": total_cluster_members_discarded,
        "genes_skipped": genes_skipped,
        "rows_dropped_with_gene": sum(g["rows_dropped"] for g in genes_skipped),
        "skipped_detail": skipped_detail,
    }
    # The ledger must reconcile, or a reader cannot tell "nothing to report" from
    # "something went unrecorded". Anything left over is a drop with no reason
    # attached, which is the one outcome this file exists to make impossible.
    unaccounted = (total_variants_seen - total_variants_processed
                   - total_skipped_validation - total_skipped_invalid
                   - total_skipped_out_of_range - total_skipped_ref_mismatch
                   - total_skipped_runtime - run_summary["rows_dropped_with_gene"])
    run_summary["variants_unaccounted"] = unaccounted
    Path(output_base).mkdir(parents=True, exist_ok=True)
    with open(Path(output_base) / "genesplicer.run_summary.json", "w") as fh:
        json.dump(run_summary, fh, indent=2)

    print(
        "Run summary: "
        f"{genes_completed}/{total_genes} genes processed, "
        f"{total_variants_processed}/{total_variants_seen} variants analysed; "
        f"skipped validation={total_skipped_validation}, "
        f"invalid={total_skipped_invalid}, "
        f"out_of_range={total_skipped_out_of_range}, "
        f"ref_mismatch={total_skipped_ref_mismatch}, "
        f"runtime={total_skipped_runtime}; "
        f"genes abandoned={len(genes_skipped)} ({run_summary['rows_dropped_with_gene']} rows); "
        f"unaccounted={unaccounted}; "
        f"non-top cluster members discarded={total_cluster_members_discarded} "
        f"-> genesplicer.run_summary.json"
    )


if __name__ == "__main__":
    main()
