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
Build the human genome-wide codon usage table used as rare_codon's null model.

WHY THIS EXISTS
---------------
cg_cotrans defines a "rare codon" relative to GENOME-WIDE codon statistics: its
README says to run calc_codon_usage.py once across many genes (`fasta` is
documented as "input MSA fasta file(s)", plural) and only then run
calc_rare_enrichment.py per gene.

BFF never does that. rare_codon_pipeline.py auto-builds a usage table when
--usage is absent, and that auto-build fails on every real input tried so far:

    SMN2:  KeyError: 'gi|556503834|ref|NC_000913.3|'   cg_cotrans' default E. coli WT GI
    F9:    KeyError: 'NNN'                             IUPAC ambiguity codons

so _build_usage_pgz_from_msa takes over and derives "rare" from the SINGLE gene
under analysis -- making the reference distribution and the test sequence the
same data. Measured against the table this script builds, that self-referential
fallback produced 7 false rare codons and missed 2 real ones out of a true set
of 4, including calling GCC rare at 9.5% when it is 39.9% of all alanine codons
in the human transcriptome.

SOURCE
------
FDA/Cancer-CoCoPUTs <https://github.com/FDA/Cancer-CoCoPUTs>, files
GRCh38_bicodons_w_ENSG_1.zip .. _10.zip. Each row is one gene: Gene (ENSG),
Transcript (ENST, primary), then 4096 bicodon columns (64x64 codon pairs).
19,322 genes.

Two upstream traps this script handles:
  * Shard 1 carries the header; shards 2-10 do NOT. Upstream intends them to be
    concatenated (`cat file1 ... file10 > combined.tsv`). Reading them as
    independent TSVs either crashes on a header mismatch or silently eats the
    first gene of every shard.
  * Their Readme names a GRCh38_codons_w_ENSG.tsv, but that file is NOT in the
    repository -- only bicodons, junction dinucleotides and categories. Hence
    the derivation below.

DERIVATION
----------
Codon counts are recovered from bicodons: for codon c,

    count(c) = sum over all 4096 bicodon columns whose FIRST codon is c

A CDS of N codons yields N-1 bicodons and codon position k is the first element
of bicodon k for k = 1..N-1, so this counts every codon except each transcript's
TERMINAL one. At a mean of 565 sense codons per gene that is a 0.177% omission,
and the omitted codon is the stop, which is excluded from the output anyway.

Relative usage is computed WITHIN each amino-acid family, matching both
calc_codon_usage.py:172 and rare_codon_pipeline.py's fallback, so the values drop
straight into the slots load_null_model reads at
overall_codon_usage[gi][c] / unweighted_codon_usage[gi][c].

VALIDATION (spot checks against published human codon usage)
    CTG Leu 39.4% (~40%)   GCC Ala 39.9% (~40%)   AGA Arg 21.6% (~20%)
    TTA Leu  7.9% (~7%)    TCG Ser  5.5% (~4%)    ATG Met 100.0%
Rare set at cg_cotrans' default threshold 0.1: CGT, CTA, TCG, TTA (4 of 61).

USAGE
    build_cocoputs_cut.py --out-dir <dir> [--keep-shards] [--force]

Writes <dir>/human_GRCh38_codon_usage.tsv. Idempotent: exits early if the table
already exists unless --force. Consumed by:
    rare_codon_pipeline.py --reference-codon-usage <dir>/human_GRCh38_codon_usage.tsv

CITATION
    Alexaki A. et al., Codon and Codon-Pair Usage Tables (CoCoPUTs), J Mol Biol
    431(13):2434-2441, 2019. Cancer-CoCoPUTs V1.1, Douglas Meyer, US FDA.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import urllib.error
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

BASE_URL = ("https://raw.githubusercontent.com/FDA/Cancer-CoCoPUTs/master/"
            "GRCh38_bicodons_w_ENSG_{n}.zip")
N_SHARDS = 10
OUT_NAME = "human_GRCh38_codon_usage.tsv"
EXPECTED_BICODON_COLS = 4096          # 64 x 64
EXPECTED_GENES_MIN = 19000            # upstream says ~19,000; guard against a truncated pull

# The 64-codon table, local so this script does not depend on an importable BFF.
_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON_TO_AA = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(
        (x, y, z) for x in _BASES for y in _BASES for z in _BASES)
}


def download_shards(work: Path, force: bool = False) -> list[Path]:
    """Fetch the 10 bicodon zips. Skips any already present unless --force."""
    zips = []
    for n in range(1, N_SHARDS + 1):
        dest = work / f"GRCh38_bicodons_w_ENSG_{n}.zip"
        if dest.exists() and dest.stat().st_size > 0 and not force:
            print(f"  EXISTS {dest.name}")
        else:
            url = BASE_URL.format(n=n)
            print(f"  GET    {dest.name}")
            try:
                urllib.request.urlretrieve(url, dest)
            except (urllib.error.URLError, urllib.error.HTTPError) as exc:
                raise RuntimeError(f"download failed for {url}: {exc}") from exc
        zips.append(dest)
    return zips


def extract_shards(zips: list[Path], work: Path) -> list[Path]:
    tsvs = []
    for z in zips:
        out_dir = work / z.stem
        out_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(z) as zf:
            for name in zf.namelist():
                target = out_dir / Path(name).name
                if not target.exists() or target.stat().st_size == 0:
                    zf.extract(name, out_dir)
                tsvs.append(out_dir / name)
    # Numeric order matters: only shard 1 carries the header.
    return sorted(tsvs, key=lambda p: int(re.search(r"_(\d+)\.tsv$", p.name).group(1)))


def accumulate_codon_counts(tsvs: list[Path]) -> tuple[dict, int]:
    header = open(tsvs[0]).readline().rstrip("\n").split("\t")
    if header[0] != "Gene" or len(header) != EXPECTED_BICODON_COLS + 2:
        raise RuntimeError(
            f"unexpected header in {tsvs[0].name}: first field {header[0]!r}, "
            f"{len(header)} columns (expected 'Gene' and {EXPECTED_BICODON_COLS + 2})")
    first_codon = [b[:3] for b in header[2:]]

    counts: dict[str, int] = defaultdict(int)
    genes = 0
    for path in tsvs:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts[0] == "Gene":       # header row, shard 1 only
                    continue
                if len(parts) != EXPECTED_BICODON_COLS + 2:
                    continue                 # malformed row; counted via `genes` shortfall
                genes += 1
                for i, v in enumerate(parts[2:]):
                    if v != "0":
                        counts[first_codon[i]] += int(v)
    if genes < EXPECTED_GENES_MIN:
        raise RuntimeError(
            f"only {genes} genes parsed (expected >= {EXPECTED_GENES_MIN}); "
            f"a shard is probably truncated or missing")
    return counts, genes


def write_table(counts: dict, out_path: Path) -> tuple[int, list[str]]:
    aa_codons: dict[str, list[str]] = defaultdict(list)
    sense = []
    for codon, aa in CODON_TO_AA.items():
        if aa != "*":
            sense.append(codon)
            aa_codons[aa].append(codon)

    rel = {}
    for aa, codons in aa_codons.items():
        denom = sum(counts[c] for c in codons)
        for c in codons:
            rel[c] = counts[c] / denom if denom else 0.0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("codon\tamino_acid\tcount\trelative_usage_within_aa\tis_rare_at_0.1\n")
        for c in sorted(sense, key=lambda x: (CODON_TO_AA[x], x)):
            f.write(f"{c}\t{CODON_TO_AA[c]}\t{counts[c]}\t{rel[c]:.6f}"
                    f"\t{1 if rel[c] <= 0.1 else 0}\n")
    return sum(counts[c] for c in sense), sorted(c for c in sense if rel[c] <= 0.1)


def validate(out_path: Path) -> bool:
    """Spot-check against published human values. Catches a silently wrong build."""
    rel = {}
    with open(out_path) as f:
        next(f)
        for line in f:
            codon, aa, count, r, _ = line.rstrip("\n").split("\t")
            rel[codon] = float(r)
    # (codon, expected, tolerance) -- wide enough to absorb assembly/annotation drift,
    # tight enough that a broken derivation (wrong column axis, missing shards) fails.
    checks = [("CTG", 0.40, 0.05), ("GCC", 0.40, 0.05), ("AGA", 0.20, 0.05),
              ("TTA", 0.08, 0.03), ("TCG", 0.05, 0.03), ("ATG", 1.00, 0.001)]
    ok = True
    for codon, expected, tol in checks:
        got = rel.get(codon, -1)
        flag = "ok" if abs(got - expected) <= tol else "FAIL"
        if flag == "FAIL":
            ok = False
        print(f"    {codon} {CODON_TO_AA[codon]}  {got:6.1%}  expected ~{expected:.0%}  {flag}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Download CoCoPUTs bicodon data and build the human "
                    "genome-wide codon usage table for rare_codon's null model.")
    ap.add_argument("--out-dir", required=True,
                    help="Directory for the table and the downloaded shards")
    ap.add_argument("--keep-shards", action="store_true",
                    help="Keep the extracted per-gene TSVs (~165 MB). They are what "
                         "calc_codon_usage.py would need for a real multi-gene run "
                         "and what the 'groups' null model requires.")
    ap.add_argument("--force", action="store_true",
                    help="Rebuild even if the table already exists")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_path = out_dir / OUT_NAME
    if out_path.exists() and out_path.stat().st_size > 0 and not args.force:
        print(f"  EXISTS {out_path} (use --force to rebuild)")
        return 0

    work = out_dir / "_cocoputs_shards"
    work.mkdir(parents=True, exist_ok=True)
    try:
        zips = download_shards(work, force=args.force)
        tsvs = extract_shards(zips, work)
        counts, genes = accumulate_codon_counts(tsvs)
        total, rare = write_table(counts, out_path)
    except Exception as exc:                       # noqa: BLE001 - reported, not swallowed
        print(f"  ERROR {exc}", file=sys.stderr)
        return 1

    print(f"  genes parsed: {genes}   sense-codon observations: {total:,}")
    print(f"  wrote {out_path}")
    print(f"  rare codons at threshold 0.1: {len(rare)} of 61  {rare}")
    print("  validation:")
    if not validate(out_path):
        print("  ERROR validation failed; the table is NOT trustworthy", file=sys.stderr)
        return 1

    if not args.keep_shards:
        for p in glob.glob(str(work / "**" / "*.tsv"), recursive=True):
            os.remove(p)
        print(f"  removed extracted TSVs (zips kept in {work}); --keep-shards to retain")
    return 0


if __name__ == "__main__":
    sys.exit(main())
