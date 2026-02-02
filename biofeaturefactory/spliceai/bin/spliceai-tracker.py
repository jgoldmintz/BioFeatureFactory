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
Live SpliceAI progress tracker.

The script inspects the Nextflow work directory, discovers every running
`run_spliceai` task, and prints processed/total mutation counts per gene.
Totals come from the mutation CSVs, while progress comes from the growing
`<gene>.spliceai.vcf` inside each task work dir.
"""
import argparse
import gzip
import pathlib
import re
import sys
import time
from datetime import datetime
from typing import Dict, Iterable, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream SpliceAI progress as processed/total mutations per gene."
    )
    parser.add_argument(
        "--mutations-path",
        type=pathlib.Path,
        required=True,
        help="Path to mutation CSV directory or single CSV.",
    )
    parser.add_argument(
        "--work-root",
        type=pathlib.Path,
        default=pathlib.Path("work"),
        help="Nextflow work directory (default: ./work).",
    )
    parser.add_argument(
        "--poll-seconds",
        type=float,
        default=5.0,
        help="Seconds between refreshes (default: 5).",
    )
    parser.add_argument(
        "--log-file",
        type=pathlib.Path,
        help="Optional log file path to append tracker output.",
    )
    return parser.parse_args()


def open_text(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def looks_like_header(line: str) -> bool:
    line = line.lower()
    tokens = ("mutant", "mutation", "gene", "chrom", "ref", "alt")
    return any(tok in line for tok in tokens)


def count_records(csv_path: Optional[pathlib.Path]) -> int:
    if not csv_path or not csv_path.exists():
        return 0
    first = None
    count = 0
    with open_text(csv_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if first is None:
                first = line
                continue
            count += 1
    if first and not looks_like_header(first):
        count += 1
    return count


def count_variants(vcf_path: pathlib.Path) -> int:
    if not vcf_path.exists():
        return 0
    with open_text(vcf_path) as fh:
        return sum(1 for line in fh if line and not line.startswith("#"))


def build_mut_index(base: pathlib.Path) -> Dict[str, pathlib.Path]:
    base = base.expanduser().resolve()
    if not base.exists():
        return {}
    index: Dict[str, pathlib.Path] = {}
    if base.is_file():
        key = base.stem.lower().replace("_mutations", "")
        index[key] = base
        return index
    for csv in base.glob("*.csv"):
        key = csv.stem.lower().replace("_mutations", "")
        index.setdefault(key, csv)
    return index


def find_mut_file(gene: str, index: Dict[str, pathlib.Path]) -> Optional[pathlib.Path]:
    if not index:
        return None
    gene_key = gene.lower()
    if gene_key in index:
        return index[gene_key]
    for key, path in index.items():
        if gene_key in key or key in gene_key:
            return path
    return None


GENE_ARG_RE = re.compile(r"--gene\s+([A-Za-z0-9_.-]+)")


def active_tasks(work_root: pathlib.Path) -> Iterable[Tuple[str, pathlib.Path]]:
    for sub in work_root.glob("*/*"):
        if not (sub / ".command.run").exists():
            continue
        if (sub / ".exitcode").exists():
            continue  # already finished
        cmd_path = sub / ".command.sh"
        if not cmd_path.exists():
            continue
        cmd = cmd_path.read_text()
        gene = None
        match = GENE_ARG_RE.search(cmd)
        if match:
            gene = match.group(1)
        else:
            marker = '-O "'
            idx = cmd.find(marker)
            if idx != -1:
                rest = cmd[idx + len(marker) :]
                gene = rest.split(".spliceai.vcf", 1)[0]
        if gene:
            yield gene, sub


def emit(lines: Iterable[str], log_fp):
    block = "\n".join(lines)
    print(block, flush=True)
    if log_fp:
        log_fp.write(block + "\n")
        log_fp.flush()


def main() -> int:
    args = parse_args()
    work_root = args.work_root.expanduser().resolve()
    mut_index = build_mut_index(args.mutations_path)
    totals_cache: Dict[str, Optional[int]] = {}
    prev_rows: Dict[str, Tuple[int, Optional[int]]] = {}
    log_fp = None
    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
        log_fp = args.log_file.open("a")

    try:
        while True:
            rows: Dict[str, Tuple[int, Optional[int]]] = {}
            for gene, workdir in active_tasks(work_root):
                processed = count_variants(workdir / f"{gene}.spliceai.vcf")
                if gene not in totals_cache:
                    mut_path = find_mut_file(gene, mut_index)
                    totals_cache[gene] = count_records(mut_path)
                rows[gene] = (processed, totals_cache.get(gene))
            if rows and rows != prev_rows:
                lines = [datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
                for gene in sorted(rows):
                    processed, total = rows[gene]
                    total_str = total if total is not None else "?"
                    lines.append(f"{gene:10s} {processed}/{total_str}")
                lines.append("")
                emit(lines, log_fp)
                prev_rows = rows
            time.sleep(max(args.poll_seconds, 0.5))
    except KeyboardInterrupt:
        return 130
    finally:
        if log_fp:
            log_fp.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
