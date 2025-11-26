#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023–2025  Jacob Goldmintz
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

"""Remove variants that already exist in a partial SpliceAI VCF."""

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple


def open_vcf(path: Path):
    if not path:
        raise ValueError("VCF path is required")
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def variant_key(line: str) -> Tuple[str, str, str, str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        raise ValueError(f"Malformed VCF line: {line!r}")
    chrom, pos, _vid, ref, alt = fields[:5]
    return chrom, pos, ref, alt


def load_processed(path: Path) -> Set[Tuple[str, str, str, str]]:
    processed: Set[Tuple[str, str, str, str]] = set()
    if not path or not path.exists():
        return processed
    with open_vcf(path) as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            raw = raw.strip()
            if not raw:
                continue
            try:
                processed.add(variant_key(raw))
            except ValueError:
                continue
    return processed


def prune_variants(source: Path, already: Path | None, dest: Path) -> Tuple[int, int]:
    processed = load_processed(already) if already else set()
    kept = 0
    total = 0
    with open_vcf(source) as src, dest.open("w") as out:
        for raw in src:
            if raw.startswith("#") or not raw.strip():
                out.write(raw)
                continue
            total += 1
            try:
                key = variant_key(raw)
            except ValueError:
                continue
            if key in processed:
                continue
            out.write(raw)
            kept += 1
    return total, kept


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Emit a VCF that excludes variants already present in another VCF."
    )
    parser.add_argument("--input", required=True, help="Input VCF (plain or gzipped)")
    parser.add_argument("--already", help="Partial/cache VCF to treat as already processed")
    parser.add_argument("--output", required=True, help="Path to write the reduced VCF")
    args = parser.parse_args()

    source = Path(args.input)
    dest = Path(args.output)
    already = Path(args.already) if args.already else None

    if not source.exists():
        print(f"[prune_vcf] Input {source} does not exist", file=sys.stderr)
        return 2

    dest.parent.mkdir(parents=True, exist_ok=True)
    total, kept = prune_variants(source, already, dest)
    print(f"[prune_vcf] total={total} kept={kept}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
