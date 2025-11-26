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

"""Merge two VCFs while deduplicating variant records."""

import argparse
import gzip
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def variant_key(line: str) -> Tuple[str, str, str, str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        raise ValueError(f"Malformed VCF line: {line!r}")
    chrom, pos, _vid, ref, alt = fields[:5]
    return chrom, pos, ref, alt


def read_header(path: Path) -> List[str]:
    lines: List[str] = []
    with open_vcf(path) as fh:
        for raw in fh:
            if not raw.startswith("#"):
                break
            lines.append(raw)
            if raw.startswith("#CHROM"):
                break
    return lines


def write_body(path: Path, out, seen: Set[Tuple[str, str, str, str]]) -> None:
    if not path.exists():
        return
    with open_vcf(path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            try:
                key = variant_key(raw)
            except ValueError:
                continue
            if key in seen:
                continue
            seen.add(key)
            out.write(raw)


def merge(existing: Optional[Path], new_path: Path, dest: Path) -> None:
    header_source = existing if existing and existing.exists() else new_path
    header_lines = read_header(header_source) if header_source and header_source.exists() else []
    dest.parent.mkdir(parents=True, exist_ok=True)
    seen: Set[Tuple[str, str, str, str]] = set()
    with dest.open("w") as out:
        for line in header_lines:
            out.write(line)
        if header_lines and not header_lines[-1].startswith("#CHROM"):
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        if existing:
            write_body(existing, out, seen)
        write_body(new_path, out, seen)


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge partial and new VCF outputs.")
    parser.add_argument("--existing", help="Existing partial VCF to include first")
    parser.add_argument("--new", required=True, help="Newly generated VCF to append")
    parser.add_argument("--output", required=True, help="Destination VCF")
    args = parser.parse_args()

    existing = Path(args.existing) if args.existing else None
    new_path = Path(args.new)
    dest = Path(args.output)

    if not new_path.exists():
        print(f"[merge_vcfs] New VCF {new_path} does not exist", file=sys.stderr)
        return 2

    merge(existing, new_path, dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
