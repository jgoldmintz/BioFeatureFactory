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

"""Wrapper for run_spliceai process with partial-cache support."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, **kwargs):
    return subprocess.run(cmd, **kwargs)


def count_variants(vcf_path: Path) -> int:
    count = 0
    if not vcf_path.exists():
        return 0
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.strip():
                count += 1
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description="Run SpliceAI with partial cache handling")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--vcf", required=True, help="Bgzipped VCF input")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--cache", required=True, help="Partial cache directory")
    parser.add_argument("--prune", required=True, help="Path to prune_vcf.py")
    parser.add_argument("--merge", required=True, help="Path to merge_vcfs.py")
    parser.add_argument("--skip-cache", action="store_true",
                        help="Bypass partial cache for this gene (used when annotation is filtered)")
    args = parser.parse_args()

    gene = args.gene
    vcf_gz = Path(args.vcf)
    reference = args.reference
    annotation = args.annotation
    cache_dir = Path(args.cache)
    prune_script = Path(args.prune)
    merge_script = Path(args.merge)

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{gene}.partial.vcf"

    input_vcf = Path(f"{gene}.input.vcf")
    remaining_vcf = Path(f"{gene}.remaining.vcf")
    remaining_vcf_gz = Path(f"{gene}.remaining.vcf.gz")
    new_vcf = Path(f"{gene}.spliceai.new.vcf")
    final_vcf = Path(f"{gene}.spliceai.vcf")

    # Check if cache should be skipped
    if args.skip_cache:
        print(f"[run_spliceai_cached] {gene}: Skipping partial cache (filtered annotation)", flush=True)
        # Decompress VCF directly to remaining_vcf
        with input_vcf.open("wb") as handle:
            run_cmd(["bgzip", "-dc", str(vcf_gz)], stdout=handle, check=True)
        shutil.copy2(input_vcf, remaining_vcf)

        # Compress and index
        with remaining_vcf_gz.open("wb") as handle:
            run_cmd(["bgzip", "-c", str(remaining_vcf)], stdout=handle, check=True)
        run_cmd(["tabix", "-p", "vcf", str(remaining_vcf_gz)], check=True)

        # Run SpliceAI directly to final output
        splice_cmd = [
            "spliceai",
            "-I", str(remaining_vcf_gz),
            "-O", str(final_vcf),
            "-R", reference,
            "-A", annotation,
        ]
        result = run_cmd(splice_cmd)
        status = result.returncode

        if status != 0:
            norm_status = 128 + abs(status) if status < 0 else status
            print(f"[run_spliceai_cached] {gene} failed (exit {status})", file=sys.stderr)
            return norm_status

        if not final_vcf.exists() or final_vcf.stat().st_size == 0:
            print(f"[run_spliceai_cached] WARNING: SpliceAI output missing for {gene}", file=sys.stderr)
            return 1

        return 0

    # Normal path with caching
    with input_vcf.open("wb") as handle:
        run_cmd(["bgzip", "-dc", str(vcf_gz)], stdout=handle, check=True)

    if cache_file.exists():
        run_cmd(["python3", str(prune_script), "--input", str(input_vcf), "--already", str(cache_file), "--output", str(remaining_vcf)], check=True)
    else:
        shutil.copy2(input_vcf, remaining_vcf)

    remaining_count = count_variants(remaining_vcf)
    if remaining_count == 0:
        print(f"[run_spliceai_cached] {gene} already complete via partial cache ({cache_file})", flush=True)
        source = cache_file if cache_file.exists() else remaining_vcf
        shutil.copy2(source, final_vcf)
        return 0

    with remaining_vcf_gz.open("wb") as handle:
        run_cmd(["bgzip", "-c", str(remaining_vcf)], stdout=handle, check=True)
    run_cmd(["tabix", "-p", "vcf", str(remaining_vcf_gz)], check=True)

    splice_cmd = [
        "spliceai",
        "-I", str(remaining_vcf_gz),
        "-O", str(new_vcf),
        "-R", reference,
        "-A", annotation,
    ]
    result = run_cmd(splice_cmd)
    status = result.returncode

    if not new_vcf.exists() or new_vcf.stat().st_size == 0:
        print(f"[run_spliceai_cached] WARNING: new SpliceAI output missing for {gene}", file=sys.stderr)

    if status != 0:
        norm_status = 128 + abs(status) if status < 0 else status
        print(f"[run_spliceai_cached] {gene} failed (exit {status}), caching partial output", file=sys.stderr)
        if new_vcf.exists() and new_vcf.stat().st_size > 0:
            shutil.copy2(new_vcf, cache_file)
        return norm_status

    merge_cmd = ["python3", str(merge_script), "--new", str(new_vcf), "--output", str(final_vcf)]
    if cache_file.exists():
        merge_cmd[2:2] = ["--existing", str(cache_file)]
    run_cmd(merge_cmd, check=True)
    if cache_file.exists():
        cache_file.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
