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
Filter SpliceAI annotation file using hybrid stratified sampling.

For genes with > max_isoforms transcripts, applies hybrid sampling:
  - Takes top N longest isoforms (default: 10)
  - Randomly samples remaining slots from the rest
  - Uses deterministic seed based on gene name for reproducibility
"""

import argparse
import hashlib
import random
import sys
from pathlib import Path
from typing import List, Tuple


def hash_gene_name(gene: str) -> int:
    """Generate deterministic seed from gene name."""
    return int(hashlib.md5(gene.encode()).hexdigest(), 16) % (2**32)


def parse_annotation_line(line: str) -> Tuple[str, int]:
    """
    Parse annotation line and return (gene_name, transcript_length).

    Format: #NAME CHROM STRAND TX_START TX_END EXON_START EXON_END
    """
    parts = line.strip().split('\t')
    if len(parts) < 5:
        return None, 0

    gene_name = parts[0]
    try:
        tx_start = int(parts[3])
        tx_end = int(parts[4])
        length = tx_end - tx_start
        return gene_name, length
    except (ValueError, IndexError):
        return gene_name, 0


def filter_annotation(
    annotation_file: Path,
    target_gene: str,
    max_isoforms: int,
    top_n: int,
    output_file: Path = None,
    log_file: Path = None
) -> None:
    """
    Filter annotation file for target gene using hybrid stratified sampling.

    Args:
        annotation_file: Path to input annotation file
        target_gene: Gene name to filter
        max_isoforms: Maximum isoforms to output
        top_n: Number of longest isoforms to always include
        output_file: Output file (None = stdout)
        log_file: Log file for filtering events (None = stderr)
    """
    # Read all lines for the target gene
    gene_lines = []
    header_line = None

    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_line = line
                continue

            gene_name, length = parse_annotation_line(line)
            if gene_name == target_gene:
                gene_lines.append((line.strip(), length))

    total_isoforms = len(gene_lines)

    # Determine if filtering is needed
    if total_isoforms <= max_isoforms:
        # No filtering needed
        log_message = f"[FILTER] {target_gene}: {total_isoforms} isoforms (no filtering needed)"
        selected_lines = [line for line, _ in gene_lines]
    else:
        # Apply hybrid stratified sampling
        # Sort by length (descending)
        gene_lines.sort(key=lambda x: x[1], reverse=True)

        # Take top N longest
        top_longest = gene_lines[:top_n]
        remaining = gene_lines[top_n:]

        # Randomly sample from remaining
        random_count = max_isoforms - top_n
        seed = hash_gene_name(target_gene)
        random.seed(seed)

        if len(remaining) <= random_count:
            # Not enough remaining, take all
            sampled_remaining = remaining
        else:
            sampled_remaining = random.sample(remaining, random_count)

        # Combine
        selected = top_longest + sampled_remaining
        selected_lines = [line for line, _ in selected]

        log_message = (
            f"[FILTER] {target_gene}: {total_isoforms} isoforms -> {len(selected_lines)} "
            f"(top {len(top_longest)} longest + {len(sampled_remaining)} random, seed={seed})"
        )

    # Write log
    if log_file:
        with open(log_file, 'a') as f:
            f.write(log_message + '\n')
    else:
        print(log_message, file=sys.stderr)

    # Write output
    out_handle = open(output_file, 'w') if output_file else sys.stdout

    try:
        if header_line:
            out_handle.write(header_line)
        for line in selected_lines:
            out_handle.write(line + '\n')
    finally:
        if output_file:
            out_handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Filter SpliceAI annotation using hybrid stratified sampling"
    )
    parser.add_argument(
        '--annotation', required=True, type=Path,
        help='Input annotation file'
    )
    parser.add_argument(
        '--gene', required=True,
        help='Gene name to filter'
    )
    parser.add_argument(
        '--max-isoforms', type=int, default=50,
        help='Maximum isoforms per gene (default: 50)'
    )
    parser.add_argument(
        '--top-n', type=int, default=10,
        help='Number of longest isoforms to always include (default: 10)'
    )
    parser.add_argument(
        '--output', type=Path,
        help='Output file (default: stdout)'
    )
    parser.add_argument(
        '--log-file', type=Path,
        help='Log file for filtering events (default: stderr)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.annotation.exists():
        print(f"ERROR: Annotation file not found: {args.annotation}", file=sys.stderr)
        sys.exit(1)

    if args.top_n >= args.max_isoforms:
        print(f"ERROR: --top-n ({args.top_n}) must be < --max-isoforms ({args.max_isoforms})", file=sys.stderr)
        sys.exit(1)

    # Run filtering
    filter_annotation(
        args.annotation,
        args.gene,
        args.max_isoforms,
        args.top_n,
        args.output,
        args.log_file
    )


if __name__ == '__main__':
    main()
