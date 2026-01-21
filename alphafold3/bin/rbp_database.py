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
POSTAR3 RBP binding site database interface.

Queries tabix-indexed POSTAR3 BED file to find RBP binding sites
overlapping with mutation positions.
"""

import gzip
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Dict, Set
import subprocess
import sys


@dataclass
class RBPBindingSite:
    """Single RBP binding site from POSTAR3."""
    chrom: str
    start: int  # 0-based
    end: int    # exclusive
    rbp_name: str
    strand: str
    method: str
    cell_line: str
    accession: str
    score: float

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2

    def distance_to(self, pos: int) -> int:
        """Distance from position to nearest edge of binding site."""
        if pos < self.start:
            return self.start - pos
        elif pos >= self.end:
            return pos - self.end + 1
        else:
            return 0  # Position is within binding site


class POSTAR3Database:
    """
    Interface to POSTAR3 RBP binding site database.

    Expects either:
    - Tabix-indexed BED file (.bed.gz + .bed.gz.tbi)
    - Plain text file (slower, loaded into memory)

    POSTAR3 format (10 columns):
    chr, start, end, ID, strand, RBP_name, method, cell_line, accession, score
    """

    def __init__(self, postar_path: str, use_tabix: bool = True):
        self.path = Path(postar_path)
        self.use_tabix = use_tabix
        self._in_memory_data: Optional[Dict[str, List[RBPBindingSite]]] = None
        self._tabix_available = False

        if not self.path.exists():
            raise FileNotFoundError(f"POSTAR3 file not found: {self.path}")

        # Check for tabix index
        if use_tabix and str(self.path).endswith('.gz'):
            tbi_path = Path(str(self.path) + '.tbi')
            if tbi_path.exists():
                self._tabix_available = self._check_tabix_binary()
            else:
                print(f"Warning: No tabix index found at {tbi_path}", file=sys.stderr)
                print("Falling back to in-memory mode (slower)", file=sys.stderr)

        if not self._tabix_available:
            self._load_into_memory()

    def _check_tabix_binary(self) -> bool:
        """Check if tabix is available."""
        try:
            result = subprocess.run(['tabix', '--version'],
                                   capture_output=True, text=True)
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def _load_into_memory(self):
        """Load entire POSTAR3 file into memory, indexed by chromosome."""
        print(f"Loading POSTAR3 into memory from {self.path}...", file=sys.stderr)
        self._in_memory_data = {}

        opener = gzip.open if str(self.path).endswith('.gz') else open
        mode = 'rt' if str(self.path).endswith('.gz') else 'r'

        count = 0
        with opener(self.path, mode) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                site = self._parse_line(line)
                if site:
                    if site.chrom not in self._in_memory_data:
                        self._in_memory_data[site.chrom] = []
                    self._in_memory_data[site.chrom].append(site)
                    count += 1

        # Sort by start position for binary search
        for chrom in self._in_memory_data:
            self._in_memory_data[chrom].sort(key=lambda x: x.start)

        print(f"Loaded {count} RBP binding sites", file=sys.stderr)

    def _parse_line(self, line: str) -> Optional[RBPBindingSite]:
        """Parse a single POSTAR3 line."""
        fields = line.split('\t')
        if len(fields) < 10:
            return None

        try:
            return RBPBindingSite(
                chrom=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                rbp_name=fields[5],
                strand=fields[4],
                method=fields[6],
                cell_line=fields[7],
                accession=fields[8],
                score=float(fields[9]) if fields[9] and fields[9] != '.' else 0.0
            )
        except (ValueError, IndexError):
            return None

    def query(self, chrom: str, start: int, end: int) -> List[RBPBindingSite]:
        """
        Query for RBP binding sites overlapping a genomic region.

        Args:
            chrom: Chromosome (e.g., 'chr1' or '1')
            start: Start position (0-based)
            end: End position (exclusive)

        Returns:
            List of overlapping RBP binding sites
        """
        # Normalize chromosome name
        chrom_variants = [chrom, f"chr{chrom}", chrom.replace('chr', '')]

        if self._tabix_available:
            return self._query_tabix(chrom, start, end)
        else:
            return self._query_memory(chrom_variants, start, end)

    def _query_tabix(self, chrom: str, start: int, end: int) -> List[RBPBindingSite]:
        """Query using tabix for fast genomic lookups."""
        results = []

        # Try with and without 'chr' prefix
        for chrom_variant in [chrom, f"chr{chrom}", chrom.replace('chr', '')]:
            try:
                region = f"{chrom_variant}:{start}-{end}"
                result = subprocess.run(
                    ['tabix', str(self.path), region],
                    capture_output=True, text=True
                )

                if result.returncode == 0 and result.stdout:
                    for line in result.stdout.strip().split('\n'):
                        if line:
                            site = self._parse_line(line)
                            if site:
                                results.append(site)
                    if results:
                        break
            except Exception:
                continue

        return results

    def _query_memory(self, chrom_variants: List[str], start: int, end: int) -> List[RBPBindingSite]:
        """Query in-memory data."""
        results = []

        for chrom in chrom_variants:
            if chrom in self._in_memory_data:
                for site in self._in_memory_data[chrom]:
                    # Check overlap
                    if site.end > start and site.start < end:
                        results.append(site)
                if results:
                    break

        return results

    def query_position(self, chrom: str, pos: int, window: int = 50) -> List[RBPBindingSite]:
        """
        Query for RBP binding sites near a single position.

        Args:
            chrom: Chromosome
            pos: Position (0-based)
            window: Window size on each side (default ±50bp)

        Returns:
            List of nearby RBP binding sites
        """
        return self.query(chrom, pos - window, pos + window + 1)

    def get_unique_rbps(self, sites: List[RBPBindingSite]) -> Set[str]:
        """Get unique RBP names from a list of binding sites."""
        return {site.rbp_name for site in sites}

    def group_by_rbp(self, sites: List[RBPBindingSite]) -> Dict[str, List[RBPBindingSite]]:
        """Group binding sites by RBP name."""
        grouped = {}
        for site in sites:
            if site.rbp_name not in grouped:
                grouped[site.rbp_name] = []
            grouped[site.rbp_name].append(site)
        return grouped


def index_postar3(input_path: str, output_path: Optional[str] = None) -> str:
    """
    Sort and index a POSTAR3 file with bgzip and tabix.

    Args:
        input_path: Path to input POSTAR3 file (plain text or gzipped)
        output_path: Output path for indexed file (default: input + .sorted.bed.gz)

    Returns:
        Path to indexed file
    """
    input_path = Path(input_path)

    if output_path is None:
        output_path = input_path.parent / (input_path.stem.replace('.gz', '') + '.sorted.bed.gz')
    else:
        output_path = Path(output_path)

    # Check for required tools
    for tool in ['sort', 'bgzip', 'tabix']:
        try:
            subprocess.run([tool, '--version'], capture_output=True)
        except FileNotFoundError:
            raise RuntimeError(f"Required tool not found: {tool}")

    print(f"Sorting and indexing {input_path}...", file=sys.stderr)

    # Decompress if needed, sort, compress with bgzip
    if str(input_path).endswith('.gz'):
        cat_cmd = ['gzcat' if sys.platform == 'darwin' else 'zcat', str(input_path)]
    else:
        cat_cmd = ['cat', str(input_path)]

    # Pipeline: cat/zcat | sort | bgzip > output
    with open(output_path, 'wb') as outfile:
        cat_proc = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
        sort_proc = subprocess.Popen(
            ['sort', '-k1,1', '-k2,2n'],
            stdin=cat_proc.stdout,
            stdout=subprocess.PIPE
        )
        bgzip_proc = subprocess.Popen(
            ['bgzip', '-c'],
            stdin=sort_proc.stdout,
            stdout=outfile
        )
        bgzip_proc.wait()

    # Index with tabix
    subprocess.run(['tabix', '-p', 'bed', str(output_path)], check=True)

    print(f"Created indexed file: {output_path}", file=sys.stderr)
    return str(output_path)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='POSTAR3 database utilities')
    subparsers = parser.add_subparsers(dest='command')

    # Index command
    index_parser = subparsers.add_parser('index', help='Sort and index POSTAR3 file')
    index_parser.add_argument('input', help='Input POSTAR3 file')
    index_parser.add_argument('-o', '--output', help='Output path')

    # Query command
    query_parser = subparsers.add_parser('query', help='Query POSTAR3 database')
    query_parser.add_argument('database', help='POSTAR3 database file')
    query_parser.add_argument('region', help='Region to query (chr:start-end)')

    args = parser.parse_args()

    if args.command == 'index':
        index_postar3(args.input, args.output)

    elif args.command == 'query':
        db = POSTAR3Database(args.database)

        # Parse region
        chrom, coords = args.region.split(':')
        start, end = map(int, coords.split('-'))

        sites = db.query(chrom, start, end)
        print(f"Found {len(sites)} binding sites:")
        for site in sites:
            print(f"  {site.rbp_name}: {site.chrom}:{site.start}-{site.end} ({site.method})")
