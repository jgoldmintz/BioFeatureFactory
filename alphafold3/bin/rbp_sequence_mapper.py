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
RBP name to protein sequence mapper.

Maps RBP gene names to UniProt IDs and retrieves protein sequences
and MSAs from local files (AlphaFold A3M format).
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Tuple, List
from dataclasses import dataclass, field


@dataclass
class RBPSequence:
    """RBP protein sequence and MSA data."""
    gene_name: str
    uniprot_id: str
    sequence: str
    length: int
    protein_name: Optional[str] = None
    msa_path: Optional[str] = None
    msa_content: Optional[str] = None
    msa_depth: int = 0

    @classmethod
    def from_fasta(cls, header: str, sequence: str, gene_name: str) -> 'RBPSequence':
        """Create from FASTA header and sequence."""
        # Parse UniProt FASTA header: >sp|P12345|GENE_HUMAN Description
        parts = header.split('|')
        if len(parts) >= 2:
            uniprot_id = parts[1]
        else:
            uniprot_id = header.split()[0].lstrip('>')

        protein_name = None
        if ' ' in header:
            protein_name = header.split(' ', 1)[1].split(' OS=')[0]

        return cls(
            gene_name=gene_name,
            uniprot_id=uniprot_id,
            sequence=sequence,
            length=len(sequence),
            protein_name=protein_name
        )

    @classmethod
    def from_a3m(cls, a3m_path: str, gene_name: str, uniprot_id: str) -> 'RBPSequence':
        """Create from A3M MSA file. Extracts query sequence (first entry)."""
        with open(a3m_path, 'r') as f:
            content = f.read()

        lines = content.strip().split('\n')

        # First sequence is the query
        query_seq = None
        msa_depth = 0

        i = 0
        while i < len(lines):
            if lines[i].startswith('>'):
                msa_depth += 1
                if query_seq is None:
                    # Collect query sequence (may span multiple lines)
                    seq_parts = []
                    i += 1
                    while i < len(lines) and not lines[i].startswith('>'):
                        seq_parts.append(lines[i].strip())
                        i += 1
                    # Remove gaps from query sequence
                    query_seq = ''.join(seq_parts).replace('-', '').replace('.', '')
                    continue
            i += 1

        if not query_seq:
            raise ValueError(f"No sequence found in {a3m_path}")

        return cls(
            gene_name=gene_name,
            uniprot_id=uniprot_id,
            sequence=query_seq,
            length=len(query_seq),
            msa_path=a3m_path,
            msa_content=content,
            msa_depth=msa_depth
        )


class RBPSequenceMapper:
    """
    Maps RBP gene names to protein sequences and MSAs.

    Uses:
    1. Gene name → UniProt ID mapping file (TSV)
    2. Local FASTA files for protein sequences OR
    3. A3M MSA files (preferred - includes sequence + alignment)
    """

    def __init__(
        self,
        mapping_file: str,
        sequence_dir: Optional[str] = None,
        sequence_fasta: Optional[str] = None,
        msa_dir: Optional[str] = None
    ):
        """
        Initialize the mapper.

        Args:
            mapping_file: Path to gene name → UniProt mapping TSV
                          (columns: Entry, Gene Names, Protein names, Length)
            sequence_dir: Directory containing protein FASTA files
            sequence_fasta: Single FASTA file with all protein sequences
            msa_dir: Directory containing A3M MSA files (AF-{uniprot}-F1-msa_v6.a3m)
        """
        self.mapping_file = Path(mapping_file)
        self.sequence_dir = Path(sequence_dir) if sequence_dir else None
        self.sequence_fasta = Path(sequence_fasta) if sequence_fasta else None
        self.msa_dir = Path(msa_dir) if msa_dir else None

        # Gene name (uppercase) -> UniProt ID
        self._gene_to_uniprot: Dict[str, str] = {}
        # UniProt ID -> full mapping info
        self._uniprot_info: Dict[str, dict] = {}
        # UniProt ID -> sequence (lazy loaded)
        self._sequences: Dict[str, str] = {}
        # Gene name (uppercase) -> sequence (direct lookup)
        self._gene_to_sequence: Dict[str, str] = {}
        # UniProt ID -> MSA file path
        self._msa_paths: Dict[str, Path] = {}
        # UniProt ID -> RBPSequence with MSA (cached)
        self._rbp_cache: Dict[str, RBPSequence] = {}

        self._load_mapping()

        # Index MSA directory first (preferred source)
        if self.msa_dir and self.msa_dir.exists():
            self._index_msa_dir()
        elif self.sequence_fasta and self.sequence_fasta.exists():
            self._load_sequences_from_fasta()
        elif self.sequence_dir and self.sequence_dir.exists():
            self._index_sequence_dir()

    def _load_mapping(self):
        """Load gene name to UniProt ID mapping."""
        if not self.mapping_file.exists():
            raise FileNotFoundError(f"Mapping file not found: {self.mapping_file}")

        print(f"Loading UniProt mapping from {self.mapping_file}...", file=sys.stderr)

        with open(self.mapping_file, 'r') as f:
            header = f.readline().strip().split('\t')

            # Find column indices
            entry_idx = header.index('Entry') if 'Entry' in header else 0
            genes_idx = header.index('Gene Names') if 'Gene Names' in header else 1
            name_idx = header.index('Protein names') if 'Protein names' in header else 2
            len_idx = header.index('Length') if 'Length' in header else 3

            for line in f:
                fields = line.strip().split('\t')
                if len(fields) <= max(entry_idx, genes_idx):
                    continue

                uniprot_id = fields[entry_idx]
                gene_names = fields[genes_idx] if genes_idx < len(fields) else ''
                protein_name = fields[name_idx] if name_idx < len(fields) else ''
                length = fields[len_idx] if len_idx < len(fields) else ''

                # Store UniProt info
                self._uniprot_info[uniprot_id] = {
                    'gene_names': gene_names,
                    'protein_name': protein_name,
                    'length': length
                }

                # Map all gene name variants to UniProt ID
                for gene in gene_names.split():
                    gene_upper = gene.upper()
                    # Prefer shorter UniProt IDs (canonical entries)
                    if gene_upper not in self._gene_to_uniprot or \
                       len(uniprot_id) < len(self._gene_to_uniprot[gene_upper]):
                        self._gene_to_uniprot[gene_upper] = uniprot_id

        print(f"Loaded {len(self._gene_to_uniprot)} gene name mappings", file=sys.stderr)

    def _load_sequences_from_fasta(self):
        """Load all sequences from a single FASTA file."""
        print(f"Loading sequences from {self.sequence_fasta}...", file=sys.stderr)

        current_header = None
        current_seq = []

        with open(self.sequence_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_header and current_seq:
                        self._store_sequence(current_header, ''.join(current_seq))
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)

            # Save last sequence
            if current_header and current_seq:
                self._store_sequence(current_header, ''.join(current_seq))

        print(f"Loaded {len(self._sequences)} protein sequences", file=sys.stderr)

    def _store_sequence(self, header: str, sequence: str):
        """Store a sequence, indexing by UniProt ID and gene names."""
        # Parse header for UniProt ID
        # Format: >sp|P12345|GENE_HUMAN ... or >P12345 ...
        parts = header.split('|')
        if len(parts) >= 2:
            uniprot_id = parts[1]
        else:
            uniprot_id = header.split()[0].lstrip('>')

        self._sequences[uniprot_id] = sequence

        # Also index by gene name if we can extract it
        if len(parts) >= 3:
            # >sp|P12345|HNRNPA1_HUMAN -> HNRNPA1
            gene_part = parts[2].split('_')[0]
            self._gene_to_sequence[gene_part.upper()] = sequence

    def _index_sequence_dir(self):
        """Index FASTA files in the sequence directory."""
        print(f"Indexing sequences from {self.sequence_dir}...", file=sys.stderr)

        fasta_extensions = ['.fasta', '.fa', '.faa']
        for fasta_file in self.sequence_dir.iterdir():
            if fasta_file.suffix.lower() in fasta_extensions:
                # Use filename as gene name hint
                gene_hint = fasta_file.stem.upper()

                with open(fasta_file, 'r') as f:
                    header = None
                    seq_parts = []
                    for line in f:
                        line = line.strip()
                        if line.startswith('>'):
                            if header and seq_parts:
                                self._store_sequence(header, ''.join(seq_parts))
                            header = line
                            seq_parts = []
                        else:
                            seq_parts.append(line)
                    if header and seq_parts:
                        self._store_sequence(header, ''.join(seq_parts))

        print(f"Indexed {len(self._sequences)} sequences from directory", file=sys.stderr)

    def _index_msa_dir(self):
        """Index A3M MSA files in the MSA directory."""
        print(f"Indexing MSAs from {self.msa_dir}...", file=sys.stderr)

        # Pattern: AF-{uniprot_id}-F1-msa_v6.a3m
        for msa_file in self.msa_dir.glob("AF-*-F1-msa*.a3m"):
            # Extract UniProt ID from filename
            # AF-P09651-F1-msa_v6.a3m -> P09651
            parts = msa_file.stem.split('-')
            if len(parts) >= 2:
                uniprot_id = parts[1]
                self._msa_paths[uniprot_id] = msa_file

        print(f"Indexed {len(self._msa_paths)} MSA files", file=sys.stderr)

    def get_uniprot_id(self, gene_name: str) -> Optional[str]:
        """Get UniProt ID for a gene name."""
        return self._gene_to_uniprot.get(gene_name.upper())

    def get_sequence(self, gene_name: str) -> Optional[str]:
        """
        Get protein sequence for an RBP by gene name.

        Args:
            gene_name: RBP gene symbol (e.g., 'HNRNPA1')

        Returns:
            Protein sequence string, or None if not found
        """
        gene_upper = gene_name.upper()

        # Try direct gene name lookup
        if gene_upper in self._gene_to_sequence:
            return self._gene_to_sequence[gene_upper]

        # Try via UniProt ID
        uniprot_id = self._gene_to_uniprot.get(gene_upper)
        if uniprot_id and uniprot_id in self._sequences:
            return self._sequences[uniprot_id]

        return None

    def get_rbp_data(self, gene_name: str, load_msa: bool = True) -> Optional[RBPSequence]:
        """
        Get full RBP data including sequence and optionally MSA.

        Args:
            gene_name: RBP gene symbol
            load_msa: Whether to load full MSA content (default True)

        Returns:
            RBPSequence object with sequence and MSA, or None if not found
        """
        gene_upper = gene_name.upper()
        uniprot_id = self.get_uniprot_id(gene_name)

        # Check cache first
        cache_key = uniprot_id or gene_upper
        if cache_key in self._rbp_cache:
            return self._rbp_cache[cache_key]

        # Try loading from MSA file (preferred - has both sequence and MSA)
        if uniprot_id and uniprot_id in self._msa_paths:
            try:
                rbp_data = RBPSequence.from_a3m(
                    str(self._msa_paths[uniprot_id]),
                    gene_name=gene_upper,
                    uniprot_id=uniprot_id
                )
                self._rbp_cache[cache_key] = rbp_data
                return rbp_data
            except Exception as e:
                print(f"Error loading MSA for {gene_name}: {e}", file=sys.stderr)

        # Fallback to sequence-only
        sequence = self.get_sequence(gene_name)
        if not sequence:
            return None

        if not uniprot_id:
            uniprot_id = 'UNKNOWN'

        info = self._uniprot_info.get(uniprot_id, {})

        rbp_data = RBPSequence(
            gene_name=gene_upper,
            uniprot_id=uniprot_id,
            sequence=sequence,
            length=len(sequence),
            protein_name=info.get('protein_name')
        )
        self._rbp_cache[cache_key] = rbp_data
        return rbp_data

    def get_msa(self, gene_name: str) -> Optional[str]:
        """
        Get MSA content for an RBP.

        Args:
            gene_name: RBP gene symbol

        Returns:
            MSA content as string (A3M format), or None if not available
        """
        rbp_data = self.get_rbp_data(gene_name, load_msa=True)
        if rbp_data:
            return rbp_data.msa_content
        return None

    def has_msa(self, gene_name: str) -> bool:
        """Check if MSA is available for an RBP."""
        uniprot_id = self.get_uniprot_id(gene_name)
        return uniprot_id is not None and uniprot_id in self._msa_paths

    def has_sequence(self, gene_name: str) -> bool:
        """Check if sequence is available for an RBP."""
        return self.get_sequence(gene_name) is not None

    def get_missing_rbps(self, rbp_names: list) -> list:
        """Get list of RBPs without available sequences."""
        return [name for name in rbp_names if not self.has_sequence(name)]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='RBP sequence mapper')
    parser.add_argument('--mapping', required=True, help='Gene-UniProt mapping TSV')
    parser.add_argument('--sequences', help='FASTA file with protein sequences')
    parser.add_argument('--seq-dir', help='Directory with FASTA files')
    parser.add_argument('--msa-dir', help='Directory with A3M MSA files')
    parser.add_argument('--query', nargs='+', help='Gene names to query')

    args = parser.parse_args()

    mapper = RBPSequenceMapper(
        mapping_file=args.mapping,
        sequence_fasta=args.sequences,
        sequence_dir=args.seq_dir,
        msa_dir=args.msa_dir
    )

    if args.query:
        for gene in args.query:
            data = mapper.get_rbp_data(gene)
            if data:
                print(f"{gene}: {data.uniprot_id} ({data.length} aa)")
                print(f"  Sequence: {data.sequence[:50]}...")
                if data.msa_content:
                    print(f"  MSA: {data.msa_depth} sequences, {len(data.msa_content)} bytes")
                else:
                    print(f"  MSA: not available")
            else:
                print(f"{gene}: NOT FOUND")
