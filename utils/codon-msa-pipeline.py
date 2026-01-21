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
BioFeatureFactory: Codon-Aware MSA Generator

Generates codon-aware multiple sequence alignments from protein alignments
and corresponding nucleotide sequences.

Workflow:
  1. Load protein MSA (from MUSCLE/MAFFT/etc.)
  2. Load corresponding nucleotide sequences for each sequence ID
  3. Back-translate protein alignment to codons
  4. Output codon MSA preserving alignment structure

Output format: FASTA with codon sequences (gaps as '---')
"""

import argparse
import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.utility import (
    read_fasta,
    write_fasta,
    codon_to_aa,
    codon_table,
)


def translate_codon(codon):
    """Translate a codon to amino acid."""
    if codon == '---' or '-' in codon:
        return '-'
    return codon_to_aa.get(codon.upper(), 'X')


def backtranslate_protein_to_codons(protein_seq, nucleotide_seq, seq_id=None):
    """
    Back-translate an aligned protein sequence to codons using the nucleotide sequence.

    Args:
        protein_seq: Aligned protein sequence (with gaps '-')
        nucleotide_seq: Unaligned nucleotide sequence (ORF)
        seq_id: Optional sequence ID for error messages

    Returns:
        str: Codon sequence with gaps as '---'

    Raises:
        ValueError: If sequences don't match
    """
    codon_seq = []
    nt_pos = 0  # Position in nucleotide sequence

    for i, aa in enumerate(protein_seq):
        if aa == '-':
            # Gap in protein alignment -> '---' in codon alignment
            codon_seq.append('---')
        else:
            # Get codon from nucleotide sequence
            if nt_pos + 3 > len(nucleotide_seq):
                raise ValueError(
                    f"Nucleotide sequence too short for protein at position {i} "
                    f"(seq: {seq_id or 'unknown'})"
                )

            codon = nucleotide_seq[nt_pos:nt_pos + 3]
            translated_aa = translate_codon(codon)

            # Verify translation matches
            if translated_aa.upper() != aa.upper() and aa.upper() != 'X':
                raise ValueError(
                    f"Translation mismatch at position {i}: "
                    f"protein has '{aa}', codon '{codon}' translates to '{translated_aa}' "
                    f"(seq: {seq_id or 'unknown'})"
                )

            codon_seq.append(codon)
            nt_pos += 3

    # Check we used all nucleotides (allowing for partial last codon)
    remaining = len(nucleotide_seq) - nt_pos
    if remaining > 3:
        print(f"Warning: {remaining} unused nucleotides at end of {seq_id or 'sequence'}",
              file=sys.stderr)

    return ''.join(codon_seq)


def load_nucleotide_sequences(nt_dir, seq_ids):
    """
    Load nucleotide sequences from a directory of FASTA files.

    Expects one FASTA file per sequence, or a single FASTA with all sequences.

    Args:
        nt_dir: Directory containing nucleotide FASTA files
        seq_ids: List of sequence IDs to load

    Returns:
        dict: {seq_id: nucleotide_sequence}
    """
    nt_seqs = {}
    nt_path = Path(nt_dir)

    # Try loading from a single file first
    combined_fastas = list(nt_path.glob('*.fasta')) + list(nt_path.glob('*.fa')) + list(nt_path.glob('*.fna'))

    for fasta_file in combined_fastas:
        seqs = read_fasta(str(fasta_file))
        for header, seq in seqs.items():
            # Match by exact ID or ID before first space/underscore
            base_id = header.split()[0].split('_')[0]
            if header in seq_ids:
                nt_seqs[header] = seq
            elif base_id in seq_ids:
                nt_seqs[base_id] = seq
            else:
                # Try fuzzy match
                for sid in seq_ids:
                    if sid in header or header in sid:
                        nt_seqs[sid] = seq
                        break

    return nt_seqs


def generate_codon_msa(protein_msa_path, nucleotide_dir, focus_seq=None, output_path=None):
    """
    Generate a codon-aware MSA from a protein MSA and nucleotide sequences.

    Args:
        protein_msa_path: Path to protein MSA FASTA
        nucleotide_dir: Directory or file containing nucleotide sequences
        focus_seq: Optional focus sequence ID to prioritize
        output_path: Output path (default: input_codon.msa.fasta)

    Returns:
        dict: {seq_id: codon_sequence}
    """
    # Load protein MSA
    protein_msa = read_fasta(protein_msa_path)
    seq_ids = list(protein_msa.keys())

    print(f"Loaded protein MSA with {len(seq_ids)} sequences")

    # Load nucleotide sequences
    if os.path.isfile(nucleotide_dir):
        nt_seqs = read_fasta(nucleotide_dir)
    else:
        nt_seqs = load_nucleotide_sequences(nucleotide_dir, seq_ids)

    print(f"Loaded {len(nt_seqs)} nucleotide sequences")

    # Process focus sequence first if specified
    if focus_seq and focus_seq in seq_ids:
        seq_ids.remove(focus_seq)
        seq_ids.insert(0, focus_seq)

    # Back-translate each sequence
    codon_msa = {}
    failed = []

    for seq_id in seq_ids:
        protein_seq = protein_msa[seq_id]

        # Find matching nucleotide sequence
        nt_seq = None
        for nt_id, seq in nt_seqs.items():
            if nt_id == seq_id or seq_id in nt_id or nt_id in seq_id:
                nt_seq = seq
                break

        if nt_seq is None:
            print(f"Warning: No nucleotide sequence found for '{seq_id}'", file=sys.stderr)
            failed.append(seq_id)
            continue

        try:
            codon_seq = backtranslate_protein_to_codons(protein_seq, nt_seq, seq_id)
            codon_msa[seq_id] = codon_seq
        except ValueError as e:
            print(f"Warning: {e}", file=sys.stderr)
            failed.append(seq_id)

    print(f"Successfully processed {len(codon_msa)} sequences, {len(failed)} failed")

    # Write output
    if output_path is None:
        base = Path(protein_msa_path).stem
        output_path = str(Path(protein_msa_path).parent / f"{base}_codon.msa.fasta")

    write_fasta(output_path, codon_msa)
    print(f"Wrote codon MSA to {output_path}")

    return codon_msa


def validate_codon_msa(codon_msa, protein_msa):
    """
    Validate that codon MSA translates back to protein MSA.

    Args:
        codon_msa: {seq_id: codon_sequence}
        protein_msa: {seq_id: protein_sequence}

    Returns:
        list: List of (seq_id, position, expected_aa, found_aa) mismatches
    """
    mismatches = []

    for seq_id in codon_msa:
        if seq_id not in protein_msa:
            continue

        codon_seq = codon_msa[seq_id]
        protein_seq = protein_msa[seq_id]

        # Should have same number of positions (codons vs AAs)
        n_codons = len(codon_seq) // 3
        if n_codons != len(protein_seq):
            mismatches.append((seq_id, -1, f"len={len(protein_seq)}", f"len={n_codons}"))
            continue

        for i in range(n_codons):
            codon = codon_seq[i*3:(i+1)*3]
            expected_aa = protein_seq[i]
            found_aa = translate_codon(codon)

            if expected_aa == '-' and codon != '---':
                mismatches.append((seq_id, i, '-', codon))
            elif expected_aa != '-' and found_aa.upper() != expected_aa.upper():
                if expected_aa.upper() != 'X':  # Allow X as wildcard
                    mismatches.append((seq_id, i, expected_aa, found_aa))

    return mismatches


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Codon-Aware MSA Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From directory of nucleotide FASTA files
  python codon-msa-pipeline.py \\
    --protein-msa /path/to/protein.msa.fasta \\
    --nucleotide-dir /path/to/nt_sequences/ \\
    --focus GENE_FOCUS_ID \\
    --output /path/to/codon.msa.fasta

  # From single nucleotide FASTA file
  python codon-msa-pipeline.py \\
    --protein-msa /path/to/protein.msa.fasta \\
    --nucleotide-file /path/to/all_sequences.fasta \\
    --output /path/to/codon.msa.fasta

Notes:
  - Protein MSA should be from MUSCLE/MAFFT/etc.
  - Nucleotide sequences should be unaligned ORFs
  - Sequence IDs must match between protein MSA and nucleotide files
"""
    )

    parser.add_argument('--protein-msa', required=True,
                        help='Path to protein MSA FASTA file')

    nt_group = parser.add_mutually_exclusive_group(required=True)
    nt_group.add_argument('--nucleotide-dir',
                          help='Directory containing nucleotide FASTA files')
    nt_group.add_argument('--nucleotide-file',
                          help='Single FASTA file with all nucleotide sequences')

    parser.add_argument('--focus', help='Focus sequence ID (placed first in output)')
    parser.add_argument('--output', '-o', help='Output FASTA path (default: input_codon.msa.fasta)')
    parser.add_argument('--validate', action='store_true',
                        help='Validate codon MSA translates correctly')

    args = parser.parse_args()

    nt_source = args.nucleotide_dir or args.nucleotide_file

    # Generate codon MSA
    codon_msa = generate_codon_msa(
        args.protein_msa,
        nt_source,
        args.focus,
        args.output
    )

    # Optional validation
    if args.validate:
        print("\nValidating codon MSA...")
        protein_msa = read_fasta(args.protein_msa)
        mismatches = validate_codon_msa(codon_msa, protein_msa)

        if mismatches:
            print(f"Found {len(mismatches)} mismatches:")
            for seq_id, pos, expected, found in mismatches[:10]:
                print(f"  {seq_id} pos {pos}: expected '{expected}', found '{found}'")
            if len(mismatches) > 10:
                print(f"  ... and {len(mismatches) - 10} more")
        else:
            print("Validation passed: all codons translate correctly")


if __name__ == "__main__":
    main()
