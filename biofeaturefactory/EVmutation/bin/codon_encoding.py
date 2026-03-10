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
Codon encoding utilities for plmc-compatible single-character codon MSA.

Maps each of the 64 DNA codons to a unique printable ASCII character,
enabling plmc to process codon-level MSAs with a 65-character alphabet
(gap + 64 codons).

Canonical lexicographic ordering (ACGT): AAA→A, AAC→B, ..., TTT→@
"""

import itertools

# Canonical lexicographic ordering of all 64 codons (ACGT base order)
_BASES = "ACGT"
_ALL_CODONS = ["".join(t) for t in itertools.product(_BASES, repeat=3)]

# 64 unique printable characters: A-Z (26) + a-z (26) + 0-9 (10) + ! @ (2)
_CHARS = (
    [chr(c) for c in range(ord('A'), ord('Z') + 1)] +  # A-Z: 26
    [chr(c) for c in range(ord('a'), ord('z') + 1)] +  # a-z: 26
    [str(d) for d in range(10)] +                       # 0-9: 10
    ['!', '@']                                          # !@:   2  → total 64
)

assert len(_ALL_CODONS) == 64
assert len(_CHARS) == 64

CODON_TO_CHAR = dict(zip(_ALL_CODONS, _CHARS))
CHAR_TO_CODON = dict(zip(_CHARS, _ALL_CODONS))

GAP_CHAR = '-'
GAP_CODON = '---'

# plmc -a alphabet string: gap first, then all 64 codon characters
CODON_ALPHABET = GAP_CHAR + "".join(_CHARS)


def encode_codon_msa(input_path, output_path):
    """
    Convert a triplet-codon FASTA MSA to single-character format for plmc.

    Each codon triplet (e.g. ATG) is mapped to a unique single character.
    Gap triplets (---) and any non-standard triplets map to the gap char '-'.

    Args:
        input_path: Path to codon MSA FASTA with triplet codons per column.
        output_path: Path to write single-character encoded FASTA.

    Returns:
        int: Number of sequences encoded.
    """
    sequences = {}
    current_id = None
    current_seq_parts = []

    with open(input_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq_parts)
                current_id = line[1:]
                current_seq_parts = []
            else:
                current_seq_parts.append(line.strip())

    if current_id is not None:
        sequences[current_id] = ''.join(current_seq_parts)

    with open(output_path, 'w') as out:
        for seq_id, seq in sequences.items():
            seq = seq.replace(' ', '').upper()

            if len(seq) % 3 != 0:
                raise ValueError(
                    f"Sequence '{seq_id}' length {len(seq)} is not divisible by 3"
                )

            encoded = []
            for i in range(0, len(seq), 3):
                triplet = seq[i:i + 3]
                if triplet == GAP_CODON:
                    encoded.append(GAP_CHAR)
                else:
                    char = CODON_TO_CHAR.get(triplet, GAP_CHAR)
                    encoded.append(char)

            out.write(f'>{seq_id}\n')
            out.write(''.join(encoded) + '\n')

    return len(sequences)
