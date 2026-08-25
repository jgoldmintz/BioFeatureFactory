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
BioFeatureFactory: shared primitives.

The lowest layer of biofeaturefactory.lib. Everything here is depended on by
two or more of the sibling modules (msa, codon_metrics, dtu_outputs, annotation)
as well as by utility itself: FASTA reading, alphabet detection, the primary-key
mint, mutation-token helpers, the codon tables and the chromosome-name map.

This module exists to make the package acyclic. The sibling modules were split
out of utility.py but still needed a handful of its primitives, which created a
cycle (utility -> msa -> utility) that only a lazy PEP 562 __getattr__ could
paper over. Hoisting the shared primitives one layer down removes the cycle
outright, so every import in the package is now plain and eager:

    {utility, core}  ->  {msa, codon_metrics, dtu_outputs, annotation}  ->  primitives

Nothing here may import from its siblings or from utility.
"""

import csv
import os
import re
import sys
import hashlib
from pathlib import Path
from typing import Dict, Optional, Tuple

chromosome_map = {
    "GRCh37": {"1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11", "4": "NC_000004.11", "5": "NC_000005.9",
               "6": "NC_000006.11", "7": "NC_000007.13", "8": "NC_000008.10", "9": "NC_000009.11", "10": "NC_000010.10",
               "11": "NC_000011.9", "12": "NC_000012.11", "13": "NC_000013.10", "14": "NC_000014.8",
               "15": "NC_000015.9", "16": "NC_000016.9", "17": "NC_000017.10", "18": "NC_000018.9", "19": "NC_000019.9",
               "20": "NC_000020.10", "21": "NC_000021.8", "22": "NC_000022.10", "X": "NC_000023.10",
               "Y": "NC_000024.9"},
    "GRCh38": {"1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12", "5": "NC_000005.10",
               "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12", "10": "NC_000010.11",
               "11": "NC_000011.10", "12": "NC_000012.12", "13": "NC_000013.11", "14": "NC_000014.9",
               "15": "NC_000015.10", "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
               "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9", "22": "NC_000022.11",
               "X": "NC_000023.11", "Y": "NC_000024.10"}}


codon_table = {'I': ['ATT', 'ATC', 'ATA'],
               'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
               'V': ['GTT', 'GTC', 'GTA', 'GTG'],
               'F': ['TTT', 'TTC'],
               'M': ['ATG'],
               'C': ['TGT', 'TGC'],
               'A': ['GCT', 'GCC', 'GCA', 'GCG'],
               'G': ['GGT', 'GGC', 'GGA', 'GGG'],
               'P': ['CCT', 'CCC', 'CCA', 'CCG'],
               'T': ['ACT', 'ACC', 'ACA', 'ACG'],
               'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
               'Y': ['TAT', 'TAC'],
               'W': ['TGG'],
               'Q': ['CAA', 'CAG'],
               'N': ['AAT', 'AAC'],
               'H': ['CAT', 'CAC'],
               'E': ['GAA', 'GAG'],
               'D': ['GAT', 'GAC'],
               'K': ['AAA', 'AAG'],
               'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
               'Stop': ['TAA', 'TAG', 'TGA'],
               '-': ['---']}


codon_to_aa = {codon: aa for aa, v in codon_table.items() for codon in v}


def read_fasta(inf, aformat="FIRST", duplicate="replace"):
    """Load sequences from a FASTA file into a name->sequence dictionary."""
    data = {}
    with open(inf, "r") as fa:
        name = ""
        for line in fa.readlines():
            if "#" in line:
                continue
            if ">" in line:
                if aformat.upper() == "NCBI":
                    name = re.search(r">[a-zA-Z]+_?\d+(\.\d+)*", line).group(0)
                elif aformat.upper() in ["FIRST", "WORD"]:
                    name = line.split()[0]
                else:
                    name = line.strip()
                name = name[1:].strip()
                if name in data.keys():
                    if duplicate.lower() in ["append", "a"]:  # simply add to existing sequence
                        pass
                    elif duplicate.lower() in ["replace", "r"]:  # reset sequence to empty
                        data[name] = ""
                    elif duplicate.lower() in ["separate", "s"]:  # add underscore+number to end of sequence name
                        matches = re.findall(r"/_\d+$/", name)
                        if matches != None and len(matches) > 0:
                            num = int(max(matches)[1:])
                            name = name[:-len(str(num))] + str(num + 1)
                            data[name] = ""
                        else:
                            name = name + "_2"
                            data[name] = ""
                else:
                    data[name] = ""
            else:
                data[name] = data[name] + line.strip()
    return data


def get_mutation_data_bioAccurate(ntposnt, is_nt):
    """Return one-based position and nucleotides for a mutation string; skips stop codons.

    is_nt is REQUIRED (no default): pass True for a nucleotide mutation token
    (e.g. G123A) or False for an amino-acid token (e.g. R213W). When True, a token
    whose ref/alt are not nucleotides (ACGTU) returns (None, None) instead of
    silently coercing an off-alphabet character into a position/ref/alt -- this
    closes the input-level validation gap (F45). Callers MUST state intent
    explicitly; a missing argument is a TypeError by design (fail-loud, no silent
    wrong default). alphafold3 opts out (is_nt=False) pending its own custom guard.
    """

    # Skip stop codons (for the case of aa)
    if 'Stop' in ntposnt or 'Sto' in ntposnt:
        return None, None

    # gDNA/intronic token. Same reasoning as get_mutation_data above: this is
    # unreachable for any currently-working input. Note the F45 alphabet guard
    # below does NOT catch it -- 'g' and the trailing base are both legal
    # nucleotide letters, so it falls through to int('d.T5000') and raises.
    # Returning (None, None) matches this function's existing contract for
    # "recognised but not scoreable here", which it already uses for stop codons
    # and off-alphabet tokens.
    #
    # This is a floor, not a fix: it converts a crash into a SILENT drop. A
    # pipeline that must report intronic tokens has to gate them at ingest with
    # warn_intronic_unsupported BEFORE reaching here.
    if is_intronic_token(ntposnt):
        return None, None

    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    # F45: reject off-alphabet tokens on nt paths (collapses the former dead if/else).
    if is_nt and (original_nt.upper() not in "ACGTU" or mutant_nt.upper() not in "ACGTU"):
        return None, None
    position = int(ntposnt[1:-1])
    return position, (original_nt, mutant_nt)


PKEY_HASH_HEX = 12


def mint_pkey(gene: str, token: str) -> str:
    """The bounded primary key for one (gene, token) pair: '{GENE}-{sha1[:12]}'.

    The key used to be '{GENE}-{token}', and the token spells out the whole REF
    allele, so key length was the length of the variant. Measured across the
    19,305 genes in grch38.txt, for a whole-gene deletion:

        mean gene span 57,992 nt; largest CNTNAP2 2,304,640 nt (a 2.3 MB key)
        99.3% exceed the 255-byte NAME_MAX that filenames are built against
        89.2% exceed PostgreSQL's ~2,704 B btree limit
        88.0% exceed MySQL InnoDB's 3,072 B index prefix
        22.9% exceed VARCHAR's 65,535 entirely

    Those are not hypothetical: miranda lost 6 runs to '[Errno 63] File name
    too long', and netMHC's workdir_stem already exists to route around the
    same wall. Length here is constant at len(gene) + 13 for every variant
    class, SNV through knockout.

    sha1 is a content digest, not a security primitive. 12 hex = 48 bits; the
    gene travels in the clear, so a collision must occur between two variants
    of the SAME gene, and the population per prefix is one gene's variants
    rather than the whole corpus.

    The digest is over the VERBATIM token -- what every existing cross-pipeline
    join already keys on -- so this changes key LENGTH and nothing else. It
    therefore inherits the property that two textual spellings of one deletion
    remain two keys. canonical_token above collapses those, but applying it
    first needs the per-space sequence for left-alignment, which does not exist
    at the point tokens are read from the mutations CSV.
    """
    digest = hashlib.sha1(token.encode("utf-8")).hexdigest()[:PKEY_HASH_HEX]
    return f"{gene.upper()}-{digest}"


PKEY_HEX_RE = re.compile(r"^[0-9a-f]{%d}$" % PKEY_HASH_HEX)


def load_pkey_map(source, gene: str = None) -> dict:
    """{pkey -> token} from pkey_mapping_<GENE>.csv (file or a directory holding it).

    The inverse of mint_pkey. mint_pkey is deterministic, so going token -> pkey
    never needs a file; this exists for the other direction, which cannot be
    computed and is what a consumer needs when the only thing its tool hands
    back is a sequence name.

    Accepts a file, a directory (searched recursively for pkey_mapping*.csv), or
    None. Returns {} rather than raising when there is nothing to read, so a
    caller can pass whatever it has and fall back to the legacy name format.
    """
    from pathlib import Path as _Path
    out: dict = {}
    if not source:
        return out
    src = _Path(source)
    if src.is_dir():
        files = sorted(src.rglob("pkey_mapping*.csv"))
        if gene:
            named = [f for f in files if gene.upper() in f.stem.upper()]
            files = named or files
    elif src.is_file():
        files = [src]
    else:
        return out
    for f in files:
        try:
            with open(f, newline="") as fh:
                for row in csv.DictReader(fh):
                    pk = (row.get("pkey") or "").strip()
                    tok = (row.get("mutant") or "").strip()
                    if pk and tok:
                        out[pk] = tok
        except Exception:
            continue
    return out


def load_token_pkey_index(source, gene: str = None) -> dict:
    """{token -> pkey}, indexed under BOTH the verbatim token and its upper-case form.

    mint_pkey hashes the token exactly as written, so a consumer holding an
    upper-cased token cannot reproduce the digest of a lower-case one:
    sha('gd.G199C') != sha('GD.G199C'), and the resulting pkey exists in no
    mapping. Several pipelines normalise early -- netNglyc's
    normalize_mutation_id upper-cases on the way in, so by the time a token
    reaches a pkey mint its original spelling is gone.

    Indexing both spellings lets those callers resolve to the SAME pkey
    variant_mapping minted from the verbatim token, without having to thread
    the original casing through. Prefer this over mint_pkey wherever the token
    in hand may already have been normalised.
    """
    out: dict = {}
    for pkey, token in load_pkey_map(source, gene).items():
        out[token] = pkey
        out[token.upper()] = pkey
    return out


def token_from_name(seq_name: str, gene_name: str = None, pkey_map: dict = None) -> str:
    """Recover the mutation token from a sequence name, in EITHER header format.

    Sequence names are the only channel back from the DTU tools -- they return
    predictions keyed by the FASTA header and nothing else -- so every consumer
    has to turn a name into a token. They each did it by hand:
    rsplit('-', 1) in extract_mutation_from_sequence_name, a gene-prefix strip
    in netMHC._process_mutant and NetSurfP3._token_from_key, split('-')[-1] in
    parse_predictions_with_mutation_filtering.

    Handles both:
        'SMN2-C330T'          legacy, token in the clear -> 'C330T'
        'SMN2-052cdd2c9c56'   hashed  -> looked up in pkey_map -> 'C330T'

    Passing gene_name strips that exact prefix, which is what makes hyphenated
    genes safe (NKX2-1, HLA-A, MT-CO1); without it the last hyphen is used,
    which is correct only because tokens contain no hyphen.

    Returns the remainder unchanged when it is not a known digest, so a caller
    wired to this keeps working against legacy headers and needs no second
    change when the mint flips.
    """
    if not isinstance(seq_name, str) or not seq_name:
        return seq_name
    name = seq_name.strip()
    if pkey_map and name in pkey_map:
        return pkey_map[name]
    if gene_name:
        prefix = f"{gene_name}-"
        rest = name[len(prefix):] if name.startswith(prefix) else name
    else:
        rest = name.rsplit("-", 1)[1] if "-" in name else name
    if pkey_map and PKEY_HEX_RE.match(rest.lower()):
        for pk, tok in pkey_map.items():
            if pk.rsplit("-", 1)[-1].lower() == rest.lower():
                return tok
    return rest


INTRONIC_PREFIX = "gd."


CHROM_PREFIX = "ch."


NON_ORF_PREFIXES = (INTRONIC_PREFIX, CHROM_PREFIX)


def is_intronic_token(token) -> bool:
    """True for any NON-ORF-space token: 'gd.T5000C' (gDNA) or 'ch.T70050267A'
    (absolute chromosomal).

    The name is historical -- it predates the chromosomal space and a gd./ch.
    token can perfectly well be exonic. What it actually tests is "is this token
    written in a coordinate space other than the ORF", which is the question
    every caller is really asking: none of them can score a token whose position
    is not an ORF offset, whichever non-ORF space it came from.

    Intronic (and other non-ORF gDNA) variants carry this prefix because a bare
    token is ambiguous between ORF and gDNA space whenever a gene's 5'UTR is
    shorter than its ORF -- SMN2's ORF spans 1-879 and its gDNA slice starts at
    18, so 'T500C' is a legal coordinate in both.

    parse_variant deliberately does NOT understand the prefix and returns None
    for it, so a pipeline that never calls this helper still fails closed. What
    it will NOT do is say anything useful: the legacy parsers raise
    ValueError("invalid literal for int() with base 10: 'd.T5000'"), which
    reads as a corrupt input file rather than as a variant class the pipeline
    cannot score.
    """
    # CASE-INSENSITIVE. Tokens are routinely upper-cased before they reach a gate
    # -- netnglyc's normalize_mutation_id does it to build the pkey, and the same
    # pattern minted 'TESTG~INTRON1' in miranda -- so a case-sensitive prefix test
    # silently stops recognising the very tokens it exists to catch, and the
    # variant resumes being reported as unparseable garbage. Widening is monotone:
    # it can only classify MORE tokens as non-ORF, never fewer, so no caller
    # becomes more permissive.
    return isinstance(token, str) and token.lower().startswith(NON_ORF_PREFIXES)


def ExtractGeneFromFASTA(file_path,count=False):
    """Extract gene name from NetNGlyc output file using read_fasta"""
    sequences = read_fasta(file_path)
    if sequences:
        first_seq_name = list(sequences.keys())[0]
        separators = ['-', '_']
        for sep in separators:
            if sep in first_seq_name:
                if count:
                    return first_seq_name.rsplit(sep, 1)[0],len(sequences)
                return first_seq_name.rsplit(sep, 1)[0]
        return first_seq_name
    return None


def extract_mutation_from_sequence_name(seq_name, gene_name=None, pkey_map=None):
    """Extract (gene, mutation_id) from a sequence name, in EITHER header format.

    'ZFP36-C330T'        -> ('ZFP36', 'C330T')
    'ZFP36-052cdd2c9c56' -> ('ZFP36', 'C330T')   when pkey_map carries it

    The bare two-argument-free call is unchanged, so existing callers keep
    working: without a pkey_map the last hyphen is used, which is correct only
    because tokens contain no hyphen. Passing gene_name strips that exact prefix
    instead and is safe for hyphenated genes (NKX2-1, HLA-A, MT-CO1); passing
    pkey_map additionally resolves a hashed name back to its token, so a caller
    wired to this needs no further edit when the header mint changes.
    """
    if not isinstance(seq_name, str) or not seq_name:
        return seq_name, None
    if gene_name:
        prefix = f"{gene_name}-"
        if seq_name.startswith(prefix):
            return gene_name, token_from_name(seq_name, gene_name, pkey_map)
        return seq_name, None
    if '-' in seq_name:
        gene, rest = seq_name.rsplit('-', 1)
        return gene, token_from_name(seq_name, gene, pkey_map)
    return seq_name, None


def should_skip_mutation(gene, mutation_id, failure_map):
    """Return True if the given gene/mutation should be filtered based on validation logs."""
    if not failure_map or not gene or not mutation_id:
        return False
    return mutation_id.upper() in failure_map.get(gene.upper(), set())


HUMAN_TAI_WEIGHTS = {
    'TTT': 0.344, 'TTC': 1.000, 'TTA': 0.051, 'TTG': 0.344,
    'TCT': 0.344, 'TCC': 0.688, 'TCA': 0.172, 'TCG': 0.086,
    'TAT': 0.344, 'TAC': 1.000, 'TAA': 0.000, 'TAG': 0.000,
    'TGT': 0.344, 'TGC': 1.000, 'TGA': 0.000, 'TGG': 1.000,
    'CTT': 0.172, 'CTC': 0.516, 'CTA': 0.086, 'CTG': 1.000,
    'CCT': 0.344, 'CCC': 0.688, 'CCA': 0.344, 'CCG': 0.172,
    'CAT': 0.344, 'CAC': 1.000, 'CAA': 0.344, 'CAG': 1.000,
    'CGT': 0.086, 'CGC': 0.344, 'CGA': 0.086, 'CGG': 0.172,
    'ATT': 0.516, 'ATC': 1.000, 'ATA': 0.086, 'ATG': 1.000,
    'ACT': 0.344, 'ACC': 1.000, 'ACA': 0.344, 'ACG': 0.172,
    'AAT': 0.344, 'AAC': 1.000, 'AAA': 0.344, 'AAG': 1.000,
    'AGT': 0.172, 'AGC': 1.000, 'AGA': 0.172, 'AGG': 0.172,
    'GTT': 0.344, 'GTC': 0.688, 'GTA': 0.172, 'GTG': 1.000,
    'GCT': 0.516, 'GCC': 1.000, 'GCA': 0.344, 'GCG': 0.172,
    'GAT': 0.344, 'GAC': 1.000, 'GAA': 0.344, 'GAG': 1.000,
    'GGT': 0.344, 'GGC': 1.000, 'GGA': 0.344, 'GGG': 0.344,
    '---': 0.000,
}


HUMAN_REFERENCE_W = {
    'TTT': 0.45, 'TTC': 1.00, 'TTA': 0.07, 'TTG': 0.13,
    'TCT': 0.44, 'TCC': 0.53, 'TCA': 0.26, 'TCG': 0.11,
    'TAT': 0.43, 'TAC': 1.00, 'TAA': 1.00, 'TAG': 0.23,
    'TGT': 0.45, 'TGC': 1.00, 'TGA': 0.47, 'TGG': 1.00,
    'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 1.00,
    'CCT': 0.52, 'CCC': 0.63, 'CCA': 0.51, 'CCG': 0.18,
    'CAT': 0.41, 'CAC': 1.00, 'CAA': 0.25, 'CAG': 1.00,
    'CGT': 0.18, 'CGC': 0.43, 'CGA': 0.14, 'CGG': 0.25,
    'ATT': 0.36, 'ATC': 1.00, 'ATA': 0.16, 'ATG': 1.00,
    'ACT': 0.37, 'ACC': 1.00, 'ACA': 0.42, 'ACG': 0.18,
    'AAT': 0.46, 'AAC': 1.00, 'AAA': 0.42, 'AAG': 1.00,
    'AGT': 0.29, 'AGC': 1.00, 'AGA': 0.45, 'AGG': 0.42,
    'GTT': 0.18, 'GTC': 0.29, 'GTA': 0.11, 'GTG': 1.00,
    'GCT': 0.45, 'GCC': 1.00, 'GCA': 0.38, 'GCG': 0.19,
    'GAT': 0.46, 'GAC': 1.00, 'GAA': 0.42, 'GAG': 1.00,
    'GGT': 0.35, 'GGC': 1.00, 'GGA': 0.46, 'GGG': 0.35,
    '---': 0.00,
}


_CODON_ENCODE_CHARS = (
    [chr(c) for c in range(ord('A'), ord('Z') + 1)]
    + [chr(c) for c in range(ord('a'), ord('z') + 1)]
    + [str(d) for d in range(10)]
    + ['!', '@']
)

CODON_ENCODED_ALPHABET = '-' + ''.join(_CODON_ENCODE_CHARS)


_CODON_ONLY_MARKERS = frozenset('0123456789!@')


def detect_alphabet(sequence):
    """Return 'codon', 'nucleotide', or 'protein' based on character composition.

    - 'codon'     : single-character codon-encoded sequence (see
                    mutation_effects/bin/codon_encoding.py). Reported only when
                    the sequence carries a codon-only marker (a digit, '!' or
                    '@') AND every non-gap character lies in the codon alphabet.
                    Evaluated case-sensitively and BEFORE the nt/protein test,
                    because the codon alphabet is case-sensitive and its digits
                    would otherwise be counted as non-nucleotide.
                    Limitation: a codon-encoded sequence that happens to use only
                    upper-case-letter codons (no digit/'!'/'@') is
                    indistinguishable from protein by composition and is reported
                    as 'protein'. A hard limit of composition-only detection.
    - 'nucleotide': >=90% of non-gap characters are ACGT + U (RNA) + N (masked
                    base). IUPAC ambiguity codes are excluded -- they collide with
                    amino-acid letters and do not appear in BFF's ORF/CDS inputs.
    - 'protein'   : otherwise.
    """
    raw = sequence.replace('-', '').replace('.', '')
    if not raw.replace('*', ''):
        raise ValueError("Empty sequence.")

    # Codon-encoded check first, case-sensitive.
    codon_alpha = set(CODON_ENCODED_ALPHABET)
    if any(c in _CODON_ONLY_MARKERS for c in raw) and all(c in codon_alpha for c in raw):
        return 'codon'

    seq = raw.upper().replace('*', '')
    # Bases + U (RNA) + N (masked base) only. IUPAC ambiguity codes
    # (RYSWKMBDHV) are intentionally EXCLUDED: 10 of them collide with amino-acid
    # one-letter codes, weakening nt-vs-protein discrimination, and they do not
    # occur in BFF's ORF/CDS inputs (the only sequences detect_alphabet runs on).
    nt_chars = set('ACGTUN')
    nt_count = sum(1 for c in seq if c in nt_chars)
    return 'nucleotide' if nt_count / len(seq) >= 0.90 else 'protein'

