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

import csv
import hashlib
import re
import os
import math
import tempfile
import subprocess
import shutil
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from urllib.parse import unquote
from Bio.Seq import Seq
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

w2n = {'zero': '0', 'one': '1', 'two': '2', 'three': '3', 'four': '4', 'five': '5', 'six': '6', 'seven': '7',
       'eight': '8', 'nine': '9', "ten": '10', "eleven": '11', "twelve": '12', 'thirteen': '13', 'fourteen': '14',
       'fifteen': '15', 'sixteen': '16', 'seventeen': '17', 'eighteen': '18', 'nineteen': '19', 'twenty': '20',
       'twenty-one': '21', 'twenty-two': '22'}

STOPWORDS = {
    "combined","merged","processed","final","updated","new","old","temp","test",
    "transcript","protein","genomic","chromosome","mapping","mutations","sequences",
    "variants","data","aa","nt","cds","orf","map","maps","reference","ref",
    "transcripts","features","tables","table","counts","count","results","out","in",
    "notes","draft","report","summary","v","ver","version"
}

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

def write_fasta(path: Path, name_to_seq: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for name, seq in name_to_seq.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def update_str(s, c, pos):
    return (s[:pos] + c + s[pos + 1:])

def subseq(seq: str, pos: int, l: int) -> str:
    """
    Return a window centered on pos.
    - pos: 0-based index into seq
    - l  : odd window length (e.g., 151)
    Truncates at sequence ends (no padding).
    """
    assert isinstance(l, int) and l > 0 and (l % 2 == 1), "l must be a positive odd integer"
    assert 0 <= pos < len(seq), "pos out of range"
    half = l // 2
    start = max(0, pos - half)
    end   = min(len(seq), pos + half + 1)  # +1 because end is exclusive
    return seq[start:end]

_FILTER_LOG_CACHE: dict[tuple[str, ...], dict[str, set[str]]] = {}
# Exon-aware validation lines look like "GENE: mutation A123G expects ..."
_LOG_MUTATION_RE = re.compile(r"^(?P<gene>[^:]+): mutation (?P<mut>[ACGT][0-9]+[ACGT])\b")


def _normalize_logs(log):
    """Return a flattened list of log file paths, expanding directories."""
    if not log:
        return []
    if isinstance(log, (str, Path)):
        items = [Path(log)]
    else:
        items = [Path(p) for p in log]
    normalized: list[Path] = []
    for item in items:
        if item is None:
            continue
        item = Path(item).expanduser()
        if item.is_dir():
            # If the user passed a directory, gather every *.log inside it
            normalized.extend(sorted(item.glob("*.log")))
        else:
            normalized.append(item)
    return normalized


def _collect_failures_from_logs(log):
    """Parse exon-aware validation logs into a mapping of failing mutations."""
    paths = _normalize_logs(log)
    if not paths:
        return {}
    cache_key = tuple(sorted(str(p) for p in paths))
    if cache_key in _FILTER_LOG_CACHE:
        return _FILTER_LOG_CACHE[cache_key]

    failures: dict[str, set[str]] = {}
    for log_path in paths:
        if not log_path.exists() or not log_path.is_file():
            continue
        try:
            with open(log_path, "r") as handle:
                for line in handle:
                    match = _LOG_MUTATION_RE.match(line.strip())
                    if not match:
                        continue
                    gene = match.group("gene").strip().upper()
                    mut = match.group("mut").strip()
                    if not gene or not mut:
                        continue
                    # Track the failing mutation for this gene in upper-case form
                    failures.setdefault(gene, set()).add(mut)
        except OSError:
            continue

    _FILTER_LOG_CACHE[cache_key] = failures
    return failures


def trim_muts(ntPosnt, log=None, gene_name=None):
    """Read a mutation CSV and optionally drop entries flagged in validation logs."""
    mut_list: list[str] = []
    with open(ntPosnt, 'r') as inf:
        for idx, line in enumerate(inf):
            if idx == 0:
                continue
            cleaned = line.replace('*', '').strip()
            if cleaned:
                # Keep order while removing stray markers (e.g., trailing '*')
                mut_list.append(cleaned)

    if not log:
        return mut_list

    gene = gene_name or extract_gene_from_filename(str(ntPosnt))
    if not gene:
        return mut_list

    failures = _collect_failures_from_logs(log)
    if not failures:
        return mut_list

    skip_set = failures.get(gene.upper())
    if not skip_set:
        return mut_list

    # Only return mutations that are not flagged in the failure log for this gene
    return [mut for mut in mut_list if mut not in skip_set]

def get_mutation_data(ntposnt):
    """Return zero-based position and nucleotides for a mutation string such as G123A."""
    # A gDNA/intronic token ('gd.T5000C') reaches int(ntposnt[1:-1]) as
    # int('d.T5000') and raises. Measured: across both legacy parsers a gd. token
    # yields 0 real parses out of 24 attempts, so this early return is unreachable
    # for any input that works today -- 'd' is not in ACGTU and '.' is not a digit,
    # so no token matching [ACGTU]+\d+[ACGTU]+ can start with the prefix.
    if is_intronic_token(ntposnt):
        return None, None
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    position = int(ntposnt[1:-1]) - 1  # Convert 1-based token position to 0-based index
    return position, (original_nt, mutant_nt)

def get_mutation_data_bioAccurate(ntposnt, is_nt):
    """Return one-based position and nucleotides for a mutation string; skips stop codons.

    is_nt is REQUIRED (no default): pass True for a nucleotide mutation token
    (e.g. G123A) or False for an amino-acid token (e.g. R213W). When True, a token
    whose ref/alt are not nucleotides (ACGTU) returns (None, None) instead of
    silently coercing an off-alphabet character into a position/ref/alt — this
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


# ---------------------------------------------------------------------------
# Length-aware variant representation (non-SNV support).
#
# Everything below is ADDITIVE. get_mutation_data (above) and
# get_mutation_data_bioAccurate (above) are the SNV path and are NOT modified:
# same signatures, same return shapes, same callers. A pipeline opts into
# non-SNV support by using parse_variant/splice_seq instead.
#
# The reason this cannot be a flag on the existing parsers: they hold REF and
# ALT as CHARACTER POSITIONS in an unstructured string (`ntposnt[0]`,
# `ntposnt[-1]`, `int(ntposnt[1:-1])`), so there is no len(REF) to compare
# against len(ALT). A record is required before the question is even askable.
# ---------------------------------------------------------------------------

_NT_ALPHABET = frozenset("ACGTU")

# Strict superset of the legacy SNV grammar. Uniquely decodable because the base
# and digit character classes are disjoint: greedy [ACGTU]+ cannot cross a digit
# and greedy \d+ cannot cross a base. Every legacy token parses identically.
_VARIANT_NT_RE = re.compile(r"^([ACGTUacgtu]+)([0-9]+)([ACGTUacgtu]+)$")
_VARIANT_AA_RE = re.compile(r"^([A-Za-z*]+)([0-9]+)([A-Za-z*]+)$")

# Multi-base sibling of _LOG_MUTATION_RE (see above). That pattern is
# [ACGT][0-9]+[ACGT] -- single base either side -- so an indel named in a
# validation log silently fails to match and is therefore never skipped.
_LOG_VARIANT_RE = re.compile(
    r"^(?P<gene>[^:]+): mutation (?P<mut>[ACGTacgtu]+[0-9]+[ACGTUacgtu]+)\b"
)

_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "u": "a", "n": "n",
}


@dataclass(frozen=True)
class Variant:
    """A length-aware mutation record.

    pos is 1-BASED and points at the first REF base, matching both the ntposnt
    token grammar and VCF. Use `pos0` for slicing. Frozen so it is hashable and
    usable as a dict key.

    orientation declares which strand ref/alt are written on. There is no safe
    default to infer for a minus-strand indel, so it is carried rather than
    guessed; 'genomic' matches the VCF convention.
    """

    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    orientation: str = "genomic"

    def __post_init__(self):
        if self.pos < 1:
            raise ValueError(f"Variant.pos is 1-based, got {self.pos}")
        if not self.ref or not self.alt:
            raise ValueError(
                "Variant.ref and Variant.alt must be non-empty; use the VCF "
                "anchor-base convention for pure insertions/deletions"
            )
        if self.orientation not in ("genomic", "transcript"):
            raise ValueError(
                f"orientation must be 'genomic' or 'transcript', got {self.orientation!r}"
            )

    @property
    def pos0(self) -> int:
        """0-based index of the first REF base, for sequence slicing."""
        return self.pos - 1

    @property
    def length_delta(self) -> int:
        """len(ALT) - len(REF). Zero for SNVs and MNVs."""
        return len(self.alt) - len(self.ref)

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1

    @property
    def kind(self) -> str:
        """One of: snv, mnv, insertion, deletion, delins."""
        if self.is_snv:
            return "snv"
        if len(self.ref) == len(self.alt):
            return "mnv"
        if len(self.ref) == 1:
            return "insertion"
        if len(self.alt) == 1:
            return "deletion"
        return "delins"

    def token(self) -> str:
        """Round-trip serialization. NOT canonical -- use canonical_token()."""
        return f"{self.ref}{self.pos}{self.alt}"


def parse_variant(token, is_nt, gene=None, orientation="genomic"):
    """Parse a mutation token into a Variant, or return None if it does not parse.

    is_nt is REQUIRED and carries the same meaning as in
    get_mutation_data_bioAccurate: True for nucleotide tokens (G123A,
    ACAA112217430A), False for amino-acid tokens (R213W). The nt grammar rejects
    off-alphabet tokens rather than coercing them, so an aa token handed to the
    nt path returns None instead of a silently wrong record.

    Returns None (never raises) for: stop-codon tokens, tokens that do not match
    the grammar, and nt tokens carrying non-ACGTU characters. This mirrors the
    (None, None) convention of get_mutation_data_bioAccurate so callers that
    already branch on a falsy parse need no new error handling.

    When you already hold ref/alt as separate fields -- a 4-column mutations
    file, a VCF record -- construct Variant(...) directly instead of formatting
    a token just to re-parse it.
    """
    if not isinstance(token, str):
        return None
    token = token.strip()
    if not token:
        return None

    # Same stop-codon guard as get_mutation_data_bioAccurate.
    if "Stop" in token or "Sto" in token:
        return None

    pattern = _VARIANT_NT_RE if is_nt else _VARIANT_AA_RE
    m = pattern.match(token)
    if not m:
        return None

    ref, pos_str, alt = m.group(1), m.group(2), m.group(3)
    if is_nt:
        ref, alt = ref.upper(), alt.upper()
        if not (set(ref) <= _NT_ALPHABET and set(alt) <= _NT_ALPHABET):
            return None

    try:
        pos = int(pos_str)
    except ValueError:
        return None
    if pos < 1:
        return None

    return Variant(pos=pos, ref=ref, alt=alt, gene=gene, orientation=orientation)


def get_variant_data(ntposnt, is_nt=True):
    """Indel-capable drop-in for get_mutation_data.

    Returns (position0, (ref, alt)) -- the SAME shape get_mutation_data returns,
    with the same 0-BASED position convention -- but ref and alt may be multiple
    bases. Returns (None, None) instead of raising when the token does not parse,
    so a caller that already guards on a falsy position needs no new handling.

    This exists so a pipeline can migrate one call site at a time:

        pos, (ref, alt) = get_mutation_data(tok)      # SNV only, raises on indel
        pos, (ref, alt) = get_variant_data(tok)       # indel-capable, never raises

    The caller must then use splice_seq(seq, pos, ref, alt) rather than
    update_str(seq, alt, pos), because only the former honours len(ref).
    Anything downstream that assumes len(ref) == len(alt) == 1 -- a codon
    extraction, a per-position join, an amino-acid projection -- is NOT made
    correct by this function alone; it only makes the variant representable.
    """
    v = parse_variant(ntposnt, is_nt=is_nt)
    if v is None:
        return None, None
    return v.pos0, (v.ref, v.alt)


def get_variant_data_bioAccurate(ntposnt, is_nt=True):
    """1-based sibling of get_variant_data, matching get_mutation_data_bioAccurate.

    Returns (position1, (ref, alt)) or (None, None). Same contract as
    get_mutation_data_bioAccurate, including the stop-codon and off-alphabet
    guards, but with multi-base ref/alt.
    """
    v = parse_variant(ntposnt, is_nt=is_nt)
    if v is None:
        return None, None
    return v.pos, (v.ref, v.alt)


def _trim_alleles(pos, ref, alt):
    """Parsimonious VCF trim: drop the shared suffix, then the shared prefix.

    Both alleles keep at least one base (the VCF anchor convention), so an SNV is
    untouched and ACAA/A at p stays ACAA/A at p rather than collapsing.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _left_align(pos, ref, alt, seq, seq_offset=1):
    """Shift an indel as far 5' as the reference allows (bcftools norm semantics).

    seq is the reference sequence ref/alt are written against; seq_offset is the
    1-based coordinate of seq[0]. Requires the reference because left-alignment
    is a property of the surrounding sequence, not of the alleles alone.
    """
    idx = pos - seq_offset
    if idx < 0 or idx + len(ref) > len(seq):
        # Cannot see the flank; leave the representation alone rather than
        # shifting against a sequence that does not cover it.
        return pos, ref, alt
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            if idx <= 0:
                break
            prev = seq[idx - 1]
            ref, alt, idx = prev + ref, prev + alt, idx - 1
        else:
            break
    # Re-anchor: a fully trimmed indel has an empty allele, which VCF forbids.
    if not ref or not alt:
        if idx <= 0:
            return idx + seq_offset, ref or seq[idx], alt or seq[idx]
        prev = seq[idx - 1]
        ref, alt, idx = prev + ref, prev + alt, idx - 1
    return idx + seq_offset, ref, alt


def canonical_token(variant, seq=None, seq_offset=1):
    """Return the canonical string form of a variant.

    This exists because BFF's entire cross-pipeline join is exact string equality
    on a concatenated {ref}{pos}{alt} token -- see the three mapping identities in
    exon_aware_mapping.py and the pkey mint in spliceai/bin/spliceai-parser.py.
    Two textual spellings of one deletion therefore become two primary keys and
    the miss is silent (the lookup returns None and the row is dropped).

    Without `seq` this normalizes representation only: uppercase plus
    parsimonious trimming. That is deterministic and sufficient when every
    producer already agrees on placement.

    With `seq` it additionally left-aligns, which is the only thing that makes
    two spellings inside a tandem repeat collapse to one key. In a CAACAA repeat
    the VCF-left and HGVS-3' placements differ by up to 3 bp -- the same
    magnitude as a 3 bp deletion -- so for repeat-adjacent indels passing `seq`
    is not optional if the token is going to be used as a join key.
    """
    pos, ref, alt = variant.pos, variant.ref.upper(), variant.alt.upper()
    if seq is not None:
        pos, ref, alt = _left_align(pos, ref, alt, seq.upper(), seq_offset)
    pos, ref, alt = _trim_alleles(pos, ref, alt)
    return f"{ref}{pos}{alt}"


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


INTRONIC_PREFIX = "gd."
CHROM_PREFIX = "ch."
# Every prefix that means "NOT ORF-relative". Kept as one tuple so the gates in
# codon_usage / rare_codon / evmutation / build_mutant_sequences_for_gene and the
# router in exon_aware_mapping cannot disagree about which spaces exist.
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


def split_intronic_tokens(tokens):
    """Partition tokens into (orf_space, gdna_space), preserving order."""
    orf_tokens, gd_tokens = [], []
    for tok in tokens:
        (gd_tokens if is_intronic_token(tok) else orf_tokens).append(tok)
    return orf_tokens, gd_tokens


def warn_intronic_unsupported(pipeline, gene, tokens, why, stream=None):
    """Emit the hard warning for intronic tokens a pipeline cannot score.

    Returns the number warned so the caller can fold it into its own accounting.

    This is deliberately loud and deliberately specific about WHY. A frame- or
    protein-dependent pipeline has no defensible output for an intronic variant:
    an intron has no reading frame and no residue, so every codon or amino-acid
    column would be undefined. Emitting a row anyway -- especially one whose
    delta columns come out 0.0 because WT and MUT are identical -- would read as
    a measured 'no effect' rather than 'not modelled', which is the failure this
    whole layer exists to prevent.
    """
    if not tokens:
        return 0
    out = stream if stream is not None else sys.stderr
    label = f"{gene}: " if gene else ""
    print(f"[{pipeline}] *** WARNING *** {label}{len(tokens)} intronic / non-ORF "
          f"(gd., ch.) token(s) were EXCLUDED.", file=out)
    print(f"[{pipeline}] These variants lie OUTSIDE the open reading frame. This "
          f"pipeline WILL NOT produce biologically meaningful results for them, "
          f"and forcing one through would yield a number that LOOKS like a "
          f"measurement and is not.", file=out)
    print(f"[{pipeline}] Reason: {why}", file=out)
    for tok in tokens:
        print(f"[{pipeline}]   excluded (NOT scored): {tok}", file=out)
    return len(tokens)


def splice_seq(seq, pos0, ref, alt, validate=True):
    """Apply a length-changing edit: seq[:pos0] + alt + seq[pos0 + len(ref):].

    The length-aware counterpart of update_str, which hardcodes a stride of
    exactly one base (`s[:pos] + c + s[pos + 1:]`) and so cannot express an indel.
    update_str is left in place and unchanged for the SNV path.

    pos0 is 0-BASED (use Variant.pos0). Raises ValueError on an out-of-range
    position or, unless validate=False, on a REF that does not match the
    sequence. The REF assertion is on by default deliberately: an unverified
    splice at a wrong coordinate produces a plausible sequence and no error,
    which is the failure mode this whole record exists to prevent.
    """
    if pos0 < 0:
        raise ValueError(f"pos0 must be >= 0, got {pos0}")
    end = pos0 + len(ref)
    if end > len(seq):
        raise ValueError(
            f"variant REF spans {pos0}..{end} (0-based) but sequence is only "
            f"{len(seq)} long"
        )
    if validate:
        observed = seq[pos0:end].upper().replace("U", "T")
        expected = ref.upper().replace("U", "T")
        if observed != expected:
            raise ValueError(
                f"REF mismatch at 0-based position {pos0}: sequence has "
                f"{seq[pos0:end]!r}, variant declares {ref!r}"
            )
    return seq[:pos0] + alt + seq[end:]


def align_wt_to_mut(wt_len, edit_offset, ref_len, alt_len):
    """WT-frame alignment across a length-changing edit.

    Returns a list of length wt_len: entry i is the MUT index corresponding to WT
    index i, or None where the edit deleted that base so no counterpart exists.

    This is the one piece of machinery that makes a per-position delta definable
    under an indel, and it is why BFF can stop emitting a positional zip that
    silently compares base i of one allele to a different base of the other.

    Three regions, all in WT coordinates:
      before the edit   identity -- unaffected by anything downstream
      inside the REF    the first min(ref_len, alt_len) bases pair up; any REF
                        base beyond that was deleted and maps to None
      after the edit    shifted by (alt_len - ref_len)

    Bases the ALT *inserted* have no WT index and therefore do not appear: the
    output is in the WT frame, matching the convention SpliceAI uses (zero-pad a
    deletion, collapse an insertion) -- the only correct alignment convention
    already present anywhere in this codebase.
    """
    delta = alt_len - ref_len
    out = []
    for i in range(wt_len):
        if i < edit_offset:
            out.append(i)
        elif i < edit_offset + ref_len:
            k = i - edit_offset
            out.append(edit_offset + k if k < alt_len else None)
        else:
            out.append(i + delta)
    return out


def infer_edit_span(wt_seq, mut_seq, frameshift=False):
    """Recover (edit_offset, ref_len, alt_len) from two sequences.

    For pipelines that hold a WT and a MUT sequence but not the variant record --
    the DTU protein tools score proteins and never see the nucleotide token. The
    edit is located by trimming the common prefix and the common suffix, which is
    exact for a single contiguous edit and needs no aligner.

    frameshift=True forces alt_len=0, meaning NOTHING after the edit aligns.
    That is required and cannot be inferred here: prefix/suffix trimming returns
    the minimal edit, so a frameshift like MKEWLTCD -> MKNG is reported as a
    6->2 replacement, which would then pair E with N and W with G. Those residues
    are not counterparts; after a frameshift every downstream residue is a
    different one. The caller knows the consequence class -- see
    protein_consequence -- and must say so.
    """
    if wt_seq == mut_seq:
        return 0, 0, 0
    n = min(len(wt_seq), len(mut_seq))
    pre = 0
    while pre < n and wt_seq[pre] == mut_seq[pre]:
        pre += 1
    suf = 0
    while suf < (n - pre) and wt_seq[len(wt_seq) - 1 - suf] == mut_seq[len(mut_seq) - 1 - suf]:
        suf += 1
    ref_len = len(wt_seq) - pre - suf
    alt_len = len(mut_seq) - pre - suf
    if frameshift:
        # Everything from the edit to the end of the WT loses its counterpart.
        return pre, len(wt_seq) - pre, 0
    return pre, ref_len, alt_len


def aligned_pairs(wt_values, mut_values, edit_offset, ref_len, alt_len):
    """Yield (wt_index, wt_value, mut_value) for WT positions with a counterpart.

    Positions deleted by the edit, or projecting outside the mutant vector, are
    skipped rather than coalesced to zero. That distinction matters: a coalesced
    0.0 reads downstream as "measured, no change", which is a fabricated
    observation, whereas a skipped position is simply absent.
    """
    mapping = align_wt_to_mut(len(wt_values), edit_offset, ref_len, alt_len)
    for i, j in enumerate(mapping):
        if j is None or not (0 <= j < len(mut_values)):
            continue
        yield i, wt_values[i], mut_values[j]


def _translate_codons(nt):
    """Translate a nucleotide string, ONE CHARACTER PER RESIDUE, stop kept as '*'.

    Deliberately NOT translate_orf_sequence: that one ends with .rstrip('*'),
    which discards exactly the stop position protein_consequence reports as
    new_stop_aa_pos. It is the right function for "give me the protein"; it is
    the wrong one for "where does translation now terminate".

    The trailing partial codon is trimmed before translating rather than after,
    because Biopython emits a BiopythonWarning for a length that is not a
    multiple of three.

    Note the repo carries two stop conventions: codon_to_aa maps a stop to the
    4-character string 'Stop' (see get_mutant_aa), while Biopython uses '*'.
    Per-residue indexing and len() require the single character, so this uses
    Biopython.
    """
    cleaned = nt.strip().upper().replace("U", "T")
    cleaned = cleaned[:len(cleaned) - len(cleaned) % 3]
    if not cleaned:
        return ""
    return str(Seq(cleaned).translate(to_stop=False))


def _trim_aa_span(aa_pos, wt_aa, mut_aa):
    """Minimal aa representation: drop residues identical on both sides.

    HGVS convention -- a 3 bp deletion whose codon span happens to start with an
    unchanged residue should be reported at the residue that actually changed.
    Unlike the nucleotide trim this may empty an allele, which is correct: a pure
    in-frame deletion has no mutant residue.
    """
    while wt_aa and mut_aa and wt_aa[-1] == mut_aa[-1]:
        wt_aa, mut_aa = wt_aa[:-1], mut_aa[:-1]
    while wt_aa and mut_aa and wt_aa[0] == mut_aa[0]:
        wt_aa, mut_aa, aa_pos = wt_aa[1:], mut_aa[1:], aa_pos + 1
    return aa_pos, wt_aa, mut_aa


def protein_consequence(variant, orf_seq, mut_orf_seq=None):
    """Codon-aware protein consequence of a nucleotide variant.

    An amino-acid token cannot express an indel's protein effect -- you need the
    ORF and the reading frame. This walks the edit through translation and
    returns the columns the DTU pipelines need:

        aa_pos          1-based residue where the change starts
        wt_aa, mut_aa   the replaced span as MULTI-CHARACTER STRINGS, minimally
                        trimmed. 'M'/'V' for an SNV -- byte-identical to the
                        current single-residue output -- 'MK'/'M' for a
                        2-residue in-frame deletion, '' for a frameshift.
                        Strings, not lists: a list written through write_tsv
                        becomes the repr "['M', 'K']" and needs literal_eval on
                        read, which then fails on the bare 'M' of an SNV row.
        n_aa_wt/n_aa_mut  len() of each, so "how many residues" is answerable
                        without parsing anything
        aa_consequence  snv | mnv | inframe_del | inframe_ins | inframe_delins |
                        frameshift | stop_gained | stop_lost | synonymous
        new_stop_aa_pos 1-based residue of the mutant's first stop, or None.
                        This is the only informative number for a frameshift,
                        where wt_aa/mut_aa are deliberately empty because
                        wt_aa[i] and mut_aa[i] no longer describe the same site.

    Returns None if the variant does not lie within the ORF.
    """
    if variant.pos0 < 0 or variant.pos0 + len(variant.ref) > len(orf_seq):
        return None
    if mut_orf_seq is None:
        mut_orf_seq = splice_seq(orf_seq, variant.pos0, variant.ref, variant.alt,
                                 validate=False)

    wt_prot_full = _translate_codons(orf_seq)
    mut_prot_full = _translate_codons(mut_orf_seq)
    new_stop = mut_prot_full.find('*')
    new_stop_aa_pos = new_stop + 1 if new_stop >= 0 else None

    delta = variant.length_delta
    first_codon = variant.pos0 // 3
    aa_pos = first_codon + 1

    if delta % 3 != 0:
        # Every residue from here to the new stop differs, and mut_aa[i] does not
        # correspond to wt_aa[i]. Emitting the two peptides as a "pair" would
        # invite exactly that false comparison, so both are left empty.
        return {
            'aa_pos': aa_pos,
            'wt_aa': '',
            'mut_aa': '',
            'n_aa_wt': 0,
            'n_aa_mut': 0,
            'aa_consequence': 'frameshift',
            'new_stop_aa_pos': new_stop_aa_pos,
        }

    last_codon = (variant.pos0 + len(variant.ref) - 1) // 3
    wt_span = orf_seq[first_codon * 3:(last_codon + 1) * 3]
    mut_span = mut_orf_seq[first_codon * 3:(last_codon + 1) * 3 + delta]
    wt_aa = _translate_codons(wt_span)
    mut_aa = _translate_codons(mut_span)
    aa_pos, wt_aa, mut_aa = _trim_aa_span(aa_pos, wt_aa, mut_aa)

    if '*' in mut_aa and '*' not in wt_aa:
        consequence = 'stop_gained'
    elif '*' in wt_aa and '*' not in mut_aa:
        consequence = 'stop_lost'
    elif not wt_aa and not mut_aa:
        consequence = 'synonymous'
    elif delta < 0:
        consequence = 'inframe_del' if not mut_aa else 'inframe_delins'
    elif delta > 0:
        consequence = 'inframe_ins' if not wt_aa else 'inframe_delins'
    elif len(wt_aa) == 1 and len(mut_aa) == 1:
        consequence = 'snv'
    else:
        consequence = 'mnv'

    return {
        'aa_pos': aa_pos,
        'wt_aa': wt_aa,
        'mut_aa': mut_aa,
        'n_aa_wt': len(wt_aa),
        'n_aa_mut': len(mut_aa),
        'aa_consequence': consequence,
        'new_stop_aa_pos': new_stop_aa_pos,
    }


def apply_variants(seq, variants, validate=True):
    """Apply SEVERAL edits to one sequence as a single haplotype.

    The Enh13 target case is a strain carrying two edits scored as one phenotype
    (`Enh13^SOX9*-3bp SRY*-12bp`). BFF's ingestion is one-token-per-line and its
    pkey is `{GENE}-{single token}`, so "apply both, score once" has no
    representation. This is the primitive for it.

    Edits are applied in DESCENDING position order. That is not a preference --
    it is required. Applying a 5' edit first shifts every downstream coordinate
    by its length delta, so a 3'-side variant parsed against the original
    sequence would then splice at the wrong place. Going 3' -> 5' leaves every
    not-yet-applied coordinate untouched.

    Overlapping edits raise. Two REF spans that touch the same base have no
    well-defined joint meaning, and silently applying one then the other would
    make the result depend on ordering.
    """
    ordered = sorted(variants, key=lambda v: v.pos, reverse=True)
    for later, earlier in zip(ordered, ordered[1:]):
        if earlier.pos0 + len(earlier.ref) > later.pos0:
            raise ValueError(
                f"overlapping edits cannot be applied as one haplotype: "
                f"{earlier.token()} spans {earlier.pos}..{earlier.pos + len(earlier.ref) - 1} "
                f"and {later.token()} starts at {later.pos}"
            )
    out = seq
    for v in ordered:
        out = splice_seq(out, v.pos0, v.ref, v.alt, validate=validate)
    return out


def canonical_haplotype_token(variants, seq=None, seq_offset=1):
    """One stable key for a set of edits applied together.

    Components are canonicalized individually then joined 5'->3' with '+', so
    the same haplotype always mints the same key regardless of the order the
    edits were listed in the input file. Pass `seq` for the same reason as in
    canonical_token: without it, repeat-adjacent indels are not left-aligned and
    two spellings of one haplotype will mint two keys.
    """
    parts = sorted(
        ((v.pos, canonical_token(v, seq=seq, seq_offset=seq_offset)) for v in variants),
        key=lambda pair: pair[0],
    )
    return "+".join(tok for _, tok in parts)


def revcomp_seq(seq):
    """Reverse-complement, length-aware and fail-loud on unknown bases.

    rc_base in exon_aware_mapping.py is a single-character dict lookup with a
    `.get(b, b)` fallback: handed a multi-base string it returns that string
    UNCOMPLEMENTED and UNREVERSED, silently. Harmless while every caller passes
    one character; wrong the moment an indel reaches a minus-strand element.
    """
    for base in seq:
        if base not in _COMPLEMENT:
            raise ValueError(f"cannot complement base {base!r} in {seq!r}")
    return str(Seq(seq).reverse_complement())


def rc_variant(variant, genomic_pos=None):
    """Flip a variant onto the opposite strand.

    Complements and REVERSES both alleles, and flips the orientation label.

    genomic_pos, when supplied, is the genomic coordinate of the variant's
    FIRST REF BASE IN TRANSCRIPT ORDER -- i.e. exactly what
    `tx_to_genome[tx_pos - 1]` yields in exon_aware_mapping.py. On the minus
    strand transcript order runs 3'->5' along the genome, so that base is the
    RIGHTMOST of the REF span and the genomic (leftmost) start is
    genomic_pos - (len(ref) - 1). Getting this wrong is invisible for an SNV,
    where the span is one base and the correction is zero.
    """
    rc_ref = revcomp_seq(variant.ref)
    rc_alt = revcomp_seq(variant.alt)
    flipped = "transcript" if variant.orientation == "genomic" else "genomic"
    if genomic_pos is None:
        return replace(variant, ref=rc_ref, alt=rc_alt, orientation=flipped)
    return replace(
        variant,
        pos=genomic_pos - (len(variant.ref) - 1),
        ref=rc_ref,
        alt=rc_alt,
        orientation=flipped,
    )


def get_mutant_aa(ntmut, ntseq, aaseq=None, index=0):
    pos_0_indexed = ntmut[0] - 1 - index

    # Check if the calculated position is valid for the sequence
    if not (0 <= pos_0_indexed < len(ntseq)):
        return None

    codon_start_pos = (pos_0_indexed // 3) * 3

    # Extract the original codon and translate it
    original_codon = ntseq[codon_start_pos: codon_start_pos + 3]
    wtaa = codon_to_aa.get(original_codon, 'X')

    mutseq = update_str(ntseq, ntmut[1][1], pos_0_indexed)

    mutated_codon = mutseq[codon_start_pos: codon_start_pos + 3]
    mutaa = codon_to_aa.get(mutated_codon, 'X')

    if aaseq is not None:
        aa_pos_0_indexed = codon_start_pos // 3

    aa_position_1_based = (codon_start_pos // 3) + 1

    return (aa_position_1_based, (wtaa, mutaa)), original_codon


def _detect_annotation_format(annotation_file):
    """Infer whether an annotation file resembles GTF, GFF3, or a custom tab format."""
    suffix = Path(annotation_file).suffix.lower()
    if suffix == ".gtf":
        return "gtf"
    if suffix in {".gff", ".gff3"}:
        return "gff3"

    try:
        with open(annotation_file, "r") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 9:
                    attr = fields[8]
                    if '"' in attr:
                        return "gtf"
                    if "=" in attr or attr.startswith("ID=") or attr.startswith("Parent="):
                        return "gff3"
                    return "gff3"
                if len(fields) >= 7:
                    return "custom"
    except FileNotFoundError:
        return "custom"
    return "custom"


def _parse_attributes(attr_field, fmt):
    """Parse the attribute column of a GTF/GFF record into a key/value dict."""
    attrs = {}
    text = attr_field.strip().strip(";")
    if not text:
        return attrs

    if fmt == "gff3":
        for part in text.split(";"):
            part = part.strip()
            if not part or "=" not in part:
                continue
            key, value = part.split("=", 1)
            key = key.strip()
            value = unquote(value.strip())
            attrs[key] = value
    else:
        for part in text.split(";"):
            part = part.strip()
            if not part:
                continue
            if " " not in part:
                attrs[part] = ""
                continue
            key, value = part.split(" ", 1)
            key = key.strip()
            value = value.strip().strip('"')
            attrs[key] = value
    return attrs


def _normalize_chrom_name(raw_name, assembly):
    """Collapse assorted chromosome labels into bare chromosome numbers for a given assembly."""
    if not raw_name:
        return raw_name

    chrom = raw_name.strip()
    rev_map = {}
    if assembly in chromosome_map:
        mapping = chromosome_map[assembly]
        rev_map = {v: k for k, v in mapping.items()}
        rev_map.update({v.split('.')[0]: k for k, v in mapping.items()})
        rev_map.update({f"chr{k}": k for k in mapping.keys()})

    if chrom in rev_map:
        return rev_map[chrom]
    chrom_core = chrom.split('.')[0]
    if chrom_core in rev_map:
        return rev_map[chrom_core]
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
        if chrom in rev_map:
            return rev_map[chrom]
        return chrom
    return chrom


def _split_multi_value(value):
    """Split a semi-colon or comma separated field into a list of trimmed strings."""
    if not value:
        return []
    parts = []
    for part in value.replace(';', ',').split(','):
        candidate = part.strip()
        if candidate:
            parts.append(candidate)
    return parts


def _infer_transcript_priority(attrs):
    """Score a transcript based on annotation tags to favour well-supported isoforms."""
    priority = 0
    biotype = (
        attrs.get("transcript_biotype")
        or attrs.get("biotype")
        or attrs.get("gene_biotype")
        or attrs.get("transcript_type")
        or attrs.get("gbkey")
        or ""
    ).lower()
    if "protein" in biotype or "coding" in biotype or "cds" in biotype:
        priority = max(priority, 2)
    elif "mrna" in biotype or "messenger" in biotype:
        priority = max(priority, 1)

    product = attrs.get("product", "").lower()
    if "protein" in product:
        priority = max(priority, 1)

    tag_field = attrs.get("tag") or ""
    if tag_field:
        tag_values = {piece.strip().lower() for piece in tag_field.replace(';', ',').split(',') if piece.strip()}
        if "mane select" in tag_values:
            priority = max(priority, 4)
        if "refseq select" in tag_values:
            priority = max(priority, 3)
        if "ccds" in tag_values:
            priority = max(priority, 2)

    return priority


def _prepare_custom_annotation(genename, annotation_file):
    """Load gene coordinates from a simple tab-delimited annotation generated by this project."""
    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue

            gene_symbol = fields[0]
            chrom = fields[1].replace('chr', '')
            strand = fields[2]
            tx_start = int(fields[3])
            tx_end = int(fields[4])

            exon_starts = [int(x) for x in fields[5].strip(',').split(',') if x.strip()]
            exon_ends = [int(x) for x in fields[6].strip(',').split(',') if x.strip()]
            exons = list(zip(exon_starts, exon_ends))

            if gene_symbol.upper() == genename.upper():
                return {
                    "chrom": chrom,
                    "strand": strand,
                    "tx_start": tx_start,
                    "tx_end": tx_end,
                    "cds_start": tx_start,
                    "cds_end": tx_end,
                    "exons": exons,
                    "transcript_id": f"{gene_symbol}_transcript",
                }
    return None


def _prepare_structured_annotation(genename, annotation_file, assembly, fmt, transcript_id=None):
    """Extract gene and transcript features from RefSeq/Ensembl-style GTF or GFF3 files.

    Args:
        genename: Gene symbol to look up.
        annotation_file: Path to annotation file.
        assembly: Reference assembly.
        fmt: Detected annotation format ('gtf' or 'gff3').
        transcript_id: Optional specific transcript ID to force selection of.
    """
    gene_upper = genename.upper()
    target_gene_ids = set()
    target_transcripts = {}
    transcript_features = {
        "transcript",
        "mrna",
        "ncrna",
        "lnc_rna",
        "primary_transcript",
        "pre_mrna",
        "mirna",
        "trna",
        "rrna",
        "snrna",
        "snorna",
        "scrna",
        "sca_rna",
    }
    gene_section_active = False

    try:
        with open(annotation_file, "r") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue

                seqname, _source, feature, start, end, _score, strand, _frame, attrs_raw = cols
                feature_lc = feature.lower()
                attrs = _parse_attributes(attrs_raw, fmt)

                attr_gene_names = {
                    attrs.get("gene_name"),
                    attrs.get("gene"),
                    attrs.get("gene_symbol"),
                    attrs.get("Name"),
                }
                attr_gene_names = {name.upper() for name in attr_gene_names if isinstance(name, str)}

                attr_gene_ids = set()
                attr_gene_ids.update(_split_multi_value(attrs.get("gene_id")))
                if feature_lc == "gene":
                    attr_gene_ids.update(_split_multi_value(attrs.get("ID")))
                parent_values = _split_multi_value(attrs.get("Parent"))
                attr_gene_ids.update(parent_values)

                matches_symbol = gene_upper in attr_gene_names
                matches_gene_id = bool(target_gene_ids.intersection(attr_gene_ids))
                matches = matches_symbol or matches_gene_id

                if feature_lc == "gene" and matches:
                    target_gene_ids.update(attr_gene_ids)
                    target_gene_ids.update(_split_multi_value(attrs.get("ID")))
                    target_gene_ids.update(_split_multi_value(attrs.get("gene_id")))
                    target_gene_ids.add(genename)
                    gene_section_active = True
                elif feature_lc == "gene" and gene_section_active:
                    break

                current_tid = None

                if feature_lc in transcript_features:
                    current_tid = attrs.get("transcript_id") or attrs.get("ID")
                elif feature_lc == "cds" and not attrs.get("transcript_id") and attrs.get("Parent"):
                    possible = [val for val in parent_values if not val.startswith("gene-")]
                    if possible:
                        current_tid = possible[0]

                if feature_lc == "gene" and matches_symbol and not target_gene_ids:
                    target_gene_ids.update(attr_gene_ids)

                if feature_lc in transcript_features and (matches or bool(set(parent_values).intersection(target_gene_ids)) or (attrs.get("gene_id") and attrs.get("gene_id") in target_gene_ids)):
                    if not current_tid:
                        continue
                    rec = target_transcripts.setdefault(current_tid, {
                        "chrom": None,
                        "strand": strand if strand in "+-" else None,
                        "exons": [],
                        "attrs": {},
                        "matched": False,
                    })
                    rec["chrom"] = _normalize_chrom_name(seqname, assembly)
                    if strand in "+-":
                        rec["strand"] = strand
                    rec["attrs"].update(attrs)
                    rec["matched"] = True
                    target_gene_ids.update(attr_gene_ids)
                    continue

                if feature_lc == "exon":
                    exon_transcript_ids = _split_multi_value(attrs.get("transcript_id"))
                    if not exon_transcript_ids:
                        exon_transcript_ids = [val for val in parent_values if val]

                    start_i, end_i = int(start), int(end)

                    for tid in exon_transcript_ids:
                        rec = target_transcripts.setdefault(tid, {
                            "chrom": None,
                            "strand": None,
                            "exons": [],
                            "attrs": {},
                            "matched": False,
                        })

                        if matches or tid in target_transcripts and target_transcripts[tid]["matched"] or bool(set(parent_values).intersection(target_gene_ids)):
                            rec["matched"] = rec["matched"] or matches or tid in target_transcripts and target_transcripts[tid]["matched"]
                            if attrs.get("gene_id"):
                                target_gene_ids.add(attrs["gene_id"])
                        elif matches_symbol:
                            rec["matched"] = True

                        if not rec["matched"]:
                            continue

                        rec["chrom"] = rec["chrom"] or _normalize_chrom_name(seqname, assembly)
                        if strand in "+-":
                            rec["strand"] = strand
                        rec["exons"].append((start_i, end_i))
                        rec["attrs"].update(attrs)

                if matches_symbol:
                    target_gene_ids.update(attr_gene_ids)

    except FileNotFoundError:
        return None

    candidates = []
    for tid, rec in target_transcripts.items():
        if not rec["matched"] or not rec["exons"]:
            continue
        if not rec["chrom"] or rec["strand"] not in "+-":
            continue

        exons_sorted = sorted(rec["exons"], key=lambda x: (x[0], x[1]))
        unique_exons = []
        for exon in exons_sorted:
            if not unique_exons or unique_exons[-1] != exon:
                unique_exons.append(exon)

        tx_start = min(s for s, _ in unique_exons)
        tx_end = max(e for _, e in unique_exons)
        exon_len = sum(e - s + 1 for s, e in unique_exons)

        priority = _infer_transcript_priority(rec["attrs"])
        score = (priority, exon_len, len(unique_exons), -tx_start)

        candidates.append((score, tid, {
            "chrom": rec["chrom"],
            "strand": rec["strand"],
            "tx_start": tx_start,
            "tx_end": tx_end,
            "exons": unique_exons,
            "attrs": rec["attrs"],
        }))

    if not candidates:
        return None

    # If a specific transcript_id is requested, try to find it among candidates
    if transcript_id:
        # Try exact match first
        exact = [c for c in candidates if c[1] == transcript_id]
        match = exact
        if not match:
            # Try matching without version suffix (e.g., NM_022162 matches NM_022162.3)
            tid_base = transcript_id.rsplit('.', 1)[0]
            match = [c for c in candidates if c[1].rsplit('.', 1)[0] == tid_base]
        if match:
            _, best_tid, best = match[0]
            # Announce a version substitution. Falling back on the bare accession
            # is the intended behaviour when the caller supplies no version
            # ('NM_022162' -> 'NM_022162.3'), but when a version WAS named and a
            # different one is returned, the caller asked a specific question and
            # got another answer. Transcript versions can differ in exon
            # structure, so every coordinate downstream is then computed against
            # a record the caller did not choose. Measured: --force-cds
            # NM_022876.3 against GRCh38.p14, which carries only NM_022876.2.
            if not exact and "." in transcript_id and best_tid != transcript_id:
                print(
                    f"WARNING: transcript '{transcript_id}' is not present for gene "
                    f"'{genename}' in this annotation; using '{best_tid}' instead "
                    f"(same accession, different version). Exon structure may differ "
                    f"from the version requested.",
                    file=sys.stderr,
                )
        else:
            available = sorted(set(c[1] for c in candidates))
            raise ValueError(
                f"Transcript '{transcript_id}' not found for gene '{genename}'. "
                f"Available transcripts: {available}"
            )
    else:
        _, best_tid, best = max(candidates, key=lambda x: x[0])

    return {
        "chrom": best["chrom"],
        "strand": best["strand"],
        "tx_start": best["tx_start"],
        "tx_end": best["tx_end"],
        "cds_start": best["tx_start"],
        "cds_end": best["tx_end"],
        "exons": best["exons"],
        "transcript_id": best_tid,
    }


def get_genome_loc(genename, annotation_file, assembly="GRCh38", transcript_id=None):
    """Return gene coordinates and exon structure from supported annotation formats.

    Args:
        genename: Gene symbol to look up.
        annotation_file: Path to annotation file (GTF/GFF3/custom).
        assembly: Reference assembly (GRCh37 or GRCh38).
        transcript_id: Optional specific transcript ID to force selection of
            (e.g., NM_022162.3). If provided and found, this transcript is
            selected instead of auto-selection based on priority scoring.

    Returns:
        dict with chrom, strand, tx_start, tx_end, exons, transcript_id, etc.
    """
    fmt = _detect_annotation_format(annotation_file)

    if fmt == "custom":
        # The custom tab-delimited format carries ONE record per gene and no
        # transcript IDs, so a requested transcript cannot be honoured here.
        # Saying so matters: the caller printed "Forcing transcript X for all
        # genes" before reaching this point, and without this line a fabricated
        # accession produced output byte-identical to an unforced run with no
        # signal anywhere that the request had been dropped.
        if transcript_id:
            print(
                f"WARNING: requested transcript '{transcript_id}' is ignored for gene "
                f"'{genename}': '{annotation_file}' is the custom tab-delimited format, "
                f"which holds a single record per gene and no transcript IDs. Supply a "
                f"GTF/GFF3 annotation to select a specific transcript.",
                file=sys.stderr,
            )
        try:
            return _prepare_custom_annotation(genename, annotation_file)
        except Exception as e:
            print(f"Error parsing annotation file {annotation_file}: {e}", file=sys.stderr)
            return None

    try:
        return _prepare_structured_annotation(genename, annotation_file, assembly, fmt, transcript_id=transcript_id)
    except Exception as e:
        print(f"Error parsing annotation file {annotation_file}: {e}", file=sys.stderr)
        return None


def load_mapping(mapping_file: str, mapType: str ='transcript') -> Dict[str, str]:
    """Load a two-column mapping CSV (mutant->mapping) using the specified column name."""

    mapping = {}
    res = validate_mapping_content(mapping_file)
    if not res:  # False on a single-column / unreadable file → graceful empty, no unpack crash (F37)
        return mapping
    tval, delim = res
    try:
        if tval:
            with open(mapping_file, 'r') as f:
                reader = csv.DictReader(f, delimiter= delim)
                for row in reader:
                    if 'mutant' in row and mapType in row:
                        ref_mutation = row['mutant']
                        mapped_mutation = row[mapType]
                        mapping[ref_mutation] = mapped_mutation

    except FileNotFoundError:
        print(f"Warning: Transcript mapping file not found: {mapping_file}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Error loading transcript mapping file {mapping_file}: {e}", file=sys.stderr)

    return mapping

def split_fasta_into_batches(fasta_file, batch_size=100, temp_dir=None):
    """Split a FASTA file into smaller batches for processing

    Args:
        fasta_file: Path to input FASTA file
        batch_size: Number of sequences per batch (default: 100)
        temp_dir: Directory to create batch files in (default: system temp directory)

    Returns:
        List of batch file paths created
    """
    try:
        # Read all sequences from FASTA file
        sequences = read_fasta(fasta_file)
        total_sequences = len(sequences)

        if total_sequences == 0:
            return []

        print(f"Splitting {total_sequences} sequences into batches of {batch_size}")

        # Calculate number of batches needed
        num_batches = math.ceil(total_sequences / batch_size)

        batch_files = []
        sequence_items = list(sequences.items())

        for i in range(num_batches):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, total_sequences)

            # Create batch sequences
            batch_sequences = dict(sequence_items[start_idx:end_idx])
            batch_count = len(batch_sequences)

            # Create temporary batch file (use temp_dir if provided, otherwise system temp)

            if temp_dir:
                # Use provided temp directory
                base_name = os.path.basename(fasta_file).replace('.fasta', f'_batch{i+1}.fasta')
                batch_filename = os.path.join(temp_dir, base_name)
            else:
                # Use system temp directory
                base_name = os.path.basename(fasta_file).replace('.fasta', f'_batch{i+1}.fasta')
                batch_filename = os.path.join(tempfile.gettempdir(), base_name)

            with open(batch_filename, 'w') as f:
                for seq_name, sequence in batch_sequences.items():
                    f.write(f">{seq_name}\n")
                    # Write sequence in lines of 80 characters
                    for j in range(0, len(sequence), 80):
                        f.write(sequence[j:j+80] + "\n")

            batch_files.append(batch_filename)
            print(f"Created batch {i+1}/{num_batches}: {batch_filename} ({batch_count} sequences)")

        return batch_files

    except Exception as e:
        print(f"Error splitting FASTA file {fasta_file}: {e}")
        return []

def combine_batch_outputs(batch_output_files, final_output_file, format_type='netnglyc', original_fasta_file=None):
    """
    Combine multiple batch outputs into a single file for backwards compatibility

    Args:
        batch_output_files: List of individual batch output files
        final_output_file: Path for combined output file
        format_type: Output format type ('netnglyc' or 'netphos')

    Returns:
        bool: True if combination successful
    """
    try:
        if not batch_output_files:
            return False

        print(f"Combining {len(batch_output_files)} batch outputs (format: {format_type})...")

        if format_type == 'netnglyc':
            return _combine_glycosylation_outputs(batch_output_files, final_output_file, original_fasta_file)
        elif format_type == 'netphos':
            return _combine_phosphorylation_outputs(batch_output_files, final_output_file)
        else:
            raise ValueError(f"Unsupported format_type: {format_type}")

    except Exception as e:
        print(f"Error combining batch outputs: {e}")
        return False

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

def _combine_glycosylation_outputs(batch_output_files, final_output_file, original_fasta_file=None):
    """Combine glycosylation prediction batch outputs (NetNGlyc format)"""
    try:
        # Count total sequences for header
        total_sequences = 0
        all_sequence_sections = []
        all_prediction_lines = []

        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()


                import os
                seperator = ['-','_']
                # Use original FASTA file if provided, otherwise fallback to counting Name: lines
                if original_fasta_file and os.path.exists(original_fasta_file):
                    try:
                        #fasta_sequences = read_fasta(original_fasta_file)
                        gene_name, total_sequences = ExtractGeneFromFASTA(original_fasta_file, count=True)
                        # Determine if this is wildtype or mutant based on sequence names
                            # For the first batch, use total sequences from original file
                        if i == 0:
                            print(f"Using original FASTA file for sequence count: {total_sequences} sequences from {gene_name}")
                            # For subsequent batches, the total was previously counted

                    except Exception as e:
                        print(f"Warning: Error reading original FASTA file {original_fasta_file}: {e}")
                        # Fallback to counting Name: lines in NetNGlyc output
                        lines = content.split('\n')
                        name_count = sum(1 for line in lines if line.startswith('Name:'))
                        total_sequences += name_count if name_count > 0 else 1

                else:
                    # Fallback: Count sequences directly from NetNGlyc output
                    lines = content.split('\n')
                    name_count = sum(1 for line in lines if line.startswith('Name:'))

                    if name_count > 0:
                        total_sequences += name_count
                    else:
                        # Parse header line for sequence count: ">debug-GENE-aa-netnglyc\t5 amino acids"
                        for line in lines:
                            if 'amino acids' in line:
                                try:
                                    parts = line.split()
                                    for j, part in enumerate(parts):
                                        if part.isdigit() and j+1 < len(parts) and 'amino' in parts[j+1]:
                                            total_sequences += int(part)
                                            break
                                    break
                                except:
                                    total_sequences += 1  # Ultimate fallback
                            else:
                                total_sequences += 1  # Fallback if no header found

                # Collect sequence display sections and prediction lines
                in_sequence_section = False
                in_prediction_section = False

                lines = content.split('\n')
                for line in lines:
                    if line.strip().startswith('>') and 'amino acids' in line:
                        continue  # Skip individual headers

                    if 'SeqName' in line and 'Position' in line:
                        in_prediction_section = True
                        continue

                    if line.startswith('    ') and len(line.strip()) > 10:
                        in_sequence_section = True
                        all_sequence_sections.append(line)
                    elif in_prediction_section and line.strip() and not line.startswith('#'):
                        all_prediction_lines.append(line)

            except Exception as e:
                print(f"Error reading batch file {batch_file}: {e}")
                continue

        # Write combined output
        with open(final_output_file, 'w') as f:
            # Write header with total sequence count
            f.write(f">{os.path.basename(final_output_file).replace('.out', '')}\t{total_sequences} amino acids\n\n")

            # Write prediction header
            f.write("SeqName                 Position  Potential  N-Glyc result  Comment\n")
            f.write("=" * 70 + "\n")

            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")

            # Write sequence sections
            if all_sequence_sections:
                f.write("\n")
                for line in all_sequence_sections:
                    f.write(line + "\n")

        print(f"Combined {len(batch_output_files)} batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")

        return True

    except Exception as e:
        print(f"Error combining glycosylation outputs: {e}")
        return False

def _combine_phosphorylation_outputs(batch_output_files, final_output_file):
    """Combine phosphorylation prediction batch outputs (NetPhos format)"""
    try:
        total_sequences = 0
        all_prediction_lines = []

        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()

                # Extract sequences count from header
                lines = content.split('\n')
                for line in lines:
                    if 'amino acids' in line and line.startswith('>'):
                        try:
                            seq_count = int(line.split()[1])
                            total_sequences += seq_count
                        except:
                            total_sequences += 1  # Fallback
                        break

                # Collect prediction lines
                for line in lines:
                    if line.startswith('# ') and len(line.split()) >= 7:
                        # This is a prediction line
                        all_prediction_lines.append(line)

            except Exception as e:
                print(f"Error reading phosphorylation batch file {batch_file}: {e}")
                continue

        # Write combined phosphorylation output
        with open(final_output_file, 'w') as f:
            # Write header
            gene_name = os.path.basename(final_output_file).replace('.out', '').replace('-netphos', '')
            f.write(f">{gene_name}\t{total_sequences} amino acids\n")
            f.write("#\n")
            f.write("#  prediction results\n")
            f.write("#\n")
            f.write("# Sequence\t\t   # x   Context     Score   Kinase    Answer\n")
            f.write("# " + "-" * 67 + "\n")

            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")

        print(f"Combined {len(batch_output_files)} phosphorylation batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")

        return True

    except Exception as e:
        print(f"Error combining phosphorylation outputs: {e}")
        return False

def extract_mutation_from_sequence_name(seq_name):
    """Extract mutation ID from sequence name (e.g., 'ZFP36-C330T' -> 'C330T')

    Args:
        seq_name: Sequence name in format 'GENE-MUTATION' or just 'GENE'

    Returns:
        tuple: (gene, mutation_id) or (gene, None) if no mutation found
    """
    if '-' in seq_name:
        parts = seq_name.rsplit('-', 1)
        if len(parts) == 2:
            return parts[0], parts[1]  # gene, mutation_id

    return seq_name, None

def load_validation_failures(log_path):
    """
    Public helper to expose validation log failures as a mapping suitable for filtering.

    Returns:
        dict[str, set[str]]: Uppercase gene -> set of mutation ids to skip.
    """
    return _collect_failures_from_logs(log_path) if log_path else {}


def should_skip_mutation(gene, mutation_id, failure_map):
    """Return True if the given gene/mutation should be filtered based on validation logs."""
    if not failure_map or not gene or not mutation_id:
        return False
    return mutation_id.upper() in failure_map.get(gene.upper(), set())


def process_single_mutation_for_sequence(seq_name, predictions, mapping_dict, is_mutant=True, tool_type='netphos', failure_map=None):
    """Process predictions for one sequence against its specific mutation

    Args:
        seq_name: Sequence name (e.g., 'ZFP36-C330T')
        predictions: List of predictions for this sequence
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant sequences (should be True for single-mutation processing)
        tool_type: 'netphos' or 'netnglyc' for tool-specific field handling

    Returns:
        list: Filtered predictions with pkeys for the specific mutation
    """
    if not is_mutant:
        raise ValueError("process_single_mutation_for_sequence should only be used for mutant sequences")

    # Extract mutation ID from sequence name
    gene, mutation_id = extract_mutation_from_sequence_name(seq_name)
    if mutation_id is None:
        return []

    # Look up this mutation in the mapping
    if mutation_id not in mapping_dict:
        return []

    if should_skip_mutation(gene, mutation_id, failure_map):
        return []

    aaposaa = mapping_dict[mutation_id]  # e.g., "Y110F"

    # Parse amino acid position and mutation info
    position_data = get_mutation_data_bioAccurate(aaposaa, is_nt=False)
    if position_data[0] is None:
        return []

    aa_pos = position_data[0]  # e.g., 110
    aa_tuple = position_data[1]  # e.g., ('Y', 'F')
    target_aa = aa_tuple[1]  # F for mutant

    # Filter predictions for this specific mutation
    results = []
    for pred in predictions:
        # Get position from prediction
        if tool_type == 'netphos':
            pred_pos = pred['pos']
            pred_aa = pred['amino_acid']
        elif tool_type == 'netnglyc':
            pred_pos = pred['position']
            pred_aa = pred['sequon'][0] if pred['sequon'] else None
        else:
            raise ValueError(f"Unsupported tool_type: {tool_type}")

        # Check if this prediction matches the mutation position and amino acid
        if pred_pos == aa_pos and pred_aa == target_aa:
            # Create pkey for this match
            pkey = f"{gene}-{mutation_id}"

            # Add pkey and fix Gene field to prediction
            result_pred = pred.copy()
            result_pred['pkey'] = pkey
            # Fix Gene field to just gene name (not gene-mutation)
            result_pred['Gene'] = gene

            # Map field names to match CSV writer expectations
            if tool_type == 'netnglyc':
                # Map NetNGlyc field names to CSV format
                if 'position' in result_pred:
                    result_pred['pos'] = result_pred.pop('position')
                if 'sequon' in result_pred:
                    result_pred['Sequon'] = result_pred.pop('sequon')
                # Remove seq_name as it's not needed in CSV
                if 'seq_name' in result_pred:
                    result_pred.pop('seq_name')

            results.append(result_pred)

    return results

def parse_predictions_with_mutation_filtering(predictions, mapping_dict, is_mutant, threshold=0.0, yes_only=False, tool_type='netphos', failure_map=None):
    """Universal prediction filtering logic for both NetPhos and NetNGlyc

    Args:
        predictions: List of all predictions
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant (True) or wildtype (False) sequences
        threshold: Score threshold for filtering
        yes_only: Only include predictions with 'YES' answer (NetPhos only)
        tool_type: 'netphos' or 'netnglyc' for tool-specific handling

    Returns:
        list: Filtered predictions with pkeys
    """
    results = []

    failure_map = failure_map or {}

    if is_mutant:
        # Group predictions by sequence name for single-mutation processing
        seq_predictions = {}
        for pred in predictions:
            seq_name = pred.get('Gene', '') if tool_type == 'netphos' else pred.get('seq_name', '')
            if seq_name not in seq_predictions:
                seq_predictions[seq_name] = []
            seq_predictions[seq_name].append(pred)

        # Process each sequence separately with its specific mutation
        for seq_name, seq_preds in seq_predictions.items():
            seq_results = process_single_mutation_for_sequence(
                seq_name, seq_preds, mapping_dict, is_mutant=True, tool_type=tool_type,
                failure_map=failure_map
            )

            # Apply additional filters
            for result in seq_results:
                # Apply threshold filter
                score_field = 'score' if tool_type == 'netphos' else 'potential'
                if result[score_field] < threshold:
                    continue

                if should_skip_mutation(result['Gene'], result.get('pkey', '').split('-')[-1] if 'pkey' in result else None, failure_map):
                    continue

                # Apply yes_only filter (NetPhos only)
                if yes_only and tool_type == 'netphos' and result.get('answer') != 'YES':
                    continue

                results.append(result)

    else:
        # Wildtype processing - use existing bulk logic (not implemented here)
        # This should use the existing wildtype processing logic from each pipeline
        raise NotImplementedError("Wildtype processing should use existing pipeline-specific logic")

    return results

def strip_all_extensions(name: str) -> str:
    # remove .csv, .csv.gz, etc.
    return re.sub(r'(\.[^.]+)+$', '', name)

def is_likely_gene_name(name: str) -> bool:
    pat = re.compile(r'^[A-Za-z0-9][A-Za-z0-9\-]{1,14}$')
    # PDB ID with chain: 1ABC_A, 7RM1_B
    pdb_chain = re.compile(r'^\d[A-Za-z0-9]{3}_[A-Za-z0-9]{1,2}$')
    return 2 <= len(name) <= 15 and (bool(pat.match(name)) or bool(pdb_chain.match(name)))

def extract_gene_from_filename(filename: str) -> str:
    """Return the most likely gene symbol from a filename."""
    name = Path(filename).name
    name = strip_all_extensions(name)

    # PDB ID with optional chain: 1ABC, 1ABC_A, 7RM1_B
    pdb = re.match(r'^(\d[A-Za-z0-9]{3}(?:_[A-Za-z0-9]{1,2})?)(?:_|$)', name)
    if pdb:
        return pdb.group(1)

    # Fast path: choose the last gene-like token bounded by underscores
    matches = re.findall(r'(?:^|_)([A-Z0-9][A-Z0-9]{1,14}(?:-[A-Z0-9]+)?)(?=$|_)', name)
    for cand in reversed(matches):
        if is_likely_gene_name(cand):
            return cand

    # Token-based cleanup for flexible prefixes/suffixes like 'transcript_mapping_'
    tokens = [t for t in name.split('_') if t]

    # strip leading stopwords
    i = 0
    while i < len(tokens) and tokens[i].lower() in STOPWORDS:
        i += 1
    # strip trailing stopwords
    j = len(tokens)
    while j > i and tokens[j-1].lower() in STOPWORDS:
        j -= 1
    core = tokens[i:j] if i < j else tokens

    # pick best candidate from remaining tokens
    candidates = [t for t in core if is_likely_gene_name(t)]
    if candidates:
        def rank(t):
            score = 0
            if t.isupper(): score += 1
            if any(ch.isdigit() for ch in t): score += 2
            if '-' in t: score += 2
            return (score, len(t))
        candidates.sort(key=rank)
        return candidates[-1]

    # Last resort: first gene-like substring anywhere
    m = re.search(r'([A-Za-z0-9][A-Za-z0-9\-]{1,14})', name)
    if m and is_likely_gene_name(m.group(1)):
        return m.group(1)

    # Fallback: basename without extensions
    return strip_all_extensions(Path(filename).stem)

def discover_mapping_files(mapping_dir):
    """Scan directory (or accept a single CSV) for mapping files and extract gene names flexibly.

    Args:
        mapping_dir: Directory path containing mapping CSV files, or path to a single CSV file

    Returns:
        dict: {gene_name: file_path} mapping
    """
    from pathlib import Path

    mapping_files = {}

    if not mapping_dir or not Path(mapping_dir).exists():
        return mapping_files

    # validate_mapping_content returns a BARE False for a single-column or
    # unreadable file and a [ok, delimiter] pair otherwise, so subscripting it
    # raises TypeError on exactly the files it means to reject -- the
    # {GENE}_mutations.csv copies sitting in the same tree. The except arm below
    # swallowed that as "Skipping", which is the right outcome reached by the
    # wrong route and hides a real read error behind the same message.
    # load_mapping already guards this (tagged F37); this did not.
    def _accepts(path) -> bool:
        res = validate_mapping_content(path)
        return bool(res[0]) if isinstance(res, (list, tuple)) else bool(res)

    p = Path(mapping_dir)
    if p.is_file():
        gene_name = extract_gene_from_filename(p.stem)
        if _accepts(p):
            mapping_files[gene_name] = str(p)
        return mapping_files

    # Scan for all CSV files
    for csv_file in Path(mapping_dir).rglob("*.csv"):
        try:
            # Extract gene name from filename
            gene_name = extract_gene_from_filename(csv_file.stem)

            # Validate CSV content structure
            if _accepts(csv_file):
                mapping_files[gene_name] = str(csv_file)

        except Exception as e:
            # Skip files that can't be processed
            print(f"Warning: Skipping {csv_file}: {e}")
            continue

    return mapping_files

def discover_fasta_files(fasta_dir):
    """Scan directory for FASTA files with flexible extensions and extract gene names

    Args:
        fasta_dir: Directory path containing FASTA files

    Returns:
        dict: {gene_name: file_path} mapping
    """
    from pathlib import Path

    fasta_files = {}

    if not fasta_dir or not Path(fasta_dir).exists():
        return fasta_files

    # Common FASTA file extensions
    fasta_extensions = ['*.fasta', '*.fa', '*.fas', '*.fna', '*.faa']

    for extension in fasta_extensions:
        for fasta_file in Path(fasta_dir).rglob(extension):
            try:
                # Extract gene name from filename
                gene_name = extract_gene_from_filename(fasta_file.stem)

                # Validate FASTA content
                if validate_fasta_content(fasta_file):
                    # If gene already found, prefer more specific naming
                    if gene_name not in fasta_files or '_aa' in fasta_file.stem:
                        fasta_files[gene_name] = str(fasta_file)

            except Exception as e:
                # Skip files that can't be processed
                print(f"Warning: Skipping {fasta_file}: {e}")
                continue

    return fasta_files

def validate_mapping_content(file_path):
    """Validate that CSV file has the expected mapping structure"""
    import csv

    try:
        with open(file_path, 'r') as f:
            # Read a sample to check for delimiters
            sample = f.read(1024)
            f.seek(0)

            # Check if file has delimiters
            has_delimiters = any(char in sample for char in ',\t;|')

            if has_delimiters:
                # Two-column file: detect delimiter
                sniffer = csv.Sniffer()
                delimiter = sniffer.sniff(sample, delimiters=',\t').delimiter
                reader = csv.DictReader(f, delimiter=delimiter)
            else:
                # Single-column file: no delimiter needed
                reader = csv.DictReader(f)

            # Check for required columns
            fieldnames = [field.lower() for field in reader.fieldnames] if reader.fieldnames else []
            has_mutation = any(col in fieldnames for col in ['mutant', 'mutation', 'nt_mutation', 'ntmutant'])
            # 'intron' and 'pre_mrna' are the value columns of the two mapping
            # CSVs exon_aware_mapping emits for intronic variants. fieldnames are
            # lower-cased above, so the pre_mRNA header arrives as 'pre_mrna'.
            # Adding to this allow-list can only ACCEPT files that were previously
            # rejected; it cannot reject anything that passes today.
            # 'pkey' is the value column of the mutant->pkey inversion table
            # exon_aware_mapping emits alongside the coordinate mappings. Without
            # it a pkey,mutant file is rejected as single-column and the four
            # pipelines that must invert a FASTA header back to a token cannot
            # load their lookup file.
            has_aa_mutation = any(col in fieldnames for col in
                                  ['aamutant', 'transcript', 'genomic', 'aa_mutation', 'amino_acid_mutation',
                                   'protein_mutation', 'chromosome', 'intron', 'pre_mrna', 'pkey'])

            # Valid formats
            return [(has_mutation and has_aa_mutation), delimiter] or [(has_mutation and not has_aa_mutation), delimiter]

    except Exception:
        return False

def validate_fasta_content(file_path):
    """Validate that file contains valid FASTA format

    Args:
        file_path: Path to FASTA file

    Returns:
        bool: True if valid FASTA file, False otherwise
    """
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()

            # Must start with '>' for FASTA format
            if not first_line.startswith('>'):
                return False

            # Check that there's sequence content
            second_line = f.readline().strip()
            if not second_line or second_line.startswith('>'):
                return False

            # Basic sequence validation (should contain valid amino acid or nucleotide characters)
            valid_chars = set('ACDEFGHIKLMNPQRSTVWYXZUOB*-')  # Amino acids + ambiguous
            if not any(char.upper() in valid_chars for char in second_line):
                return False

            return True

    except Exception:
        return False

def load_wt_sequences(input_dir: str, wt_header: str = "transcript") -> Dict[str, str]:
    """
    Load WT sequences (configured by wt_header) into memory keyed by gene symbol.
    Accepts a single FASTA file or a directory of FASTA files.
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"WT fasta path not found: {input_dir}")

    sequences: Dict[str, str] = {}
    if input_path.is_file():
        fasta_files = [input_path] if input_path.suffix.lower() in (".fa", ".fasta", ".fna") else []
    else:
        fasta_files = sorted([f for f in input_path.iterdir() if f.is_file() and f.suffix.lower() in (".fa", ".fasta", ".fna")])
    print(f"[WT] Scanning {len(fasta_files)} WT FASTA files")
    for fasta_file in fasta_files:
        data = read_fasta(str(fasta_file))
        if not data:
            continue

        seq = None
        if wt_header in data and data[wt_header] and data[wt_header].strip():
            seq = data[wt_header].strip()
        else:
            non_empty = [(h, s) for h, s in data.items() if s and s.strip()]
            if len(non_empty) == 1:
                # Single record: no ambiguity — use it regardless of header label.
                seq = non_empty[0][1].strip()
            elif len(non_empty) > 1:
                # F39: >1 record and none matches wt_header. The old max-by-length
                # pick silently seeded WT/MUT from an arbitrary isoform. Refuse to
                # guess — skip this file loudly so the gene is absent, not wrong.
                print(f"[WT][SKIP] {fasta_file.name}: header '{wt_header}' not found among "
                      f"{len(non_empty)} records ({', '.join(h for h, _ in non_empty)}); "
                      f"cannot disambiguate WT isoform. Pass --wt-header to select one. Skipping.")
                continue

        if not seq:
            continue

        gene = extract_gene_from_filename(fasta_file.name) or fasta_file.stem
        sequences[gene.upper()] = seq

    header_label = wt_header.upper() if wt_header else "SEQUENCES"
    print(f"[WT] Loaded {len(sequences)} WT {header_label} into memory")
    return sequences


def resolve_output_base(output_arg, input_arg, tool_name):
    """Resolve the output base path as a nested per-gene/per-tool directory.

    Always treats *output_arg* as a base directory and constructs:
        {output_arg}/{gene}/{tool_name}/{gene}

    The nested directory is created automatically.

    Args:
        output_arg: The raw output argument from argparse (always a directory).
        input_arg:  The raw input argument (FASTA path / directory).
        tool_name:  Tool subdirectory name (e.g. 'NetPhos', 'NetMHC').

    Returns:
        str: A resolved output base path suitable for appending '.tsv' etc.
    """
    out = Path(output_arg)
    gene = extract_gene_from_filename(Path(input_arg).stem) or Path(input_arg).stem
    nested = out / gene / tool_name
    nested.mkdir(parents=True, exist_ok=True)
    return str(nested / gene)


# =============================================================================
# Shared FASTA Synthesis Functions
# =============================================================================

def translate_orf_sequence(nt_sequence: str) -> str:
    """Translate a nucleotide ORF into an amino acid sequence, trimming trailing stops."""
    if not nt_sequence:
        return ""
    cleaned = nt_sequence.strip().upper().replace("U", "T")
    if not cleaned:
        return ""
    aa_seq = str(Seq(cleaned).translate(to_stop=False))
    return aa_seq.rstrip('*').strip()


def load_wt_sequence_map(input_path: str, wt_header: str = "ORF"):
    """
    Load WT nucleotide sequences from a file or directory using shared utility helpers.

    Returns:
        tuple(dict, tempfile.TemporaryDirectory | None): (gene -> nt sequence, temp dir holder)
    """
    import tempfile as _tmpmod
    source = Path(input_path)
    temp_dir = None

    if source.is_dir():
        load_path = str(source)
    elif source.is_file():
        temp_dir = _tmpmod.TemporaryDirectory(prefix="wt_src_")
        shutil.copy2(str(source), os.path.join(temp_dir.name, source.name))
        load_path = temp_dir.name
    else:
        raise FileNotFoundError(f"Input FASTA path not found: {input_path}")

    sequences = load_wt_sequences(load_path, wt_header=wt_header)
    return sequences, temp_dir


def infer_aamutation_from_nt(mutant_id: str, nt_sequence: str):
    """Infer the amino-acid mutation string (e.g., K543E) from a nucleotide mutation."""
    nt_info = get_mutation_data_bioAccurate(mutant_id, is_nt=True)
    if nt_info[0] is None:
        return None
    aa_info = get_mutant_aa(nt_info, nt_sequence)
    if not aa_info:
        return None
    (aa_pos, (wt_aa, mut_aa)), _ = aa_info
    # Skip nonsense (stop-gain) variants: codon_to_aa maps stop codons to the
    # 4-char string 'Stop', which must NOT be spliced in as a residue. Matches
    # the stop-skip already done for aa tokens in get_mutation_data_bioAccurate.
    if not wt_aa or not mut_aa or mut_aa == 'Stop':
        return None
    aa_pos = int(aa_pos)
    return aa_pos, wt_aa, mut_aa


def infer_aavariant_from_nt(mutant_id: str, nt_sequence: str):
    """Length-aware sibling of infer_aamutation_from_nt.

    Same inputs (nt token + WT ORF) and the same first three return values, but
    wt_aa/mut_aa may be MULTI-RESIDUE strings and a fourth element carries the
    consequence class:

        (aa_pos, wt_aa, mut_aa, aa_consequence)   or   None

    infer_aamutation_from_nt is the SNV-only original and is left untouched: it
    routes through get_mutation_data_bioAccurate, which raises on any multi-base
    token, and it returns None for stop-gain. Three tests pin that behaviour
    (test_utils_sequence.py::TestInferAamutationFromNt), so this is additive.

    This is the shared chokepoint for indel support in the four protein
    pipelines: gitnexus reports infer_aamutation_from_nt as CRITICAL with 22
    impacted symbols, reached twice per pipeline -- once via
    build_mutant_sequences_for_gene when synthesizing the mutant protein, and
    once when a mapping CSV lacks an `aamutant` column and the token has to be
    derived (netnglyc_pipeline.py:1233). An indel returns None from the original,
    and the netNglyc caller then drops the mutation outright at :1237-1238.

    Unlike the original this does NOT drop stop-gain -- a premature stop is a
    real consequence and is reported as 'stop_gained'.
    """
    variant = parse_variant(mutant_id, is_nt=True)
    if variant is None or not nt_sequence:
        return None
    cons = protein_consequence(variant, nt_sequence)
    if cons is None:
        return None
    aa_pos, wt_aa, mut_aa = cons['aa_pos'], cons['wt_aa'], cons['mut_aa']

    if cons['aa_consequence'] == 'synonymous':
        # _trim_aa_span empties BOTH alleles here because nothing changed, not
        # because an allele is absent. Letting that fall through to the indel
        # anchoring below re-anchored on the PRECEDING residue and returned
        # aa_pos - 1, so every synonymous variant was reported one codon early
        # (nt 13 is codon 5 'CTG'/L and came back as codon 4 'TGG'/W). The
        # position from protein_consequence is already right; only the residue
        # has to be filled back in.
        prot_full = _translate_codons(nt_sequence)
        idx = aa_pos - 1
        residue = prot_full[idx] if 0 <= idx < len(prot_full) else ''
        return (aa_pos, residue, residue, 'synonymous')

    # protein_consequence returns the MINIMAL (trimmed) form, so a pure deletion
    # has mut_aa == '' and a pure insertion has wt_aa == ''. An empty allele
    # cannot be rendered as a token, and a caller that receives one drops the
    # mutation -- the exact silent drop this replaces. Re-anchor on the preceding
    # residue, the VCF convention already used for nucleotides, so both alleles
    # are non-empty and parse_variant(is_nt=False) can read the token back.
    if (not wt_aa or not mut_aa) and cons['aa_consequence'] != 'frameshift':
        prot = _translate_codons(nt_sequence).split('*')[0]
        idx = aa_pos - 1
        if idx - 1 >= 0 and idx - 1 < len(prot):          # anchor 5' (preferred)
            anchor = prot[idx - 1]
            aa_pos, wt_aa, mut_aa = aa_pos - 1, anchor + wt_aa, anchor + mut_aa
        elif idx < len(prot):                              # at residue 1: anchor 3'
            anchor = prot[idx + len(cons['wt_aa'])] if idx + len(cons['wt_aa']) < len(prot) else ''
            if anchor:
                wt_aa, mut_aa = wt_aa + anchor, mut_aa + anchor
    return (aa_pos, wt_aa, mut_aa, cons['aa_consequence'])


def format_aa_token(aa_pos, wt_aa, mut_aa, aa_consequence=None):
    """Render an aa-level change as a mapping-CSV `aamutant` token.

    'K541E' for a substitution -- byte-identical to what
    f"{wt_aa}{aa_pos}{mut_aa}" produces today, so existing rows are unchanged.
    'KEW541R' for a multi-residue in-frame change, which parse_variant(is_nt=False)
    reads back correctly and get_mutation_data_bioAccurate does not.

    A frameshift has no bounded wt_aa/mut_aa (see protein_consequence), so it
    renders in the HGVS style 'K541fs'. That form is deliberately NOT parseable by
    parse_variant: there is no residue pair to recover, and inventing one would be
    the fabrication this whole layer exists to prevent. Callers must branch on
    aa_consequence rather than trying to re-parse it.
    """
    if aa_consequence == 'frameshift':
        # No residue pair exists; carry the position and the class only.
        return f"{(wt_aa or '')[:1]}{aa_pos}fs"
    if not wt_aa or not mut_aa:
        # infer_aavariant_from_nt anchors these, so an empty allele reaching here
        # means the caller built the tuple by hand without a WT protein. Refuse
        # rather than emit '' -- an empty token is dropped downstream, silently.
        raise ValueError(
            f"cannot format an aa token with an empty allele "
            f"(pos={aa_pos}, wt={wt_aa!r}, alt={mut_aa!r}); anchor it first")
    return f"{wt_aa}{aa_pos}{mut_aa}"


def _non_snv_mutant_protein(gene_name, token, nt_sequence, non_snp):
    """Build the mutant PROTEIN for a non-SNV token. Returns (header, seq) or None.

    Returns None when the caller should fall through to the existing SNV path:
    the token is an SNV, it does not parse, --non-snp is off, or there is no
    nucleotide sequence to splice.

    An indel's protein consequence cannot be produced by substituting a residue.
    update_str(aa_sequence, mut_aa, idx) replaces exactly one amino acid, which
    is right for a missense SNV and meaningless for a length change. The mutant
    protein is the translation of the edited ORF, so this splices at the
    NUCLEOTIDE level and retranslates.

    A frameshift is built, not refused: translation to the new stop is
    well-defined and the DTU tools can score the resulting protein. What is not
    well-defined is any position-wise WT<->MUT comparison, and that is the
    consumer's problem to suppress -- see protein_consequence's aa_consequence.

    Raises ValueError on a REF that disagrees with the ORF, so the caller's
    per-token handler records it by name instead of dropping it silently.
    """
    v = parse_variant(token, is_nt=True)
    if v is None or v.is_snv or not non_snp or not nt_sequence:
        return None
    if v.pos0 + len(v.ref) > len(nt_sequence):
        raise ValueError(
            f"{token}: REF spans past the ORF ({v.pos0 + len(v.ref)} > {len(nt_sequence)})")
    observed = nt_sequence[v.pos0:v.pos0 + len(v.ref)].upper()
    if observed != v.ref.upper():
        raise ValueError(
            f"{token}: REF mismatch, ORF has {observed!r} at 0-based {v.pos0}")
    mut_nt = splice_seq(nt_sequence, v.pos0, v.ref, v.alt, validate=False)
    # Truncate at the FIRST stop, not just trailing ones. translate_orf_sequence
    # uses to_stop=False and .rstrip('*'), which is fine for an SNV but leaves an
    # internal '*' embedded after a frameshift -- a real run produced 'MKNG*PVI'.
    # Handing that to a DTU tool scores four residues that translation never
    # reaches. A frameshift ends at its new stop.
    #
    # Trim to a whole number of codons first: a frameshift leaves a partial
    # trailing codon, which makes Biopython emit a BiopythonWarning.
    mut_nt = mut_nt[:len(mut_nt) - len(mut_nt) % 3]
    mut_aa_seq = translate_orf_sequence(mut_nt).split('*')[0]
    if not mut_aa_seq:
        raise ValueError(
            f"{token}: mutant ORF translates to nothing before the first stop")
    return f"{gene_name}-{token}", mut_aa_seq


def build_mutant_sequences_for_gene(
    gene_name: str,
    nt_sequence,
    aa_sequence: str,
    mapping_file,
    log_path,
    failure_map,
    input_type: str = 'nt',
    non_snp: bool = False,
):
    """
    Return a dict of {header: sequence} for all mutants of a given gene.

    Handles both CSV format (with headers) and single-column format.

    Args:
        gene_name: Gene identifier
        nt_sequence: Wild-type nucleotide sequence (None if input_type='aa')
        aa_sequence: Wild-type amino acid sequence
        mapping_file: Path to mutation file (single-column or CSV)
        log_path: Optional validation log path
        failure_map: Optional map of failed mutations to skip
        input_type: 'nt' for nucleotide mutations (e.g., A1002T), 'aa' for amino acid mutations (e.g., M334V)
    """
    if not mapping_file or not os.path.exists(mapping_file):
        return {}

    allowed_mutations = None
    if log_path:
        try:
            allowed_mutations = {
                entry.split(',')[0].strip().upper()
                for entry in trim_muts(mapping_file, log=log_path, gene_name=gene_name)
                if entry
            }
        except Exception:
            allowed_mutations = None

    mutant_sequences = {}
    # Per-token failures are collected, not raised. A malformed token must cost
    # its own row and nothing else -- the outer try below still returns {} for
    # file-level failures (unreadable file, malformed CSV header), which is a
    # genuinely gene-wide problem.
    skipped_tokens = []
    # Intronic tokens are gated HERE, at the single chokepoint through which
    # netNglyc, netphos, netMHC and NetSurfP3 all obtain their mutant proteins
    # (gitnexus: direct callers NetSurfP3.main, netMHC.main, synthesize_gene_fastas;
    # netNglyc calls it from _synthesize_gene_fastas_non_snv). Gating once here
    # rather than four times in the pipelines keeps the four from drifting, and it
    # sits well before any DTU/Docker invocation -- a warning emitted after the
    # binary has already been handed a sequence is not a gate.
    intronic_tokens = []
    try:
        with open(mapping_file, 'r') as handle:
            lines = handle.readlines()

        # Detect format: single-column or CSV
        is_single_column = True
        if lines and ',' in lines[0]:
            first_line_lower = lines[0].lower()
            if any(keyword in first_line_lower for keyword in ['mutant', 'mutation', 'aamutant']):
                is_single_column = False

        mutant_keys = ['mutant', 'mutation', 'nt_mutation', 'ntmutant']
        aa_keys = ['aamutant', 'aa_mutation', 'amino_acid_mutation', 'protein_mutation']

        if is_single_column:
            for line in lines:
                mutant_id = line.strip()
                if not mutant_id or mutant_id.lower() == 'mutant':
                    continue

                try:
                    mutant_clean = mutant_id.replace(" ", "")
                    if is_intronic_token(mutant_clean):
                        intronic_tokens.append(mutant_clean)
                        continue
                    if allowed_mutations and mutant_clean.upper() not in allowed_mutations:
                        continue
                    if should_skip_mutation(gene_name, mutant_clean, failure_map):
                        continue

                    built = _non_snv_mutant_protein(gene_name, mutant_clean,
                                                    nt_sequence, non_snp)
                    if built is not None:
                        mutant_sequences[built[0]] = built[1]
                        continue

                    pos = None
                    wt_aa = mut_aa = None

                    if input_type == 'aa':
                        aa_info = get_mutation_data_bioAccurate(mutant_clean, is_nt=False)
                        if aa_info[0] is not None and aa_info[1]:
                            pos = aa_info[0]
                            wt_aa, mut_aa = aa_info[1]
                    else:
                        inferred = infer_aamutation_from_nt(mutant_clean, nt_sequence)
                        if inferred is not None:
                            pos, wt_aa, mut_aa = inferred

                    if pos is None or not wt_aa or not mut_aa:
                        continue

                    idx = int(pos) - 1
                    if idx < 0 or idx >= len(aa_sequence):
                        continue
                    if wt_aa and aa_sequence[idx].upper() != wt_aa.upper():
                        continue

                    header = f"{gene_name}-{mutant_clean}"
                    mutant_sequences[header] = update_str(aa_sequence, mut_aa, idx)
                except Exception as exc:
                    skipped_tokens.append((mutant_id, f"{type(exc).__name__}: {exc}"))
                    continue

        else:
            with open(mapping_file, 'r') as handle:
                reader = csv.DictReader(handle)

                for row in reader:
                    mutant_id = ""
                    for key in mutant_keys:
                        if key in row and row[key]:
                            mutant_id = row[key].strip()
                            break
                    if not mutant_id:
                        continue

                    try:
                        mutant_clean = mutant_id.replace(" ", "")
                        if is_intronic_token(mutant_clean):
                            intronic_tokens.append(mutant_clean)
                            continue
                        if allowed_mutations and mutant_clean.upper() not in allowed_mutations:
                            continue
                        if should_skip_mutation(gene_name, mutant_clean, failure_map):
                            continue

                        built = _non_snv_mutant_protein(gene_name, mutant_clean,
                                                        nt_sequence, non_snp)
                        if built is not None:
                            mutant_sequences[built[0]] = built[1]
                            continue

                        aa_string = ""
                        for key in aa_keys:
                            if key in row and row[key]:
                                aa_string = row[key].strip()
                                break

                        pos = None
                        wt_aa = mut_aa = None
                        if aa_string:
                            pos, nts = get_mutation_data_bioAccurate(aa_string, is_nt=False)
                            if pos is not None and nts:
                                wt_aa, mut_aa = nts

                        if pos is None or not wt_aa or not mut_aa:
                            if input_type == 'aa':
                                aa_info = get_mutation_data_bioAccurate(mutant_clean, is_nt=False)
                                if aa_info[0] is not None and aa_info[1]:
                                    pos = aa_info[0]
                                    wt_aa, mut_aa = aa_info[1]
                            else:
                                inferred = infer_aamutation_from_nt(mutant_clean, nt_sequence)
                                if inferred is not None:
                                    pos, wt_aa, mut_aa = inferred

                        if pos is None or not wt_aa or not mut_aa:
                            continue

                        idx = int(pos) - 1
                        if idx < 0 or idx >= len(aa_sequence):
                            continue
                        if wt_aa and aa_sequence[idx].upper() != wt_aa.upper():
                            continue

                        header = f"{gene_name}-{mutant_clean}"
                        mutant_sequences[header] = update_str(aa_sequence, mut_aa, idx)
                    except Exception as exc:
                        skipped_tokens.append((mutant_id, f"{type(exc).__name__}: {exc}"))
                        continue

    except Exception as exc:
        # File-level failure only. Per-token failures never reach here.
        print(f"Warning: Failed to synthesize mutants for {gene_name} ({mapping_file}): {exc}")
        return {}

    warn_intronic_unsupported(
        'protein_synthesis', gene_name, intronic_tokens,
        "An intron has no reading frame and no residue, so no mutant protein "
        "exists to score. This excludes them from netNglyc, netphos, netMHC, "
        "NetSurfP3 and EVmutation. Score intronic variants with RNAfold, miranda, "
        "genesplicer or AlphaFold3 instead.")

    if skipped_tokens:
        shown = ", ".join(f"{tok} ({why})" for tok, why in skipped_tokens[:5])
        more = f" (+{len(skipped_tokens) - 5} more)" if len(skipped_tokens) > 5 else ""
        print(
            f"Warning: {gene_name}: skipped {len(skipped_tokens)} unparseable "
            f"mutation(s), kept {len(mutant_sequences)}: {shown}{more}"
        )

    return mutant_sequences


def synthesize_gene_fastas(wt_sequences, mapping_lookup, sequence_root, log_path=None, failure_map=None):
    """Create WT and mutant FASTAs for each gene, returning the directories and summary."""
    sequence_root = Path(sequence_root)
    wt_dir = sequence_root / "wt"
    mut_dir = sequence_root / "mut"
    wt_dir.mkdir(parents=True, exist_ok=True)
    mut_dir.mkdir(parents=True, exist_ok=True)

    summary = []
    for gene_name, wt_seq in wt_sequences.items():
        gene_name = gene_name.upper()
        seq_upper = wt_seq.strip().upper()

        # Auto-detect the WT alphabet: nucleotide -> translate to protein;
        # protein -> use directly (no translation); codon-encoded -> skip.
        # Prevents the unconditional-translate crash (TranslationError: invalid
        # codon) that amino-acid input previously caused here.
        try:
            detected = detect_alphabet(seq_upper)
        except ValueError:
            print(f"Skipping {gene_name}: empty sequence")
            continue
        if detected == 'nucleotide':
            nt_for_build = seq_upper
            aa_seq = translate_orf_sequence(seq_upper)
            build_input_type = 'nt'
            if not aa_seq:
                print(f"Skipping {gene_name}: unable to translate ORF")
                continue
        elif detected == 'protein':
            nt_for_build = None
            aa_seq = seq_upper
            build_input_type = 'aa'
        else:  # codon-encoded
            print(f"Skipping {gene_name}: codon-encoded input, expected nt or aa")
            continue

        wt_header = f"{gene_name}-wt"
        wt_path = wt_dir / f"{gene_name}-wt.fasta"
        write_fasta(wt_path, {wt_header: aa_seq})

        mapping_file = mapping_lookup.get(gene_name.upper())
        # non_snp=True is REQUIRED here, not optional. This is the only route by
        # which netNglyc, netphos, netMHC and NetSurfP3 obtain their mutant
        # proteins, and without it build_mutant_sequences_for_gene never reaches
        # _non_snv_mutant_protein -- every indel dies at synthesis, before any of
        # those pipelines' own non-SNV handling can see it. An SNV token is
        # unaffected: the non-SNV branch returns None for it and the original
        # path runs unchanged.
        mutant_sequences = build_mutant_sequences_for_gene(
            gene_name,
            nt_for_build,
            aa_seq,
            mapping_file,
            log_path,
            failure_map,
            input_type=build_input_type,
            non_snp=True,
        )

        mut_path = None
        if mutant_sequences:
            mut_path = mut_dir / f"{gene_name}_aa.fasta"
            write_fasta(mut_path, mutant_sequences)

        summary.append({
            "gene": gene_name,
            "wt_path": str(wt_path),
            "mut_path": str(mut_path) if mut_path else None,
            "mutant_count": len(mutant_sequences),
        })

    return wt_dir, mut_dir, summary


# =============================================================================
# Codon Usage Functions
# =============================================================================

def get_codon_counts(seq):
    """
    Compute codon and codon-pair statistics from a nucleotide sequence.

    Returns:
        tuple: (codondata, codonpairdata) dictionaries containing:
            - codondata['counts']: Raw codon counts
            - codondata['RSCU']: Relative Synonymous Codon Usage
            - codondata['W']: Relative adaptiveness (codon frequency / max synonymous frequency)
            - codonpairdata['counts']: Raw bicodon counts
            - codonpairdata['RSCPU']: Relative Synonymous Codon Pair Usage
            - codonpairdata['CPS']: Codon Pair Score (ln of observed/expected)
            - codonpairdata['noln CPS']: CPS without natural log
            - codonpairdata['W_CP']: Relative adaptiveness for codon pairs
    """
    import numpy as np

    codondata = {
        "counts": {codon: 0 for codon in codon_to_aa.keys()},
        "RSCU": {codon: 0 for codon in codon_to_aa.keys()}
    }
    codonpairdata = {
        "counts": {codon1 + codon2: 0 for codon1 in codon_to_aa.keys() for codon2 in codon_to_aa.keys()},
        "RSCPU": {codon1 + codon2: 0 for codon1 in codon_to_aa.keys() for codon2 in codon_to_aa.keys()},
        "noln CPS": {},
        "CPS": {}
    }

    # Ensure sequence length is a multiple of 3
    seq_len_multiple_of_3 = (len(seq) // 3) * 3

    for i in range(0, seq_len_multiple_of_3, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and codon in codondata["counts"]:
            codondata["counts"][codon] += 1

        if i + 6 <= seq_len_multiple_of_3:
            bicodon = seq[i:i + 6]
            if bicodon in codonpairdata["counts"]:
                codonpairdata["counts"][bicodon] += 1

    # Calculate RSCU for each codon
    for codon1 in codondata["counts"].keys():
        aa = codon_to_aa.get(codon1)
        if not aa or aa == 'Stop' or aa == '-':
            codondata["RSCU"][codon1] = np.nan
            continue
        syn_codons = codon_table.get(aa, [])
        numsyn = sum(codondata["counts"].get(c, 0) for c in syn_codons)
        try:
            codondata["RSCU"][codon1] = codondata["counts"][codon1] / (numsyn / len(syn_codons))
        except ZeroDivisionError:
            codondata["RSCU"][codon1] = np.nan

    # Calculate W (relative adaptiveness) for codons
    codondata["W"] = {}
    for codon1 in codondata["RSCU"].keys():
        aa = codon_to_aa.get(codon1)
        if not aa or aa == 'Stop' or aa == '-' or np.isnan(codondata["RSCU"].get(codon1, np.nan)):
            codondata["W"][codon1] = np.nan
            continue
        syn_codons = codon_table.get(aa, [])
        max_rscu = max([codondata["RSCU"].get(c, 0) for c in syn_codons] or [1])
        codondata["W"][codon1] = codondata["RSCU"][codon1] / max_rscu if max_rscu > 0 else np.nan

    # Calculate RSCPU for codon pairs
    for cp1 in codonpairdata["counts"].keys():
        codon_a, codon_b = cp1[:3], cp1[3:]
        aa1 = codon_to_aa.get(codon_a)
        aa2 = codon_to_aa.get(codon_b)

        if not aa1 or not aa2 or aa1 in ('Stop', '-') or aa2 in ('Stop', '-'):
            codonpairdata["RSCPU"][cp1] = np.nan
            continue

        syn_cp_codons1 = codon_table.get(aa1, [])
        syn_cp_codons2 = codon_table.get(aa2, [])

        numsyn = sum(codonpairdata["counts"].get(cpa + cpb, 0)
                     for cpa in syn_cp_codons1 for cpb in syn_cp_codons2)
        syn_cp_count = len(syn_cp_codons1) * len(syn_cp_codons2)

        try:
            codonpairdata["RSCPU"][cp1] = codonpairdata["counts"][cp1] / (numsyn / syn_cp_count)
        except ZeroDivisionError:
            codonpairdata["RSCPU"][cp1] = np.nan

        # Calculate CPS
        try:
            expected_count = codondata["counts"].get(codon_a, 0) * codondata["counts"].get(codon_b, 0)
            if expected_count == 0:
                raise ZeroDivisionError
            noln_cps = codonpairdata["counts"][cp1] / expected_count
            codonpairdata["noln CPS"][cp1] = noln_cps
            codonpairdata["CPS"][cp1] = math.log(noln_cps)
        except (ZeroDivisionError, ValueError):
            codonpairdata["noln CPS"][cp1] = np.nan
            codonpairdata["CPS"][cp1] = np.nan

    # Calculate W_CP for codon pairs
    codonpairdata["W_CP"] = {}
    for cp1 in codonpairdata["RSCPU"].keys():
        if np.isnan(codonpairdata["RSCPU"].get(cp1, np.nan)):
            codonpairdata["W_CP"][cp1] = np.nan
            continue
        codon_a, codon_b = cp1[:3], cp1[3:]
        aa1 = codon_to_aa.get(codon_a)
        aa2 = codon_to_aa.get(codon_b)
        if not aa1 or not aa2:
            codonpairdata["W_CP"][cp1] = np.nan
            continue
        syn_cp_codons1 = codon_table.get(aa1, [])
        syn_cp_codons2 = codon_table.get(aa2, [])
        max_rscpu = max([codonpairdata["RSCPU"].get(cpa + cpb, 0)
                         for cpa in syn_cp_codons1 for cpb in syn_cp_codons2] or [1])
        codonpairdata["W_CP"][cp1] = codonpairdata["RSCPU"][cp1] / max_rscpu if max_rscpu > 0 else np.nan

    return codondata, codonpairdata


# Human tRNA adaptation weights (tAI)
# Based on tRNA gene copy numbers and wobble pairing efficiency
# Sources: dos Reis et al. 2004, Tuller et al. 2010
# Format: codon -> tAI weight (0-1 scale, normalized)
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

# Human reference W values for CAI calculation
# Based on highly expressed genes (Sharp & Li 1987, adapted for human)
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


def compute_cai(seq, w_values=None):
    """
    Compute Codon Adaptation Index (CAI) for a sequence.

    CAI = geometric mean of W values across all codons
    CAI = exp((1/L) * sum(ln(W_i)))

    Args:
        seq: Nucleotide sequence (in-frame ORF)
        w_values: Dict of codon -> W values. If None, uses HUMAN_REFERENCE_W

    Returns:
        float: CAI value (0-1), or None if cannot compute
    """
    import numpy as np

    if w_values is None:
        w_values = HUMAN_REFERENCE_W

    seq = seq.upper().replace('U', 'T')
    seq_len = (len(seq) // 3) * 3

    log_w_sum = 0.0
    codon_count = 0

    for i in range(0, seq_len, 3):
        codon = seq[i:i+3]
        if codon in ('TAA', 'TAG', 'TGA', '---'):  # Skip stops and gaps
            continue
        w = w_values.get(codon)
        if w is None or w <= 0:
            continue
        log_w_sum += np.log(w)
        codon_count += 1

    if codon_count == 0:
        return None

    return np.exp(log_w_sum / codon_count)


def compute_tai(seq, tai_weights=None):
    """
    Compute tRNA Adaptation Index (tAI) for a sequence.

    tAI = geometric mean of tRNA adaptation weights across all codons.

    Args:
        seq: Nucleotide sequence (in-frame ORF)
        tai_weights: Dict of codon -> tAI weights. If None, uses HUMAN_TAI_WEIGHTS

    Returns:
        float: tAI value (0-1), or None if cannot compute
    """
    import numpy as np

    if tai_weights is None:
        tai_weights = HUMAN_TAI_WEIGHTS

    seq = seq.upper().replace('U', 'T')
    seq_len = (len(seq) // 3) * 3

    log_tai_sum = 0.0
    codon_count = 0

    for i in range(0, seq_len, 3):
        codon = seq[i:i+3]
        if codon in ('TAA', 'TAG', 'TGA', '---'):  # Skip stops and gaps
            continue
        w = tai_weights.get(codon)
        if w is None or w <= 0:
            continue
        log_tai_sum += np.log(w)
        codon_count += 1

    if codon_count == 0:
        return None

    return np.exp(log_tai_sum / codon_count)


def get_codon_tai(codon, tai_weights=None):
    """Get tAI weight for a single codon."""
    if tai_weights is None:
        tai_weights = HUMAN_TAI_WEIGHTS
    return tai_weights.get(codon.upper().replace('U', 'T'), None)


def get_codon_cai_w(codon, w_values=None):
    """Get CAI W value (reference adaptiveness) for a single codon."""
    if w_values is None:
        w_values = HUMAN_REFERENCE_W
    return w_values.get(codon.upper().replace('U', 'T'), None)


def extract_codon_with_bicodons(ntposnt, seq):
    """
    Extract codon and bicodons for a given SNP, respecting biological constraints.

    Biology rules:
    - First codon: only forward bicodon possible (codon1 + codon2)
    - Last codon: only reverse bicodon possible (codon_n-1 + codon_n)
    - Middle codons: both forward and reverse bicodons possible

    Args:
        ntposnt (str): SNP notation (e.g., "A123G") - 1-based position
        seq (str): DNA sequence

    Returns:
        tuple: (original_codon, forward_bicodon, reverse_bicodon, pos_in_codon, pos, codon_number)
               where bicodons may be empty strings if not biologically possible
    """
    pos, mut = get_mutation_data_bioAccurate(ntposnt, is_nt=True)
    if pos is None:
        return None, "", "", 0, 0, 0

    # Convert to 0-based indexing
    pos_0_indexed = pos - 1

    # F46: apply the ALT base before slicing. Previously `mut` was bound and never
    # used, so every codon/bicodon (and every RSCU/W/CAI/tAI derived from them) was
    # sliced from the unmodified WT ORF while the caller labeled it 'mutated_codon'
    # — the mutation was invisible on the only CLI-reachable path. Only the mutated
    # codon changes; neighbouring codons in the bicodons stay WT, which is correct.
    if pos_0_indexed < 0 or pos_0_indexed >= len(seq):
        return None, "", "", 0, 0, 0
    ref_nt, alt_nt = mut
    # F46 REF GUARD: splicing ALT onto a base that is not REF invents a codon present
    # in NEITHER the WT nor the true mutant sequence. Measured on the production corpus
    # (59 genes / 35,089 tokens vs Bio_DBs/fastas ORF records): 1,537 tokens (4.38%)
    # have a REF that disagrees with the ORF base — concentrated in MECP2 295/423,
    # NOD2 271/387, FGFR2 215/290, MPZ 195/270 — i.e. an ORF/isoform mismatch, not noise.
    # Reject rather than fabricate; the caller turns None into a skipped row.
    if ref_nt and seq[pos_0_indexed].upper() != ref_nt.upper():
        return None, "", "", 0, 0, 0
    seq = seq[:pos_0_indexed] + alt_nt.upper() + seq[pos_0_indexed + 1:]

    pos_in_codon = pos_0_indexed % 3

    # Find the codon start position and codon number
    codon_start_pos = (pos_0_indexed // 3) * 3
    codon_number = (pos_0_indexed // 3) + 1  # 1-based codon numbering

    # Calculate total number of complete codons in sequence
    total_codons = len(seq) // 3

    # Extract the original codon containing the mutation
    original_codon = seq[codon_start_pos:codon_start_pos + 3]

    # Initialize bicodons
    forward_bicodon = ""
    reverse_bicodon = ""

    # Determine which bicodons are biologically possible
    is_first_codon = (codon_number == 1)
    is_last_codon = (codon_number == total_codons)

    # Forward bicodon (current codon + following codon)
    # Possible for first and middle codons, but not last codon
    if not is_last_codon and codon_start_pos + 6 <= len(seq):
        following_codon = seq[codon_start_pos + 3:codon_start_pos + 6]
        if len(following_codon) == 3:
            forward_bicodon = original_codon + following_codon

    # Reverse bicodon (preceding codon + current codon)
    # Possible for middle and last codons, but not first codon
    if not is_first_codon and codon_start_pos >= 3:
        preceding_codon = seq[codon_start_pos - 3:codon_start_pos]
        if len(preceding_codon) == 3:
            reverse_bicodon = preceding_codon + original_codon

    return original_codon, forward_bicodon, reverse_bicodon, pos_in_codon, pos, codon_number


# =============================================================================
# MSA Generation and Processing Utilities
# =============================================================================

# Single-character codon-encoded alphabet. MUST stay in lockstep with
# mutation_effects/bin/codon_encoding.py: 64 codons -> A-Z (26) + a-z (26) +
# 0-9 (10) + '!' '@' (2), plus '-' gap = 65 symbols. Case-sensitive by design
# ('A' encodes codon AAA, 'a' a different codon), which is why the codon check
# below runs BEFORE any upper-casing.
_CODON_ENCODE_CHARS = (
    [chr(c) for c in range(ord('A'), ord('Z') + 1)]
    + [chr(c) for c in range(ord('a'), ord('z') + 1)]
    + [str(d) for d in range(10)]
    + ['!', '@']
)
CODON_ENCODED_ALPHABET = '-' + ''.join(_CODON_ENCODE_CHARS)
# Symbols that occur ONLY in the codon-encoded alphabet — never in nucleotide
# (ACGTU/IUPAC) or standard protein (20 aa + BXZUO*) sequences. Their presence
# is an unambiguous codon-encoding signal.
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
                    base). IUPAC ambiguity codes are excluded — they collide with
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


def prepare_protein_query(query_fasta):
    """Return a protein FASTA path suitable for jackhmmer.

    If the query sequence is nucleotide, translates to protein (to stop codon)
    and writes the result to a temporary file.

    Returns:
        (protein_fasta_path, tmp_path_or_None)
        Caller must delete tmp_path when finished if it is not None.
    """
    import tempfile
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    record = next(SeqIO.parse(query_fasta, 'fasta'))
    alphabet = detect_alphabet(str(record.seq))

    if alphabet == 'protein':
        return query_fasta, None

    print(f"Nucleotide query detected for '{record.id}' — translating to protein.")
    aa_seq = record.seq.translate(to_stop=True)
    if not aa_seq:
        raise ValueError(f"Translation of '{record.id}' produced an empty protein sequence.")

    aa_record = SeqRecord(aa_seq, id=record.id, description='translated')
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    SeqIO.write([aa_record], tmp, 'fasta')
    tmp.close()
    print(f"Translated protein length: {len(aa_seq)} aa")
    return tmp.name, tmp.name


def run_jackhmmer(query_fasta, database, output_sto, jackhmmer_binary,
                  iterations=5, evalue_inclusion=1e-3, threads=4, max_seqs=10000):
    """
    Run jackhmmer iterative search against a sequence database.

    Args:
        query_fasta: Path to query protein sequence (FASTA)
        database: Path to UniRef90 or similar database
        output_sto: Path for Stockholm output
        jackhmmer_binary: Path to jackhmmer executable
        iterations: Number of search iterations
        evalue_inclusion: E-value threshold for inclusion
        threads: Number of CPU threads
        max_seqs: Maximum number of sequences to include in the alignment (default: 10000).
                  Prevents profile drift for genes in large superfamilies.

    Returns:
        Path to Stockholm output file
    """
    cmd = [
        jackhmmer_binary,
        '-N', str(iterations),
        '--incE', str(evalue_inclusion),
        '-A', output_sto,
        '--cpu', str(threads),
        '--noali',
        query_fasta,
        database
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"jackhmmer failed: {result.stderr}")

    return Path(output_sto)


def parse_stockholm(stockholm_file):
    """
    Parse Stockholm format MSA file.

    Args:
        stockholm_file: Path to Stockholm file

    Returns:
        tuple: (dict {seq_id: sequence}, rf_annotation str or None)
            rf_annotation is the concatenated #=GC RF string where 'x' marks
            match columns and '.' marks insert columns. None if absent.
    """
    from collections import defaultdict
    current_seqs = defaultdict(str)
    rf_parts = []

    with open(stockholm_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            # Capture RF column annotation before skipping other # lines
            if line.startswith('#=GC RF'):
                parts = line.split()
                if len(parts) >= 3:
                    rf_parts.append(parts[2])
                continue

            # Skip remaining comments and empty lines
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                current_seqs[seq_id] += seq

    rf_annotation = ''.join(rf_parts) if rf_parts else None
    return dict(current_seqs), rf_annotation


def stockholm_to_a2m(msa, focus_seq_id, rf_annotation=None):
    """
    Convert Stockholm MSA to A2M format.

    A2M format:
    - Uppercase: match states (aligned to query)
    - Lowercase: insertions relative to query
    - '-': deletions (gaps in sequence, not in query)
    - '.': gaps in query (insertions in other sequences)

    Args:
        msa: dict {seq_id: sequence}
        focus_seq_id: ID of the focus/query sequence
        rf_annotation: #=GC RF string from parse_stockholm where 'x' = match
            column and '.' = insert column. When provided, this is used as the
            authoritative match/insert column definition. Falls back to focus
            sequence non-gap positions when None.

    Returns:
        dict: {seq_id: a2m_sequence}
    """
    if focus_seq_id not in msa:
        for seq_id in msa:
            if focus_seq_id in seq_id or seq_id in focus_seq_id:
                focus_seq_id = seq_id
                break
        else:
            raise ValueError(f"Focus sequence '{focus_seq_id}' not found in MSA")

    # Identify match columns using RF annotation when available
    if rf_annotation is not None:
        match_columns = {i for i, c in enumerate(rf_annotation) if c == 'x'}
    else:
        focus_seq = msa[focus_seq_id]
        match_columns = {i for i, c in enumerate(focus_seq) if c not in '-.'}

    a2m_msa = {}
    for seq_id, seq in msa.items():
        a2m_seq = []
        for i, c in enumerate(seq):
            if i in match_columns:
                if c in '-.':
                    a2m_seq.append('-')
                else:
                    a2m_seq.append(c.upper())
            else:
                if c in '-.':
                    a2m_seq.append('.')
                else:
                    a2m_seq.append(c.lower())
        a2m_msa[seq_id] = ''.join(a2m_seq)

    return a2m_msa


def filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=0.6, a2m_format=False, focus_id=None):
    """
    Filter MSA by removing gappy sequences and columns.

    Args:
        msa: dict {seq_id: sequence}
        max_seq_gaps: Maximum fraction of gaps allowed per sequence
        max_col_gaps: Maximum fraction of gaps allowed per column
        a2m_format: When True, gap fractions are computed over match-state
            columns only (uppercase + '-'). Insert columns ('.' and lowercase)
            are structural in A2M and excluded from gap accounting. Column
            filtering also operates on match-state columns only.
        focus_id: When provided, columns where the focus sequence has a residue
            (non-gap match state) are always retained regardless of gap fraction.

    Returns:
        dict: Filtered MSA
    """
    if not msa:
        return {}

    if a2m_format:
        # Identify match-state column indices (uppercase or '-' in any sequence)
        sample = next(iter(msa.values()))
        match_cols = [i for i, c in enumerate(sample) if c == '-' or c.isupper()]

        # Remove sequences with too many deletions in match-state columns
        filtered_seqs = {}
        for seq_id, seq in msa.items():
            match_states = [seq[i] for i in match_cols]
            gap_count = sum(1 for c in match_states if c == '-')
            gap_frac = gap_count / len(match_states) if match_states else 1.0
            if gap_frac <= max_seq_gaps:
                filtered_seqs[seq_id] = seq

        if not filtered_seqs:
            return {}

        focus_seq = filtered_seqs.get(focus_id) if focus_id else None

        # Remove match-state columns with too many gaps; always retain focus residue columns
        seq_list = list(filtered_seqs.values())
        n_seqs = len(seq_list)
        cols_to_keep = []
        for i in match_cols:
            if focus_seq is not None and focus_seq[i] != '-':
                cols_to_keep.append(i)
                continue
            col = [s[i] for s in seq_list]
            gap_frac = sum(1 for c in col if c == '-') / n_seqs
            if gap_frac <= max_col_gaps:
                cols_to_keep.append(i)

        # Reconstruct full sequences keeping only retained match columns;
        # insert columns are dropped since downstream tools use match states only
        final_msa = {}
        for seq_id, seq in filtered_seqs.items():
            final_msa[seq_id] = ''.join(seq[i] for i in cols_to_keep)

        return final_msa

    # Stockholm / non-A2M path
    filtered_seqs = {}
    for seq_id, seq in msa.items():
        gap_count = seq.count('-') + seq.count('.') + seq.count('!')
        gap_frac = gap_count / len(seq) if len(seq) > 0 else 1.0
        if gap_frac <= max_seq_gaps:
            filtered_seqs[seq_id] = seq

    if not filtered_seqs:
        return {}

    seq_list = list(filtered_seqs.values())
    seq_len = len(seq_list[0])
    n_seqs = len(seq_list)

    cols_to_keep = []
    for i in range(seq_len):
        col = [s[i] for s in seq_list]
        gap_count = sum(1 for c in col if c in '-.')
        gap_frac = gap_count / n_seqs
        if gap_frac <= max_col_gaps:
            cols_to_keep.append(i)

    final_msa = {}
    for seq_id, seq in filtered_seqs.items():
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        final_msa[seq_id] = new_seq

    return final_msa


def _chunk_codons(seq):
    """Split a nucleotide sequence into codon triplets.

    Any triplet containing a gap character ('-' or '.') is normalized to '---'
    so the gap_tokens filter in compute_sequence_weights treats partial-gap
    codons as gaps rather than as informative residues. Without this, a column
    drop in filter_msa_by_gaps that breaks codon alignment leaks partial-gap
    codons (e.g. 'A-G') into the identity calculation, inflating both the
    match and aligned counts and skewing N_eff.

    Trailing 1-2 nt that don't complete a codon are dropped (codon MSAs are
    expected to have lengths that are multiples of 3).
    """
    codons = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        c = seq[i:i+3]
        codons.append('---' if any(ch in '-.' for ch in c) else c)
    return codons


def compute_sequence_weights(msa, identity_threshold=0.8, codon_mode=False):
    """
    Compute sequence weights based on clustering at identity threshold.

    Used for N_eff calculation. Weight = 1 / number of neighbors.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold (0.8 = 80% identity)
        codon_mode: If True, compare codon triplets instead of single characters.
                    Use for codon MSAs where each unit is a 3-nt codon.

    Returns:
        dict: {seq_id: weight}
    """
    seq_ids = list(msa.keys())
    n_seqs = len(seq_ids)

    if n_seqs == 0:
        return {}

    if codon_mode:
        seqs = [_chunk_codons(msa[sid]) for sid in seq_ids]
        gap_tokens = {'---', '...'}
    else:
        seqs = [msa[sid] for sid in seq_ids]
        gap_tokens = {'-', '.'}

    neighbor_counts = [1] * n_seqs

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            matches = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a == b and a not in gap_tokens)
            aligned = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a not in gap_tokens and b not in gap_tokens)
            if aligned > 0:
                identity = matches / aligned
                if identity >= identity_threshold:
                    neighbor_counts[i] += 1
                    neighbor_counts[j] += 1

    weights = {seq_ids[i]: 1.0 / neighbor_counts[i] for i in range(n_seqs)}
    return weights


def compute_neff(msa, identity_threshold=0.8, codon_mode=False):
    """
    Compute effective number of sequences (N_eff).

    N_eff = sum of sequence weights, where weight = 1/n_neighbors.
    Higher N_eff indicates more diverse MSA with better evolutionary signal.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold
        codon_mode: If True, compare codon triplets instead of single characters.

    Returns:
        float: N_eff value
    """
    weights = compute_sequence_weights(msa, identity_threshold, codon_mode=codon_mode)
    return sum(weights.values())



def write_tsv(rows, path, fieldnames=None, *, extrasaction='ignore', mkdir=True):
    """Write a list of row dicts to a tab-separated file.

    Args:
        rows: list of dicts to write.
        path: output file path (str or Path).
        fieldnames: column names. If None, derived from rows[0].keys().
            When rows is empty and fieldnames is provided, a header-only
            file is written; when both are empty, no file is created.
        extrasaction: passed to csv.DictWriter ('ignore' or 'raise').
        mkdir: create parent directories if they don't exist.
    """
    path = Path(path)
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        return
    if fieldnames is None:
        if not rows:
            return
        fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction=extrasaction)
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def mutation_class(wt_aa, mut_aa):
    """Classify a mutation as SYNONYMOUS, MISSENSE, STOP_GAIN, STOP_LOSS, or UNKNOWN."""
    if wt_aa == "X" or mut_aa == "X":
        return "UNKNOWN"
    if wt_aa == mut_aa:
        return "SYNONYMOUS"
    if mut_aa == "Stop" and wt_aa != "Stop":
        return "STOP_GAIN"
    if wt_aa == "Stop" and mut_aa != "Stop":
        return "STOP_LOSS"
    return "MISSENSE"


def write_a2m(msa, output_path, focus_seq_id=None):
    """
    Write MSA to A2M format file.

    Args:
        msa: dict {seq_id: sequence}
        output_path: Output file path
        focus_seq_id: If provided, write focus sequence first
    """
    with open(output_path, 'w') as f:
        if focus_seq_id and focus_seq_id in msa:
            f.write(f">{focus_seq_id}\n{msa[focus_seq_id]}\n")
        for seq_id, seq in msa.items():
            if seq_id != focus_seq_id:
                f.write(f">{seq_id}\n{seq}\n")
