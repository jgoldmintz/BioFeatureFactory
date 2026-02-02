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
import inspect
import warnings
import re
import os
import math
import tempfile
import subprocess
import shutil
import requests
import time
import traceback
import sys
import json
import logging
from collections import Counter
from pathlib import Path
from urllib.parse import unquote
from Bio import Entrez
from typing import Dict, List, Tuple

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
                    name = re.search(">[a-zA-Z]+_?\d+(\.\d+)*", line).group(0)
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
                        matches = re.findall("/_\d+$/", name)
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
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1]) - 1  # Convert to 0-based index
    return position, (original_nt, mutant_nt)

def get_mutation_data_bioAccurate(ntposnt):
    """Return one-based position and nucleotides for a mutation string; skips stop codons."""

    # Skip stop codons (for the case of aa)
    if 'Stop' in ntposnt or 'Sto' in ntposnt:
        return None, None

    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1])
    return position, (original_nt, mutant_nt)

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

def convert_position(seq1, seq2, position1, space="-"):
    """Project a one-based index from seq1 onto seq2 while accounting for gap characters."""
    error = None

    if position1 == 0:
        warnings.warn("\033[93m" + "Position given is not 1-indexed in function " + str(inspect.stack()[1].function) + "\033[0m")
        error = "Position given is not 1-indexed\n" + str(inspect.stack()[1].function)

    i1 = 0
    i2 = 0
    increment = 0
    while i1 < int(position1) and increment < len(seq1) and increment < len(seq2):
        if seq1[increment] != space:
            i1 += 1
        if seq2[increment] != space:
            i2 += 1
        increment += 1
    if not seq1[increment - 1] == seq1.replace(space, "")[i1 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq1[increment - 1] + " " + seq1.replace(space, "")[i1 - 1] + "\033[0m")
    elif seq2[increment - 1] != space and not seq2[increment - 1] == seq2.replace(space, "")[i2 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq2[increment - 1] + " " + seq2.replace(space, "")[i2 - 1] + "\033[0m")
    if seq2[increment - 1] == space:
        error = "Sequence 1 position aligns with a gap in sequence 2\n" + str(inspect.stack()[1].function)
    return (i2, error)

# retries a requests library function
# assumes a response object is returned (which can be raised for HTML status)
def retry_request(func, positional_arguments=[], keyword_arguments={}, lim=10, wait=2):
    """Retry a requests-style call, returning the response or None after repeated failures."""
    for i in range(lim):
        try:
            response = func(*positional_arguments, **keyword_arguments)
            response.raise_for_status()
            return (response)
        except:
            time.sleep(wait * i)
    warnings.warn("\033[93m" + str(func.__name__) + " failed after " + str(lim) + " tries." + "\033[0m")
    return (None)

def retry_func(func, positional_arguments=[], keyword_arguments={}, lim=10, wait=2):
    """Retry a callable that may raise exceptions, returning its result or None."""
    for i in range(lim):
        try:
            return (func(*positional_arguments, **keyword_arguments))
        except Exception as e:
            time.sleep(wait * i)
    warnings.warn("\033[93m" + str(func.__name__) + " failed after " + str(lim) + " tries." + "\033[0m")
    return (None)

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
        match = [c for c in candidates if c[1] == transcript_id]
        if not match:
            # Try matching without version suffix (e.g., NM_022162 matches NM_022162.3)
            tid_base = transcript_id.rsplit('.', 1)[0]
            match = [c for c in candidates if c[1].rsplit('.', 1)[0] == tid_base]
        if match:
            _, best_tid, best = match[0]
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

def detect_reference_build_from_gid(gid):
    """Return reference build label for a RefSeq genomic accession."""
    if not gid:
        return None
    for build, mapping in chromosome_map.items():
        if gid in mapping.values():
            return build
    return None


def collect_gene_reference_builds(mutation_files, annotation_file, assembly="GRCh38"):
    """Infer reference builds for mutation files using local annotation data."""
    if not annotation_file:
        raise ValueError("annotation_file is required to collect reference builds")

    gene_builds, failed_genes, build_counts = {}, {}, Counter()
    for mutation_file in mutation_files:
        try:
            gene_name = extract_gene_from_filename(str(mutation_file))
            print(f"Detecting reference build for {gene_name}...")

            genome_loc_result = get_genome_loc(gene_name, annotation_file, assembly=assembly)
            if not genome_loc_result:
                failed_genes[gene_name] = "Could not retrieve genome location"
                continue

            chromosome = genome_loc_result["chrom"]
            coordinates = [genome_loc_result["tx_start"], genome_loc_result["tx_end"]]
            strand = genome_loc_result["strand"]
            gid = chromosome_map.get(assembly, {}).get(str(chromosome))
            build = detect_reference_build_from_gid(gid)
            if not build:
                failed_genes[gene_name] = f"Unknown reference build for gid: {gid}"
                continue

            gene_builds[gene_name] = {
                "build": build,
                "gid": gid,
                "chromosome": chromosome,
                "coordinates": coordinates,
                "strand": strand,
            }
            build_counts[build] += 1
        except Exception as e:
            failed_genes[gene_name] = str(e)
    return build_counts, gene_builds, failed_genes


def determine_consensus_build(build_counts):
    """Return the most frequent reference build observed."""
    return build_counts.most_common(1)[0][0] if build_counts else None


def get_reference_cache_path(cache_dir="./reference_cache"):
    """Return path to reference metadata cache file, ensuring directory exists."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(exist_ok=True)
    return cache_path / "reference_metadata.json"


def load_reference_cache(cache_dir="./reference_cache"):
    """Load reference metadata cache from disk if present."""
    cache_file = get_reference_cache_path(cache_dir)
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not load reference cache: {e}", file=sys.stderr)
    return {}


def save_reference_cache(cache_data, cache_dir="./reference_cache"):
    """Persist reference metadata cache to disk."""
    cache_file = get_reference_cache_path(cache_dir)
    try:
        with open(cache_file, "w") as f:
            json.dump(cache_data, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not save reference cache: {e}", file=sys.stderr)


def download_reference_genome(build, cache_dir="./reference_cache"):
    """Download and decompress a reference genome FASTA for the requested build."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(exist_ok=True)
    urls = {
        "GRCh38": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz",
        "GRCh37": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz",
    }
    if build not in urls:
        print(f"Error: Unknown reference build: {build}", file=sys.stderr)
        return None

    url = urls[build]
    output_file = cache_path / f"{build}_reference.fna.gz"
    uncompressed_file = cache_path / f"{build}_reference.fna"

    print(f"Downloading {build} reference genome...")

    if shutil.which("aria2c"):
        cmd = [
            "aria2c",
            "--max-connection-per-server=8",
            "--split=8",
            "--min-split-size=1M",
            "--continue=true",
            "--dir",
            str(cache_path),
            "--out",
            output_file.name,
            url,
        ]
    elif shutil.which("wget"):
        cmd = ["wget", "-O", str(output_file), url]
    elif shutil.which("curl"):
        cmd = ["curl", "-L", "-o", str(output_file), url]
    else:
        print("Error: No download tool available (aria2c, wget, or curl required)", file=sys.stderr)
        return None

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error downloading reference: {result.stderr}", file=sys.stderr)
            return None

        print(f"Decompressing {build} reference genome...")
        result = subprocess.run(["gunzip", str(output_file)], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error decompressing reference: {result.stderr}", file=sys.stderr)
            return None

        print(f"Successfully downloaded {build} reference genome")
        return str(uncompressed_file)
    except Exception as e:
        print(f"Error downloading reference genome: {e}", file=sys.stderr)
        return None

def load_mapping(mapping_file: str, mapType: str ='transcript') -> Dict[str, str]:
    """Load a two-column mapping CSV (mutant->mapping) using the specified column name."""

    mapping = {}
    tval, delim = validate_mapping_content(mapping_file)
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
    position_data = get_mutation_data_bioAccurate(aaposaa)
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
    pat = re.compile(r'^[A-Za-z][A-Za-z0-9\-]{1,14}$')
    return 2 <= len(name) <= 15 and bool(pat.match(name))

def extract_gene_from_filename(filename: str) -> str:
    """Return the most likely gene symbol from a filename."""
    name = Path(filename).name
    name = strip_all_extensions(name)

    # Fast path: choose the last gene-like token bounded by underscores
    matches = re.findall(r'(?:^|_)([A-Z][A-Z0-9]{1,14}(?:-[A-Z0-9]+)?)(?=$|_)', name)
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
    m = re.search(r'([A-Za-z][A-Za-z0-9\-]{1,14})', name)
    if m and is_likely_gene_name(m.group(1)):
        return m.group(1)

    # Fallback: basename without extensions
    return strip_all_extensions(Path(filename).stem)

def discover_mapping_files(mapping_dir):
    """Scan directory for CSV mapping files and extract gene names flexibly

    Args:
        mapping_dir: Directory path containing mapping CSV files

    Returns:
        dict: {gene_name: file_path} mapping
    """
    from pathlib import Path

    mapping_files = {}

    if not mapping_dir or not Path(mapping_dir).exists():
        return mapping_files

    # Scan for all CSV files
    for csv_file in Path(mapping_dir).glob("*.csv"):
        try:
            # Extract gene name from filename
            gene_name = extract_gene_from_filename(csv_file.stem)

            # Validate CSV content structure
            if validate_mapping_content(csv_file)[0]:
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
        for fasta_file in Path(fasta_dir).glob(extension):
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
            has_aa_mutation = any(col in fieldnames for col in
                                  ['aamutant', 'transcript', 'genomic', 'aa_mutation', 'amino_acid_mutation',
                                   'protein_mutation', 'chromosome'])

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
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"WT fasta directory not found: {input_dir}")

    sequences: Dict[str, str] = {}
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
            if non_empty:
                seq = max(non_empty, key=lambda kv: len(kv[1]))[1].strip()

        if not seq:
            continue

        gene = extract_gene_from_filename(fasta_file.name) or fasta_file.stem
        sequences[gene.upper()] = seq

    header_label = wt_header.upper() if wt_header else "SEQUENCES"
    print(f"[WT] Loaded {len(sequences)} WT {header_label} into memory")
    return sequences


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
    pos, mut = get_mutation_data_bioAccurate(ntposnt)
    if pos is None:
        return None, "", "", 0, 0, 0

    # Convert to 0-based indexing
    pos_0_indexed = pos - 1
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

def run_jackhmmer(query_fasta, database, output_sto, jackhmmer_binary,
                  iterations=5, evalue_inclusion=1e-3, threads=4):
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
        dict: {seq_id: sequence} mapping
    """
    from collections import defaultdict
    current_seqs = defaultdict(str)

    with open(stockholm_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                current_seqs[seq_id] += seq

    return dict(current_seqs)


def stockholm_to_a2m(msa, focus_seq_id):
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

    focus_seq = msa[focus_seq_id]
    match_columns = [i for i, c in enumerate(focus_seq) if c not in '-.']

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


def filter_msa_by_gaps(msa, max_seq_gaps=0.4, max_col_gaps=0.6):
    """
    Filter MSA by removing gappy sequences and columns.

    Args:
        msa: dict {seq_id: sequence}
        max_seq_gaps: Maximum fraction of gaps allowed per sequence
        max_col_gaps: Maximum fraction of gaps allowed per column

    Returns:
        dict: Filtered MSA
    """
    if not msa:
        return {}

    # Remove sequences with too many gaps
    filtered_seqs = {}
    for seq_id, seq in msa.items():
        gap_count = seq.count('-') + seq.count('.') + seq.count('!')
        gap_frac = gap_count / len(seq) if len(seq) > 0 else 1.0
        if gap_frac <= max_seq_gaps:
            filtered_seqs[seq_id] = seq

    if not filtered_seqs:
        return {}

    # Identify columns with too many gaps
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

    # Remove gappy columns
    final_msa = {}
    for seq_id, seq in filtered_seqs.items():
        new_seq = ''.join(seq[i] for i in cols_to_keep)
        final_msa[seq_id] = new_seq

    return final_msa


def compute_sequence_weights(msa, identity_threshold=0.8):
    """
    Compute sequence weights based on clustering at identity threshold.

    Used for N_eff calculation. Weight = 1 / number of neighbors.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold (0.8 = 80% identity)

    Returns:
        dict: {seq_id: weight}
    """
    seq_ids = list(msa.keys())
    n_seqs = len(seq_ids)

    if n_seqs == 0:
        return {}

    seqs = [msa[sid] for sid in seq_ids]
    neighbor_counts = [1] * n_seqs

    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            matches = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a == b and a not in '-.')
            aligned = sum(1 for a, b in zip(seqs[i], seqs[j])
                         if a not in '-.' and b not in '-.')
            if aligned > 0:
                identity = matches / aligned
                if identity >= identity_threshold:
                    neighbor_counts[i] += 1
                    neighbor_counts[j] += 1

    weights = {seq_ids[i]: 1.0 / neighbor_counts[i] for i in range(n_seqs)}
    return weights


def compute_neff(msa, identity_threshold=0.8):
    """
    Compute effective number of sequences (N_eff).

    N_eff = sum of sequence weights, where weight = 1/n_neighbors.
    Higher N_eff indicates more diverse MSA with better evolutionary signal.

    Args:
        msa: dict {seq_id: sequence}
        identity_threshold: Clustering threshold

    Returns:
        float: N_eff value
    """
    weights = compute_sequence_weights(msa, identity_threshold)
    return sum(weights.values())


def validate_msa_quality(msa, min_neff_ratio=10, query_length=None, focus_seq_id=None):
    """
    Validate MSA quality for EVmutation analysis.

    Args:
        msa: dict {seq_id: sequence}
        min_neff_ratio: Minimum N_eff / L ratio (default: 10)
        query_length: Length of query sequence (computed if not provided)
        focus_seq_id: ID of focus sequence (for length calculation)

    Returns:
        dict: Quality metrics including pass/fail status
    """
    if not msa:
        return {'pass': False, 'error': 'Empty MSA'}

    # Get query length
    if query_length is None:
        if focus_seq_id and focus_seq_id in msa:
            query_length = len(msa[focus_seq_id].replace('-', '').replace('.', ''))
        else:
            first_seq = next(iter(msa.values()))
            query_length = len(first_seq.replace('-', '').replace('.', ''))

    # Compute N_eff
    neff = compute_neff(msa)
    neff_ratio = neff / query_length if query_length > 0 else 0
    min_neff = min_neff_ratio * query_length

    return {
        'pass': neff >= min_neff,
        'n_sequences': len(msa),
        'n_eff': round(neff, 1),
        'n_eff_ratio': round(neff_ratio, 2),
        'query_length': query_length,
        'min_neff_required': min_neff
    }


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
