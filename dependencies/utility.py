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
from collections import Counter
from pathlib import Path
from urllib.parse import unquote
from Bio import Entrez
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


def read_fasta(inf, aformat="FIRST", duplicate="replace"):
    """Load sequences from a FASTA file into a name→sequence dictionary."""
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


def _prepare_structured_annotation(genename, annotation_file, assembly, fmt):
    """Extract gene and transcript features from RefSeq/Ensembl-style GTF or GFF3 files."""
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

                transcript_id = None

                if feature_lc in transcript_features:
                    transcript_id = attrs.get("transcript_id") or attrs.get("ID")
                elif feature_lc == "cds" and not attrs.get("transcript_id") and attrs.get("Parent"):
                    possible = [val for val in parent_values if not val.startswith("gene-")]
                    if possible:
                        transcript_id = possible[0]

                if feature_lc == "gene" and matches_symbol and not target_gene_ids:
                    target_gene_ids.update(attr_gene_ids)

                if feature_lc in transcript_features and (matches or bool(set(parent_values).intersection(target_gene_ids)) or (attrs.get("gene_id") and attrs.get("gene_id") in target_gene_ids)):
                    if not transcript_id:
                        continue
                    rec = target_transcripts.setdefault(transcript_id, {
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


def get_genome_loc(genename, annotation_file, assembly="GRCh38"):
    """Return gene coordinates and exon structure from supported annotation formats."""
    fmt = _detect_annotation_format(annotation_file)

    if fmt == "custom":
        try:
            return _prepare_custom_annotation(genename, annotation_file)
        except Exception as e:
            print(f"Error parsing annotation file {annotation_file}: {e}", file=sys.stderr)
            return None

    try:
        return _prepare_structured_annotation(genename, annotation_file, assembly, fmt)
    except Exception as e:
        print(f"Error parsing annotation file {annotation_file}: {e}", file=sys.stderr)
        return None


def get_genome_loc3(genename, annotation_file, assembly="GRCh38"):
    """
    Get gene location using local annotation file.

    Args:
        genename (str): Gene symbol (e.g., "CD44")
        annotation_file (str): Path to local annotation file
        assembly (str): Reference assembly ("GRCh37" or "GRCh38")

    Returns:
        tuple: (chromosome, [start, end], strand) or None
    """
    try:
        with open(annotation_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 7:
                    continue

                # Format: NAME, CHROM, STRAND, TX_START, TX_END, EXON_START, EXON_END
                gene_symbol = fields[0]
                chromosome = fields[1].replace('chr', '')  # Remove 'chr' prefix if present
                strand = fields[2]
                tx_start = int(fields[3])
                tx_end = int(fields[4])

                if gene_symbol.upper() == genename.upper():
                    return (chromosome, [tx_start, tx_end], strand)

    except Exception as e:
        print(f"Error parsing annotation file {annotation_file}: {e}", file=sys.stderr)

    return None

def get_genome_loc2(genename, assembly="GRCh38", organism="homo sapiens"):
    """Query NCBI Gene for genomic location data when local annotations are unavailable."""
    organism_formatted = re.sub(r'\s+', '+', organism, flags=re.MULTILINE)
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={genename}[gene]+AND+{organism_formatted}[organism]"
    r = retry_request(requests.get, [search_url])

    if not r:
        return (None)

    ids = re.findall("<Id>(\d+)</Id>", r.text)
    for accid in ids:
        try:
            def fetch_gene_record():
                return Entrez.efetch(db="gene", id=accid, retmode="xml")

            handle = retry_func(fetch_gene_record, [], {})
            if not handle:
                continue

            record = str(handle.read())

            # Case-insensitive gene name matching
            if genename.lower() not in record.lower():
                continue

            # NEW APPROACH: Extract chromosome from genomic Gene-commentary sections
            gene_commentaries = record.split("<Gene-commentary>")

            # Look for genomic commentary sections that might contain coordinate information
            for section in gene_commentaries:
                # Skip if this section doesn't look like a genomic section
                if 'type value="genomic"' not in section and 'genomic' not in section.lower():
                    continue

                # Try to extract chromosome from this specific genomic section
                chromosome = None

                # Method 1: Look for chromosome information with prioritized patterns
                # More specific patterns first to avoid false matches
                chromosome_patterns = [
                    r"Chromosome (\d+) Reference",  # Most specific: "Chromosome 17 Reference"
                    r"Chromosome ([1-9]|1[0-9]|2[0-2]|X|Y) Reference",  # "Chromosome X Reference"
                    r"Chromosome ([1-9]|1[0-9]|2[0-2]|X|Y) Alternate",  # "Chromosome 17 Alternate"
                    r"Chromosome ([1-9]|1[0-9]|2[0-2]|X|Y) Primary",     # "Chromosome 17 Primary"
                    r">Chromosome ([1-9]|1[0-9]|2[0-2]|X|Y)<",           # XML tags
                    r"Chr ([1-9]|1[0-9]|2[0-2]|X|Y) ",                   # "Chr 17 "
                    r"Chromosome ([1-9]|1[0-9]|2[0-2]|X|Y)",             # Generic fallback
                ]

                for pattern in chromosome_patterns:
                    matches = re.findall(pattern, section, re.IGNORECASE)
                    if matches:
                        chromosome = matches[0]
                        break

                if not chromosome:
                    continue

                # Convert word-based chromosome names to numbers if needed
                if chromosome.lower() in w2n:
                    chromosome = w2n[chromosome.lower()]

                # Check if this chromosome exists in our mapping
                if str(chromosome) not in chromosome_map[assembly]:
                    continue

                gid = chromosome_map[assembly][str(chromosome)]

                # Look for sections containing this RefSeq accession
                refseq_sections = [x for x in gene_commentaries if gid.split(".")[0] in x]
                version_sections = [x for x in refseq_sections if re.search(
                    rf"<Gene-commentary_version>{gid.split('.')[-1]}</Gene-commentary_version>", x)]

                if len(version_sections) > 0:
                    s = version_sections[0]

                    # Extract coordinates
                    from_matches = re.findall(r"<Seq-interval_from>(\d+)</Seq-interval_from>", s)
                    to_matches = re.findall(r"<Seq-interval_to>(\d+)</Seq-interval_to>", s)

                    if not from_matches or not to_matches:
                        continue

                    r1 = int(max(from_matches, key=len))
                    r2 = int(max(to_matches, key=len))

                    # Extract strand information
                    strand = "unknown"

                    # Look for Na-strand value format
                    na_strand_matches = re.findall(r'<Na-strand value="([^"]+)"/>', s)
                    if na_strand_matches:
                        na_strand = na_strand_matches[0].lower()
                        if na_strand == "plus":
                            strand = "+"
                        elif na_strand == "minus":
                            strand = "-"
                        else:
                            strand = na_strand
                    else:
                        # Fallback: Look for numeric strand format
                        strand_matches = re.findall(r"<Seq-interval_strand>(\d+)</Seq-interval_strand>", s)
                        if strand_matches:
                            strand_num = int(strand_matches[0])
                            if strand_num == 1:
                                strand = "+"
                            elif strand_num == 2:
                                strand = "-"
                            else:
                                strand = f"unknown({strand_num})"

                    return chromosome, [r1, r2], strand, gid

        except KeyError as e:
            traceback.print_exc()
        except Exception as e:
            print(traceback.format_exc())
            print(e)
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

def load_mapping(mapping_file, mapType='transcript'):
    """Load a two-column mapping CSV (mutant→mapping) using the specified column name."""

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

def run_docker_command(docker_image, fasta_file, command_template, output_file, timeout=300):
    """Generic Docker execution wrapper

    Args:
        docker_image: Docker image name
        fasta_file: Input FASTA file path
        command_template: Command template (e.g., "./ape -m netphos {input}")
        output_file: Output file path
        timeout: Command timeout in seconds

    Returns:
        tuple: (success, output_content, error_message)
    """
    work_dir = tempfile.mkdtemp()

    try:
        # Copy the FASTA file to the work directory
        docker_input = os.path.join(work_dir, "input.fasta")
        shutil.copy2(fasta_file, docker_input)

        # Build Docker command
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{work_dir}:/data",
            docker_image
        ]

        # Add command from template
        command = command_template.format(input="/data/input.fasta")
        docker_cmd.extend(command.split())

        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )

        if result.returncode == 0:
            # Save output to file
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            return True, result.stdout, None
        else:
            return False, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        return False, "", f"Docker command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)
    finally:
        # Clean up work directory
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir, ignore_errors=True)

def extract_mutation_from_sequence_name(seq_name):
    """Extract mutation ID from sequence name (e.g., 'ZFP36-C330T' → 'C330T')

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

def process_single_mutation_for_sequence(seq_name, predictions, mapping_dict, is_mutant=True, tool_type='netphos'):
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

def parse_predictions_with_mutation_filtering(predictions, mapping_dict, is_mutant, threshold=0.0, yes_only=False, tool_type='netphos'):
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
                seq_name, seq_preds, mapping_dict, is_mutant=True, tool_type=tool_type
            )

            # Apply additional filters
            for result in seq_results:
                # Apply threshold filter
                score_field = 'score' if tool_type == 'netphos' else 'potential'
                if result[score_field] < threshold:
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

def extract_gene_from_filename(filename):
    """Extract gene name from filename using intelligent pattern matching

    Args:
        filename: Filename with or without extension (e.g., 'BRCA1_aa.csv', 'combined_ACADM.csv', 'CYBB_transcript_mutations.csv')

    Returns:
        str: Extracted gene name (e.g., 'BRCA1', 'ACADM', 'CYBB')
    """
    import re
    from pathlib import Path

    # Remove file extension if present
    if '.' in filename:
        filename = Path(filename).stem

    # Method 1: Try to directly extract gene name patterns
    # Look for the pattern at the beginning of the string (before the first underscore with a suffix)
    gene_pattern = re.compile(r'^([A-Z][A-Z0-9]{1,14}(?:[-][A-Z0-9]+)?)(?:_|$)')
    match = gene_pattern.match(filename)

    if match and is_likely_gene_name(match.group(1)):
        return match.group(1)

    # Method 2: Strip known prefixes and suffixes
    gene_name = filename

    # Common prefixes to remove
    prefixes = [
        r'^combined_',
        r'^merged_',
        r'^processed_',
        r'^final_',
        r'^updated_',
        r'^new_',
        r'^old_',
        r'^temp_',
        r'^test_',
    ]

    # Common suffixes to remove (order matters - more specific first)
    suffixes = [
        r'_transcript_mutations$',  # Remove '_transcript_mutations' suffix
        r'_protein_mutations$',  # Remove '_protein_mutations' suffix
        r'_nt_to_aa_mapping$',  # Remove '_nt_to_aa_mapping' suffix
        r'_mutations$',  # Remove '_mutations' suffix
        r'_transcript$',  # Remove '_transcript' suffix  <-- ADD THIS
        r'_protein$',  # Remove '_protein' suffix
        r'_mapping$',  # Remove '_mapping' suffix
        r'_sequences$',  # Remove '_sequences' suffix
        r'_variants$',  # Remove '_variants' suffix
        r'_data$',  # Remove '_data' suffix
        r'_aa$',  # Remove '_aa' suffix
        r'_nt$',  # Remove '_nt' suffix
        r'_cds$',  # Remove '_cds' suffix
        r'_orf$',  # Remove '_orf' suffix
    ]

    # Remove prefixes
    for pattern in prefixes:
        gene_name = re.sub(pattern, '', gene_name, flags=re.IGNORECASE)

    # Remove suffixes
    for pattern in suffixes:
        gene_name = re.sub(pattern, '', gene_name, flags=re.IGNORECASE)

    # Clean up underscores
    gene_name = gene_name.strip('_')

    # Validate and return
    if gene_name and is_likely_gene_name(gene_name):
        return gene_name

    # Fallback to cleaned name or original if nothing worked
    return gene_name if gene_name else Path(filename).stem

def is_likely_gene_name(name):
    """Check if a string looks like a gene name

    Args:
        name: String to check

    Returns:
        bool: True if likely a gene name
    """
    import re

    # Gene names are typically:
    # - 2-15 characters long
    # - Start with a letter
    # - Contain only letters, numbers, and occasionally hyphens

    gene_pattern = re.compile(r'^[A-Za-z][A-Za-z0-9\-]{1,14}$')

    # Additional check: common gene name patterns
    if 2 <= len(name) <= 15 and gene_pattern.match(name):
        return True

    return False

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
