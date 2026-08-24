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
BioFeatureFactory: Annotation

Genome-annotation loading: GTF/GFF3 and the custom tab-delimited format,
transcript selection, and chromosome-name normalisation. Consumers:
variant_mapping and vcf_converter.

Split out of utility.py, which had grown to 92 symbols. utility.py re-exports
every name here lazily, so existing `from ...utility import X` callers are
unaffected.
"""

import csv
import os
import re
import sys
import math
import shutil
import tempfile
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Dict, Optional, Tuple
from urllib.parse import unquote  # GFF3 attribute values are percent-encoded

from Bio.Seq import Seq

from biofeaturefactory.lib.primitives import (
    chromosome_map,
)


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

