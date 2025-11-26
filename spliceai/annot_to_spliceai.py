#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023–2025  Jacob Goldmintz
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

"""Convert a GTF file into SpliceAI-compatible annotation format."""

import argparse
from collections import defaultdict
from pathlib import Path

HEADER = "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n"

SIMPLE_TO_NCBI = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

NCBI_TO_SIMPLE = {value.split(".")[0]: key for key, value in SIMPLE_TO_NCBI.items()}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert GTF transcript annotations to SpliceAI tab format",
    )
    parser.add_argument("gtf", help="Input GTF file (e.g. converted from RefSeq GFF3)")
    parser.add_argument(
        "-o",
        "--out",
        help="Output annotation file (default: <input>.spliceai.txt)",
    )
    parser.add_argument(
        "--chromosome-format",
        choices=["simple", "ncbi", "ucsc"],
        default="simple",
        help="How to label chromosomes in the output (default: simple numbers)",
    )
    return parser.parse_args()


def _format_chromosome(raw: str, fmt: str) -> str | None:
    chrom = raw.strip()
    stripped = chrom[3:] if chrom.lower().startswith("chr") else chrom
    base = stripped.split(".")[0]

    if fmt == "simple":
        if stripped.startswith("NC_"):
            return NCBI_TO_SIMPLE.get(base, stripped)
        return stripped

    if fmt == "ncbi":
        if stripped.startswith("NC_"):
            return stripped
        simple = stripped.upper()
        return SIMPLE_TO_NCBI.get(simple, stripped)

    if fmt == "ucsc":
        simple = _format_chromosome(raw, "simple")
        if simple is None:
            return None
        return simple if simple.lower().startswith("chr") else f"chr{simple}"

    return stripped


def collect_transcripts(gtf_path: Path, chrom_format: str):
    transcripts = defaultdict(lambda: {"chrom": None, "gene": None, "strand": None, "exons": []})

    with gtf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = parts
            if feature.lower() != "exon":
                continue

            chrom_name = _format_chromosome(chrom, chrom_format)
            if chrom_name is None:
                continue

            attr_pairs = {}
            for field in attrs.split(";"):
                field = field.strip()
                if not field or " " not in field:
                    continue
                key, value = field.split(" ", 1)
                attr_pairs[key.strip()] = value.strip().strip('"')

            transcript_id = attr_pairs.get("transcript_id")
            if not transcript_id:
                continue

            gene_name = (
                attr_pairs.get("gene_name")
                or attr_pairs.get("gene")
                or attr_pairs.get("gene_id")
                or transcript_id
            )
            if gene_name.startswith('#NAME'):
                continue

            data = transcripts[transcript_id]
            data["chrom"] = chrom_name
            data["gene"] = gene_name
            data["strand"] = strand
            data["exons"].append((int(start), int(end)))

    return transcripts


def write_spliceai(transcripts, out_path: Path):
    with out_path.open("w") as out:
        out.write(HEADER)
        for tx, info in transcripts.items():
            exons = info["exons"]
            if not exons:
                continue

            sorted_exons = sorted(exons, key=lambda e: e[0])
            starts = ",".join(str(s) for s, _ in sorted_exons) + ","
            ends = ",".join(str(e) for _, e in sorted_exons) + ","
            tx_start = sorted_exons[0][0]
            tx_end = sorted_exons[-1][1]

            out.write(
                f"{info['gene']}\t{info['chrom']}\t{info['strand']}\t"
                f"{tx_start}\t{tx_end}\t{starts}\t{ends}\n"
            )


def main():
    args = parse_args()
    gtf_path = Path(args.gtf).resolve()
    if not gtf_path.exists():
        raise SystemExit(f"Input GTF '{gtf_path}' does not exist")

    transcripts = collect_transcripts(gtf_path, args.chromosome_format)

    out_path = (
        Path(args.out)
        if args.out
        else gtf_path.with_suffix("").with_suffix(".spliceai.txt")
    )
    write_spliceai(transcripts, out_path.resolve())
    print(f"Wrote SpliceAI annotation: {out_path}")


if __name__ == "__main__":
    main()
