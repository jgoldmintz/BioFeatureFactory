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

"""Preflight: can this VCF join to this chromosome mapping?

Runs BEFORE run_spliceai, not after. The join is what spliceai-parser.py does at
the END of the pipeline, so a mapping that does not describe the VCF is currently
discovered only once SpliceAI has already scored every variant -- minutes to hours
of GPU per gene, thrown away. This is the same question asked in seconds, with no
model involved.

The check cannot disagree with the parser, because it calls the parser's own
resolve_pkey: same three tiers, same order (exact string equality, then
canonical_token representation-normalisation, then reference left-alignment).
spliceai-parser.py is loaded by path rather than imported by name because its
filename carries a hyphen.

Exit codes:
    0  every record joined, or enough did (see --min-match-fraction)
    1  the VCF carries no data records, or too few joined

WHEN THIS FIRES IN A NORMAL RUN: it should not. With generate_vcfs in the loop
main.nf hands the SAME chromosome mapping to the converter (main.nf:266) and to
the parser (main.nf:299), so a VCF built by core/vcf_converter.py cannot disagree
with the mapping it was built from. A failure here means one of:
  * --input_vcf_path / --skip_vcf_generation supplied a VCF with independent
    provenance -- built against a different reference build, a different
    chromosome naming convention, or an older version of the mappings;
  * the mappings were regenerated after the VCF was written;
  * the VCF is empty because every token was refused during conversion.
"""

import argparse
import importlib.util
import os
import sys
from pathlib import Path

from biofeaturefactory.lib.utility import load_mapping

_HERE = Path(__file__).resolve().parent


def _load_parser_module():
    """Load spliceai-parser.py by path (its filename is not a valid module name)."""
    path = _HERE / "spliceai-parser.py"
    if not path.exists():
        sys.exit(f"[check_vcf_mapping] ERROR: {path} not found; it ships beside this script.")
    spec = importlib.util.spec_from_file_location("_spliceai_parser", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_vcf_records(vcf_file):
    """(chrom, pos, ref, alt) per data line. ALT is left as written, commas and all."""
    records = []
    malformed = 0
    with open(vcf_file) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 8:
                malformed += 1
                continue
            records.append((fields[0], fields[1], fields[3], fields[4]))
    return records, malformed


def check(vcf_file, mapping_file, reference=None, min_fraction=0.0, gene=None):
    """Return (ok, message). ok is False when the run should not proceed."""
    P = _load_parser_module()

    label = gene or Path(vcf_file).stem

    if not os.path.exists(mapping_file):
        return False, f"chromosome mapping not found: {mapping_file}"

    chromosome_mapping = load_mapping(mapping_file, mapType='chromosome')
    if not chromosome_mapping:
        return False, (f"{label}: chromosome mapping '{mapping_file}' holds no usable entries. "
                       f"Regenerate it with core/variant_mapping.py.")

    records, malformed = read_vcf_records(vcf_file)
    if not records:
        return False, (f"{label}: '{vcf_file}' carries NO data records. SpliceAI would score "
                       f"nothing and every downstream table would be header-only. If this VCF "
                       f"came from core/vcf_converter.py, check its stderr -- every token was "
                       f"likely refused (commonly no_chromosome_mapping).")

    # Same three tiers the parser uses, in the same order, via its own function.
    notation_index = P.build_notation_index(chromosome_mapping)
    refbases = P.ReferenceBases(reference) if reference else None

    matched = 0
    modes = {}
    unmatched_examples = []
    left_align_index = None
    left_align_chrom = None

    for chrom, pos, ref, alt in records:
        if refbases is not None and chrom != left_align_chrom:
            left_align_index = P.build_left_align_index(chromosome_mapping, chrom, refbases)
            left_align_chrom = chrom
        # One allele is enough to call the record joinable; a multi-allelic record
        # is scored per allele by SpliceAI and joined per allele by the parser.
        hit_mode = None
        for allele in alt.split(','):
            # gene_context is a placeholder: this asks whether the genomic notation
            # RESOLVES, not what its pkey is. resolve_pkey returns 'no_gene_context'
            # only after a successful match, so a non-empty string keeps the two
            # failure modes distinguishable.
            _pkey, mode = P.resolve_pkey(
                pos, ref, allele, "PREFLIGHT", chromosome_mapping, {},
                None, notation_index, chrom, refbases, left_align_index,
            )
            if mode in ('exact', 'canonical', 'left_aligned'):
                hit_mode = mode
                break
        if hit_mode:
            matched += 1
            modes[hit_mode] = modes.get(hit_mode, 0) + 1
        elif len(unmatched_examples) < 5:
            unmatched_examples.append(f"{ref}{pos}{alt}")

    total = len(records)
    fraction = matched / total if total else 0.0
    detail = ", ".join(f"{k}={v}" for k, v in sorted(modes.items())) or "none"
    summary = (f"{label}: {matched}/{total} VCF records join the chromosome mapping "
               f"({fraction:.1%}; {detail})")
    if malformed:
        summary += f"; {malformed} line(s) refused for having <8 columns"

    if matched == 0:
        return False, (
            f"{summary}.\n"
            f"       The mapping does not describe this VCF. Nothing would be written after "
            f"SpliceAI runs.\n"
            f"       unmatched e.g.: {', '.join(unmatched_examples)}\n"
            f"       mapping e.g.  : "
            f"{', '.join(list(chromosome_mapping.values())[:5])}\n"
            f"       Likely causes: the VCF was built against a different reference build or "
            f"chromosome naming convention, or the mappings were regenerated after it."
        )
    if fraction < min_fraction:
        return False, (f"{summary}, below --min-match-fraction {min_fraction:.1%}.\n"
                       f"       unmatched e.g.: {', '.join(unmatched_examples)}")
    return True, summary


def main():
    parser = argparse.ArgumentParser(
        description="Preflight check that a VCF and a chromosome mapping describe the same variants."
    )
    parser.add_argument("-i", "--vcf", required=True, help="VCF to check (annotated or not)")
    parser.add_argument("-c", "--chromosome-mapping", required=True,
                        help="chr_mapping_<GENE>.csv from core/variant_mapping.py")
    parser.add_argument("-r", "--reference",
                        help="Reference FASTA (needs a .fai). Enables the same indel "
                             "left-alignment tier the parser uses; without it, repeat-adjacent "
                             "indels can read as unmatched here and still join later.")
    parser.add_argument("-g", "--gene", help="Gene label for messages")
    parser.add_argument("--min-match-fraction", type=float, default=0.0,
                        help="Fail when the matched fraction is below this (default 0.0: fail "
                             "only when NOTHING matches)")
    parser.add_argument("--warn-only", action="store_true",
                        help="Report and exit 0 regardless. Use to survey an existing tree.")
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        print(f"[check_vcf_mapping] ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    ok, message = check(args.vcf, args.chromosome_mapping, args.reference,
                        args.min_match_fraction, args.gene)
    if ok:
        print(f"[check_vcf_mapping] OK {message}")
        return 0
    print(f"[check_vcf_mapping] ERROR {message}", file=sys.stderr)
    if args.warn_only:
        print("[check_vcf_mapping] --warn-only: continuing anyway.", file=sys.stderr)
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
