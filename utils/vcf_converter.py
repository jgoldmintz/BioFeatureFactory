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
A robust script to convert gene-specific SNP files to standard VCF files.
Updated to support RefSeq chromosome naming for SpliceAI compatibility,
with exon-aware mapping and optional sequence verification.
"""

import argparse
import csv
import datetime
import json
import os
import sys
from pathlib import Path

import pysam

from utils.utility import (
    get_mutation_data_bioAccurate,
    trim_muts,
    extract_gene_from_filename,
    chromosome_map,
    get_genome_loc,
    load_mapping,
    read_fasta,
)

# Sequence verification

def verify_sequences_against_reference(verify_path, reference_path):
    try:
        verify_path = Path(verify_path)

        # Collect all FASTA files to verify
        fasta_files = []
        if verify_path.is_file():
            fasta_files = [verify_path]
        elif verify_path.is_dir():
            # Find all FASTA files in directory
            fasta_extensions = ['*.fasta', '*.fa', '*.fas', '*.fna', '*.faa']
            for ext in fasta_extensions:
                fasta_files.extend(verify_path.glob(ext))
            if not fasta_files:
                return {"success": False, "error": f"No FASTA files found in directory {verify_path}"}
        else:
            return {"success": False, "error": f"Path {verify_path} does not exist"}

        print(f"Verifying {len(fasta_files)} FASTA file(s) against reference...")

        ref_fasta = pysam.FastaFile(reference_path)
        results = {
            "success": True,
            "verified_sequences": {},
            "failed_sequences": {},
            "total_sequences": 0,
            "files_processed": 0,
        }

        for fasta_file in fasta_files:
            print(f"\nProcessing file: {fasta_file.name}")

            verify_sequences = read_fasta(str(fasta_file))
            if not verify_sequences:
                print(f"  Warning: No sequences found in {fasta_file.name}")
                continue

            results["files_processed"] += 1
            results["total_sequences"] += len(verify_sequences)

            for seq_name, sequence in verify_sequences.items():
                print(f"  Verifying {seq_name}...")
                found, found_locations = False, []

                for ref_name in ref_fasta.references:
                    ref_sequence = ref_fasta.fetch(ref_name)
                    pos = ref_sequence.upper().find(sequence.upper())
                    if pos != -1:
                        found = True
                        found_locations.append(
                            {"chromosome": ref_name, "position": pos + 1, "length": len(sequence)}
                        )

                if found:
                    # Include file info in the results
                    results["verified_sequences"][f"{fasta_file.name}:{seq_name}"] = {
                        "sequence_length": len(sequence),
                        "locations": found_locations,
                        "source_file": str(fasta_file)
                    }
                    print(f"    [OK] Found at {len(found_locations)} location(s)")
                else:
                    results["failed_sequences"][f"{fasta_file.name}:{seq_name}"] = {
                        "sequence_length": len(sequence),
                        "error": "Sequence not found in reference genome",
                        "source_file": str(fasta_file)
                    }
                    print("    [X] Not found in reference genome")

        ref_fasta.close()

        verified_count = len(results["verified_sequences"])
        failed_count = len(results["failed_sequences"])

        print(f"\nVerification summary:")
        print(f"  Files processed: {results['files_processed']}")
        print(f"  Total sequences: {results['total_sequences']}")
        print(f"  Verified: {verified_count}")
        print(f"  Failed: {failed_count}")

        if failed_count > 0:
            print(f"WARNING: {failed_count} sequences not found in reference genome")

        return results
    except Exception as e:
        return {"success": False, "error": str(e)}



# Helper functions

def chromosome_to_refseq(chromosome, version="GRCh38"):
    return chromosome_map.get(version, {}).get(str(chromosome), str(chromosome))


def format_chromosome(chromosome, format_type="refseq"):
    if format_type == "refseq":
        return chromosome_to_refseq(chromosome)
    elif format_type == "ucsc":
        return f"chr{chromosome}"
    return str(chromosome)


def reverse_complement(nucleotide):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return comp.get(nucleotide.upper(), nucleotide)


def get_reference_nucleotide(fasta_file, chromosome, position):
    try:
        fasta = pysam.FastaFile(fasta_file)
        candidates = [str(chromosome), f"chr{chromosome}"]
        for build in ["GRCh38", "GRCh37"]:
            if str(chromosome) in chromosome_map.get(build, {}):
                candidates.append(chromosome_map[build][str(chromosome)])
        for ref_name in fasta.references:
            for build in ["GRCh38", "GRCh37"]:
                if str(chromosome) in chromosome_map.get(build, {}):
                    base = chromosome_map[build][str(chromosome)].split(".")[0]
                    if ref_name.startswith(base):
                        candidates.append(ref_name)
        for cand in dict.fromkeys(candidates):  # preserve order, unique
            if cand in fasta.references:
                return fasta.fetch(cand, position - 1, position).upper()
        return None
    except Exception as e:
        print(f"Error extracting reference nucleotide: {e}", file=sys.stderr)
        return None


def _resolve_mapping_source(chromosome_mapping_input, gene_name):
    """Return the mapping CSV path for a gene (if one should exist)."""
    if not chromosome_mapping_input:
        return None
    candidate = Path(chromosome_mapping_input)
    if candidate.is_dir():
        return candidate / f"combined_{gene_name}.csv"
    return candidate


def _collect_log_metadata(log_path):
    """Collect resolved log paths and their mtimes for cache invalidation."""
    paths = []
    mtimes = {}
    if not log_path:
        return paths, mtimes

    candidate = Path(log_path)
    if candidate.is_dir():
        resolved_dir = str(candidate.resolve())
        try:
            mtimes[resolved_dir] = candidate.stat().st_mtime
        except OSError:
            mtimes[resolved_dir] = None
        for log_file in sorted(candidate.glob("*.log")):
            resolved = str(log_file.resolve())
            paths.append(resolved)
            try:
                mtimes[resolved] = log_file.stat().st_mtime
            except OSError:
                mtimes[resolved] = None
    else:
        resolved = str(candidate.resolve())
        paths.append(resolved)
        try:
            mtimes[resolved] = candidate.stat().st_mtime
        except OSError:
            mtimes[resolved] = None

    return paths, mtimes


def _load_cache(cache_path):
    """Load the JSON cache if present, otherwise return an empty dict."""
    if cache_path.exists():
        try:
            with open(cache_path, "r") as handle:
                return json.load(handle)
        except (OSError, json.JSONDecodeError):
            return {}
    return {}


def _save_cache(cache_path, cache):
    """Persist the cache dictionary to disk (best-effort)."""
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as handle:
            json.dump(cache, handle, indent=2, sort_keys=True)
    except OSError as exc:
        print(f"Warning: Unable to write cache file {cache_path}: {exc}", file=sys.stderr)


def _is_cache_valid(entry, info):
    """Return True if the cached entry still matches the current input state."""
    try:
        if entry.get("source_mtime") != info.get("source_mtime"):
            return False

        output_path = entry.get("output_path")
        if not output_path or not Path(output_path).exists():
            return False
        if entry.get("output_mtime") != Path(output_path).stat().st_mtime:
            return False

        if entry.get("mapping_mtime") != info.get("mapping_mtime"):
            return False

        if entry.get("log_paths") != info.get("log_paths"):
            return False

        if entry.get("log_mtimes") != info.get("log_mtimes"):
            return False

        return True
    except OSError:
        return False


# Core file processing

def process_single_file(
    mutations_file,
    outdir,
    chromosome_format="refseq",
    reference_fasta=None,
    annotation_file=None,
    chromosome_mapping_input=None,
    validate_mapping=False,
    log_path=None,
):
    try:
        gene_name = extract_gene_from_filename(mutations_file)
        print(f"Processing gene '{gene_name}' from '{mutations_file}'...")

        gene_info = get_genome_loc(gene_name, annotation_file)
        if not gene_info:
            return False, 0, f"Could not retrieve genome location for {gene_name}", None

        chromosome = gene_info["chrom"]
        gene_start, gene_end, strand = gene_info["tx_start"], gene_info["tx_end"], gene_info["strand"]
        exons = gene_info["exons"]

        print(
            f"Located {gene_name} on chr{chromosome}:{gene_start}-{gene_end} ({strand}), "
            f"{len(exons)} exon(s)"
        )

        # Load mutations
        mut_list = trim_muts(mutations_file, log=log_path, gene_name=gene_name)

        if not mut_list:
            return False, 0, f"No mutations loaded for {gene_name}", None

        formatted_chr = format_chromosome(chromosome, chromosome_format)
        variants, processed_count = [], 0

        mapping_lookup = {}
        mapping_source = _resolve_mapping_source(chromosome_mapping_input, gene_name)
        if mapping_source:
            if mapping_source.exists():
                mapping_lookup = load_mapping(str(mapping_source), mapType="chromosome")
                if not mapping_lookup:
                    print(
                        f"  Warning: Mapping file '{mapping_source}' contained no usable chromosome entries.",
                        file=sys.stderr,
                    )
            else:
                print(
                    f"  Warning: No chromosome mapping file found for {gene_name} at {mapping_source}.",
                    file=sys.stderr,
                )

        for snp_string in mut_list:
            try:
                relative_pos, nts = get_mutation_data_bioAccurate(snp_string)
                if relative_pos is None or not nts:
                    continue
                wt_nt, mut_nt = nts
                mapping_entry = mapping_lookup.get(snp_string) if mapping_lookup else None
                if mapping_entry:
                    m_ref = mapping_entry[0]
                    m_alt = mapping_entry[-1]
                    try:
                        m_pos = int(mapping_entry[1:-1])
                    except ValueError:
                        print(
                            f"Skipping mapping for {snp_string}: could not parse position in '{mapping_entry}'.",
                            file=sys.stderr,
                        )
                        mapping_entry = None
                    else:
                        abs_pos = m_pos
                        ref, alt = m_ref.upper(), m_alt.upper()
                        if validate_mapping and reference_fasta:
                            actual_ref = get_reference_nucleotide(reference_fasta, chromosome, abs_pos)
                            if actual_ref and actual_ref != ref:
                                print(
                                    f"  Warning: Mapping REF mismatch for {snp_string} (mapping {ref}, reference {actual_ref}). Recomputing from annotation.",
                                    file=sys.stderr,
                                )
                                mapping_entry = None
                        if mapping_entry:
                            variants.append((abs_pos, formatted_chr, ref, alt))
                            processed_count += 1
                            continue
                # fallback to annotation-driven computation
                if relative_pos is None:
                    continue
                if strand == "+":
                    abs_pos = gene_start + relative_pos - 1
                    ref, alt = wt_nt, mut_nt
                else:
                    abs_pos = gene_end - relative_pos + 1
                    ref, alt = reverse_complement(wt_nt), reverse_complement(mut_nt)
                if reference_fasta:
                    actual_ref = get_reference_nucleotide(reference_fasta, chromosome, abs_pos)
                    if actual_ref:
                        ref = actual_ref
                variants.append((abs_pos, formatted_chr, ref, alt))
                processed_count += 1
            except Exception as e:
                print(f"Skipping mutation {snp_string}: {e}", file=sys.stderr)

        # Sort and write VCF
        variants.sort(key=lambda x: x[0])
        out_vcf = Path(outdir) / f"{gene_name}.vcf"
        with open(out_vcf, "w") as f:
            f.write("##fileformat=VCFv4.3\n")
            f.write(f"##fileDate={datetime.datetime.now():%Y%m%d}\n")
            f.write("##source=SnpToVcfConverter\n")
            f.write(f"##gene_context={gene_name}\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for pos, chrom, ref, alt in variants:
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

        return True, processed_count, None, str(out_vcf)
    except Exception as e:
        return False, 0, str(e), None



# Main

def main():
    parser = argparse.ArgumentParser(
        description="Convert gene-specific SNP CSVs to VCF with RefSeq support and exon-aware mapping"
    )
    parser.add_argument("-m", "--mutation", required=True, help="Path to mutation CSV or dir")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--chromosome-format", choices=["simple", "refseq", "ucsc"], default="refseq"
    )
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA")
    parser.add_argument("-a", "--annotation", required=True, help="Gene annotation file")
    parser.add_argument(
        "--chromosome-mapping-input",
        help="Path to a chromosome mapping CSV or directory of combined_<GENE>.csv files",
    )
    parser.add_argument(
        "--validate-mapping",
        action="store_true",
        help="Cross-check supplied chromosome mappings against the reference FASTA",
    )
    parser.add_argument(
        "--log",
        help="Validation log (or directory of logs) listing mutations to skip during conversion.",
    )
    parser.add_argument(
        "--clear-cache",
        action="store_true",
        help="Remove cached VCF metadata before processing (forces regeneration).",
    )
    parser.add_argument("--verify-sequences", help="FASTA file or directory containing FASTA files to verify against reference")
    args = parser.parse_args()

    mutation_path = Path(args.mutation)
    outdir_path = Path(args.outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    reference_fasta, annotation_file = args.reference, args.annotation
    failed_files, successful_files, total_processed = [], 0, 0

    if args.verify_sequences:
        results = verify_sequences_against_reference(args.verify_sequences, reference_fasta)
        if not results["success"]:
            print(f"Sequence verification failed: {results['error']}", file=sys.stderr)
            sys.exit(1)

    files = [mutation_path] if mutation_path.is_file() else list(mutation_path.glob("*.csv"))

    cache_path = outdir_path / ".vcf_converter_cache.json"
    if args.clear_cache and cache_path.exists():
        try:
            cache_path.unlink()
        except OSError as exc:
            print(f"Warning: Unable to clear cache {cache_path}: {exc}", file=sys.stderr)
    cache = _load_cache(cache_path)
    cache_dirty = False

    log_paths_global, log_mtimes_global = _collect_log_metadata(args.log)
    reference_resolved = str(Path(reference_fasta).resolve())
    annotation_resolved = str(Path(annotation_file).resolve())

    for f in files:
        gene_name = extract_gene_from_filename(str(f))
        mutation_file = Path(f).resolve()
        mapping_source = _resolve_mapping_source(args.chromosome_mapping_input, gene_name)
        mapping_source_str = None
        mapping_mtime = None
        if mapping_source:
            try:
                mapping_source_str = str(mapping_source.resolve())
            except OSError:
                mapping_source_str = str(mapping_source)
            if mapping_source.exists():
                try:
                    mapping_mtime = mapping_source.stat().st_mtime
                except OSError:
                    mapping_mtime = None

        key_payload = {
            "gene": gene_name,
            "mutation_file": str(mutation_file),
            "reference": reference_resolved,
            "annotation": annotation_resolved,
            "chromosome_format": args.chromosome_format,
            "chromosome_mapping_input": str(args.chromosome_mapping_input) if args.chromosome_mapping_input else None,
            "mapping_source": mapping_source_str,
            "validate_mapping": bool(args.validate_mapping),
            "log_input": str(args.log) if args.log else None,
            "log_paths": list(log_paths_global),
        }
        cache_key = json.dumps(key_payload, sort_keys=True)

        try:
            source_mtime = mutation_file.stat().st_mtime
        except OSError:
            source_mtime = None

        info = {
            "source_mtime": source_mtime,
            "mapping_mtime": mapping_mtime,
            "log_paths": list(log_paths_global),
            "log_mtimes": dict(log_mtimes_global),
        }

        cached_entry = cache.get(cache_key)
        if cached_entry and _is_cache_valid(cached_entry, info):
            successful_files += 1
            total_processed += cached_entry.get("variant_count", 0)
            print(
                f"[OK] {Path(f).name} -> {cached_entry.get('output_path')} (cached {cached_entry.get('variant_count', 0)} variants)"
            )
            continue

        success, count, err, vcf = process_single_file(
            str(f),
            str(outdir_path),
            chromosome_format=args.chromosome_format,
            reference_fasta=reference_fasta,
            annotation_file=annotation_file,
            chromosome_mapping_input=args.chromosome_mapping_input,
            validate_mapping=args.validate_mapping,
            log_path=args.log,
        )
        if success:
            successful_files += 1
            total_processed += count
            print(f"[OK] {Path(f).name} -> {vcf} ({count} variants)")
            try:
                output_path = Path(vcf).resolve()
                output_mtime = output_path.stat().st_mtime
            except OSError:
                output_path = Path(vcf)
                output_mtime = None
            cache[cache_key] = {
                "source_mtime": source_mtime,
                "output_path": str(output_path),
                "output_mtime": output_mtime,
                "mapping_mtime": mapping_mtime,
                "log_paths": list(log_paths_global),
                "log_mtimes": dict(log_mtimes_global),
                "variant_count": count,
            }
            cache_dirty = True
        else:
            failed_files.append((Path(f).name, err))
            print(f"[X] {Path(f).name}: {err}", file=sys.stderr)
            if cache_key in cache:
                cache_dirty = True
                cache.pop(cache_key, None)

    print(
        f"\nSummary: {successful_files}/{len(files)} files successful "
        f"({total_processed} variants processed)"
    )
    if failed_files:
        print(f"Failed files: {failed_files}")

    if cache_dirty:
        _save_cache(cache_path, cache)


if __name__ == "__main__":
    main()
