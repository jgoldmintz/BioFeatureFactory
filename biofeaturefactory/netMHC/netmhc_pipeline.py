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
NetMHC Pipeline for MHC Binding Prediction

Predicts MHC class I and II binding peptides for WT and mutant protein sequences.
Supports both Docker and native execution modes.
Generates ensemble TSV outputs with WT vs MUT comparisons.

Key features:
- MHC class I and II binding prediction
- Multiple HLA allele support
- WT vs mutant comparison with delta scores
- Batch processing for large FASTA files
- Integration with mutation mapping CSVs
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
import sys
import time
import platform
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Dict, List, Tuple


# Import utility functions
from biofeaturefactory.utils.utility import (
    write_fasta,
    write_tsv,
    combine_batch_outputs,
    discover_fasta_files,
    ExtractGeneFromFASTA,
    parse_predictions_with_mutation_filtering,
    load_validation_failures,
    should_skip_mutation,
    load_wt_sequences,
    load_wt_sequence_map,
    extract_gene_from_filename,
    extract_mutation_from_sequence_name,
    translate_orf_sequence,
    build_mutant_sequences_for_gene,
    detect_alphabet,
)


def is_linux_host():
    """Return True when running on a Linux kernel."""
    return platform.system().lower() == "linux"


def resolve_native_netmhc_path(user_path=None, tool_version="netMHC"):
    """
    Resolve a usable native NetMHC executable when available.

    Search order:
    1. Explicit --native-netmhc-path value
    2. $NETMHC_PATH or $NETMHCPAN_PATH environment variable
    3. Common install locations

    Args:
        user_path: User-specified path to NetMHC binary
        tool_version: Which NetMHC tool to use. Only netMHC-4.0 is supported —
            the invocation and the 15-column/'HLA'-header parser are hardcoded to
            it (see --netmhc-tool, choices=['netMHC']). The default is 'netMHC' so
            a programmatic caller that omits it cannot silently resolve netMHCpan,
            whose output this parser reads as zero predictions (F5).
    """
    candidates = []

    def _add(path):
        if path:
            candidates.append(os.path.expanduser(path))

    _add(user_path)
    _add(os.environ.get("NETMHC_PATH"))

    netmhc_home = os.environ.get("NETMHC_HOME")
    if netmhc_home:
        _add(os.path.join(netmhc_home, tool_version))

    home = Path.home()
    # netMHC-4.0 only (see --netmhc-tool): do NOT auto-discover netMHCpan — its
    # output format is unreadable by this parser and would silently yield zero
    # predictions.
    common_roots = [
        home / "netMHC" / tool_version,
        Path(f"/opt/netMHC/{tool_version}"),
        Path(f"/usr/local/bin/{tool_version}"),
    ]

    for candidate in common_roots:
        _add(str(candidate))

    for path in candidates:
        if not path:
            continue
        # If given a directory, look for the binary inside it
        if os.path.isdir(path):
            for subpath in [
                os.path.join(path, "bin", tool_version),
                os.path.join(path, tool_version),
            ]:
                if os.path.isfile(subpath) and os.access(subpath, os.X_OK):
                    return os.path.abspath(subpath)
        elif os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)

    return None


def build_netmhc_executor(args, parser):
    """
    Resolve the native NetMHC binary and return a callable executor.

    Returns a callable executor and a description string.
    """
    native_netmhc = resolve_native_netmhc_path(
        getattr(args, "native_netmhc_path", None),
        getattr(args, "netmhc_tool", "netMHC")   # F5: netMHC-4.0 is the only supported tool
    )
    if not native_netmhc:
        parser.error(
            "Native NetMHC binary not found. Provide --native-netmhc-path, "
            "or set NETMHC_PATH / NETMHC_HOME environment variable."
        )

    verbose_flag = getattr(args, "verbose", True)
    if verbose_flag:
        print(f"Execution mode: native NetMHC ({native_netmhc})")

    def _runner(fasta_file, output_file, timeout, alleles):
        return _run_native_netmhc(fasta_file, output_file, timeout, native_netmhc, alleles)
    return _runner, f"native:{native_netmhc}"


def _run_native_netmhc(fasta_file, output_file, timeout, netmhc_path, alleles=None):
    """
    Execute NetMHC using native installation.

    Args:
        fasta_file: Input FASTA file path
        output_file: Output file path
        timeout: Command timeout in seconds
        netmhc_path: Path to NetMHC executable
        alleles: List of HLA alleles to predict

    Returns:
        tuple: (success, output_content, error_message)
    """
    if not alleles:
        alleles = ["HLA-A0201"]  # Default allele (netMHC-4.0 format)

    all_outputs = []

    try:
        # NetMHC processes one allele at a time, so run for each allele
        for allele in alleles:
            # Build native NetMHC command
            # Format: netMHC -a HLA-A*02:01 -f input.fasta
            cmd = [netmhc_path, "-a", allele, "-f", fasta_file]

            # netMHC binary expects $NETMHC pointing to the platform dir
            # (the parent of bin/), e.g. netMHC-4.0/Darwin_x86_64/
            env = os.environ.copy()
            netmhc_bin_dir = os.path.dirname(os.path.abspath(netmhc_path))
            env["NETMHC"] = os.path.dirname(netmhc_bin_dir) if os.path.basename(netmhc_bin_dir) == "bin" else netmhc_bin_dir

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                env=env,
            )

            # netMHC-4.0 returns exit code 1 even on success;
            # check stderr and stdout for actual error indicators
            if result.stderr and "error" in result.stderr.lower():
                return False, result.stdout, result.stderr
            if "cannot be found" in result.stdout or (not result.stdout.strip()):
                return False, result.stdout, result.stderr or result.stdout

            all_outputs.append(result.stdout)

        # Combine all outputs
        combined_output = "\n".join(all_outputs)

        # Write to output file
        with open(output_file, 'w') as f:
            f.write(combined_output)

        return True, combined_output, None

    except subprocess.TimeoutExpired:
        return False, "", f"NetMHC command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)


def parse_netmhc_output(output_file):
    """
    Parse NetMHC output file and extract binding predictions.

    NetMHC output format (space-separated):
    pos HLA peptide Core Offset I_pos I_len D_pos D_len iCore Identity 1-log50k(aff) Affinity(nM) %Rank BindLevel

    Returns:
        list: List of prediction dictionaries with keys:
              pos, mhc_allele, peptide, core, affinity, rank, bind_level, identity
    """
    predictions = []

    try:
        with open(output_file, 'r') as f:
            in_prediction_section = False
            current_identity = ""

            for line in f:
                original_line = line
                line = line.strip()

                # Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue

                # Skip separator lines (all dashes)
                if line.startswith('---'):
                    continue

                # Detect prediction section header
                if 'pos' in line.lower() and 'peptide' in line.lower() and 'HLA' in line:
                    in_prediction_section = True
                    continue

                # Stop at summary line
                if line.startswith('Protein ') and 'Allele' in line:
                    in_prediction_section = False
                    continue

                if in_prediction_section:
                    # Parse prediction line
                    # Format: "  0  HLA-A0201  TMDKSELVQ  ...  28676.59  43.00"
                    # Or:     "219  HLA-A0201  QLLRDNLTL  ...   167.10   1.50 <= WB"

                    fields = line.split()

                    # Need at least 14 fields for valid prediction
                    if len(fields) < 14:
                        continue

                    try:
                        # Extract identity (sequence name) from field 10
                        identity = fields[10] if len(fields) > 10 else ""

                        # BindLevel is optional and appears as "<= WB" or "<= SB" (2 tokens)
                        bind_level = ""
                        if len(fields) >= 16 and fields[14] == "<=":
                            bind_level = fields[15]  # WB or SB
                        elif len(fields) == 15 and fields[14] not in ["<=", ""]:
                            # Sometimes just "WB" or "SB" without <=
                            bind_level = fields[14]

                        prediction = {
                            'pos': int(fields[0]),
                            'mhc_allele': fields[1],
                            'peptide': fields[2],
                            'core': fields[3],
                            'affinity': float(fields[12]),  # Affinity(nM)
                            'rank': float(fields[13]),      # %Rank
                            'bind_level': bind_level,
                            'identity': identity,
                        }
                        predictions.append(prediction)
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue

    except Exception as e:
        print(f"Error parsing NetMHC output {output_file}: {e}")
        return []

    return predictions


def compare_wt_mut_predictions(gene_name, mutation, wt_preds, mut_preds, threshold=0.5):
    """
    Compare WT and MUT predictions to classify epitope changes.

    Classification logic:
    - Gained: wt_rank > 2.0 AND mut_rank <= threshold (new strong binder)
    - Lost: wt_rank <= threshold AND mut_rank > 2.0 (lost strong binder)
    - Strengthened: both bind AND delta_rank < -5 (improved binding)
    - Weakened: both bind AND delta_rank > 5 (reduced binding)
    - Stable: both bind AND abs(delta_rank) <= 5 (no significant change)

    Args:
        gene_name: Gene name
        mutation: Mutation identifier
        wt_preds: List of WT predictions
        mut_preds: List of MUT predictions
        threshold: Rank threshold for strong binder (default 0.5)

    Returns:
        List of event dictionaries
    """
    events = []

    # Build lookup maps keyed on the REGISTER (allele, pos) — NOT the peptide.
    # A window covering the mutated residue has different WT vs MUT peptide
    # strings by construction; keying on the sequence would give them different
    # keys so they would never be paired (the delta would be computed against the
    # non-binder sentinel instead of the matched register). pos aligns 1:1 because
    # SNV/missense preserves length.
    wt_map = {}
    for pred in wt_preds:
        key = (pred['mhc_allele'], pred['pos'])
        wt_map[key] = pred

    mut_map = {}
    for pred in mut_preds:
        key = (pred['mhc_allele'], pred['pos'])
        mut_map[key] = pred

    # Every (allele, pos) register seen in either allele.
    all_keys = set(wt_map.keys()) | set(mut_map.keys())

    for key in all_keys:
        allele, pos = key
        wt_pred = wt_map.get(key)
        mut_pred = mut_map.get(key)
        wt_peptide = wt_pred['peptide'] if wt_pred else ''
        mut_peptide = mut_pred['peptide'] if mut_pred else ''

        # Skip if both missing (shouldn't happen)
        if not wt_pred and not mut_pred:
            continue

        # Extract ranks and affinities
        wt_rank = wt_pred['rank'] if wt_pred else 100.0  # Non-binder
        mut_rank = mut_pred['rank'] if mut_pred else 100.0
        wt_affinity = wt_pred['affinity'] if wt_pred else 50000.0
        mut_affinity = mut_pred['affinity'] if mut_pred else 50000.0
        wt_bind = wt_pred['bind_level'] if wt_pred else ''
        mut_bind = mut_pred['bind_level'] if mut_pred else ''

        # Calculate delta
        delta_rank = mut_rank - wt_rank
        delta_affinity = mut_affinity - wt_affinity

        # Classify by NetMHC binding band (per the tool's own definitions):
        #   2 = SB (strong binder, %Rank < threshold, default 0.5)
        #   1 = WB (weak binder,   threshold <= %Rank < 2.0)
        #   0 = NB (non-binder,    %Rank >= 2.0)
        # Strength order SB(2) > WB(1) > NB(0). Classification is a band transition;
        # delta_rank/delta_affinity are retained as reported magnitudes only.
        # (The former delta_rank <>±5 test was dead: inside the both-bind branch
        # delta_rank is bounded to [-2, 2] and could never reach ±5.)
        wt_band = 2 if wt_rank < threshold else (1 if wt_rank < 2.0 else 0)
        mut_band = 2 if mut_rank < threshold else (1 if mut_rank < 2.0 else 0)

        if wt_band == 0 and mut_band > 0:
            classification = "gained"
            classification_code = 2
        elif wt_band > 0 and mut_band == 0:
            classification = "lost"
            classification_code = -2
        elif mut_band > wt_band:
            classification = "strengthened"
            classification_code = 1
        elif mut_band < wt_band:
            classification = "weakened"
            classification_code = -1
        else:
            classification = "stable"
            classification_code = 0

        event = {
            'Gene': gene_name,
            'pkey': f"{gene_name}-{mutation}",
            'mutation': mutation,
            'wt_peptide': wt_peptide,
            'mut_peptide': mut_peptide,
            'pos': pos,
            'mhc_allele': allele,
            'wt_rank': wt_rank,
            'mut_rank': mut_rank,
            'delta_rank': delta_rank,
            'wt_affinity': wt_affinity,
            'mut_affinity': mut_affinity,
            'delta_affinity': delta_affinity,
            'bind_level_wt': wt_bind,
            'bind_level_mut': mut_bind,
            'classification': classification,
            'classification_code': classification_code,
        }

        events.append(event)

    return events


def summarize_epitope_changes(gene_name, mutation, events):
    """
    Generate summary statistics for a mutation's epitope changes.

    Args:
        gene_name: Gene name
        mutation: Mutation identifier
        events: List of epitope event dictionaries

    Returns:
        Dictionary with summary statistics
    """
    # Count events by classification
    count_gained = sum(1 for e in events if e['classification'] == 'gained')
    count_lost = sum(1 for e in events if e['classification'] == 'lost')
    count_strengthened = sum(1 for e in events if e['classification'] == 'strengthened')
    count_weakened = sum(1 for e in events if e['classification'] == 'weakened')
    count_stable = sum(1 for e in events if e['classification'] == 'stable')

    # Count epitopes (rank <= 0.5)
    n_epitopes_wt = sum(1 for e in events if e['wt_rank'] <= 0.5)
    n_epitopes_mut = sum(1 for e in events if e['mut_rank'] <= 0.5)

    # Find maximum absolute delta_rank and sum
    abs_deltas = [abs(e['delta_rank']) for e in events]
    max_abs_delta = max(abs_deltas) if abs_deltas else 0.0
    sum_abs_delta = sum(abs_deltas)

    # Find top event (largest absolute delta)
    top_event = None
    if events:
        top_event = max(events, key=lambda e: abs(e['delta_rank']))

    # QC flags
    qc_flags = []
    if n_epitopes_wt == 0:
        qc_flags.append("no_wt_epitopes")
    if n_epitopes_mut == 0:
        qc_flags.append("no_mut_epitopes")
    if not events:
        qc_flags.append("no_predictions")

    summary = {
        'pkey': f"{gene_name}-{mutation}",
        'Gene': gene_name,
        'mutation': mutation,
        'n_epitopes_wt': n_epitopes_wt,
        'n_epitopes_mut': n_epitopes_mut,
        'count_gained': count_gained,
        'count_lost': count_lost,
        'count_strengthened': count_strengthened,
        'count_weakened': count_weakened,
        'count_stable': count_stable,
        'max_abs_delta_rank': max_abs_delta,
        'sum_abs_delta_rank': sum_abs_delta,
        'top_event_type': top_event['classification'] if top_event else '',
        'top_event_allele': top_event['mhc_allele'] if top_event else '',
        'top_event_wt_peptide': top_event['wt_peptide'] if top_event else '',
        'top_event_mut_peptide': top_event['mut_peptide'] if top_event else '',
        'top_event_delta_rank': top_event['delta_rank'] if top_event else 0.0,
        'qc_flags': ';'.join(qc_flags) if qc_flags else '',
    }

    return summary


def write_summary_tsv(summary_rows, output_file):
    """Write summary TSV file."""
    if not summary_rows:
        print(f"Warning: No summary data to write")
        return

    fieldnames = [
        'pkey', 'Gene', 'mutation',
        'n_epitopes_wt', 'n_epitopes_mut',
        'count_gained', 'count_lost', 'count_strengthened', 'count_weakened', 'count_stable',
        'max_abs_delta_rank', 'sum_abs_delta_rank',
        'top_event_type', 'top_event_allele', 'top_event_wt_peptide', 'top_event_mut_peptide', 'top_event_delta_rank',
        'qc_flags'
    ]

    write_tsv(summary_rows, output_file, fieldnames, extrasaction='ignore')
    print(f"Wrote {len(summary_rows)} summary rows to {output_file}")


def write_events_tsv(events, output_file):
    """Write events TSV file."""
    if not events:
        print(f"Warning: No event data to write")
        return

    fieldnames = [
        'pkey', 'Gene', 'mutation', 'wt_peptide', 'mut_peptide', 'pos', 'mhc_allele',
        'wt_rank', 'mut_rank', 'delta_rank',
        'wt_affinity', 'mut_affinity', 'delta_affinity',
        'bind_level_wt', 'bind_level_mut',
        'classification', 'classification_code'
    ]

    write_tsv(events, output_file, fieldnames, extrasaction='ignore')
    print(f"Wrote {len(events)} events to {output_file}")


def write_sites_tsv(sites, output_file):
    """Write sites TSV file."""
    if not sites:
        print(f"Warning: No site data to write")
        return

    fieldnames = [
        'pkey', 'Gene', 'sequence_type',
        'pos', 'mhc_allele', 'peptide', 'core',
        'affinity', 'rank', 'bind_level', 'identity'
    ]

    write_tsv(sites, output_file, fieldnames, extrasaction='ignore')
    print(f"Wrote {len(sites)} site predictions to {output_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="NetMHC pipeline for MHC binding prediction with WT/mutant comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('input', nargs='?',
                       help='Input: WT FASTA file/directory')
    parser.add_argument('output', nargs='?',
                       help='Output base directory')

    # MHC-specific options
    parser.add_argument('--alleles', nargs='+',
                       help='HLA alleles to predict (e.g., HLA-A*02:01 HLA-B*07:02). If not specified, uses default set.')
    # Only netMHC-4.0 is supported: the invocation (allele format, `-a`/`-f`) and
    # the parser (15-column layout, 'HLA' header) are hardcoded to it. netMHCpan /
    # netMHCII use different output layouts the parser cannot read and would
    # silently yield zero predictions, so they are not offered.
    parser.add_argument('--netmhc-tool', choices=['netMHC'],
                       default='netMHC',
                       help='NetMHC tool to use (only netMHC-4.0 is supported)')
    # Execution backend
    parser.add_argument('--native-netmhc-path',
                       help='Path to native NetMHC executable')

    # Processing options
    parser.add_argument('--mutations', '-m',
                       help='Mutation file or directory of mutation CSVs (single-column NT mutations)')
    parser.add_argument('--log',
                       help='Validation log file to skip failed mutations')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Binding rank threshold for strong binders (default: 0.5)')
    parser.add_argument('--batch-size', type=int, default=100,
                       help='(deprecated; unused — netMHC input is no longer batched)')
    parser.add_argument('--timeout', type=int, default=600,
                       help='Command timeout in seconds (default: 600)')
    parser.add_argument('--max-workers', type=int, default=4,
                       help='Concurrent netMHC runs per gene (default: 4; set 1 for serial)')

    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Validate arguments
    if not args.input or not args.output:
        parser.error("input and output arguments are required")

    if not args.mutations:
        parser.error("--mutations is REQUIRED for full-pipeline mode")

    # Build NetMHC executor
    executor, exec_desc = build_netmhc_executor(args, parser)

    # Load validation failures if provided
    failure_map = None
    if args.log:
        failure_map = load_validation_failures(args.log)
        if args.verbose:
            print(f"Loaded validation failures from {args.log}")

    # Load WT sequences - handles both files and directories
    wt_sequences, temp_holder = load_wt_sequence_map(args.input, wt_header='ORF')
    if not wt_sequences:
        print(f"Error: No FASTA sequences found in {args.input}")
        return 1

    if args.verbose:
        print(f"Loaded {len(wt_sequences)} WT sequences")

    # Discover mutation files
    mutation_files = {}
    if args.mutations:
        mut_path = Path(args.mutations)
        if mut_path.is_file():
            mutation_files[extract_gene_from_filename(mut_path.name)] = str(mut_path)
        elif mut_path.is_dir():
            for csv_file in mut_path.glob("*.csv"):
                mutation_files[extract_gene_from_filename(csv_file.name)] = str(csv_file)

    # Process each gene
    for gene_name, wt_seq in wt_sequences.items():
        gene_summary_rows = []
        gene_events = []
        gene_sites = []
        if args.verbose:
            print(f"\nProcessing gene: {gene_name}")

        # Auto-detect the WT alphabet: nucleotide -> translate; protein -> use
        # directly; codon-encoded -> skip. Prevents the unconditional-translate
        # crash (TranslationError) that amino-acid input previously caused.
        try:
            detected = detect_alphabet(wt_seq)
        except ValueError:
            print(f"Warning: {gene_name} WT sequence is empty, skipping")
            continue
        if detected == 'nucleotide':
            wt_nt_seq = wt_seq
            wt_aa_seq = translate_orf_sequence(wt_nt_seq)
            build_input_type = 'nt'
            if not wt_aa_seq:
                print(f"Warning: Could not translate {gene_name}, skipping")
                continue
        elif detected == 'protein':
            wt_nt_seq = None
            wt_aa_seq = wt_seq
            build_input_type = 'aa'
        else:  # codon-encoded
            print(f"Warning: {gene_name} WT looks codon-encoded; netMHC needs nt or aa, skipping")
            continue

        # Build mutant sequences
        mapping_file = mutation_files.get(gene_name)
        mutant_seqs = build_mutant_sequences_for_gene(
            gene_name, wt_nt_seq, wt_aa_seq, mapping_file, args.log, failure_map,
            input_type=build_input_type
        )

        if args.verbose:
            print(f"  Generated {len(mutant_seqs)} mutant sequences")

        # Create temp directory for this gene
        gene_workdir = tempfile.mkdtemp(prefix=f"netmhc_{gene_name}_")

        try:
            # Write WT FASTA
            wt_fasta = os.path.join(gene_workdir, f"{gene_name}_wt.fasta")
            write_fasta(Path(wt_fasta), {f"{gene_name}_WT": wt_aa_seq})

            # Run netMHC on the single-record WT FASTA. netMHC digests the
            # sequence into all k-mer windows itself, so no input batching is
            # needed (audit F8 — the former split-by-record batching was a no-op).
            wt_output = os.path.join(gene_workdir, f"{gene_name}_wt.out")
            success, _, error = executor(wt_fasta, wt_output, args.timeout, args.alleles)
            if not success:
                print(f"Warning: NetMHC failed for {gene_name} WT: {error}")
                continue
            wt_predictions = parse_netmhc_output(wt_output)
            if args.verbose:
                print(f"  WT: {len(wt_predictions)} predictions")

            # Add to sites output (all WT predictions)
            for pred in wt_predictions:
                gene_sites.append({
                    'Gene': gene_name,
                    'pkey': f"{gene_name}-WT",
                    'sequence_type': 'wt',
                    **pred
                })

            # Process each mutant in parallel. Every netMHC run is a separate OS
            # process (subprocess), so threads give real concurrency (the GIL is
            # released during subprocess.run) while keeping ONE netMHC invocation
            # per mutant — unambiguous output attribution and per-mutant failure
            # isolation. Genes with thousands of mutations no longer pay a fully
            # serial process-spawn cost.
            def _process_mutant(item):
                mut_header, mut_aa_seq = item
                # Key is f"{gene_name}-{mutant_clean}" (utility.py:1544); strip the
                # exact gene-name prefix so hyphenated genes (NKX2-1, HLA-A, MT-CO1)
                # parse correctly instead of split('-',1) corrupting the token.
                mutation = mut_header[len(gene_name) + 1:] if mut_header.startswith(f"{gene_name}-") else mut_header
                try:
                    mut_fasta = os.path.join(gene_workdir, f"{gene_name}_{mutation}.fasta")
                    write_fasta(Path(mut_fasta), {mut_header: mut_aa_seq})
                    mut_output = os.path.join(gene_workdir, f"{gene_name}_{mutation}.out")
                    success, _, error = executor(mut_fasta, mut_output, args.timeout, args.alleles)
                    if not success:
                        print(f"Warning: NetMHC failed for {gene_name} {mutation}: {error}")
                        return None
                    mut_predictions = parse_netmhc_output(mut_output)
                    if args.verbose:
                        print(f"  {mutation}: {len(mut_predictions)} predictions")
                    site_rows = [
                        {'Gene': gene_name, 'pkey': f"{gene_name}-{mutation}",
                         'sequence_type': 'mut', **pred}
                        for pred in mut_predictions
                    ]
                    events = compare_wt_mut_predictions(
                        gene_name, mutation, wt_predictions, mut_predictions, args.threshold
                    )
                    summary = summarize_epitope_changes(gene_name, mutation, events)
                    return site_rows, events, summary
                except Exception as e:
                    print(f"Warning: netMHC mutant {gene_name} {mutation} failed: {e}")
                    return None

            # map() preserves input order (deterministic output) and runs concurrently.
            with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as pool:
                for result in pool.map(_process_mutant, list(mutant_seqs.items())):
                    if result is None:
                        continue
                    site_rows, events, summary = result
                    gene_sites.extend(site_rows)
                    gene_events.extend(events)
                    gene_summary_rows.append(summary)

        finally:
            # Clean up temp directory
            if os.path.exists(gene_workdir):
                shutil.rmtree(gene_workdir, ignore_errors=True)

        gene_out = Path(args.output) / gene_name / "NetMHC"
        gene_out.mkdir(parents=True, exist_ok=True)
        write_summary_tsv(gene_summary_rows, str(gene_out / f"{gene_name}.tsv"))
        write_events_tsv(gene_events, str(gene_out / f"{gene_name}.events.tsv"))
        write_sites_tsv(gene_sites, str(gene_out / f"{gene_name}.sites.tsv"))

    if args.verbose:
        print(f"\nPipeline complete!")

    # Cleanup temp holder if used
    if temp_holder:
        temp_holder.cleanup()

    return 0


if __name__ == '__main__':
    sys.exit(main())
