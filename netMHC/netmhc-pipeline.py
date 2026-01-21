#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023–2026  Jacob Goldmintz
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
import csv
import subprocess
import tempfile
import shutil
from pathlib import Path
import json
import hashlib
import sys
import time
import platform
from datetime import datetime
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional, Dict, List, Tuple
from Bio.Seq import Seq


# Import utility functions
sys.path.append(os.path.join(os.path.dirname(__file__), '../utils'))

from utility import (
    read_fasta,
    get_mutation_data_bioAccurate,
    write_fasta,
    split_fasta_into_batches,
    combine_batch_outputs,
    discover_mapping_files,
    discover_fasta_files,
    ExtractGeneFromFASTA,
    parse_predictions_with_mutation_filtering,
    load_validation_failures,
    should_skip_mutation,
    load_wt_sequences,
    trim_muts,
    get_mutant_aa,
    update_str,
    extract_gene_from_filename,
    extract_mutation_from_sequence_name,
)


def translate_orf_sequence(nt_sequence: str) -> str:
    """Translate a nucleotide ORF into an amino acid sequence, trimming trailing stops."""
    if not nt_sequence:
        return ""
    cleaned = nt_sequence.strip().upper().replace("U", "T")
    if not cleaned:
        return ""
    aa_seq = str(Seq(cleaned).translate(to_stop=False))
    return aa_seq.rstrip('*').strip()


def is_linux_host():
    """Return True when running on a Linux kernel."""
    return platform.system().lower() == "linux"


def resolve_native_netmhc_path(user_path=None, tool_version="netMHCpan"):
    """
    Resolve a usable native NetMHC executable when available.

    Search order:
    1. Explicit --native-netmhc-path value
    2. $NETMHC_PATH or $NETMHCPAN_PATH environment variable
    3. Common install locations

    Args:
        user_path: User-specified path to NetMHC binary
        tool_version: Which NetMHC tool to use (netMHCpan, netMHCII, netMHC)
    """
    candidates = []

    def _add(path):
        if path:
            candidates.append(os.path.expanduser(path))

    _add(user_path)
    _add(os.environ.get("NETMHC_PATH"))
    _add(os.environ.get("NETMHCPAN_PATH"))

    netmhc_home = os.environ.get("NETMHC_HOME")
    if netmhc_home:
        _add(os.path.join(netmhc_home, tool_version))

    home = Path.home()
    common_roots = [
        home / "netMHCpan-4.1" / "netMHCpan",
        home / "netMHCpan" / "netMHCpan",
        home / "netMHC" / tool_version,
        Path(f"/opt/netMHC/{tool_version}"),
        Path(f"/usr/local/bin/{tool_version}"),
    ]

    for candidate in common_roots:
        _add(str(candidate))

    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)

    return None


def build_netmhc_executor(args, parser):
    """
    Decide whether to use Docker or a native NetMHC installation.

    Returns a callable executor and a description string.
    """
    if getattr(args, "force_native", False) and getattr(args, "force_docker", False):
        parser.error("--force-native and --force-docker cannot be used together")

    native_netmhc = resolve_native_netmhc_path(
        getattr(args, "native_netmhc_path", None),
        getattr(args, "netmhc_tool", "netMHCpan")
    )
    linux_host = is_linux_host()
    use_native = False

    if getattr(args, "force_native", False):
        if not native_netmhc:
            parser.error("--force-native requires a valid --native-netmhc-path or NETMHC_PATH")
        use_native = True
    elif getattr(args, "force_docker", False):
        use_native = False
    elif native_netmhc and (linux_host or os.environ.get("NETMHC_ALLOW_NATIVE", "0") == "1"):
        use_native = True

    verbose_flag = getattr(args, "verbose", True)

    if use_native:
        if verbose_flag:
            print(f"Execution mode: native NetMHC ({native_netmhc})")
        def _runner(fasta_file, output_file, timeout, alleles):
            return _run_native_netmhc(fasta_file, output_file, timeout, native_netmhc, alleles)
        return _runner, f"native:{native_netmhc}"

    if verbose_flag:
        reason = "forced" if getattr(args, "force_docker", False) else (
            "native binary not found" if not native_netmhc else
            "non-Linux host" if not linux_host else "default")
        print(f"Execution mode: Docker NetMHC ({reason})")
    return _run_docker_netmhc, "docker"


def _run_docker_netmhc(fasta_file, output_file, timeout=600, alleles=None):
    """
    Execute NetMHC via Docker container.

    Args:
        fasta_file: Input FASTA file path
        output_file: Output file path
        timeout: Command timeout in seconds
        alleles: List of HLA alleles to predict (e.g., ['HLA-A*02:01', 'HLA-B*07:02'])

    Returns:
        tuple: (success, output_content, error_message)
    """
    if not alleles:
        alleles = ["HLA-A*02:01"]  # Default allele

    work_dir = tempfile.mkdtemp(prefix="netmhc_docker_")
    all_outputs = []

    try:
        # Copy FASTA to work directory
        docker_input = os.path.join(work_dir, "input.fasta")
        shutil.copy2(fasta_file, docker_input)

        # NetMHC processes one allele at a time, so run for each allele
        for allele in alleles:
            # Build Docker command
            # Format: docker run --rm --platform linux/386 -v /path:/data biofeaturefactory:latest netMHC -a ALLELE -f /data/input.fasta
            docker_cmd = [
                "docker", "run", "--rm",
                "--platform", "linux/386",
                "-v", f"{work_dir}:/data",
                "biofeaturefactory:latest",
                "netMHC",
                "-a", allele,
                "-f", "/data/input.fasta"
            ]

            result = subprocess.run(
                docker_cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                return False, result.stdout, result.stderr

            all_outputs.append(result.stdout)

        # Combine all outputs
        combined_output = "\n".join(all_outputs)

        # Write to output file
        with open(output_file, 'w') as f:
            f.write(combined_output)

        return True, combined_output, None

    except subprocess.TimeoutExpired:
        return False, "", f"Docker command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)
    finally:
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir, ignore_errors=True)


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
        alleles = ["HLA-A*02:01"]  # Default allele

    all_outputs = []

    try:
        # NetMHC processes one allele at a time, so run for each allele
        for allele in alleles:
            # Build native NetMHC command
            # Format: netMHC -a HLA-A*02:01 -f input.fasta
            cmd = [netmhc_path, "-a", allele, "-f", fasta_file]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=os.path.dirname(netmhc_path)
            )

            if result.returncode != 0:
                return False, result.stdout, result.stderr

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


def write_netmhc_tsv(predictions, output_file, include_pkey=False):
    """
    Write NetMHC predictions to TSV format.

    Args:
        predictions: List of prediction dictionaries
        output_file: Output TSV file path
        include_pkey: Whether to include pkey column
    """
    if not predictions:
        print(f"Warning: No predictions to write to {output_file}")
        with open(output_file, 'w') as f:
            f.write("# No predictions generated\n")
        return

    fieldnames = ['Gene', 'pos', 'mhc_allele', 'peptide', 'core', 'score', 'rank', 'bind_level']
    if include_pkey:
        fieldnames.insert(0, 'pkey')

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(predictions)

    print(f"Wrote {len(predictions)} NetMHC predictions to {output_file}")


def build_mutant_sequences_for_gene(
    gene_name: str,
    nt_sequence: str,
    aa_sequence: str,
    mapping_file: Optional[str],
    log_path: Optional[str],
    failure_map: Optional[dict],
):
    """
    Build mutant amino acid sequences from nucleotide mutations.

    Returns a dict of {header: sequence} for all mutants of a given gene.
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
    try:
        with open(mapping_file, 'r') as handle:
            reader = csv.DictReader(handle)
            mutant_keys = ['mutant', 'mutation', 'nt_mutation', 'ntmutant']
            aa_keys = ['aamutant', 'aa_mutation', 'amino_acid_mutation', 'protein_mutation']

            for row in reader:
                mutant_id = ""
                for key in mutant_keys:
                    if key in row and row[key]:
                        mutant_id = row[key].strip()
                        break
                if not mutant_id:
                    continue

                mutant_clean = mutant_id.replace(" ", "")
                if allowed_mutations and mutant_clean.upper() not in allowed_mutations:
                    continue
                if should_skip_mutation(gene_name, mutant_clean, failure_map):
                    continue

                aa_string = ""
                for key in aa_keys:
                    if key in row and row[key]:
                        aa_string = row[key].strip()
                        break

                pos = None
                wt_aa = mut_aa = None
                if aa_string:
                    pos, nts = get_mutation_data_bioAccurate(aa_string)
                    if pos is not None and nts:
                        wt_aa, mut_aa = nts

                if pos is None or not wt_aa or not mut_aa:
                    # Try to infer from nucleotide mutation
                    nt_info = get_mutation_data_bioAccurate(mutant_clean)
                    if nt_info[0] is not None:
                        aa_info = get_mutant_aa(nt_info, nt_sequence)
                        if aa_info:
                            (pos, (wt_aa, mut_aa)), _ = aa_info

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
        print(f"Warning: Failed to synthesize mutants for {gene_name} ({mapping_file}): {exc}")
        return {}

    return mutant_sequences


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

    # Build lookup maps: (peptide, allele, pos) -> prediction
    wt_map = {}
    for pred in wt_preds:
        key = (pred['peptide'], pred['mhc_allele'], pred['pos'])
        wt_map[key] = pred

    mut_map = {}
    for pred in mut_preds:
        key = (pred['peptide'], pred['mhc_allele'], pred['pos'])
        mut_map[key] = pred

    # Find all unique peptide/allele/position combinations
    all_keys = set(wt_map.keys()) | set(mut_map.keys())

    for key in all_keys:
        peptide, allele, pos = key
        wt_pred = wt_map.get(key)
        mut_pred = mut_map.get(key)

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

        # Classify event
        classification = "stable"
        classification_code = 0

        if wt_rank > 2.0 and mut_rank <= threshold:
            classification = "gained"
            classification_code = 2
        elif wt_rank <= threshold and mut_rank > 2.0:
            classification = "lost"
            classification_code = -2
        elif wt_rank <= 2.0 and mut_rank <= 2.0:
            # Both are binders (weak or strong)
            if delta_rank < -5:
                classification = "strengthened"
                classification_code = 1
            elif delta_rank > 5:
                classification = "weakened"
                classification_code = -1
            else:
                classification = "stable"
                classification_code = 0

        event = {
            'Gene': gene_name,
            'pkey': f"{gene_name}-{mutation}",
            'mutation': mutation,
            'peptide': peptide,
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
        'top_event_peptide': top_event['peptide'] if top_event else '',
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
        'top_event_type', 'top_event_allele', 'top_event_peptide', 'top_event_delta_rank',
        'qc_flags'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote {len(summary_rows)} summary rows to {output_file}")


def write_events_tsv(events, output_file):
    """Write events TSV file."""
    if not events:
        print(f"Warning: No event data to write")
        return

    fieldnames = [
        'pkey', 'Gene', 'mutation', 'peptide', 'pos', 'mhc_allele',
        'wt_rank', 'mut_rank', 'delta_rank',
        'wt_affinity', 'mut_affinity', 'delta_affinity',
        'bind_level_wt', 'bind_level_mut',
        'classification', 'classification_code'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(events)

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

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sites)

    print(f"Wrote {len(sites)} site predictions to {output_file}")


# ============================================================================
# Caching Functions
# ============================================================================

def get_file_cache_key(fasta_file, alleles=None):
    """Generate cache key based on file content and alleles."""
    try:
        file_path = Path(fasta_file).resolve()
        stat = file_path.stat()
        # Include alleles in cache key since different alleles = different results
        allele_str = ",".join(sorted(alleles)) if alleles else "default"
        cache_data = f"{file_path}:{stat.st_size}:{stat.st_mtime}:{allele_str}"
        return hashlib.md5(cache_data.encode()).hexdigest()
    except Exception:
        return None


def get_cache_dir(custom_dir=None):
    """Get or create cache directory."""
    if custom_dir:
        cache_dir = Path(custom_dir)
    else:
        cache_dir = Path.home() / ".netmhc_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def get_cached_result(fasta_file, alleles=None, cache_dir=None):
    """Check if cached result exists for file and alleles."""
    cache_key = get_file_cache_key(fasta_file, alleles)
    if not cache_key:
        return None, None

    cache_dir = get_cache_dir(cache_dir)
    cache_file = cache_dir / f"{cache_key}_netmhc.out"
    metadata_file = cache_dir / f"{cache_key}_netmhc.json"

    if cache_file.exists() and metadata_file.exists():
        try:
            with open(metadata_file) as f:
                metadata = json.load(f)
            return str(cache_file), metadata
        except (json.JSONDecodeError, IOError):
            # Remove corrupted cache files
            for f in [cache_file, metadata_file]:
                if f.exists():
                    f.unlink()
    return None, None


def save_to_cache(fasta_file, output_file, alleles=None, cache_dir=None, metadata=None):
    """Save result to cache."""
    cache_key = get_file_cache_key(fasta_file, alleles)
    if not cache_key or not os.path.exists(output_file):
        return

    try:
        cache_dir = get_cache_dir(cache_dir)
        cache_file = cache_dir / f"{cache_key}_netmhc.out"
        metadata_file = cache_dir / f"{cache_key}_netmhc.json"

        # Copy output to cache
        shutil.copy2(output_file, cache_file)

        # Save metadata
        cache_metadata = {
            "source_file": str(fasta_file),
            "cache_key": cache_key,
            "alleles": alleles,
            "cached_at": datetime.now().isoformat(),
        }
        if metadata:
            cache_metadata.update(metadata)

        with open(metadata_file, 'w') as f:
            json.dump(cache_metadata, f, indent=2)

    except Exception as e:
        print(f"Warning: Failed to save cache: {e}")


def clear_cache(cache_dir=None):
    """Clear all cached NetMHC results."""
    cache_dir = get_cache_dir(cache_dir)
    count = 0
    for f in cache_dir.glob("*_netmhc.*"):
        f.unlink()
        count += 1
    print(f"Cleared {count} cached files from {cache_dir}")


def run_netmhc_with_cache(fasta_file, output_file, alleles=None, timeout=600, executor=None, use_cache=True, cache_dir=None):
    """
    Run NetMHC with caching support.

    Args:
        fasta_file: Input FASTA file
        output_file: Output file path
        alleles: List of HLA alleles
        timeout: Command timeout
        executor: NetMHC executor function
        use_cache: Whether to use cached results
        cache_dir: Custom cache directory

    Returns:
        bool: True if successful
    """
    if executor is None:
        executor = _run_docker_netmhc

    # Check cache first
    if use_cache:
        cached_file, cached_metadata = get_cached_result(fasta_file, alleles, cache_dir)
        if cached_file:
            print(f"Using cached result from {cached_metadata.get('cached_at', 'unknown')}")
            shutil.copy2(cached_file, output_file)
            return True

    # Run NetMHC
    success, output, error = executor(fasta_file, output_file, timeout, alleles)

    if not success:
        print(f"NetMHC execution failed: {error}")
        return False

    # Save to cache
    if use_cache:
        save_to_cache(fasta_file, output_file, alleles, cache_dir)

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="NetMHC pipeline for MHC binding prediction with WT/mutant comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('input', nargs='?',
                       help='Input: WT FASTA file/directory (full-pipeline mode) or NetMHC output directory (parse mode)')
    parser.add_argument('output', nargs='?',
                       help='Output TSV file')

    # Mode selection
    parser.add_argument('--mode', choices=['full-pipeline', 'netmhc-only', 'parse-only'],
                       default='full-pipeline',
                       help='Processing mode: full pipeline (NetMHC + parse), NetMHC-only (no parsing), or parse-only')

    # MHC-specific options
    parser.add_argument('--alleles', nargs='+',
                       help='HLA alleles to predict (e.g., HLA-A*02:01 HLA-B*07:02). If not specified, uses default set.')
    parser.add_argument('--netmhc-tool', choices=['netMHCpan', 'netMHC', 'netMHCII'],
                       default='netMHCpan',
                       help='Which NetMHC tool to use (default: netMHCpan)')
    parser.add_argument('--peptide-lengths', nargs='+', type=int,
                       help='Peptide lengths to predict (e.g., 8 9 10 11). Default depends on tool.')

    # Execution backend
    parser.add_argument('--native-netmhc-path',
                       help='Path to native NetMHC executable')
    parser.add_argument('--force-native', action='store_true',
                       help='Force native NetMHC execution')
    parser.add_argument('--force-docker', action='store_true',
                       help='Force Docker NetMHC execution')
    parser.add_argument('--docker-image', default='biofeaturefactory:latest',
                       help='Docker image name (default: biofeaturefactory:latest)')

    # Processing options
    parser.add_argument('--mapping-dir',
                       help='Directory containing mutation mapping CSV files (REQUIRED for full-pipeline mode)')
    parser.add_argument('--log',
                       help='Validation log file to skip failed mutations')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Binding rank threshold for strong binders (default: 0.5)')
    parser.add_argument('--batch-size', type=int, default=100,
                       help='Batch size for large FASTA files')
    parser.add_argument('--timeout', type=int, default=600,
                       help='Command timeout in seconds (default: 600)')

    # Cache and output options
    parser.add_argument('--cache-dir',
                       help='Custom cache directory for NetMHC results')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable result caching')
    parser.add_argument('--keep-intermediates', action='store_true',
                       help='Keep intermediate files for debugging')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Validate arguments
    if not args.input or not args.output:
        parser.error("input and output arguments are required")

    if args.mode == 'full-pipeline' and not args.mapping_dir:
        parser.error("--mapping-dir is REQUIRED for full-pipeline mode")

    # Build NetMHC executor
    executor, exec_desc = build_netmhc_executor(args, parser)

    # Load validation failures if provided
    failure_map = None
    if args.log:
        failure_map = load_validation_failures(args.log)
        if args.verbose:
            print(f"Loaded validation failures from {args.log}")

    # Initialize result collections
    all_summary_rows = []
    all_events = []
    all_sites = []

    if args.mode == 'full-pipeline':
        # Discover FASTA files
        fasta_files = discover_fasta_files(args.input)
        if not fasta_files:
            print(f"Error: No FASTA files found in {args.input}")
            return 1

        if args.verbose:
            print(f"Found {len(fasta_files)} FASTA files to process")

        # Discover mapping files
        mapping_files = discover_mapping_files(args.mapping_dir) if args.mapping_dir else {}

        # Process each gene
        for gene_name, fasta_path in fasta_files.items():
            # gene_name already extracted by discover_fasta_files
            if args.verbose:
                print(f"\nProcessing gene: {gene_name}")

            # Load WT sequence
            wt_sequences = read_fasta(fasta_path)
            if not wt_sequences:
                print(f"Warning: No sequences in {fasta_path}, skipping")
                continue
            if len(wt_sequences) > 1:
                if 'ORF' in wt_sequences:
                    wt_header, wt_nt_seq = 'ORF', wt_sequences['ORF']
                else:
                    wt_header, wt_nt_seq = next(iter(wt_sequences.items()))
            else:
                wt_header, wt_nt_seq = next(iter(wt_sequences.items()))


            # Translate to amino acids
            wt_aa_seq = translate_orf_sequence(wt_nt_seq)
            if not wt_aa_seq:
                print(f"Warning: Could not translate {gene_name}, skipping")
                continue

            # Build mutant sequences
            mapping_file = mapping_files.get(gene_name)
            mutant_seqs = build_mutant_sequences_for_gene(
                gene_name, wt_nt_seq, wt_aa_seq, mapping_file, args.log, failure_map
            )

            if args.verbose:
                print(f"  Generated {len(mutant_seqs)} mutant sequences")

            # Create temp directory for this gene
            gene_workdir = tempfile.mkdtemp(prefix=f"netmhc_{gene_name}_")

            try:
                # Write WT FASTA
                wt_fasta = os.path.join(gene_workdir, f"{gene_name}_wt.fasta")
                write_fasta(Path(wt_fasta), {f"{gene_name}_WT": wt_aa_seq})

                # Check if batching is needed for WT
                wt_predictions = []
                if len(wt_aa_seq) > args.batch_size and args.batch_size > 0:
                    # Split into batches
                    if args.verbose:
                        print(f"  Batching WT sequence ({len(wt_aa_seq)} AA) into chunks of {args.batch_size}")

                    batch_dir = os.path.join(gene_workdir, "wt_batches")
                    os.makedirs(batch_dir, exist_ok=True)

                    batch_files = split_fasta_into_batches(wt_fasta, batch_dir, args.batch_size)

                    # Run NetMHC on each batch
                    for batch_file in batch_files:
                        batch_output = batch_file.replace('.fasta', '.out')
                        success, _, error = executor(batch_file, batch_output, args.timeout, args.alleles)

                        if not success:
                            print(f"Warning: NetMHC failed for {gene_name} WT batch {batch_file}: {error}")
                            continue

                        # Parse batch predictions
                        batch_preds = parse_netmhc_output(batch_output)
                        wt_predictions.extend(batch_preds)
                else:
                    # Run NetMHC on full WT sequence
                    wt_output = os.path.join(gene_workdir, f"{gene_name}_wt.out")
                    success, _, error = executor(wt_fasta, wt_output, args.timeout, args.alleles)

                    if not success:
                        print(f"Warning: NetMHC failed for {gene_name} WT: {error}")
                        continue

                    # Parse WT predictions
                    wt_predictions = parse_netmhc_output(wt_output)
                if args.verbose:
                    print(f"  WT: {len(wt_predictions)} predictions")

                # Add to sites output (all WT predictions)
                for pred in wt_predictions:
                    all_sites.append({
                        'Gene': gene_name,
                        'pkey': f"{gene_name}-WT",
                        'sequence_type': 'wt',
                        **pred
                    })

                # Process each mutant
                for mut_header, mut_aa_seq in mutant_seqs.items():
                    mutation = mut_header.split('-', 1)[1] if '-' in mut_header else mut_header

                    # Write mutant FASTA
                    mut_fasta = os.path.join(gene_workdir, f"{gene_name}_{mutation}.fasta")
                    write_fasta(Path(mut_fasta), {mut_header: mut_aa_seq})

                    # Check if batching is needed for mutant
                    mut_predictions = []
                    if len(mut_aa_seq) > args.batch_size and args.batch_size > 0:
                        # Split into batches
                        if args.verbose:
                            print(f"  Batching {mutation} sequence ({len(mut_aa_seq)} AA)")

                        batch_dir = os.path.join(gene_workdir, f"{mutation}_batches")
                        os.makedirs(batch_dir, exist_ok=True)

                        batch_files = split_fasta_into_batches(mut_fasta, batch_dir, args.batch_size)

                        # Run NetMHC on each batch
                        for batch_file in batch_files:
                            batch_output = batch_file.replace('.fasta', '.out')
                            success, _, error = executor(batch_file, batch_output, args.timeout, args.alleles)

                            if not success:
                                print(f"Warning: NetMHC failed for {gene_name} {mutation} batch {batch_file}: {error}")
                                continue

                            # Parse batch predictions
                            batch_preds = parse_netmhc_output(batch_output)
                            mut_predictions.extend(batch_preds)
                    else:
                        # Run NetMHC on full mutant sequence
                        mut_output = os.path.join(gene_workdir, f"{gene_name}_{mutation}.out")
                        success, _, error = executor(mut_fasta, mut_output, args.timeout, args.alleles)

                        if not success:
                            print(f"Warning: NetMHC failed for {gene_name} {mutation}: {error}")
                            continue

                        # Parse mutant predictions
                        mut_predictions = parse_netmhc_output(mut_output)
                    if args.verbose:
                        print(f"  {mutation}: {len(mut_predictions)} predictions")

                    # Add to sites output (all mutant predictions)
                    for pred in mut_predictions:
                        all_sites.append({
                            'Gene': gene_name,
                            'pkey': f"{gene_name}-{mutation}",
                            'sequence_type': 'mut',
                            **pred
                        })

                    # Compare WT vs MUT and classify epitope changes
                    epitope_events = compare_wt_mut_predictions(
                        gene_name, mutation, wt_predictions, mut_predictions, args.threshold
                    )

                    # Add events to collection
                    all_events.extend(epitope_events)

                    # Generate summary for this mutation
                    summary = summarize_epitope_changes(gene_name, mutation, epitope_events)
                    all_summary_rows.append(summary)

            finally:
                # Clean up temp directory
                if not args.keep_intermediates and os.path.exists(gene_workdir):
                    shutil.rmtree(gene_workdir, ignore_errors=True)

        # Write output files
        write_summary_tsv(all_summary_rows, args.output)
        write_events_tsv(all_events, args.output.replace('.tsv', '.events.tsv'))
        write_sites_tsv(all_sites, args.output.replace('.tsv', '.sites.tsv'))

        if args.verbose:
            print(f"\n Pipeline complete!")
            print(f"  Summary: {args.output}")
            print(f"  Events: {args.output.replace('.tsv', '.events.tsv')}")
            print(f"  Sites: {args.output.replace('.tsv', '.sites.tsv')}")

    else:
        print(f"Mode {args.mode} not yet implemented")
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
