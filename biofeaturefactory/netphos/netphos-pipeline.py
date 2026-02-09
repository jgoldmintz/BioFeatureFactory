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
NetPhos Pipeline

Synthesizes WT + mutant AA FASTAs, runs NetPhos via native APE, parses raw
output, and builds unified ensemble comparison tables (summary / events / sites
TSVs) with per-mutation deltas.

Modes:
  full-pipeline  – synthesize FASTAs → run NetPhos → parse → ensemble
  netphos-only   – run NetPhos on supplied FASTAs (no parsing)
  parse-only     – parse existing NetPhos output directories → ensemble
"""

import re
import sys
import argparse
import os
import csv
import subprocess
import tempfile
import shutil
import hashlib
import json
import platform
from pathlib import Path

from biofeaturefactory.utils.utility import (
    split_fasta_into_batches,
    combine_batch_outputs,
    discover_mapping_files,
    discover_fasta_files,
    load_validation_failures,
    load_wt_sequence_map,
    translate_orf_sequence,
    build_mutant_sequences_for_gene,
    synthesize_gene_fastas,
    extract_mutation_from_sequence_name,
    extract_gene_from_filename,
    read_fasta,
    resolve_output_base,
)


# ---------------------------------------------------------------------------
# Platform helpers
# ---------------------------------------------------------------------------

def is_linux_host():
    """Return True when running on a Linux kernel."""
    return platform.system().lower() == "linux"


def resolve_native_ape_path(user_path=None):
    """
    Resolve path to native APE binary.
    Checks: 1) user-provided path, 2) environment variables, 3) common locations.
    Returns absolute path if found and executable, else None.
    """
    candidates = []

    def _add(path):
        if path:
            candidates.append(os.path.expanduser(path))

    _add(user_path)
    _add(os.environ.get("NETPHOS_APE_PATH"))
    netphos_home = os.environ.get("NETPHOS_HOME")
    if netphos_home:
        _add(os.path.join(netphos_home, "ape-1.0", "ape"))
    netnglyc_home = os.environ.get("NETNGLYC_HOME")
    if netnglyc_home:
        _add(os.path.join(netnglyc_home, "ape-1.0", "ape"))

    home = Path.home()
    common_roots = [
        home / "ape-1.0" / "ape",
        home / "netphos" / "ape-1.0" / "ape",
        home / "netNglyc" / "ape-1.0" / "ape",
        Path("/opt/netphos/ape-1.0/ape"),
        Path("/opt/netnglyc/ape-1.0/ape"),
        Path("/usr/local/bin/ape"),
    ]
    for candidate in common_roots:
        _add(str(candidate))

    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)
    return None


# ---------------------------------------------------------------------------
# NetPhos raw output parser
# ---------------------------------------------------------------------------

def parse_netphos_file(input_file):
    """Parse NetPhos output file and extract phosphorylation predictions."""
    predictions = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Skip header lines
            if line.startswith('# Sequence') or line.startswith('# ---'):
                continue

            # Parse data lines starting with #
            if line.startswith('#'):
                clean_line = line[1:].strip()  # Remove leading #
                parts = clean_line.split()

                # Expect: seq_name pos aa context score kinase answer
                if len(parts) >= 6 and parts[2] in ['S', 'T', 'Y']:
                    try:
                        predictions.append({
                            'seq_name': parts[0],
                            'pos': int(parts[1]),
                            'amino_acid': parts[2],
                            'context': parts[3],
                            'score': float(parts[4]),
                            'kinase': parts[5],
                            'answer': parts[6] if len(parts) > 6 else '.',
                        })
                    except (ValueError, IndexError):
                        continue

    return predictions


# ---------------------------------------------------------------------------
# NetPhos execution
# ---------------------------------------------------------------------------

def _run_native_netphos(fasta_file, output_file, timeout=300, ape_bin=None):
    """Run NetPhos using native APE binary."""
    if not ape_bin:
        return False, "No APE binary path provided"

    try:
        ape_dir = os.path.dirname(ape_bin)
        result = subprocess.run(
            [ape_bin, "-m", "netphos", fasta_file],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=ape_dir
        )

        if result.returncode == 0:
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            return True, result.stdout
        else:
            error_msg = f"APE failed with return code {result.returncode}\n"
            error_msg += f"STDERR: {result.stderr}"
            return False, error_msg
    except subprocess.TimeoutExpired:
        return False, f"Native APE command timed out after {timeout} seconds"
    except Exception as e:
        return False, str(e)


def process_netphos_batched(fasta_file, output_file, batch_size=100, timeout=300, executor_fn=None, ape_bin=None):
    """Process large FASTA files using batching to prevent segmentation faults."""
    if executor_fn is None:
        executor_fn = _run_native_netphos

    try:
        batch_files = split_fasta_into_batches(fasta_file, batch_size)

        if not batch_files:
            print("No sequences found in FASTA file")
            return False

        print(f"Processing {len(batch_files)} batches...")
        batch_outputs = []

        for i, batch_file in enumerate(batch_files):
            batch_output = output_file.replace('.out', f'-batch-{i+1}.out')
            print(f"Processing batch {i+1}/{len(batch_files)}...")

            success, error = executor_fn(batch_file, batch_output, timeout, ape_bin)

            if success:
                batch_outputs.append(batch_output)
                print(f"Batch {i+1} completed: {batch_output}")
            else:
                print(f"Batch {i+1} failed: {error}")

        if batch_outputs:
            if len(batch_outputs) == 1:
                single_batch_output = batch_outputs[0]
                try:
                    shutil.move(single_batch_output, output_file)
                    print(f"Single batch completed: {output_file}")
                    return True
                except Exception as e:
                    print(f"Failed to move single batch output: {e}")
                    return False
            else:
                success = combine_batch_outputs(batch_outputs, output_file, format_type='netphos')
                if success:
                    print(f"Combined {len(batch_outputs)} batch outputs into {output_file}")
                    return True
                else:
                    print("Failed to combine batch outputs")
                    return False
        else:
            print("No successful batches to combine")
            return False

    except Exception as e:
        print(f"Batch processing failed: {e}")
        return False


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def get_file_cache_key(fasta_file):
    """Generate cache key based on file path and modification time."""
    try:
        file_path = str(fasta_file)
        stat = os.stat(file_path)
        cache_data = f"{file_path}:{stat.st_size}:{stat.st_mtime}"
        return hashlib.md5(cache_data.encode()).hexdigest()
    except Exception:
        return None


def get_cache_dir():
    """Get or create cache directory."""
    cache_dir = os.path.expanduser("~/.netphos_cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


def get_cached_result(fasta_file, cache_type="netphos"):
    """Check if cached result exists for file."""
    cache_key = get_file_cache_key(fasta_file)
    if not cache_key:
        return None, None

    cache_dir = get_cache_dir()
    cache_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.out")
    metadata_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.json")

    if os.path.exists(cache_file) and os.path.exists(metadata_file):
        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            return cache_file, metadata
        except Exception:
            for f in [cache_file, metadata_file]:
                if os.path.exists(f):
                    os.remove(f)

    return None, None


def save_to_cache(fasta_file, output_file, cache_type="netphos", metadata=None):
    """Save result to cache."""
    cache_key = get_file_cache_key(fasta_file)
    if not cache_key or not os.path.exists(output_file):
        return False

    try:
        cache_dir = get_cache_dir()
        cache_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.out")
        metadata_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.json")

        shutil.copy2(output_file, cache_file)

        cache_metadata = {
            "original_file": str(fasta_file),
            "cache_key": cache_key,
            "cache_type": cache_type,
            "cached_at": str(subprocess.run(["date"], capture_output=True, text=True).stdout.strip()),
            "file_size": os.path.getsize(fasta_file)
        }
        if metadata:
            cache_metadata.update(metadata)

        with open(metadata_file, 'w') as f:
            json.dump(cache_metadata, f, indent=2)

        return True
    except Exception as e:
        print(f"Warning: Failed to save cache: {e}")
        return False


def clear_cache():
    """Clear all cached NetPhos results."""
    cache_dir = get_cache_dir()
    try:
        if os.path.exists(cache_dir):
            cache_files = list(Path(cache_dir).glob("*"))
            for cache_file in cache_files:
                os.remove(cache_file)
            print(f"Cleared {len(cache_files)} cached files from {cache_dir}")
            return len(cache_files)
        else:
            print("No cache directory found")
            return 0
    except Exception as e:
        print(f"Error clearing cache: {e}")
        return 0


def count_fasta_sequences(fasta_file):
    """Count sequences in FASTA file."""
    count = 0
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        return 1
    return count


def run_netphos_with_fasta(fasta_file, output_file, batch_size=None, timeout=300, use_cache=True,
                           executor_fn=None, ape_bin=None):
    """Run NetPhos on FASTA file with intelligent processing strategy selection and caching."""
    if executor_fn is None:
        executor_fn = _run_native_netphos

    mode_label = "native"

    def _run(fasta, output, tmo):
        return executor_fn(fasta, output, tmo, ape_bin)

    if use_cache:
        cached_file, cached_metadata = get_cached_result(fasta_file, "netphos")
        if cached_file:
            print(f"Using cached NetPhos result for {fasta_file}")
            if cached_metadata:
                print(f"  Cached at: {cached_metadata.get('cached_at', 'Unknown')}")
            shutil.copy2(cached_file, output_file)
            return True

    if batch_size:
        print(f"Using batch processing (batch_size: {batch_size})...")
        result = process_netphos_batched(fasta_file, output_file, batch_size, timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        if result and use_cache:
            seq_count = count_fasta_sequences(fasta_file)
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": batch_size, "sequence_count": seq_count})
        return result

    seq_count = count_fasta_sequences(fasta_file)
    print(f"Processing {seq_count} sequence(s)...")

    if seq_count == 1:
        print(f"Single sequence detected - using single {mode_label} run...")
        success, error = _run(fasta_file, output_file, timeout)

        if success:
            print(f"NetPhos completed: {output_file}")
            if use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "single", "sequence_count": seq_count})
            return True
        else:
            print(f"Single run failed: {error}")
            return False

    elif seq_count <= 10:
        print(f"Small sequence set - attempting single {mode_label} run...")
        success, error = _run(fasta_file, output_file, timeout)

        if success:
            print(f"NetPhos completed: {output_file}")
            if use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "single", "sequence_count": seq_count})
            return True
        else:
            print(f"Single run failed: {error}")
            print("Falling back to batch processing for small sequence set...")
            result = process_netphos_batched(fasta_file, output_file, batch_size=10, timeout=timeout,
                                             executor_fn=executor_fn, ape_bin=ape_bin)
            if result and use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "batch_fallback", "batch_size": 10, "sequence_count": seq_count})
            return result

    elif seq_count <= 100:
        print("Medium sequence set - using batch processing...")
        result = process_netphos_batched(fasta_file, output_file, batch_size=25, timeout=timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        if result and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 25, "sequence_count": seq_count})
        return result

    else:
        print("Large sequence set - using batch processing...")
        result = process_netphos_batched(fasta_file, output_file, batch_size=50, timeout=timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        if result and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 50, "sequence_count": seq_count})
        return result


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------

def _classify_netphos_event(wt_score, mut_score, threshold, delta_threshold=0.05):
    """Classify a single (position, kinase) pair between WT and MUT.

    Returns (classification, classification_code, delta).
    """
    wt_above = wt_score is not None and wt_score >= threshold
    mut_above = mut_score is not None and mut_score >= threshold

    wt_val = wt_score if wt_score is not None else 0.0
    mut_val = mut_score if mut_score is not None else 0.0
    delta = mut_val - wt_val

    if wt_above and not mut_above:
        return "lost", -2, delta
    if not wt_above and mut_above:
        return "gained", 2, delta
    if wt_above and mut_above:
        if delta >= delta_threshold:
            return "strengthened", 1, delta
        if delta <= -delta_threshold:
            return "weakened", -1, delta
        return "stable", 0, delta
    return "subthreshold", -3, delta


# ---------------------------------------------------------------------------
# Ensemble builder
# ---------------------------------------------------------------------------

def _collect_output_files(directory):
    """Collect NetPhos output files from a directory."""
    output_files = []
    for ext in ['*.out', '*.txt']:
        output_files.extend(Path(directory).glob(ext))
    return sorted(output_files)


def _parse_directory_predictions(directory):
    """Parse all NetPhos output files in a directory, returning a dict keyed by seq_name."""
    all_preds = {}
    output_files = _collect_output_files(directory)
    for fpath in output_files:
        preds = parse_netphos_file(str(fpath))
        for pred in preds:
            sn = pred['seq_name']
            all_preds.setdefault(sn, []).append(pred)
    return all_preds


def build_netphos_ensemble(wt_preds_by_gene, mut_preds_by_mutation, mapping_lookup,
                           threshold=0.5, delta_threshold=0.05):
    """Build ensemble comparison tables from parsed WT and MUT predictions.

    Args:
        wt_preds_by_gene: dict  gene_name -> list[pred_dict]
        mut_preds_by_mutation: dict  (gene, nt_mutation) -> list[pred_dict]
        mapping_lookup: dict  gene -> mapping CSV path
        threshold: score threshold for YES/NO
        delta_threshold: minimum absolute delta for strengthened/weakened

    Returns:
        (summary_rows, events_rows, sites_rows)
    """
    summary_rows = []
    events_rows = []
    sites_rows = []

    for (gene, nt_mutation), mut_preds in mut_preds_by_mutation.items():
        pkey = f"{gene}-{nt_mutation}"
        wt_preds = wt_preds_by_gene.get(gene, [])

        # Build lookup: (position, kinase) -> best score
        wt_map = {}
        for p in wt_preds:
            key = (p['pos'], p['kinase'])
            if key not in wt_map or p['score'] > wt_map[key]['score']:
                wt_map[key] = p

        mut_map = {}
        for p in mut_preds:
            key = (p['pos'], p['kinase'])
            if key not in mut_map or p['score'] > mut_map[key]['score']:
                mut_map[key] = p

        all_keys = set(wt_map.keys()) | set(mut_map.keys())

        # --- sites rows (raw predictions tagged with allele) ---
        for p in wt_preds:
            sites_rows.append({
                'pkey': pkey,
                'Gene': gene,
                'allele': 'WT',
                'seq_name': p['seq_name'],
                'position': p['pos'],
                'amino_acid': p['amino_acid'],
                'context': p['context'],
                'score': p['score'],
                'kinase': p['kinase'],
                'answer': p['answer'],
            })
        for p in mut_preds:
            sites_rows.append({
                'pkey': pkey,
                'Gene': gene,
                'allele': 'MUT',
                'seq_name': p['seq_name'],
                'position': p['pos'],
                'amino_acid': p['amino_acid'],
                'context': p['context'],
                'score': p['score'],
                'kinase': p['kinase'],
                'answer': p['answer'],
            })

        # --- events rows ---
        mutation_events = []
        for pos, kinase in sorted(all_keys):
            wp = wt_map.get((pos, kinase))
            mp = mut_map.get((pos, kinase))

            wt_score = wp['score'] if wp else None
            mut_score = mp['score'] if mp else None
            wt_aa = wp['amino_acid'] if wp else '.'
            mut_aa = mp['amino_acid'] if mp else '.'
            wt_answer = wp['answer'] if wp else '.'
            mut_answer = mp['answer'] if mp else '.'

            classification, code, delta = _classify_netphos_event(
                wt_score, mut_score, threshold, delta_threshold)

            event = {
                'pkey': pkey,
                'Gene': gene,
                'position': pos,
                'amino_acid_wt': wt_aa,
                'amino_acid_mut': mut_aa,
                'kinase': kinase,
                'wt_score': wt_score if wt_score is not None else '',
                'mut_score': mut_score if mut_score is not None else '',
                'delta': round(delta, 6),
                'wt_answer': wt_answer,
                'mut_answer': mut_answer,
                'classification': classification,
                'classification_code': code,
            }
            events_rows.append(event)
            mutation_events.append(event)

        # --- summary row ---
        n_sites_wt = sum(1 for k in all_keys
                         if (w := wt_map.get(k)) and w['score'] >= threshold)
        n_sites_mut = sum(1 for k in all_keys
                          if (m := mut_map.get(k)) and m['score'] >= threshold)

        count_gained = sum(1 for e in mutation_events if e['classification'] == 'gained')
        count_lost = sum(1 for e in mutation_events if e['classification'] == 'lost')
        count_strengthened = sum(1 for e in mutation_events if e['classification'] == 'strengthened')
        count_weakened = sum(1 for e in mutation_events if e['classification'] == 'weakened')
        count_stable = sum(1 for e in mutation_events if e['classification'] == 'stable')

        abs_deltas = [abs(e['delta']) for e in mutation_events]
        max_abs_delta = max(abs_deltas) if abs_deltas else 0.0
        sum_abs_delta = sum(abs_deltas)

        non_stable = [e for e in mutation_events if e['classification'] not in ('stable', 'subthreshold')]
        n_kinases_affected = len({e['kinase'] for e in non_stable})

        top_event = max(mutation_events, key=lambda e: abs(e['delta'])) if mutation_events else None

        qc_flags = []
        if not wt_preds:
            qc_flags.append("missing_wt")
        if not mut_preds:
            qc_flags.append("missing_mut")
        if max_abs_delta == 0.0 and mutation_events:
            qc_flags.append("no_delta")

        summary_rows.append({
            'pkey': pkey,
            'Gene': gene,
            'n_sites_wt': n_sites_wt,
            'n_sites_mut': n_sites_mut,
            'count_gained': count_gained,
            'count_lost': count_lost,
            'count_strengthened': count_strengthened,
            'count_weakened': count_weakened,
            'count_stable': count_stable,
            'max_abs_delta': round(max_abs_delta, 6),
            'sum_abs_delta': round(sum_abs_delta, 6),
            'n_kinases_affected': n_kinases_affected,
            'top_event_type': top_event['classification'] if top_event else '',
            'top_event_delta': round(top_event['delta'], 6) if top_event else 0.0,
            'top_event_position': top_event['position'] if top_event else '',
            'top_event_kinase': top_event['kinase'] if top_event else '',
            'top_event_classification_code': top_event['classification_code'] if top_event else '',
            'qc_flags': '|'.join(qc_flags) if qc_flags else '',
        })

    return summary_rows, events_rows, sites_rows


def write_ensemble_outputs(output_base, summary_rows, events_rows, sites_rows):
    """Write the three ensemble TSV files."""
    summary_path = f"{output_base}.tsv"
    events_path = f"{output_base}.events.tsv"
    sites_path = f"{output_base}.sites.tsv"

    # Summary
    summary_fields = [
        'pkey', 'Gene', 'n_sites_wt', 'n_sites_mut',
        'count_gained', 'count_lost', 'count_strengthened', 'count_weakened', 'count_stable',
        'max_abs_delta', 'sum_abs_delta', 'n_kinases_affected',
        'top_event_type', 'top_event_delta', 'top_event_position', 'top_event_kinase',
        'top_event_classification_code', 'qc_flags',
    ]
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {len(summary_rows)} summary rows to {summary_path}")

    # Events
    events_fields = [
        'pkey', 'Gene', 'position', 'amino_acid_wt', 'amino_acid_mut',
        'kinase', 'wt_score', 'mut_score', 'delta',
        'wt_answer', 'mut_answer', 'classification', 'classification_code',
    ]
    with open(events_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=events_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(events_rows)
    print(f"Wrote {len(events_rows)} events to {events_path}")

    # Sites
    sites_fields = [
        'pkey', 'Gene', 'allele', 'seq_name', 'position',
        'amino_acid', 'context', 'score', 'kinase', 'answer',
    ]
    with open(sites_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=sites_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sites_rows)
    print(f"Wrote {len(sites_rows)} site predictions to {sites_path}")


# ---------------------------------------------------------------------------
# Mode runners
# ---------------------------------------------------------------------------

def _run_netphos_on_directory(fasta_dir, output_dir, args, executor_fn, ape_bin):
    """Run NetPhos on all FASTAs in fasta_dir, writing outputs to output_dir."""
    discovered = discover_fasta_files(str(fasta_dir))
    fasta_files = list(discovered.values())

    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        return []

    print(f"Found {len(fasta_files)} FASTA files in {fasta_dir}")
    outputs = []
    use_cache = not args.no_cache

    for fasta_path in fasta_files:
        input_basename = os.path.basename(fasta_path)
        for ext in ['.fasta', '.fa', '.fas', '.fna']:
            if input_basename.lower().endswith(ext):
                input_basename = input_basename[:-len(ext)]
                break
        netphos_output = os.path.join(str(output_dir), f"{input_basename}-netphos.out")

        seq_count = count_fasta_sequences(fasta_path)
        print(f"Processing {fasta_path} ({seq_count} sequences)...")
        success = run_netphos_with_fasta(fasta_path, netphos_output, args.batch_size, args.timeout, use_cache,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        if success:
            outputs.append(netphos_output)
            print(f"NetPhos completed: {netphos_output}")
        else:
            print(f"NetPhos failed for: {fasta_path}")

    return outputs


def _pair_predictions_with_mutations(wt_dir, mut_dir, mapping_lookup):
    """Parse WT and MUT output directories and pair predictions by gene/mutation.

    Returns (wt_preds_by_gene, mut_preds_by_mutation).
    """
    # Parse WT outputs - keyed by gene
    wt_preds_by_gene = {}
    wt_raw = _parse_directory_predictions(wt_dir)
    for seq_name, preds in wt_raw.items():
        gene, _ = extract_mutation_from_sequence_name(seq_name)
        gene = gene.upper().replace('_AA', '').replace('-WT', '').replace('_WT', '')
        gene_clean = extract_gene_from_filename(gene) or gene
        wt_preds_by_gene.setdefault(gene_clean.upper(), []).extend(preds)

    # Parse MUT outputs - keyed by (gene, nt_mutation)
    mut_preds_by_mutation = {}
    mut_raw = _parse_directory_predictions(mut_dir)
    for seq_name, preds in mut_raw.items():
        gene, mutation = extract_mutation_from_sequence_name(seq_name)
        gene = gene.upper().replace('_AA', '')
        gene_clean = extract_gene_from_filename(gene) or gene
        if mutation:
            mut_preds_by_mutation.setdefault((gene_clean.upper(), mutation), []).extend(preds)

    return wt_preds_by_gene, mut_preds_by_mutation


def run_full_pipeline_mode(args, executor_fn, ape_bin):
    """Synthesize FASTAs, run NetPhos, parse, build ensemble."""
    if not args.mapping_dir:
        print("ERROR: --mapping-dir is required for full-pipeline mode")
        return 1

    args.output = resolve_output_base(args.output, args.input, "netphos")

    wt_header = getattr(args, 'wt_header', 'ORF')
    verbose = getattr(args, 'verbose', False)
    keep = getattr(args, 'keep_intermediates', False)

    # Load WT sequences
    wt_sequences, temp_holder = load_wt_sequence_map(args.input, wt_header=wt_header)
    if not wt_sequences:
        print("ERROR: No WT sequences found")
        return 1

    if verbose:
        print(f"Loaded {len(wt_sequences)} WT sequences")

    # Discover mapping files
    mapping_lookup = discover_mapping_files(args.mapping_dir)
    if not mapping_lookup:
        print(f"ERROR: No mapping files found in {args.mapping_dir}")
        return 1

    if verbose:
        print(f"Found mappings for genes: {list(mapping_lookup.keys())}")

    # Load failure map
    failure_map = load_validation_failures(args.log) if args.log else {}

    # Synthesize FASTAs
    work_dir = tempfile.mkdtemp(prefix="netphos_pipeline_")
    seq_root = os.path.join(work_dir, "sequences")
    wt_dir, mut_dir, synth_summary = synthesize_gene_fastas(
        wt_sequences, mapping_lookup, seq_root,
        log_path=args.log, failure_map=failure_map,
    )

    total_mutants = sum(s['mutant_count'] for s in synth_summary)
    print(f"Synthesized FASTAs: {len(synth_summary)} genes, {total_mutants} mutants")

    # Run NetPhos on WT and MUT directories
    wt_output_dir = Path(work_dir) / "wt_outputs"
    mut_output_dir = Path(work_dir) / "mut_outputs"
    wt_output_dir.mkdir(parents=True, exist_ok=True)
    mut_output_dir.mkdir(parents=True, exist_ok=True)

    print("Running NetPhos on WT sequences...")
    wt_outputs = _run_netphos_on_directory(wt_dir, wt_output_dir, args, executor_fn, ape_bin)
    print("Running NetPhos on MUT sequences...")
    mut_outputs = _run_netphos_on_directory(mut_dir, mut_output_dir, args, executor_fn, ape_bin)

    if not wt_outputs and not mut_outputs:
        print("ERROR: NetPhos execution failed for all files")
        return 1

    print(f"NetPhos completed: {len(wt_outputs)} WT, {len(mut_outputs)} MUT output files")

    # Parse and pair predictions
    wt_preds_by_gene, mut_preds_by_mutation = _pair_predictions_with_mutations(
        str(wt_output_dir), str(mut_output_dir), mapping_lookup)

    # Build ensemble
    summary, events, sites = build_netphos_ensemble(
        wt_preds_by_gene, mut_preds_by_mutation, mapping_lookup,
        threshold=args.threshold, delta_threshold=0.05)

    # Strip .tsv from output if user included it, since write_ensemble_outputs adds it
    output_base = args.output
    if output_base.endswith('.tsv'):
        output_base = output_base[:-4]

    write_ensemble_outputs(output_base, summary, events, sites)

    # Cleanup
    if not keep:
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        print(f"Intermediate files preserved at: {work_dir}")

    if temp_holder:
        temp_holder.cleanup()

    return 0


def run_netphos_only_mode(args, executor_fn, ape_bin):
    """Run NetPhos on supplied FASTAs without parsing."""
    temp_output_dir = tempfile.mkdtemp(prefix="netphos_outputs_")

    if os.path.isfile(args.input):
        fasta_file = args.input
        netphos_output = os.path.join(temp_output_dir,
                                      os.path.basename(args.input).replace('.fasta', '-netphos.out'))

        print(f"Running NetPhos on {fasta_file}...")
        use_cache = not args.no_cache
        success = run_netphos_with_fasta(fasta_file, netphos_output, args.batch_size, args.timeout, use_cache,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        if not success:
            print("ERROR: NetPhos execution failed")
            return 1

        print(f"NetPhos completed: {netphos_output}")

    elif os.path.isdir(args.input):
        outputs = _run_netphos_on_directory(args.input, temp_output_dir, args, executor_fn, ape_bin)
        if not outputs:
            print("ERROR: NetPhos execution failed for all files")
            return 1
        print(f"NetPhos completed for {len(outputs)} files")

    else:
        print(f"ERROR: Input must be a FASTA file or directory: {args.input}")
        return 1

    print(f"NetPhos outputs in: {temp_output_dir}")
    print("Use --mode parse-only to parse these outputs.")
    return 0


def run_parse_mode(args):
    """Parse existing NetPhos output directories and build ensemble tables."""
    if not args.mapping_dir:
        print("ERROR: --mapping-dir is required for parse-only mode")
        return 1

    args.output = resolve_output_base(args.output, args.input, "netphos")

    if not os.path.exists(args.mapping_dir):
        print(f"ERROR: Mapping path not found: {args.mapping_dir}")
        return 1

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: Input path not found: {args.input}")
        return 1

    mapping_lookup = discover_mapping_files(args.mapping_dir)
    if not mapping_lookup:
        print(f"ERROR: No mapping files found in {args.mapping_dir}")
        return 1

    # Expect input directory with wt/ and mut/ subdirs, or a flat directory
    wt_dir = input_path / "wt"
    mut_dir = input_path / "mut"

    if wt_dir.is_dir() and mut_dir.is_dir():
        wt_preds_by_gene, mut_preds_by_mutation = _pair_predictions_with_mutations(
            str(wt_dir), str(mut_dir), mapping_lookup)
    else:
        # Flat directory: treat all outputs as from a single allele context
        print("Warning: No wt/ and mut/ subdirectories found. Attempting flat directory parsing.")
        all_preds = _parse_directory_predictions(str(input_path))

        wt_preds_by_gene = {}
        mut_preds_by_mutation = {}
        for seq_name, preds in all_preds.items():
            gene, mutation = extract_mutation_from_sequence_name(seq_name)
            gene = gene.upper().replace('_AA', '')
            gene_clean = extract_gene_from_filename(gene) or gene
            if mutation and mutation.upper() not in ('WT',):
                mut_preds_by_mutation.setdefault((gene_clean.upper(), mutation), []).extend(preds)
            else:
                wt_preds_by_gene.setdefault(gene_clean.upper(), []).extend(preds)

    summary, events, sites = build_netphos_ensemble(
        wt_preds_by_gene, mut_preds_by_mutation, mapping_lookup,
        threshold=args.threshold, delta_threshold=0.05)

    output_base = args.output
    if output_base.endswith('.tsv'):
        output_base = output_base[:-4]

    write_ensemble_outputs(output_base, summary, events, sites)
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='NetPhos pipeline: synthesize FASTAs, run NetPhos, build ensemble comparison tables',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input', nargs='?',
                        help='Input: WT FASTA file/directory (full-pipeline), FASTA file/directory (netphos-only), or NetPhos output directory (parse-only)')
    parser.add_argument('output', nargs='?',
                        help='Output base path (produces .tsv, .events.tsv, .sites.tsv)')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Phosphorylation score threshold (default: 0.5)')
    parser.add_argument('--yes-only', action='store_true',
                        help='Only include predictions marked as YES')
    parser.add_argument('--mapping-dir',
                        help='Directory or single CSV file containing mutation mapping(s)')
    parser.add_argument('--log',
                        help='Validation log file or directory to skip failed mutations')
    parser.add_argument('--wt-header', default='ORF',
                        help='FASTA header identifying WT sequence (default: ORF)')
    parser.add_argument('--keep-intermediates', action='store_true',
                        help='Keep temp directories for debugging')
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose output')

    # Mode selection
    parser.add_argument('--mode', choices=['full-pipeline', 'netphos-only', 'parse-only'], default='full-pipeline',
                        help='Processing mode (default: full-pipeline)')
    parser.add_argument('--batch-size', type=int,
                        help='Batch size for large FASTA files')
    parser.add_argument('--timeout', type=int, default=300,
                        help='Command timeout in seconds (default: 300)')

    # Cache options
    parser.add_argument('--no-cache', action='store_true',
                        help='Disable result caching')
    parser.add_argument('--clear-cache', action='store_true',
                        help='Clear all cached NetPhos results and exit')

    # Native execution options
    parser.add_argument('--native-ape-path',
                        help='Path to native APE binary')

    args = parser.parse_args()

    # Handle cache clearing
    if args.clear_cache:
        cleared_count = clear_cache()
        if cleared_count > 0:
            print(f"Cache cleared successfully ({cleared_count} files removed)")
        else:
            print("No cached files found")
        return 0

    # Validate required arguments
    if not args.input or not args.output:
        parser.error("input and output arguments are required (unless using --clear-cache)")

    if not os.path.exists(args.input):
        parser.error(f"Input path not found: {args.input}")

    # Validate threshold/yes-only
    if args.yes_only and args.threshold > 0.0:
        print("Warning: --yes-only specified, ignoring --threshold")
        args.threshold = 0.0

    # Resolve APE binary (needed for execution modes)
    native_ape = None
    executor_fn = None
    ape_bin = None

    if args.mode in ['full-pipeline', 'netphos-only']:
        native_ape = resolve_native_ape_path(getattr(args, 'native_ape_path', None))
        if not native_ape:
            parser.error(
                "Native APE binary not found. Provide --native-ape-path, or set "
                "NETPHOS_APE_PATH / NETPHOS_HOME environment variable."
            )
        executor_fn = _run_native_netphos
        ape_bin = native_ape
        if args.verbose:
            print(f"Execution mode: native APE ({native_ape})")

    # Dispatch to mode runner
    if args.mode == 'full-pipeline':
        return run_full_pipeline_mode(args, executor_fn, ape_bin)
    elif args.mode == 'netphos-only':
        return run_netphos_only_mode(args, executor_fn, ape_bin)
    elif args.mode == 'parse-only':
        return run_parse_mode(args)

    return 0


if __name__ == '__main__':
    exit(main())
