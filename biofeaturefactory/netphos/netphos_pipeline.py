#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# Licensed under the PolyForm Noncommercial License 1.0.0.
# You may use this software for any purpose other than commercial use.
# See LICENSE, or https://polyformproject.org/licenses/noncommercial/1.0.0/
#
# Required Notice: Copyright (C) 2023-2026 Jacob Goldmintz
# (https://github.com/jgoldmintz)

"""
NetPhos Pipeline

Synthesizes WT + mutant AA FASTAs, runs NetPhos via native APE, parses raw
output, and builds unified ensemble comparison tables (summary / events / sites
TSVs) with per-mutation deltas.

"""

import re
import sys
import argparse
import os
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
    resolve_output_base,
    write_tsv,
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
    """Process large FASTA files in batches (APE caps at ~2000 seqs/run — see
    software/ape-1.0/ape 'maxn').

    Returns (ok, complete):
      ok       -- True if at least one batch produced output (combined result usable)
      complete -- True only if EVERY batch succeeded (no mutants dropped)

    Each successful batch's output is content-cached, so a rerun reuses the good
    batches (APE skipped) and re-runs only the failed batch(es). Mutants in failed
    batches are written to <output_file>.dropped and printed — a partial run is
    never reported as a silent success.
    """
    if executor_fn is None:
        executor_fn = _run_native_netphos

    try:
        batch_files = split_fasta_into_batches(fasta_file, batch_size)

        if not batch_files:
            print("No sequences found in FASTA file")
            return False, False

        print(f"Processing {len(batch_files)} batches...")
        batch_outputs = []
        failed_batches = []

        for i, batch_file in enumerate(batch_files):
            batch_output = output_file.replace('.out', f'-batch-{i+1}.out')

            cached = _get_cached_batch(batch_file)
            if cached:
                shutil.copy2(cached, batch_output)
                batch_outputs.append(batch_output)
                print(f"Batch {i+1}/{len(batch_files)}: cache hit (APE skipped)")
                continue

            print(f"Processing batch {i+1}/{len(batch_files)}...")
            success, error = executor_fn(batch_file, batch_output, timeout, ape_bin)

            if success:
                batch_outputs.append(batch_output)
                _save_batch_cache(batch_file, batch_output)
                print(f"Batch {i+1} completed: {batch_output}")
            else:
                failed_batches.append(batch_file)
                print(f"Batch {i+1} FAILED: {error}")

        # Surface dropped mutants (the sequence headers of the failed batches).
        dropped = []
        for bf in failed_batches:
            try:
                with open(bf) as fh:
                    dropped.extend(
                        line[1:].strip().split()[0]
                        for line in fh if line.startswith('>') and line[1:].strip()
                    )
            except Exception:
                pass
        if dropped:
            drop_path = output_file + '.dropped'
            try:
                with open(drop_path, 'w') as fh:
                    fh.write('\n'.join(dropped) + '\n')
            except Exception:
                drop_path = '(unwritten)'
            print(f"WARNING: {len(dropped)} mutant(s) DROPPED from {len(failed_batches)} "
                  f"failed batch(es); listed in {drop_path}. Re-run to retry — successful "
                  f"batches are cached, so only the failed batch(es) re-run.")

        complete = not failed_batches

        if not batch_outputs:
            print("No successful batches to combine")
            return False, False

        if len(batch_outputs) == 1:
            try:
                shutil.move(batch_outputs[0], output_file)
                print(f"Single batch completed: {output_file}")
            except Exception as e:
                print(f"Failed to move single batch output: {e}")
                return False, False
        else:
            if not combine_batch_outputs(batch_outputs, output_file, format_type='netphos'):
                print("Failed to combine batch outputs")
                return False, False
            print(f"Combined {len(batch_outputs)} batch output(s) into {output_file}"
                  + ("" if complete else f" ({len(dropped)} mutant(s) dropped)"))

        return True, complete

    except Exception as e:
        print(f"Batch processing failed: {e}")
        return False, False


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


def _batch_content_cache_key(fasta_file):
    """Content-based md5 of a FASTA file (independent of path/mtime), so an
    identical batch reuses its cached result across reruns."""
    try:
        with open(fasta_file, 'rb') as fh:
            return hashlib.md5(fh.read()).hexdigest()
    except Exception:
        return None


def _get_cached_batch(fasta_file):
    """Return the cached .out path for a batch's content, or None."""
    key = _batch_content_cache_key(fasta_file)
    if not key:
        return None
    cache_file = os.path.join(get_cache_dir(), f"{key}_netphos_batch.out")
    return cache_file if os.path.exists(cache_file) else None


def _save_batch_cache(fasta_file, output_file):
    """Cache a successful batch's output keyed by the batch's content."""
    key = _batch_content_cache_key(fasta_file)
    if not key or not os.path.exists(output_file):
        return
    try:
        shutil.copy2(output_file, os.path.join(get_cache_dir(), f"{key}_netphos_batch.out"))
    except Exception as e:
        print(f"Warning: failed to cache batch: {e}")


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
        ok, complete = process_netphos_batched(fasta_file, output_file, batch_size, timeout,
                                               executor_fn=executor_fn, ape_bin=ape_bin)
        # Cache the whole-gene result only when COMPLETE — never cache a partial
        # (dropped-batch) output, or reruns would serve the incomplete result.
        if ok and complete and use_cache:
            seq_count = count_fasta_sequences(fasta_file)
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": batch_size, "sequence_count": seq_count})
        return ok

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
            ok, complete = process_netphos_batched(fasta_file, output_file, batch_size=10, timeout=timeout,
                                                   executor_fn=executor_fn, ape_bin=ape_bin)
            if ok and complete and use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "batch_fallback", "batch_size": 10, "sequence_count": seq_count})
            return ok

    elif seq_count <= 100:
        print("Medium sequence set - using batch processing...")
        ok, complete = process_netphos_batched(fasta_file, output_file, batch_size=25, timeout=timeout,
                                               executor_fn=executor_fn, ape_bin=ape_bin)
        if ok and complete and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 25, "sequence_count": seq_count})
        return ok

    else:
        print("Large sequence set - using batch processing...")
        ok, complete = process_netphos_batched(fasta_file, output_file, batch_size=50, timeout=timeout,
                                               executor_fn=executor_fn, ape_bin=ape_bin)
        if ok and complete and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 50, "sequence_count": seq_count})
        return ok


# ---------------------------------------------------------------------------
# Classification logic
# ---------------------------------------------------------------------------

def _classify_netphos_event(wt_score, mut_score, wt_above, mut_above, delta_threshold=0.05):
    """Classify a single (position, kinase) pair between WT and MUT.

    `wt_above`/`mut_above` are the site-membership booleans decided by the caller:
    score >= threshold, or (under --yes-only) NetPhos answer == 'YES'.
    Returns (classification, classification_code, delta).
    """

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
    """Collect NetPhos output files from a directory, EXCLUDING per-batch
    intermediates (`*-batch-N.out`). Those are already folded into the combined
    `{GENE}-netphos.out`, so including them would parse every mutant twice and
    double-count rows in sites.tsv."""
    output_files = []
    for ext in ['*.out', '*.txt']:
        output_files.extend(
            p for p in Path(directory).glob(ext) if '-batch-' not in p.name
        )
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
                           threshold=0.5, delta_threshold=0.05, yes_only=False,
                           expected_mutations=None):
    """Build ensemble comparison tables from parsed WT and MUT predictions.

    Args:
        wt_preds_by_gene: dict  gene_name -> list[pred_dict]
        mut_preds_by_mutation: dict  (gene, nt_mutation) -> list[pred_dict]
        mapping_lookup: dict  gene -> mapping CSV path
        threshold: score threshold for YES/NO
        delta_threshold: minimum absolute delta for strengthened/weakened
        expected_mutations: optional iterable of (gene, nt_mutation) that were
            submitted to NetPhos. Any submitted mutation with no parsed
            predictions (dropped batch, failed run) gets a zeroed summary row
            flagged missing_mut instead of vanishing from the output.

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

            if yes_only:
                wt_above = (wt_answer == 'YES')
                mut_above = (mut_answer == 'YES')
            else:
                wt_above = wt_score is not None and wt_score >= threshold
                mut_above = mut_score is not None and mut_score >= threshold
            classification, code, delta = _classify_netphos_event(
                wt_score, mut_score, wt_above, mut_above, delta_threshold)

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
        if yes_only:
            n_sites_wt = sum(1 for k in all_keys
                             if (w := wt_map.get(k)) and w['answer'] == 'YES')
            n_sites_mut = sum(1 for k in all_keys
                              if (m := mut_map.get(k)) and m['answer'] == 'YES')
        else:
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

    # Submitted mutations with no parsed predictions never enter the loop above,
    # so a dropped batch would silently shrink the table. Emit a flagged row.
    seen_keys = set(mut_preds_by_mutation.keys())
    for key in (expected_mutations or []):
        gene, nt_mutation = tuple(key)
        if (gene, nt_mutation) in seen_keys:
            continue
        seen_keys.add((gene, nt_mutation))

        wt_preds = wt_preds_by_gene.get(gene, [])
        wt_map = {}
        for p in wt_preds:
            wkey = (p['pos'], p['kinase'])
            if wkey not in wt_map or p['score'] > wt_map[wkey]['score']:
                wt_map[wkey] = p

        if yes_only:
            n_sites_wt = sum(1 for p in wt_map.values() if p['answer'] == 'YES')
        else:
            n_sites_wt = sum(1 for p in wt_map.values() if p['score'] >= threshold)

        qc_flags = []
        if not wt_preds:
            qc_flags.append("missing_wt")
        qc_flags.append("missing_mut")

        summary_rows.append({
            'pkey': f"{gene}-{nt_mutation}",
            'Gene': gene,
            'n_sites_wt': n_sites_wt,
            'n_sites_mut': 0,
            'count_gained': 0,
            'count_lost': 0,
            'count_strengthened': 0,
            'count_weakened': 0,
            'count_stable': 0,
            'max_abs_delta': 0.0,
            'sum_abs_delta': 0.0,
            'n_kinases_affected': 0,
            'top_event_type': '',
            'top_event_delta': 0.0,
            'top_event_position': '',
            'top_event_kinase': '',
            'top_event_classification_code': '',
            'qc_flags': '|'.join(qc_flags),
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
    write_tsv(summary_rows, summary_path, summary_fields, extrasaction='ignore')
    print(f"Wrote {len(summary_rows)} summary rows to {summary_path}")

    # Events
    events_fields = [
        'pkey', 'Gene', 'position', 'amino_acid_wt', 'amino_acid_mut',
        'kinase', 'wt_score', 'mut_score', 'delta',
        'wt_answer', 'mut_answer', 'classification', 'classification_code',
    ]
    write_tsv(events_rows, events_path, events_fields, extrasaction='ignore')
    print(f"Wrote {len(events_rows)} events to {events_path}")

    # Sites
    sites_fields = [
        'pkey', 'Gene', 'allele', 'seq_name', 'position',
        'amino_acid', 'context', 'score', 'kinase', 'answer',
    ]
    write_tsv(sites_rows, sites_path, sites_fields, extrasaction='ignore')
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


def _expected_mutations_from_synthesis(synth_summary):
    """Recover the (gene, mutation) keys submitted to NetPhos from the synthesized
    mutant FASTAs. Keys are normalized exactly as _pair_predictions_with_mutations
    normalizes parsed output names, so the two sets are comparable.
    """
    expected = []
    for entry in synth_summary:
        mut_path = entry.get('mut_path')
        if not mut_path:
            continue
        gene = str(entry.get('gene', '')).upper()
        if not gene:
            continue
        gene_key = (extract_gene_from_filename(gene) or gene).upper()
        prefix = f"{gene}-"
        try:
            with open(mut_path) as fh:
                for line in fh:
                    if not line.startswith('>'):
                        continue
                    header = line[1:].strip()
                    if not header:
                        continue
                    header = header.split()[0]
                    if not header.upper().startswith(prefix):
                        continue
                    mutation = header[len(prefix):]
                    if mutation:
                        expected.append((gene_key, mutation))
        except OSError as exc:
            print(f"Warning: could not read synthesized mutants for {gene}: {exc}")
    return expected


def run_full_pipeline_mode(args, executor_fn, ape_bin):
    """Synthesize FASTAs, run NetPhos, parse, build ensemble."""
    if not args.mapping_dir:
        print("ERROR: --mapping-dir is required for full-pipeline mode")
        return 1

    wt_header = getattr(args, 'wt_header', 'ORF')
    verbose = getattr(args, 'verbose', False)

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

    expected_mutations = _expected_mutations_from_synthesis(synth_summary)

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
        threshold=args.threshold, delta_threshold=0.05, yes_only=args.yes_only,
        expected_mutations=expected_mutations)

    n_missing = sum(1 for r in summary if 'missing_mut' in r['qc_flags'])
    if n_missing:
        print(f"WARNING: {n_missing} submitted mutation(s) produced no NetPhos "
              f"predictions; written as summary rows flagged missing_mut")

    genes_in_results = set()
    for row in summary:
        genes_in_results.add(row.get('Gene', ''))
    for row in events:
        genes_in_results.add(row.get('Gene', ''))
    for row in sites:
        genes_in_results.add(row.get('Gene', ''))
    genes_in_results.discard('')

    for gene in genes_in_results:
        gene_summary = [r for r in summary if r.get('Gene') == gene]
        gene_events = [r for r in events if r.get('Gene') == gene]
        gene_sites = [r for r in sites if r.get('Gene') == gene]
        gene_out = Path(args.output) / gene / "NetPhos"
        gene_out.mkdir(parents=True, exist_ok=True)
        output_base = str(gene_out / gene)
        write_ensemble_outputs(output_base, gene_summary, gene_events, gene_sites)

    if not genes_in_results:
        gene_out = Path(args.output) / "NetPhos"
        gene_out.mkdir(parents=True, exist_ok=True)
        write_ensemble_outputs(str(gene_out / "output"), summary, events, sites)

    # Cleanup
    shutil.rmtree(work_dir, ignore_errors=True)

    if temp_holder:
        temp_holder.cleanup()

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
                        help='Input: WT FASTA file/directory')
    parser.add_argument('output', nargs='?',
                        help='Output base directory (writes {GENE}/NetPhos/{GENE}.tsv, .events.tsv, .sites.tsv)')
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
    parser.add_argument('--verbose', action='store_true',
                        help='Verbose output')

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

    # --yes-only redefines a "site" as a NetPhos YES call inside the ensemble
    # (see build_netphos_ensemble); it no longer zeroes the threshold.

    # Resolve APE binary (needed for execution modes)
    native_ape = None
    executor_fn = None
    ape_bin = None

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

    return run_full_pipeline_mode(args, executor_fn, ape_bin)


if __name__ == '__main__':
    exit(main())
