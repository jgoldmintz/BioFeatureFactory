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

"""

import argparse
import csv
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
    extract_mutation_from_sequence_name,
    extract_gene_from_filename,
    write_tsv,
    write_fasta,
    detect_alphabet,
    trim_muts,
    should_skip_mutation,
    parse_variant,
    protein_consequence,
    infer_edit_span,
    align_wt_to_mut,
    splice_seq,
    is_intronic_token,
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

# Classification codes. The original six occupy -3..2 (subthreshold -3, lost -2,
# weakened -1, stable 0, strengthened 1, gained 2); -4/+4 are the residue-level
# outcomes that only exist once alleles can differ in length, placed outside that
# range so an old consumer cannot mistake one for a comparison.
_CODE_DELETED_RESIDUE = -4
_CODE_INSERTED_RESIDUE = 4


def _classify_netphos_event(wt_score, mut_score, wt_above, mut_above, delta_threshold=0.05):
    """Classify a single ALIGNED (position, kinase) pair between WT and MUT.

    ALIGNED is a precondition: the residue exists in BOTH alleles. Callers must
    route a position with no counterpart -- one the edit deleted, or one the edit
    inserted -- to the deleted/inserted rows instead, because the arithmetic below
    would fabricate a delta against a residue that does not exist.

    `wt_above`/`mut_above` are the site-membership booleans decided by the caller:
    score >= threshold, or (under --yes-only) NetPhos answer == 'YES'.
    Returns (classification, classification_code, delta).

    The None -> 0.0 coalescence below is deliberate and stays. Here it means "the
    residue is present in this allele and NetPhos reported no site for this kinase
    on it", i.e. below the tool's own reporting cut -- which is a measurement, and
    is exactly what the gained/lost classification is built on. That is a
    different fact from "the residue is not there at all", which is what the
    caller now handles separately and reports as an EMPTY delta rather than a
    zero. Coalescing the second case would claim a measurement that was never
    made; coalescing the first is the pipeline's existing definition of "not a
    site" and changing it would silently redefine every gained/lost call.
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
# Synthesis with per-token accounting and residue-frame alignment
#
# netphos does NOT call synthesize_gene_fastas. Three reasons, in order of
# weight:
#
#   1. synthesize_gene_fastas (utility.py:2380) calls
#      build_mutant_sequences_for_gene without passing `non_snp`, so it always
#      takes the default False and the shared non-SNV protein builder
#      (_non_snv_mutant_protein) is unreachable through it. Every indel token
#      dies at synthesis.
#   2. Comparing two alleles across a length change needs the WT and MUT PROTEIN
#      STRINGS in memory to build the residue projection. synthesize_gene_fastas
#      returns only paths and a count.
#   3. Naming a rejected token needs the REQUESTED token set. Recovering it from
#      the synthesized mutant FASTA -- what this pipeline used to do -- cannot
#      see a token that failed synthesis, which is precisely the set that needs
#      naming.
#
# build_mutant_sequences_for_gene is called directly instead. That is an already
# exported shared function, already called directly this way by netsurfp3 and
# netmhc; nothing is added to or changed in utils/utility.py.
# ---------------------------------------------------------------------------

# Mirrors build_mutant_sequences_for_gene's own format detection
# (utility.py:2242-2250) so the requested-token set is exactly the set that
# function iterates. Kept in sync by construction: same keys, same precedence.
_MUTANT_COLUMN_KEYS = ['mutant', 'mutation', 'nt_mutation', 'ntmutant']


def _requested_tokens(mapping_file, log_path, gene_name, failure_map):
    """Return the mutation tokens this gene's mapping file asks for, in order.

    This is the denominator for every accounting statement the pipeline makes. It
    comes from the mapping file -- the source of truth -- not from the synthesized
    FASTA, so a token that failed synthesis is still counted and can still be
    named.

    The validation-log and failure-map filters are applied here for the same
    reason build_mutant_sequences_for_gene applies them: a mutation deliberately
    excluded upstream is not a rejection and must not be reported as one.
    """
    if not mapping_file or not os.path.exists(mapping_file):
        return []
    try:
        with open(mapping_file, 'r') as handle:
            lines = handle.readlines()
    except OSError as exc:
        print(f"Warning: could not read mapping file for {gene_name}: {exc}")
        return []
    if not lines:
        return []

    is_single_column = True
    if ',' in lines[0]:
        first_line_lower = lines[0].lower()
        if any(k in first_line_lower for k in ['mutant', 'mutation', 'aamutant']):
            is_single_column = False

    raw = []
    if is_single_column:
        for line in lines:
            tok = line.strip()
            if tok and tok.lower() != 'mutant':
                raw.append(tok)
    else:
        reader = csv.DictReader(lines)
        for row in reader:
            for key in _MUTANT_COLUMN_KEYS:
                if row.get(key):
                    raw.append(row[key].strip())
                    break

    # Same allow-list the builder derives from the validation log.
    allowed = None
    if log_path:
        try:
            allowed = {
                entry.split(',')[0].strip().upper()
                for entry in trim_muts(mapping_file, log=log_path, gene_name=gene_name)
                if entry
            }
        except Exception:
            allowed = None

    out = []
    for tok in raw:
        clean = tok.replace(" ", "")
        if not clean:
            continue
        if allowed is not None and clean.upper() not in allowed:
            continue
        if should_skip_mutation(gene_name, clean, failure_map):
            continue
        out.append(clean)
    return out


def _aa_level_non_snv_protein(token, wt_prot):
    """Mutant protein for a NON-SNV AMINO-ACID token. Returns the sequence or None.

    Only reached when the WT input is a protein, so there is no ORF and no
    nucleotide token. The shared builder cannot serve this case: its aa path goes
    through get_mutation_data_bioAccurate, which does int(token[1:-1]) and so
    raises on any multi-residue token, and _non_snv_mutant_protein deliberately
    works at the nucleotide level because an indel's protein effect normally
    requires retranslation.

    When the token is ALREADY written in residues there is nothing to translate:
    the edit is exactly splice_seq one level up, with the same whole-span REF
    guard. An aa token cannot express a frameshift -- it names a bounded residue
    span on both sides -- so the result is always a well-defined in-frame protein.

    Returns None for an SNV (the untouched shared path handles it), for a token
    that does not parse, and for one whose REF does not match the WT protein;
    the caller names each of those instead.
    """
    variant = parse_variant(token, is_nt=False)
    if variant is None or variant.is_snv or not wt_prot:
        return None
    try:
        # validate=True: the whole REF span must match, not just its first
        # residue. A multi-residue REF whose first character happens to line up
        # would otherwise splice at a wrong coordinate and produce a plausible,
        # wrong protein.
        return splice_seq(wt_prot, variant.pos0, variant.ref, variant.alt,
                          validate=True)
    except ValueError:
        return None


def _rejection_reason(token, orf_nt, is_nt_input):
    """Name why a requested token produced no mutant protein.

    Returning a reason string is the whole point: the alternative -- and the
    previous behaviour -- is that the token vanishes, which downstream is
    indistinguishable from "this gene never had that mutation".

    Reason codes follow the shared convention: UPPER_SNAKE, optionally with a
    `:detail` tail carrying the observed value.
    """
    # An out-of-ORF token is WELL FORMED and was declined on scope, not syntax.
    # UNPARSEABLE_TOKEN gives it the same reason as actual garbage, which defeats
    # the purpose of naming reasons: build_mutant_sequences_for_gene already
    # identifies it correctly on stderr, and this string is the only record that
    # reaches the TSV.
    if is_intronic_token(token):
        return 'NON_ORF_TOKEN:no_reading_frame_at_protein_level'
    variant = parse_variant(token, is_nt=is_nt_input)
    if variant is None:
        return 'UNPARSEABLE_TOKEN'
    if not is_nt_input or not orf_nt:
        # Amino-acid input: there is no ORF to check the token against. A non-SNV
        # aa token only reaches here after _aa_level_non_snv_protein declined it,
        # which it does exactly when the REF span disagrees with the WT protein.
        if not variant.is_snv:
            return 'REF_MISMATCH:wt_protein_span_differs'
        return 'SYNTHESIS_SKIPPED:aa_input'
    end = variant.pos0 + len(variant.ref)
    if end > len(orf_nt):
        return f'REF_SPANS_PAST_ORF:{end}>{len(orf_nt)}'
    observed = orf_nt[variant.pos0:end].upper()
    if observed != variant.ref.upper():
        return f'REF_MISMATCH:orf_has_{observed}'
    cons = protein_consequence(variant, orf_nt)
    if cons is None:
        return 'NO_PROTEIN_CONSEQUENCE'
    # An SNV that creates a premature stop is refused by infer_aamutation_from_nt
    # (utility.py:2061) because codon_to_aa renders a stop as the 4-character
    # string 'Stop', which must not be spliced in as a residue. That refusal is
    # correct and is the SNV path's long-standing behaviour; what was wrong is
    # that the token then disappeared without a trace.
    if variant.is_snv and cons['aa_consequence'] == 'stop_gained':
        return 'STOP_GAINED:snv_path_does_not_synthesize'
    return f"SYNTHESIS_SKIPPED:{cons['aa_consequence']}"


def _alignment_context(token, orf_nt, wt_prot, mut_prot, is_nt_input):
    """Residue-frame projection inputs for one synthesized mutant.

    Returns the dict build_netphos_ensemble needs to map a WT residue position
    onto the mutant residue it actually became.

    The span is recovered from the two PROTEIN strings that were actually handed
    to NetPhos, not from the nucleotide token's idealised codon span, because the
    mutant protein is truncated at its first stop (utility.py:2186-2187) and the
    codon span does not describe that. infer_edit_span's docstring names this
    exact use: pipelines that hold a WT and a MUT sequence but not the variant
    record.

    TWO cases must be taken away from infer_edit_span, because its prefix/suffix
    trimming reports the MINIMAL edit and both of these are larger than minimal:

      frameshift -- every residue from the edit onward is a different one.
          infer_edit_span cannot see this and documents that it must be told;
          protein_consequence supplies the fact.

      re-termination -- the edit moved the stop codon, so the mutant's
          C-terminus is not the WT's C-terminus and the two tails are unrelated
          sequence. Trimming a common suffix off them pairs residues that are not
          counterparts: a truncating delins whose mutant happens to end in the
          same residues as the full protein gets its last WT residues mapped onto
          the mutant's last residues, fabricating gained/lost events while the
          residues' true counterparts are separately reported deleted. The
          accounting in _variant_qc_flags is identical either way, so nothing
          downstream can detect it.

    Re-termination is detected from lengths, not from the consequence string: for
    an in-frame edit the protein length changes by exactly the codon delta unless
    translation now stops somewhere else. When it does, only the common PREFIX is
    trusted and the edit is declared to run to the end of both alleles -- which is
    what an early stop (nothing after it exists) or a read-through (everything
    after the old stop is new) actually is.
    """
    consequence = ''
    new_stop_aa_pos = None
    frameshift = False

    # Whether the TOKEN is an SNV, decided at the level the token is written at.
    # This is not the same question as the protein consequence and must not be
    # inferred from it in either direction: an MNV such as TAT13TTT changes one
    # residue, so its protein consequence is 'snv' while the token is not; and a
    # synonymous SNV changes none, so its consequence is 'synonymous' while the
    # token is an ordinary SNV whose output must not move.
    variant = parse_variant(token, is_nt=bool(is_nt_input))
    is_snv_token = variant.is_snv if variant is not None else False

    if is_nt_input and orf_nt and variant is not None:
        cons = protein_consequence(variant, orf_nt)
        if cons is not None:
            consequence = cons['aa_consequence']
            new_stop_aa_pos = cons['new_stop_aa_pos']
            frameshift = consequence == 'frameshift'

    if frameshift:
        offset, ref_len, alt_len = infer_edit_span(wt_prot, mut_prot, frameshift=True)
    else:
        # In frame, the protein length changes by exactly the codon delta -- unless
        # translation now terminates somewhere else. An aa-level token retranslates
        # nothing, so its expected delta is whatever was observed and this can never
        # fire for one; an SNV's is 0 with equal lengths, so the SNV path is
        # untouched.
        observed_delta = len(mut_prot) - len(wt_prot)
        if is_nt_input and variant is not None:
            expected_delta = variant.length_delta // 3
        else:
            expected_delta = observed_delta
        if observed_delta == expected_delta:
            offset, ref_len, alt_len = infer_edit_span(wt_prot, mut_prot)
        else:
            offset = 0
            while (offset < min(len(wt_prot), len(mut_prot))
                   and wt_prot[offset] == mut_prot[offset]):
                offset += 1
            ref_len = len(wt_prot) - offset
            alt_len = len(mut_prot) - offset

    if not consequence:
        # aa-token input, or a token whose consequence could not be derived.
        # Name the class from the observed protein lengths rather than leaving
        # the column empty -- the lengths are a fact of the two sequences.
        if ref_len == alt_len == 0:
            consequence = 'synonymous'
        elif ref_len == alt_len:
            consequence = 'snv' if ref_len == 1 else 'mnv'
        elif alt_len == 0:
            consequence = 'inframe_del'
        elif ref_len == 0:
            consequence = 'inframe_ins'
        else:
            consequence = 'inframe_delins'

    return {
        'aa_offset': offset,
        'aa_ref_len': ref_len,
        'aa_alt_len': alt_len,
        'aa_consequence': consequence,
        'is_snv_token': is_snv_token,
        'new_stop_aa_pos': new_stop_aa_pos,
        'n_aa_wt': len(wt_prot),
        'n_aa_mut': len(mut_prot),
    }


def synthesize_with_context(wt_sequences, mapping_lookup, sequence_root,
                            log_path=None, failure_map=None):
    """Write WT/MUT protein FASTAs and return the accounting needed to compare them.

    Returns (wt_dir, mut_dir, summary, variant_ctx, rejected):
      summary      per-gene dict, same shape synthesize_gene_fastas returned
      variant_ctx  (gene, token) -> _alignment_context dict, for every mutant built
      rejected     (gene, token) -> reason string, for every requested token that
                   was not built

    The gene half of both keys is normalized exactly as
    _pair_predictions_with_mutations normalizes gene names parsed back out of the
    NetPhos output, so the three dictionaries join.
    """
    sequence_root = Path(sequence_root)
    wt_dir = sequence_root / "wt"
    mut_dir = sequence_root / "mut"
    wt_dir.mkdir(parents=True, exist_ok=True)
    mut_dir.mkdir(parents=True, exist_ok=True)

    summary = []
    variant_ctx = {}
    rejected = {}

    for gene_name, wt_seq in wt_sequences.items():
        gene_name = gene_name.upper()
        gene_key = (extract_gene_from_filename(gene_name) or gene_name).upper()
        seq_upper = wt_seq.strip().upper()
        mapping_file = mapping_lookup.get(gene_name)
        # Read the requested set BEFORE the alphabet dispatch. A gene the dispatch
        # cannot use still had mutations asked of it, and skipping the gene used to
        # take all of them with it -- no ctx entry, no rejection, so no row at all
        # and a stdout line as the only trace. Naming the reason per token is the
        # same rule the per-token path already follows; the gene-level branches
        # were the one hole left in it.
        requested = _requested_tokens(mapping_file, log_path, gene_name, failure_map)

        def _decline_gene(reason):
            for tok in requested:
                rejected[(gene_key, tok)] = reason

        # Same alphabet dispatch as synthesize_gene_fastas (utility.py:2393-2415):
        # nucleotide -> translate; protein -> use as-is; codon-encoded -> skip.
        try:
            detected = detect_alphabet(seq_upper)
        except ValueError:
            print(f"Skipping {gene_name}: empty sequence")
            _decline_gene('GENE_SKIPPED:empty_wt_sequence')
            continue
        if detected == 'nucleotide':
            nt_for_build = seq_upper
            aa_seq = translate_orf_sequence(seq_upper)
            build_input_type = 'nt'
            if not aa_seq:
                print(f"Skipping {gene_name}: unable to translate ORF")
                _decline_gene('GENE_SKIPPED:orf_does_not_translate')
                continue
        elif detected == 'protein':
            nt_for_build = None
            aa_seq = seq_upper
            build_input_type = 'aa'
        else:
            print(f"Skipping {gene_name}: codon-encoded input, expected nt or aa")
            _decline_gene(f'GENE_SKIPPED:alphabet_is_{detected}')
            continue

        wt_path = wt_dir / f"{gene_name}-wt.fasta"
        write_fasta(wt_path, {f"{gene_name}-wt": aa_seq})

        # non_snp=True is not a user preference and gets no flag. It routes each
        # token through _non_snv_mutant_protein, which returns None for an SNV and
        # lets it fall through to the untouched SNV path; whether a token is
        # length-changing is a fact of the record, decided per token.
        mutant_sequences = build_mutant_sequences_for_gene(
            gene_name, nt_for_build, aa_seq, mapping_file, log_path, failure_map,
            input_type=build_input_type, non_snp=True,
        )

        # Protein WT input only: recover the non-SNV aa tokens the shared builder
        # structurally cannot express (see _aa_level_non_snv_protein). Done here
        # rather than in utils/utility.py because the shared builder is not
        # modified by this change; an SNV never enters this branch, so the
        # existing aa path is untouched.
        if build_input_type == 'aa':
            for token in requested:
                header = f"{gene_name}-{token}"
                if header in mutant_sequences:
                    continue
                built = _aa_level_non_snv_protein(token, aa_seq)
                if built:
                    mutant_sequences[header] = built

        mut_path = None
        if mutant_sequences:
            mut_path = mut_dir / f"{gene_name}_aa.fasta"
            write_fasta(mut_path, mutant_sequences)

        prefix = f"{gene_name}-"
        for header, mut_prot in mutant_sequences.items():
            token = header[len(prefix):] if header.startswith(prefix) else header
            variant_ctx[(gene_key, token)] = _alignment_context(
                token, nt_for_build, aa_seq, mut_prot, build_input_type == 'nt')

        for token in requested:
            if f"{gene_name}-{token}" in mutant_sequences:
                continue
            rejected[(gene_key, token)] = _rejection_reason(
                token, nt_for_build, build_input_type == 'nt')

        summary.append({
            "gene": gene_name,
            "wt_path": str(wt_path),
            "mut_path": str(mut_path) if mut_path else None,
            "mutant_count": len(mutant_sequences),
        })

    return wt_dir, mut_dir, summary, variant_ctx, rejected


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


def _best_by_position_kinase(preds):
    """Collapse predictions to one per (position, kinase), keeping the best score."""
    best = {}
    for p in preds:
        key = (p['pos'], p['kinase'])
        if key not in best or p['score'] > best[key]['score']:
            best[key] = p
    return best


def _residue_projection(ctx):
    """Return (wt_to_mut, mut_to_wt) 0-based residue index maps, or (None, None).

    (None, None) means "no context available -- treat WT position i and MUT
    position i as the same residue", which is the identity assumption and the
    only correct one for an equal-length pair.

    For any variant whose two alleles are the same length -- every SNV, every MNV
    -- align_wt_to_mut returns the identity, so routing the SNV path through this
    projection changes nothing about its output. That is what makes the join
    correct under an indel without disturbing the existing behaviour.
    """
    if not ctx or not ctx.get('n_aa_wt'):
        return None, None
    wt_to_mut = align_wt_to_mut(ctx['n_aa_wt'], ctx['aa_offset'],
                                ctx['aa_ref_len'], ctx['aa_alt_len'])
    n_mut = ctx['n_aa_mut']
    mut_to_wt = {}
    for i, j in enumerate(wt_to_mut):
        # The range test is a bounds INVARIANT, not a case handler. Truncation is
        # already carried by aa_alt_len -- _alignment_context sets it from the
        # mutant's real length -- so a span produced there cannot project past the
        # mutant (verified over 400k random allele pairs: 0 out-of-range). It
        # stays because build_netphos_ensemble accepts a caller-supplied ctx, and
        # a hand-built span with a wrong alt_len must degrade to "no counterpart"
        # rather than index a residue that does not exist.
        if j is None or not (0 <= j < n_mut):
            wt_to_mut[i] = None
            continue
        mut_to_wt[j] = i
    return wt_to_mut, mut_to_wt


def build_netphos_ensemble(wt_preds_by_gene, mut_preds_by_mutation, variant_ctx,
                           threshold=0.5, delta_threshold=0.05, yes_only=False,
                           expected_mutations=None, rejected=None):
    """Build ensemble comparison tables from parsed WT and MUT predictions.

    Args:
        wt_preds_by_gene: dict  gene_name -> list[pred_dict]
        mut_preds_by_mutation: dict  (gene, nt_mutation) -> list[pred_dict]
        variant_ctx: dict  (gene, nt_mutation) -> alignment context from
            _alignment_context. Supplies the residue projection, so a WT position
            is compared with the residue it actually became rather than with
            whatever slid into its index. A missing entry falls back to the
            identity projection.
        threshold: score threshold for YES/NO
        delta_threshold: minimum absolute delta for strengthened/weakened
        expected_mutations: optional iterable of (gene, nt_mutation) that were
            submitted to NetPhos. Any submitted mutation with no parsed
            predictions (dropped batch, failed run) gets a zeroed summary row
            flagged missing_mut instead of vanishing from the output.
        rejected: optional dict (gene, nt_mutation) -> reason, for tokens that
            never reached NetPhos because synthesis declined them. Each gets a
            summary row carrying the reason, so a rejected token is named rather
            than absent.

    The `mapping_lookup` parameter this function used to take was never read in
    the body; it is replaced by variant_ctx rather than kept alongside it.

    Returns:
        (summary_rows, events_rows, sites_rows)
    """
    summary_rows = []
    events_rows = []
    sites_rows = []
    rejected = rejected or {}

    def _eventless_summary_row(pkey, gene, n_sites_wt, qc_flags, mut_measured):
        """Summary row for a variant that produced no events at all.

        Both no-event cases go through here -- a mutation submitted to NetPhos
        that came back with no predictions, and a token synthesis declined so it
        was never submitted -- so the two cannot drift apart in shape.

        mut_measured says whether a mutant protein was ever scored, and it is the
        only thing that decides between 0 and EMPTY on the mutant-side columns.
        A rejected token has no mutant protein at all, so `n_sites_mut: 0` and
        `count_lost: 0` would be fabricated observations -- "the mutant had no
        phosphosites and lost none" is a measurement, and none was made. They are
        emitted empty, with the reason already in qc_flags.

        missing_mut keeps its zeros: there the mutant WAS submitted and the run
        came back empty, and that row shape is pre-existing output. Changing it
        would move bytes on an SNV run whose batch was dropped.
        """
        n = 0 if mut_measured else ''
        f = 0.0 if mut_measured else ''
        return {
            'pkey': pkey, 'Gene': gene,
            'n_sites_wt': n_sites_wt, 'n_sites_mut': n,
            'count_gained': n, 'count_lost': n, 'count_strengthened': n,
            'count_weakened': n, 'count_stable': n,
            'max_abs_delta': f, 'sum_abs_delta': f, 'n_kinases_affected': n,
            'top_event_type': '', 'top_event_delta': f,
            'top_event_position': '', 'top_event_kinase': '',
            'top_event_classification_code': '',
            'qc_flags': '|'.join(qc_flags),
        }

    def _wt_site_count(wt_map):
        if yes_only:
            return sum(1 for p in wt_map.values() if p['answer'] == 'YES')
        return sum(1 for p in wt_map.values() if p['score'] >= threshold)

    for (gene, nt_mutation), mut_preds in mut_preds_by_mutation.items():
        pkey = f"{gene}-{nt_mutation}"
        wt_preds = wt_preds_by_gene.get(gene, [])
        ctx = variant_ctx.get((gene, nt_mutation)) if variant_ctx else None
        wt_to_mut, mut_to_wt = _residue_projection(ctx)

        wt_map = _best_by_position_kinase(wt_preds)
        mut_map = _best_by_position_kinase(mut_preds)

        def _project(wt_pos):
            """WT residue position (1-based) -> (MUT position or None, align_status).

            The range check is not redundant and is not a deletion: wt_pos comes
            from parsed NetPhos output, not from the protein string, so a gene key
            that aggregated predictions from two WT files, or a malformed output
            line, can carry a position the WT protein does not have. That has no
            counterpart, but calling it 'deleted' would blame the variant for a
            data problem -- so it gets its own name. An IndexError here would take
            the whole run down instead.
            """
            if wt_to_mut is None:
                return wt_pos, 'aligned'
            idx = wt_pos - 1
            if not (0 <= idx < len(wt_to_mut)):
                return None, 'wt_position_unmapped'
            j = wt_to_mut[idx]
            return (None, 'deleted') if j is None else (j + 1, 'aligned')

        def _to_wt_pos(mut_pos):
            """MUT residue position (1-based) -> WT residue position, or None."""
            if mut_to_wt is None:
                return mut_pos
            j = mut_to_wt.get(mut_pos - 1)
            return None if j is None else j + 1

        # --- sites rows (raw predictions tagged with allele) ---
        # WT rows are emitted once per pkey, i.e. repeated across the gene's
        # mutations. That is denormalization, not duplication by accident: the
        # pkey IS the WT/MUT pair, so `sites[sites.pkey == k]` has to return both
        # alleles of that comparison. position_wt_frame is what makes the two
        # halves joinable once an indel has made `position` allele-specific.
        for p in wt_preds:
            sites_rows.append({
                'pkey': pkey, 'Gene': gene, 'allele': 'WT',
                'seq_name': p['seq_name'], 'position': p['pos'],
                'position_wt_frame': p['pos'],
                'amino_acid': p['amino_acid'], 'context': p['context'],
                'score': p['score'], 'kinase': p['kinase'], 'answer': p['answer'],
            })
        for p in mut_preds:
            origin = _to_wt_pos(p['pos'])
            sites_rows.append({
                'pkey': pkey, 'Gene': gene, 'allele': 'MUT',
                'seq_name': p['seq_name'], 'position': p['pos'],
                'position_wt_frame': origin if origin is not None else '',
                'amino_acid': p['amino_acid'], 'context': p['context'],
                'score': p['score'], 'kinase': p['kinase'], 'answer': p['answer'],
            })

        # --- events rows ---
        # Keyed on the WT residue where one exists. Three cases, and every one
        # produces a row:
        #   aligned  -- the residue exists in both alleles
        #   deleted  -- a WT residue the edit removed (or truncated away)
        #   inserted -- a mutant residue with no WT origin
        # (The aligned_N/M figure in qc_flags is counted over RESIDUES, not over
        # these rows -- see _variant_qc_flags for why.)
        wt_keyed = {}          # (wt_pos, kinase) -> (mut_pos or None, align_status)
        for wt_pos, kinase in wt_map:
            wt_keyed[(wt_pos, kinase)] = _project(wt_pos)
        inserted_keys = []     # (mut_pos, kinase) with no WT origin
        for mut_pos, kinase in mut_map:
            origin = _to_wt_pos(mut_pos)
            if origin is None:
                inserted_keys.append((mut_pos, kinase))
            else:
                wt_keyed.setdefault((origin, kinase), (mut_pos, 'aligned'))

        mutation_events = []
        for (wt_pos, kinase) in sorted(wt_keyed):
            mut_pos, status = wt_keyed[(wt_pos, kinase)]
            wp = wt_map.get((wt_pos, kinase))
            mp = mut_map.get((mut_pos, kinase)) if mut_pos is not None else None

            if mut_pos is None:
                # No counterpart residue. There is nothing to subtract, so the
                # delta is EMPTY -- a 0.0 here would read downstream as
                # "measured, no change", a fabricated observation.
                #
                # wp is indexed directly, not guarded: a key can only carry a
                # None mut_pos if it came from wt_map, because the other source
                # (a mutant position's WT origin) exists precisely when the
                # projection is not None. A guard here would be unreachable, and
                # if the invariant ever broke a KeyError is the right outcome.
                event = {
                    'pkey': pkey, 'Gene': gene,
                    'position': wt_pos, 'position_mut': '',
                    'align_status': status,
                    'amino_acid_wt': wp['amino_acid'],
                    'amino_acid_mut': '',
                    'kinase': kinase,
                    'wt_score': wp['score'],
                    'mut_score': '',
                    'delta': '',
                    'wt_answer': wp['answer'],
                    'mut_answer': '',
                    'classification': ('deleted_residue' if status == 'deleted'
                                       else 'unmapped_residue'),
                    'classification_code': _CODE_DELETED_RESIDUE,
                }
                events_rows.append(event)
                mutation_events.append(event)
                continue

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
                'position': wt_pos,
                'position_mut': mut_pos,
                'align_status': status,
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

        for (mut_pos, kinase) in sorted(inserted_keys):
            mp = mut_map[(mut_pos, kinase)]
            event = {
                'pkey': pkey, 'Gene': gene,
                'position': '', 'position_mut': mut_pos,
                'align_status': 'inserted',
                'amino_acid_wt': '', 'amino_acid_mut': mp['amino_acid'],
                'kinase': kinase,
                'wt_score': '', 'mut_score': mp['score'],
                'delta': '',
                'wt_answer': '', 'mut_answer': mp['answer'],
                'classification': 'inserted_residue',
                'classification_code': _CODE_INSERTED_RESIDUE,
            }
            events_rows.append(event)
            mutation_events.append(event)

        # --- summary row ---
        n_sites_wt = _wt_site_count(wt_map)
        if yes_only:
            n_sites_mut = sum(1 for p in mut_map.values() if p['answer'] == 'YES')
        else:
            n_sites_mut = sum(1 for p in mut_map.values() if p['score'] >= threshold)

        count_gained = sum(1 for e in mutation_events if e['classification'] == 'gained')
        count_lost = sum(1 for e in mutation_events if e['classification'] == 'lost')
        count_strengthened = sum(1 for e in mutation_events if e['classification'] == 'strengthened')
        count_weakened = sum(1 for e in mutation_events if e['classification'] == 'weakened')
        count_stable = sum(1 for e in mutation_events if e['classification'] == 'stable')

        # Only aligned events carry a numeric delta; deleted/inserted ones carry
        # '' and must not be swept into the magnitude statistics as zeros.
        scored = [e for e in mutation_events if e['delta'] != '']
        abs_deltas = [abs(e['delta']) for e in scored]
        max_abs_delta = max(abs_deltas) if abs_deltas else 0.0
        # float() so the column stays float-typed when nothing was comparable:
        # sum([]) is int 0 and would render as "0" in a column that is "0.0"
        # everywhere else. A no-op on a non-empty list of floats.
        sum_abs_delta = float(sum(abs_deltas))

        non_stable = [e for e in mutation_events if e['classification'] not in ('stable', 'subthreshold')]
        n_kinases_affected = len({e['kinase'] for e in non_stable})

        top_event = max(scored, key=lambda e: abs(e['delta'])) if scored else None

        qc_flags = []
        if not wt_preds:
            qc_flags.append("missing_wt")
        if not mut_preds:
            qc_flags.append("missing_mut")
        if scored and max_abs_delta == 0.0:
            qc_flags.append("no_delta")
        elif mutation_events and not scored:
            # Events exist but none of them is a comparison: every site either
            # lost its residue or gained one. "no_delta" would claim the deltas
            # were measured and came out zero, which is the fabrication this
            # whole pass removes. (For an SNV `scored` is always the full event
            # list, so this branch cannot fire and the flag above is unchanged.)
            qc_flags.append("no_comparable_sites")
        qc_flags.extend(_variant_qc_flags(ctx, wt_to_mut))

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
            # EMPTY, not 0.0, when NOTHING was comparable. `scored` is empty only
            # when every event is a deleted or inserted residue -- a frameshift
            # that removes all of a protein's phosphosites reported max/sum 0.0,
            # indistinguishable downstream from a variant that genuinely changed
            # nothing. no_comparable_sites already flags the row, but no mean or
            # ranking reads a qc column. For an SNV `scored` is always the full
            # event list, so this cannot fire on one: SNV output is byte-identical.
            'max_abs_delta': round(max_abs_delta, 6) if scored else '',
            'sum_abs_delta': round(sum_abs_delta, 6) if scored else '',
            'n_kinases_affected': n_kinases_affected,
            'top_event_type': top_event['classification'] if top_event else '',
            # '' not 0.0: with no comparable event there is no top event, and 0.0
            # names a measured magnitude. _eventless_summary_row already uses ''.
            'top_event_delta': round(top_event['delta'], 6) if top_event else '',
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
        wt_map = _best_by_position_kinase(wt_preds)

        qc_flags = []
        if not wt_preds:
            qc_flags.append("missing_wt")
        # A token synthesis declined never reached NetPhos, so "missing_mut" on
        # its own would blame the run for a decision made before the run. Carry
        # the reason instead.
        reason = rejected.get((gene, nt_mutation))
        qc_flags.append("rejected_token" if reason else "missing_mut")
        if reason:
            qc_flags.append(reason)

        summary_rows.append(_eventless_summary_row(
            f"{gene}-{nt_mutation}", gene, _wt_site_count(wt_map), qc_flags,
            mut_measured=reason is None))

    return summary_rows, events_rows, sites_rows


def _variant_qc_flags(ctx, wt_to_mut):
    """Named qc tokens describing the variant class and the alignment accounting.

    Follows the RNAfold/codon_usage convention: the consequence class and the
    union accounting live in the qc column rather than in new columns, so an SNV
    row's qc_flags stays empty and every other column keeps its meaning across
    variant classes.

    The counts are over RESIDUES, not over scored events. That distinction is the
    whole point of the accounting: an insertion of a residue no kinase scores
    produces no unmatched event, so an event-level count would report it as
    aligned_N/N -- "fully aligned" -- while the protein demonstrably gained a
    residue. Counting the residue union states the sequence fact regardless of
    whether NetPhos happened to score the inserted residue.
    """
    if not ctx:
        return []
    flags = []
    consequence = ctx.get('aa_consequence')
    # Gated on the TOKEN being non-SNV, not on the protein consequence. An SNV's
    # qc_flags must stay exactly as it was, including for a synonymous one whose
    # consequence string is not 'snv'; and every non-SNV token must be named,
    # including an MNV whose one-residue effect makes its consequence read 'snv'.
    if consequence and not ctx.get('is_snv_token'):
        flags.append(f"aa_consequence:{consequence}")
    if consequence == 'frameshift' and ctx.get('new_stop_aa_pos'):
        flags.append(f"new_stop_aa:{ctx['new_stop_aa_pos']}")

    n_wt, n_mut = ctx.get('n_aa_wt', 0), ctx.get('n_aa_mut', 0)
    if n_wt != n_mut:
        n_aligned = sum(1 for j in wt_to_mut if j is not None) if wt_to_mut else n_wt
        n_deleted = n_wt - n_aligned
        # Mutant residues with no WT origin. Denominator is the UNION of both
        # alleles: n_wt alone would report an insertion as fully aligned however
        # large it is, because every WT residue does keep a counterpart -- the
        # residues without one are all on the mutant side.
        n_inserted = n_mut - n_aligned
        flags.append(
            f"length_changed:{n_mut - n_wt:+d}aa;"
            f"aligned_{n_aligned}/{n_wt + n_inserted};"
            f"deleted_{n_deleted};inserted_{n_inserted}")
    return flags


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

    # Events. position/position_mut are the WT and MUT residue coordinates of the
    # same site; they differ once an indel shifts the frame, and align_status says
    # which of the two exists. A deleted or inserted residue carries an EMPTY
    # delta, never 0.0.
    events_fields = [
        'pkey', 'Gene', 'position', 'position_mut', 'align_status',
        'amino_acid_wt', 'amino_acid_mut',
        'kinase', 'wt_score', 'mut_score', 'delta',
        'wt_answer', 'mut_answer', 'classification', 'classification_code',
    ]
    write_tsv(events_rows, events_path, events_fields, extrasaction='ignore')
    print(f"Wrote {len(events_rows)} events to {events_path}")

    # Sites. `position` is in the row's own allele frame; position_wt_frame is the
    # WT coordinate it corresponds to (empty for a residue the mutant inserted),
    # which is what makes the WT and MUT halves of a pkey joinable under an indel.
    sites_fields = [
        'pkey', 'Gene', 'allele', 'seq_name', 'position', 'position_wt_frame',
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


def _pair_predictions_with_mutations(wt_dir, mut_dir):
    """Parse WT and MUT output directories and pair predictions by gene/mutation.

    Gene keys are normalized here and in synthesize_with_context by the same
    expression, so variant_ctx / rejected join onto these keys.

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

    # Synthesize FASTAs. Non-SNV tokens are handled here BY DEFAULT: there is no
    # flag, because whether a token changes length is a fact of the record, not a
    # user preference, and the token grammar is uniquely decodable so no parser
    # mode has to be selected.
    work_dir = tempfile.mkdtemp(prefix="netphos_pipeline_")
    seq_root = os.path.join(work_dir, "sequences")
    wt_dir, mut_dir, synth_summary, variant_ctx, rejected = synthesize_with_context(
        wt_sequences, mapping_lookup, seq_root,
        log_path=args.log, failure_map=failure_map,
    )

    total_mutants = sum(s['mutant_count'] for s in synth_summary)
    print(f"Synthesized FASTAs: {len(synth_summary)} genes, {total_mutants} mutants")
    if rejected:
        print(f"{len(rejected)} requested mutation(s) declined at synthesis; each "
              f"gets a summary row naming the reason")

    # The submitted set is the set the MAPPING FILES asked for -- built mutants
    # plus declined ones. Recovering it from the synthesized mutant FASTA (the
    # previous approach) could not see a token that failed synthesis, so exactly
    # the tokens that needed naming were the ones that vanished.
    expected_mutations = list(variant_ctx.keys()) + list(rejected.keys())

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
        str(wt_output_dir), str(mut_output_dir))

    # Build ensemble
    summary, events, sites = build_netphos_ensemble(
        wt_preds_by_gene, mut_preds_by_mutation, variant_ctx,
        threshold=args.threshold, delta_threshold=0.05, yes_only=args.yes_only,
        expected_mutations=expected_mutations, rejected=rejected)

    n_missing = sum(1 for r in summary if 'missing_mut' in r['qc_flags'])
    if n_missing:
        print(f"WARNING: {n_missing} submitted mutation(s) produced no NetPhos "
              f"predictions; written as summary rows flagged missing_mut")
    n_rejected = sum(1 for r in summary if 'rejected_token' in r['qc_flags'])
    if n_rejected:
        print(f"WARNING: {n_rejected} requested mutation(s) never reached NetPhos; "
              f"written as summary rows flagged rejected_token with a named reason")

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
