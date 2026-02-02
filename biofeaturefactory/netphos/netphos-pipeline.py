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
NetPhos Output Parser

Parses NetPhos output files and converts them to TSV format.
Input: NetPhos text output with phosphorylation predictions
Output: TSV with columns: Gene, pos, amino_acid, context, score, kinase, answer
With --mapping-dir: Output includes pkey column and filters by mutation positions
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

# Import shared batch processing utilities
from biofeaturefactory.utils.utility import (
    split_fasta_into_batches,
    combine_batch_outputs,
    get_mutation_data_bioAccurate_unified,
    discover_mapping_files,
    discover_fasta_files,
    parse_predictions_with_mutation_filtering,
    load_validation_failures,
    should_skip_mutation,
)

# Use the unified function from utility.py
def get_mutation_data_bioAccurate(aaposaa):
    """Extract amino acid position from mutation notation (e.g., 'D333E' -> 333)"""
    return get_mutation_data_bioAccurate_unified(aaposaa)


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
    # Check netnglyc home as well since APE ships with netnglyc
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

def load_nt_to_aa_mapping(mapping_file):
    """Load NT->AA mapping from CSV file"""
    mapping_dict = {}
    try:
        with open(mapping_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mutation_id = row['mutant']  # e.g., "A1000G"
                aaposaa = row['aamutant']    # e.g., "D333E"
                mapping_dict[mutation_id] = aaposaa
    except Exception as e:
        print(f"Error reading mapping file {mapping_file}: {e}")
    return mapping_dict

def parse_netphos_file(input_file):
    """Parse NetPhos output file and extract phosphorylation predictions"""
    predictions = []
    current_gene_base = None
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Extract gene name from header line
            if line.startswith('>') and 'amino acids' in line:
                # e.g., ">CHRNE_aa     497 amino acids" -> "CHRNE"
                match = re.match(r'>(\S+)', line)
                if match:
                    header_gene = match.group(1)
                    # Extract base gene name (remove _aa suffix if present)
                    current_gene_base = header_gene.replace('_aa', '')
                continue
            
            # Skip empty lines
            if not line:
                continue
                
            # Check if this is a prediction line (starts with # followed by gene name)
            if line.startswith('#') and current_gene_base:
                # Remove the "# " prefix for parsing
                clean_line = line[2:].strip()
                parts = clean_line.split()
                
                # Check if this looks like a prediction line
                if (len(parts) >= 7 and 
                    parts[0].startswith(current_gene_base) and 
                    parts[2] in ['S', 'T', 'Y']):
                    
                    try:
                        gene_name = parts[0]
                        position = int(parts[1])
                        amino_acid = parts[2]
                        context = parts[3]
                        score = float(parts[4])
                        kinase = parts[5]
                        answer = parts[6] if len(parts) > 6 else '.'
                        
                        predictions.append({
                            'Gene': gene_name,
                            'pos': position,
                            'amino_acid': amino_acid,
                            'context': context,
                            'score': score,
                            'kinase': kinase,
                            'answer': answer
                        })
                        
                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue
    
    return predictions

def write_tsv(predictions, output_file, include_pkey=False):
    """Write predictions to TSV file"""
    with open(output_file, 'w') as f:
        # Write header
        if include_pkey:
            f.write("pkey\tGene\tpos\tamino_acid\tcontext\tscore\tkinase\tanswer\n")
        else:
            f.write("Gene\tpos\tamino_acid\tcontext\tscore\tkinase\tanswer\n")
        
        # Write predictions
        for pred in predictions:
            if include_pkey:
                f.write(f"{pred['pkey']}\t{pred['Gene']}\t{pred['pos']}\t{pred['amino_acid']}\t{pred['context']}\t{pred['score']}\t{pred['kinase']}\t{pred['answer']}\n")
            else:
                f.write(f"{pred['Gene']}\t{pred['pos']}\t{pred['amino_acid']}\t{pred['context']}\t{pred['score']}\t{pred['kinase']}\t{pred['answer']}\n")

def parse_with_mutation_filtering_directory(input_path, mapping_dir, threshold=0.0, yes_only=False, is_mutant=False, failure_map=None):
    """Parse NetPhos files from directory or single file and filter by mutation positions"""
    all_results = []
    
    if os.path.isfile(input_path):
        # Single file processing
        results = parse_with_mutation_filtering_single_file(
            input_path, mapping_dir, threshold, yes_only, is_mutant, failure_map=failure_map
        )
        all_results.extend(results)
    elif os.path.isdir(input_path):
        # Directory processing - find all NetPhos output files
        netphos_files = []
        for ext in ['*.out', '*.txt']:
            netphos_files.extend(Path(input_path).glob(ext))
        
        if not netphos_files:
            print(f"Warning: No NetPhos output files (.out, .txt) found in directory: {input_path}")
            return []
        
        print(f"Found {len(netphos_files)} NetPhos output files in {input_path}")
        
        # Process each file
        for file_path in netphos_files:
            print(f"Processing {file_path}...")
            results = parse_with_mutation_filtering_single_file(
                str(file_path), mapping_dir, threshold, yes_only, is_mutant, failure_map=failure_map
            )
            all_results.extend(results)
    
    return all_results

def parse_with_mutation_filtering_single_file(input_file, mapping_dir, threshold=0.0, yes_only=False, is_mutant=False, failure_map=None):
    """Parse NetPhos file and filter by mutation positions, generating pkeys"""

    failure_map = failure_map or {}
    
    # First, parse all predictions from the file
    all_predictions = parse_netphos_file(input_file)
    
    if not all_predictions:
        return []
    
    # Extract gene name from predictions (assume all predictions are from same gene)
    raw_gene = all_predictions[0]['Gene']
    # Extract base gene name (e.g., "CHRNE-C915T" -> "CHRNE")
    gene = raw_gene.split('-')[0]
    
    # Load NT->AA mapping for this gene using flexible discovery
    discovered_mappings = discover_mapping_files(mapping_dir)
    if gene not in discovered_mappings:
        print(f"Warning: No mapping file found for {gene} in {mapping_dir}")
        print(f"Available mappings: {list(discovered_mappings.keys())}")
        return []
    
    mapping_file = discovered_mappings[gene]
    
    mapping_dict = load_nt_to_aa_mapping(mapping_file)
    print(f"Loaded {len(mapping_dict)} mutations for {gene}")
    
    # Use unified filtering logic
    if is_mutant:
        # Use new single-mutation processing for mutant sequences
        try:
            results = parse_predictions_with_mutation_filtering(
                all_predictions,
                mapping_dict,
                is_mutant=True,
                threshold=threshold,
                yes_only=yes_only,
                tool_type='netphos',
                failure_map=failure_map,
            )
            print(f"Predictions at mutation sites: {len(results)}")
            return results
        except Exception as e:
            print(f"Error in single-mutation processing: {e}")
            # Fall back to original logic if needed
            pass
    
    # Wildtype processing or fallback - use original bulk logic
    results = []
    matches_found = 0
    stop_codons_skipped = 0
    
    # Pre-process mapping data to create lookup structures
    position_mutations = {}  # pos -> list of (mutation_id, aaposaa, target_aa)
    
    for mutation_id, aaposaa in mapping_dict.items():
        # Skip stop codons and get position + amino acid info
        position_data = get_mutation_data_bioAccurate(aaposaa)
        if position_data[0] is None:
            stop_codons_skipped += 1
            continue
        
        aa_pos = position_data[0]  # e.g., 541
        aa_tuple = position_data[1]  # e.g., ('K', 'E') for K541E
        
        # Determine which amino acid to match based on mutant flag
        if is_mutant:
            target_aa = aa_tuple[1]  # Mutant amino acid (E for mutant)
        else:
            target_aa = aa_tuple[0]  # Original amino acid (K for wildtype)
        
        # Skip logged failures during mutant parsing
        if is_mutant and should_skip_mutation(gene, mutation_id, failure_map):
            continue

        # Group mutations by position for efficient lookup
        if aa_pos not in position_mutations:
            position_mutations[aa_pos] = []
        position_mutations[aa_pos].append((mutation_id, aaposaa, target_aa))
    
    # Now iterate through predictions and find all matching mutations
    for pred in all_predictions:
        pred_pos = pred['pos']
        pred_aa = pred['amino_acid']
        
        # Find all mutations that match this prediction's position
        if pred_pos in position_mutations:
            for mutation_id, aaposaa, target_aa in position_mutations[pred_pos]:
                # Check if amino acid matches
                if pred_aa == target_aa:
                    # Apply filters
                    if pred['score'] < threshold:
                        continue
                    if yes_only and pred['answer'] != 'YES':
                        continue
                    
                    matches_found += 1
                    pkey = f"{gene}-{mutation_id}"  # e.g., "PMS2-A1622G"
                    
                    # Create result entry with pkey
                    results.append({
                        'pkey': pkey,
                        'Gene': pred['Gene'],
                        'pos': pred['pos'],
                        'amino_acid': pred['amino_acid'],
                        'context': pred['context'],
                        'score': pred['score'],
                        'kinase': pred['kinase'],
                        'answer': pred['answer']
                    })
    
    print(f"Mutation positions checked: {len(mapping_dict) - stop_codons_skipped}")
    print(f"Stop codons skipped: {stop_codons_skipped}")
    print(f"Predictions at mutation sites: {matches_found}")
    
    return results

def _run_native_netphos(fasta_file, output_file, timeout=300, ape_bin=None):
    """Run NetPhos using native APE binary (Linux only)."""
    if not ape_bin:
        return False, "No APE binary path provided"

    try:
        # APE needs to run from its installation directory
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
    """Process large FASTA files using batching to prevent segmentation faults"""
    if executor_fn is None:
        executor_fn = _run_native_netphos

    try:
        # Split FASTA into batches
        batch_files = split_fasta_into_batches(fasta_file, batch_size)

        if not batch_files:
            print("No sequences found in FASTA file")
            return False

        print(f"Processing {len(batch_files)} batches...")
        batch_outputs = []

        # Process each batch
        for i, batch_file in enumerate(batch_files):
            batch_output = output_file.replace('.out', f'-batch-{i+1}.out')
            print(f"Processing batch {i+1}/{len(batch_files)}...")

            success, error = executor_fn(batch_file, batch_output, timeout, ape_bin)
            
            if success:
                batch_outputs.append(batch_output)
                print(f"Batch {i+1} completed: {batch_output}")
            else:
                print(f"Batch {i+1} failed: {error}")
                # Continue with other batches
        
        # Handle batch outputs
        if batch_outputs:
            if len(batch_outputs) == 1:
                # Single batch: just rename/move the file instead of combining
                single_batch_output = batch_outputs[0]
                try:
                    shutil.move(single_batch_output, output_file)
                    print(f"Single batch completed: {output_file}")
                    return True
                except Exception as e:
                    print(f"Failed to move single batch output: {e}")
                    return False
            else:
                # Multiple batches: combine them
                success = combine_batch_outputs(batch_outputs, output_file, format_type='netphos')
                
                if success:
                    print(f"Combined {len(batch_outputs)} batch outputs into {output_file}")
                    
                    # Preserve individual batch files for parsing
                    print(f"Preserved {len(batch_outputs)} batch files for parsing:")
                    for batch_output in batch_outputs:
                        print(f"  {batch_output}")
                    
                    return True
                else:
                    print("Failed to combine batch outputs")
                    return False
        else:
            print("No successful batches to combine")
            return False
            
        # Clean up batch FASTA files
        for batch_file in batch_files:
            if os.path.exists(batch_file):
                os.remove(batch_file)
                
    except Exception as e:
        print(f"Batch processing failed: {e}")
        return False

def get_file_cache_key(fasta_file):
    """Generate cache key based on file path and modification time"""
    try:
        file_path = str(fasta_file)
        stat = os.stat(file_path)
        # Use file path, size, and mtime for cache key
        cache_data = f"{file_path}:{stat.st_size}:{stat.st_mtime}"
        return hashlib.md5(cache_data.encode()).hexdigest()
    except Exception:
        return None

def get_cache_dir():
    """Get or create cache directory"""
    cache_dir = os.path.expanduser("~/.netphos_cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir

def get_cached_result(fasta_file, cache_type="netphos"):
    """Check if cached result exists for file"""
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
            # Remove corrupted cache files
            for f in [cache_file, metadata_file]:
                if os.path.exists(f):
                    os.remove(f)
    
    return None, None

def save_to_cache(fasta_file, output_file, cache_type="netphos", metadata=None):
    """Save result to cache"""
    cache_key = get_file_cache_key(fasta_file)
    if not cache_key or not os.path.exists(output_file):
        return False
    
    try:
        cache_dir = get_cache_dir()
        cache_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.out")
        metadata_file = os.path.join(cache_dir, f"{cache_key}_{cache_type}.json")
        
        # Copy output to cache
        shutil.copy2(output_file, cache_file)
        
        # Save metadata
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
    """Clear all cached NetPhos results"""
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
    """Count sequences in FASTA file"""
    count = 0
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        return 1  # Fallback assumption
    return count

def run_netphos_with_fasta(fasta_file, output_file, batch_size=None, timeout=300, use_cache=True,
                           executor_fn=None, ape_bin=None):
    """Run NetPhos on FASTA file with intelligent processing strategy selection and caching"""
    if executor_fn is None:
        executor_fn = _run_native_netphos

    mode_label = "native"

    def _run(fasta, output, tmo):
        return executor_fn(fasta, output, tmo, ape_bin)

    # Check cache first if enabled
    if use_cache:
        cached_file, cached_metadata = get_cached_result(fasta_file, "netphos")
        if cached_file:
            print(f"Using cached NetPhos result for {fasta_file}")
            if cached_metadata:
                print(f"  Cached at: {cached_metadata.get('cached_at', 'Unknown')}")
            shutil.copy2(cached_file, output_file)
            return True

    # If batch_size is explicitly specified, use batching
    if batch_size:
        print(f"Using batch processing (batch_size: {batch_size})...")
        result = process_netphos_batched(fasta_file, output_file, batch_size, timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)

        # Cache result if successful
        if result and use_cache:
            seq_count = count_fasta_sequences(fasta_file)
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": batch_size, "sequence_count": seq_count})
        return result

    # Count sequences to determine optimal strategy
    seq_count = count_fasta_sequences(fasta_file)
    print(f"Processing {seq_count} sequence(s)...")

    # Intelligent strategy selection based on sequence count
    if seq_count == 1:
        # Single sequence: only try single run, no fallback
        print(f"Single sequence detected - using single {mode_label} run...")
        success, error = _run(fasta_file, output_file, timeout)

        if success:
            print(f"NetPhos completed: {output_file}")
            # Cache result if successful
            if use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "single", "sequence_count": seq_count})
            return True
        else:
            print(f"Single run failed: {error}")
            print(f"NetPhos failed for single sequence. Check {mode_label} setup and APE system availability.")
            return False

    elif seq_count <= 10:
        # Small batch (2-10 sequences): try single run first, then batch if needed
        print(f"Small sequence set - attempting single {mode_label} run...")
        success, error = _run(fasta_file, output_file, timeout)

        if success:
            print(f" NetPhos completed: {output_file}")
            # Cache result if successful
            if use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "single", "sequence_count": seq_count})
            return True
        else:
            print(f" Single run failed: {error}")
            print("Falling back to batch processing for small sequence set...")
            result = process_netphos_batched(fasta_file, output_file, batch_size=10, timeout=timeout,
                                             executor_fn=executor_fn, ape_bin=ape_bin)
            # Cache result if successful
            if result and use_cache:
                save_to_cache(fasta_file, output_file, "netphos",
                              {"processing_mode": "batch_fallback", "batch_size": 10, "sequence_count": seq_count})
            return result

    elif seq_count <= 100:
        # Medium batch (11-100 sequences): use batch processing with size 25
        print("Medium sequence set - using batch processing...")
        result = process_netphos_batched(fasta_file, output_file, batch_size=25, timeout=timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        # Cache result if successful
        if result and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 25, "sequence_count": seq_count})
        return result

    else:
        # Large batch (100+ sequences): use batch processing with size 50
        print("Large sequence set - using batch processing...")
        result = process_netphos_batched(fasta_file, output_file, batch_size=50, timeout=timeout,
                                         executor_fn=executor_fn, ape_bin=ape_bin)
        # Cache result if successful
        if result and use_cache:
            save_to_cache(fasta_file, output_file, "netphos",
                          {"processing_mode": "batch", "batch_size": 50, "sequence_count": seq_count})
        return result

def main():
    parser = argparse.ArgumentParser(description='Parse NetPhos output to TSV format')
    parser.add_argument('input', nargs='?', help='Input: FASTA file/directory (full-pipeline/netphos-only modes) or NetPhos output file/directory (parse-only mode)')
    parser.add_argument('output', nargs='?', help='Output TSV file')
    parser.add_argument('--threshold', type=float, default=0.0, 
                       help='Minimum score threshold (default: 0.0 - include all). Ignored if --yes-only is used.')
    parser.add_argument('--yes-only', action='store_true',
                       help='Only include predictions marked as YES (high-confidence predictions, overrides --threshold)')
    parser.add_argument('--mapping-dir', 
                       help='Directory containing mutation mapping CSV files (enables pkey generation and mutation filtering)')
    parser.add_argument('--is-mutant', action='store_true',
                       help='Process mutant sequences (matches mutant amino acid). Default: wildtype (matches original amino acid)')
    parser.add_argument('--log',
                       help='Validation log file or directory to skip failed mutations (mutant mode only)')
    
    # Mode selection
    parser.add_argument('--mode', choices=['full-pipeline', 'netphos-only', 'parse-only'], default='parse-only',
                       help='Processing mode: "full-pipeline" (run NetPhos + parse), "netphos-only" (run NetPhos only), "parse-only" (parse existing outputs). Default: parse-only')
    parser.add_argument('--batch-size', type=int, 
                       help='Batch size for large FASTA files (triggers batch processing)')
    parser.add_argument('--timeout', type=int, default=300,
                       help='Docker command timeout in seconds (default: 300)')
    
    # Cache options
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable result caching (always reprocess files)')
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
    
    # Validate required arguments when not clearing cache
    if not args.input or not args.output:
        parser.error("input and output arguments are required (unless using --clear-cache)")
    
    # Mode-specific validation
    if args.mode in ['full-pipeline', 'netphos-only']:
        # Validate inputs for NetPhos execution modes
        if not args.input or not args.output:
            parser.error(f"For {args.mode} mode: both input (FASTA file/directory) and output (TSV file) are required")
        
        if not os.path.exists(args.input):
            parser.error(f"Input FASTA file/directory not found: {args.input}")
            
        # Accept either single FASTA file or directory
        if os.path.isfile(args.input):
            if not args.input.endswith(('.fasta', '.fa', '.fas')):
                parser.error(f"For {args.mode} mode: input file must be a FASTA file (.fasta, .fa, .fas)")
        elif not os.path.isdir(args.input):
            parser.error(f"Input must be a FASTA file or directory: {args.input}")
    
    elif args.mode == 'parse-only':
        # Validate inputs for parse-only mode
        if not args.input or not args.output:
            parser.error(f"For {args.mode} mode: both input (NetPhos output file/directory) and output (TSV file) are required")
        
        if not os.path.exists(args.input):
            parser.error(f"Input NetPhos output file/directory not found: {args.input}")
        
        # Accept either single file or directory
        if os.path.isfile(args.input):
            if not args.input.endswith(('.out', '.txt')):
                print(f"Warning:  Warning: Input file {args.input} doesn't have typical NetPhos extension (.out, .txt)")
        elif not os.path.isdir(args.input):
            parser.error(f"Input must be a file or directory: {args.input}")
    
    # Validate threshold/yes-only logic
    if args.yes_only and args.threshold > 0.0:
        print("Warning:  Warning: --yes-only specified, ignoring --threshold value (YES answers are already high-confidence)")
        args.threshold = 0.0  # Reset threshold when yes-only is used
    
    failure_map = load_validation_failures(args.log) if (args.log and args.is_mutant) else {}

    # Resolve native APE binary
    native_ape = resolve_native_ape_path(getattr(args, 'native_ape_path', None))
    if not native_ape:
        parser.error(
            "Native APE binary not found. Provide --native-ape-path, or set "
            "NETPHOS_APE_PATH / NETPHOS_HOME environment variable."
        )
    executor_fn = _run_native_netphos
    ape_bin = native_ape
    print(f"Execution mode: native APE ({native_ape})")

    # Handle different processing modes
    if args.mode in ['full-pipeline', 'netphos-only']:
        # Input can be a FASTA file or directory containing FASTA files
        # Create temporary directory for NetPhos outputs (like NetNGlyc does)
        import tempfile
        temp_output_dir = tempfile.mkdtemp(prefix="netphos_outputs_")
        
        if os.path.isfile(args.input):
            # Single file processing - create output in temp directory
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
            
            print(f" NetPhos completed: {netphos_output}")
            
            # For netphos-only mode, we're done
            if args.mode == 'netphos-only':
                print("NetPhos execution completed. Use --mode parse-only to parse the output.")
                return 0
            
            # For full-pipeline mode, continue to parsing
            args.input = temp_output_dir
            print(f"Proceeding to parse outputs from {temp_output_dir}")
            
        elif os.path.isdir(args.input):
            # Directory processing - find all FASTA files using flexible discovery
            discovered_fastas = discover_fasta_files(args.input)
            fasta_files = [Path(file_path) for file_path in discovered_fastas.values()]
            
            if not fasta_files:
                print(f"Warning: No FASTA files found in directory: {args.input}")
                print(f"Searched for extensions: .fasta, .fa, .fas, .fna")
                return 1
            
            print(f"Found {len(fasta_files)} FASTA files in {args.input}")
            
            # Process each FASTA file individually - outputs go directly to temp directory
            netphos_outputs = []
            for fasta_file in fasta_files:
                fasta_path = str(fasta_file)
                seq_count = count_fasta_sequences(fasta_path)
                # Create output in temp directory with flexible extension handling
                input_basename = os.path.basename(fasta_path)
                # Remove any FASTA extension and add -netphos.out
                for ext in ['.fasta', '.fa', '.fas', '.fna']:
                    if input_basename.lower().endswith(ext):
                        input_basename = input_basename[:-len(ext)]
                        break
                netphos_output = os.path.join(temp_output_dir, f"{input_basename}-netphos.out")
                
                print(f"Processing {fasta_path} ({seq_count} sequences)...")
                use_cache = not args.no_cache
                # Let each file determine its own optimal processing strategy
                success = run_netphos_with_fasta(fasta_path, netphos_output, args.batch_size, args.timeout, use_cache,
                                                 executor_fn=executor_fn, ape_bin=ape_bin)

                if success:
                    netphos_outputs.append(netphos_output)
                    print(f" NetPhos completed: {netphos_output}")
                else:
                    print(f" NetPhos failed for: {fasta_path}")
            
            if not netphos_outputs:
                print("ERROR: NetPhos execution failed for all files")
                return 1
            
            print(f" NetPhos completed for {len(netphos_outputs)}/{len(fasta_files)} files")
            
            # For netphos-only mode, we're done
            if args.mode == 'netphos-only':
                print("NetPhos execution completed. Use --mode parse-only to parse the outputs.")
                return 0
            
            # For full-pipeline mode, continue to parsing the directory of outputs
            args.input = temp_output_dir
            print(f"Proceeding to parse outputs from {temp_output_dir}")
        else:
            print(f"ERROR: Input must be a FASTA file or directory: {args.input}")
            return 1
    
    # Parse mode (either full-pipeline continuing or parse-only)
    if args.mode in ['full-pipeline', 'parse-only']:
        # Validate mapping directory is required for parsing modes
        if not args.mapping_dir:
            parser.error(f"For {args.mode} mode: --mapping-dir is REQUIRED (directory containing mutation mapping CSV files)")
        
        if not os.path.exists(args.mapping_dir):
            print(f"ERROR: Mapping directory not found: {args.mapping_dir}")
            return 1
        
        # Use mutation filtering mode - generates pkeys (always required for parsing modes)
        sequence_type = "mutant" if args.is_mutant else "wildtype"
        print(f"Using mutation filtering mode with pkey generation (sequence_type: {sequence_type})...")
        filtered_predictions = parse_with_mutation_filtering_directory(
            args.input,
            args.mapping_dir,
            args.threshold,
            args.yes_only,
            args.is_mutant,
            failure_map=failure_map if args.is_mutant else None,
        )
        
        # Write output with pkeys
        write_tsv(filtered_predictions, args.output, include_pkey=True)
        print(f"Wrote {len(filtered_predictions)} predictions with pkeys to {args.output}")
    
    return 0

if __name__ == '__main__':
    exit(main())
