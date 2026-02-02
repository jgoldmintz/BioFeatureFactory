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
NetNGlyc Pipeline with Host SignalP 6 Integration
Runs SignalP 6 on host, then executes native NetNGlyc binary.
NetNGlyc analyzes FULL sequence with SignalP context.
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
from typing import Optional
from Bio.Seq import Seq


# Import utility functions
from biofeaturefactory.utils.utility import (
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


# ---------------------------------------------------------------------------
# Native execution support
# ---------------------------------------------------------------------------

def is_linux_host():
    """Return True when running on a Linux kernel."""
    return platform.system().lower() == "linux"


def resolve_native_netnglyc_path(user_path=None):
    """
    Resolve a usable native NetNGlyc executable when available.

    Search order:
    1. Explicit --native-netnglyc-bin value
    2. $NETNGLYC_PATH environment variable
    3. $NETNGLYC_HOME/netNglyc
    4. Common install locations
    """
    candidates = []

    def _add(path):
        if path:
            candidates.append(os.path.expanduser(path))

    _add(user_path)
    _add(os.environ.get("NETNGLYC_PATH"))

    netnglyc_home = os.environ.get("NETNGLYC_HOME")
    if netnglyc_home:
        _add(os.path.join(netnglyc_home, "netNglyc"))

    home = Path.home()
    common_roots = [
        home / "netNglyc" / "netNglyc",
        Path("/opt/netnglyc/netNglyc"),
        Path("/usr/local/bin/netNglyc"),
    ]

    for candidate in common_roots:
        _add(str(candidate))

    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)

    return None


def translate_orf_sequence(nt_sequence: str) -> str:
    """Translate a nucleotide ORF into an amino acid sequence, trimming trailing stops."""
    if not nt_sequence:
        return ""
    cleaned = nt_sequence.strip().upper().replace("U", "T")
    if not cleaned:
        return ""
    aa_seq = str(Seq(cleaned).translate(to_stop=False))
    return aa_seq.rstrip('*').strip()


def load_wt_sequence_map(input_path: str, wt_header: str = "ORF"):
    """
    Load WT nucleotide sequences from a file or directory using shared utility helpers.

    Returns:
        tuple(dict, tempfile.TemporaryDirectory | None): (gene -> nt sequence, temp dir holder)
    """
    source = Path(input_path)
    temp_dir = None

    if source.is_dir():
        load_path = str(source)
    elif source.is_file():
        temp_dir = tempfile.TemporaryDirectory(prefix="netnglyc_wt_src_")
        shutil.copy2(str(source), os.path.join(temp_dir.name, source.name))
        load_path = temp_dir.name
    else:
        raise FileNotFoundError(f"Input FASTA path not found: {input_path}")

    sequences = load_wt_sequences(load_path, wt_header=wt_header)
    return sequences, temp_dir


def infer_aamutation_from_nt(mutant_id: str, nt_sequence: str):
    """Infer the amino-acid mutation string (e.g., K543E) from a nucleotide mutation."""
    nt_info = get_mutation_data_bioAccurate(mutant_id)
    if nt_info[0] is None:
        return None
    aa_info = get_mutant_aa(nt_info, nt_sequence)
    if not aa_info:
        return None
    (aa_pos, (wt_aa, mut_aa)), _ = aa_info
    if not wt_aa or not mut_aa:
        return None
    aa_pos = int(aa_pos)
    return aa_pos, wt_aa, mut_aa


def build_mutant_sequences_for_gene(
    gene_name: str,
    nt_sequence: str,
    aa_sequence: str,
    mapping_file: Optional[str],
    log_path: Optional[str],
    failure_map: Optional[dict],
):
    """Return a dict of {header: sequence} for all mutants of a given gene."""
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
                    inferred = infer_aamutation_from_nt(mutant_clean, nt_sequence)
                    if inferred is None:
                        continue
                    pos, wt_aa, mut_aa = inferred

                idx = int(pos) - 1
                if idx < 0 or idx >= len(aa_sequence):
                    continue
                if wt_aa and aa_sequence[idx].upper() != wt_aa.upper():
                    # Skip if reference AA does not match translated ORF
                    continue

                header = f"{gene_name}-{mutant_clean}"
                mutant_sequences[header] = update_str(aa_sequence, mut_aa, idx)
    except Exception as exc:
        print(f"Warning: Failed to synthesize mutants for {gene_name} ({mapping_file}): {exc}")
        return {}

    return mutant_sequences


def synthesize_gene_fastas(wt_sequences, mapping_lookup, sequence_root, log_path=None, failure_map=None):
    """Create WT and mutant FASTAs for each gene, returning the directories and summary."""
    sequence_root = Path(sequence_root)
    wt_dir = sequence_root / "wt"
    mut_dir = sequence_root / "mut"
    wt_dir.mkdir(parents=True, exist_ok=True)
    mut_dir.mkdir(parents=True, exist_ok=True)

    summary = []
    for gene_name, nt_seq in wt_sequences.items():
        gene_name = gene_name.upper()
        nt_seq_upper = nt_seq.strip().upper()
        aa_seq = translate_orf_sequence(nt_seq_upper)
        if not aa_seq:
            print(f"Skipping {gene_name}: unable to translate ORF")
            continue

        wt_header = f"{gene_name}-wt"
        wt_path = wt_dir / f"{gene_name}-wt.fasta"
        write_fasta(wt_path, {wt_header: aa_seq})

        mapping_file = mapping_lookup.get(gene_name.upper())
        mutant_sequences = build_mutant_sequences_for_gene(
            gene_name,
            nt_seq_upper,
            aa_seq,
            mapping_file,
            log_path,
            failure_map,
        )

        mut_path = None
        if mutant_sequences:
            mut_path = mut_dir / f"{gene_name}_aa.fasta"
            write_fasta(mut_path, mutant_sequences)

        summary.append({
            "gene": gene_name,
            "wt_path": str(wt_path),
            "mut_path": str(mut_path) if mut_path else None,
            "mutant_count": len(mutant_sequences),
        })

    return wt_dir, mut_dir, summary


def normalize_mutation_id(mut_id: str) -> str:
    """Normalize mutation identifiers to a consistent format used in headers/pkeys."""
    return (mut_id or "").replace(" ", "").strip().upper()


def parse_signalp_summary(file_path: str):
    """Extract per-sequence SignalP info from NetNGlyc output footer."""
    summary = {}
    try:
        with open(file_path, "r") as handle:
            for line in handle:
                line = line.strip()
                if not line.startswith("# ") or "Signal peptide" not in line:
                    continue
                payload = line[2:]
                if ":" not in payload:
                    continue
                seq_id, desc = payload.split(":", 1)
                seq_id = seq_id.strip()
                desc = desc.strip()
                has_signal = "Signal peptide detected" in desc
                cleavage = None
                probability = None
                if "probability" in desc:
                    try:
                        prob_str = desc.split("probability:")[-1].strip(") ")
                        probability = float(prob_str)
                    except ValueError:
                        probability = None
                if "cleavage at position" in desc:
                    try:
                        cleavage_token = desc.split("cleavage at position", 1)[1]
                        cleavage = int(cleavage_token.split()[0])
                    except (ValueError, IndexError):
                        cleavage = None
                summary[seq_id] = {
                    "has_signal": has_signal,
                    "probability": probability,
                    "cleavage_site": cleavage,
                }
    except FileNotFoundError:
        pass
    return summary


def load_mutation_index(mapping_lookup):
    """Build quick lookup of expected mutations per gene from mapping CSV files."""
    mutation_index: dict[str, set[str]] = {}
    for gene, mapping_file in mapping_lookup.items():
        gene_key = gene.upper()
        mutation_index.setdefault(gene_key, set())
        try:
            with open(mapping_file, "r") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    mutant = (
                        row.get("mutant")
                        or row.get("mutation")
                        or row.get("nt_mutation")
                        or row.get("ntmutant")
                        or ""
                    ).strip()
                    if mutant:
                        mutation_index[gene_key].add(normalize_mutation_id(mutant))
        except OSError:
            continue
    return mutation_index


def normalize_gene_symbol(name: str) -> str:
    """Normalize a candidate gene token to uppercase canonical form."""
    cleaned = name.strip()
    cleaned = cleaned.replace("_wt", "").replace("-wt", "").replace("_WT", "").replace("-WT", "")
    gene = extract_gene_from_filename(cleaned) or cleaned
    return gene.upper()


def iter_netnglyc_output_files(directory: Path):
    """Yield NetNGlyc output files (combined and per-sequence) from a directory."""
    patterns = ["*-netnglyc.out", "seq_*_output.out"]
    for pattern in patterns:
        for file_path in directory.glob(pattern):
            if file_path.is_file():
                yield file_path


def load_signalp_cache(cache_root):
    """Load cached SignalP summaries from json files under the cache root."""
    cache_data = {}
    bases = []
    if cache_root:
        bases.append(Path(cache_root))
        bases.append(Path(cache_root) / ".signalp6_cache")
    bases.append(Path(os.path.expanduser("~")) / ".signalp6_cache")
    seen_dirs: set[Path] = set()
    for base in bases:
        if not base or not base.exists():
            continue
        for pred_file in base.glob("*_sp6_output/prediction_results.txt"):
            parent = pred_file.parent
            if parent in seen_dirs:
                continue
            seen_dirs.add(parent)
            try:
                with open(pred_file, "r") as handle:
                    for line in handle:
                        if not line or line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) < 4:
                            parts = line.strip().split()
                        if len(parts) < 4:
                            continue
                        seq_id = parts[0].strip()
                        prediction = parts[1].strip().upper()
                        sp_prob = None
                        try:
                            sp_prob = float(parts[3]) if parts[3] else None
                        except (ValueError, IndexError):
                            sp_prob = None
                        cleavage = None
                        if len(parts) > 4 and parts[4]:
                            token = parts[4].strip()
                            try:
                                if "-" in token:
                                    cleavage = int(token.split("-")[0])
                                else:
                                    cleavage = int(token)
                            except ValueError:
                                cleavage = None
                        cache_data[seq_id] = {
                            "has_signal": prediction.startswith("SP"),
                            "probability": sp_prob,
                            "cleavage_site": cleavage,
                        }
            except OSError:
                continue
    return cache_data

def _process_single_sequence_worker(args):
    """
    Worker function for parallel processing - must be at module level for pickling
    """
    (temp_fasta, temp_output, seq_name, use_signalp, cache_dir,
     docker_timeout, verbose, native_path) = args

    worker_processor = RobustDockerNetNGlyc(
        use_signalp=use_signalp,
        max_workers=1,
        cache_dir=cache_dir,
        docker_timeout=docker_timeout,
        verbose=verbose,
        native_bin=native_path,
    )

    try:
        success, output, error = worker_processor.process_single_fasta(
            temp_fasta, temp_output, 0
        )
        return success, temp_output, error, seq_name
    except Exception as e:
        return False, temp_output, str(e), seq_name
    finally:
        # Clean up worker processor
        if hasattr(worker_processor, 'temp_dir') and os.path.exists(worker_processor.temp_dir):
            shutil.rmtree(worker_processor.temp_dir, ignore_errors=True)


class SignalP6Handler:
    """Handle SignalP 6 predictions on the host before NetNGlyc execution"""

    def __init__(self, cache_dir=None, verbose=False, logger=None):
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".signalp6_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        self.last_output_dir = None  # Store the output directory path
        self.verbose = verbose
        self.logger = logger
        # Check SignalP availability once at initialization
        self.signalp6_available = self.check_signalp6_available()

    def check_signalp6_available(self):
        """Check if SignalP 6 is available on the host"""
        try:
            result = subprocess.run(
                ["signalp6", "--version"],
                capture_output=True,
                timeout=50,
                text=True
            )
            if result.returncode == 0:
                if self.verbose:
                    print(f"Using SignalP 6.0: {result.stdout.strip()}")
                return True
        except:
            pass
        return False

    def get_cache_key(self, fasta_file):
        """Generate cache key for SignalP results"""
        hasher = hashlib.md5()
        with open(fasta_file, 'rb') as f:
            hasher.update(f.read())
        return hasher.hexdigest()

    def clear_cache(self):
        """Clear all cached SignalP results"""
        if os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
            os.makedirs(self.cache_dir, exist_ok=True)
            return True
        return False

    def run_signalp6(self, fasta_file, output_dir=None):
        """
        Run SignalP 6 on the host
        Returns: (dict of results, path to output directory)
        """
        # Check cache first
        cache_key = self.get_cache_key(fasta_file)
        cache_file = os.path.join(self.cache_dir, f"{cache_key}_sp6.json")
        cache_dir = os.path.join(self.cache_dir, f"{cache_key}_sp6_output")

        if os.path.exists(cache_file) and os.path.exists(cache_dir):
            if self.verbose:
                print(f"Using cached SignalP 6 results for {os.path.basename(fasta_file)}")
            with open(cache_file, 'r') as f:
                return json.load(f), cache_dir

        if not self.signalp6_available:
            error_msg = "ERROR: SignalP 6.0 is required but not available. Please install SignalP 6.0 and ensure it's in your PATH."
            if self.logger:
                self.logger.error(error_msg)
            print(error_msg, file=sys.stderr)
            raise Exception("SignalP 6.0 not available - required for NetNGlyc processing with SignalP integration")

        if self.verbose:
            print(f"Running SignalP 6.0 on {os.path.basename(fasta_file)}")

        # Use provided output_dir or create temp
        if output_dir:
            signalp_output_dir = output_dir
            temp_dir = None
        else:
            temp_dir = tempfile.mkdtemp()
            signalp_output_dir = temp_dir

        results = {}

        try:
            # Run SignalP 6
            cmd = [
                "signalp6",
                "--fastafile", fasta_file,
                "--output_dir", signalp_output_dir,
                "--organism", "eukarya",
                "--format", "txt",
                "--mode", "fast"
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode == 0:
                # Parse prediction_results.txt
                pred_file = os.path.join(signalp_output_dir, "prediction_results.txt")
                if os.path.exists(pred_file):
                    with open(pred_file, 'r') as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) >= 9:
                                seq_id = parts[0]
                                prediction = parts[1]  # 'SP' or 'OTHER'
                                sp_prob = float(parts[3]) if parts[3] else 0.0

                                # Parse CS Position (column 9)
                                cs_pos = None
                                if len(parts) > 8 and parts[8] and parts[8].strip():
                                    cs_info = parts[8].strip()
                                    try:
                                        if '-' in cs_info:
                                            cs_pos = int(cs_info.split('-')[0])
                                        elif cs_info.isdigit():
                                            cs_pos = int(cs_info)
                                    except:
                                        cs_pos = 25

                                has_sp = prediction == 'SP'
                                if has_sp and not cs_pos:
                                    cs_pos = 25  # Default

                                results[seq_id] = {
                                    'has_signal': has_sp,
                                    'cleavage_site': cs_pos if has_sp else None,
                                    'probability': sp_prob
                                }

                # Cache results and output directory
                with open(cache_file, 'w') as f:
                    json.dump(results, f)

                # Copy output directory to cache
                if os.path.exists(cache_dir):
                    shutil.rmtree(cache_dir)
                shutil.copytree(signalp_output_dir, cache_dir)

                self.last_output_dir = signalp_output_dir
                return results, signalp_output_dir
            else:
                print(f"SignalP 6 failed: {result.stderr[:200]}")

        except Exception as e:
            if self.logger:
                self.logger.error(f"SignalP 6 error: {e}")
            else:
                print(f"SignalP 6 error: {e}")
        finally:
            # Clean up temp dir if one was created
            if temp_dir and temp_dir != signalp_output_dir:
                shutil.rmtree(temp_dir)

        return results, None


class RobustDockerNetNGlyc:
    """
    NetNGlyc processor with integrated SignalP 6 preprocessing on host.

    Requires:
    - NetNGlyc binary installed on the host (native execution)
    - SignalP 6.0 installed on host system (when use_signalp=True)
    """

    def __init__(self, use_signalp=True,
                 max_workers=4, cache_dir=None, docker_timeout=600, keep_intermediates=False,
                 verbose=False, native_bin=None, **_ignored):
        self.max_workers = max_workers
        self.temp_dir = tempfile.mkdtemp(prefix="netnglyc_")
        self.use_signalp = use_signalp
        self.docker_timeout = docker_timeout
        self.keep_intermediates = keep_intermediates
        self.verbose = verbose
        self.native_bin = native_bin

        # Initialize SignalP handler
        if use_signalp:
            self.signalp_handler = SignalP6Handler(
                cache_dir=cache_dir, verbose=self.verbose, logger=getattr(self, "error_logger", None)
            )
        else:
            self.signalp_handler = None

        # Initialize cache
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".netnglyc_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        # Setup error logging with date-named log file
        self._setup_error_logging()

        # Resolve native binary path
        self._setup_execution_mode()

    def _setup_execution_mode(self):
        """Resolve native NetNGlyc binary path."""
        native_path = resolve_native_netnglyc_path(self.native_bin)
        if not native_path:
            raise ValueError(
                "Native NetNGlyc binary not found. Provide native_bin argument, "
                "or set NETNGLYC_PATH / NETNGLYC_HOME environment variable."
            )
        self.use_native = True
        self.native_path = native_path

        if self.verbose:
            print(f"Execution mode: native ({self.native_path})")

    def _setup_error_logging(self):
        """Setup error logging with date-named log file"""
        # Create error log with current date
        today = datetime.now().strftime("%Y-%m-%d")
        log_filename = f"netnglyc_errors_{today}.log"
        log_path = os.path.join(self.cache_dir, log_filename)
        
        # Setup logger
        self.error_logger = logging.getLogger(f'netnglyc_errors_{id(self)}')
        self.error_logger.setLevel(logging.ERROR)
        
        # Remove existing handlers to avoid duplicates
        for handler in self.error_logger.handlers[:]:
            self.error_logger.removeHandler(handler)
        
        # Create file handler
        handler = logging.FileHandler(log_path)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.error_logger.addHandler(handler)
        
        # Don't propagate to root logger
        self.error_logger.propagate = False

    def clear_cache(self):
        """Clear all cached NetNGlyc results"""
        if os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
            os.makedirs(self.cache_dir, exist_ok=True)
            return True
        return False

    def _run_native_netnglyc(self, fasta_file):
        """
        Run NetNGlyc using native Linux binary.

        Args:
            fasta_file: Path to input FASTA file

        Returns:
            NetNGlyc output as string

        Raises:
            Exception: If native execution fails
        """
        try:
            result = subprocess.run(
                [self.native_path, fasta_file],
                capture_output=True,
                text=True,
                timeout=self.docker_timeout,
                cwd=os.path.dirname(self.native_path)
            )

            if result.returncode == 0:
                if "Predictions for N-Glycosylation sites" in result.stdout:
                    return result.stdout
                elif "No Asparagines in the input sequences" in result.stdout:
                    return result.stdout

            raise Exception(f"Native NetNGlyc failed with return code {result.returncode}: {result.stderr}")

        except subprocess.TimeoutExpired:
            error_msg = f"Native NetNGlyc timeout ({self.docker_timeout}s) for file: {os.path.basename(fasta_file)}"
            print(f"ERROR: {error_msg}")
            self.error_logger.error(error_msg)
            raise Exception(f"Native NetNGlyc timeout after {self.docker_timeout} seconds")
        except Exception as e:
            error_msg = f"Native NetNGlyc failed for file {os.path.basename(fasta_file)}: {e}"
            print(f"ERROR: {error_msg}")
            self.error_logger.error(error_msg)
            raise Exception(f"Native NetNGlyc failed: {e}")

    def process_single_fasta(self, fasta_file, output_file, worker_id=0):
        """
        Process a single FASTA file with SignalP 6 preprocessing and licensed NetNGlyc

        1. Run SignalP 6 on host (if enabled)
        2. Run native NetNGlyc on ORIGINAL FASTA
        3. NetNGlyc analyzes full sequence with SignalP context
        """
        #print(f"Worker {worker_id}: Processing {os.path.basename(fasta_file)}")

        # Check cache
        cache_key = hashlib.md5(open(fasta_file, 'rb').read()).hexdigest()[:16]
        cache_file = os.path.join(self.cache_dir, f"{cache_key}_netnglyc.out")

        if os.path.exists(cache_file):
            #print(f"Worker {worker_id}: Using cached result")
            shutil.copy(cache_file, output_file)
            return True, output_file, None

        try:
            # Step 1: Run SignalP 6 on host (if enabled)
            signalp_results = {}
            signalp_output_dir = None

            if self.use_signalp and self.signalp_handler:
                # Create a temp directory for SignalP output
                sp_temp_dir = os.path.join(self.temp_dir, f"signalp_{worker_id}_{os.getpid()}")
                os.makedirs(sp_temp_dir, exist_ok=True)

                # Run SignalP 6 and get output directory
                signalp_results, signalp_output_dir = self.signalp_handler.run_signalp6(
                    fasta_file, sp_temp_dir
                )

                if signalp_results:
                    for seq_id, info in signalp_results.items():
                        if info['has_signal']:
                            pass  # SignalP info processed elsewhere
    
            # Step 2: Run NetNGlyc with ORIGINAL FASTA
            netnglyc_output = self._run_native_netnglyc(fasta_file)

            if netnglyc_output:
                # Step 3: Save results
                with open(output_file, 'w') as f:
                    f.write(netnglyc_output)

                    # Add SignalP information summary at the end
                    if signalp_results:
                        f.write("\n" + "=" * 70 + "\n")
                        f.write("# SignalP 6 predictions summary:\n")
                        for seq_id, info in signalp_results.items():
                            if info['has_signal']:
                                f.write(f"# {seq_id}: Signal peptide detected, "
                                        f"cleavage at position {info['cleavage_site']} "
                                        f"(probability: {info['probability']:.3f})\n")
                            else:
                                f.write(f"# {seq_id}: No signal peptide detected "
                                        f"(probability: {info['probability']:.6f})\n")

                # Cache the result
                shutil.copy(output_file, cache_file)

                # Parse and report glycosylation sites
                sites = self._parse_netnglyc_output(netnglyc_output)
                #if sites:
                    #print(f"Worker {worker_id}: Found {len(sites)} glycosylation site(s)")
                #else:
                    #print(f"Worker {worker_id}: No glycosylation sites found")

                return True, output_file, None
            else:
                return False, output_file, "NetNGlyc failed to produce output"

        except Exception as e:
            return False, output_file, str(e)
        finally:
            # Clean up SignalP output directory if it exists in temp
            if signalp_output_dir and signalp_output_dir.startswith(self.temp_dir):
                if os.path.exists(signalp_output_dir):
                    shutil.rmtree(signalp_output_dir)

    def _parse_netnglyc_output(self, output):
        """Parse NetNGlyc output to extract predicted sites"""
        sites = []
        in_results = False

        for line in output.split('\n'):
            if "SeqName" in line and "Position" in line:
                in_results = True
                continue

            if in_results and line.strip() and not line.startswith('#') and not line.startswith('-'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        # NetNGlyc output format: SeqName Position Motif Potential ...
                        position = int(parts[1]) if parts[1].isdigit() else int(parts[1].split()[0])
                        potential = 0.0

                        # Find the potential value (usually 4th column but can vary)
                        for part in parts[2:]:
                            try:
                                val = float(part)
                                if 0 <= val <= 1:
                                    potential = val
                                    break
                            except:
                                continue

                        if potential > 0.5:  # Threshold
                            sites.append({
                                'position': position,
                                'sequence': parts[2] if len(parts) > 2 and not parts[2].replace('.',
                                                                                                '').isdigit() else '',
                                'potential': potential
                            })
                    except (ValueError, IndexError):
                        pass

        return sites

    def load_nt_to_aa_mapping(self, mapping_file):
        """Load NT to AA mutation mapping from CSV file"""
        import csv
        mapping = {}
        try:
            with open(mapping_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    mutant = row['mutant']  # e.g., A1002T
                    aamutant = row['aamutant']  # e.g., V334V
                    
                    # Extract AA position from aamutant (e.g., V334V -> 334)
                    import re
                    pos_match = re.search(r'(\d+)', aamutant)
                    if pos_match:
                        aa_pos = int(pos_match.group(1))
                        mapping[mutant] = aa_pos
        except Exception as e:
            print(f"Error loading mapping file {mapping_file}: {e}")
        return mapping
    
    def load_nt_to_aa_mapping_enhanced(self, mapping_file):
        """Load NT to AA mutation mapping with amino acid data for 1:many logic"""
        import csv
        mapping = {}
        try:
            with open(mapping_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    mutant = row['mutant']  # e.g., A1002T
                    aamutant = row['aamutant']  # e.g., K541E
                    
                    # Parse amino acid mutation to get position and amino acids
                    position_data = get_mutation_data_bioAccurate(aamutant)
                    if position_data[0] is not None:
                        aa_pos = position_data[0]  # e.g., 541
                        aa_tuple = position_data[1]  # e.g., ('K', 'E')
                        mapping[mutant] = {
                            'aa_pos': aa_pos,
                            'original_aa': aa_tuple[0],
                            'mutant_aa': aa_tuple[1],
                            'aamutant': aamutant
                        }
        except Exception as e:
            print(f"Error loading enhanced mapping file {mapping_file}: {e}")
        return mapping

    def parse_netnglyc_output(self, file_path, threshold=0.5):
        """Parse NetNGlyc output file and extract all predictions"""
        predictions = []
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
                
            in_results = False
            skip_next_line = False
            found_first_separator = False
            for line in content.split('\n'):
                # Find the results table - handle both correct and malformed formats
                if ("SeqName" in line and "Position" in line and "Potential" in line) or \
                   ("agreement" in line and "result" in line and len(line.split()) <= 3):
                    in_results = True
                    # Only skip next line if this is the main header (not the subheader)
                    if "SeqName" in line:
                        skip_next_line = True  # Skip the second header line
                    continue
                
                # Skip the second header line  
                if skip_next_line:
                    skip_next_line = False
                    continue
                    
                # Handle separator lines - first one starts data, second one ends it
                # Support both --- and === separators for different NetNGlyc formats
                if in_results and (line.strip().startswith('---') or line.strip().startswith('===')):
                    if not found_first_separator:
                        found_first_separator = True
                        continue  # Continue after first separator to read data
                    else:
                        break  # Stop at second separator
                    
                if in_results and line.strip() and not line.startswith('#'):
                    # Parse NetNGlyc table format - handle tab-separated data
                    # Replace multiple tabs/spaces with single space for easier parsing
                    cleaned_line = ' '.join(line.split())
                    parts = cleaned_line.split()
                    if len(parts) >= 5:
                        try:
                            seq_name = parts[0]
                            position = int(parts[1])
                            sequon = parts[2]  # e.g., NKSE
                            potential = float(parts[3])
                            jury = parts[4] if len(parts) > 4 else ""
                            result = parts[5] if len(parts) > 5 else ""
                            
                            if potential >= threshold:
                                predictions.append({
                                    'seq_name': seq_name,
                                    'position': position,
                                    'sequon': sequon,
                                    'potential': potential,
                                    'jury_agreement': jury,
                                    'n_glyc_result': result
                                })
                        except (ValueError, IndexError) as e:
                            continue
                            
        except Exception as e:
            print(f"Error parsing NetNGlyc output {file_path}: {e}")
            
        return predictions


    def process_fasta_batched(self, fasta_file, output_base, batch_size=100, worker_id=0):
        """Process a FASTA file in batches to handle large numbers of sequences"""
        try:
            # Split FASTA into batches using shared utility (save to outputs directory for troubleshooting)
            outputs_dir = os.path.dirname(output_base)
            batch_files = split_fasta_into_batches(fasta_file, batch_size, temp_dir=outputs_dir)
            
            if not batch_files:
                return False, None, "Failed to create batch files"
            
            # If only one batch, use regular processing
            if len(batch_files) == 1:
                #print(f"Worker {worker_id}: Single batch processing")
                success, output_file, error = self.process_single_fasta(
                    batch_files[0], output_base, worker_id
                )
                # Clean up batch file (only if not keeping intermediates)
                if not self.keep_intermediates and os.path.exists(batch_files[0]):
                    os.remove(batch_files[0])
                return success, output_file, error
            
            # Process multiple batches
            #print(f"Worker {worker_id}: Processing {len(batch_files)} batches")
            batch_outputs = []
            
            for i, batch_file in enumerate(batch_files):
                # Create batch-specific output file name for parsing
                batch_output = output_base.replace('.out', f'-batch-{i+1}.out')
                
                #print(f"Worker {worker_id}: Processing batch {i+1}/{len(batch_files)}")
                
                try:
                    # Process this batch
                    success, batch_result, error = self.process_single_fasta(
                        batch_file, batch_output, worker_id
                    )
                    
                    if success:
                        batch_outputs.append(batch_output)
                        if self.verbose:
                            print(f"   Batch {i+1} saved: {os.path.basename(batch_output)}")
                    else:
                        if self.verbose:
                            print(f"Batch {i+1} failed: {error}")
                        # Clean up batch files (only if not keeping intermediates)
                        if not self.keep_intermediates:
                            for bf in batch_files:
                                if os.path.exists(bf):
                                    os.remove(bf)
                        return False, output_base, f"Batch {i+1} failed: {error}"
                        
                finally:
                    # Clean up this batch file (only if not keeping intermediates)
                    if not self.keep_intermediates and os.path.exists(batch_file):
                        os.remove(batch_file)
            
            # Create combined output for backwards compatibility (optional)
            #print(f"Worker {worker_id}: Creating combined output for compatibility")
            combined_output = combine_batch_outputs(batch_outputs, output_base, format_type='netnglyc', original_fasta_file=fasta_file)
            
            # PRESERVE all batch files for parsing - do NOT delete them
            #print(f"   Preserved {len(batch_outputs)} batch files for parsing:")
            for batch_output in batch_outputs:
                if os.path.exists(batch_output):
                    print(f"      - {os.path.basename(batch_output)}")
            if combined_output:
                return True, output_base, None
            else:
                return False, output_base, "Failed to combine batch outputs"
                
        except Exception as e:
            # Clean up any remaining batch files (only if not keeping intermediates)
            try:
                if 'batch_files' in locals() and not self.keep_intermediates:
                    for bf in batch_files:
                        if os.path.exists(bf):
                            os.remove(bf)
            except:
                pass
            return False, output_base, f"Batch processing error: {e}"


    def parse_netnglyc_multisequence_output(self, file_path, threshold=0.5):
        """Parse NetNGlyc output file containing multiple sequences"""
        predictions_by_sequence = {}
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            
            current_sequence = None
            in_results = False
            skip_next_line = False
            found_first_separator = False
            found_results_header = False
            line_count = 0
            results_lines_found = 0
            
            lines = content.split('\n')
            
            for i, line in enumerate(lines):
                line = line.strip()
                line_count += 1
                
                # Detect sequence names in the output
                if line.startswith('Name:') and len(line.split()) >= 2:
                    current_sequence = line.split()[1]
                    predictions_by_sequence[current_sequence] = []
                    in_results = False
                    found_first_separator = False
                    continue
                
                # Find the results table headers
                if "SeqName" in line and "Position" in line and "Potential" in line:
                    in_results = True
                    found_results_header = True
                    skip_next_line = True  # Skip the second header line
                    continue
                
                # Skip the second header line  
                if skip_next_line:
                    skip_next_line = False
                    continue
                
                # Handle separator lines - first one starts data, second one ends it
                if in_results and line.startswith('---'):
                    if not found_first_separator:
                        found_first_separator = True
                        continue  # Continue after first separator to read data
                    else:
                        in_results = False  # Stop at second separator
                        found_first_separator = False
                        continue
                
                # Parse prediction lines
                if in_results and found_first_separator and line and not line.startswith('#'):
                    results_lines_found += 1
                    
                    # Check for "No sites predicted" message
                    if "No sites predicted" in line:
                        continue
                    
                    # Parse NetNGlyc table format
                    cleaned_line = ' '.join(line.split())
                    parts = cleaned_line.split()
                    
                    if len(parts) >= 5:
                        try:
                            seq_name = parts[0]
                            if seq_name not in predictions_by_sequence:
                                predictions_by_sequence[seq_name] = []
                            position = int(parts[1])
                            sequon = parts[2]  # e.g., NKSE
                            potential = float(parts[3])
                            jury = parts[4] if len(parts) > 4 else ""
                            result = parts[5] if len(parts) > 5 else ""
                            
                            
                            if potential >= threshold:
                                # If the current sequence has yet to be detected, use seq_name
                                if current_sequence is None:
                                    current_sequence = seq_name
                                    predictions_by_sequence[current_sequence] = []
                                
                                prediction = {
                                    'seq_name': seq_name,
                                    'position': position,
                                    'sequon': sequon,
                                    'potential': potential,
                                    'jury_agreement': jury,
                                    'n_glyc_result': result
                                }
                                
                                # Add to the appropriate sequence's predictions
                                if seq_name in predictions_by_sequence:
                                    predictions_by_sequence[seq_name].append(prediction)
                                else:
                                    predictions_by_sequence[seq_name] = [prediction]
                            else:
                                pass  # Prediction below threshold
                        except (ValueError, IndexError) as e:
                            continue
            
            
            # Parsing completed
            pass
        except Exception as e:
            print(f"Error parsing multi-sequence NetNGlyc output {file_path}: {e}")
            
        return predictions_by_sequence

    def parse_simple_netnglyc_output(self, file_path, threshold=0.5):
        """Simple parser for NetNGlyc output - extracts predictions from results table"""
        predictions = []
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Find results table section
            lines = content.split('\n')
            in_results = False
            found_header = False
            
            for line in lines:
                line = line.strip()
                
                # Find results table header
                if "SeqName" in line and "Position" in line and "Potential" in line:
                    found_header = True
                    continue
                
                # Skip header continuation and separators
                if found_header and (line.startswith('---') or 'agreement' in line):
                    if line.startswith('---') and not in_results:
                        in_results = True  # Start reading predictions
                    elif line.startswith('---') and in_results:
                        break  # End of predictions
                    continue
                
                # Parse prediction lines
                if in_results and line and not line.startswith('#'):
                    if "No sites predicted" in line:
                        break
                        
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            seq_name = parts[0]
                            position = int(parts[1])
                            sequon = parts[2] if len(parts) > 2 else ""
                            potential = float(parts[3])
                            jury = parts[4] if len(parts) > 4 else ""
                            result = parts[5] if len(parts) > 5 else ""
                            
                            if potential >= threshold:
                                predictions.append({
                                    'seq_name': seq_name,
                                    'position': position,
                                    'sequon': sequon,
                                    'potential': potential,
                                    'jury_agreement': jury,
                                    'n_glyc_result': result
                                })
                        except (ValueError, IndexError):
                            continue
            
        except Exception as e:
            print(f"Error in simple NetNGlyc parsing: {e}")
        
        return predictions

    def parse_mutant_files(self, input_dir, threshold=0.5, mapping_dir=None, fasta_dir=None, failure_map=None):
        """Parse mutant NetNGlyc files using sequence names from predictions"""
        import os
        import csv
        import logging
        from datetime import datetime
        from pathlib import Path
        
        results = []
        failure_map = failure_map or {}
        wt_sequences = {}
        wt_temp_holder = None
        
        if fasta_dir:
            try:
                wt_sequences, wt_temp_holder = load_wt_sequence_map(fasta_dir, wt_header="ORF")
            except Exception as exc:
                if self.verbose:
                    print(f"Warning: unable to load WT sequences for mutant inference: {exc}")
                wt_sequences = {}
        
        if mapping_dir is None:
            print("ERROR: mapping_dir is required for mutant parsing")
            if wt_temp_holder:
                wt_temp_holder.cleanup()
            return results
        
        # Find NetNGlyc output files: both regular and batch files
        # Pattern 1: {GENE}_aa-netnglyc.out (regular files)
        regular_files = list(Path(input_dir).glob("*_aa-netnglyc.out"))
        # Pattern 2: {GENE}_aa-netnglyc-batch-*.out (batch files)
        batch_files = list(Path(input_dir).glob("*_aa-netnglyc-batch-*.out"))
        
        # Group files by gene for batch combination
        gene_files = {}
        
        # Build gene->files map by parsing names directly (avoid reading .out as FASTA)
        for file_path in regular_files:
            filename = file_path.stem  # Remove .out extension
            if '_aa-netnglyc' in filename:
                gene = extract_gene_from_filename(filename)
                if gene:
                    gene = gene.upper()
                    gene_files.setdefault(gene, []).append(file_path)
        
        # Process batch files using the same filename-based gene extraction
        for file_path in batch_files:
            filename = file_path.stem  # Remove .out extension
            if '_aa-netnglyc-batch-' in filename:
                gene = extract_gene_from_filename(filename)
                if gene:
                    gene = gene.upper()
                    gene_files.setdefault(gene, []).append(file_path)
        
        #print(f"Found files for {len(gene_files)} genes:")
        for gene, files in gene_files.items():
            file_types = []
            regular_count = sum(1 for f in files if '-batch-' not in f.stem)
            batch_count = sum(1 for f in files if '-batch-' in f.stem)
            if regular_count > 0:
                file_types.append(f"{regular_count} regular")
            if batch_count > 0:
                file_types.append(f"{batch_count} batch")
            print(f"  {gene}: {', '.join(file_types)} files")
        
        # Simplified processing: combine batches then process single file per gene
        gene_count = 0
        total_genes = len(gene_files)
        
        for gene, gene_file_list in gene_files.items():
            gene_count += 1
            
            if self.verbose:
                print(f"\nProcessing {gene} using simplified single-file approach ({len(gene_file_list)} source files)")
            else:
                print(f"Processing {gene}: {gene_count}/{total_genes}")
            
            # Step 1: If there are batch files, combine them into a single file
            combined_file_path = None
            if len(gene_file_list) > 1:
                # Check batch files needs combining
                batch_files_for_gene = [f for f in gene_file_list if '-batch-' in f.stem]
                if batch_files_for_gene:
                    if self.verbose:
                        print(f"  Combining {len(batch_files_for_gene)} batch files into single file")
                    combined_file_path = os.path.join(input_dir, f"{gene}_aa-netnglyc-combined.out")
                    try:
                        # Try to find the original FASTA file for this gene
                        possible_fasta_files = [
                            os.path.join(input_dir, f"{gene}_aa.fasta"),
                            os.path.join(input_dir, f"debug-{gene}-aa.fasta"),
                            os.path.join(input_dir, f"{gene}.fasta"),
                        ]
                        original_fasta_file = None
                        for fasta_path in possible_fasta_files:
                            if os.path.exists(fasta_path):
                                original_fasta_file = fasta_path
                                break
                        
                        success = combine_batch_outputs(
                            [str(f) for f in batch_files_for_gene], 
                            combined_file_path, 
                            format_type='netnglyc',
                            original_fasta_file=original_fasta_file
                        )
                        if success:
                            if self.verbose:
                                print(f"  Combined file created: {combined_file_path}")
                        else:
                            if self.verbose:
                                print(f"  Failed to combine batch files, using first file")
                            combined_file_path = str(gene_file_list[0])
                    except Exception as e:
                        if self.verbose:
                            print(f"  Error combining batch files: {e}, using first file")
                        combined_file_path = str(gene_file_list[0])
                else:
                    # Use the single regular file
                    combined_file_path = str(gene_file_list[0])
            else:
                # Single file, use as-is
                combined_file_path = str(gene_file_list[0])
            
            # Step 2: Parse the single combined file using simplified logic
            if self.verbose:
                print(f"  Parsing combined file: {os.path.basename(combined_file_path)}")
            
            # Load mapping file for this gene using flexible discovery
            discovered_mappings = discover_mapping_files(mapping_dir)
            if gene not in discovered_mappings:
                if self.verbose:
                    print(f"  Warning: No mapping file found for {gene} in {mapping_dir}")
                    print(f"  Available mappings: {list(discovered_mappings.keys())}")
                continue
            
            mapping_file = discovered_mappings[gene]
            
            # Read mapping file
            mapping_dict = {}
            try:
                wt_nt_seq = wt_sequences.get(gene.upper(), "").upper() if wt_sequences else ""
                with open(mapping_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        ntposnt = (
                            row.get('mutant')
                            or row.get('mutation')
                            or row.get('nt_mutation')
                            or row.get('ntmutant')
                            or ""
                        ).strip()
                        if not ntposnt:
                            continue
                        aaposaa = (
                            row.get('aamutant')
                            or row.get('aa_mutation')
                            or row.get('amino_acid_mutation')
                            or row.get('protein_mutation')
                            or ""
                        ).strip()
                        if not aaposaa and wt_nt_seq:
                            inferred = infer_aamutation_from_nt(ntposnt, wt_nt_seq)
                            if inferred:
                                aa_pos, wt_aa, mut_aa = inferred
                                aaposaa = f"{wt_aa}{aa_pos}{mut_aa}"
                        if not aaposaa:
                            continue
                        mapping_dict[ntposnt] = aaposaa
                if failure_map:
                    skip_set = {m.upper() for m in failure_map.get(gene.upper(), set())}
                    if skip_set:
                        before = len(mapping_dict)
                        mapping_dict = {
                            mut_id: aa_mut
                            for mut_id, aa_mut in mapping_dict.items()
                            if mut_id.upper() not in skip_set
                        }
                        skipped = before - len(mapping_dict)
                        if skipped and self.verbose:
                            print(f"  Skipping {skipped} mutations flagged in validation logs")
                if self.verbose:
                    print(f"  Loaded mapping for {len(mapping_dict)} mutations")
            except Exception as e:
                print(f"  Error reading mapping file {mapping_file}: {e}")
                continue
            
            if not mapping_dict:
                if self.verbose:
                    print("  No remaining mutations to evaluate after filtering; skipping gene")
                continue
            
            # Step 3: Parse predictions from the combined file
            try:
                all_predictions = self.parse_simple_netnglyc_output(combined_file_path, threshold)
                if self.verbose:
                    print(f"  Parsed {len(all_predictions)} predictions from combined file")
                
                if not all_predictions:
                    continue
                
                # Step 4: Use shared utility function for single-mutation processing
                try:
                    gene_results = parse_predictions_with_mutation_filtering(
                        all_predictions,
                        mapping_dict,
                        is_mutant=True,
                        threshold=threshold,
                        yes_only=False,
                        tool_type='netnglyc',
                        failure_map=failure_map,
                    )
                    if self.verbose:
                        print(f"  Found {len(gene_results)} predictions at mutation sites")
                    results.extend(gene_results)
                    
                except Exception as e:
                    print(f"  Error in unified mutation processing: {e}")
                    # Fall back to original logic if needed
                    print(f"  Falling back to original processing logic")
                    continue
                    
            except Exception as e:
                print(f"  Error parsing NetNGlyc output {combined_file_path}: {e}")
                continue
        
        if wt_temp_holder:
            wt_temp_holder.cleanup()
        return results

    def parse_wildtype_files(self, input_dir, mapping_dir, threshold=0.5):
        """Parse wildtype NetNGlyc files: {GENE}_aa-netnglyc.out or {GENE}_wt-netnglyc.out"""
        import os
        from pathlib import Path
        
        results = []
        # Look for both _aa and _wt patterns to handle different naming conventions
        output_files = list(Path(input_dir).glob("*-netnglyc.out"))
        
        for file_path in output_files:
            filename = file_path.stem  # Remove .out extension
            
            # Extract gene from file content using robust extraction
            if '_aa-netnglyc' in filename or '_wt-netnglyc' in filename:
                gene = ExtractGeneFromFASTA(str(file_path))
                if not gene:
                    # Skip files where gene extraction failed
                    continue
            else:
                # Skip files that don't match expected patterns
                continue
            
            # Load enhanced mapping for this gene with amino acid data using flexible discovery
            discovered_mappings = discover_mapping_files(mapping_dir)
            if gene in discovered_mappings:
                mapping_file = discovered_mappings[gene]
                mapping = self.load_nt_to_aa_mapping_enhanced(mapping_file)
                
                # Parse NetNGlyc predictions once
                predictions = self.parse_netnglyc_output(file_path, threshold)
                
                # Pre-process mapping data to create lookup structures for 1:many logic
                position_mutations = {}  # pos -> list of (mutation_id, mutation_data)
                stop_codons_skipped = 0
                
                for mutation_id, mutation_data in mapping.items():
                    # Skip if parsing failed
                    if not mutation_data:
                        stop_codons_skipped += 1
                        continue
                    
                    aa_pos = mutation_data['aa_pos']
                    # For wildtype processing, match against original amino acid
                    target_aa = mutation_data['original_aa']
                    
                    # Group mutations by position for efficient lookup
                    if aa_pos not in position_mutations:
                        position_mutations[aa_pos] = []
                    position_mutations[aa_pos].append((mutation_id, mutation_data, target_aa))
                
                # Now iterate through predictions and find all matching mutations (1:many logic)
                for pred in predictions:
                    pred_pos = pred['position']
                    pred_aa = pred['sequon'][0]  # First letter of sequon is the amino acid
                    
                    # Find all mutations that match this prediction's position AND amino acid
                    if pred_pos in position_mutations:
                        for mutation_id, mutation_data, target_aa in position_mutations[pred_pos]:
                            # Check if amino acid matches (critical fix!)
                            if pred_aa == target_aa:
                                pkey = f"{gene}-{mutation_id}"
                                results.append({
                                    'pkey': pkey,
                                    'Gene': gene,
                                    'pos': pred['position'],
                                    'Sequon': pred['sequon'],
                                    'potential': pred['potential'],
                                    'jury_agreement': pred['jury_agreement'],
                                    'n_glyc_result': pred['n_glyc_result']
                                })
                
                if stop_codons_skipped > 0:
                    print(f"Skipped {stop_codons_skipped} mutations with stop codons for {gene}")
        
        return results

    def collect_wt_sites(self, directories, threshold, signalp_cache=None):
        """Collect WT site predictions across provided directories."""
        wt_sites_by_gene: dict[str, list[dict]] = {}
        wt_signalp: dict[str, dict] = {}
        site_rows = []
        directories = [Path(d) for d in directories if d]
        signalp_cache = signalp_cache or {}
        for directory in directories:
            if not directory.exists():
                continue
            for file_path in iter_netnglyc_output_files(directory):
                predictions_by_seq = self.parse_netnglyc_multisequence_output(
                    str(file_path), threshold=0.0
                )
                signalp_map = parse_signalp_summary(str(file_path))
                for seq_name, predictions in predictions_by_seq.items():
                    gene_key = normalize_gene_symbol(seq_name)
                    seq_signalp = (
                        signalp_map.get(seq_name)
                        or signalp_cache.get(seq_name)
                        or {}
                    )
                    for pred in predictions:
                        jury_score = self._parse_jury_score(pred["jury_agreement"])
                        row = {
                            "pkey": f"{gene_key}-WT",
                            "Gene": gene_key,
                            "allele": "WT",
                            "seq_name": seq_name,
                            "position": pred["position"],
                            "sequon": pred["sequon"],
                            "potential": pred["potential"],
                            "jury_agreement": pred["jury_agreement"],
                            "jury_agreement_score": jury_score,
                            "n_glyc_result": pred["n_glyc_result"],
                            "n_glyc_result_code": self._encode_n_glyc_result(
                                pred["n_glyc_result"]
                            ),
                            "signalp_has_signal": self._encode_optional_bool(
                                seq_signalp.get("has_signal")
                            ),
                            "signalp_probability": seq_signalp.get("probability"),
                            "signalp_cleavage": seq_signalp.get("cleavage_site"),
                            "above_threshold": self._encode_bool(
                                pred["potential"] >= threshold
                            ),
                        }
                        site_rows.append(row)
                        wt_sites_by_gene.setdefault(gene_key, []).append(row)
                    if seq_name.lower().endswith("-wt") and seq_signalp:
                        wt_signalp[gene_key] = seq_signalp
                    elif gene_key not in wt_signalp and seq_signalp:
                        wt_signalp[gene_key] = seq_signalp
        for gene, sites in wt_sites_by_gene.items():
            if not sites and self.verbose:
                print(f"[WARN] {gene}: missing WT predictions in provided outputs")
        for gene, sites in wt_sites_by_gene.items():
            if not sites and hasattr(self, "error_logger"):
                self.error_logger.error(f"WT parser: {gene} missing predictions in provided outputs")
        return wt_sites_by_gene, wt_signalp, site_rows

    def collect_mutant_sites(self, directories, threshold, mutation_index=None, signalp_cache=None):
        """Collect mutant site predictions keyed by pkey."""
        mut_sites_by_pkey: dict[str, list[dict]] = {}
        mut_signalp: dict[str, dict] = {}
        pkey_gene_map: dict[str, str] = {}
        site_rows = []
        directories = [Path(d) for d in directories if d]
        mutation_index = mutation_index or {}
        signalp_cache = signalp_cache or {}
        for directory in directories:
            if not directory.exists():
                continue
            for file_path in iter_netnglyc_output_files(directory):
                stem = file_path.stem
                if "-wt-netnglyc" in stem:
                    continue
                predictions_by_seq = self.parse_netnglyc_multisequence_output(
                    str(file_path), threshold=0.0
                )
                signalp_map = parse_signalp_summary(str(file_path))
                for seq_name, predictions in predictions_by_seq.items():
                    gene_name, mutation_id = extract_mutation_from_sequence_name(seq_name)
                    if not mutation_id:
                        continue
                    gene_key = (gene_name or "").upper()
                    mut_id_norm = normalize_mutation_id(mutation_id)
                    if mutation_index and gene_key in mutation_index:
                        if mut_id_norm not in mutation_index[gene_key]:
                            continue
                    pkey = f"{gene_key}-{mut_id_norm}"
                    seq_signalp = signalp_map.get(seq_name, {})
                    if not seq_signalp:
                        seq_signalp = signalp_cache.get(seq_name, {})
                    pkey_gene_map[pkey] = gene_key
                    mut_signalp[pkey] = seq_signalp
                    for pred in predictions:
                        jury_score = self._parse_jury_score(pred["jury_agreement"])
                        row = {
                            "pkey": pkey,
                            "Gene": gene_key,
                            "allele": "MUT",
                            "seq_name": seq_name,
                            "position": pred["position"],
                            "sequon": pred["sequon"],
                            "potential": pred["potential"],
                            "jury_agreement": pred["jury_agreement"],
                            "jury_agreement_score": jury_score,
                            "n_glyc_result": pred["n_glyc_result"],
                            "n_glyc_result_code": self._encode_n_glyc_result(
                                pred["n_glyc_result"]
                            ),
                            "signalp_has_signal": self._encode_optional_bool(
                                seq_signalp.get("has_signal")
                            ),
                            "signalp_probability": seq_signalp.get("probability"),
                            "signalp_cleavage": seq_signalp.get("cleavage_site"),
                            "above_threshold": self._encode_bool(
                                pred["potential"] >= threshold
                            ),
                        }
                        site_rows.append(row)
                        mut_sites_by_pkey.setdefault(pkey, []).append(row)
        for pkey, sites in mut_sites_by_pkey.items():
            if not sites and self.verbose:
                print(f"[WARN] {pkey}: missing mutant predictions in provided outputs")
        for pkey, sites in mut_sites_by_pkey.items():
            if not sites and hasattr(self, "error_logger"):
                self.error_logger.error(f"MUT parser: {pkey} missing predictions in provided outputs")
        return mut_sites_by_pkey, mut_signalp, site_rows, pkey_gene_map

    def _classify_netnglyc_event(self, wt_site, mut_site, threshold, delta_threshold=0.05):
        wt_potential = wt_site["potential"] if wt_site else None
        mut_potential = mut_site["potential"] if mut_site else None
        above_wt = wt_potential is not None and wt_potential >= threshold
        above_mut = mut_potential is not None and mut_potential >= threshold
        if above_wt and above_mut:
            delta = mut_potential - wt_potential
            if delta >= delta_threshold:
                return "strengthened"
            if delta <= -delta_threshold:
                return "weakened"
            return "stable"
        if above_wt and not above_mut:
            return "lost"
        if not above_wt and above_mut:
            return "gained"
        if wt_site or mut_site:
            return "subthreshold"
        return "no_site"

    def _encode_bool(self, value):
        return 1 if value else 0

    def _encode_optional_bool(self, value):
        if value is None:
            return None
        return 1 if value else 0

    def _parse_jury_score(self, jury_string):
        if not jury_string or "/" not in jury_string:
            return None
        cleaned = jury_string.strip().strip("()")
        try:
            num, den = cleaned.split("/")
            num = float(num)
            den = float(den)
            if den == 0:
                return None
            return num / den
        except (ValueError, ZeroDivisionError):
            return None

    def _encode_n_glyc_result(self, result):
        scale = {
            "+++": 3,
            "++": 2,
            "+": 1,
            "-": -1,
            "--": -2,
            "---": -3,
        }
        return scale.get((result or "").strip(), 0)

    def _encode_classification(self, classification):
        mapping = {
            "gained": 2,
            "strengthened": 1,
            "stable": 0,
            "subthreshold": -1,
            "lost": -2,
            "weakened": -1,
            "no_site": 0,
        }
        return mapping.get(classification, 0)

    def _summarize_netnglyc_variant(
        self,
        gene,
        pkey,
        wt_sites,
        mut_sites,
        wt_sig,
        mut_sig,
        threshold,
    ):
        events = []
        wt_by_pos = {site["position"]: site for site in wt_sites}
        mut_by_pos = {site["position"]: site for site in mut_sites}
        all_positions = sorted(set(wt_by_pos.keys()) | set(mut_by_pos.keys()))
        for pos in all_positions:
            wt_site = wt_by_pos.get(pos)
            mut_site = mut_by_pos.get(pos)
            wt_potential = wt_site["potential"] if wt_site else None
            mut_potential = mut_site["potential"] if mut_site else None
            delta = None
            if wt_potential is not None or mut_potential is not None:
                delta = (mut_potential or 0.0) - (wt_potential or 0.0)
            classification = self._classify_netnglyc_event(wt_site, mut_site, threshold)
            cleavage = wt_sig.get("cleavage_site") if wt_sig else None
            if cleavage is None and mut_sig:
                cleavage = mut_sig.get("cleavage_site")
            post_cleavage = None
            if cleavage is not None:
                post_cleavage = pos >= cleavage
            events.append(
                {
                    "pkey": pkey,
                    "Gene": gene,
                    "position": pos,
                    "wt_potential": wt_potential,
                    "mut_potential": mut_potential,
                    "delta": delta,
                    "classification": classification,
                    "wt_sequon": wt_site["sequon"] if wt_site else None,
                    "mut_sequon": mut_site["sequon"] if mut_site else None,
                    "wt_above_threshold": self._encode_bool(
                        wt_potential is not None and wt_potential >= threshold
                    ),
                    "mut_above_threshold": self._encode_bool(
                        mut_potential is not None and mut_potential >= threshold
                    ),
                    "classification_code": self._encode_classification(classification),
                    "post_cleavage": self._encode_bool(post_cleavage),
                }
            )
        n_sites_wt = sum(1 for site in wt_sites if site["potential"] >= threshold)
        n_sites_mut = sum(1 for site in mut_sites if site["potential"] >= threshold)
        classification_counts = {
            "gained": 0,
            "lost": 0,
            "strengthened": 0,
            "weakened": 0,
            "stable": 0,
        }
        deltas = []
        post_cleavage_delta = 0.0
        for event in events:
            cls = event["classification"]
            if cls in classification_counts:
                classification_counts[cls] += 1
            if event["delta"] is not None:
                abs_delta = abs(event["delta"])
                deltas.append(abs_delta)
                if event["post_cleavage"]:
                    post_cleavage_delta += abs_delta
        sum_abs_delta = sum(deltas)
        max_abs_delta = max(deltas) if deltas else 0.0
        top_event = None
        if events:
            top_event = max(
                events,
                key=lambda ev: (abs(ev["delta"] or 0), ev["classification"] == "gained"),
            )
        frac_post_cleavage = (
            (post_cleavage_delta / sum_abs_delta) if sum_abs_delta > 0 else 0.0
        )
        qc_flags = []
        if not wt_sites:
            qc_flags.append("missing_wt")
        if not mut_sites:
            qc_flags.append("missing_mut")
        if sum_abs_delta == 0 and (n_sites_wt or n_sites_mut):
            qc_flags.append("no_delta")
        if wt_sig is None and mut_sig is None:
            qc_flags.append("no_signalp")
        summary = {
            "pkey": pkey,
            "Gene": gene,
            "n_sites_wt": n_sites_wt,
            "n_sites_mut": n_sites_mut,
            "count_gained": classification_counts["gained"],
            "count_lost": classification_counts["lost"],
            "count_strengthened": classification_counts["strengthened"],
            "count_weakened": classification_counts["weakened"],
            "count_stable": classification_counts["stable"],
            "max_abs_delta": max_abs_delta,
            "sum_abs_delta": sum_abs_delta,
            "top_event_type": top_event["classification"] if top_event else None,
            "top_event_delta": top_event["delta"] if top_event else None,
            "top_event_position": top_event["position"] if top_event else None,
            "top_event_classification_code": self._encode_classification(
                top_event["classification"]
            )
            if top_event
            else None,
            "wt_signalp_has_signal": self._encode_optional_bool(
                wt_sig.get("has_signal") if wt_sig else None
            ),
            "wt_signalp_probability": None
            if wt_sig is None
            else wt_sig.get("probability"),
            "wt_signalp_cleavage": None
            if wt_sig is None
            else wt_sig.get("cleavage_site"),
            "mut_signalp_has_signal": self._encode_optional_bool(
                mut_sig.get("has_signal") if mut_sig else None
            ),
            "mut_signalp_probability": None
            if mut_sig is None
            else mut_sig.get("probability"),
            "mut_signalp_cleavage": None
            if mut_sig is None
            else mut_sig.get("cleavage_site"),
            "frac_effect_post_cleavage": frac_post_cleavage,
            "qc_flags": ",".join(qc_flags) if qc_flags else "",
        }
        for event in events:
            event["wt_above_threshold"] = self._encode_bool(
                event["wt_above_threshold"]
            )
            event["mut_above_threshold"] = self._encode_bool(
                event["mut_above_threshold"]
            )
            event["post_cleavage"] = self._encode_optional_bool(event["post_cleavage"])
            event["classification_code"] = self._encode_classification(
                event["classification"]
            )
        return summary, events

    def _write_tsv(self, path, fieldnames, rows):
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def build_netnglyc_ensemble(
        self,
        wt_dirs,
        mut_dirs,
        mapping_lookup,
        threshold,
        summary_path,
        signalp_cache=None,
    ):
        mutation_index = load_mutation_index(mapping_lookup)
        signalp_cache = signalp_cache or load_signalp_cache(self.cache_dir)
        wt_sites_by_gene, wt_signalp, wt_site_rows = self.collect_wt_sites(
            wt_dirs, threshold, signalp_cache=signalp_cache
        )
        (
            mut_sites_by_pkey,
            mut_signalp,
            mut_site_rows,
            pkey_gene_map,
        ) = self.collect_mutant_sites(
            mut_dirs, threshold, mutation_index, signalp_cache=signalp_cache
        )
        sites_rows = wt_site_rows + mut_site_rows
        expected_pkeys = set(mut_sites_by_pkey.keys())
        for gene, mutation_ids in mutation_index.items():
            for mut_id in mutation_ids:
                expected_pkeys.add(f"{gene}-{mut_id}")
                pkey_gene_map.setdefault(f"{gene}-{mut_id}", gene)
        summary_rows = []
        events_rows = []
        for pkey in sorted(expected_pkeys):
            gene = pkey_gene_map.get(pkey)
            if not gene:
                gene = pkey.split("-", 1)[0]
            wt_sites = wt_sites_by_gene.get(gene, [])
            mut_sites = mut_sites_by_pkey.get(pkey, [])
            wt_sig = wt_signalp.get(gene)
            mut_sig = mut_signalp.get(pkey)
            summary, events = self._summarize_netnglyc_variant(
                gene, pkey, wt_sites, mut_sites, wt_sig or {}, mut_sig or {}, threshold
            )
            summary_rows.append(summary)
            events_rows.extend(events)
        events_fields = [
            "pkey",
            "Gene",
            "position",
            "wt_potential",
            "mut_potential",
            "delta",
            "classification",
            "classification_code",
            "wt_sequon",
            "mut_sequon",
            "wt_above_threshold",
            "mut_above_threshold",
            "post_cleavage",
        ]
        summary_fields = [
            "pkey",
            "Gene",
            "n_sites_wt",
            "n_sites_mut",
            "count_gained",
            "count_lost",
            "count_strengthened",
            "count_weakened",
            "count_stable",
            "max_abs_delta",
            "sum_abs_delta",
            "top_event_type",
            "top_event_classification_code",
            "top_event_delta",
            "top_event_position",
            "wt_signalp_has_signal",
            "wt_signalp_probability",
            "wt_signalp_cleavage",
            "mut_signalp_has_signal",
            "mut_signalp_probability",
            "mut_signalp_cleavage",
            "frac_effect_post_cleavage",
            "qc_flags",
        ]
        sites_fields = [
            "pkey",
            "Gene",
            "allele",
            "seq_name",
            "position",
            "sequon",
            "potential",
            "jury_agreement",
            "jury_agreement_score",
            "n_glyc_result",
            "n_glyc_result_code",
            "signalp_has_signal",
            "signalp_probability",
            "signalp_cleavage",
            "above_threshold",
        ]
        summary_path = Path(summary_path)
        events_path = summary_path.with_suffix(".events.tsv")
        sites_path = summary_path.with_suffix(".sites.tsv")
        self._write_tsv(str(summary_path), summary_fields, summary_rows)
        self._write_tsv(str(events_path), events_fields, events_rows)
        self._write_tsv(str(sites_path), sites_fields, sites_rows)

    def save_parsed_results(self, results, output_file):
        """Save parsed results to TSV file"""
        import csv
        
        if not results:
            print("No results to save")
            return
            
        with open(output_file, 'w', newline='') as f:
            fieldnames = ['pkey', 'Gene', 'pos', 'Sequon', 'potential', 'jury_agreement', 'n_glyc_result']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
            
        print(f"Saved {len(results)} parsed results to {output_file}")

    def process_directory(self, input_dir, output_dir, pattern=None, processing_mode="auto"):
        """
        Process all FASTA files in directory with intelligent processing mode detection
        Applies intelligent processing strategy to each file based on sequence count
        """
        os.makedirs(output_dir, exist_ok=True)

        # Find all FASTA files using flexible discovery
        if pattern:
            # Legacy support: use provided pattern
            input_files = list(Path(input_dir).glob(pattern))
        else:
            # New flexible discovery: scan for FASTA files with various extensions
            discovered_fastas = discover_fasta_files(input_dir)
            input_files = [Path(file_path) for file_path in discovered_fastas.values()]
            
        if not input_files:
            if pattern:
                print(f"No files matching {pattern} found")
            else:
                print(f"No FASTA files found in {input_dir}")
                if self.verbose:
                    print(f"  Searched for extensions: .fasta, .fa, .fas, .fna")
            return {"total": 0, "success": 0, "failed": 0}

        if self.verbose:
            print(f"Processing {len(input_files)} files from directory: {input_dir}")
            if not pattern:
                # Show discovered extensions for flexible discovery
                extensions_found = set()
                for file_path in input_files:
                    ext = file_path.suffix.lower()
                    if ext:
                        extensions_found.add(ext)
                if extensions_found:
                    print(f"  FASTA file extensions found: {', '.join(sorted(extensions_found))}")
        else:
            print(f"Processing {len(input_files)} files from directory: {input_dir}")
        if self.use_signalp and self.signalp_handler.signalp6_available:
            print("SignalP 6 will be used for signal peptide prediction")

        # Analyze each file and determine processing strategy
        file_strategies = []
        total_sequences = 0
        
        print("\n=== Analyzing files for optimal processing strategies ===")
        for fasta_file in input_files:
            seq_count = self.count_sequences_in_fasta(str(fasta_file))
            total_sequences += seq_count
            
            # Determine processing mode for this specific file
            file_processing_mode = self.detect_processing_mode(str(fasta_file), processing_mode)
            
            output_file = os.path.join(output_dir, fasta_file.stem + "-netnglyc.out")
            
            file_strategies.append({
                'file_path': str(fasta_file),
                'output_file': output_file,
                'seq_count': seq_count,
                'processing_mode': file_processing_mode,
                'file_name': fasta_file.name
            })
            
            print(f"{fasta_file.name}: {seq_count} sequences -> {file_processing_mode} mode")

        print(f"\nTotal: {len(input_files)} files, {total_sequences} sequences")
        print("=" * 60)

        # Process files with their determined strategies
        results = {"total": len(input_files), "success": 0, "failed": 0, "errors": [], "processing_summary": {}}
        start_time = time.time()

        # Track processing modes used
        mode_counts = {"single": 0, "parallel": 0, "batch": 0, "sequential": 0}
        
        # Use ProcessPoolExecutor for true parallelism in directory processing
        with ProcessPoolExecutor(max_workers=min(self.max_workers, len(input_files))) as executor:
            futures = {}
            
            # Submit jobs based on each file's optimal processing strategy
            for strategy in file_strategies:
                future = executor.submit(
                    self._process_single_file_with_strategy,
                    strategy
                )
                futures[future] = strategy
                mode_counts[strategy['processing_mode']] += 1

            # Process completed tasks
            for i, future in enumerate(as_completed(futures)):
                strategy = futures[future]
                try:
                    success, output_file, error, processing_info = future.result()
                    
                    if success:
                        results["success"] += 1
                        # Validate the output
                        if os.path.exists(output_file):
                            with open(output_file, 'r') as f:
                                content = f.read()
                            is_valid, warnings, prediction_count = self.validate_netnglyc_output(
                                content, strategy['seq_count']
                            )
                            
                            processing_info.update({
                                'predictions_found': prediction_count,
                                'validation_warnings': warnings
                            })
                            
                        
                        results["processing_summary"][strategy['file_name']] = processing_info
                        print(f"SUCCESS: {strategy['file_name']}: {processing_info.get('predictions_found', 0)} predictions")
                    else:
                        results["failed"] += 1
                        results["errors"].append({
                            "file": strategy['file_path'],
                            "error": error
                        })
                        print(f"ERROR: {strategy['file_name']}: {error}")
                        
                except Exception as e:
                    results["failed"] += 1
                    results["errors"].append({
                        "file": strategy['file_path'],
                        "error": f"Processing exception: {e}"
                    })
                    print(f"EXCEPTION: {strategy['file_name']}: Exception - {e}")

                # Progress update
                done = results["success"] + results["failed"]
                if done % 5 == 0 or done == results["total"]:
                    elapsed = time.time() - start_time
                    rate = done / elapsed if elapsed > 0 else 0
                    print(f"Progress: {done}/{results['total']} files ({rate:.1f} files/sec)")

        elapsed_total = time.time() - start_time
        
        # Print comprehensive summary
        print(f"\n{'=' * 60}")
        print(f"Directory Processing Summary:")
        print(f"   - Total files: {results['total']}")
        print(f"   - Successful: {results['success']}")
        print(f"   - Failed: {results['failed']}")
        print(f"   - Processing time: {elapsed_total / 60:.1f} minutes")
        print(f"   - Processing modes used:")
        for mode, count in mode_counts.items():
            if count > 0:
                print(f"     - {mode}: {count} files")
        
        # Show total predictions found
        total_predictions = sum(info.get('predictions_found', 0) 
                              for info in results["processing_summary"].values())
        print(f"   - Total predictions found: {total_predictions}")
        
        if results["errors"]:
            print(f"   - Errors encountered: {len(results['errors'])}")
        
        print("=" * 60)

        return results

    def _process_single_file_with_strategy(self, strategy):
        """
        Worker function for processing a single file with its determined strategy
        This runs in a separate process for true parallelism
        """
        file_path = strategy['file_path']
        output_file = strategy['output_file']
        processing_mode = strategy['processing_mode']
        seq_count = strategy['seq_count']
        
        # Create a new processor instance for this worker process
        worker_processor = RobustDockerNetNGlyc(
            use_signalp=self.use_signalp,
            max_workers=1,  # Each worker handles one file
            cache_dir=self.cache_dir,
            docker_timeout=self.docker_timeout,
            keep_intermediates=self.keep_intermediates,
            verbose=self.verbose,
            native_bin=getattr(self, 'native_path', None),
        )
        
        processing_info = {
            'processing_mode': processing_mode,
            'sequence_count': seq_count,
            'start_time': time.time()
        }
        
        try:
            # Apply the determined processing strategy
            if processing_mode == "single":
                success, output, error = worker_processor.process_single_fasta(file_path, output_file, 0)
            elif processing_mode == "parallel":
                # For directory processing,  a smaller worker count per file is used
                parallel_workers = min(4, max(1, seq_count // 50))  # Adaptive worker count
                success, output, error = worker_processor.process_parallel_docker(
                    file_path, output_file, parallel_workers
                )
                processing_info['parallel_workers_used'] = parallel_workers
            elif processing_mode == "batch":
                # Use adaptive batch size based on sequence count
                batch_size = min(100, max(20, seq_count // 10))
                success, output, error = worker_processor.process_fasta_batched(
                    file_path, output_file, batch_size, 0
                )
                processing_info['batch_size_used'] = batch_size
            elif processing_mode == "sequential":
                # Sequential is same as single for now
                success, output, error = worker_processor.process_single_fasta(file_path, output_file, 0)
            else:
                return False, output_file, f"Unknown processing mode: {processing_mode}", processing_info
            
            processing_info['end_time'] = time.time()
            processing_info['duration'] = processing_info['end_time'] - processing_info['start_time']
            processing_info['success'] = success
            
            return success, output_file, error, processing_info
            
        except Exception as e:
            processing_info['end_time'] = time.time()
            processing_info['duration'] = processing_info['end_time'] - processing_info['start_time']
            processing_info['success'] = False
            return False, output_file, f"Worker exception: {e}", processing_info
        finally:
            # Clean up worker processor
            if hasattr(worker_processor, 'temp_dir') and os.path.exists(worker_processor.temp_dir):
                shutil.rmtree(worker_processor.temp_dir, ignore_errors=True)

    def count_sequences_in_fasta(self, fasta_file):
        """Count the number of sequences in a FASTA file efficiently"""
        try:
            count = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count
        except Exception as e:
            print(f"Error counting sequences in {fasta_file}: {e}")
            return 0

    def detect_processing_mode(self, fasta_file, user_mode="auto"):
        """
        Intelligently detect the best processing mode based on sequence count
        - 1-50 sequences: single
        - 51-500 sequences: parallel
        - 501+ sequences: batch (split into batches)
        """
        if user_mode != "auto":
            print(f"Using user-specified processing mode: {user_mode}")
            return user_mode
        
        seq_count = self.count_sequences_in_fasta(fasta_file)
        
        if seq_count == 0:
            print("Warning: No sequences found in FASTA file")
            return "single"
        elif seq_count <= 50:
            mode = "single"
            print(f"Auto-detected mode: {mode} ({seq_count} sequences)")
        elif seq_count <= 500:
            mode = "parallel"
            print(f"Auto-detected mode: {mode} ({seq_count} sequences)")
        else:
            mode = "batch"
            print(f"Auto-detected mode: {mode} ({seq_count} sequences - requires batching to avoid NetNGlyc limits)")
        
        return mode

    def validate_netnglyc_output(self, output_content, expected_sequences=None):
        """
        Validate NetNGlyc output for quality and completeness
        Returns: (is_valid, warning_messages, prediction_count)
        """
        warnings = []
        prediction_count = 0
        
        if not output_content or len(output_content.strip()) < 100:
            return False, ["Output too short or empty"], 0
        
        # Check for successful execution - support multiple NetNGlyc output formats
        valid_headers = [
            "# Predictions for N-Glycosylation sites",  # Standard single-file format
            ">"  # Batch/combined file format (starts with >GENE_aa-netnglyc)
        ]
        
        has_valid_header = any(header in output_content for header in valid_headers)
        if not has_valid_header:
            return False, ["Missing NetNGlyc header - execution may have failed"], 0
            
        # Additional content validation - ensure it contains NetNGlyc-specific content
        netnglyc_indicators = [
            "Name:",  # Sequence entries
            "(Threshold=0.5)",  # Threshold indicators  
            "N-Glyc result",  # Result headers
            "netNglyc: SignalP"  # SignalP integration messages
        ]
        
        content_score = sum(1 for indicator in netnglyc_indicators if indicator in output_content)
        if content_score < 2:  # Require at least 2 indicators present
            return False, ["Output does not contain expected NetNGlyc content"], 0
        
        # Count actual sequences processed - handle both single and batch formats
        sequence_names_standard = output_content.count("Name:")  # Standard format
        sequence_names_batch = len([line for line in output_content.split('\n') if line.startswith('>')]) # Batch format
        sequence_names = max(sequence_names_standard, sequence_names_batch)
        
        if expected_sequences and sequence_names != expected_sequences:
            warnings.append(f"Expected {expected_sequences} sequences, found {sequence_names}")
            if self.verbose:
                warnings.append(f"Standard format count: {sequence_names_standard}, Batch format count: {sequence_names_batch}")
        
        # Count predictions more accurately
        lines = output_content.split('\n')
        in_results = False
        separators_seen = 0
        
        for line in lines:
            line = line.strip()
            
            # Find results table start
            if "SeqName" in line and "Position" in line:
                in_results = True
                continue
            
            # Handle separators in results section
            if in_results and line.startswith('---'):
                separators_seen += 1
                if separators_seen == 1:
                    continue  # Start reading predictions after first separator
                elif separators_seen == 2:
                    break  # End of results section
            
            # Count valid prediction lines
            if in_results and separators_seen == 1 and line and not line.startswith('#'):
                # Skip header continuation lines
                if 'agreement' in line and 'result' in line and len(line.split()) <= 2:
                    continue
                
                # Check for valid prediction line format
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        # Validate it's a real prediction (position should be a number)
                        int(parts[1])
                        float(parts[3])  # Potential should be a float
                        prediction_count += 1
                    except (ValueError, IndexError):
                        continue
        
        # Check for common issues
        if "No sites predicted" in output_content and prediction_count == 0:
            warnings.append("No glycosylation sites predicted - this may be biologically accurate")
        
        if "error" in output_content.lower() or "failed" in output_content.lower():
            warnings.append("Potential processing errors detected in output")
        
        # Final validation - file is valid if it has proper header and content
        is_valid = has_valid_header and content_score >= 2
        return is_valid, warnings, prediction_count

    def process_parallel_docker(self, fasta_file, output_file, max_workers=4):
        """
        Process FASTA file by splitting into individual sequences and running
        them in parallel.
        """
        
        try:
            # Read all sequences
            sequences = read_fasta(fasta_file)
            total_sequences = len(sequences)
            
            if total_sequences == 0:
                return False, output_file, "No sequences found in FASTA file"
            
            #print(f"Processing {total_sequences} sequences with {max_workers} parallel Docker containers")
            
            # Create temporary files for individual sequences (save to outputs directory for troubleshooting)
            outputs_dir = os.path.dirname(output_file)
            temp_files = []
            output_files = []
            
            for i, (seq_name, sequence) in enumerate(sequences.items()):
                # Create individual FASTA file (save to outputs directory for troubleshooting)
                temp_fasta = os.path.join(outputs_dir, f"seq_{i}_{seq_name.replace('/', '_')}.fasta")
                temp_output = os.path.join(outputs_dir, f"seq_{i}_output.out")
                
                with open(temp_fasta, 'w') as f:
                    f.write(f">{seq_name}\n")
                    # Write sequence in lines of 80 characters
                    for j in range(0, len(sequence), 80):
                        f.write(sequence[j:j+80] + "\n")
                
                temp_files.append((temp_fasta, temp_output, seq_name))
                output_files.append(temp_output)
            
            # Process sequences in parallel using ProcessPoolExecutor
            #print(f"Starting parallel processing with {max_workers} workers...")
            successful_outputs = []
            failed_count = 0
            
            # Prepare arguments for worker function (add instance parameters)
            worker_args = []
            for temp_fasta, temp_output, seq_name in temp_files:
                args = (temp_fasta, temp_output, seq_name,
                        self.use_signalp, self.cache_dir, self.docker_timeout, self.verbose,
                        getattr(self, 'native_path', None))
                worker_args.append((args, seq_name))
            
            # Use ProcessPoolExecutor for true parallelism
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(_process_single_sequence_worker, args): seq_name 
                          for args, seq_name in worker_args}
                
                for future in as_completed(futures):
                    seq_name = futures[future]
                    try:
                        success, temp_output, error, _ = future.result()
                        if success and os.path.exists(temp_output):
                            successful_outputs.append(temp_output)
                            print(f"[OK] Completed: {seq_name}")
                        else:
                            failed_count += 1
                            print(f"[X] Failed: {seq_name} - {error}")
                    except Exception as e:
                        failed_count += 1
                        print(f"[X] Exception processing {seq_name}: {e}")
            
            print(f"Parallel processing completed: {len(successful_outputs)} successful, {failed_count} failed")
            
            # Combine all successful outputs
            if successful_outputs:
                combined_success = self.combine_parallel_outputs(successful_outputs, output_file, total_sequences)
                if combined_success:
                    #print(f"Successfully combined outputs into {output_file}")
                    return True, output_file, None
                else:
                    return False, output_file, "Failed to combine parallel outputs"
            else:
                return False, output_file, "All parallel processes failed"
                
        except Exception as e:
            return False, output_file, f"Parallel processing error: {e}"
        finally:
            # Clean up temporary files (only if not keeping intermediates)
            if not self.keep_intermediates:
                for temp_fasta, temp_output, _ in temp_files:
                    for temp_file in [temp_fasta, temp_output]:
                        if os.path.exists(temp_file):
                            try:
                                os.remove(temp_file)
                            except:
                                pass

    def combine_parallel_outputs(self, output_files, final_output_file, total_sequences):
        """Combine multiple parallel NetNGlyc outputs into a single file"""
        try:
            if not output_files:
                return False
            
            # Start building combined output
            combined_lines = []
            
            # Create header section (SignalP warning removed as it's misleading - SignalP 6.0 is working)
            combined_lines.extend([
                f"# Predictions for N-Glycosylation sites in {total_sequences} sequences",
                ""
            ])
            
            # Collect all sequence information and results
            all_predictions = []
            sequence_info = []
            
            for output_file in output_files:
                try:
                    with open(output_file, 'r') as f:
                        content = f.read()
                    
                    # Extract sequence information (Name: lines and sequence display)
                    lines = content.split('\n')
                    in_sequence_display = False
                    current_seq_lines = []
                    
                    for line in lines:
                        if line.startswith('Name:'):
                            if current_seq_lines:
                                sequence_info.extend(current_seq_lines)
                                current_seq_lines = []
                            sequence_info.append(line)
                            in_sequence_display = True
                        elif in_sequence_display and (line and not line.startswith('#') and 'SeqName' not in line):
                            if line.strip().startswith('---'):
                                in_sequence_display = False
                                current_seq_lines.append("")  # Add spacing
                            else:
                                current_seq_lines.append(line)
                    
                    if current_seq_lines:
                        sequence_info.extend(current_seq_lines)
                    
                    # Extract prediction results  
                    in_results = False
                    found_first_separator = False
                    
                    for line in lines:
                        if "SeqName" in line and "Position" in line:
                            in_results = True
                            continue
                        elif line.strip().startswith('---') and in_results:
                            if not found_first_separator:
                                found_first_separator = True
                                continue
                            else:
                                break  # End of results
                        elif in_results and found_first_separator and line.strip() and not line.startswith('#'):
                            all_predictions.append(line.strip())
                
                except Exception as e:
                    print(f"Error reading parallel output {output_file}: {e}")
                    continue
            
            # Add sequence information
            combined_lines.extend(sequence_info)
            
            # Add results table
            if all_predictions:
                combined_lines.extend([
                    "(Threshold=0.5)",
                    "----------------------------------------------------------------------",
                    "SeqName        Position      Potential          Jury          N-Glyc",
                    "\t\t\t\t              agreement       result", 
                    "----------------------------------------------------------------------"
                ])
                combined_lines.extend(all_predictions)
                combined_lines.append("----------------------------------------------------------------------")
            else:
                combined_lines.extend([
                    "(Threshold=0.5)",
                    "----------------------------------------------------------------------", 
                    "No sites predicted",
                    "----------------------------------------------------------------------"
                ])
            
            combined_lines.append("")
            
            # Write combined output
            with open(final_output_file, 'w') as f:
                f.write('\n'.join(combined_lines))
            
            print(f"Combined {len(output_files)} parallel outputs into {final_output_file}")
            return True
            
        except Exception as e:
            print(f"Error combining parallel outputs: {e}")
            return False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)


def clear_all_caches(cache_dir=None):
    """Clear all caches (both SignalP and NetNGlyc)"""
    cleared = []

    # Clear SignalP cache
    sp_cache = os.path.join(os.path.expanduser("~"), ".signalp6_cache")
    if cache_dir:
        sp_cache = os.path.join(cache_dir, ".signalp6_cache")

    if os.path.exists(sp_cache):
        shutil.rmtree(sp_cache)
        cleared.append(f"SignalP cache: {sp_cache}")

    # Clear NetNGlyc cache
    ng_cache = os.path.join(os.path.expanduser("~"), ".netnglyc_cache")
    if cache_dir:
        ng_cache = os.path.join(cache_dir, ".netnglyc_cache")

    if os.path.exists(ng_cache):
        shutil.rmtree(ng_cache)
        cleared.append(f"NetNGlyc cache: {ng_cache}")

    return cleared


def run_full_pipeline_mode(args, failure_map, parser):
    """Execute the full NetNGlyc pipeline: synthesize FASTAs, run NetNGlyc, parse outputs."""
    if not args.input or not args.output:
        parser.error("For full-pipeline mode: both input (FASTA file/directory) and output (TSV file) are required")

    if not args.mapping_dir:
        parser.error("For full-pipeline mode: --mapping-dir is REQUIRED (directory containing mutation mapping CSV files)")

    temp_output_dir = tempfile.mkdtemp(prefix="netnglyc_outputs_")
    temp_sequence_dir = tempfile.mkdtemp(prefix="netnglyc_sequences_")
    wt_temp_holder = None

    try:
        wt_sequences, wt_temp_holder = load_wt_sequence_map(args.input, wt_header="ORF")
        if not wt_sequences:
            parser.error(f"No WT ORF sequences were found in {args.input}")

        mapping_lookup = {
            gene.upper(): path for gene, path in discover_mapping_files(args.mapping_dir).items()
        }
        if not mapping_lookup:
            parser.error(f"No mapping CSV files found in {args.mapping_dir}")

        wt_fasta_dir, mut_fasta_dir, build_summary = synthesize_gene_fastas(
            wt_sequences,
            mapping_lookup,
            temp_sequence_dir,
            log_path=args.log,
            failure_map=failure_map,
        )

        wt_fastas = sorted(Path(wt_fasta_dir).glob("*.fasta"))
        if not wt_fastas:
            parser.error("Failed to synthesize any WT FASTA files for NetNGlyc")
        mut_fastas = sorted(Path(mut_fasta_dir).glob("*.fasta"))

        print("\n=== FASTA synthesis summary ===")
        for entry in build_summary:
            mut_msg = f"{entry['mutant_count']} mutants" if entry['mutant_count'] else "no mutants"
            print(f"  {entry['gene']}: WT -> {entry['wt_path']} | Mutants -> {mut_msg}")
        print("=" * 60)

        wt_outputs_dir = Path(temp_output_dir) / "wt"
        mut_outputs_dir = Path(temp_output_dir) / "mut"
        wt_outputs_dir.mkdir(parents=True, exist_ok=True)
        mut_outputs_dir.mkdir(parents=True, exist_ok=True)

        with RobustDockerNetNGlyc(
            use_signalp=True,
            max_workers=args.workers,
            cache_dir=args.cache_dir,
            docker_timeout=args.batch_timeout,
            keep_intermediates=args.keep_intermediates,
            verbose=args.verbose,
            native_bin=args.native_netnglyc_bin,
        ) as processor:

            print("\nRunning NetNGlyc on WT amino acid FASTAs...")
            wt_results_summary = processor.process_directory(
                str(wt_fasta_dir),
                str(wt_outputs_dir),
            )
            print(f"WT NetNGlyc runs completed: {wt_results_summary['success']}/{wt_results_summary['total']}")

            if mut_fastas:
                print("\nRunning NetNGlyc on mutant amino acid FASTAs...")
                mut_results_summary = processor.process_directory(
                    str(mut_fasta_dir),
                    str(mut_outputs_dir),
                )
                print(f"Mutant NetNGlyc runs completed: {mut_results_summary['success']}/{mut_results_summary['total']}")
            else:
                print("No mutant FASTA files were generated; skipping mutant NetNGlyc execution")

            print("\nBuilding NetNGlyc ensemble outputs...")
            signalp_cache = load_signalp_cache(args.cache_dir)
            processor.build_netnglyc_ensemble(
                wt_dirs=[str(wt_outputs_dir)],
                mut_dirs=[str(mut_outputs_dir)],
                mapping_lookup=mapping_lookup,
                threshold=args.threshold,
                summary_path=args.output,
                signalp_cache=signalp_cache,
            )
            print(f"NetNGlyc ensemble outputs written to {args.output}")

    finally:
        if wt_temp_holder:
            wt_temp_holder.cleanup()
        if args.keep_intermediates:
            print(f"WT/Mutant FASTAs saved in {temp_sequence_dir}")
            print(f"Raw NetNGlyc outputs saved in {temp_output_dir}")
        else:
            shutil.rmtree(temp_sequence_dir, ignore_errors=True)
            shutil.rmtree(temp_output_dir, ignore_errors=True)

    return 0


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Licensed NetNGlyc pipeline with SignalP 6 integration and intelligent parallel processing. Requires licensed NetNGlyc 1.0 native binary and SignalP 6.0 installation.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", nargs='?', help="Input FASTA file or directory (required for process/full-pipeline modes)")
    parser.add_argument("output", nargs='?', help="Output file or directory (required for all modes except --test/--clear-cache)")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (used with --processing-mode parallel, default: 4)")
    parser.add_argument("--cache-dir", help="Custom cache directory for SignalP/NetNGlyc results")
    parser.add_argument("--test", action="store_true",
                        help="Run test with ABCB1 sequence (no other args required)")
    parser.add_argument("--clear-cache", action="store_true",
                        help="Clear all cached results and exit (no other args required)")
    parser.add_argument("--mode", choices=["process", "parse", "full-pipeline"], default="process",
                        help="Processing mode: 'process' (run NetNGlyc only), 'parse' (parse existing outputs), 'full-pipeline' (process + parse)")
    parser.add_argument("--mapping-dir",
                        help="Directory containing mutation mapping CSV files (REQUIRED for parsing modes)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Minimum glycosylation potential threshold for predictions (default: 0.5)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="Number of sequences per batch for batch processing mode (default: 100)")
    parser.add_argument("--batch-timeout", type=int, default=5000,
                        help="Timeout in seconds for NetNGlyc execution (default: 5000s)")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep intermediate output files for debugging (don't clean up temporary directories)")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose output showing detailed processing information")
    parser.add_argument("--log",
                        help="Validation log file or directory to skip failed mutations (mutant modes only)")

    # Native execution options
    parser.add_argument("--native-netnglyc-bin",
                        help="Path to native NetNGlyc binary")

    args = parser.parse_args()
    failure_map = load_validation_failures(args.log) if args.log else {}

    # Handle cache clearing
    if args.clear_cache:
        print("Clearing all caches...")
        cleared = clear_all_caches(args.cache_dir)
        if cleared:
            for cache in cleared:
                print(f"  Cleared: {cache}")
            print("All caches cleared successfully")
        else:
            print("No caches found to clear")
        return 0

    if args.test:
        print("Running test mode...")
        test_root = tempfile.mkdtemp(prefix="netnglyc_test_")
        fasta_dir = Path(test_root) / "wt"
        mapping_dir = Path(test_root) / "mapping"
        fasta_dir.mkdir(parents=True, exist_ok=True)
        mapping_dir.mkdir(parents=True, exist_ok=True)

        test_fasta_path = fasta_dir / "ABCB1_nt.fasta"
        test_mapping_path = mapping_dir / "ABCB1_mapping.csv"
        test_output_tsv = Path("test_abcb1_results.tsv")
        if test_output_tsv.exists():
            test_output_tsv.unlink()

        test_nt_sequence = (
            "ATGGATCTGGAAGGTGATCGTAATGGTGGTGCTAAAAAAAAAAATTTTTTTAAACTGAAT"
        )

        with open(test_fasta_path, 'w') as handle:
            handle.write(">orf\n")
            handle.write(test_nt_sequence + "\n")

        with open(test_mapping_path, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(["mutant", "aamutant"])
            writer.writerow(["G10A", ""])
            writer.writerow(["A22G", "N8D"])

        test_args = argparse.Namespace(**vars(args))
        test_args.input = str(fasta_dir)
        test_args.output = str(test_output_tsv)
        test_args.mapping_dir = str(mapping_dir)
        test_args.mode = "full-pipeline"
        test_args.test = False

        try:
            run_full_pipeline_mode(test_args, {}, parser)
            print(f"\nTest pipeline completed successfully. Parsed TSV: {test_output_tsv}")
            return 0
        finally:
            shutil.rmtree(test_root, ignore_errors=True)

    # Handle parsing mode
    if args.mode == "parse":
        if not args.input or not args.output:
            parser.error("For parse mode: both input (NetNGlyc output directory) and output (TSV file) are required")
        
        # Validate parsing requirements
        if not args.mapping_dir:
            parser.error("For parse mode: --mapping-dir is REQUIRED (directory containing mutation mapping CSV files)")
        
        mapping_lookup = discover_mapping_files(args.mapping_dir)
        processor = RobustDockerNetNGlyc(
            use_signalp=False,
            cache_dir=args.cache_dir,
            docker_timeout=args.batch_timeout,
            verbose=args.verbose,
            native_bin=args.native_netnglyc_bin,
        )
        input_path = Path(args.input)
        wt_dirs = []
        mut_dirs = []
        if (input_path / "wt").exists():
            wt_dirs.append(str(input_path / "wt"))
        wt_dirs.append(str(input_path))
        if (input_path / "mut").exists():
            mut_dirs.append(str(input_path / "mut"))
        mut_dirs.append(str(input_path))

        signalp_cache = load_signalp_cache(args.cache_dir)
        processor.build_netnglyc_ensemble(
            wt_dirs=wt_dirs,
            mut_dirs=mut_dirs,
            mapping_lookup=mapping_lookup,
            threshold=args.threshold,
            summary_path=args.output,
            signalp_cache=signalp_cache,
        )
        print(f"Wrote NetNGlyc ensemble summary to {args.output}")
        return 0
        
    elif args.mode == "full-pipeline":
        return run_full_pipeline_mode(args, failure_map, parser)

    # Normal processing mode - validate required arguments  
    if not args.input or not args.output:
        parser.error("For process mode: both input (FASTA file/directory) and output (NetNGlyc output file/directory) are required")

    with RobustDockerNetNGlyc(
            use_signalp=True,
            max_workers=args.workers,
            cache_dir=args.cache_dir,
            docker_timeout=args.batch_timeout,
            verbose=args.verbose,
            native_bin=args.native_netnglyc_bin,
    ) as processor:

        if os.path.isfile(args.input):
            # Single file
            success, output, error = processor.process_single_fasta(
                args.input, args.output, 0
            )
            if success:
                print(f"Successfully processed {args.input}")
            else:
                print(f"Failed: {error}")
                sys.exit(1)
        else:
            # Directory processing with intelligent mode detection
            results = processor.process_directory(args.input, args.output)
            print(f"\n=== Processing Summary ===")
            print(f"Total files: {results['total']}")
            print(f"Successful: {results['success']}")
            print(f"Failed: {results['failed']}")

            if results['errors']:
                print("\nErrors:")
                for err in results['errors'][:5]:
                    print(f"  {err['file']}: {err['error']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
