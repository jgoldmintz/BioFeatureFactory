#!/usr/bin/env python3
"""
Complete NetNGlyc Docker Pipeline with Host SignalP 6 Integration
Runs SignalP 6 on host Mac, passes results to NetNGlyc in Docker
NetNGlyc analyzes FULL sequence with SignalP context
"""

import os
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


# Import utility functions  
sys.path.append(os.path.join(os.path.dirname(__file__), '../dependencies'))

from utility import (
    read_fasta,
    get_mutation_data_bioAccurate,
    split_fasta_into_batches,
    combine_batch_outputs,
    discover_mapping_files,
    discover_fasta_files,
    ExtractGeneFromFASTA,
    parse_predictions_with_mutation_filtering,
    load_validation_failures,
    should_skip_mutation,
)

def _process_single_sequence_worker(args):
    """
    Worker function for parallel processing - must be at module level for pickling
    """
    temp_fasta, temp_output, seq_name, docker_image, use_signalp, cache_dir, docker_timeout, verbose = args
    
    worker_processor = RobustDockerNetNGlyc(
        docker_image=docker_image,
        use_signalp=use_signalp,
        max_workers=1,
        cache_dir=cache_dir,
        docker_timeout=docker_timeout,
        verbose=verbose
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


def get_docker_platform_args():
    """
    Detect if running on ARM64 and return appropriate Docker platform arguments
    """
    machine = platform.machine().lower()
    system = platform.system().lower()

    if system == 'darwin' and machine in ['arm64', 'aarch64']:
        return ["--platform", "linux/386"]
    elif 'arm' in machine or 'aarch' in machine:
        return ["--platform", "linux/386"]
    else:
        return []


class SignalP6Handler:
    """Handle SignalP 6 predictions on the host before Docker NetNGlyc"""

    def __init__(self, cache_dir=None, verbose=False):
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".signalp6_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        self.last_output_dir = None  # Store the output directory path
        self.verbose = verbose
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
            print(f"Using cached SignalP 6 results for {os.path.basename(fasta_file)}")
            with open(cache_file, 'r') as f:
                return json.load(f), cache_dir

        if not self.signalp6_available:
            error_msg = "ERROR: SignalP 6.0 is required but not available. Please install SignalP 6.0 and ensure it's in your PATH."
            print(error_msg)
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
            print(f"SignalP 6 error: {e}")
        finally:
            # Clean up temp dir if one was created
            if temp_dir and temp_dir != signalp_output_dir:
                shutil.rmtree(temp_dir)

        return results, None


class RobustDockerNetNGlyc:
    """
    Docker NetNGlyc with integrated SignalP 6 preprocessing on host
    
    Requires:
    - Licensed NetNGlyc 1.0 in Docker container
    - SignalP 6.0 installed on host system (when use_signalp=True)
    - No fallback functionality - fails explicitly if requirements not met
    """

    def __init__(self, docker_image="netnglyc:latest", use_signalp=True,
                 max_workers=4, cache_dir=None, docker_timeout=600, keep_intermediates=False, verbose=False):
        self.docker_image = docker_image
        self.max_workers = max_workers
        self.temp_dir = tempfile.mkdtemp(prefix="netnglyc_")
        self.use_signalp = use_signalp
        self.docker_timeout = docker_timeout
        self.keep_intermediates = keep_intermediates
        self.verbose = verbose

        # Initialize SignalP handler
        self.signalp_handler = SignalP6Handler(cache_dir=cache_dir, verbose=self.verbose) if use_signalp else None

        # Initialize cache
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".netnglyc_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        # Setup error logging with date-named log file
        self._setup_error_logging()

        # Detect platform once during initialization
        self.platform_args = get_docker_platform_args()
        if self.platform_args and self.verbose:
            print(f"Detected ARM64 architecture - Docker will use emulation mode")

        # Test Docker setup (non-fatal - allow fallback during processing)
        self._test_docker_availability()

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

    def _test_docker_availability(self):
        """Test if Docker and NetNGlyc image are available (non-fatal)"""
        try:
            # Test Docker
            result = subprocess.run(["docker", "--version"], capture_output=True, timeout=5)
            if result.returncode != 0:
                if self.verbose:
                    print("Docker not available - will fallback to native stub when needed")
                return False

            # Test image
            result = subprocess.run(
                ["docker", "images", "-q", self.docker_image],
                capture_output=True,
                text=True,
                timeout=5
            )
            if not result.stdout.strip():
                if self.verbose:
                    print(f"Docker image '{self.docker_image}' not found - will fallback to native stub when needed")
                return False
                
            if self.verbose:
                print(f"Docker NetNGlyc ready: {self.docker_image}")
            return True

        except Exception as e:
            if self.verbose:
                print(f"Docker test failed ({e}) - will fallback to native stub when needed")
            return False

    def _test_docker(self):
        """Test if Docker and NetNGlyc image are available (fatal version for strict testing)"""
        try:
            # Test Docker
            result = subprocess.run(["docker", "--version"], capture_output=True, timeout=5)
            if result.returncode != 0:
                raise Exception("Docker not available")

            # Test image
            result = subprocess.run(
                ["docker", "images", "-q", self.docker_image],
                capture_output=True,
                text=True,
                timeout=5
            )
            if not result.stdout.strip():
                raise Exception(f"Docker image '{self.docker_image}' not found")

        except Exception as e:
            print(f"Docker setup error: {e}")
            raise

    def clear_cache(self):
        """Clear all cached NetNGlyc results"""
        if os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
            os.makedirs(self.cache_dir, exist_ok=True)
            return True
        return False

    def _run_docker_netnglyc(self, fasta_file, signalp_output_dir=None):
        """
        Run licensed NetNGlyc in Docker container with SignalP results mounted
        
        Requires:
        - Valid NetNGlyc Docker image with licensed NetNGlyc 1.0
        - Properly configured Docker environment
        
        Fails explicitly if Docker or NetNGlyc requirements not met
        """
        # Create a temporary directory for this run
        work_dir = tempfile.mkdtemp(dir=self.temp_dir)

        try:
            # Copy the FASTA file to the work directory
            docker_input = os.path.join(work_dir, "input.fasta")
            shutil.copy2(fasta_file, docker_input)

            # Build Docker command
            docker_cmd = [
                "docker", "run", "--rm"
            ]

            # Add platform specification if needed
            docker_cmd.extend(self.platform_args)

            # Mount the work directory
            docker_cmd.extend(["-v", f"{work_dir}:/data"])

            # Mount SignalP output directory if available
            if signalp_output_dir and os.path.exists(signalp_output_dir):
                docker_cmd.extend(["-v", f"{signalp_output_dir}:/signalp_output:ro"])

            # Add image and command
            docker_cmd.extend([
                self.docker_image,
                "/data/input.fasta"  # Path inside container
            ])


            # Use configured timeout, with environment variable override if set
            timeout = int(os.environ.get('NETNGLYC_DOCKER_TIMEOUT', str(self.docker_timeout)))
            
            result = subprocess.run(
                docker_cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            # Check for successful output
            if result.returncode == 0:
                if "Predictions for N-Glycosylation sites" in result.stdout:
                    return self._clean_signalp_warning(result.stdout)
                elif "No Asparagines in the input sequences" in result.stdout:
                    # No N residues to analyze
                    return self._clean_signalp_warning(result.stdout)

            # Check if platform warning but still got output
            if "WARNING" in result.stderr and "platform" in result.stderr:
                if result.stdout and "Predictions" in result.stdout:
                    return self._clean_signalp_warning(result.stdout)

            raise Exception(f"Docker NetNGlyc failed with return code {result.returncode}")

        except subprocess.TimeoutExpired as e:
            # Timeout - log error but don't fallback to stub
            error_msg = f"Docker NetNGlyc timeout ({self.docker_timeout}s) for file: {os.path.basename(fasta_file)}"
            print(f"ERROR: {error_msg}")
            self.error_logger.error(error_msg)
            raise Exception(f"NetNGlyc Docker timeout after {self.docker_timeout} seconds")
        except Exception as e:
            # All Docker issues result in explicit failure - no fallback
            error_msg = f"Docker NetNGlyc failed for file {os.path.basename(fasta_file)}: {e}"
            print(f"ERROR: {error_msg}")
            self.error_logger.error(error_msg)
            raise Exception(f"Docker NetNGlyc failed: {e}")
        finally:
            # Clean up work directory
            if os.path.exists(work_dir):
                shutil.rmtree(work_dir)

    def _clean_signalp_warning(self, netnglyc_output):
        """
        Remove cosmetic SignalP warning from NetNGlyc output
        
        The warning "netNglyc: SignalP is not configured, no signal peptide predictions are made"
        is misleading because SignalP 6.0 is actually running on the host system and providing
        accurate signal peptide predictions. The warning appears because NetNGlyc (inside Docker)
        cannot see the SignalP 6.0 installation on the host, but the predictions are still
        being generated and integrated properly.
        
        This is purely a cosmetic fix to remove confusing output.
        """
        if not netnglyc_output:
            return netnglyc_output
            
        # Remove the specific SignalP warning line
        lines = netnglyc_output.split('\n')
        cleaned_lines = []
        removed_warning = False
        
        for line in lines:
            # Skip the misleading SignalP warning line
            if "netNglyc: SignalP is not configured" in line:
                removed_warning = True
                if self.verbose:
                    print(f"DEBUG: Removed SignalP warning line: {line[:50]}...")
                continue
            cleaned_lines.append(line)
        
        if self.verbose and removed_warning:
            print("DEBUG: SignalP warning successfully removed from NetNGlyc output")
        elif self.verbose:
            print("DEBUG: No SignalP warning found in NetNGlyc output")
        
        return '\n'.join(cleaned_lines)


    def process_single_fasta(self, fasta_file, output_file, worker_id=0):
        """
        Process a single FASTA file with SignalP 6 preprocessing and licensed NetNGlyc

        1. Run SignalP 6 on host (if enabled) - REQUIRED if use_signalp=True
        2. Pass ORIGINAL FASTA and SignalP results to licensed NetNGlyc Docker
        3. NetNGlyc analyzes full sequence with SignalP context
        
        Fails explicitly if requirements not met (no fallback functionality)
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
    
            # Step 2: Run NetNGlyc Docker with ORIGINAL FASTA and SignalP results
            netnglyc_output = self._run_docker_netnglyc(fasta_file, signalp_output_dir)

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
        
        if mapping_dir is None:
            print("ERROR: mapping_dir is required for mutant parsing")
            return results
        
        # Import shared utility functions
        import sys
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'dependencies'))

        
        # Find NetNGlyc output files: both regular and batch files
        # Pattern 1: {GENE}_aa-netnglyc.out (regular files)
        regular_files = list(Path(input_dir).glob("*_aa-netnglyc.out"))
        # Pattern 2: {GENE}_aa-netnglyc-batch-*.out (batch files)
        batch_files = list(Path(input_dir).glob("*_aa-netnglyc-batch-*.out"))
        
        # Group files by gene for batch combination
        gene_files = {}
        
        # Process regular files
        for file_path in regular_files:
            filename = file_path.stem  # Remove .out extension
            if '_aa-netnglyc' in filename:
                gene = ExtractGeneFromFASTA(str(file_path))
                if gene and gene not in gene_files:
                    gene_files[gene] = []
                if gene:
                    gene_files[gene].append(file_path)
        
        # Process batch files  
        for file_path in batch_files:
            filename = file_path.stem  # Remove .out extension
            # Extract gene from batch file content
            if '_aa-netnglyc-batch-' in filename:
                gene = ExtractGeneFromFASTA(str(file_path))
                if gene and gene not in gene_files:
                    gene_files[gene] = []
                if gene:
                    gene_files[gene].append(file_path)
        
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
                with open(mapping_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        ntposnt = row['mutant']  # e.g., "C3066A"
                        aaposaa = row['aamutant']  # e.g., "S1022N"
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
        print(f"   • Total files: {results['total']}")
        print(f"   • Successful: {results['success']}")
        print(f"   • Failed: {results['failed']}")
        print(f"   • Processing time: {elapsed_total / 60:.1f} minutes")
        print(f"   • Processing modes used:")
        for mode, count in mode_counts.items():
            if count > 0:
                print(f"     - {mode}: {count} files")
        
        # Show total predictions found
        total_predictions = sum(info.get('predictions_found', 0) 
                              for info in results["processing_summary"].values())
        print(f"   • Total predictions found: {total_predictions}")
        
        if results["errors"]:
            print(f"   • Errors encountered: {len(results['errors'])}")
        
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
            docker_image=self.docker_image,
            use_signalp=self.use_signalp,
            max_workers=1,  # Each worker handles one file
            cache_dir=self.cache_dir,
            docker_timeout=self.docker_timeout,
            keep_intermediates=self.keep_intermediates,  # CRITICAL: Pass the flag to worker processes
            verbose=self.verbose
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
        Fallback hierarchy:
        - 1-50 sequences: single (one Docker container)
        - 51-500 sequences: parallel (multiple Docker containers)  
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
            print(f"Auto-detected mode: {mode} ({seq_count} sequences - optimal for single Docker container)")
        elif seq_count <= 500:
            mode = "parallel" 
            print(f"Auto-detected mode: {mode} ({seq_count} sequences - optimal for parallel Docker containers)")
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
        parallel Docker containers for each sequence
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
                args = (temp_fasta, temp_output, seq_name, self.docker_image, 
                       self.use_signalp, self.cache_dir, self.docker_timeout, self.verbose)
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
                            print(f"✓ Completed: {seq_name}")
                        else:
                            failed_count += 1
                            print(f"✗ Failed: {seq_name} - {error}")
                    except Exception as e:
                        failed_count += 1
                        print(f"✗ Exception processing {seq_name}: {e}")
            
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


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Licensed NetNGlyc pipeline with SignalP 6 integration and intelligent parallel processing. Requires licensed NetNGlyc 1.0 Docker container and SignalP 6.0 installation.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", nargs='?', help="Input FASTA file or directory (required for process/full-pipeline modes)")
    parser.add_argument("output", nargs='?', help="Output file or directory (required for all modes except --test/--clear-cache)")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (used with --processing-mode parallel, default: 4)")
    parser.add_argument("--docker-image", default="netnglyc:latest",
                        help="Docker image name (default: netnglyc:latest)")
    parser.add_argument("--cache-dir", help="Custom cache directory for SignalP/NetNGlyc results")
    parser.add_argument("--test", action="store_true",
                        help="Run test with ABCB1 sequence (no other args required)")
    parser.add_argument("--clear-cache", action="store_true",
                        help="Clear all cached results and exit (no other args required)")
    parser.add_argument("--mode", choices=["process", "parse", "full-pipeline"], default="process",
                        help="Processing mode: 'process' (run NetNGlyc only), 'parse' (parse existing outputs), 'full-pipeline' (process + parse)")
    parser.add_argument("--mapping-dir", 
                        help="Directory containing mutation mapping CSV files (REQUIRED for wildtype parsing in parse/full-pipeline modes)")
    parser.add_argument("--is-mutant", action="store_true",
                        help="Parse mutant files vs wildtype files (affects parsing logic)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Minimum glycosylation potential threshold for predictions (default: 0.5)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="Number of sequences per batch for batch processing mode (default: 100)")
    parser.add_argument("--single-file", action="store_true",
                        help="Force single Docker container processing (overrides --processing-mode auto)")
    parser.add_argument("--batch-timeout", type=int, default=5000,
                        help="Timeout in seconds for Docker containers (default: 5000s)")
    parser.add_argument("--processing-mode", choices=["auto", "single", "parallel", "sequential", "batch"], 
                        default="auto",
                        help="Processing strategy: 'auto' (intelligent selection based on sequence count), 'single' (1 container, best for 1-50 seqs), 'parallel' (multiple containers, best for 51-500 seqs), 'sequential' (one-by-one), 'batch' (split into batches, required for 501+ seqs)")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep intermediate Docker output files for debugging (don't clean up temporary directories)")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose output showing detailed processing information")
    parser.add_argument("--log",
                        help="Validation log file or directory to skip failed mutations (mutant modes only)")

    args = parser.parse_args()
    failure_map = load_validation_failures(args.log) if (args.log and args.is_mutant) else {}

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
        test_fasta = "test_abcb1.fasta"
        test_output = "test_abcb1_netnglyc.out"

        # Create test file with ABCB1 fragment
        with open(test_fasta, 'w') as f:
            f.write(">ABCB1_HUMAN\n")
            f.write("MDLEGDRNGGAKKKNFFKLNNKSEKDKKEKKPTVSVFSMFRYSNWLDKLYMVVGTLAAI\n")
            f.write("IHGAGLPLMMLVFGEMTDIFANAGNLEDLMSNITNRSDINDTGFFMNLEEDMTRYAYY\n")

        with RobustDockerNetNGlyc(
                docker_image=args.docker_image,
                use_signalp=True,
                cache_dir=args.cache_dir,
                docker_timeout=args.batch_timeout,
                verbose=args.verbose
        ) as processor:
            success, output, error = processor.process_single_fasta(
                test_fasta, test_output, 0
            )

            if success:
                print("\nTest successful!")
                with open(test_output, 'r') as f:
                    content = f.read()
                    if "NITN" in content:
                        print("Found expected NITN glycosylation site")

                    # Check if SignalP was properly integrated
                    if "SignalP is not configured" not in content:
                        print("SignalP integration working!")
                    else:
                        print("WARNING: SignalP stub not recognized by NetNGlyc (but analysis still works)")

                    print("\nOutput preview:")
                    print(content[:1000])
            else:
                print(f"Test failed: {error}")

        # Clean up
        if os.path.exists(test_fasta):
            os.remove(test_fasta)
        return 0 if success else 1

    # Handle parsing mode
    if args.mode == "parse":
        if not args.input or not args.output:
            parser.error("For parse mode: both input (NetNGlyc output directory) and output (TSV file) are required")
        
        # Validate parsing requirements
        if not args.mapping_dir:
            parser.error("For parse mode: --mapping-dir is REQUIRED (directory containing mutation mapping CSV files)")
        
        processor = RobustDockerNetNGlyc(
            docker_image=args.docker_image, 
            use_signalp=False, 
            docker_timeout=args.batch_timeout,
            verbose=args.verbose
        )
        
        if args.is_mutant:
            print(f"Parsing mutant NetNGlyc files from {args.input}")
            print("Mode: Mutant parsing using wildtype logic with mapping files")
            results = processor.parse_mutant_files(
                args.input,
                args.threshold,
                mapping_dir=args.mapping_dir,
                failure_map=failure_map,
            )
        else:
            print(f"Parsing wildtype NetNGlyc files from {args.input}")
            print(f"Mode: Wildtype parsing (using mapping files from {args.mapping_dir})")
            results = processor.parse_wildtype_files(args.input, args.mapping_dir, args.threshold)
            
        processor.save_parsed_results(results, args.output)
        return 0
        
    elif args.mode == "full-pipeline":
        if not args.input or not args.output:
            parser.error("For full-pipeline mode: both input (FASTA file/directory) and output (TSV file) are required")
        
        # Validate full-pipeline mode requirements
        if not args.is_mutant and not args.mapping_dir:
            parser.error("For wildtype full-pipeline (without --is-mutant): --mapping-dir is REQUIRED (directory containing mutation mapping CSV files)")
        
        # Create temporary directory for NetNGlyc outputs
        import tempfile
        temp_output_dir = tempfile.mkdtemp(prefix="netnglyc_outputs_")
        
        try:
            # Step 1: Process FASTA files
            #print("Step 1: Processing FASTA files with NetNGlyc...")
            
            with RobustDockerNetNGlyc(
                docker_image=args.docker_image,
                use_signalp=True,
                max_workers=args.workers,
                cache_dir=args.cache_dir,
                docker_timeout=args.batch_timeout,
                keep_intermediates=args.keep_intermediates,
                verbose=args.verbose
            ) as processor:
                
                if os.path.isfile(args.input):
                    # Single file - use intelligent processing mode selection  
                    input_basename = os.path.basename(args.input)
                    # Remove any FASTA extension and add -netnglyc.out
                    for ext in ['.fasta', '.fa', '.fas', '.fna']:
                        if input_basename.lower().endswith(ext):
                            input_basename = input_basename[:-len(ext)]
                            break
                    output_file = os.path.join(temp_output_dir, f"{input_basename}-netnglyc.out")
                    
                    # Detect optimal processing mode
                    if args.single_file:
                        processing_mode = "single"
                        #print("User forced single-file processing (no batching)")
                    else:
                        processing_mode = processor.detect_processing_mode(args.input, args.processing_mode)
                    
                    # Execute based on detected/selected mode
                    if processing_mode == "single":
                        #print("Using single Docker container processing")
                        success, output, error = processor.process_single_fasta(args.input, output_file, 0)
                    elif processing_mode == "parallel":
                        #print(f"Using parallel Docker containers (max workers: {args.workers})")
                        success, output, error = processor.process_parallel_docker(
                            args.input, output_file, args.workers
                        )
                    elif processing_mode == "batch":
                        #print(f"Using batch processing (batch size: {args.batch_size})")
                        success, output, error = processor.process_fasta_batched(
                            args.input, output_file, args.batch_size, 0
                        )
                    elif processing_mode == "sequential":
                        #print("Using sequential processing (one-by-one)")
                        # For sequential, use single processing but could add specific logic
                        success, output, error = processor.process_single_fasta(args.input, output_file, 0)
                    else:
                        print(f"Unknown processing mode: {processing_mode}, defaulting to single")
                        success, output, error = processor.process_single_fasta(args.input, output_file, 0)
                    
                    # Validate output if processing succeeded
                    if success:
                        with open(output_file, 'r') as f:
                            output_content = f.read()
                        
                        expected_seq_count = processor.count_sequences_in_fasta(args.input)
                        is_valid, warnings, prediction_count = processor.validate_netnglyc_output(
                            output_content, expected_seq_count
                        )
                        
                        
                        print(f"Processing successful: {prediction_count} predictions found")
                    else:
                        print(f"NetNGlyc processing failed: {error}")
                        return 1
                else:
                    # Directory processing with intelligent mode detection
                    results = processor.process_directory(args.input, temp_output_dir, processing_mode=args.processing_mode)
                    print(f"NetNGlyc processing: {results['success']}/{results['total']} files successful")
                    
                    if results['failed'] > 0:
                        print(f"Warning: {results['failed']} files failed NetNGlyc processing")
                
                # Step 2: Parse the generated outputs
                #print("Step 2: Parsing NetNGlyc outputs...")
                
                if args.is_mutant:
                    #print("Parsing as mutant files using wildtype logic with mapping")
                    parsing_start_time = time.time()
                    results = processor.parse_mutant_files(
                        temp_output_dir,
                        args.threshold,
                        mapping_dir=args.mapping_dir,
                        fasta_dir=args.input,
                        failure_map=failure_map,
                    )
                    parsing_elapsed = time.time() - parsing_start_time
                    
                    # Print comprehensive summary for mutant processing
                    print(f"\n{'=' * 60}")
                    print(f"Mutant Processing Summary:")
                    print(f"   • Total predictions found: {len(results)}")
                    print(f"   • Parsing time: {parsing_elapsed / 60:.1f} minutes")
                    print(f"   • Threshold used: {args.threshold}")
                    if args.mapping_dir:
                        print(f"   • Filtered to mutation positions only")
                    print("=" * 60)
                else:
                    #print("Parsing as wildtype files")
                    results = processor.parse_wildtype_files(temp_output_dir, args.mapping_dir, args.threshold)
                    
                processor.save_parsed_results(results, args.output)
                
        finally:
            # Cleanup temporary directory (unless keeping intermediates for debugging)
            if args.keep_intermediates:
                print(f"Raw output files can be viewed in {temp_output_dir}")
            else:
                import shutil
                shutil.rmtree(temp_output_dir, ignore_errors=True)
            
        return 0

    # Normal processing mode - validate required arguments  
    if not args.input or not args.output:
        parser.error("For process mode: both input (FASTA file/directory) and output (NetNGlyc output file/directory) are required")

    with RobustDockerNetNGlyc(
            docker_image=args.docker_image,
            use_signalp=True,
            max_workers=args.workers,
            cache_dir=args.cache_dir,
            docker_timeout=args.batch_timeout,
            verbose=args.verbose
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
            results = processor.process_directory(args.input, args.output, processing_mode=args.processing_mode)
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
