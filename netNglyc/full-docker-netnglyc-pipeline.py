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
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from dependencies.utility import read_fasta, get_mutation_data_bioAccurate

# Import utility functions  
sys.path.append(os.path.join(os.path.dirname(__file__), 'dependencies'))


def _process_single_sequence_worker(args):
    """
    Worker function for parallel processing - must be at module level for pickling
    """
    temp_fasta, temp_output, seq_name, docker_image, use_signalp, cache_dir, docker_timeout = args
    
    worker_processor = RobustDockerNetNGlyc(
        docker_image=docker_image,
        use_signalp=use_signalp,
        max_workers=1,
        cache_dir=cache_dir,
        docker_timeout=docker_timeout
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

    def __init__(self, cache_dir=None):
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".signalp6_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        self.last_output_dir = None  # Store the output directory path

    def check_signalp6_available(self):
        """Check if SignalP 6 is available on the host"""
        try:
            result = subprocess.run(
                ["signalp6", "--version"],
                capture_output=True,
                timeout=5,
                text=True
            )
            if result.returncode == 0:
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

        if not self.check_signalp6_available():
            print("SignalP 6.0 not available - will use compatible stub within Docker container")
            return {}, None

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
            # Clean up temp dir if we created one
            if temp_dir and temp_dir != signalp_output_dir:
                shutil.rmtree(temp_dir)

        return results, None


class RobustDockerNetNGlyc:
    """
    Docker NetNGlyc with integrated SignalP 6 preprocessing on host
    """

    def __init__(self, docker_image="netnglyc:latest", use_signalp=True,
                 max_workers=4, cache_dir=None, docker_timeout=600, keep_intermediates=False):
        self.docker_image = docker_image
        self.max_workers = max_workers
        self.temp_dir = tempfile.mkdtemp(prefix="netnglyc_")
        self.use_signalp = use_signalp
        self.docker_timeout = docker_timeout
        self.keep_intermediates = keep_intermediates

        # Initialize SignalP handler
        self.signalp_handler = SignalP6Handler(cache_dir=cache_dir) if use_signalp else None

        # Initialize cache
        if cache_dir is None:
            cache_dir = os.path.join(os.path.expanduser("~"), ".netnglyc_cache")
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        # Detect platform once during initialization
        self.platform_args = get_docker_platform_args()
        if self.platform_args:
            print(f"Detected ARM64 architecture - Docker will use emulation mode")

        # Test Docker setup (non-fatal - allow fallback during processing)
        self._test_docker_availability()

    def _test_docker_availability(self):
        """Test if Docker and NetNGlyc image are available (non-fatal)"""
        try:
            # Test Docker
            result = subprocess.run(["docker", "--version"], capture_output=True, timeout=5)
            if result.returncode != 0:
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
                print(f"Docker image '{self.docker_image}' not found - will fallback to native stub when needed")
                return False
                
            print(f"Docker NetNGlyc ready: {self.docker_image}")
            return True

        except Exception as e:
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
        Run NetNGlyc in Docker container with SignalP results mounted
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
                    return result.stdout
                elif "No Asparagines in the input sequences" in result.stdout:
                    # No N residues to analyze
                    return result.stdout

            # Check if platform warning but still got output
            if "WARNING" in result.stderr and "platform" in result.stderr:
                if result.stdout and "Predictions" in result.stdout:
                    return result.stdout

            raise Exception(f"Docker NetNGlyc failed with return code {result.returncode}")

        except Exception as e:
            print(f"Docker NetNGlyc failed ({e}), attempting native stub fallback...")
            return self._run_native_netnglyc_stub(fasta_file, signalp_results=None)
        finally:
            # Clean up work directory
            if os.path.exists(work_dir):
                shutil.rmtree(work_dir)

    def _run_native_netnglyc_stub(self, fasta_file, signalp_results=None):
        """
        Run netnglyc_stub directly on host system (license-free alternative)
        
        Args:
            fasta_file: Path to FASTA file to analyze
            signalp_results: Dict of SignalP 6.0 results (optional)
            
        Returns:
            str: NetNGlyc-compatible output string, or None if failed
        """
        try:
            # Find netnglyc_stub script
            script_dir = os.path.dirname(os.path.abspath(__file__))
            stub_path = os.path.join(script_dir, "netnglyc_stub")
            
            if not os.path.exists(stub_path):
                print(f"NetNGlyc stub not found at: {stub_path}")
                return None
            
            print(f"Running native NetNGlyc stub on {os.path.basename(fasta_file)}")
            
            # Execute netnglyc_stub
            cmd = ["python3", stub_path, fasta_file]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode == 0 and result.stdout:
                stub_output = result.stdout
                
                # Enhance output with SignalP information if available
                if signalp_results:
                    enhanced_output = self._integrate_signalp_with_stub_output(
                        stub_output, signalp_results
                    )
                    return enhanced_output
                else:
                    return stub_output
            else:
                print(f"NetNGlyc stub failed: {result.stderr[:200] if result.stderr else 'No error message'}")
                return None
                
        except subprocess.TimeoutExpired:
            print("NetNGlyc stub execution timed out")
            return None
        except Exception as e:
            print(f"NetNGlyc stub execution error: {e}")
            return None

    def _integrate_signalp_with_stub_output(self, stub_output, signalp_results):
        """
        Integrate SignalP 6.0 results into netnglyc_stub output for enhanced accuracy
        
        Args:
            stub_output: Original netnglyc_stub output
            signalp_results: Dict of SignalP predictions per sequence
            
        Returns:
            str: Enhanced output with SignalP information
        """
        try:
            lines = stub_output.split('\n')
            enhanced_lines = []
            
            # Add SignalP header information
            enhanced_lines.append("##############################################################################")
            enhanced_lines.append("# NetNGlyc stub with SignalP 6.0 integration - Enhanced N-glycosylation prediction")
            enhanced_lines.append("##############################################################################")
            enhanced_lines.append("")
            
            # Process each line and enhance with SignalP context
            in_header = True
            for line in lines:
                if line.startswith("# Predictions for N-Glycosylation sites"):
                    enhanced_lines.append(line)
                    enhanced_lines.append("")
                    
                    # Add SignalP summary
                    if signalp_results:
                        enhanced_lines.append("# SignalP 6.0 predictions:")
                        for seq_id, info in signalp_results.items():
                            if info['has_signal']:
                                enhanced_lines.append(f"# {seq_id}: Signal peptide detected, "
                                                    f"cleavage at position {info['cleavage_site']} "
                                                    f"(probability: {info['probability']:.3f})")
                            else:
                                enhanced_lines.append(f"# {seq_id}: No signal peptide detected "
                                                    f"(probability: {info['probability']:.6f})")
                        enhanced_lines.append("")
                    in_header = False
                elif not line.startswith("##############################################################################") and not (in_header and line.startswith("#")):
                    enhanced_lines.append(line)
            
            return '\n'.join(enhanced_lines)
            
        except Exception as e:
            print(f"Warning: Failed to integrate SignalP results: {e}")
            return stub_output

    def process_single_fasta(self, fasta_file, output_file, worker_id=0):
        """
        Process a single FASTA file with optional SignalP 6 preprocessing

        1. Run SignalP 6 on host (if enabled)
        2. Pass ORIGINAL FASTA and SignalP results to Docker
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
    
            # Step 2: Run NetNGlyc (Docker or native stub) with ORIGINAL FASTA and SignalP results
            if getattr(self, 'use_native_stub', False):
                #print(f"Worker {worker_id}: Running native NetNGlyc stub...")
                netnglyc_output = self._run_native_netnglyc_stub(fasta_file, signalp_results)
            else:
                #print(f"Worker {worker_id}: Running NetNGlyc in Docker...")
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
                # Find the results table
                if "SeqName" in line and "Position" in line and "Potential" in line:
                    in_results = True
                    skip_next_line = True  # Skip the second header line
                    continue
                
                # Skip the second header line  
                if skip_next_line:
                    skip_next_line = False
                    continue
                    
                # Handle separator lines - first one starts data, second one ends it
                if in_results and line.strip().startswith('---'):
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

    def split_fasta_into_batches(self, fasta_file, batch_size=100):
        """Split a FASTA file into smaller batches for processing"""
        import math
        
        try:
            # Read all sequences from FASTA file
            sequences = read_fasta(fasta_file)
            total_sequences = len(sequences)
            
            if total_sequences == 0:
                return []
            
            print(f"Splitting {total_sequences} sequences into batches of {batch_size}")
            
            # Calculate number of batches needed
            num_batches = math.ceil(total_sequences / batch_size)
            
            batch_files = []
            sequence_items = list(sequences.items())
            
            for i in range(num_batches):
                start_idx = i * batch_size
                end_idx = min((i + 1) * batch_size, total_sequences)
                
                # Create batch sequences
                batch_sequences = dict(sequence_items[start_idx:end_idx])
                batch_count = len(batch_sequences)
                
                # Create temporary batch file
                batch_filename = fasta_file.replace('.fasta', f'_batch{i+1}.fasta')
                
                with open(batch_filename, 'w') as f:
                    for seq_name, sequence in batch_sequences.items():
                        f.write(f">{seq_name}\n")
                        # Write sequence in lines of 80 characters
                        for j in range(0, len(sequence), 80):
                            f.write(sequence[j:j+80] + "\n")
                
                batch_files.append(batch_filename)
                print(f"Created batch {i+1}/{num_batches}: {batch_filename} ({batch_count} sequences)")
            
            return batch_files
            
        except Exception as e:
            print(f"Error splitting FASTA file {fasta_file}: {e}")
            return []

    def process_fasta_batched(self, fasta_file, output_base, batch_size=100, worker_id=0):
        """Process a FASTA file in batches to handle large numbers of sequences"""
        try:
            # Split FASTA into batches
            batch_files = self.split_fasta_into_batches(fasta_file, batch_size)
            
            if not batch_files:
                return False, None, "Failed to create batch files"
            
            # If only one batch, use regular processing
            if len(batch_files) == 1:
                #print(f"Worker {worker_id}: Single batch processing")
                success, output_file, error = self.process_single_fasta(
                    batch_files[0], output_base, worker_id
                )
                # Clean up batch file
                if os.path.exists(batch_files[0]):
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
                        print(f"   Batch {i+1} saved: {os.path.basename(batch_output)}")
                    else:
                        print(f"Batch {i+1} failed: {error}")
                        # Clean up batch files
                        for bf in batch_files:
                            if os.path.exists(bf):
                                os.remove(bf)
                        return False, output_base, f"Batch {i+1} failed: {error}"
                        
                finally:
                    # Clean up this batch file
                    if os.path.exists(batch_file):
                        os.remove(batch_file)
            
            # Create combined output for backwards compatibility (optional)
            #print(f"Worker {worker_id}: Creating combined output for compatibility")
            combined_output = self.combine_batch_outputs(batch_outputs, output_base)
            
            # PRESERVE all batch files for parsing - do NOT delete them
            #print(f"   Preserved {len(batch_outputs)} batch files for parsing:")
            for batch_output in batch_outputs:
                if os.path.exists(batch_output):
                    print(f"      - {os.path.basename(batch_output)}")
            
            # No longer cleaning up batch files - they're needed for parsing!
            
            if combined_output:
                return True, output_base, None
            else:
                return False, output_base, "Failed to combine batch outputs"
                
        except Exception as e:
            # Clean up any remaining batch files
            try:
                if 'batch_files' in locals():
                    for bf in batch_files:
                        if os.path.exists(bf):
                            os.remove(bf)
            except:
                pass
            return False, output_base, f"Batch processing error: {e}"

    def combine_batch_outputs(self, batch_output_files, final_output_file):
        """
        Combine multiple NetNGlyc batch outputs into a single file for backwards compatibility
        
        IMPORTANT: This creates a combined file but does NOT delete the individual batch files.
        The individual batch files are preserved for parsing since they contain all predictions.
        
        Args:
            batch_output_files: List of individual batch NetNGlyc output files  
            final_output_file: Path for combined output file
            
        Returns:
            bool: True if combination successful
        """
        try:
            if not batch_output_files:
                return False
            
            print(f"Combining {len(batch_output_files)} batch outputs...")
            
            # Count total sequences for header
            total_sequences = 0
            all_sequence_sections = []
            all_prediction_lines = []
            total_predictions = 0
            
            for i, batch_file in enumerate(batch_output_files):
                try:
                    with open(batch_file, 'r') as f:
                        content = f.read()
                    
                    batch_seq_count = content.count('Name:')
                    total_sequences += batch_seq_count
                    print(f"   {os.path.basename(batch_file)}: {batch_seq_count} sequences")
                    
                    # Extract complete sequence sections (from Name: to first ---)
                    lines = content.split('\n')
                    sequence_section = []
                    prediction_section = []
                    
                    collecting_sequences = False
                    collecting_predictions = False
                    found_results_header = False
                    separators_seen = 0
                    
                    for line in lines:
                        # Collect sequence information
                        if line.startswith('Name:'):
                            collecting_sequences = True
                            sequence_section.append(line)
                        elif collecting_sequences:
                            # Check if we've hit the results section
                            if "SeqName" in line and "Position" in line and "Potential" in line:
                                # Hit results header - stop collecting sequences
                                collecting_sequences = False
                            elif line.strip().startswith('---') and not collecting_predictions:
                                # Hit separator before results - stop collecting sequences
                                collecting_sequences = False
                            else:
                                sequence_section.append(line)
                        
                        # Collect prediction results
                        if "SeqName" in line and "Position" in line and "Potential" in line:
                            found_results_header = True
                            continue
                        
                        if found_results_header and line.strip().startswith('---'):
                            separators_seen += 1
                            if separators_seen == 1:
                                collecting_predictions = True  # Start collecting after first separator
                                continue
                            elif separators_seen == 2:
                                collecting_predictions = False  # Stop at second separator
                                break
                        
                        if collecting_predictions and line.strip():
                            # Skip header continuation lines
                            if 'agreement' in line and 'result' in line and len(line.split()) <= 2:
                                continue
                            # Check for valid prediction line
                            parts = line.strip().split()
                            if len(parts) >= 4 and not line.startswith('#'):
                                try:
                                    # Validate it's a real prediction (position should be a number)
                                    int(parts[1])
                                    prediction_section.append(line.strip())
                                    total_predictions += 1
                                except (ValueError, IndexError):
                                    continue
                    
                    # Add this batch's sections to the combined collections
                    if sequence_section:
                        all_sequence_sections.extend(sequence_section)
                        all_sequence_sections.append("")  # Add spacing between batches
                    
                    all_prediction_lines.extend(prediction_section)
                    
                    batch_predictions = len(prediction_section)
                    print(f"   Batch {i+1}: {batch_predictions} predictions extracted")
                    
                except Exception as e:
                    print(f"   Error processing batch {i+1} ({batch_file}): {e}")
                    continue
            
            print(f"   Total sequences: {total_sequences}, Total predictions: {total_predictions}")
            
            # Build the final combined output
            combined_content = []
            
            # 1. Header
            combined_content.extend([
                f"# Predictions for N-Glycosylation sites in {total_sequences} sequences",
                "",
                "netNglyc: SignalP is not configured, no signal peptide predictions are made",
                ""
            ])
            
            # 2. All sequence information
            combined_content.extend(all_sequence_sections)
            
            # 3. Results table
            if all_prediction_lines:
                combined_content.extend([
                    "(Threshold=0.5)",
                    "----------------------------------------------------------------------",
                    "SeqName        Position      Potential          Jury          N-Glyc",
                    "\t\t\t\t              agreement       result",
                    "----------------------------------------------------------------------"
                ])
                combined_content.extend(all_prediction_lines)
                combined_content.append("----------------------------------------------------------------------")
                print(f"   Added {len(all_prediction_lines)} prediction lines to results table")
            else:
                combined_content.extend([
                    "(Threshold=0.5)",
                    "----------------------------------------------------------------------",
                    "No sites predicted",
                    "----------------------------------------------------------------------"
                ])
                print(f"   WARNING: No predictions found - adding 'No sites predicted'")
            
            combined_content.append("")  # Final newline
            
            # Write the combined file
            try:
                with open(final_output_file, 'w') as f:
                    f.write('\n'.join(combined_content))
                
                print(f"   Successfully wrote combined output to {os.path.basename(final_output_file)}")
                print(f"   Final stats: {total_sequences} sequences, {total_predictions} predictions")
                
                # Simple validation
                with open(final_output_file, 'r') as f:
                    validation_content = f.read()
                
                actual_seq_count = validation_content.count('Name:')
                print(f"   Validation: {actual_seq_count} sequences in final file (expected: {total_sequences})")
                
                if actual_seq_count != total_sequences:
                    print(f"   WARNING: Sequence count mismatch! Expected {total_sequences}, got {actual_seq_count}")
                
                return True
                
            except Exception as e:
                print(f"   Error writing combined output: {e}")
                return False
            
        except Exception as e:
            print(f"Error combining batch outputs: {e}")
            return False

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
                                # If we haven't detected the current sequence yet, use seq_name
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

    def parse_mutant_files(self, input_dir, threshold=0.5, mapping_dir=None, fasta_dir=None):
        """Parse mutant NetNGlyc files using sequence names from predictions"""
        import os
        import csv
        import logging
        from datetime import datetime
        from pathlib import Path
        
        results = []
        
        if mapping_dir is None:
            print("ERROR: mapping_dir is required for mutant parsing")
            return results
        
        # Find NetNGlyc output files: both regular and batch files
        # Pattern 1: {GENE}_aa-netnglyc.out (regular files)
        regular_files = list(Path(input_dir).glob("*_aa-netnglyc.out"))
        # Pattern 2: {GENE}_aa-netnglyc-batch-*.out (batch files)
        batch_files = list(Path(input_dir).glob("*_aa-netnglyc-batch-*.out"))
        
        # Group files by gene
        gene_files = {}
        
        # Process regular files
        for file_path in regular_files:
            filename = file_path.stem  # Remove .out extension
            if '_aa-netnglyc' in filename:
                gene = filename.replace('_aa-netnglyc', '')
                if gene not in gene_files:
                    gene_files[gene] = []
                gene_files[gene].append(file_path)
        
        # Process batch files  
        for file_path in batch_files:
            filename = file_path.stem  # Remove .out extension
            # Extract gene from batch filename: ABCB1_aa-netnglyc-batch-1 -> ABCB1
            if '_aa-netnglyc-batch-' in filename:
                gene = filename.split('_aa-netnglyc-batch-')[0]
                if gene not in gene_files:
                    gene_files[gene] = []
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
        
        # Process each gene's files
        for gene, gene_file_list in gene_files.items():
            print(f"\nProcessing {gene} using sequence-based mutant parsing ({len(gene_file_list)} files)")
            
            # Set up debug logging for this gene (once per gene)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            log_file = os.path.join(input_dir, f"netnglyc_debug_{gene}_{timestamp}.log")
            
            # Configure logging
            logger = logging.getLogger(f"netnglyc_{gene}")
            logger.setLevel(logging.DEBUG)
            # Clear any existing handlers
            logger.handlers.clear()
            
            # File handler for detailed logging
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(message)s')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            
            logger.info(f"=== NetNGlyc Parsing Debug Log for {gene} ===")
            logger.info(f"Threshold: {threshold}")
            logger.info(f"Processing {len(gene_file_list)} files for this gene")
            
            # Load mapping file for this gene (once per gene)
            mapping_file = os.path.join(mapping_dir, f"{gene}_nt_to_aa_mapping.csv")
            if not os.path.exists(mapping_file):
                print(f"Warning: Mapping file not found for {gene}: {mapping_file}")
                continue
            
            # Read mapping file into dictionary: {ntposnt: aaposaa}
            mapping_dict = {}
            try:
                with open(mapping_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        ntposnt = row['mutant']  # e.g., "C3066A"
                        aaposaa = row['aamutant']  # e.g., "S1022N"
                        mapping_dict[ntposnt] = aaposaa
                
                print(f"Loaded mapping for {len(mapping_dict)} mutations")
                
                # DEBUG: Show a few mapping examples
                sample_mappings = list(mapping_dict.items())[:3]
                
            except Exception as e:
                print(f"Error reading mapping file {mapping_file}: {e}")
                continue
            
            # Process all files for this gene (both regular and batch)
            gene_predictions = []
            for file_index, file_path in enumerate(gene_file_list):
                file_type = "batch" if "-batch-" in file_path.stem else "regular"
                print(f"  Processing file {file_index+1}/{len(gene_file_list)}: {file_path.name} ({file_type})")
                logger.info(f"Processing {file_type} file: {file_path.name}")
                
                # Get list of mutation sequence names from the FASTA file if provided
                mutant_sequences = []
                if fasta_dir:
                    fasta_file = os.path.join(fasta_dir, f"{gene}_aa.fasta")
                    if os.path.exists(fasta_file):
                        try:
                            fasta_data = read_fasta(fasta_file)
                            mutant_sequences = list(fasta_data.keys())
                            print(f"Found {len(mutant_sequences)} sequences in FASTA file")
                        except Exception as e:
                            print(f"Error reading FASTA file {fasta_file}: {e}")
                
                # Parse NetNGlyc predictions from this file
                try:
                    file_predictions = self.parse_simple_netnglyc_output(file_path, threshold)
                    print(f"    Parsed {len(file_predictions)} predictions from {file_path.name}")
                    
                    # Add to gene's total predictions
                    gene_predictions.extend(file_predictions)
                    
                    # DEBUG: Show some sample predictions
                    if file_predictions:
                        sample_preds = file_predictions[:2]
                        positions = [f"pos={p['position']}" for p in sample_preds]
                        print(f"    Sample predictions: {positions}")
                    
                except Exception as e:
                    print(f"    Error parsing NetNGlyc output {file_path}: {e}")
                    logger.error(f"Failed to parse {file_path.name}: {e}")
                    continue
            
            # Now process all predictions for this gene together
            print(f"  Total predictions for {gene}: {len(gene_predictions)}")
            logger.info(f"Total predictions collected from all files: {len(gene_predictions)}")
            
            # Process each prediction individually using its sequence name
            processed_mutations = set()
            total_processed = 0
            matches_found = 0
            stop_codons_skipped = 0
            below_threshold = 0
            position_mismatches = 0
            
            logger.info(f"Processing {len(gene_predictions)} predictions")
            
            for pred in gene_predictions:
                total_processed += 1
                try:
                    seq_name = pred['seq_name']  # e.g., "ABCB1-A1002T"
                    
                    # Skip stop codons in sequence name
                    if 'Sto' in seq_name:
                        continue
                    
                    # Extract mutation ID from sequence name: ABCB1-A1002T -> A1002T
                    if '-' in seq_name:
                        mutation_id = seq_name.rsplit('-', 1)[1]  # A1002T
                    else:
                        continue
                    
                    # Look up this mutation in the mapping
                    if mutation_id not in mapping_dict:
                        continue
                    
                    aaposaa = mapping_dict[mutation_id]  # e.g., "V334V" or "790Sto"
                    
                    # Skip stop codons in mapping result
                    if 'Sto' in aaposaa:
                        stop_codons_skipped += 1
                        logger.debug(f"SKIP - Stop codon: {mutation_id} -> {aaposaa}")
                        continue
                    
                    # Extract amino acid position from aaposaa
                    aa_pos, _ = get_mutation_data_bioAccurate(aaposaa)  # e.g., 334
                    
                    # Log mapping details to file
                    logger.debug(f"MAPPING: {mutation_id} -> {aaposaa} -> position {aa_pos}")
                    
                    # Get all prediction positions for this sequence
                    seq_predictions = [p['position'] for p in gene_predictions if p['seq_name'] == seq_name]
                    logger.debug(f"POSITIONS: {seq_name} has predictions at {seq_predictions}")
                    
                    # Check if this prediction is at the mutation position (RESTORE ORIGINAL LOGIC)
                    if pred['position'] == aa_pos:
                        matches_found += 1
                        pkey = f"{gene}-{mutation_id}"  # e.g., "ABCB1-A1002T"
                        
                        # Create result entry for this prediction
                        results.append({
                            'pkey': pkey,
                            'Gene': gene,
                            'pos': pred['position'],
                            'Sequon': pred['sequon'],
                            'potential': pred['potential'],
                            'jury_agreement': pred.get('jury_agreement', ''),
                            'n_glyc_result': pred.get('n_glyc_result', '')
                        })
                        
                        # Track that we found this mutation
                        processed_mutations.add(mutation_id)
                        
                        # Log match to file and stdout
                        logger.info(f"MATCH: {seq_name} at pos {pred['position']} (potential={pred['potential']:.4f})")
                        #print(f"   MATCH: {seq_name} at position {pred['position']}")
                        
                    else:
                        position_mismatches += 1
                        logger.debug(f"NO_MATCH: {seq_name} prediction_pos={pred['position']} != mutation_pos={aa_pos}")
                        
                        # Log detailed comparison for sample mutations
                        if mutation_id in ['A1002T', 'C3066A', 'T3435C', 'A2895T', 'T1758C']:
                            logger.info(f"SAMPLE: {mutation_id} mutation_pos={aa_pos}, predictions={seq_predictions} -> NO MATCH at {pred['position']}")
                
                except Exception as e:
                    print(f"Error processing prediction for {pred.get('seq_name', 'unknown')}: {e}")
                    continue
            
            # Summary statistics to stdout
            print(f"   Processed {total_processed} predictions from {len(set([p['seq_name'].rsplit('-', 1)[1] for p in gene_predictions if '-' in p['seq_name']]))} unique mutations")
            print(f"   Found {matches_found} predictions at mutation positions")
            print(f"   Skipped {stop_codons_skipped} stop codons")
            print(f"   📝 Debug log: {log_file}")
            
            # Final summary to log file
            logger.info(f"=== FINAL SUMMARY ===")
            logger.info(f"Total predictions processed: {total_processed}")
            logger.info(f"Matches found at mutation positions: {matches_found}")
            logger.info(f"Position mismatches: {position_mismatches}")
            logger.info(f"Stop codons skipped: {stop_codons_skipped}")
            logger.info(f"Unique mutations with matches: {len(processed_mutations)}")
            
            # Close log handler
            file_handler.close()
            logger.removeHandler(file_handler)
        
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
            
            # Extract gene from filename: ABCB1_aa-netnglyc or ABCB1_wt-netnglyc
            gene = None
            if '_aa-netnglyc' in filename:
                gene = filename.replace('_aa-netnglyc', '')
            elif '_wt-netnglyc' in filename:
                gene = filename.replace('_wt-netnglyc', '')
            else:
                # Skip files that don't match expected patterns
                continue
            
            # Load mapping for this gene
            mapping_file = os.path.join(mapping_dir, f"{gene}_nt_to_aa_mapping.csv")
            if os.path.exists(mapping_file):
                mapping = self.load_nt_to_aa_mapping(mapping_file)
                
                # Parse NetNGlyc predictions once
                predictions = self.parse_netnglyc_output(file_path, threshold)
                
                # For each mutation, filter predictions at that position
                for ntposnt, mutation_aa_pos in mapping.items():
                    for pred in predictions:
                        if pred['position'] == mutation_aa_pos:
                            pkey = f"{gene}-{ntposnt}"
                            results.append({
                                'pkey': pkey,
                                'Gene': gene,
                                'pos': pred['position'],
                                'Sequon': pred['sequon'],
                                'potential': pred['potential'],
                                'jury_agreement': pred['jury_agreement'],
                                'n_glyc_result': pred['n_glyc_result']
                            })
        
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

    def process_directory(self, input_dir, output_dir, pattern="*.fasta", processing_mode="auto"):
        """
        Process all FASTA files in directory with intelligent processing mode detection
        Applies intelligent processing strategy to each file based on sequence count
        """
        os.makedirs(output_dir, exist_ok=True)

        # Find all FASTA files
        input_files = list(Path(input_dir).glob(pattern))
        if not input_files:
            print(f"No files matching {pattern} found")
            return {"total": 0, "success": 0, "failed": 0}

        print(f"Processing {len(input_files)} files from directory: {input_dir}")
        if self.use_signalp and self.signalp_handler.check_signalp6_available():
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
                            
                            if warnings:
                                print(f"WARNING: {strategy['file_name']}: {len(warnings)} validation warnings")
                                for warning in warnings:
                                    print(f"    • {warning}")
                        
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
            docker_timeout=self.docker_timeout
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
                # For directory processing, we use a smaller worker count per file
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
        
        # Check for successful execution
        if "Predictions for N-Glycosylation sites" not in output_content:
            return False, ["Missing NetNGlyc header - execution may have failed"], 0
        
        # Count actual sequences processed
        sequence_names = output_content.count("Name:")
        if expected_sequences and sequence_names != expected_sequences:
            warnings.append(f"Expected {expected_sequences} sequences, found {sequence_names}")
        
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
        
        is_valid = "Predictions for N-Glycosylation sites" in output_content
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
            
            # Create temporary files for individual sequences
            temp_files = []
            output_files = []
            
            for i, (seq_name, sequence) in enumerate(sequences.items()):
                # Create individual FASTA file
                temp_fasta = os.path.join(self.temp_dir, f"seq_{i}_{seq_name.replace('/', '_')}.fasta")
                temp_output = os.path.join(self.temp_dir, f"seq_{i}_output.out")
                
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
                       self.use_signalp, self.cache_dir, self.docker_timeout)
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
            # Clean up temporary files
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
            
            # Create header section
            combined_lines.extend([
                f"# Predictions for N-Glycosylation sites in {total_sequences} sequences",
                "",
                "netNglyc: SignalP is not configured, no signal peptide predictions are made",
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
        description="Complete NetNGlyc pipeline with SignalP 6 integration and intelligent parallel processing",
        epilog="""
USAGE EXAMPLES:

1. Test mode (no other args required):
   %(prog)s --test

2. Process single FASTA file (automatic mode selection):
   %(prog)s input.fasta output-netnglyc.out

3. Full pipeline with automatic processing mode (wildtype sequences):
   %(prog)s --mode full-pipeline --mapping-dir ./mutations/combined/aa/ input.fasta results.tsv

4. Full pipeline with automatic processing mode (mutant sequences):
   %(prog)s --mode full-pipeline --is-mutant input.fasta results.tsv

5. Force parallel processing with 8 workers:
   %(prog)s --mode full-pipeline --processing-mode parallel --workers 8 input.fasta results.tsv

6. Parse existing NetNGlyc outputs (mutant mode):
   %(prog)s --mode parse --is-mutant netnglyc-outputs/ parsed_results.tsv

7. Parse existing NetNGlyc outputs (wildtype mode - requires mapping):
   %(prog)s --mode parse --mapping-dir ./mutations/combined/aa/ netnglyc-outputs/ parsed_results.tsv

PROCESSING MODES:
- auto: Intelligent selection (1-50 seqs→single, 51-500→parallel, 501+→batch)
- single: One Docker container (optimal for small datasets)
- parallel: Multiple Docker containers (fast for medium datasets) 
- batch: Sequential batching (required for very large datasets)

REQUIRED ARGUMENTS BY MODE:
- process: input, output
- parse (mutant): input, output, --is-mutant
- parse (wildtype): input, output, --mapping-dir  
- full-pipeline (mutant): input, output, --is-mutant
- full-pipeline (wildtype): input, output, --mapping-dir
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", nargs='?', help="Input FASTA file or directory (required for process/full-pipeline modes)")
    parser.add_argument("output", nargs='?', help="Output file or directory (required for all modes except --test/--clear-cache)")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (used with --processing-mode parallel, default: 4)")
    parser.add_argument("--no-signalp", action="store_true",
                        help="Skip SignalP 6 preprocessing (faster but less accurate)")
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
                        help="Directory containing nt_to_aa_mapping.csv files (REQUIRED for wildtype parsing in parse/full-pipeline modes)")
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

    args = parser.parse_args()

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
                use_signalp=not args.no_signalp,
                cache_dir=args.cache_dir,
                docker_timeout=args.batch_timeout
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
            parser.error("For parse mode: --mapping-dir is REQUIRED (directory containing {GENE}_nt_to_aa_mapping.csv files)")
        
        processor = RobustDockerNetNGlyc(
            docker_image=args.docker_image, 
            use_signalp=False, 
            docker_timeout=args.batch_timeout
        )
        
        if args.is_mutant:
            print(f"Parsing mutant NetNGlyc files from {args.input}")
            print("Mode: Mutant parsing using wildtype logic with mapping files")
            results = processor.parse_mutant_files(args.input, args.threshold, mapping_dir=args.mapping_dir)
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
            parser.error("For wildtype full-pipeline (without --is-mutant): --mapping-dir is REQUIRED (directory containing {GENE}_nt_to_aa_mapping.csv files)")
        
        # Create temporary directory for NetNGlyc outputs
        import tempfile
        temp_output_dir = tempfile.mkdtemp(prefix="netnglyc_outputs_")
        
        try:
            # Step 1: Process FASTA files
            #print("Step 1: Processing FASTA files with NetNGlyc...")
            
            with RobustDockerNetNGlyc(
                docker_image=args.docker_image,
                use_signalp=not args.no_signalp,
                max_workers=args.workers,
                cache_dir=args.cache_dir,
                docker_timeout=args.batch_timeout,
                keep_intermediates=args.keep_intermediates
            ) as processor:
                
                if os.path.isfile(args.input):
                    # Single file - use intelligent processing mode selection
                    output_file = os.path.join(temp_output_dir, 
                                             os.path.basename(args.input).replace('.fasta', '-netnglyc.out'))
                    
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
                        # For sequential, we'll use single processing but could add specific logic
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
                        
                        if warnings:
                            print("Output validation warnings:")
                            for warning in warnings:
                                print(f"  WARNING: {warning}")
                        
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
                    results = processor.parse_mutant_files(temp_output_dir, args.threshold, mapping_dir=args.mapping_dir, fasta_dir=args.input)
                else:
                    #print("Parsing as wildtype files")
                    results = processor.parse_wildtype_files(temp_output_dir, args.mapping_dir, args.threshold)
                    
                processor.save_parsed_results(results, args.output)
                
        finally:
            # Cleanup temporary directory (unless keeping intermediates for debugging)
            if args.keep_intermediates:
                #print(f"DEBUG: Keeping intermediate Docker outputs in: {temp_output_dir}")
                print(f"    You can examine individual NetNGlyc output files for debugging")
                print(f"    Remember to clean up manually when done: rm -rf {temp_output_dir}")
            else:
                import shutil
                shutil.rmtree(temp_output_dir, ignore_errors=True)
            
        return 0

    # Normal processing mode - validate required arguments  
    if not args.input or not args.output:
        parser.error("For process mode: both input (FASTA file/directory) and output (NetNGlyc output file/directory) are required")

    with RobustDockerNetNGlyc(
            docker_image=args.docker_image,
            use_signalp=not args.no_signalp,
            max_workers=args.workers,
            cache_dir=args.cache_dir,
            docker_timeout=args.batch_timeout
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