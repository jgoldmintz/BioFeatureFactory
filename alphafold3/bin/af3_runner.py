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
AlphaFold3 execution wrapper.

Generates AF3 input JSON files and executes predictions in:
- local: Direct subprocess execution (requires GPU)
- batch: Generate SLURM submission scripts
- cloud: Generate GCP Batch configuration
"""

import json
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from enum import Enum
import hashlib


class ExecutionMode(Enum):
    LOCAL = "local"
    BATCH = "batch"
    CLOUD = "cloud"


class PredictionBackend(Enum):
    AF3 = "af3"
    BOLTZ = "boltz"


@dataclass
class AF3Input:
    """AlphaFold3 input specification for RNA-protein complex."""
    name: str
    rna_sequence: str
    protein_sequence: str
    rna_chain_id: str = "A"
    protein_chain_id: str = "B"
    protein_msa: Optional[str] = None  # A3M format MSA content

    def to_json_dict(self) -> Dict[str, Any]:
        """Convert to AF3 input JSON format."""
        # RNA chain
        rna_entry = {
            "rna": {
                "id": self.rna_chain_id,
                "sequence": self.rna_sequence.replace('T', 'U')  # Ensure RNA
            }
        }

        # Protein chain with optional MSA
        protein_entry = {
            "protein": {
                "id": self.protein_chain_id,
                "sequence": self.protein_sequence
            }
        }

        # Add MSA if provided
        if self.protein_msa:
            protein_entry["protein"]["unpairedMsa"] = self.protein_msa

        return {
            "name": self.name,
            "modelSeeds": [1],
            "sequences": [rna_entry, protein_entry]
        }

    def get_hash(self) -> str:
        """Get unique hash for this input (for caching)."""
        content = f"{self.rna_sequence}_{self.protein_sequence}"
        return hashlib.md5(content.encode()).hexdigest()[:12]


@dataclass
class AF3Job:
    """AF3 job configuration."""
    job_id: str
    input_json_path: Path
    output_dir: Path
    wt_input: AF3Input
    mut_input: Optional[AF3Input] = None
    status: str = "pending"
    result_path: Optional[Path] = None


@dataclass
class AF3RunnerConfig:
    """Configuration for structure prediction runner."""
    output_base_dir: str = "./af3_outputs"
    execution_mode: ExecutionMode = ExecutionMode.LOCAL

    # Backend selection
    backend: PredictionBackend = PredictionBackend.BOLTZ

    # AF3-specific config
    af3_binary: str = "alphafold3"
    model_dir: Optional[str] = None
    use_docker: bool = True
    docker_image: str = "alphafold3"
    docker_gpu_flag: str = "--gpus all"

    # Boltz-specific config
    boltz_binary: str = "boltz"
    use_msa: bool = True

    # Batch/SLURM config
    slurm_partition: str = "gpu"
    slurm_time: str = "4:00:00"
    slurm_gpus: int = 1
    slurm_mem: str = "64G"

    # Cloud/GCP config
    gcp_project: Optional[str] = None
    gcp_region: str = "us-central1"
    gcp_machine_type: str = "a2-highgpu-1g"


class AF3Runner:
    """
    AlphaFold3 execution manager.

    Handles input generation, job submission, and output collection
    for RNA-protein complex predictions.
    """

    def __init__(self, config: AF3RunnerConfig):
        self.config = config
        self.output_dir = Path(config.output_base_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self._job_cache: Dict[str, AF3Job] = {}

    def create_input(
        self,
        name: str,
        rna_sequence: str,
        protein_sequence: str,
        rna_chain_id: str = "A",
        protein_chain_id: str = "B"
    ) -> AF3Input:
        """Create an AF3 input specification."""
        return AF3Input(
            name=name,
            rna_sequence=rna_sequence,
            protein_sequence=protein_sequence,
            rna_chain_id=rna_chain_id,
            protein_chain_id=protein_chain_id
        )

    def _write_input_json(self, af3_input: AF3Input, job_dir: Path) -> Path:
        """Write AF3 input JSON file."""
        json_path = job_dir / f"{af3_input.name}.json"
        with open(json_path, 'w') as f:
            json.dump(af3_input.to_json_dict(), f, indent=2)
        return json_path

    def _get_job_dir(self, job_id: str) -> Path:
        """Get output directory for a job."""
        job_dir = self.output_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    def _check_cache(self, af3_input: AF3Input) -> Optional[Path]:
        """Check if result exists in cache."""
        input_hash = af3_input.get_hash()
        cache_dir = self.output_dir / "cache" / input_hash

        if cache_dir.exists():
            # Look for output files
            for suffix in ['_model.cif', '_confidences.json', '_summary.json']:
                results = list(cache_dir.glob(f"*{suffix}"))
                if results:
                    return cache_dir

        return None

    def submit_job(
        self,
        wt_input: AF3Input,
        mut_input: Optional[AF3Input] = None,
        job_id: Optional[str] = None,
        use_cache: bool = True
    ) -> AF3Job:
        """
        Submit an AF3 prediction job.

        Args:
            wt_input: Wild-type RNA-protein complex input
            mut_input: Optional mutant input (for paired analysis)
            job_id: Optional job identifier
            use_cache: Whether to check cache first

        Returns:
            AF3Job with status and paths
        """
        if job_id is None:
            job_id = f"{wt_input.name}_{wt_input.get_hash()}"

        job_dir = self._get_job_dir(job_id)

        # Check cache
        if use_cache:
            cached = self._check_cache(wt_input)
            if cached:
                print(f"Using cached result for {wt_input.name}", file=sys.stderr)
                return AF3Job(
                    job_id=job_id,
                    input_json_path=cached / f"{wt_input.name}.json",
                    output_dir=cached,
                    wt_input=wt_input,
                    mut_input=mut_input,
                    status="completed",
                    result_path=cached
                )

        # Write input JSON(s)
        wt_json = self._write_input_json(wt_input, job_dir)

        if mut_input:
            self._write_input_json(mut_input, job_dir)

        job = AF3Job(
            job_id=job_id,
            input_json_path=wt_json,
            output_dir=job_dir,
            wt_input=wt_input,
            mut_input=mut_input,
            status="pending"
        )

        # Execute based on mode
        if self.config.execution_mode == ExecutionMode.LOCAL:
            self._run_local(job)
        elif self.config.execution_mode == ExecutionMode.BATCH:
            self._generate_slurm_script(job)
        elif self.config.execution_mode == ExecutionMode.CLOUD:
            self._generate_cloud_config(job)

        self._job_cache[job_id] = job
        return job

    def _run_local(self, job: AF3Job):
        """Run prediction locally."""
        if self.config.backend == PredictionBackend.BOLTZ:
            self._run_boltz(job)
        elif self.config.use_docker:
            self._run_docker(job)
        else:
            self._run_native(job)

    def _run_docker(self, job: AF3Job):
        """Run AF3 via Docker container."""
        # Resolve absolute paths for volume mounts
        input_dir = job.input_json_path.parent.resolve()
        output_dir = job.output_dir.resolve()
        input_filename = job.input_json_path.name

        cmd = [
            "docker", "run", "--rm",
        ]

        # Add GPU flag
        if self.config.docker_gpu_flag:
            cmd.extend(self.config.docker_gpu_flag.split())

        # Volume mounts
        cmd.extend([
            "-v", f"{input_dir}:/input",
            "-v", f"{output_dir}:/output",
        ])

        # Model directory mount if specified
        if self.config.model_dir:
            model_dir = Path(self.config.model_dir).resolve()
            cmd.extend(["-v", f"{model_dir}:/models"])

        # Image and arguments
        cmd.append(self.config.docker_image)
        cmd.extend([
            f"--json_path=/input/{input_filename}",
            "--output_dir=/output"
        ])

        if self.config.model_dir:
            cmd.append("--model_dir=/models")

        print(f"Running AF3 via Docker: {' '.join(cmd)}", file=sys.stderr)

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )

            if result.returncode == 0:
                job.status = "completed"
                job.result_path = job.output_dir
            else:
                job.status = "failed"
                print(f"AF3 Docker failed: {result.stderr}", file=sys.stderr)

        except subprocess.TimeoutExpired:
            job.status = "timeout"
            print(f"AF3 timed out for {job.job_id}", file=sys.stderr)
        except FileNotFoundError:
            job.status = "failed"
            print("Docker not found. Install Docker or set use_docker=False", file=sys.stderr)

    def _run_native(self, job: AF3Job):
        """Run AF3 as native binary (non-Docker)."""
        cmd = [
            self.config.af3_binary,
            f"--json_path={job.input_json_path}",
            f"--output_dir={job.output_dir}"
        ]

        if self.config.model_dir:
            cmd.append(f"--model_dir={self.config.model_dir}")

        print(f"Running AF3 native: {' '.join(cmd)}", file=sys.stderr)

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )

            if result.returncode == 0:
                job.status = "completed"
                job.result_path = job.output_dir
            else:
                job.status = "failed"
                print(f"AF3 failed: {result.stderr}", file=sys.stderr)

        except subprocess.TimeoutExpired:
            job.status = "timeout"
            print(f"AF3 timed out for {job.job_id}", file=sys.stderr)
        except FileNotFoundError:
            job.status = "failed"
            print(f"AF3 binary not found: {self.config.af3_binary}", file=sys.stderr)

    def _write_boltz_input(self, af3_input: AF3Input, job_dir: Path) -> Path:
        """Write Boltz input FASTA file for RNA-protein complex."""
        fasta_path = job_dir / f"{af3_input.name}.fasta"

        with open(fasta_path, 'w') as f:
            # Boltz uses special headers for multi-chain complexes
            # Format: >chain_id|entity_type
            f.write(f">{af3_input.protein_chain_id}|protein\n")
            f.write(f"{af3_input.protein_sequence}\n")
            f.write(f">{af3_input.rna_chain_id}|rna\n")
            f.write(f"{af3_input.rna_sequence.replace('T', 'U')}\n")

        # Write MSA if provided
        if af3_input.protein_msa:
            msa_path = job_dir / f"{af3_input.name}_msa.a3m"
            with open(msa_path, 'w') as f:
                f.write(af3_input.protein_msa)

        return fasta_path

    def _run_boltz(self, job: AF3Job):
        """Run Boltz-1 for structure prediction."""
        # Write Boltz input format
        fasta_path = self._write_boltz_input(job.wt_input, job.output_dir)

        cmd = [
            self.config.boltz_binary,
            "predict",
            str(fasta_path),
            "--out_dir", str(job.output_dir),
            "--recycling_steps", "3",
            "--sampling_steps", "20"
        ]

        # Add MSA flag if available
        if job.wt_input.protein_msa and self.config.use_msa:
            msa_path = job.output_dir / f"{job.wt_input.name}_msa.a3m"
            cmd.extend(["--msa_file", str(msa_path)])

        print(f"Running Boltz: {' '.join(cmd)}", file=sys.stderr)

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )

            if result.returncode == 0:
                job.status = "completed"
                job.result_path = job.output_dir
                print(f"Boltz completed for {job.job_id}", file=sys.stderr)
            else:
                job.status = "failed"
                print(f"Boltz failed: {result.stderr}", file=sys.stderr)

        except subprocess.TimeoutExpired:
            job.status = "timeout"
            print(f"Boltz timed out for {job.job_id}", file=sys.stderr)
        except FileNotFoundError:
            job.status = "failed"
            print(f"Boltz not found. Install with: pip install boltz", file=sys.stderr)

    def _generate_slurm_script(self, job: AF3Job):
        """Generate SLURM submission script for batch execution."""
        script_path = job.output_dir / "submit.slurm"

        script_content = f"""#!/bin/bash
#SBATCH --job-name=af3_{job.job_id}
#SBATCH --partition={self.config.slurm_partition}
#SBATCH --time={self.config.slurm_time}
#SBATCH --gpus={self.config.slurm_gpus}
#SBATCH --mem={self.config.slurm_mem}
#SBATCH --output={job.output_dir}/slurm_%j.out
#SBATCH --error={job.output_dir}/slurm_%j.err

# Load required modules (adjust for your cluster)
# module load cuda/12.0
# module load alphafold3

cd {job.output_dir}

{self.config.af3_binary} \\
    --json_path={job.input_json_path} \\
    --output_dir={job.output_dir}
"""
        if self.config.model_dir:
            script_content = script_content.rstrip() + f" \\\n    --model_dir={self.config.model_dir}\n"

        with open(script_path, 'w') as f:
            f.write(script_content)

        job.status = "script_generated"
        print(f"SLURM script written to: {script_path}", file=sys.stderr)
        print(f"Submit with: sbatch {script_path}", file=sys.stderr)

    def _generate_cloud_config(self, job: AF3Job):
        """Generate GCP Batch configuration."""
        config_path = job.output_dir / "gcp_batch.json"

        gcp_config = {
            "taskGroups": [{
                "taskSpec": {
                    "runnables": [{
                        "container": {
                            "imageUri": "gcr.io/deepmind-alphafold/alphafold3:latest",
                            "commands": [
                                f"--json_path=/input/{job.input_json_path.name}",
                                f"--output_dir=/output"
                            ]
                        }
                    }],
                    "computeResource": {
                        "cpuMilli": 8000,
                        "memoryMib": 65536
                    },
                    "volumes": [{
                        "gcs": {"remotePath": f"gs://{self.config.gcp_project}/af3_jobs/{job.job_id}/input"},
                        "mountPath": "/input"
                    }, {
                        "gcs": {"remotePath": f"gs://{self.config.gcp_project}/af3_jobs/{job.job_id}/output"},
                        "mountPath": "/output"
                    }]
                },
                "taskCount": 1
            }],
            "allocationPolicy": {
                "instances": [{
                    "policy": {
                        "machineType": self.config.gcp_machine_type,
                        "accelerators": [{
                            "type": "nvidia-tesla-a100",
                            "count": 1
                        }]
                    }
                }],
                "location": {
                    "allowedLocations": [f"regions/{self.config.gcp_region}"]
                }
            },
            "logsPolicy": {
                "destination": "CLOUD_LOGGING"
            }
        }

        with open(config_path, 'w') as f:
            json.dump(gcp_config, f, indent=2)

        job.status = "config_generated"
        print(f"GCP Batch config written to: {config_path}", file=sys.stderr)
        print(f"Upload input and submit with: gcloud batch jobs submit {job.job_id} --config={config_path}", file=sys.stderr)

    def get_job_status(self, job_id: str) -> Optional[str]:
        """Get status of a submitted job."""
        if job_id in self._job_cache:
            return self._job_cache[job_id].status
        return None

    def collect_results(self, job: AF3Job) -> Optional[Dict[str, Path]]:
        """
        Collect output files from a completed job.

        Returns:
            Dict mapping output type to file path
        """
        if job.status != "completed" or not job.result_path:
            return None

        results = {}
        output_dir = job.result_path

        # Standard AF3 output files
        patterns = {
            'model': '*_model.cif',
            'confidences': '*_confidences.json',
            'summary': '*_summary_confidences.json',
            'ranking': '*ranking_scores.csv'
        }

        for name, pattern in patterns.items():
            matches = list(output_dir.glob(pattern))
            if matches:
                results[name] = matches[0]

        return results if results else None


def create_rna_protein_input(
    job_name: str,
    rna_seq: str,
    protein_seq: str,
    protein_msa: Optional[str] = None,
    mutation_pos: Optional[int] = None
) -> AF3Input:
    """
    Helper to create AF3 input for RNA-protein binding prediction.

    Args:
        job_name: Identifier for the job
        rna_seq: RNA sequence (will be converted T->U)
        protein_seq: Protein sequence
        protein_msa: Optional pre-computed MSA in A3M format
        mutation_pos: Optional position to mark (for naming)

    Returns:
        AF3Input ready for submission
    """
    return AF3Input(
        name=job_name,
        rna_sequence=rna_seq.upper().replace('T', 'U'),
        protein_sequence=protein_seq.upper(),
        rna_chain_id="R",
        protein_chain_id="P",
        protein_msa=protein_msa
    )


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Structure Prediction Runner')
    parser.add_argument('--mode', choices=['local', 'batch', 'cloud'], default='local')
    parser.add_argument('--backend', choices=['boltz', 'af3'], default='boltz',
                       help='Prediction backend (default: boltz)')
    parser.add_argument('--rna', required=True, help='RNA sequence')
    parser.add_argument('--protein', required=True, help='Protein sequence')
    parser.add_argument('--name', default='test', help='Job name')
    parser.add_argument('--output-dir', default='./af3_outputs')
    parser.add_argument('--af3-binary', default='alphafold3')
    parser.add_argument('--docker-image', default='alphafold3',
                       help='Docker image for AF3')
    parser.add_argument('--no-docker', action='store_true',
                       help='Run AF3 as native binary instead of Docker')
    parser.add_argument('--model-dir', help='Path to AF3 model weights')

    args = parser.parse_args()

    config = AF3RunnerConfig(
        execution_mode=ExecutionMode(args.mode),
        backend=PredictionBackend(args.backend),
        output_base_dir=args.output_dir,
        af3_binary=args.af3_binary,
        docker_image=args.docker_image,
        use_docker=not args.no_docker,
        model_dir=args.model_dir
    )

    runner = AF3Runner(config)

    af3_input = create_rna_protein_input(
        job_name=args.name,
        rna_seq=args.rna,
        protein_seq=args.protein
    )

    job = runner.submit_job(af3_input)
    print(f"Job {job.job_id}: {job.status}")
