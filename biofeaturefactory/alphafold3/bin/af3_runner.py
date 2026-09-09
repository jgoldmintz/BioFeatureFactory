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
AlphaFold3 execution wrapper.

Generates AF3 input JSON files and executes predictions in:
- local: Direct subprocess execution (requires GPU)

ExecutionMode.BATCH and ExecutionMode.CLOUD raise. Their script generators were
removed because both emitted artifacts that were never ingested. SLURM is handled
by biofeaturefactory.alphafold3.burst (manifest-driven array + L1 cache); cloud
burst (AWS Batch / GCP Batch) is out of scope.
"""

import json
import fcntl
import os
import queue
import re
import shutil
import stat
import string
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, Future
from pathlib import Path
from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Tuple, Iterator
from enum import Enum
import hashlib

from biofeaturefactory.alphafold3.bin.burst_manifest import is_result_complete


def detect_gpu_count() -> int:
    """Detect number of NVIDIA GPUs via nvidia-smi. Returns 1 if detection fails."""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            lines = [l.strip() for l in result.stdout.strip().splitlines() if l.strip()]
            return max(len(lines), 1)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return 1


class GPUPool:
    """Thread-safe pool of GPU device IDs."""

    def __init__(self, gpu_ids: List[int]):
        self._pool: queue.Queue = queue.Queue()
        for gid in gpu_ids:
            self._pool.put(gid)
        self.num_gpus = len(gpu_ids)

    def acquire(self) -> int:
        return self._pool.get()

    def release(self, gpu_id: int):
        self._pool.put(gpu_id)


class ExecutionMode(Enum):
    LOCAL = "local"
    BATCH = "batch"
    CLOUD = "cloud"


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
        # RNA chain (unpairedMsa required when --norun_data_pipeline)
        rna_entry = {
            "rna": {
                "id": self.rna_chain_id,
                "sequence": self.rna_sequence.replace('T', 'U'),
                "unpairedMsa": ""
            }
        }

        # Protein chain with MSA + templates (all required when --norun_data_pipeline)
        protein_entry = {
            "protein": {
                "id": self.protein_chain_id,
                "sequence": self.protein_sequence,
                "unpairedMsa": self.protein_msa if self.protein_msa else "",
                "pairedMsa": "",
                "templates": []
            }
        }

        return {
            "dialect": "alphafold3",
            "version": 2,
            "name": self.name,
            "modelSeeds": [1],
            "sequences": [rna_entry, protein_entry]
        }

    def get_hash(self) -> str:
        """Return the versioned content identity used by AF3 caches."""
        prediction_input = self.to_json_dict()
        prediction_input.pop("name")
        canonical_input = json.dumps(
            {
                "hash_schema_version": 1,
                "prediction_input": prediction_input,
            },
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        return hashlib.sha256(canonical_input.encode("utf-8")).hexdigest()

    def sanitised_name(self) -> str:
        """Return the filename-safe name used by AlphaFold3."""
        spaceless_name = self.name.replace(' ', '_')
        allowed_chars = set(string.ascii_letters + string.digits + '_-.')
        return ''.join(char for char in spaceless_name if char in allowed_chars)


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
class _DeferredJob:
    job: AF3Job
    future: Future


@dataclass
class AF3RunnerConfig:
    """Configuration for structure prediction runner.

    LOCAL/Docker is the only execution mode wired through this runner. SLURM
    and cloud bursts are handled by biofeaturefactory.alphafold3.burst, which
    has its own argparse for cluster-specific options.
    """
    output_base_dir: str = "./af3_outputs"
    execution_mode: ExecutionMode = ExecutionMode.LOCAL

    # AF3 config
    af3_binary: str = "alphafold3"
    model_dir: str = ""
    docker_image: str = "alphafold3"
    docker_gpu_flag: str = "--gpus all"

    # Parallelism
    max_gpus: Optional[int] = None  # None = auto-detect
    batch_size: int = 16

    # Resume and process-local acceleration
    resume: bool = True
    adopt_legacy_results: bool = False
    jax_cache_dir: Optional[str] = None
    timeout_per_job: int = 7200


class AF3Runner:
    """
    AlphaFold3 execution manager.

    Handles input generation, job submission, and output collection
    for RNA-protein complex predictions.
    """

    _PROVENANCE_SUFFIX = ".bff-af3-provenance.json"
    _PROVENANCE_SCHEMA_VERSION = 1
    _MODEL_FILE_PATTERNS = tuple(re.compile(pattern) for pattern in (
        r".*\.[0-9]+\.bin\.zst",
        r".*\.bin\.zst\.[0-9]+",
        r".*\.[0-9]+\.bin",
        r".*\.bin\.[0-9]+",
        r".*\.bin\.zst",
        r".*\.bin",
    ))

    def __init__(self, config: AF3RunnerConfig):
        if config.batch_size < 1:
            raise ValueError("batch_size must be at least 1")
        if config.timeout_per_job < 1:
            raise ValueError("timeout_per_job must be at least 1 second")
        if config.model_dir and not Path(config.model_dir).expanduser().is_dir():
            raise ValueError(
                "model_dir must be the directory containing AF3 model weights, "
                f"not a weights file: {config.model_dir}"
            )

        self.config = config
        self.output_dir = Path(config.output_base_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self._state_dir = self.output_dir.parent / ".af3_local"
        self._work_dir = self._state_dir / "work"
        self._log_dir = self._state_dir / "logs"
        self._lock_dir = self._state_dir / "locks"
        self._jax_cache_dir = (
            Path(config.jax_cache_dir).expanduser()
            if config.jax_cache_dir
            else self.output_dir.parent / ".cache" / "af3-jax"
        )
        self._model_dir = (
            Path(config.model_dir).expanduser() if config.model_dir else None
        )
        self._work_dir.mkdir(parents=True, exist_ok=True)
        self._log_dir.mkdir(parents=True, exist_ok=True)
        self._lock_dir.mkdir(parents=True, exist_ok=True)
        self._jax_cache_dir.mkdir(parents=True, exist_ok=True)

        self._job_cache: Dict[str, AF3Job] = {}
        self._deferred_jobs: Optional[List[_DeferredJob]] = None
        self._inflight_job_ids = set()
        self._job_lock_handles: Dict[str, Any] = {}
        self._job_lock = threading.Lock()
        self._docker_failed: bool = False
        self._consecutive_failures: int = 0
        self._last_error: Optional[str] = None
        self._failure_lock = threading.Lock()
        self._provenance_lock = threading.Lock()
        self._runtime_provenance: Optional[Dict[str, Any]] = None

        # GPU pool and thread pool for parallel execution
        num_gpus = config.max_gpus if config.max_gpus else detect_gpu_count()
        if num_gpus < 1:
            raise ValueError("max_gpus must be at least 1")
        self._num_gpus = num_gpus
        self._gpu_pool = GPUPool(list(range(num_gpus)))
        self._executor = ThreadPoolExecutor(max_workers=num_gpus)

        # Eager runtime preflight. _current_runtime_provenance inspects the
        # Docker image and hashes every weights file, and it raises RuntimeError
        # when the image, the weights or the docker binary is absent. Before
        # this it was reached only from the first job or a resume check, i.e.
        # AFTER POSTAR3 and the sequence mapper had loaded and once PER MUTATION,
        # each raise landing in the pipeline's per-mutation handler as another
        # FAILED row. Calling it here turns all three into one startup error.
        # The result is memoised on self._runtime_provenance, so the later call
        # sites pay nothing.
        self._preflight_runtime()

    def _preflight_runtime(self) -> None:
        """Fail fast when the container or the weights are not usable.

        Kept separate from __init__ so a caller that genuinely wants a runner
        without a live Docker daemon can subclass or patch it.
        """
        try:
            self._current_runtime_provenance()
        except RuntimeError as exc:
            raise RuntimeError(
                f"AF3 runtime preflight failed: {exc}\n"
                f"  docker image : {self.config.docker_image}\n"
                f"  model_dir    : {self.config.model_dir}\n"
                "Nothing has been submitted. Fix the above and re-run."
            ) from exc

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
        """Atomically write an AF3 input JSON file."""
        json_path = job_dir / f"{af3_input.name}.json"
        tmp_path = json_path.with_name(
            f"{json_path.name}.{uuid.uuid4().hex}.tmp"
        )
        try:
            with open(tmp_path, 'w') as output_handle:
                json.dump(af3_input.to_json_dict(), output_handle, indent=2)
                output_handle.flush()
                os.fsync(output_handle.fileno())
            os.replace(tmp_path, json_path)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()
        return json_path

    def _get_job_dir(self, job_id: str) -> Path:
        """Get output directory for a job."""
        job_dir = self.output_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    @staticmethod
    def _sha256_file(path: Path) -> str:
        """Hash one model file without loading it into memory."""
        digest = hashlib.sha256()
        with open(path, "rb") as input_handle:
            while chunk := input_handle.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()

    def _model_files(self) -> List[Path]:
        """Return direct model files matching AF3's supported naming forms."""
        if self._model_dir is None:
            raise RuntimeError("model_dir is required for local execution")
        model_files = sorted(
            path for path in self._model_dir.iterdir()
            if path.is_file() and any(
                pattern.fullmatch(path.name)
                for pattern in self._MODEL_FILE_PATTERNS
            )
        )
        if not model_files:
            raise RuntimeError(
                f"No AF3 .bin or .bin.zst model files found in {self._model_dir}"
            )
        return model_files

    def _current_runtime_provenance(self) -> Dict[str, Any]:
        """Fingerprint the container and weights used for local inference."""
        with self._provenance_lock:
            if self._runtime_provenance is not None:
                return self._runtime_provenance

            try:
                image_result = subprocess.run(
                    [
                        "docker", "image", "inspect",
                        "--format={{.Id}}", self.config.docker_image,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
            except (OSError, subprocess.TimeoutExpired) as error:
                raise RuntimeError(
                    f"Unable to inspect Docker image {self.config.docker_image}: "
                    f"{error}"
                ) from error
            image_id = image_result.stdout.strip()
            if image_result.returncode != 0 or not image_id:
                diagnostic = image_result.stderr.strip() or "image not found"
                raise RuntimeError(
                    f"Unable to inspect Docker image {self.config.docker_image}: "
                    f"{diagnostic}"
                )

            model_files = [
                {
                    "name": path.name,
                    "size": path.stat().st_size,
                    "sha256": self._sha256_file(path),
                }
                for path in self._model_files()
            ]
            self._runtime_provenance = {
                "schema_version": self._PROVENANCE_SCHEMA_VERSION,
                "docker_image_id": image_id,
                "model_files": model_files,
                "inference_mode": "norun_data_pipeline",
            }
            return self._runtime_provenance

    def _expected_job_provenance(self, af3_input: AF3Input) -> Dict[str, Any]:
        """Combine runtime identity with the exact named AF3 payload."""
        canonical_input = json.dumps(
            af3_input.to_json_dict(),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        return {
            "runtime": self._current_runtime_provenance(),
            "input_sha256": hashlib.sha256(
                canonical_input.encode()
            ).hexdigest(),
        }

    def _result_provenance_path(self, result_dir: Path) -> Path:
        """Place metadata beside Docker-owned output directories."""
        return result_dir.parent / f".{result_dir.name}{self._PROVENANCE_SUFFIX}"

    def _read_result_provenance(
        self,
        result_dir: Path
    ) -> Optional[Dict[str, Any]]:
        """Read a result provenance sidecar, rejecting malformed metadata."""
        provenance_path = self._result_provenance_path(result_dir)
        try:
            with open(provenance_path, "r") as provenance_handle:
                provenance = json.load(provenance_handle)
        except (OSError, json.JSONDecodeError):
            return None
        return provenance if isinstance(provenance, dict) else None

    def _write_result_provenance(
        self,
        result_dir: Path,
        af3_input: AF3Input
    ) -> None:
        """Atomically persist the scientific identity beside one result."""
        provenance_path = self._result_provenance_path(result_dir)
        temporary_path = provenance_path.with_name(
            f"{provenance_path.name}.{uuid.uuid4().hex}.tmp"
        )
        try:
            with open(temporary_path, "w") as provenance_handle:
                json.dump(
                    self._expected_job_provenance(af3_input),
                    provenance_handle,
                    indent=2,
                    sort_keys=True,
                )
                provenance_handle.flush()
                os.fsync(provenance_handle.fileno())
            os.replace(temporary_path, provenance_path)
        finally:
            if temporary_path.exists():
                temporary_path.unlink()

    def _provenance_matches(
        self,
        result_dir: Path,
        af3_input: AF3Input
    ) -> bool:
        """Return whether saved result provenance matches this exact run."""
        return self._read_result_provenance(result_dir) == (
            self._expected_job_provenance(af3_input)
        )

    def _reserve_job_id(self, job_id: str) -> None:
        """Lock one mutable job path across threads and runner processes."""
        with self._job_lock:
            if job_id in self._inflight_job_ids:
                raise RuntimeError(f"AF3 job_id is already in flight: {job_id}")
            lock_name = hashlib.sha256(job_id.encode()).hexdigest()
            lock_path = self._lock_dir / f"{lock_name}.lock"
            lock_handle = open(lock_path, "a+")
            try:
                fcntl.flock(
                    lock_handle.fileno(),
                    fcntl.LOCK_EX | fcntl.LOCK_NB,
                )
            except BlockingIOError:
                lock_handle.close()
                raise RuntimeError(
                    f"AF3 job_id is already running in another process: {job_id}"
                ) from None
            except BaseException:
                lock_handle.close()
                raise
            lock_handle.seek(0)
            lock_handle.truncate()
            lock_handle.write(f"pid={os.getpid()} job_id={job_id}\n")
            lock_handle.flush()
            self._inflight_job_ids.add(job_id)
            self._job_lock_handles[job_id] = lock_handle

    def _release_job_id(self, job_id: str) -> None:
        """Release a completed or aborted job path reservation."""
        with self._job_lock:
            self._inflight_job_ids.discard(job_id)
            lock_handle = self._job_lock_handles.pop(job_id, None)
            if lock_handle is not None:
                try:
                    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
                finally:
                    lock_handle.close()

    @staticmethod
    def _input_json_matches(json_path: Path, af3_input: AF3Input) -> bool:
        """Return whether a saved input is exactly the requested AF3 input."""
        try:
            with open(json_path, 'r') as input_handle:
                saved_input = json.load(input_handle)
        except (OSError, json.JSONDecodeError):
            return False
        return saved_input == af3_input.to_json_dict()

    @staticmethod
    def _is_complete_result(result_dir: Path, output_name: str) -> bool:
        """Require AF3's three nonempty top-ranked files for one exact job."""
        return is_result_complete(result_dir, output_name)

    def _find_resumable_result(
        self,
        af3_input: AF3Input,
        job_dir: Path
    ) -> Optional[Path]:
        """Find a complete result whose saved input exactly matches this job."""
        input_json_path = job_dir / f"{af3_input.name}.json"
        if not self._input_json_matches(input_json_path, af3_input):
            return None

        output_root = job_dir / "output"
        output_name = af3_input.sanitised_name()
        exact_result = output_root / output_name
        candidates = [exact_result]
        if output_root.is_dir():
            timestamped = sorted(
                (
                    path for path in output_root.glob(f"{output_name}_*")
                    if path.is_dir()
                ),
                key=lambda path: path.stat().st_mtime,
                reverse=True
            )
            candidates.extend(timestamped)

        for candidate in candidates:
            if self._is_complete_result(candidate, output_name):
                return candidate
        return None

    @staticmethod
    def _quarantine_path(path: Path) -> Optional[Path]:
        """Move stale or partial output aside without deleting it."""
        if not path.exists():
            return None
        suffix = f"stale-{time.strftime('%Y%m%d-%H%M%S')}-{uuid.uuid4().hex[:8]}"
        quarantined = path.with_name(f"{path.name}.{suffix}")
        os.replace(path, quarantined)
        return quarantined

    def _quarantine_result(self, result_dir: Path) -> Optional[Path]:
        """Preserve a result together with its candidate-specific provenance."""
        provenance_path = self._result_provenance_path(result_dir)
        quarantined = self._quarantine_path(result_dir)
        if quarantined is None:
            return None
        if provenance_path.exists():
            os.replace(
                provenance_path,
                self._result_provenance_path(quarantined),
            )
        return quarantined

    def _prepare_job(
        self,
        af3_input: AF3Input,
        mut_input: Optional[AF3Input],
        job_id: str,
        use_cache: bool
    ) -> AF3Job:
        """Create one job, adopting a matching complete result when available."""
        self._reserve_job_id(job_id)
        try:
            job_dir = self._get_job_dir(job_id)
            input_json_path = job_dir / f"{af3_input.name}.json"
            input_matches = self._input_json_matches(input_json_path, af3_input)
            should_resume = use_cache and self.config.resume
            rejected_provenance = False

            if should_resume and input_matches:
                resumed_result = self._find_resumable_result(af3_input, job_dir)
                if resumed_result is not None:
                    provenance_path = self._result_provenance_path(
                        resumed_result
                    )
                    if (not provenance_path.exists()
                            and not self.config.adopt_legacy_results):
                        raise RuntimeError(
                            f"Complete legacy AF3 output for {af3_input.name} "
                            "has no model/image provenance. Verify that it was "
                            "created with the current weights and Docker image, "
                            "then rerun once with --adopt-legacy-results."
                        )
                    saved_provenance = self._read_result_provenance(
                        resumed_result
                    )
                    expected_provenance = self._expected_job_provenance(af3_input)
                    if saved_provenance == expected_provenance:
                        pass
                    elif (not provenance_path.exists()
                          and self.config.adopt_legacy_results):
                        self._write_result_provenance(
                            resumed_result, af3_input
                        )
                        print(
                            f"Adopted legacy AF3 result for {af3_input.name}",
                            file=sys.stderr
                        )
                    else:
                        rejected_provenance = True

                    if not rejected_provenance:
                        print(
                            f"Resuming completed AF3 result for {af3_input.name}",
                            file=sys.stderr
                        )
                        job = AF3Job(
                            job_id=job_id,
                            input_json_path=input_json_path,
                            output_dir=job_dir,
                            wt_input=af3_input,
                            mut_input=mut_input,
                            status="completed",
                            result_path=resumed_result
                        )
                        self._job_cache[job_id] = job
                        self._release_job_id(job_id)
                        return job

            if (job_dir / "output").exists() and (
                not input_matches or not should_resume or rejected_provenance
            ):
                quarantined = self._quarantine_path(job_dir / "output")
                if quarantined is not None:
                    print(
                        f"Preserved non-resumable output for {af3_input.name} at "
                        f"{quarantined}",
                        file=sys.stderr
                    )

            written_json = self._write_input_json(af3_input, job_dir)
            if mut_input:
                self._write_input_json(mut_input, job_dir)

            job = AF3Job(
                job_id=job_id,
                input_json_path=written_json,
                output_dir=job_dir,
                wt_input=af3_input,
                mut_input=mut_input,
                status="pending"
            )
            self._job_cache[job_id] = job
            return job
        except BaseException:
            self._release_job_id(job_id)
            raise

    @staticmethod
    def _resolved_future(job: AF3Job) -> Future:
        """Return an already-resolved Future for a resumed job."""
        future: Future = Future()
        future.set_result(job)
        return future

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

        job = self._prepare_job(wt_input, mut_input, job_id, use_cache)
        if job.status == "completed":
            return job

        # Execute based on mode
        if self.config.execution_mode == ExecutionMode.LOCAL:
            self._run_local(job)
        elif self.config.execution_mode == ExecutionMode.BATCH:
            self._release_job_id(job.job_id)
            raise RuntimeError(
                "ExecutionMode.BATCH is no longer supported in af3_runner. "
                "The prior stub generated a broken script and never ingested "
                "results. Use the burst driver instead: "
                "`python -m biofeaturefactory.alphafold3.burst submit ...`"
            )
        elif self.config.execution_mode == ExecutionMode.CLOUD:
            self._release_job_id(job.job_id)
            raise RuntimeError(
                "ExecutionMode.CLOUD is not implemented; the prior stub "
                "never ingested results. A cloud-burst driver (AWS Batch / "
                "GCP Batch) is out of scope."
            )
        else:
            self._release_job_id(job.job_id)
            raise ValueError(
                f"Unknown execution mode: {self.config.execution_mode}"
            )

        return job

    def submit_job_async(
        self,
        wt_input: AF3Input,
        mut_input: Optional[AF3Input] = None,
        job_id: Optional[str] = None,
        use_cache: bool = True
    ) -> Future:
        """
        Non-blocking job submission. Returns a Future[AF3Job].
        Cache hits resolve immediately without acquiring a GPU.
        """
        if job_id is None:
            job_id = f"{wt_input.name}_{wt_input.get_hash()}"

        job = self._prepare_job(wt_input, mut_input, job_id, use_cache)
        if job.status == "completed":
            return self._resolved_future(job)

        if self.config.execution_mode == ExecutionMode.LOCAL:
            if self._deferred_jobs is not None:
                future: Future = Future()
                self._deferred_jobs.append(_DeferredJob(job=job, future=future))
                return future
            try:
                worker_future = self._executor.submit(
                    self._run_local_with_gpu, job
                )
            except BaseException:
                self._release_job_id(job.job_id)
                raise
            worker_future.add_done_callback(
                lambda completed, submitted_job=job:
                self._release_cancelled_job(completed, submitted_job)
            )
            return worker_future
        elif self.config.execution_mode == ExecutionMode.BATCH:
            self._release_job_id(job.job_id)
            raise RuntimeError(
                "ExecutionMode.BATCH is no longer supported in af3_runner. "
                "The prior stub generated a broken script and never ingested "
                "results. Use the burst driver instead: "
                "`python -m biofeaturefactory.alphafold3.burst submit ...`"
            )
        elif self.config.execution_mode == ExecutionMode.CLOUD:
            self._release_job_id(job.job_id)
            raise RuntimeError(
                "ExecutionMode.CLOUD is not implemented; the prior stub "
                "never ingested results. A cloud-burst driver (AWS Batch / "
                "GCP Batch) is out of scope."
            )
        else:
            self._release_job_id(job.job_id)
            raise ValueError(f"Unknown execution mode: {self.config.execution_mode}")

    @contextmanager
    def batch_submissions(self) -> Iterator[None]:
        """Defer asynchronous submissions and dispatch bounded AF3 batches."""
        if self.config.batch_size == 1:
            yield
            return
        if self._deferred_jobs is not None:
            raise RuntimeError("Nested AF3 batch submission scopes are not supported")

        deferred_jobs: List[_DeferredJob] = []
        self._deferred_jobs = deferred_jobs
        try:
            yield
        except BaseException:
            for deferred in deferred_jobs:
                deferred.future.cancel()
                self._release_job_id(deferred.job.job_id)
            raise
        else:
            try:
                self._dispatch_deferred_jobs(deferred_jobs)
            except BaseException:
                for deferred in deferred_jobs:
                    if not deferred.future.done():
                        deferred.future.cancel()
                    self._release_job_id(deferred.job.job_id)
                raise
        finally:
            self._deferred_jobs = None

    def _dispatch_deferred_jobs(self, deferred_jobs: List[_DeferredJob]) -> None:
        """Submit bounded, balanced groups while keeping every GPU eligible."""
        runnable_jobs = []
        for deferred in deferred_jobs:
            if deferred.future.set_running_or_notify_cancel():
                runnable_jobs.append(deferred)
            else:
                deferred.job.status = "cancelled"
                self._job_cache[deferred.job.job_id] = deferred.job
                self._release_job_id(deferred.job.job_id)

        job_count = len(runnable_jobs)
        if job_count == 0:
            return
        size_limited_batches = (
            job_count + self.config.batch_size - 1
        ) // self.config.batch_size
        batch_count = max(
            size_limited_batches,
            min(self._num_gpus, job_count),
        )
        base_size, extra_jobs = divmod(job_count, batch_count)
        start = 0
        for batch_index in range(batch_count):
            batch_length = base_size + (batch_index < extra_jobs)
            batch = runnable_jobs[start:start + batch_length]
            start += batch_length
            batch_future = self._executor.submit(
                self._run_local_batch_with_gpu,
                [deferred.job for deferred in batch]
            )
            batch_future.add_done_callback(
                lambda completed, items=batch: self._resolve_batch_futures(
                    completed, items
                )
            )

    def _release_cancelled_job(
        self,
        completed_future: Future,
        job: AF3Job
    ) -> None:
        """Release a queued executor job when cancellation prevents startup."""
        if not completed_future.cancelled():
            return
        job.status = "cancelled"
        self._job_cache[job.job_id] = job
        self._release_job_id(job.job_id)

    @staticmethod
    def _resolve_batch_futures(
        batch_future: Future,
        deferred_jobs: List[_DeferredJob]
    ) -> None:
        """Propagate one worker result to every per-job Future."""
        try:
            completed_jobs = batch_future.result()
        except BaseException as error:
            for deferred in deferred_jobs:
                if not deferred.future.done():
                    deferred.future.set_exception(error)
            return

        for deferred, job in zip(deferred_jobs, completed_jobs):
            if not deferred.future.done():
                deferred.future.set_result(job)

    def _run_local_batch_with_gpu(self, jobs: List[AF3Job]) -> List[AF3Job]:
        """Acquire one GPU for a bounded batch and resolve every job status."""
        try:
            with self._failure_lock:
                docker_failed = self._docker_failed
            if docker_failed:
                for job in jobs:
                    job.status = "failed"
                return jobs

            if not self.config.model_dir:
                with self._failure_lock:
                    if not self._docker_failed:
                        self._docker_failed = True
                        print(
                            "Error: --model-dir is required for local execution.",
                            file=sys.stderr
                        )
                for job in jobs:
                    job.status = "failed"
                return jobs

            gpu_id = self._gpu_pool.acquire()
            try:
                self._run_docker_batch(jobs, gpu_id=gpu_id)
            finally:
                self._gpu_pool.release(gpu_id)
            return jobs
        finally:
            for job in jobs:
                self._job_cache[job.job_id] = job
                self._release_job_id(job.job_id)

    def _run_local_with_gpu(self, job: AF3Job) -> AF3Job:
        """Acquire a GPU, run Docker pinned to it, release. Called from thread pool."""
        try:
            with self._failure_lock:
                docker_failed = self._docker_failed
            if docker_failed:
                job.status = "failed"
                return job

            if not self.config.model_dir:
                job.status = "failed"
                with self._failure_lock:
                    if not self._docker_failed:
                        self._docker_failed = True
                        print(
                            "Error: --model-dir is required for local execution.",
                            file=sys.stderr
                        )
                return job

            gpu_id = self._gpu_pool.acquire()
            try:
                self._run_docker(job, gpu_id=gpu_id)
            finally:
                self._gpu_pool.release(gpu_id)
            return job
        finally:
            self._job_cache[job.job_id] = job
            self._release_job_id(job.job_id)

    def _run_local(self, job: AF3Job):
        """Run prediction locally via Docker (synchronous, backward-compatible)."""
        try:
            if not self.config.model_dir:
                job.status = "failed"
                with self._failure_lock:
                    if not self._docker_failed:
                        self._docker_failed = True
                        print(
                            "Error: --model-dir is required for local execution. "
                            "Point it to the directory containing AF3 model weights.",
                            file=sys.stderr
                        )
                return
            gpu_id = self._gpu_pool.acquire()
            try:
                self._run_docker(job, gpu_id=gpu_id)
            finally:
                self._gpu_pool.release(gpu_id)
        finally:
            self._job_cache[job.job_id] = job
            self._release_job_id(job.job_id)

    def _build_docker_batch_command(
        self,
        input_dir: Path,
        output_dir: Path,
        cidfile: Path,
        container_name: str,
        gpu_id: Optional[int]
    ) -> List[str]:
        """Build one native AF3 directory-mode Docker invocation."""
        if self._model_dir is None:
            raise ValueError("model_dir is required for local execution")
        model_dir = self._model_dir.resolve()
        docker_image_id = self._current_runtime_provenance()["docker_image_id"]
        command = [
            "docker", "run", "--rm",
            "--name", container_name,
            "--cidfile", str(cidfile),
        ]
        if gpu_id is not None:
            command.extend(["--gpus", f"device={gpu_id}"])
        elif self.config.docker_gpu_flag:
            command.extend(self.config.docker_gpu_flag.split())

        command.extend([
            "-v", f"{input_dir.resolve()}:/root/af_input:ro",
            "-v", f"{output_dir.resolve()}:/root/af_output",
            "-v", f"{model_dir}:/root/models:ro",
            "-v", f"{self._jax_cache_dir.resolve()}:/root/jax_cache",
            docker_image_id,
            "python", "run_alphafold.py",
            "--input_dir=/root/af_input",
            "--model_dir=/root/models",
            "--output_dir=/root/af_output",
            "--jax_compilation_cache_dir=/root/jax_cache",
            "--norun_data_pipeline",
        ])
        return command

    @staticmethod
    def _read_log_tail(log_path: Path, max_bytes: int = 12000) -> str:
        """Read only the bounded tail needed for diagnostics."""
        try:
            with open(log_path, 'rb') as log_handle:
                log_handle.seek(0, os.SEEK_END)
                size = log_handle.tell()
                log_handle.seek(max(0, size - max_bytes))
                return log_handle.read().decode(errors='replace').strip()
        except OSError:
            return ""

    @staticmethod
    def _is_infrastructure_error(error_text: str) -> bool:
        """Return whether retrying biological inputs cannot fix this error."""
        lowered = error_text.lower()
        return any(message in lowered for message in (
            "cannot connect to the docker daemon",
            "error response from daemon",
            "pull access denied",
            "unable to find image",
            "failed to create shim task",
        ))

    @staticmethod
    def _stop_container(cidfile: Path, container_name: str) -> bool:
        """Stop a timed-out container and confirm it is no longer present."""
        try:
            container_id = cidfile.read_text().strip()
        except OSError:
            container_id = ""
        target = container_id or container_name
        try:
            remove_result = subprocess.run(
                ["docker", "rm", "-f", target],
                capture_output=True,
                text=True,
                timeout=30
            )
            if remove_result.returncode == 0:
                return True
            inspect_result = subprocess.run(
                ["docker", "inspect", target],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if inspect_result.returncode == 0:
                return False
            diagnostic = (
                f"{inspect_result.stdout}\n{inspect_result.stderr}"
            ).lower()
            return (
                "no such object" in diagnostic
                or "no such container" in diagnostic
            )
        except (OSError, subprocess.TimeoutExpired):
            return False

    def _publish_completed_outputs(
        self,
        jobs: List[AF3Job],
        staged_output_dir: Path
    ) -> List[AF3Job]:
        """Atomically publish independently complete children from a batch."""
        completed_jobs = []
        for job in jobs:
            output_name = job.wt_input.sanitised_name()
            staged_result = staged_output_dir / output_name
            if not self._is_complete_result(staged_result, output_name):
                continue
            canonical_root = job.output_dir / "output"
            canonical_root.mkdir(parents=True, exist_ok=True)
            canonical_result = canonical_root / output_name
            if canonical_result.exists():
                if (self._is_complete_result(canonical_result, output_name)
                        and self._provenance_matches(
                            canonical_result, job.wt_input
                        )):
                    job.status = "completed"
                    job.result_path = canonical_result
                    completed_jobs.append(job)
                    continue
                self._quarantine_result(canonical_result)

            os.replace(staged_result, canonical_result)
            self._write_result_provenance(canonical_result, job.wt_input)
            job.status = "completed"
            job.result_path = canonical_result
            completed_jobs.append(job)
        return completed_jobs

    def _record_single_failure(
        self,
        job: AF3Job,
        error_text: str,
        timed_out: bool = False
    ) -> None:
        """Apply the existing per-job failure breaker after fallback isolation."""
        if timed_out:
            job.status = "timeout"
            print(f"AF3 timed out for {job.job_id}", file=sys.stderr)
            return

        job.status = "failed"
        lines = error_text.splitlines()
        error_key = lines[-1] if lines else error_text
        with self._failure_lock:
            if error_key == self._last_error:
                self._consecutive_failures += 1
            else:
                self._consecutive_failures = 1
                self._last_error = error_key
            if self._consecutive_failures >= 3:
                self._docker_failed = True
                print(
                    "AF3 failed 3 consecutive times with the same isolated "
                    f"error, aborting: {error_text[-300:]}",
                    file=sys.stderr
                )
            else:
                print(f"AF3 failed for {job.job_id}: {error_text}",
                      file=sys.stderr)

    def _preserve_staged_outputs(
        self,
        output_dir: Path,
        jobs: List[AF3Job],
        batch_id: str
    ) -> None:
        """Rename Docker-owned incomplete trees out of temporary storage."""
        if not output_dir.is_dir():
            return
        jobs_by_name = {
            job.wt_input.sanitised_name(): job for job in jobs
        }
        fallback_root = self._state_dir / "incomplete" / batch_id
        for staged_path in list(output_dir.iterdir()):
            if staged_path.is_dir() and not any(staged_path.iterdir()):
                staged_path.rmdir()
                continue
            job = jobs_by_name.get(staged_path.name)
            if job is None:
                destination_root = fallback_root
            else:
                destination_root = job.output_dir / "output"
            destination_root.mkdir(parents=True, exist_ok=True)
            destination = destination_root / (
                f"{staged_path.name}.incomplete-{batch_id}"
            )
            if destination.exists():
                destination = destination.with_name(
                    f"{destination.name}-{uuid.uuid4().hex[:8]}"
                )
            try:
                os.replace(staged_path, destination)
            except PermissionError:
                original_mode = stat.S_IMODE(staged_path.stat().st_mode)
                staged_path.chmod(original_mode | stat.S_IWUSR)
                os.replace(staged_path, destination)
                destination.chmod(original_mode)
            print(
                f"Preserved incomplete AF3 output at {destination}",
                file=sys.stderr,
            )

    def _run_docker_batch(
        self,
        jobs: List[AF3Job],
        gpu_id: Optional[int] = None
    ) -> None:
        """Run a bounded AF3 input directory and retry only missing children."""
        if not jobs:
            return
        with self._failure_lock:
            docker_failed = self._docker_failed
        if docker_failed:
            for job in jobs:
                job.status = "failed"
            return

        output_names = [job.wt_input.sanitised_name() for job in jobs]
        if any(not name for name in output_names):
            raise ValueError("AF3 job names must contain filename-safe characters")
        if len(set(output_names)) != len(output_names):
            raise ValueError("AF3 batch contains colliding sanitised job names")

        batch_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{uuid.uuid4().hex[:8]}"
        log_path = self._log_dir / f"batch-{batch_id}.log"
        stderr_path = self._log_dir / f"batch-{batch_id}.stderr.log"
        return_code: Optional[int] = None
        timed_out = False
        invocation_error = ""

        temporary_path = Path(tempfile.mkdtemp(
            prefix=f"batch-{batch_id}-",
            dir=self._work_dir
        ))
        cleanup_safe = True
        output_dir = temporary_path / "outputs"
        try:
            input_dir = temporary_path / "inputs"
            input_dir.mkdir()
            output_dir.mkdir()
            for output_name in output_names:
                (output_dir / output_name).mkdir()
            cidfile = temporary_path / "container.cid"
            container_name = f"bff-af3-{uuid.uuid4().hex[:20]}"

            for index, job in enumerate(jobs):
                staged_json = input_dir / f"{index:04d}.json"
                with open(staged_json, "w") as staged_handle:
                    json.dump(job.wt_input.to_json_dict(), staged_handle)

            command = self._build_docker_batch_command(
                input_dir, output_dir, cidfile, container_name, gpu_id
            )
            print(
                f"Running AF3 batch ({len(jobs)} inputs) [GPU {gpu_id}]: "
                f"{' '.join(command)}",
                file=sys.stderr
            )

            try:
                with open(log_path, 'w') as log_handle, open(stderr_path, 'w') as stderr_handle:
                    result = subprocess.run(
                        command,
                        stdout=log_handle,
                        stderr=stderr_handle,
                        text=True,
                        timeout=self.config.timeout_per_job * len(jobs)
                    )
                return_code = result.returncode
            except subprocess.TimeoutExpired:
                timed_out = True
                if not self._stop_container(cidfile, container_name):
                    cleanup_safe = False
                    invocation_error = (
                        "AF3 timed out and container cleanup could not be "
                        f"confirmed for {container_name}"
                    )
            except OSError as error:
                invocation_error = str(error)

            log_tail = self._read_log_tail(log_path)
            stderr_tail = self._read_log_tail(stderr_path)
            log_tail = '\n'.join(part for part in (log_tail, stderr_tail) if part)
            if invocation_error:
                log_tail = invocation_error
            completed_jobs = (
                [] if not cleanup_safe or invocation_error else
                self._publish_completed_outputs(jobs, output_dir)
            )
            completed_ids = {job.job_id for job in completed_jobs}
            incomplete_jobs = [
                job for job in jobs if job.job_id not in completed_ids
            ]
        finally:
            if cleanup_safe:
                try:
                    self._preserve_staged_outputs(output_dir, jobs, batch_id)
                    shutil.rmtree(temporary_path)
                except OSError as error:
                    cleanup_safe = False
                    invocation_error = (
                        "Could not safely preserve or remove Docker-owned AF3 "
                        f"temporary output at {temporary_path}: {error}"
                    )
            if not cleanup_safe:
                print(
                    f"Left AF3 work directory intact at {temporary_path}",
                    file=sys.stderr,
                )

        if invocation_error:
            log_tail = invocation_error

        if not incomplete_jobs:
            with self._failure_lock:
                self._consecutive_failures = 0
                self._last_error = None
            return

        if invocation_error or (return_code in (125, 126, 127)
                                and self._is_infrastructure_error(stderr_tail)):
            with self._failure_lock:
                self._docker_failed = True
            for job in incomplete_jobs:
                job.status = "failed"
            diagnostic = log_tail or "Docker invocation failed"
            print(
                f"Docker error (aborting remaining jobs): {diagnostic}",
                file=sys.stderr
            )
            return

        if len(jobs) > 1:
            print(
                f"AF3 batch left {len(incomplete_jobs)} incomplete input(s); "
                "retrying those inputs individually",
                file=sys.stderr
            )
            for job in incomplete_jobs:
                self._run_docker_batch([job], gpu_id=gpu_id)
            return

        if return_code == 0:
            diagnostic = "AF3 exited successfully without complete primary outputs"
            if log_tail:
                diagnostic += f"; see {log_path}\n{log_tail}"
            log_tail = diagnostic
        elif not log_tail:
            log_tail = f"AF3 exited with status {return_code}"
        self._record_single_failure(
            incomplete_jobs[0],
            log_tail,
            timed_out=timed_out
        )

    def _run_docker(self, job: AF3Job, gpu_id: Optional[int] = None):
        """Run one AF3 job through the shared directory-mode implementation."""
        self._run_docker_batch([job], gpu_id=gpu_id)

    # SLURM/BATCH and GCP/CLOUD script generators were removed; both produced
    # broken artifacts and had no result-ingestion path. SLURM is now handled
    # by biofeaturefactory.alphafold3.burst (separate submit + ingest driver
    # with a manifest-driven SLURM array). Cloud-burst (AWS Batch / GCP Batch)
    # is out of scope.

    def shutdown(self):
        """Shut down the thread pool. Call after all jobs complete."""
        self._executor.shutdown(wait=True)

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

        # Standard AF3 output files (may be in a timestamped subdirectory)
        patterns = {
            'model': '**/*_model.cif',
            'confidences': '**/*_confidences.json',
            'summary': '**/*_summary_confidences.json',
            'ranking': '**/*ranking_scores.csv'
        }

        for name, pattern in patterns.items():
            matches = list(output_dir.glob(pattern))
            # Prefer top-level ranked files over per-sample files
            top_level = [f for f in matches if 'seed-' not in f.name]
            if top_level:
                results[name] = top_level[0]
            elif matches:
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
    parser.add_argument('--rna', required=True, help='RNA sequence')
    parser.add_argument('--protein', required=True, help='Protein sequence')
    parser.add_argument('--name', default='test', help='Job name')
    parser.add_argument('--output-dir', default='./af3_outputs')
    parser.add_argument('--af3-binary', default='alphafold3')
    parser.add_argument('--docker-image', default='alphafold3',
                       help='Docker image for AF3')
    parser.add_argument('--model-dir', help='Path to AF3 model weights')
    parser.add_argument('--no-resume', action='store_false', dest='resume')
    parser.add_argument('--adopt-legacy-results', action='store_true')
    parser.add_argument('--jax-cache-dir')

    args = parser.parse_args()

    config = AF3RunnerConfig(
        execution_mode=ExecutionMode(args.mode),
        output_base_dir=args.output_dir,
        af3_binary=args.af3_binary,
        docker_image=args.docker_image,
        model_dir=args.model_dir,
        resume=args.resume,
        adopt_legacy_results=args.adopt_legacy_results,
        jax_cache_dir=args.jax_cache_dir
    )

    runner = AF3Runner(config)

    af3_input = create_rna_protein_input(
        job_name=args.name,
        rna_seq=args.rna,
        protein_seq=args.protein
    )

    job = runner.submit_job(af3_input)
    print(f"Job {job.job_id}: {job.status}")
