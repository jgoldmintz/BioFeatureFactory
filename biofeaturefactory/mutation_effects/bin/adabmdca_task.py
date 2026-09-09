"""Run one planned adabmDCA task and publish verified, complete artifacts."""

from __future__ import annotations

import argparse
import errno
import fcntl
import hashlib
import json
import math
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager, nullcontext
from pathlib import Path


CUDA_OOM = re.compile(
    r"(?:torch\.cuda\.OutOfMemoryError\s*:|"
    r"(?:RuntimeError|torch\.OutOfMemoryError)\s*:\s*"
    r"CUDA(?: out of memory| error:\s*out of memory))",
    re.IGNORECASE,
)
THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "BFF_TASK_THREADS",
)


def _gpu_free_memory() -> dict[str, float]:
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=uuid,memory.free", "--format=csv,noheader,nounits"],
            check=True, capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError(f"GPU resource probe failed: {error}") from error
    available = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = [field.strip() for field in line.split(",")]
        try:
            if len(fields) != 2 or not fields[0].startswith("GPU-"):
                raise ValueError("Expected GPU UUID and free memory")
            free_gib = float(fields[1]) / 1024
            if not math.isfinite(free_gib) or free_gib < 0:
                raise ValueError("Invalid free memory")
        except ValueError as error:
            raise ValueError(f"Invalid GPU resource probe response: {line!r}") from error
        available[fields[0]] = free_gib
    return available


@contextmanager
def lease_gpu(plan: dict, probe=None, poll_interval: float = 0.25):
    """Hold a process-shared GPU lock, including inherited backend descriptors."""
    eligible = plan.get("eligible_gpu_uuids")
    if not isinstance(eligible, list) or not eligible or any(
        not isinstance(uuid, str) or not re.fullmatch(r"GPU-[A-Za-z0-9-]+", uuid)
        for uuid in eligible
    ):
        raise ValueError("CUDA tasks require explicit eligible physical GPU UUIDs")
    required = float(plan["gpu_memory_gib"])
    timeout = float(plan.get("gpu_wait_timeout", 600))
    if not math.isfinite(required) or required <= 0:
        raise ValueError("GPU memory requirement must be positive and finite")
    if not math.isfinite(timeout) or timeout < 0 or poll_interval <= 0:
        raise ValueError("GPU wait timeout must be finite and nonnegative")
    lease_directory = Path(plan.get("lease_dir") or Path.home() / ".cache" / "biofeaturefactory" / "gpu-leases")
    lease_directory.mkdir(parents=True, exist_ok=True)
    memory_probe = probe or _gpu_free_memory
    deadline = time.monotonic() + timeout
    selected_descriptor = None
    while selected_descriptor is None:
        for uuid in dict.fromkeys(eligible):
            descriptor = os.open(lease_directory / f"{uuid}.lock", os.O_CREAT | os.O_RDWR, 0o600)
            try:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except OSError as error:
                    if error.errno not in {errno.EACCES, errno.EAGAIN}:
                        raise
                    continue
                available = memory_probe()
                if not set(eligible).intersection(available):
                    raise ValueError("GPU resource probe found none of the planned GPU UUIDs")
                if available.get(uuid, -1) >= required:
                    selected_uuid, selected_descriptor = uuid, descriptor
                    descriptor = None
                    break
            finally:
                if descriptor is not None:
                    os.close(descriptor)
        if selected_descriptor is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise ValueError(
                    f"GPU resource wait timed out after {timeout:g}s: no unlocked eligible "
                    f"GPU has {required:g} GiB free"
                )
            time.sleep(min(poll_interval, remaining))
    try:
        yield selected_uuid, selected_descriptor
    finally:
        os.close(selected_descriptor)


def _file_record(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def _task_names(gene: str, side: str) -> tuple[str, str, str]:
    if not isinstance(gene, str) or not gene or gene in {".", ".."}:
        raise ValueError("Plan requires a valid gene name")
    if Path(gene).name != gene or "\\" in gene or "\x00" in gene:
        raise ValueError("Gene name must not contain a path")
    if side not in {"protein", "codon"}:
        raise ValueError("Plan side must be protein or codon")
    stem = f"{gene}.{side}"
    return f"{stem}.tsv", f"{stem}_adabm_params", f"{stem}.complete.json"


def verify_completion(
    path: str | Path, fingerprint: str,
    artifact_paths: dict[str, str | Path] | None = None,
) -> bool:
    """Check a completion manifest and every named artifact against its hashes."""
    try:
        manifest_path = Path(path)
        manifest = json.loads(manifest_path.read_text())
        tsv_name, params_name, manifest_name = _task_names(
            manifest["gene"], manifest["side"],
        )
        if manifest_path.name != manifest_name or manifest["success"] is not True:
            return False
        if not fingerprint or manifest["fingerprint"] != fingerprint:
            return False
        if manifest["device"] not in {"cuda", "cpu"}:
            return False
        if set(manifest["files"]) != {tsv_name, params_name}:
            return False
        if artifact_paths is not None and set(artifact_paths) != {tsv_name, params_name}:
            return False
        for filename, expected in manifest["files"].items():
            artifact = (
                Path(artifact_paths[filename]) if artifact_paths is not None
                else manifest_path.parent / filename
            )
            if artifact.name != filename or not artifact.is_file() or artifact.stat().st_size == 0:
                return False
            if _file_record(artifact) != expected:
                return False
        return True
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return False


def _write_json(path: Path, contents: dict) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as handle:
        json.dump(contents, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _prepare_command(
    backend_script: Path, backend_args: list[str], plan: dict,
    device: str, attempt: Path, task_directory: Path,
) -> tuple[list[str], Path]:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("--fasta", "-f")
    parser.add_argument("--mutations", "-m")
    parser.add_argument("--msa")
    parser.add_argument("--codon-msa", "-cm")
    parser.add_argument("--protein-params", "-pp")
    parser.add_argument("--codon-params", "-cp")
    parser.add_argument("--validation-log", "-vl")
    parser.add_argument("--output", "-o")
    parser.add_argument("--adabmdca-device", "-ad")
    parser.add_argument("--gene", "-g")
    arguments, remaining = parser.parse_known_args(backend_args)
    gene, side = plan["gene"], plan["side"]
    command = [sys.executable, str(backend_script), *remaining]
    for option in ("fasta", "mutations", "validation_log"):
        value = getattr(arguments, option)
        if value:
            source = (task_directory / value).resolve()
            if option == "fasta" and not source.is_file():
                raise ValueError("The task wrapper requires a single-gene FASTA file")
            command.extend([f"--{option.replace('_', '-')}", str(source)])
    msa_option = "msa" if side == "protein" else "codon_msa"
    msa_value = getattr(arguments, msa_option)
    if msa_value:
        source = (task_directory / msa_value).resolve()
        local_msa = attempt / "inputs" / source.name
        local_msa.parent.mkdir()
        shutil.copyfile(source, local_msa)
        command.extend([f"--{msa_option.replace('_', '-')}", str(local_msa)])
    params_name = f"{gene}.{side}_adabm_params"
    params_directory = attempt / "params"
    params_directory.mkdir()
    params_path = params_directory / params_name
    supplied_params = plan.get("params") if "params" in plan else getattr(
        arguments, f"{side}_params",
    )
    if supplied_params:
        source = (task_directory / supplied_params).resolve()
        if source.is_dir():
            source = source / params_name
        if source.is_file():
            if source.stat().st_size == 0:
                raise ValueError(f"Explicit params are empty: {source}")
            shutil.copyfile(source, params_path)
        elif "params" in plan:
            raise ValueError(f"Explicit params do not exist: {source}")
    params_argument = params_path if params_path.is_file() else params_directory
    command.extend([
        f"--{side}-params", str(params_argument),
        "--output", str(attempt / "output"),
        "--adabmdca-device", device, "--gene", gene,
    ])
    return command, params_path


def _task_environment(device: str, threads: int, gpu_uuid: str | None) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update({variable: str(threads) for variable in THREAD_VARIABLES})
    environment["PYTHONUNBUFFERED"] = "1"
    environment["CUDA_VISIBLE_DEVICES"] = gpu_uuid if device == "cuda" else ""
    return environment


def _run_backend(command: list[str], attempt: Path, environment: dict, descriptor: int | None):
    cuda_oom = False
    with (attempt / "backend.log").open("w") as log:
        process = subprocess.Popen(
            command, cwd=attempt, env=environment, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, errors="replace", bufsize=1,
            pass_fds=(descriptor,) if descriptor is not None else (),
        )
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                sys.stdout.write(line)
                sys.stdout.flush()
                cuda_oom = cuda_oom or bool(CUDA_OOM.search(line))
            return process.wait(), cuda_oom
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            process.stdout.close()


def run_task(
    plan: dict, device: str, backend_script: str | Path, threads: int,
    backend_args: list[str], task_directory: str | Path | None = None,
) -> int:
    if threads < 1:
        raise ValueError("Threads must be positive")
    if device not in {"cuda", "cpu"}:
        raise ValueError("Device must be cuda or cpu")
    tsv_name, params_name, manifest_name = _task_names(plan["gene"], plan["side"])
    if not isinstance(plan.get("fingerprint"), str) or not plan["fingerprint"]:
        raise ValueError("Plan requires a nonempty fingerprint")
    destination = Path(task_directory or Path.cwd()).resolve()
    attempt_root = destination / ".adabmdca-attempts"
    attempt_root.mkdir(exist_ok=True)
    attempt = Path(tempfile.mkdtemp(
        prefix=f"{plan['gene']}.{plan['side']}.{device}.", dir=attempt_root,
    ))
    command, params_path = _prepare_command(
        Path(backend_script).resolve(), backend_args, plan, device, attempt, destination,
    )
    oom_name = f"{plan['gene']}.{plan['side']}.gpu_oom.json"
    for filename in (tsv_name, params_name, manifest_name, oom_name):
        previous = destination / filename
        if previous.exists() or previous.is_symlink():
            previous_directory = attempt / "previous"
            previous_directory.mkdir(exist_ok=True)
            os.replace(previous, previous_directory / filename)
    allocation = lease_gpu(plan) if device == "cuda" else nullcontext((None, None))
    with allocation as (gpu_uuid, descriptor):
        environment = _task_environment(device, threads, gpu_uuid)
        resource_report = attempt / "backend_resources.json"
        environment["BFF_TASK_RESOURCE_REPORT"] = str(resource_report)
        metadata = {
            **plan, "device": device, "threads": threads,
            "attempt": attempt.name, "command": command,
            "cuda_visible_devices": environment["CUDA_VISIBLE_DEVICES"],
        }
        _write_json(attempt / "attempt.json", metadata)
        returncode, cuda_oom = _run_backend(command, attempt, environment, descriptor)
        metadata["backend_returncode"] = returncode
        memory_scale = 1024 ** 3 if sys.platform == "darwin" else 1024 ** 2
        metadata["subprocess_peak_rss_gib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / memory_scale
        if resource_report.is_file():
            try:
                backend_resources = json.loads(resource_report.read_text())
            except (OSError, ValueError):
                backend_resources = None
            if isinstance(backend_resources, dict):
                metadata["backend_resources"] = backend_resources
        _write_json(attempt / "attempt.json", metadata)
    if returncode != 0:
        if device == "cuda" and cuda_oom:
            _write_json(destination / oom_name, {
                **metadata, "success": False, "reason": "cuda_oom",
                "backend_returncode": returncode,
            })
            print(f"CUDA OOM for {plan['gene']}.{plan['side']}; CPU retry required", file=sys.stderr)
            return 0
        print(f"adabmDCA backend failed with exit code {returncode}; see {attempt / 'backend.log'}", file=sys.stderr)
        return returncode if returncode > 0 else 128 - returncode
    tsv_path = attempt / "output" / plan["gene"] / "adabmDCA" / tsv_name
    artifacts = {tsv_name: tsv_path, params_name: params_path}
    for artifact in artifacts.values():
        if not artifact.is_file() or artifact.stat().st_size == 0:
            print(f"adabmDCA backend did not produce a nonempty artifact: {artifact}", file=sys.stderr)
            return 1
    files = {name: _file_record(artifact) for name, artifact in artifacts.items()}
    for filename, artifact in artifacts.items():
        os.replace(artifact, destination / filename)
    _write_json(destination / manifest_name, {
        **metadata, "success": True, "files": files,
    })
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--device", required=True, choices=("cuda", "cpu"))
    parser.add_argument("--backend-script", type=Path, default=Path(__file__).resolve().parents[1] / "adabmdca_pipeline.py")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("backend_args", nargs=argparse.REMAINDER)
    arguments = parser.parse_args(argv)
    backend_args = arguments.backend_args
    if backend_args[:1] == ["--"]:
        backend_args = backend_args[1:]
    try:
        plan = json.loads(arguments.plan.read_text())
        return run_task(
            plan, arguments.device, arguments.backend_script,
            arguments.threads, backend_args,
        )
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"adabmDCA task failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
