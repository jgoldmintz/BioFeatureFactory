"""Resource estimates and deterministic placement for mutation-effects tasks."""

from __future__ import annotations

import argparse
import csv
import copy
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import subprocess
import sys


GIB = 1024 ** 3
SCHEMA_VERSION = 1
MODEL_TENSORS = {"bmDCA": 8, "eaDCA": 10, "edDCA": 10, "pseudoDCA": 6}


def positive(value, name):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def positive_integer(value, name):
    number = positive(value, name)
    if not number.is_integer():
        raise ValueError(f"{name} must be an integer")
    return int(number)


def automatic_threads(hardware, tasks):
    """Greedily admit jobs by RAM then ID, matching each GPU job to one UUID.

    Fixed CPU requests and one core per automatic job must fit the usable CPU
    budget. Remaining cores are shared by admitted automatic jobs. This static,
    deterministic estimate is not an optimal packing or dynamic rebalancing.
    """
    cpu_budget = positive_integer(hardware["cpus"], "CPU budget")
    memory_budget = positive(hardware["memory_gib"], "RAM budget")
    available_gpus = {gpu["uuid"] for gpu in hardware.get("gpus", [])}
    candidates = []
    for task in tasks:
        fixed_threads = task.get("threads")
        if fixed_threads is not None:
            fixed_threads = positive_integer(fixed_threads, "threads")
            if fixed_threads > cpu_budget:
                raise ValueError(f"{task['id']}: {fixed_threads} threads exceed CPU budget {cpu_budget}")
        requested_gpus = set(task.get("eligible_gpu_uuids", []))
        candidates.append({
            "id": str(task["id"]),
            "memory_gib": positive(task["memory_gib"], "task RAM"),
            "threads": fixed_threads,
            "gpu_required": bool(requested_gpus),
            "eligible_gpus": sorted(requested_gpus & available_gpus),
        })
    candidates.sort(key=lambda task: (task["memory_gib"], task["id"]))
    gpu_assignments = {}

    def assign_gpu(candidate_index, visited):
        for uuid in candidates[candidate_index]["eligible_gpus"]:
            if uuid in visited:
                continue
            visited.add(uuid)
            previous = gpu_assignments.get(uuid)
            if previous is None or assign_gpu(previous, visited):
                gpu_assignments[uuid] = candidate_index
                return True
        return False

    memory_used = 0.0
    minimum_cpus = 0
    fixed_cpus = 0
    automatic_jobs = 0
    concurrent_jobs = 0
    for candidate_index, task in enumerate(candidates):
        required_cpus = task["threads"] if task["threads"] is not None else 1
        if memory_used + task["memory_gib"] > memory_budget or minimum_cpus + required_cpus > cpu_budget:
            continue
        if task["gpu_required"] and not assign_gpu(candidate_index, set()):
            continue
        memory_used += task["memory_gib"]
        minimum_cpus += required_cpus
        concurrent_jobs += 1
        if task["threads"] is None:
            automatic_jobs += 1
        else:
            fixed_cpus += task["threads"]
    return {
        "threads": max(1, (cpu_budget - fixed_cpus) // max(1, automatic_jobs)),
        "concurrent_jobs": concurrent_jobs,
    }


def file_digest(filename):
    digest = hashlib.sha256()
    with open(filename, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def available_memory_gib():
    try:
        import psutil
        available = psutil.virtual_memory().available
    except ImportError:
        if sys.platform.startswith("linux"):
            entries = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
            available = int(entries["MemAvailable"].split()[0]) * 1024
        else:
            available = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    for limit_path, usage_path in (
        ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
        ("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes"),
    ):
        try:
            limit = int(Path(limit_path).read_text())
            usage = int(Path(usage_path).read_text())
            available = min(available, max(0, limit - usage))
        except (OSError, ValueError):
            pass
    return available / GIB


def parse_gpus(output, visible=None, headroom=0.9):
    devices = []
    for row in csv.reader(output.splitlines(), skipinitialspace=True):
        if len(row) != 4:
            raise ValueError("Unexpected nvidia-smi GPU inventory")
        index, uuid, total, free = (part.strip() for part in row)
        free_memory = float(free)
        if not math.isfinite(free_memory) or free_memory < 0:
            raise ValueError("GPU free memory must be finite and nonnegative")
        capacity = positive(total, "GPU memory") / 1024
        devices.append({"index": index, "uuid": uuid, "memory_gib": capacity * headroom})
    if visible is not None:
        requested = [item.strip() for item in visible.split(",") if item.strip()]
        if not requested or requested == ["-1"]:
            return []
        selected = []
        for identifier in requested:
            matches = [device for device in devices if identifier == device["index"] or device["uuid"].startswith(identifier)]
            if len(matches) != 1:
                raise ValueError(f"Cannot resolve visible GPU {identifier!r}; supply an allocated hardware profile")
            if matches[0] not in selected:
                selected.append(matches[0])
        devices = selected
    return devices


def detect_hardware(device="auto", memory_gib=None, cpus=None, headroom=0.9):
    if not 0 < headroom < 1:
        raise ValueError("Resource headroom must be between zero and one")
    cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota != "max":
            cores = min(cores, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        pass
    for controller_path in ("/sys/fs/cgroup/cpu", "/sys/fs/cgroup/cpu,cpuacct"):
        try:
            quota = int((Path(controller_path) / "cpu.cfs_quota_us").read_text())
            period = int((Path(controller_path) / "cpu.cfs_period_us").read_text())
            if quota > 0 and period > 0:
                cores = min(cores, max(1, quota // period))
        except (OSError, ValueError):
            pass
    if os.environ.get("SLURM_CPUS_PER_TASK"):
        cores = min(cores, int(os.environ["SLURM_CPUS_PER_TASK"]))
    cores = max(1, cores - 1)
    usable = available_memory_gib() * headroom
    if os.environ.get("SLURM_MEM_PER_NODE"):
        usable = min(usable, float(os.environ["SLURM_MEM_PER_NODE"]) / 1024 * headroom)
    if memory_gib is not None:
        usable = min(usable, positive(memory_gib, "RAM budget"))
    if cpus is not None:
        cores = min(cores, int(positive(cpus, "CPU budget")))
    devices = []
    if device != "cpu" and os.environ.get("CUDA_VISIBLE_DEVICES") not in ("", "-1"):
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=index,uuid,memory.total,memory.free", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=10, check=True,
            )
        except FileNotFoundError:
            if device.startswith("cuda"):
                raise ValueError("CUDA requested but nvidia-smi is unavailable")
        except (subprocess.SubprocessError, OSError) as error:
            raise ValueError(f"Cannot inventory GPU resources: {error}") from error
        else:
            devices = parse_gpus(result.stdout, os.environ.get("CUDA_VISIBLE_DEVICES"), headroom)
    if device.startswith("cuda:"):
        index = int(device.split(":", 1)[1])
        if index < 0 or index >= len(devices):
            raise ValueError(f"Requested device {device} is outside visible GPU allocation")
        devices = [devices[index]]
    if device.startswith("cuda") and not devices:
        raise ValueError("CUDA requested but no allocated GPU is available")
    return {"cpus": cores, "memory_gib": usable, "gpus": devices}


def alignment_shape(filename, codon=False):
    width = None
    current = 0
    count = 0
    with open(filename) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count:
                    if current != width:
                        raise ValueError(f"Unequal alignment widths in {filename}")
                count += 1
                current = 0
            else:
                if not count:
                    raise ValueError(f"Missing FASTA header in {filename}")
                current += len(line) if codon else sum(not char.islower() and char != "." for char in line)
                if count == 1:
                    width = current
    if not count or not width or current != width:
        raise ValueError(f"Empty or unaligned MSA: {filename}")
    if codon and width % 3:
        raise ValueError(f"Codon alignment width is not divisible by three: {filename}")
    return count, width // 3 if codon else width


def estimate_memory(sequences, sites, states, model, dtype, chains, margin=1.15, prebuilt=False):
    """Conservative, uncalibrated tensor envelope; measured overrides take precedence."""
    scalar = 8 if dtype == "float64" else 4
    matrix = sites * sites * states * states * scalar / GIB
    data = sequences * sites * states * scalar / GIB
    weighting = sequences * sequences * scalar / GIB
    chain_memory = chains * sites * states * scalar / GIB if model != "pseudoDCA" else 0
    scoring = 2 * sites * sites * states * states * 4 / GIB + 2
    if prebuilt:
        return {"gpu_memory_gib": 0.0, "gpu_host_memory_gib": scoring * margin, "cpu_memory_gib": scoring * margin}
    multiplier = MODEL_TENSORS[model]
    mask = sites * sites * states * states / GIB
    training = multiplier * matrix + mask + 2 * data + weighting + 2 * chain_memory + 2
    return {
        "gpu_memory_gib": training * margin,
        "gpu_host_memory_gib": (scoring + data + 2) * margin,
        "cpu_memory_gib": max(training, scoring) * margin,
    }


def task_settings(config, gene, side):
    settings = dict(config["settings"])
    route = config.get("routing", {}).get(gene)
    if route is not None:
        settings["score_missense_codon"] = side == "codon" and route.get("score_missense_codon", False)
        if side == "protein":
            settings["skip_codon"] = settings.get("skip_codon", False) or not route["codon"]
    return settings


def task_fingerprint(gene, side, fasta, msa, mutations, config, params=None):
    inputs = {"fasta": file_digest(fasta), "msa": file_digest(msa), "mutations": file_digest(mutations)}
    if params:
        inputs["params"] = file_digest(params)
    payload = {
        "schema": SCHEMA_VERSION, "gene": gene, "side": side,
        "inputs": inputs, "settings": task_settings(config, gene, side),
        "implementation": config.get("implementation", {}),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def plan_task(gene, side, fasta, msa, mutations, config, params=None):
    if side not in {"protein", "codon"}:
        raise ValueError("Task side must be protein or codon")
    hardware = config["hardware"]
    settings = task_settings(config, gene, side)
    if params:
        from biofeaturefactory.lib.utility import find_gene_file
        params = find_gene_file(params, gene, [f"*.{side}_adabm_params", f"*.{side}.dat", "*.dat"])
        if params:
            params = str(Path(params).resolve())
    sequences, sites = alignment_shape(msa, codon=side == "codon")
    states = 65 if side == "codon" else len(settings.get("protein_alphabet", "-ACDEFGHIKLMNPQRSTVWY"))
    estimates = estimate_memory(sequences, sites, states, settings["model"], settings["dtype"], settings["nchains"], config.get("memory_margin", 1.15), bool(params))
    override = config.get("overrides", {}).get(f"{gene}.{side}", {})
    for key in estimates:
        if key in override:
            estimates[key] = positive(override[key], key)
    threads = positive_integer(override.get("threads", config["threads"]), "threads")
    if threads > hardware["cpus"]:
        raise ValueError(f"{gene}/{side}: {threads} threads exceed CPU budget {hardware['cpus']}")
    eligible = [gpu for gpu in hardware["gpus"] if gpu["memory_gib"] >= estimates["gpu_memory_gib"]]
    mode = config["device"]
    cpu_fits = estimates["cpu_memory_gib"] <= hardware["memory_gib"]
    gpu_fits = eligible and estimates["gpu_host_memory_gib"] <= hardware["memory_gib"]
    if params or mode == "cpu":
        device = "cpu"
    elif gpu_fits:
        device = "cuda"
    elif mode == "auto":
        device = "cpu"
    else:
        raise ValueError(f"{gene}/{side}: CUDA request cannot fit any visible GPU (estimate {estimates['gpu_memory_gib']:.2f} GiB)")
    if device == "cpu" and not cpu_fits:
        raise ValueError(f"{gene}/{side}: unschedulable; CPU estimate {estimates['cpu_memory_gib']:.2f} GiB exceeds RAM budget {hardware['memory_gib']:.2f} GiB")
    return {
        "gene": gene, "side": side, "device": device, "threads": threads,
        "settings": dict(settings),
        "sites": sites, "sequences": sequences, "states": states, **estimates,
        "params": params,
        "eligible_gpu_uuids": [item["uuid"] for item in eligible],
        "lease_dir": config.get("lease_dir"),
        "gpu_wait_timeout": config.get("gpu_wait_timeout", 600),
        "cpu_fallback_allowed": mode == "auto" and cpu_fits,
        "estimate_basis": "explicit override" if override else "uncalibrated conservative tensor envelope",
        "fingerprint": task_fingerprint(gene, side, fasta, msa, mutations, config, params),
    }


def validate_hardware(hardware):
    hardware["cpus"] = positive_integer(hardware["cpus"], "CPU budget")
    hardware["memory_gib"] = positive(hardware["memory_gib"], "RAM budget")
    if hardware["cpus"] < 1:
        raise ValueError("At least one CPU is required")
    seen = set()
    for gpu in hardware["gpus"]:
        if not isinstance(gpu["uuid"], str) or not gpu["uuid"] or gpu["uuid"] in seen:
            raise ValueError("GPU UUIDs must be unique and nonempty")
        seen.add(gpu["uuid"])
        gpu["memory_gib"] = float(gpu["memory_gib"])
        if not math.isfinite(gpu["memory_gib"]) or gpu["memory_gib"] < 0:
            raise ValueError("GPU budget must be finite and nonnegative")
    return hardware


def make_config(args, hardware=None):
    settings = {name.removeprefix("adabmdca_"): getattr(args, name) for name in (
        "adabmdca_model", "adabmdca_dtype", "adabmdca_nepochs", "adabmdca_tol", "adabmdca_patience",
        "adabmdca_check_every", "adabmdca_target", "adabmdca_lr", "adabmdca_nchains", "adabmdca_nsweeps", "adabmdca_seed",
    )}
    settings["skip_codon"] = args.skip_codon_adabmdca
    settings["validation_log"] = file_digest(args.validation_log) if args.validation_log else None
    implementation = {name: file_digest(Path(__file__).parent.parent / name) for name in ("adabmdca_pipeline.py", "bin/codon_encoding.py")}
    for package in ("adabmDCA", "torch", "numpy"):
        try:
            implementation[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            implementation[package] = "unavailable"
    declared_hardware = hardware is not None or args.resource_hardware is not None
    if hardware is None:
        hardware = json.loads(args.resource_hardware.read_text()) if args.resource_hardware else detect_hardware(args.adabmdca_device, args.resource_memory_gib, args.resource_cpus, args.resource_headroom)
    hardware = validate_hardware(copy.deepcopy(hardware))
    if args.resource_cpus is not None:
        hardware["cpus"] = min(hardware["cpus"], positive_integer(args.resource_cpus, "CPU budget"))
    if args.resource_memory_gib is not None:
        hardware["memory_gib"] = min(hardware["memory_gib"], positive(args.resource_memory_gib, "RAM budget"))
    if args.adabmdca_device == "cpu":
        hardware["gpus"] = []
    elif declared_hardware and args.adabmdca_device.startswith("cuda:"):
        index = int(args.adabmdca_device.split(":", 1)[1])
        if index < 0 or index >= len(hardware["gpus"]):
            raise ValueError(f"{args.adabmdca_device} is outside the declared GPU allocation")
        hardware["gpus"] = [hardware["gpus"][index]]
    if positive(args.resource_memory_margin, "memory margin") < 1:
        raise ValueError("Memory margin must be at least one")
    overrides = json.loads(args.resource_overrides.read_text()) if args.resource_overrides else {}
    if not isinstance(overrides, dict):
        raise ValueError("Resource overrides must be a gene.side mapping")
    for key, values in overrides.items():
        if not isinstance(values, dict) or set(values) - {"gpu_memory_gib", "cpu_memory_gib", "gpu_host_memory_gib", "threads"}:
            raise ValueError(f"Invalid resource override for {key}")
        for name, value in values.items():
            (positive_integer if name == "threads" else positive)(value, name)
    return {
        "schema": SCHEMA_VERSION, "hardware": hardware,
        "device": args.adabmdca_device,
        "threads": 1 if args.threads is None else positive_integer(args.threads, "threads"),
        "threads_explicit": args.threads is not None,
        "settings": settings, "implementation": implementation,
        "memory_margin": args.resource_memory_margin,
        "overrides": overrides,
        "lease_dir": str(args.gpu_lease_dir.resolve()),
        "gpu_wait_timeout": args.gpu_wait_timeout,
    }


def resource_error(gene, side, backend, error):
    return {
        "gene": gene, "side": side, "backend": backend,
        "message": str(error), "stage": "resource_planning",
    }


def blocked_plan(gene, side, backend, error):
    return {
        "gene": gene, "side": side, "backend": backend, "device": "blocked",
        "resource_error": resource_error(gene, side, backend, error),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan"])
    for name in ("gene", "side", "msa", "fasta", "mutations", "config", "output"):
        parser.add_argument(f"--{name}", required=True, choices=("protein", "codon") if name == "side" else None)
    parser.add_argument("--params")
    parser.add_argument("--skip-codon", action="store_true")
    parser.add_argument("--defer-errors", action="store_true")
    args = parser.parse_args(argv)
    try:
        config = json.loads(Path(args.config).read_text())
        if not isinstance(config, dict):
            raise ValueError("Resource config must be a mapping")
        validate_hardware(config["hardware"])
        positive_integer(config["threads"], "threads")
        if positive(config.get("memory_margin", 1.15), "memory margin") < 1:
            raise ValueError("Memory margin must be at least one")
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Resource planning failed: {error}\n")
    try:
        plan = plan_task(args.gene, args.side, args.fasta, args.msa, args.mutations, config, args.params)
    except (OSError, ValueError, KeyError) as error:
        if not args.defer_errors:
            parser.exit(1, f"Resource planning failed: {error}\n")
        plan = blocked_plan(args.gene, args.side, "adabmdca", error)
    Path(args.output).write_text(json.dumps(plan, indent=2) + "\n")
    if "resource_error" not in plan:
        print(f"[resources] {args.gene}/{args.side}: {plan['device']}; GPU {plan['gpu_memory_gib']:.2f} GiB, CPU {plan['cpu_memory_gib']:.2f} GiB, host-on-GPU {plan['gpu_host_memory_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
