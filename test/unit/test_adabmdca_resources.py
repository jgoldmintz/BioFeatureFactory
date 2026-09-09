"""Resource planning with simulated allocations; no GPU execution is claimed."""

from argparse import Namespace
from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from biofeaturefactory.mutation_effects.bin import resource_planner as planner


GPU_INVENTORY = "0, GPU-small-0001, 40960, 32768\n1, GPU-large-0002, 81920, 71680\n"


@pytest.fixture
def hardware():
    return {
        "cpus": 32,
        "memory_gib": 220,
        "gpus": [
            {"index": "0", "uuid": "GPU-first", "memory_gib": 80},
            {"index": "1", "uuid": "GPU-second", "memory_gib": 80},
        ],
    }


@pytest.fixture
def config(hardware, tmp_path):
    return {
        "hardware": hardware,
        "device": "auto",
        "threads": 4,
        "settings": {"model": "bmDCA", "dtype": "float32", "nchains": 64},
        "implementation": {"pipeline": "test-version"},
        "memory_margin": 1.15,
        "overrides": {},
        "lease_dir": str(tmp_path / "leases"),
    }


@pytest.fixture
def inputs(tmp_path):
    files = {
        "fasta": tmp_path / "GENE.fasta",
        "msa": tmp_path / "GENE.msa.a2m",
        "mutations": tmp_path / "GENE_mutations.csv",
    }
    files["fasta"].write_text(">ORF\nATGGCT\n")
    files["msa"].write_text(">ORF\nMA\n>homolog\nM-\n")
    files["mutations"].write_text("mutant\nG4A\n")
    return files


@pytest.fixture
def detected_environment(monkeypatch):
    virtual_files = {}
    original_read = Path.read_text

    def read_virtual(path, *args, **kwargs):
        filename = str(path)
        if filename in virtual_files:
            return virtual_files[filename]
        if filename.startswith("/sys/fs/cgroup/"):
            raise FileNotFoundError(filename)
        return original_read(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read_virtual)
    monkeypatch.setattr(planner.os, "sched_getaffinity", lambda process: set(range(32)), raising=False)
    monkeypatch.setitem(sys.modules, "psutil", Namespace(
        virtual_memory=lambda: Namespace(available=256 * planner.GIB),
    ))
    monkeypatch.setattr(planner.subprocess, "run", lambda *args, **kwargs: Namespace(stdout=GPU_INVENTORY))
    for name in ("CUDA_VISIBLE_DEVICES", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(name, raising=False)
    return virtual_files


def config_args(tmp_path, **overrides):
    values = {
        "adabmdca_model": "bmDCA", "adabmdca_dtype": "float32",
        "adabmdca_nepochs": 50, "adabmdca_tol": 0.001,
        "adabmdca_patience": 3, "adabmdca_check_every": 10,
        "adabmdca_target": 0.95, "adabmdca_lr": 0.01,
        "adabmdca_nchains": 64, "adabmdca_nsweeps": 10, "adabmdca_seed": 0,
        "skip_codon_adabmdca": False, "validation_log": None,
        "resource_hardware": None, "resource_memory_gib": None,
        "resource_cpus": None, "resource_headroom": 0.9,
        "resource_memory_margin": 1.15, "resource_overrides": None,
        "adabmdca_device": "auto", "threads": 4,
        "gpu_lease_dir": tmp_path / "leases", "gpu_wait_timeout": 600,
    }
    values.update(overrides)
    return Namespace(**values)


@pytest.mark.parametrize("visible,expected", [
    (None, ["GPU-small-0001", "GPU-large-0002"]),
    ("1,0", ["GPU-large-0002", "GPU-small-0001"]),
    ("GPU-large", ["GPU-large-0002"]),
    ("0,GPU-small-0001", ["GPU-small-0001"]),
    ("", []), ("-1", []),
])
def test_gpu_visibility_respects_index_uuid_order_and_empty(visible, expected):
    devices = planner.parse_gpus(GPU_INVENTORY, visible, headroom=0.9)
    assert [device["uuid"] for device in devices] == expected
    for device in devices:
        assert device["memory_gib"] == pytest.approx(36 if device["index"] == "0" else 72)


@pytest.mark.parametrize("visible", ["9", "GPU-missing", "GPU-"])
def test_gpu_visibility_rejects_unresolved_or_ambiguous_ids(visible):
    with pytest.raises(ValueError, match="resolve visible GPU"):
        planner.parse_gpus(GPU_INVENTORY, visible)


@pytest.mark.parametrize("free", ["-1", "nan", "inf"])
def test_gpu_inventory_rejects_invalid_free_memory(free):
    with pytest.raises(ValueError):
        planner.parse_gpus(f"0, GPU-first, 81920, {free}\n")


@pytest.mark.parametrize("mode", ["auto", "cuda"])
@pytest.mark.parametrize("free_mib", [0, 20480, 81920])
def test_busy_gpu_is_planned_by_capacity_for_runtime_waiting(config, inputs, mode, free_mib):
    config["device"] = mode
    config["hardware"]["gpus"] = planner.parse_gpus(f"0, GPU-busy, 81920, {free_mib}\n")
    config["overrides"]["GENE.protein"] = {"gpu_memory_gib": 60}
    planner.validate_hardware(config["hardware"])
    plan = planner.plan_task("GENE", "protein", config=config, **inputs)
    assert config["hardware"]["gpus"][0]["memory_gib"] == pytest.approx(72)
    assert plan["device"] == "cuda"
    assert plan["eligible_gpu_uuids"] == ["GPU-busy"]


def test_zero_declared_gpu_capacity_uses_cpu(config, inputs):
    config["hardware"]["gpus"] = [{"uuid": "GPU-disabled", "memory_gib": 0}]
    planner.validate_hardware(config["hardware"])
    plan = planner.plan_task("GENE", "protein", config=config, **inputs)
    assert plan["device"] == "cpu"
    assert plan["eligible_gpu_uuids"] == []


def test_cpu_detection_avoids_gpu_probe(detected_environment, monkeypatch):
    monkeypatch.setattr(planner.subprocess, "run", lambda *args, **kwargs: pytest.fail("CPU mode queried GPUs"))
    hardware = planner.detect_hardware("cpu")
    assert hardware["gpus"] == []
    assert hardware["cpus"] == 31
    assert hardware["memory_gib"] == pytest.approx(256 * 0.9)


@pytest.mark.parametrize("visible", ["", "-1"])
def test_empty_visibility_hides_gpus(detected_environment, monkeypatch, visible):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", visible)
    assert planner.detect_hardware("auto")["gpus"] == []
    with pytest.raises(ValueError, match="no allocated GPU"):
        planner.detect_hardware("cuda")


def test_no_gpu_binary_allows_auto_and_rejects_cuda(detected_environment, monkeypatch):
    def absent_binary(*args, **kwargs):
        raise FileNotFoundError("nvidia-smi")

    monkeypatch.setattr(planner.subprocess, "run", absent_binary)
    assert planner.detect_hardware("auto")["gpus"] == []
    with pytest.raises(ValueError, match="nvidia-smi is unavailable"):
        planner.detect_hardware("cuda")


def test_cgroup_quota_memory_usage_and_explicit_caps(detected_environment):
    detected_environment.update({
        "/sys/fs/cgroup/cpu.max": "800000 100000",
        "/sys/fs/cgroup/memory.max": str(100 * planner.GIB),
        "/sys/fs/cgroup/memory.current": str(20 * planner.GIB),
    })
    uncapped = planner.detect_hardware("cpu")
    capped = planner.detect_hardware("cpu", memory_gib=50, cpus=6)
    assert uncapped["cpus"] == 7
    assert uncapped["memory_gib"] == pytest.approx(72)
    assert capped["cpus"] == 6
    assert capped["memory_gib"] == 50


def test_cgroup_v1_memory_and_slurm_limits(detected_environment, monkeypatch):
    detected_environment.update({
        "/sys/fs/cgroup/memory/memory.limit_in_bytes": str(100 * planner.GIB),
        "/sys/fs/cgroup/memory/memory.usage_in_bytes": str(20 * planner.GIB),
    })
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "32768")
    hardware = planner.detect_hardware("cpu")
    assert hardware["cpus"] == 7
    assert hardware["memory_gib"] == pytest.approx(28.8)


@pytest.mark.parametrize("controller_path", ["/sys/fs/cgroup/cpu", "/sys/fs/cgroup/cpu,cpuacct"])
@pytest.mark.parametrize("quota,expected_cores", [(800000, 7), (150000, 1), (-1, 31)])
def test_cgroup_v1_cpu_quotas_preserve_reserved_core(detected_environment, controller_path, quota, expected_cores):
    detected_environment.update({
        f"{controller_path}/cpu.cfs_quota_us": str(quota),
        f"{controller_path}/cpu.cfs_period_us": "100000",
    })
    assert planner.detect_hardware("cpu")["cpus"] == expected_cores


def test_cgroup_v1_does_not_relax_more_restrictive_v2_quota(detected_environment):
    detected_environment.update({
        "/sys/fs/cgroup/cpu.max": "400000 100000",
        "/sys/fs/cgroup/cpu/cpu.cfs_quota_us": "800000",
        "/sys/fs/cgroup/cpu/cpu.cfs_period_us": "100000",
    })
    assert planner.detect_hardware("cpu")["cpus"] == 3


def test_cuda_index_addresses_visible_order(detected_environment, monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "GPU-large,GPU-small")
    hardware = planner.detect_hardware("cuda:1")
    assert [gpu["uuid"] for gpu in hardware["gpus"]] == ["GPU-small-0001"]
    with pytest.raises(ValueError, match="outside visible GPU"):
        planner.detect_hardware("cuda:2")


@pytest.mark.parametrize("source", ["argument", "profile"])
def test_declared_hardware_honors_cuda_index(tmp_path, hardware, source):
    args = config_args(tmp_path, adabmdca_device="cuda:1")
    if source == "profile":
        args.resource_hardware = tmp_path / "hardware.json"
        args.resource_hardware.write_text(json.dumps(hardware))
        result = planner.make_config(args)
    else:
        result = planner.make_config(args, hardware=hardware)
    assert [gpu["uuid"] for gpu in result["hardware"]["gpus"]] == ["GPU-second"]


def test_declared_hardware_honors_cpu_and_ram_caps(tmp_path, hardware):
    args = config_args(tmp_path, resource_cpus=8, resource_memory_gib=100, adabmdca_device="cpu")
    result = planner.make_config(args, hardware=hardware)
    assert result["hardware"] == {"cpus": 8, "memory_gib": 100, "gpus": []}


@pytest.mark.parametrize("margin", [-1, 0, 0.5, float("nan"), float("inf")])
def test_invalid_memory_margins_are_rejected(tmp_path, hardware, margin):
    args = config_args(tmp_path, resource_memory_margin=margin)
    with pytest.raises(ValueError):
        planner.make_config(args, hardware=hardware)


@pytest.mark.parametrize("override", [
    {"GENE.protein": {"cpu_mem_gib": 20}},
    {"GENE.protein": [20]},
    [20],
])
def test_malformed_override_profiles_are_rejected(tmp_path, hardware, override):
    profile = tmp_path / "overrides.json"
    profile.write_text(json.dumps(override))
    args = config_args(tmp_path, resource_overrides=profile)
    with pytest.raises(ValueError):
        planner.make_config(args, hardware=hardware)


@pytest.mark.parametrize("contents,codon,shape", [
    (">first\nAcD.E-F\n>second\nAD-EF\n", False, (2, 5)),
    (">first\nAA\nCG\n>second\nA-CG\n", False, (2, 4)),
    (">first\nATGgct---\n>second\nATGGCT---\n", True, (2, 3)),
])
def test_alignment_shape_protein_insertions_and_codon_columns(tmp_path, contents, codon, shape):
    msa = tmp_path / "alignment.fasta"
    msa.write_text(contents)
    assert planner.alignment_shape(msa, codon=codon) == shape


@pytest.mark.parametrize("contents,codon", [
    ("", False), ("ACGT\n", False), (">one\n", False),
    (">one\nAC\n>two\nA\n", False), (">one\nACGT\n", True),
])
def test_invalid_alignments_are_rejected(tmp_path, contents, codon):
    msa = tmp_path / "bad.fasta"
    msa.write_text(contents)
    with pytest.raises(ValueError):
        planner.alignment_shape(msa, codon=codon)


@pytest.mark.parametrize("model", ["bmDCA", "eaDCA", "edDCA", "pseudoDCA"])
@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_all_models_and_dtypes_have_positive_conservative_estimates(model, dtype):
    estimate = planner.estimate_memory(819, 758, 65, model, dtype, 64)
    assert all(value > 0 for value in estimate.values())
    assert estimate["cpu_memory_gib"] >= estimate["gpu_host_memory_gib"]
    if dtype == "float64":
        single = planner.estimate_memory(819, 758, 65, model, "float32", 64)
        assert estimate["gpu_memory_gib"] > single["gpu_memory_gib"]


def test_npm1_estimate_exceeds_observed_memory_lower_bound():
    estimate = planner.estimate_memory(819, 758, 65, "bmDCA", "float32", 64)
    assert estimate["gpu_memory_gib"] > 64.20 + 9.04


@pytest.mark.parametrize("side,states", [("protein", 21), ("codon", 65)])
def test_planning_uses_side_specific_states(config, inputs, side, states):
    if side == "codon":
        inputs["msa"].write_text(">ORF\nATGGCT\n>homolog\nATG---\n")
    plan = planner.plan_task("GENE", side, config=config, **inputs)
    assert (plan["sequences"], plan["sites"], plan["states"]) == (2, 2, states)


def test_heterogeneous_gpus_choose_only_fitting_devices(config, inputs):
    config["hardware"]["gpus"][0]["memory_gib"] = 40
    config["overrides"]["GENE.protein"] = {"gpu_memory_gib": 60}
    plan = planner.plan_task("GENE", "protein", config=config, **inputs)
    assert plan["device"] == "cuda"
    assert plan["eligible_gpu_uuids"] == ["GPU-second"]


@pytest.mark.parametrize("mode,expected", [("auto", "cpu"), ("cpu", "cpu"), ("cuda", None)])
def test_gpu_oversize_obeys_cpu_fallback_and_explicit_device(config, inputs, mode, expected):
    config["device"] = mode
    config["overrides"]["GENE.protein"] = {"gpu_memory_gib": 100, "cpu_memory_gib": 100}
    if expected is None:
        with pytest.raises(ValueError, match="CUDA request cannot fit"):
            planner.plan_task("GENE", "protein", config=config, **inputs)
    else:
        plan = planner.plan_task("GENE", "protein", config=config, **inputs)
        assert plan["device"] == expected
        assert plan["cpu_fallback_allowed"] == (mode == "auto")


def test_gpu_host_ram_shortfall_falls_back_to_fitting_cpu(config, inputs):
    config["overrides"]["GENE.protein"] = {"gpu_host_memory_gib": 250, "cpu_memory_gib": 100}
    assert planner.plan_task("GENE", "protein", config=config, **inputs)["device"] == "cpu"


def test_no_gpu_allocation_uses_cpu(config, inputs):
    config["hardware"]["gpus"] = []
    assert planner.plan_task("GENE", "protein", config=config, **inputs)["device"] == "cpu"


def test_unschedulable_task_is_rejected(config, inputs):
    config["overrides"]["GENE.protein"] = {"gpu_memory_gib": 100, "cpu_memory_gib": 250}
    with pytest.raises(ValueError, match="unschedulable"):
        planner.plan_task("GENE", "protein", config=config, **inputs)


def test_prebuilt_params_are_planned_for_cpu_even_with_forced_cuda(config, inputs, tmp_path):
    params = tmp_path / "GENE.protein_adabm_params"
    params.write_text("h 0 A 0\nh 1 A 0\n")
    config["device"] = "cuda"
    plan = planner.plan_task("GENE", "protein", config=config, params=params, **inputs)
    assert plan["device"] == "cpu"
    assert plan["gpu_memory_gib"] == 0
    assert plan["params"] == str(params.resolve())


def test_seven_task_simulated_allocation_has_two_gpu_and_three_cpu_slots(config, inputs):
    plans = []
    for task_index in range(7):
        gene = f"GENE{task_index}"
        config["overrides"][f"{gene}.protein"] = {
            "gpu_memory_gib": 64 if task_index < 3 else 96,
            "gpu_host_memory_gib": 10,
            "cpu_memory_gib": 100 if task_index < 3 else 64,
        }
        plans.append(planner.plan_task(gene, "protein", config=config, **inputs))

    assert [plan["device"] for plan in plans] == ["cuda"] * 3 + ["cpu"] * 4
    admitted_gpu = plans[:2]
    admitted_cpu = plans[3:6]
    used_ram = sum(plan["gpu_host_memory_gib"] for plan in admitted_gpu)
    used_ram += sum(plan["cpu_memory_gib"] for plan in admitted_cpu)
    used_threads = sum(plan["threads"] for plan in admitted_gpu + admitted_cpu)
    assert used_ram == 212 <= config["hardware"]["memory_gib"]
    assert used_threads == 20 <= config["hardware"]["cpus"]
    assert plans[6]["cpu_memory_gib"] > config["hardware"]["memory_gib"] - used_ram
    for gpu, active in zip(config["hardware"]["gpus"], admitted_gpu):
        assert active["gpu_memory_gib"] <= gpu["memory_gib"]
        assert plans[2]["gpu_memory_gib"] > gpu["memory_gib"] - active["gpu_memory_gib"]


@pytest.mark.parametrize("changed", ["fasta", "msa", "mutations", "params", "model", "implementation"])
def test_fingerprint_changes_with_scientific_inputs(config, inputs, tmp_path, changed):
    params = tmp_path / "params.dat"
    params.write_text("h 0 A 0\n")
    before = planner.task_fingerprint("GENE", "protein", config=config, params=params, **inputs)
    if changed == "model":
        config["settings"]["model"] = "pseudoDCA"
    elif changed == "implementation":
        config["implementation"]["pipeline"] = "another-version"
    else:
        path = params if changed == "params" else inputs[changed]
        path.write_text(path.read_text() + "changed\n")
    after = planner.task_fingerprint("GENE", "protein", config=config, params=params, **inputs)
    assert after != before


def test_fingerprint_ignores_allocation_and_input_paths(config, inputs, tmp_path):
    params = tmp_path / "original_params.dat"
    params.write_text("h 0 A 0\n")
    before = planner.task_fingerprint("GENE", "protein", config=config, params=params, **inputs)
    relocated = {}
    for kind, path in inputs.items():
        target = tmp_path / f"relocated_{path.name}"
        target.write_bytes(path.read_bytes())
        relocated[kind] = target
    alternate = deepcopy(config)
    alternate.update({"device": "cpu", "threads": 8, "lease_dir": "/another/lease/path"})
    alternate["hardware"] = {"cpus": 16, "memory_gib": 128, "gpus": []}
    alternate["overrides"] = {"GENE.protein": {"cpu_memory_gib": 20}}
    relocated_params = tmp_path / "relocated_params.dat"
    relocated_params.write_bytes(params.read_bytes())
    after = planner.task_fingerprint("GENE", "protein", config=alternate, params=relocated_params, **relocated)
    assert after == before


@pytest.mark.parametrize("key,value", [
    ("threads", 0), ("threads", -1), ("threads", 0.5), ("threads", 100),
    ("cpu_memory_gib", 0), ("gpu_memory_gib", -1),
    ("gpu_host_memory_gib", float("nan")), ("cpu_memory_gib", float("inf")),
])
def test_invalid_task_overrides_are_rejected(config, inputs, key, value):
    config["overrides"]["GENE.protein"] = {key: value}
    with pytest.raises(ValueError):
        planner.plan_task("GENE", "protein", config=config, **inputs)


@pytest.mark.parametrize("changes", [
    {"cpus": 0}, {"cpus": -1}, {"cpus": 1.5}, {"memory_gib": 0}, {"memory_gib": float("nan")},
    {"gpus": [{"uuid": "", "memory_gib": 80}]},
    {"gpus": [{"uuid": "GPU-first", "memory_gib": -1}]},
    {"gpus": [{"uuid": "GPU-first", "memory_gib": 80}, {"uuid": "GPU-first", "memory_gib": 80}]},
])
def test_invalid_hardware_profiles_are_rejected(hardware, changes):
    hardware.update(changes)
    with pytest.raises(ValueError):
        planner.validate_hardware(hardware)


@pytest.mark.parametrize("defer_errors", [False, True])
@pytest.mark.parametrize("failure_kind", ["oversize", "invalid_msa", "missing_msa"])
def test_cli_defers_only_task_planning_errors(tmp_path, config, inputs, capsys, defer_errors, failure_kind):
    if failure_kind == "oversize":
        config["overrides"]["GENE.protein"] = {"gpu_memory_gib": 100, "cpu_memory_gib": 250}
    elif failure_kind == "invalid_msa":
        inputs["msa"].write_text("")
    else:
        inputs["msa"] = tmp_path / "absent.fasta"
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    arguments = ["plan", "--gene", "GENE", "--side", "protein", "--config", str(config_file),
                 "--output", str(output)]
    for name, path in inputs.items():
        arguments.extend([f"--{name}", str(path)])
    if not defer_errors:
        with pytest.raises(SystemExit) as failure:
            planner.main(arguments)
        assert failure.value.code == 1
        assert not output.exists()
        assert "Resource planning failed" in capsys.readouterr().err
        return
    planner.main([*arguments, "--defer-errors"])
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""
    plan = json.loads(output.read_text())
    assert plan["device"] == "blocked"
    assert "threads" not in plan
    assert plan["resource_error"] == {
        "gene": "GENE", "side": "protein", "backend": "adabmdca",
        "message": plan["resource_error"]["message"], "stage": "resource_planning",
    }
    assert plan["resource_error"]["message"]


def test_cli_deferred_mode_keeps_valid_plan_unchanged(tmp_path, config, inputs):
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    arguments = ["plan", "--gene", "GENE", "--side", "protein", "--config", str(config_file),
                 "--output", str(output), "--defer-errors"]
    for name, path in inputs.items():
        arguments.extend([f"--{name}", str(path)])
    planner.main(arguments)
    assert json.loads(output.read_text()) == planner.plan_task("GENE", "protein", config=config, **inputs)


@pytest.mark.parametrize("failure_kind", ["config", "threads", "hardware"])
def test_cli_deferred_mode_does_not_hide_global_input_errors(tmp_path, config, inputs, failure_kind):
    config_file = tmp_path / "config.json"
    if failure_kind == "threads":
        config["threads"] = 0
    elif failure_kind == "hardware":
        config["hardware"]["memory_gib"] = 0
    if failure_kind != "config":
        config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    arguments = ["plan", "--gene", "GENE", "--side", "protein", "--config", str(config_file),
                 "--output", str(output), "--defer-errors"]
    for name, path in inputs.items():
        arguments.extend([f"--{name}", str(path)])
    with pytest.raises(SystemExit) as failure:
        planner.main(arguments)
    assert failure.value.code == 1
    assert not output.exists()


@pytest.mark.parametrize("blocked", [False, True])
def test_cli_deferred_mode_does_not_hide_output_write_errors(tmp_path, config, inputs, blocked):
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    if blocked:
        inputs["msa"] = tmp_path / "missing.fasta"
    arguments = ["plan", "--gene", "GENE", "--side", "protein", "--config", str(config_file),
                 "--output", str(tmp_path), "--defer-errors"]
    for name, path in inputs.items():
        arguments.extend([f"--{name}", str(path)])
    with pytest.raises(OSError):
        planner.main(arguments)
