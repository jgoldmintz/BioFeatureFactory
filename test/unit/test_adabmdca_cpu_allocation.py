"""Static CPU sharing against simulated RAM and GPU allocations."""

from argparse import Namespace
from copy import deepcopy
import json

import pytest

from biofeaturefactory.mutation_effects.bin import resource_planner as planner


@pytest.fixture
def hardware():
    return {
        "cpus": 20,
        "memory_gib": 220,
        "gpus": [
            {"index": "0", "uuid": "GPU-first", "memory_gib": 80},
            {"index": "1", "uuid": "GPU-second", "memory_gib": 80},
        ],
    }


def task(identifier, memory=10, gpus=(), threads=None):
    return {
        "id": identifier, "memory_gib": memory,
        "eligible_gpu_uuids": list(gpus), "threads": threads,
    }


def config_args(tmp_path, **overrides):
    values = {
        "adabmdca_model": "pseudoDCA", "adabmdca_dtype": "float32",
        "adabmdca_nepochs": 50, "adabmdca_tol": 0.001,
        "adabmdca_patience": 3, "adabmdca_check_every": 10,
        "adabmdca_target": 0.95, "adabmdca_lr": 0.01,
        "adabmdca_nchains": 64, "adabmdca_nsweeps": 10, "adabmdca_seed": 0,
        "skip_codon_adabmdca": False, "validation_log": None,
        "resource_hardware": None, "resource_memory_gib": None,
        "resource_cpus": None, "resource_headroom": 0.9,
        "resource_memory_margin": 1.15, "resource_overrides": None,
        "adabmdca_device": "auto", "threads": None,
        "gpu_lease_dir": tmp_path / "leases", "gpu_wait_timeout": 600,
    }
    values.update(overrides)
    return Namespace(**values)


def test_two_cpu_jobs_share_already_reserved_twenty_cores(hardware):
    result = planner.automatic_threads(hardware, [task("first"), task("second")])
    assert result == {"threads": 10, "concurrent_jobs": 2}


def test_three_gpu_jobs_share_only_two_gpu_slots(hardware):
    jobs = [task(str(index), gpus=("GPU-first", "GPU-second")) for index in range(3)]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 10, "concurrent_jobs": 2}


def test_ram_limited_single_job_gets_all_usable_cores(hardware):
    hardware["memory_gib"] = 100
    jobs = [task("first", memory=64), task("second", memory=64)]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 20, "concurrent_jobs": 1}


def test_gpu_matching_reassigns_flexible_job_for_constrained_job(hardware):
    jobs = [
        task("a-flexible", gpus=("GPU-first", "GPU-second")),
        task("b-constrained", gpus=("GPU-first",)),
    ]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 10, "concurrent_jobs": 2}


def test_gpu_matching_does_not_count_ineligible_second_gpu(hardware):
    jobs = [task("first", gpus=("GPU-second",)), task("second", gpus=("GPU-second",))]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 20, "concurrent_jobs": 1}


def test_unallocated_gpu_is_not_treated_as_cpu_job(hardware):
    jobs = [task("outside", gpus=("GPU-other",)), task("cpu")]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 20, "concurrent_jobs": 1}


def test_seven_jobs_admit_two_gpu_and_three_ram_limited_cpu_jobs(hardware):
    jobs = [task(f"gpu-{index}", gpus=("GPU-first", "GPU-second")) for index in range(3)]
    jobs += [task(f"cpu-{index}", memory=64) for index in range(4)]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 4, "concurrent_jobs": 5}


@pytest.mark.parametrize("fixed,automatic_count,expected", [
    (10, 1, {"threads": 10, "concurrent_jobs": 2}),
    (16, 2, {"threads": 2, "concurrent_jobs": 3}),
    (20, 2, {"threads": 1, "concurrent_jobs": 1}),
])
def test_fixed_override_reserves_cores_before_sharing(hardware, fixed, automatic_count, expected):
    jobs = [task("fixed", memory=1, threads=fixed)]
    jobs += [task(f"automatic-{index}") for index in range(automatic_count)]
    assert planner.automatic_threads(hardware, jobs) == expected


def test_multiple_fixed_overrides_leave_only_unreserved_cores(hardware):
    jobs = [task("fixed-1", threads=6), task("fixed-2", threads=8), task("auto")]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 6, "concurrent_jobs": 3}


def test_unschedulable_fixed_override_rejected(hardware):
    with pytest.raises(ValueError, match="exceed CPU budget"):
        planner.automatic_threads(hardware, [task("too-many", threads=21)])


def test_empty_pending_tasks_use_cpu_budget(hardware):
    assert planner.automatic_threads(hardware, []) == {"threads": 20, "concurrent_jobs": 0}


def test_no_ram_fitting_jobs_have_zero_concurrency(hardware):
    assert planner.automatic_threads(hardware, [task("too-large", memory=221)]) == {
        "threads": 20, "concurrent_jobs": 0,
    }


def test_cpu_budget_caps_minimum_one_thread_per_task(hardware):
    hardware["cpus"] = 3
    jobs = [task(str(index)) for index in range(7)]
    assert planner.automatic_threads(hardware, jobs) == {"threads": 1, "concurrent_jobs": 3}


def test_single_cpu_budget_is_not_reserved_again(hardware):
    hardware["cpus"] = 1
    assert planner.automatic_threads(hardware, [task("only")]) == {"threads": 1, "concurrent_jobs": 1}


def test_ram_sorted_greedy_policy_is_deterministic_and_does_not_mutate_inputs(hardware):
    jobs = [task("large", memory=180), task("small-2", memory=80), task("small-1", memory=80)]
    original_hardware = deepcopy(hardware)
    original_jobs = deepcopy(jobs)
    expected = {"threads": 10, "concurrent_jobs": 2}
    assert planner.automatic_threads(hardware, jobs) == expected
    assert planner.automatic_threads(hardware, list(reversed(jobs))) == expected
    assert jobs == original_jobs
    assert hardware == original_hardware


@pytest.mark.parametrize("value", [0, -1, 1.5, float("nan"), float("inf")])
def test_invalid_fixed_threads_rejected(hardware, value):
    with pytest.raises(ValueError, match="threads"):
        planner.automatic_threads(hardware, [task("invalid", threads=value)])


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_task_memory_rejected(hardware, value):
    with pytest.raises(ValueError, match="task RAM"):
        planner.automatic_threads(hardware, [task("invalid", memory=value)])


@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_pseudodca_memory_is_independent_of_unused_mcmc_chains(dtype):
    estimates = [planner.estimate_memory(819, 758, 65, "pseudoDCA", dtype, chains)
                 for chains in (1, 10000, 1000000)]
    assert estimates[0] == estimates[1] == estimates[2]


@pytest.mark.parametrize("model", ["bmDCA", "eaDCA", "edDCA"])
def test_mcmc_models_retain_chain_memory_estimates(model):
    small = planner.estimate_memory(819, 758, 65, model, "float32", 1)
    large = planner.estimate_memory(819, 758, 65, model, "float32", 1000000)
    assert large["gpu_memory_gib"] > small["gpu_memory_gib"]
    assert large["cpu_memory_gib"] > small["cpu_memory_gib"]


@pytest.mark.parametrize("dtype,scalar", [("float32", 4), ("float64", 8)])
def test_pseudodca_memory_accounts_for_dense_boolean_mask(dtype, scalar):
    sequences, sites, states = 819, 758, 65
    elements = sites * sites * states * states
    expected = (
        planner.MODEL_TENSORS["pseudoDCA"] * elements * scalar
        + elements
        + 2 * sequences * sites * states * scalar
        + sequences * sequences * scalar
    ) / planner.GIB + 2
    estimate = planner.estimate_memory(sequences, sites, states, "pseudoDCA", dtype, 64, margin=1)
    assert estimate["gpu_memory_gib"] == pytest.approx(expected)
    assert estimate["cpu_memory_gib"] == pytest.approx(expected)


def test_config_uses_provisional_one_thread_for_automatic_allocation(hardware, tmp_path):
    config = planner.make_config(config_args(tmp_path), hardware=hardware)
    assert config["threads"] == 1
    assert config["threads_explicit"] is False


def test_config_preserves_explicit_global_threads(hardware, tmp_path):
    config = planner.make_config(config_args(tmp_path, threads=8), hardware=hardware)
    assert config["threads"] == 8
    assert config["threads_explicit"] is True


@pytest.mark.parametrize("global_threads", [None, 8])
def test_task_override_remains_authoritative_with_automatic_or_explicit_global_threads(hardware, tmp_path, global_threads):
    overrides = tmp_path / "overrides.json"
    overrides.write_text(json.dumps({"GENE.protein": {"threads": 12}}))
    config = planner.make_config(config_args(tmp_path, threads=global_threads, resource_overrides=overrides), hardware=hardware)
    fasta = tmp_path / "GENE.fasta"
    msa = tmp_path / "GENE.a2m"
    mutations = tmp_path / "GENE.csv"
    fasta.write_text(">ORF\nATGGCT\n")
    msa.write_text(">ORF\nMA\n>homolog\nM-\n")
    mutations.write_text("mutant\nG4A\n")
    plan = planner.plan_task("GENE", "protein", fasta, msa, mutations, config)
    assert plan["threads"] == 12
