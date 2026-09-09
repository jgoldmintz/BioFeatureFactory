"""Controller defaults and resource-aware thread sharing without model training."""

import json
import sys

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller
from biofeaturefactory.mutation_effects.bin.adabmdca_task import _file_record


def write_file(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def parse_inputs(tmp_path, monkeypatch, *extra, genes=2, gpu_count=2, backend="--adabmdca-only"):
    source = tmp_path / "source"
    for index in range(genes):
        gene = f"GENE{index}"
        write_file(source / gene / "fastas" / f"{gene}.fasta", ">ORF\nATGGCT\n")
        write_file(source / gene / "mappings" / "mutations" / f"{gene}_mutations.csv", "mutant\nG4A\n")
        write_file(source / gene / "MSA" / f"{gene}.msa.a2m", ">ORF\nMA\n>homolog\nMT\n")
    hardware = write_file(tmp_path / "hardware.json", json.dumps({
        "cpus": 20, "memory_gib": 220,
        "gpus": [{"uuid": f"GPU-{index}", "memory_gib": 80} for index in range(gpu_count)],
    }))
    monkeypatch.setattr(sys, "argv", [
        "mutEffects_controller", "--fasta", str(source), "--output", str(tmp_path / "output"),
        backend, "--resource-hardware", str(hardware), *map(str, extra),
    ])
    return controller.parse_args()


def prepare(args):
    genes = controller._resolve_genes(args.fasta)
    manifest = controller.build_manifest(genes, args)
    config, plans = controller.prepare_resource_plans(args, manifest, genes)
    return config, plans, manifest, genes


def test_default_two_gpu_jobs_share_twenty_usable_cpus(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch)
    config, plans, manifest, genes = prepare(args)
    assert args.threads is None
    assert args.adabmdca_model == "pseudoDCA"
    assert args.adabmdca_nepochs is None
    assert config["cpu_allocation"] == {"threads": 10, "concurrent_jobs": 2}
    assert {plan["threads"] for plan in plans} == {10}
    assert {plan["settings"]["model"] for plan in plans} == {"pseudoDCA"}
    assert all(plan["settings"]["nepochs"] is None for plan in plans)


@pytest.mark.parametrize("genes,gpu_count,threads,concurrency", [
    (2, 1, 20, 1), (3, 2, 10, 2), (5, 0, 4, 5),
])
def test_concurrency_not_total_gene_count_controls_share(tmp_path, monkeypatch, genes, gpu_count, threads, concurrency):
    args = parse_inputs(tmp_path, monkeypatch, genes=genes, gpu_count=gpu_count)
    config, plans, manifest, genes = prepare(args)
    assert config["cpu_allocation"] == {"threads": threads, "concurrent_jobs": concurrency}
    assert {plan["threads"] for plan in plans} == {threads}


def test_mixed_cpu_gpu_queues_share_cpu_budget(tmp_path, monkeypatch):
    overrides = write_file(tmp_path / "overrides.json", json.dumps({
        f"GENE{index}.protein": {"gpu_memory_gib": 100, "cpu_memory_gib": 40}
        for index in range(2, 5)
    }))
    args = parse_inputs(tmp_path, monkeypatch, "--resource-overrides", overrides, genes=5)
    config, plans, manifest, genes = prepare(args)
    assert [plan["device"] for plan in plans].count("cuda") == 2
    assert [plan["device"] for plan in plans].count("cpu") == 3
    assert config["cpu_allocation"] == {"threads": 4, "concurrent_jobs": 5}
    assert {plan["threads"] for plan in plans} == {4}


def test_explicit_threads_and_training_settings_are_preserved(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch, "--threads", 3, "--adabmdca-model", "bmDCA", "--adabmdca-nepochs", 17)
    config, plans, manifest, genes = prepare(args)
    assert config["threads_explicit"] is True
    assert config["cpu_allocation"] == {"threads": 3, "concurrent_jobs": None}
    assert {plan["threads"] for plan in plans} == {3}
    assert all(plan["settings"]["model"] == "bmDCA" and plan["settings"]["nepochs"] == 17 for plan in plans)


def test_thread_override_reserves_its_share(tmp_path, monkeypatch):
    overrides = write_file(tmp_path / "overrides.json", json.dumps({"GENE0.protein": {"threads": 12}}))
    args = parse_inputs(tmp_path, monkeypatch, "--resource-overrides", overrides)
    config, plans, manifest, genes = prepare(args)
    assert config["threads"] == 8
    assert {plan["gene"]: plan["threads"] for plan in plans} == {"GENE0": 12, "GENE1": 8}


@pytest.mark.parametrize("budget,threads", [(6, 3), (19, 9)])
def test_explicit_cpu_budget_is_shared_after_clamping(tmp_path, monkeypatch, budget, threads):
    args = parse_inputs(tmp_path, monkeypatch, "--resource-cpus", budget)
    config, plans, manifest, genes = prepare(args)
    assert config["hardware"]["cpus"] == budget
    assert config["threads"] == threads
    assert {plan["threads"] for plan in plans} == {threads}


def test_missing_alignment_counts_once_not_each_blocked_backend(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch, genes=1)
    args.run_evmutation = True
    (args.fasta / "GENE0" / "MSA" / "GENE0.msa.a2m").unlink()
    config, plans, manifest, genes = prepare(args)
    tasks = controller.pending_resource_tasks(args, manifest, genes, config, plans)
    assert tasks == [{"id": "GENE0.protein.msa", "memory_gib": 8, "eligible_gpu_uuids": [], "threads": None}]
    assert config["cpu_allocation"] == {"threads": 20, "concurrent_jobs": 1}


def test_ready_scoring_shares_with_missing_alignment_generation(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch)
    (args.fasta / "GENE1" / "MSA" / "GENE1.msa.a2m").unlink()
    config, plans, manifest, genes = prepare(args)
    tasks = controller.pending_resource_tasks(args, manifest, genes, config, plans)
    assert {task["id"] for task in tasks} == {"GENE0.protein.adabmdca", "GENE1.protein.msa"}
    assert config["cpu_allocation"] == {"threads": 10, "concurrent_jobs": 2}
    assert plans[0]["threads"] == 10


def test_evmutation_only_also_uses_automatic_share(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch, backend="--evmutation-only")
    config, plans, manifest, genes = prepare(args)
    assert plans == []
    assert config["cpu_allocation"] == {"threads": 10, "concurrent_jobs": 2}


def test_completed_tasks_do_not_consume_thread_share(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch)
    config, plans, manifest, genes = prepare(args)
    completed = plans[0]
    gene, side = completed["gene"], completed["side"]
    directory = args.output / gene / "adabmDCA"
    tsv = write_file(directory / f"{gene}.{side}.tsv", "pkey\tscore\nmutation\t1\n")
    params = write_file(args.output / f"adabmdca_{side}_params" / f"{gene}.{side}_adabm_params", "h 0 A 1\n")
    write_file(directory / f"{gene}.{side}.complete.json", json.dumps({
        "gene": gene, "side": side, "device": completed["device"],
        "success": True, "fingerprint": completed["fingerprint"],
        "files": {path.name: _file_record(path) for path in (tsv, params)},
    }))
    config, current, manifest, genes = prepare(args)
    assert config["cpu_allocation"] == {"threads": 20, "concurrent_jobs": 1}
    assert current[0]["device"] == "cached"
    assert current[1]["threads"] == 20
    assert current[1]["fingerprint"] == plans[1]["fingerprint"]


def test_nextflow_receives_resolved_share_for_all_task_families(tmp_path, monkeypatch):
    args = parse_inputs(tmp_path, monkeypatch)
    config, plans, manifest, genes = prepare(args)
    args.resource_config = config
    args.resource_config_path = tmp_path / "config.json"
    command = controller.build_nextflow_cmd(args, str(tmp_path / "manifest.json"))
    for option in ("--threads", "--msa_cpus", "--evmutation_cpus"):
        assert command[command.index(option) + 1] == "10"
    assert command[command.index("--adabmdca_model") + 1] == "pseudoDCA"
