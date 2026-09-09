"""Native-allocation planning tests without running plmc or loading dense models."""

from copy import deepcopy
import json
import struct
import tracemalloc

import pytest

from biofeaturefactory.mutation_effects.bin import plmc_resources as planner


def write_msa(tmp_path, contents):
    filename = tmp_path / "alignment.fasta"
    filename.write_text(contents)
    return filename


def write_params(tmp_path, sites=2, states=20, sequences=2):
    filename = tmp_path / "model.params"
    pairs = sites * (sites - 1) // 2
    size = 40 + states + 4 * sequences + 5 * sites + 8 * sites * states + 8 * pairs * states * states
    header = struct.pack("=5i5f", sites, states, sequences, 0, 500, 0.2, 0.01, 16.2, 0.0, float(sequences))
    filename.write_bytes(header + bytes(size - len(header)))
    return filename


@pytest.fixture
def config():
    return {"threads": 4, "memory_margin": 1.15,
            "hardware": {"cpus": 20, "memory_gib": 220},
            "evmutation": {"memory_gib": None}}


def test_protein_focus_counts_actual_alignment_columns_not_orf_or_first_row(tmp_path):
    msa = write_msa(tmp_path, ">homolog\nMAAA--\n>GENE/12-17 description\nM-A--A\n")
    shape = planner.alignment_shape(msa, "protein", "GENE")
    assert shape["sites"] == 3
    assert shape["alignment_sites"] == 6
    assert shape["states"] == 20
    assert shape["focus_header"] == "GENE/12-17 description"
    assert shape["focus_mode"] == "prefix_non_gap"


def test_missing_focus_keeps_full_width_including_gaps(tmp_path):
    msa = write_msa(tmp_path, ">ORF\nMA--AA\n>homolog\nM-A-AA\n")
    shape = planner.alignment_shape(msa, "protein", "GENE")
    assert shape["sites"] == 6
    assert shape["focus_found"] is False
    assert shape["focus_mode"] == "missing_focus_full_alignment"


def test_first_prefix_focus_wins_for_protein(tmp_path):
    msa = write_msa(tmp_path, ">GENE-first\nM---\n>GENE-second\nMAAA\n")
    assert planner.alignment_shape(msa, "protein", "GENE")["sites"] == 1


def test_explicit_custom_protein_alphabet_does_not_strip_lowercase_columns(tmp_path):
    msa = write_msa(tmp_path, ">GENE\nMAaa--\n>homolog\nMAAA--\n")
    shape = planner.alignment_shape(msa, "protein", "GENE")
    assert shape["sites"] == 4
    assert shape["valid_sequences"] == 1
    assert shape["sequences"] == 2


def test_codon_encoding_preserves_lowercase_encoded_symbols_and_coerces_invalid_triplets(tmp_path):
    msa = write_msa(tmp_path, ">ORF\nAAACGG---NNNATG\n>homolog\nAAACGGAAAAAAATG\n")
    shape = planner.alignment_shape(msa, "codon", "ORF")
    assert planner.CODON_TO_CHAR["CGG"].islower()
    assert shape["sites"] == 3
    assert shape["alignment_sites"] == 5
    assert shape["states"] == 64


def test_codon_duplicate_header_matches_shared_encoder_last_record_wins(tmp_path):
    msa = write_msa(tmp_path, ">ORF\nAAA------\n>homolog\nAAAAAAAAA\n>ORF\nAAAAAAAAA\n")
    shape = planner.alignment_shape(msa, "codon", "ORF")
    assert shape["sites"] == 3
    assert shape["sequences"] == 3


def test_codon_header_trailing_spaces_do_not_merge_distinct_encoder_records(tmp_path):
    msa = write_msa(tmp_path, ">ORF\nAAAAAAAAA\n>ORF \nAAA------\n")
    shape = planner.alignment_shape(msa, "codon", "ORF")
    assert shape["sites"] == 3
    assert shape["focus_header"] == "ORF"


@pytest.mark.parametrize("contents,side,error", [
    ("", "protein", "Empty alignment"),
    ("MA\n", "protein", "Missing FASTA"),
    (">GENE\n", "protein", "Empty alignment record"),
    (">GENE\nMA\n>other\nMAA\n", "protein", "Unequal alignment"),
    (">GENE\nZZZ\n", "protein", "no valid"),
    (">GENE\n---\n", "protein", "no modeled"),
    (">ORF\nATGA\n", "codon", "divisible by three"),
])
def test_invalid_alignments_fail_before_training(tmp_path, contents, side, error):
    msa = write_msa(tmp_path, contents)
    with pytest.raises(ValueError, match=error):
        planner.alignment_shape(msa, side, "ORF" if side == "codon" else "GENE")


def test_native_optimizer_counts_nineteen_parameter_arrays_and_thread_blocks():
    estimate = planner.estimate_memory(10, 40, 20, workspace_threads=5, margin=1)
    parameters = 40 * 20 + (40 * 39 // 2) * 20 * 20
    assert estimate["parameter_count"] == parameters
    assert estimate["optimizer_bytes"] == 18 * parameters * 8
    assert estimate["marginal_bytes"] == parameters * 8
    assert estimate["thread_workspace_bytes"] == 5 * (2 * 40 * 20 * 20 + 2 * 20) * 8
    assert estimate["lbfgs_history"] == 6
    assert estimate["native_scalar_bytes"] == 8
    assert estimate["assumed_native_precision"] == "float64"


def test_scoring_counts_four_full_float64_arrays_even_for_prebuilt_model():
    estimate = planner.estimate_memory(10, 40, 64, prebuilt=True, margin=1)
    assert estimate["scoring_dense_bytes"] == 4 * 40 * 40 * 64 * 64 * 8
    assert estimate["training_memory_gib"] == 0
    assert estimate["optimizer_bytes"] == 0
    assert estimate["estimated_memory_gib"] == estimate["scoring_memory_gib"]


def test_memory_uses_larger_phase_not_sum_of_sequential_training_and_scoring():
    estimate = planner.estimate_memory(20, 300, 64)
    assert estimate["estimated_memory_gib"] == max(estimate["training_memory_gib"], estimate["scoring_memory_gib"])
    assert estimate["estimated_memory_gib"] < estimate["training_memory_gib"] + estimate["scoring_memory_gib"]


def test_sequence_count_and_raw_msa_size_contribute_to_input_allowance():
    small = planner.estimate_memory(10, 40, 20, alignment_sites=50, msa_bytes=500)
    large = planner.estimate_memory(10000, 40, 20, alignment_sites=50, msa_bytes=500000)
    assert large["input_allowance_bytes"] > small["input_allowance_bytes"]
    assert large["estimated_memory_gib"] > small["estimated_memory_gib"]
    assert large["parameter_count"] == small["parameter_count"]


def test_native_parameter_integer_overflow_is_rejected_not_reported_schedulable():
    with pytest.raises(ValueError, match="native signed 32-bit parameter count"):
        planner.estimate_memory(819, 1025, 64)


def test_prebuilt_scoring_does_not_apply_native_training_integer_limit():
    estimate = planner.estimate_memory(819, 1025, 64, prebuilt=True)
    assert estimate["training_memory_gib"] == 0
    assert estimate["scoring_memory_gib"] > 100


def test_native_alignment_integer_overflow_is_rejected():
    with pytest.raises(ValueError, match="sequence indexing overflows"):
        planner.estimate_memory(30000000, 20, 20, alignment_sites=100)


def test_prebuilt_header_dimensions_drive_scoring_not_current_msa(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\nMA\n>homolog\nM-\n")
    params = write_params(tmp_path, sites=7, states=65)
    plan = planner.plan_evmutation_task("GENE", "protein", msa, config, params=params)
    assert plan["sites"] == 7
    assert plan["states"] == 65
    assert plan["prebuilt"] is True
    assert plan["training_memory_gib"] == 0
    assert plan["scoring_dense_bytes"] == 4 * 7 * 7 * 65 * 65 * 8


@pytest.mark.parametrize("contents,error", [
    (b"", "header"),
    (b"stub params", "header"),
    (struct.pack("=5i", -1, 20, 2, 0, 500), "dimensions"),
    (struct.pack("=5i", 2, 0, 2, 0, 500), "dimensions"),
    (struct.pack("=5i", 2, 20, -1, 0, 500), "dimensions"),
    (struct.pack("=5i", 2, 20, 2, 0, 500), "payload"),
])
def test_invalid_or_partial_prebuilt_params_rejected(tmp_path, contents, error):
    params = tmp_path / "partial.params"
    params.write_bytes(contents)
    with pytest.raises(ValueError, match=error):
        planner.params_shape(params)


def test_missing_prebuilt_params_not_silently_treated_as_training(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    with pytest.raises(FileNotFoundError):
        planner.plan_evmutation_task("GENE", "protein", msa, config, params=tmp_path / "absent")


def test_full_cpu_budget_reserves_workspace_without_changing_task_threads(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\n" + "A" * 40 + "\n")
    config["threads"] = 1
    original = deepcopy(config)
    plan = planner.plan_evmutation_task("GENE", "protein", msa, config)
    assert plan["threads"] == 1
    assert plan["workspace_threads"] == 20
    assert config == original


@pytest.mark.parametrize("floor", [0.01, 10])
def test_explicit_ram_is_a_minimum_not_estimator_bypass(tmp_path, config, floor):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    config["evmutation"]["memory_gib"] = floor
    plan = planner.plan_evmutation_task("GENE", "protein", msa, config)
    assert plan["memory_gib"] == max(floor, plan["estimated_memory_gib"])


def test_ram_oversize_job_fails_preflight(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\n" + "A" * 1200 + "\n")
    config["hardware"]["memory_gib"] = 10
    with pytest.raises(ValueError, match="EVmutation unschedulable"):
        planner.plan_evmutation_task("GENE", "protein", msa, config)


@pytest.mark.parametrize("margin", [0, 0.9, float("nan"), float("inf")])
def test_invalid_memory_margin_rejected(margin):
    with pytest.raises(ValueError, match="margin"):
        planner.estimate_memory(10, 20, 20, margin=margin)


def test_planner_cli_without_config_uses_explicit_fallback_threads_and_ram_floor(tmp_path):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    output = tmp_path / "plan.json"
    planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                  "--threads", "3", "--memory-gib", "4", "--output", str(output)])
    plan = json.loads(output.read_text())
    assert plan["memory_gib"] == 4
    assert plan["threads"] == 3


def test_planner_cli_config_takes_precedence_over_legacy_fallbacks(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                  "--config", str(config_file), "--threads", "9", "--memory-gib", "80", "--output", str(output)])
    plan = json.loads(output.read_text())
    assert plan["threads"] == 4
    assert plan["memory_gib"] == 80


def test_cli_ram_floor_cannot_lower_configured_minimum(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    config["evmutation"]["memory_gib"] = 20
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                  "--config", str(config_file), "--memory-gib", "2", "--output", str(output)])
    plan = json.loads(output.read_text())
    assert plan["memory_gib"] == 20


@pytest.mark.parametrize("cli_floor,config_floor", [(0, 20), (-1, 20), (float("nan"), 20), (20, 0), (20, -1)])
def test_cli_validates_both_ram_floors_when_config_is_supplied(tmp_path, config, cli_floor, config_floor):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    config["evmutation"]["memory_gib"] = config_floor
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    with pytest.raises(SystemExit) as failure:
        planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                      "--config", str(config_file), "--memory-gib", str(cli_floor), "--output", str(tmp_path / "plan.json")])
    assert failure.value.code == 1


def test_planner_cli_failure_does_not_write_success_plan(tmp_path):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    output = tmp_path / "plan.json"
    with pytest.raises(SystemExit) as failure:
        planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                      "--params", str(tmp_path / "absent"), "--output", str(output)])
    assert failure.value.code == 1
    assert not output.exists()


def test_alignment_inspection_memory_is_bounded_by_record_not_sequence_count(tmp_path):
    msa = tmp_path / "many.fasta"
    with msa.open("w") as handle:
        for index in range(10000):
            handle.write(f">sequence-{index}\n{'A' * 100}\n")
    tracemalloc.start()
    try:
        shape = planner.alignment_shape(msa, "protein", "missing")
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert shape["sequences"] == 10000
    assert peak < 300000


@pytest.mark.parametrize("defer_errors", [False, True])
@pytest.mark.parametrize("failure_kind", ["oversize", "invalid_msa", "missing_msa", "missing_params"])
def test_cli_defers_only_task_planning_errors(tmp_path, config, capsys, defer_errors, failure_kind):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    extra_args = []
    if failure_kind == "oversize":
        config["hardware"]["memory_gib"] = 0.1
    elif failure_kind == "invalid_msa":
        msa.write_text(">GENE\nMA\n>other\nM\n")
    elif failure_kind == "missing_msa":
        msa = tmp_path / "absent.fasta"
    else:
        extra_args = ["--params", str(tmp_path / "absent.params")]
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    arguments = ["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                 "--config", str(config_file), "--output", str(output), *extra_args]
    if not defer_errors:
        with pytest.raises(SystemExit) as failure:
            planner.main(arguments)
        assert failure.value.code == 1
        assert not output.exists()
        assert "resource planning failed" in capsys.readouterr().err
        return
    planner.main([*arguments, "--defer-errors"])
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""
    plan = json.loads(output.read_text())
    assert plan["device"] == "blocked"
    assert "memory_gib" not in plan
    assert plan["resource_error"] == {
        "gene": "GENE", "side": "protein", "backend": "evmutation",
        "message": plan["resource_error"]["message"], "stage": "resource_planning",
    }
    assert plan["resource_error"]["message"]


def test_cli_deferred_mode_keeps_valid_plan_unchanged(tmp_path, config):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config))
    output = tmp_path / "plan.json"
    planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                  "--config", str(config_file), "--output", str(output), "--defer-errors"])
    assert json.loads(output.read_text()) == planner.plan_evmutation_task("GENE", "protein", msa, config)


@pytest.mark.parametrize("failure_kind", ["config", "threads", "memory"])
def test_cli_deferred_mode_does_not_hide_global_input_errors(tmp_path, failure_kind):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    output = tmp_path / "plan.json"
    extra_args = {
        "config": ["--config", str(tmp_path / "missing.json")],
        "threads": ["--threads", "0"],
        "memory": ["--memory-gib", "0"],
    }[failure_kind]
    with pytest.raises(SystemExit) as failure:
        planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                      "--output", str(output), "--defer-errors", *extra_args])
    assert failure.value.code == 1
    assert not output.exists()


@pytest.mark.parametrize("blocked", [False, True])
def test_cli_deferred_mode_does_not_hide_output_write_errors(tmp_path, blocked):
    msa = write_msa(tmp_path, ">GENE\nMA\n")
    if blocked:
        msa = tmp_path / "missing.fasta"
    with pytest.raises(SystemExit) as failure:
        planner.main(["plan", "--gene", "GENE", "--side", "protein", "--msa", str(msa),
                      "--output", str(tmp_path), "--defer-errors"])
    assert failure.value.code == 1
