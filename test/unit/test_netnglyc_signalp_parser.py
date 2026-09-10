"""SignalP table formats and authoritative raw-output cache parsing."""

import json
import runpy
import subprocess
from pathlib import Path

import pytest

from biofeaturefactory.netNglyc import netnglyc_pipeline as pipeline


EUKARYA_COLUMNS = ("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "CS Position")
OTHER_COLUMNS = (
    "ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)",
    "TAT(Tat/SPI)", "TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS Position",
)
EXPECTED = {
    "SIGNAL": {"has_signal": True, "probability": 0.9, "cleavage_site": 3},
    "NEGATIVE": {"has_signal": False, "probability": 0.1, "cleavage_site": None},
}


def prediction_table(organism="eukarya", header=True, columns=None):
    selected_columns = columns or (EUKARYA_COLUMNS if organism == "eukarya" else OTHER_COLUMNS)
    rows = []
    for identifier, prediction, other_prob, signal_prob, cleavage in (
        ("SIGNAL", "SP", "0.1", "0.9", "CS pos: 3-4. Pr: 0.8500"),
        ("NEGATIVE", "OTHER", "0.9", "0.1", ""),
    ):
        values = {name: "0.0" for name in OTHER_COLUMNS}
        values.update({
            "ID": identifier, "Prediction": prediction, "OTHER": other_prob,
            "SP(Sec/SPI)": signal_prob, "CS Position": cleavage,
        })
        rows.append("\t".join(values[name] for name in selected_columns))
    prefix = (
        f"# SignalP-6.0\tOrganism: {organism.capitalize()}\tTimestamp: 20260909120000\n"
        "# " + "\t".join(selected_columns) + "\n"
    ) if header else ""
    return prefix + "\n".join(rows) + "\n"


def write_predictions(tmp_path, content):
    destination = tmp_path / "prediction_results.txt"
    destination.write_text(content)
    return destination


@pytest.mark.parametrize("organism", ("eukarya", "other"))
@pytest.mark.parametrize("header", (False, True))
def test_supported_formats_keep_signal_and_negative_rows(tmp_path, organism, header):
    source = write_predictions(tmp_path, prediction_table(organism, header))
    assert pipeline.parse_signalp_predictions(source, expected_ids=set(EXPECTED)) == EXPECTED


@pytest.mark.parametrize("columns", (
    ("CS Position", "OTHER", "ID", "SP(Sec/SPI)", "Prediction"),
    tuple(reversed(OTHER_COLUMNS)),
))
def test_header_names_control_reordered_columns(tmp_path, columns):
    source = write_predictions(tmp_path, prediction_table(columns=columns))
    assert pipeline.parse_signalp_predictions(source, expected_ids=set(EXPECTED)) == EXPECTED


@pytest.mark.parametrize("cleavage", ("CS pos: 3-4. Pr: 0.8500", "3-4", "3"))
@pytest.mark.parametrize("column_count", (5, 9))
def test_legacy_cleavage_spellings_remain_supported(tmp_path, cleavage, column_count):
    fields = ["SIGNAL", "SP", "0.1", "0.9", *(["0"] * (column_count - 5)), cleavage]
    source = write_predictions(tmp_path, "\t".join(fields) + "\n")
    assert pipeline.parse_signalp_predictions(source) == {"SIGNAL": EXPECTED["SIGNAL"]}


@pytest.mark.parametrize("probability", ("nan", "NaN", "inf", "-inf", "-0.01", "1.01", "broken", ""))
def test_invalid_signal_probabilities_are_rejected(tmp_path, probability):
    content = prediction_table().replace("SIGNAL\tSP\t0.1\t0.9\t", f"SIGNAL\tSP\t0.1\t{probability}\t")
    source = write_predictions(tmp_path, content)
    with pytest.raises(ValueError):
        pipeline.parse_signalp_predictions(source)


@pytest.mark.parametrize("content", (
    "SIGNAL\tSP\n",
    "NEGATIVE\tOTHER\t0.9\t0.1\n",
    "SIGNAL\tSP\t0.1\t0.9\t0\tCS pos: 3-4. Pr: 0.85\n",
    "SIGNAL\tSP\t0.1\t0.9\t\n",
    "SIGNAL\tSP\t0.1\t0.9\tgarbage\n",
    "\tSP\t0.1\t0.9\t3\n",
    "SIGNAL\tUNKNOWN\t0.1\t0.9\t3\n",
    "# ID\tPrediction\tOTHER\tCS Position\nSIGNAL\tSP\t0.1\t3\n",
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nSIGNAL\tSP\t0.1\t0.9\n",
))
def test_malformed_rows_are_rejected(tmp_path, content):
    source = write_predictions(tmp_path, content)
    with pytest.raises(ValueError):
        pipeline.parse_signalp_predictions(source)


def test_duplicate_sequence_ids_are_rejected(tmp_path):
    content = prediction_table() + "SIGNAL\tSP\t0.1\t0.9\t3\n"
    source = write_predictions(tmp_path, content)
    with pytest.raises(ValueError):
        pipeline.parse_signalp_predictions(source)


@pytest.mark.parametrize("expected_ids", ({"SIGNAL", "NEGATIVE", "MISSING"}, {"SIGNAL"}))
def test_expected_ids_must_match_the_complete_output(tmp_path, expected_ids):
    source = write_predictions(tmp_path, prediction_table())
    with pytest.raises(ValueError):
        pipeline.parse_signalp_predictions(source, expected_ids=expected_ids)


@pytest.mark.parametrize("content", ("", "# SignalP-6.0\tOrganism: Eukarya\n"))
def test_empty_output_cannot_satisfy_expected_inputs(tmp_path, content):
    source = write_predictions(tmp_path, content)
    with pytest.raises(ValueError):
        pipeline.parse_signalp_predictions(source, expected_ids={"SIGNAL"})


@pytest.mark.parametrize("organism", ("eukarya", "other"))
def test_cache_loader_recovers_cleavage_from_raw_prediction_tables(tmp_path, monkeypatch, organism):
    monkeypatch.setenv("HOME", str(tmp_path / "unused-home"))
    output = tmp_path / "cache" / "fixture_sp6_output"
    output.mkdir(parents=True)
    write_predictions(output, prediction_table(organism))
    (output.parent / "fixture_sp6.json").write_text(json.dumps({
        "SIGNAL": {"has_signal": True, "probability": 0.9, "cleavage_site": None},
    }))
    assert pipeline.load_signalp_cache(output.parent) == EXPECTED


def test_cache_loader_does_not_publish_partial_rows_from_malformed_file(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "unused-home"))
    broken = tmp_path / "cache" / "broken_sp6_output"
    valid = tmp_path / "cache" / "valid_sp6_output"
    broken.mkdir(parents=True)
    valid.mkdir()
    write_predictions(broken, prediction_table().replace("NEGATIVE\tOTHER\t0.9\t0.1\t", "NEGATIVE\tOTHER"))
    write_predictions(valid, "VALID\tOTHER\t0.9\t0.1\t\n")
    assert pipeline.load_signalp_cache(broken.parent) == {
        "VALID": {"has_signal": False, "probability": 0.1, "cleavage_site": None},
    }


@pytest.fixture
def handler(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "unused-home"))
    monkeypatch.setattr(pipeline, "resolve_signalp6_path", lambda path: "/unused/signalp6")
    monkeypatch.setattr(pipeline.SignalP6Handler, "check_signalp6_available", lambda self: True)
    return pipeline.SignalP6Handler(str(tmp_path / "cache"), signalp_bin="/unused/signalp6")


def make_input(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">SIGNAL\nMKTE\n>NEGATIVE\nMKTE\n")
    return fasta


def seed_cache(handler, fasta, raw_content, json_content):
    checksum = handler.get_cache_key(str(fasta))
    cache = Path(handler.cache_dir)
    raw_directory = cache / f"{checksum}_sp6_output"
    raw_directory.mkdir()
    write_predictions(raw_directory, raw_content)
    (cache / f"{checksum}_sp6.json").write_text(json_content)
    return raw_directory


@pytest.mark.parametrize("organism", ("eukarya", "other"))
@pytest.mark.parametrize("json_content", ("{bad json", json.dumps({
    "SIGNAL": {"has_signal": True, "probability": 0.5, "cleavage_site": None},
    "NEGATIVE": {"has_signal": False, "probability": 0.5, "cleavage_site": None},
})))
def test_runtime_cache_reparses_raw_predictions_without_inference(
    tmp_path, monkeypatch, handler, organism, json_content,
):
    fasta = make_input(tmp_path)
    raw_directory = seed_cache(handler, fasta, prediction_table(organism), json_content)
    handler.signalp6_available = False
    monkeypatch.setattr(pipeline.subprocess, "run", lambda *args, **kwargs: pytest.fail("Unexpected inference"))
    results, directory = handler.run_signalp6(str(fasta))
    assert results == EXPECTED
    assert Path(directory) == raw_directory


@pytest.mark.parametrize("raw_content", (
    "SIGNAL\tSP\n",
    "SIGNAL\tSP\t0.1\t0.9\t3\n",
))
def test_invalid_raw_cache_cannot_be_hidden_by_valid_json(tmp_path, monkeypatch, handler, raw_content):
    fasta = make_input(tmp_path)
    seed_cache(handler, fasta, raw_content, json.dumps(EXPECTED))
    calls = []

    def invoke(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 7, "", "simulated inference failure")

    monkeypatch.setattr(pipeline.subprocess, "run", invoke)
    with pytest.raises(RuntimeError, match="SignalP 6 prediction failed"):
        handler.run_signalp6(str(fasta), str(tmp_path / "fresh"))
    assert len(calls) == 1


def test_fresh_eukarya_prediction_then_cached_reparse(tmp_path, monkeypatch, handler):
    fasta = make_input(tmp_path)
    calls = []

    def invoke(command, **kwargs):
        calls.append(command)
        assert command[command.index("--organism") + 1] == "eukarya"
        assert command[command.index("--mode") + 1] == "fast"
        output = Path(command[command.index("--output_dir") + 1])
        output.mkdir(parents=True, exist_ok=True)
        write_predictions(output, prediction_table())
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(pipeline.subprocess, "run", invoke)
    fresh, _ = handler.run_signalp6(str(fasta), str(tmp_path / "fresh"))
    assert fresh == EXPECTED
    cached, _ = handler.run_signalp6(str(fasta))
    assert cached == fresh
    assert len(calls) == 1


@pytest.mark.parametrize("organism", ("eukarya", "other"))
def test_legacy_adapter_already_accepts_both_table_formats(tmp_path, monkeypatch, capsys, organism):
    write_predictions(tmp_path, prediction_table(organism))
    monkeypatch.setenv("SIGNALP6_RESULTS_DIR", str(tmp_path))
    adapter_path = Path(pipeline.__file__).parent / "bin" / "signalp6_adapter"
    adapter = runpy.run_path(str(adapter_path), run_name="signalp_adapter_test")
    adapter["main"]()
    rows = [line.split() for line in capsys.readouterr().out.splitlines()]
    assert [(row[0], len(row), row[13]) for row in rows] == [
        ("SIGNAL", 14, "Y"), ("NEGATIVE", 14, "N"),
    ]
