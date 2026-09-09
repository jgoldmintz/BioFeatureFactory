"""SignalP failures and explicit executable propagation through real workers."""

import csv
import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

from biofeaturefactory.netNglyc import netnglyc_pipeline as pipeline


@pytest.fixture
def executables(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    tools_dir = tmp_path / "licensed tools"
    tools_dir.mkdir()
    signalp = tools_dir / "explicit-signalp"
    signalp.write_text(f'''#!{sys.executable}
import pathlib
import sys
if "--version" in sys.argv:
    print("SignalP 6 stub")
    raise SystemExit(0)
fasta = pathlib.Path(sys.argv[sys.argv.index("--fastafile") + 1])
headers = [line[1:].strip() for line in fasta.read_text().splitlines() if line.startswith(">")]
if any("BAD" in header for header in headers):
    print("simulated missing weights", file=sys.stderr)
    raise SystemExit(7)
output = pathlib.Path(sys.argv[sys.argv.index("--output_dir") + 1])
output.mkdir(parents=True, exist_ok=True)
rows = ["\\t".join([header, "OTHER", "0.9", "0.1", "0", "0", "0", "0", ""]) for header in headers]
(output / "prediction_results.txt").write_text("\\n".join(rows) + "\\n")
''')
    signalp.chmod(0o755)
    netnglyc = tools_dir / "explicit-netnglyc"
    netnglyc.write_text(f'''#!{sys.executable}
import os
import pathlib
import sys
assert (pathlib.Path(os.environ["SIGNALP6_RESULTS_DIR"]) / "prediction_results.txt").is_file()
headers = [line[1:].strip() for line in pathlib.Path(sys.argv[1]).read_text().splitlines() if line.startswith(">")]
print("# Predictions for N-Glycosylation sites")
for header in headers:
    print("Name:", header, "Length: 4")
    print("MKTE")
    print("----------------------------------------------------------------------")
print("SeqName Position Potential Jury N-Glyc")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
''')
    netnglyc.chmod(0o755)
    return signalp, netnglyc


def _processor(tmp_path, executables):
    signalp, netnglyc = executables
    return pipeline.RobustDockerNetNGlyc(
        cache_dir=str(tmp_path / "cache"), max_workers=2,
        native_bin=str(netnglyc), signalp_bin=str(signalp))


def test_valid_negative_and_signalp_cache_are_not_failure(tmp_path, executables):
    fasta = tmp_path / "good.fasta"
    fasta.write_text(">GOOD\nMKTE\n")
    with _processor(tmp_path, executables) as processor:
        output = tmp_path / "good.out"
        success, _, error = processor.process_single_fasta(str(fasta), str(output))
        assert success, error
        assert pipeline.parse_signalp_summary(output)["GOOD"] == {
            "has_signal": False, "probability": 0.1, "cleavage_site": None}
        executables[0].unlink()
        assert processor.process_single_fasta(str(fasta), str(tmp_path / "cached.out"))[0]


@pytest.mark.parametrize("mode", ["nonzero", "timeout", "missing", "malformed", "empty", "nan"])
def test_signalp_failure_never_returns_negative_or_publishes_cache(tmp_path, executables, monkeypatch, mode):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">SEQ\nMKTE\n")
    handler = pipeline.SignalP6Handler(str(tmp_path / "cache"), signalp_bin=str(executables[0]))

    def invoke(command, **kwargs):
        if mode == "timeout":
            raise subprocess.TimeoutExpired(command, 300)
        if mode == "nonzero":
            return subprocess.CompletedProcess(command, 7, "", "missing weights")
        output_dir = Path(command[command.index("--output_dir") + 1])
        output_dir.mkdir(parents=True, exist_ok=True)
        if mode != "missing":
            content = {"malformed": "SEQ\tOTHER\n", "empty": "# header\n",
                       "nan": "SEQ\tOTHER\t0\tnan\t0\t0\t0\t0\t\n"}[mode]
            (output_dir / "prediction_results.txt").write_text(content)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(pipeline.subprocess, "run", invoke)
    with pytest.raises(RuntimeError, match="SignalP 6 prediction failed"):
        handler.run_signalp6(str(fasta), str(tmp_path / "run"))
    assert not list((tmp_path / "cache").glob("*_sp6.json"))


def test_legacy_empty_cache_and_output_without_signalp_are_not_reused(tmp_path, executables):
    fasta = tmp_path / "bad.fasta"
    fasta.write_text(">BAD\nMKTE\n")
    with _processor(tmp_path, executables) as processor:
        checksum = hashlib.md5(fasta.read_bytes()).hexdigest()
        cache = Path(processor.cache_dir)
        (cache / f"{checksum}_sp6.json").write_text("{}")
        (cache / f"{checksum}_sp6_output").mkdir()
        (cache / f"{checksum[:16]}_netnglyc.out").write_text("Name: BAD\n")
        output = tmp_path / "bad.out"
        success, _, error = processor.process_single_fasta(str(fasta), str(output))
    assert not success
    assert "missing weights" in error
    assert not output.exists()


@pytest.mark.parametrize("cleavage", ["CS pos: 3-4. Pr: 0.8", "3-4", "3"])
def test_valid_positive_cleavage_is_preserved(tmp_path, executables, monkeypatch, cleavage):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">SEQ\nMKTE\n")
    handler = pipeline.SignalP6Handler(str(tmp_path / "cache"), signalp_bin=str(executables[0]))

    def invoke(command, **kwargs):
        output = Path(command[command.index("--output_dir") + 1])
        output.mkdir(parents=True, exist_ok=True)
        (output / "prediction_results.txt").write_text(
            f"SEQ\tSP\t0.1\t0.9\t0\t0\t0\t0\t{cleavage}\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(pipeline.subprocess, "run", invoke)
    results, directory = handler.run_signalp6(str(fasta), str(tmp_path / "run"))
    assert results["SEQ"] == {"has_signal": True, "probability": 0.9, "cleavage_site": 3}
    assert Path(directory, "prediction_results.txt").is_file()


@pytest.mark.parametrize("mode", ["single", "parallel", "batch", "sequential"])
def test_explicit_executable_survives_directory_worker_strategies(tmp_path, executables, mode):
    inputs = tmp_path / "inputs"
    outputs = tmp_path / "outputs"
    inputs.mkdir()
    (inputs / "GENE.fasta").write_text(">GOOD1\nMKTE\n>GOOD2\nMKTE\n")
    with _processor(tmp_path, executables) as processor:
        summary = processor.process_directory(str(inputs), str(outputs), processing_mode=mode)
    assert summary["success"] == 1, summary
    assert summary["failed"] == 0
    assert "GOOD1" in (outputs / "GENE-netnglyc.out").read_text()
    assert "GOOD2" in (outputs / "GENE-netnglyc.out").read_text()


@pytest.mark.parametrize("mode", ["parallel", "batch"])
def test_partial_sequence_failure_keeps_good_outputs_and_returns_failure(tmp_path, executables, mode):
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(">GOOD1\nMKTE\n>BAD\nMKTE\n>GOOD2\nMKTE\n")
    output = tmp_path / "mixed-netnglyc.out"
    with _processor(tmp_path, executables) as processor:
        if mode == "parallel":
            success, _, error = processor.process_parallel_docker(str(fasta), str(output), 2)
        else:
            success, _, error = processor.process_fasta_batched(str(fasta), str(output), batch_size=1)
    assert not success
    assert "failed" in error
    assert output.exists()
    assert "GOOD1" in output.read_text()
    assert "GOOD2" in output.read_text()
    assert "Name: BAD" not in output.read_text()


def test_cli_partial_gene_failure_publishes_tables_and_exits_nonzero(tmp_path, executables):
    input_root = tmp_path / "inputs"
    output_root = tmp_path / "outputs"
    for gene in ["GOOD", "BAD"]:
        fasta_dir = input_root / gene / "fastas"
        mapping_dir = input_root / gene / "mappings" / "aa"
        fasta_dir.mkdir(parents=True)
        mapping_dir.mkdir(parents=True)
        (fasta_dir / f"{gene}.fasta").write_text(">ORF\nATGAAAACCTAA\n")
        (mapping_dir / f"{gene}_aa_mapping.csv").write_text("mutant,aamutant\nA4G,K2E\n")
        pkey_dir = mapping_dir.parent / "pkey"
        pkey_dir.mkdir()
        (pkey_dir / f"pkey_mapping_{gene}.csv").write_text(
            f"pkey,mutant\n{pipeline.mint_pkey(gene, 'A4G')},A4G\n")
    command = [sys.executable, "-m", "biofeaturefactory.netNglyc.netnglyc_pipeline",
               "-i", str(input_root), "-md", str(input_root), "-o", str(output_root), "-w", "2",
               "-cd", str(tmp_path / "cache"), "-snp", str(executables[0]),
               "-nnb", str(executables[1])]
    result = subprocess.run(command, text=True, capture_output=True, timeout=90)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "execution(s) FAILED" in result.stderr
    for gene in ["GOOD", "BAD"]:
        path = output_root / gene / "NetNglyc" / f"{gene}.tsv"
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        assert ("NOT_SCORED" in rows[0]["qc_flags"]) == (gene == "BAD")
        if gene == "GOOD":
            assert rows[0]["wt_signalp_has_signal"] == "0"
            assert rows[0]["mut_signalp_has_signal"] == "0"
