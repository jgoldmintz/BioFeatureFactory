"""Run protein CLIs through publication with external predictor fixtures only."""

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[2]
MUTATIONS = ("A4G", "T9C", "T9TTCT")
PROTEINS = {"NPM1": "MKSTNATEF", "PAM": "MKSTNATEFG"}
PIPELINES = {
    "NetSurfP3": "NetSurfP3/netsurfp3_pipeline.py",
    "NetMHC": "netMHC/netmhc_pipeline.py",
    "NetNglyc": "netNglyc/netnglyc_pipeline.py",
    "NetPhos": "netphos/netphos_pipeline.py",
}


def _write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _pkey(gene, token):
    return f"{gene}-{hashlib.sha1(token.encode()).hexdigest()[:12]}"


def _inputs(directory, extension):
    inputs = {}
    for gene in PROTEINS:
        gene_root = directory / gene
        nucleotide = "ATGAAATCTACTAATGCTACTGAATTT" + ("GGT" if gene == "PAM" else "") + "TAA"
        fasta = _write(gene_root / "fastas" / f"{gene}{extension}",
                       f">transcript\nCCCC{nucleotide}CCCC\n>ORF\n{nucleotide}\n")
        mutations = _write(gene_root / "mappings/mutations" / f"{gene}_mutations.csv",
                           "mutant\n" + "\n".join(MUTATIONS) + "\n")
        mapping = _write(gene_root / "mappings/aa" / f"{gene}_aa_mapping.csv",
                         "mutant,aamutant\nA4G,K2E\nT9C,S3S\nT9TTCT,S3SS\n")
        _write(gene_root / "mappings/transcript" / f"transcript_mapping_{gene}.csv",
               "mutant,transcript\nA4G,A8G\nT9C,T13C\nT9TTCT,T13TTCT\n")
        _write(gene_root / "mappings/pkey" / f"pkey_mapping_{gene}.csv",
               "pkey,mutant\n" + "".join(f"{_pkey(gene, token)},{token}\n" for token in MUTATIONS))
        _write(gene_root / "CodonMSA" / f"{gene}.codon.msa.fasta", ">ORF\nAAACCC\n>decoy\nGGGTTT\n")
        inputs[gene] = (fasta, mutations, mapping)
    return inputs


PREDICTOR_SUPPORT = '''
import json
import os
from pathlib import Path
import uuid

def records(path):
    sequences = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].split()[0]
            sequences[header] = ""
        elif line.strip():
            sequences[header] += line.strip()
    return sequences

def record(tool, sequences):
    destination = Path(os.environ["BFF_TEST_PREDICTOR_INPUTS"])
    destination.mkdir(parents=True, exist_ok=True)
    (destination / (uuid.uuid4().hex + ".json")).write_text(
        json.dumps({"tool": tool, "sequences": sequences}))
'''


def _external_predictors(directory):
    support = _write(directory / "predictor_support.py", PREDICTOR_SUPPORT)
    _write(directory / "sitecustomize.py",
           "import os\nif os.environ.get('BFF_TEST_NSP3_STUB') == '1':\n"
           "    import nsp3\n    from nsp3 import main, cli\n")
    native = '''
import sys
from predictor_support import records, record

tool = Path(sys.argv[0]).name
if tool == "signalp6" and "--version" in sys.argv:
    print("SignalP 6 fixture")
    raise SystemExit(0)
if tool == "signalp6":
    sequences = records(sys.argv[sys.argv.index("--fastafile") + 1])
    record(tool, sequences)
    destination = Path(sys.argv[sys.argv.index("--output_dir") + 1])
    destination.mkdir(parents=True, exist_ok=True)
    content = "# ID\\tPrediction\\tOTHER\\tSP(Sec/SPI)\\tCS Position\\n"
    content += "".join(header + "\\tSP\\t0.1\\t0.9\\tCS pos: 2-3. Pr: 0.8\\n" for header in sequences)
    (destination / "prediction_results.txt").write_text(content)
elif tool == "netnglyc":
    sequences = records(sys.argv[1])
    record(tool, sequences)
    assert Path(os.environ["SIGNALP6_RESULTS_DIR"], "prediction_results.txt").is_file()
    print("# Predictions for N-Glycosylation sites")
    for header, sequence in sequences.items():
        print("Name:", header, "Length:", len(sequence))
        print(sequence)
    print("SeqName Position Potential Jury N-Glyc")
    print("------------------------------------------------")
    for header, sequence in sequences.items():
        position = sequence.index("N") + 1
        print(header, position, sequence[position - 1:position + 3], 0.8, "9/9", "+")
elif tool == "ape":
    sequences = records(sys.argv[-1])
    record(tool, sequences)
    for header, sequence in sequences.items():
        for position, residue in enumerate(sequence, 1):
            if residue in "STY":
                print("#", header, position, residue, "XXX" + residue + "XXX", "0.800", "PKA", "YES")
elif tool == "netMHC":
    sequences = records(sys.argv[sys.argv.index("-f") + 1])
    record(tool, sequences)
    print("Pos HLA Peptide Core Identity Affinity(nM) %Rank BindLevel")
    for header, sequence in sequences.items():
        for offset in range(len(sequence) - 8):
            peptide = sequence[offset:offset + 9]
            print(offset, "HLA-A0201", peptide, peptide, header, "50", "0.2", "<= SB")
'''
    executables = {}
    for name in ("signalp6", "netnglyc", "ape", "netMHC"):
        executable = _write(directory / name, f"#!{sys.executable}\nfrom pathlib import Path\nimport os\n" + native)
        executable.chmod(0o755)
        executables[name] = executable
    _write(directory / "nsp3/__init__.py", "")
    _write(directory / "nsp3/cli.py", "def load_config(path):\n    return {'arch': {'args': {}}}\n")
    _write(directory / "nsp3/main.py", '''
from types import SimpleNamespace
import numpy as np
import torch
from predictor_support import records, record

class SecondaryFeatures:
    def __init__(self, model, checkpoint):
        self.model = model

    def __call__(self, fasta_path):
        sequences = records(fasta_path)
        record("nsp3", sequences)
        identifiers, residues, outputs = [], [], []
        for header, sequence in sequences.items():
            tensors = [np.full((1, len(sequence), width), 1 / width)
                       for width in (8, 3, 2, 1, 1, 1)]
            tensors[3] *= 0.2
            identifiers.append([header])
            residues.append([sequence])
            outputs.append(tensors)
        return identifiers, residues, outputs

def get_instance(module, name, config):
    return torch.nn.Identity()

module_arch = SimpleNamespace()
module_pred = SimpleNamespace(SecondaryFeatures=SecondaryFeatures)
''')
    return executables, support.parent


def _run(pipeline, run_root, input_path, companion, executables, module_path):
    output = run_root / "output"
    recording = run_root / "predictor_inputs"
    home = run_root / "home"
    home.mkdir(parents=True)
    command = [sys.executable, str(REPOSITORY / "biofeaturefactory" / PIPELINES[pipeline]),
               "-i", str(input_path), "-o", str(output)]
    if companion:
        command += ["-md" if pipeline in ("NetNglyc", "NetPhos") else "-m", str(companion)]
    if pipeline == "NetSurfP3":
        command += ["-M", "fixture.pth", "-c", "fixture.yml", "-bs", "2"]
    elif pipeline == "NetMHC":
        command += ["-nnp", str(executables["netMHC"]), "--max-workers", "1"]
    elif pipeline == "NetNglyc":
        command += ["-nnb", str(executables["netnglyc"]), "-snp", str(executables["signalp6"]),
                    "-cd", str(run_root / "cache"), "-w", "1"]
    else:
        command += ["-nap", str(executables["ape"]), "-nc"]
    environment = dict(os.environ, HOME=str(home), PYTHONDONTWRITEBYTECODE="1",
                       PYTHONPATH=os.pathsep.join([str(module_path), str(REPOSITORY)]),
                       BFF_TEST_NSP3_STUB="1" if pipeline == "NetSurfP3" else "0",
                       BFF_TEST_PREDICTOR_INPUTS=str(recording))
    result = subprocess.run(command, cwd=REPOSITORY, env=environment, text=True,
                            capture_output=True, timeout=120)
    _write(run_root / "command.json", json.dumps(command))
    _write(run_root / "stdout.log", result.stdout)
    _write(run_root / "stderr.log", result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    tables = {}
    for path in output.rglob("*.tsv"):
        with path.open() as handle:
            tables[str(path.relative_to(output))] = list(csv.DictReader(handle, delimiter="\t"))
    calls = {}
    for path in recording.glob("*.json"):
        call = json.loads(path.read_text())
        calls.setdefault(call["tool"], {}).update(call["sequences"])
    return tables, calls


@pytest.mark.parametrize("pipeline", PIPELINES)
@pytest.mark.parametrize("extension", (".fasta", ".fa"))
@pytest.mark.parametrize("detached", (False, True), ids=("tree_files", "standalone_files"))
def test_protein_cli_file_and_root_directory_publish_equivalent_complete_tables(tmp_path, pipeline, extension, detached):
    input_root = tmp_path / "input tree"
    inputs = _inputs(input_root, extension)
    executables, module_path = _external_predictors(tmp_path / "external tools")
    directory_tables, directory_calls = _run(
        pipeline, tmp_path / "directory", input_root, None, executables, module_path)
    assert len(directory_tables) == 6
    expected_calls = {}
    for gene, protein in PROTEINS.items():
        fasta, mutations, mapping = inputs[gene]
        companion = mapping if pipeline in ("NetNglyc", "NetPhos") else mutations
        if detached:
            standalone = tmp_path / "standalone files" / gene
            fasta = _write(standalone / fasta.name, fasta.read_text())
            companion = _write(standalone / companion.name, companion.read_text())
        file_tables, file_calls = _run(
            pipeline, tmp_path / gene, fasta, companion, executables, module_path)
        assert len(file_tables) == 3
        assert all(rows for rows in file_tables.values()), file_tables
        for name, rows in file_tables.items():
            assert directory_tables[name] == rows
        summary_name = (f"{gene}/NetSurfP3/{gene}.netsurfp3.summary.tsv"
                        if pipeline == "NetSurfP3" else f"{gene}/{pipeline}/{gene}.tsv")
        summaries = file_tables[summary_name]
        assert len(summaries) == 3
        assert {row["pkey"] for row in summaries} == {_pkey(gene, token) for token in MUTATIONS}
        assert all(not any(flag in row.get("qc_flags", "")
                           for flag in ("NOT_SCORED", "missing_mut", "FAILED", "REF_MISMATCH"))
                   for row in summaries), summaries
        sites_name = (f"{gene}/NetSurfP3/{gene}.netsurfp3.residues.tsv"
                      if pipeline == "NetSurfP3" else f"{gene}/{pipeline}/{gene}.sites.tsv")
        site_rows = file_tables[sites_name]
        if pipeline == "NetSurfP3":
            assert len(site_rows) == 6 * len(protein) + 1
            assert all(float(row["rsa"]) == 0.2 for row in site_rows)
        elif pipeline == "NetMHC":
            assert len(site_rows) == 4 * (len(protein) - 8) + 1
            assert all(float(row["affinity"]) == 50.0 and float(row["rank"]) == 0.2 for row in site_rows)
        elif pipeline == "NetNglyc":
            assert len(site_rows) == 4
            assert all(float(row["potential"]) == 0.8 and row["signalp_cleavage"] == "2" for row in site_rows)
        else:
            assert len(site_rows) == 19
            assert all(float(row["score"]) == 0.8 for row in site_rows)
        expected_sequences = {
            _pkey(gene, "A4G"): "ME" + protein[2:],
            _pkey(gene, "T9C"): protein,
            _pkey(gene, "T9TTCT"): protein[:3] + "S" + protein[3:],
        }
        expected_tools = {"NetSurfP3": {"nsp3"}, "NetMHC": {"netMHC"},
                          "NetNglyc": {"netnglyc", "signalp6"}, "NetPhos": {"ape"}}
        assert set(file_calls) == expected_tools[pipeline]
        for tool, sequences in file_calls.items():
            assert all(sequences.get(header) == sequence for header, sequence in expected_sequences.items()), sequences
            wildtypes = {header: sequence for header, sequence in sequences.items() if header not in expected_sequences}
            assert len(wildtypes) == 1 and set(wildtypes.values()) == {protein}, sequences
            expected_calls.setdefault(tool, {}).update(sequences)
    assert directory_calls == expected_calls
