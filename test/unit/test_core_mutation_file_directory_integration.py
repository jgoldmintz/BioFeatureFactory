"""Real BFF CLI paths; only search/alignment executables are synthetic.

Optional licensed/upstream dependencies are supplied via PYTHONPATH. The
controller check additionally requires BFF_TEST_NEXTFLOW=1 and runs an unedited
copy of the BFF controller/workflow in a temporary directory.
"""

import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import pickle
import shutil
import struct
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from biofeaturefactory.lib.utility import mint_pkey
from biofeaturefactory.mutation_effects.bin.codon_encoding import CODON_ALPHABET, CODON_TO_CHAR


REPO = Path(__file__).resolve().parents[2]
PACKAGE = REPO / "biofeaturefactory"
GENES = ("NPM1", "PAM")
ORF = "ATGGCTGAATTCCAGTAA"
PROTEIN = "MAEFQ"
MUTATIONS = ("G4A", "T6C", "C13T")
PROTEIN_ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"


def write_file(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def executable(path, source):
    write_file(path, f"#!{sys.executable}\n" + textwrap.dedent(source))
    path.chmod(0o755)
    return path


def run_cli(script, arguments, workdir, environment=None):
    runtime_env = dict(os.environ)
    runtime_env.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": os.pathsep.join(filter(None, [str(REPO), os.environ.get("PYTHONPATH")])),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    runtime_env.update(environment or {})
    result = subprocess.run(
        [sys.executable, str(script), *map(str, arguments)], cwd=workdir,
        env=runtime_env, text=True, capture_output=True, timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Error in rare codon analysis" not in result.stderr, result.stderr
    return result


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture
def gene_tree(tmp_path):
    root = tmp_path / "input root with spaces"
    files = {}
    for gene in GENES:
        gene_root = root / gene
        files[gene] = {
            "fasta": write_file(gene_root / "fastas" / f"{gene}.fasta", f">ORF\n{ORF}\n>transcript\nCCC{ORF}CCC\n>genomic\nAAA{ORF}AAA\n"),
            "mutations": write_file(gene_root / "mappings" / "mutations" / f"{gene}_mutations.csv", "mutant\n" + "\n".join(MUTATIONS) + "\n"),
            "msa": write_file(gene_root / "MSA" / f"{gene}.msa.a2m", f">{gene}\n{PROTEIN}\n>homolog\nMTEFQ\n"),
            "codon_msa": write_file(gene_root / "CodonMSA" / f"{gene}.codon.msa.fasta", ">ORF\nATGGCTGAATTCCAG\n>homolog\nATGGCCGAATTCCAG\n"),
        }
        write_file(gene_root / "mappings" / "aa" / f"{gene}_aa_mapping.csv", "pkey,mutant\nwrong-key,NOT_A_MUTATION\n")
    return root, files


def assert_rows_equal(file_output, directory_output, suffix, expected_tokens=MUTATIONS):
    for gene in GENES:
        file_rows = rows(file_output / gene / suffix.format(gene=gene))
        directory_rows = rows(directory_output / gene / suffix.format(gene=gene))
        assert file_rows == directory_rows
        assert {row["pkey"] for row in file_rows} == {mint_pkey(gene, token) for token in expected_tokens}
        assert len(file_rows) == len(expected_tokens)
    return file_rows


def detached_files(files, tmp_path):
    detached = tmp_path / "detached inputs"
    detached.mkdir(exist_ok=True)
    for paths in files.values():
        copied = {}
        for artifact, source in paths.items():
            target = detached / source.name
            shutil.copy2(source, target)
            copied[artifact] = target
        yield copied


def test_codon_usage_cli_file_directory_predictions(gene_tree, tmp_path):
    root, files = gene_tree
    script = PACKAGE / "codon_usage" / "codon_usage_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    run_cli(script, ["-f", root, "-o", directory_output], tmp_path)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-f", paths["fasta"], "-m", paths["mutations"], "-o", file_output], tmp_path)
    assert_rows_equal(file_output, directory_output, "CodonUsage/{gene}.codon_usage.tsv")
    for gene in GENES:
        predictions = {row["pkey"]: row for row in rows(file_output / gene / "CodonUsage" / f"{gene}.codon_usage.tsv")}
        for token, expected_codon in (("G4A", "ACT"), ("T6C", "GCC"), ("C13T", "TAG")):
            prediction = predictions[mint_pkey(gene, token)]
            assert prediction["codon_mut"] == expected_codon
            assert prediction["codon_wt"] == ("CAG" if token == "C13T" else "GCT")
            assert prediction["RSCU_wt"] != ""


def test_protein_msa_cli_file_directory_processing(gene_tree, tmp_path):
    root, files = gene_tree
    invocation_log = tmp_path / "jackhmmer-inputs.jsonl"
    jackhmmer = executable(tmp_path / "tool bin" / "jackhmmer", '''
        import json
        import os
        from pathlib import Path
        import sys

        arguments = sys.argv[1:]
        query = Path(arguments[-2]).read_text()
        records = query.splitlines()
        assert records[1] == 'MAEFQ', query
        with open(os.environ['BFF_SEARCH_LOG'], 'a') as handle:
            handle.write(json.dumps(query) + '\\n')
        Path(arguments[arguments.index('-A') + 1]).write_text(
            '# STOCKHOLM 1.0\\n' + records[0][1:].split()[0] + ' MAEFQ\\nhomolog MTEFQ\\n#=GC RF xxxxx\\n//\\n'
        )
    ''')
    database = write_file(tmp_path / "database.faa", ">homolog\nMTEFQ\n")
    script = PACKAGE / "core" / "msa_generation_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["-d", database, "-j", jackhmmer, "--threads", "1"]
    environment = {"BFF_SEARCH_LOG": str(invocation_log)}
    run_cli(script, ["-i", root, "-o", directory_output, *common], tmp_path, environment)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-i", paths["fasta"], "-o", file_output, *common], tmp_path, environment)
    assert len(invocation_log.read_text().splitlines()) == 4
    for gene in GENES:
        relative = Path(gene) / "MSA" / f"{gene}.msa.a2m"
        assert (file_output / relative).read_text() == (directory_output / relative).read_text()
        assert (file_output / relative).read_text().startswith(f">{gene}\nMAEFQ\n")
        stats = Path(gene) / "MSA" / f"{gene}.msa.stats.json"
        assert json.loads((file_output / stats).read_text()) == json.loads((directory_output / stats).read_text())
        assert json.loads((file_output / stats).read_text())["query_length"] == 5


@pytest.mark.parametrize("compressed_assembly", [True, False], ids=["gzipped-cds", "plain-cds"])
def test_codon_msa_cli_file_directory_backtranslation(gene_tree, tmp_path, compressed_assembly):
    root, files = gene_tree
    tool_bin = tmp_path / "tools"
    mmseqs = executable(tool_bin / "mmseqs", '''
        import json
        import os
        from pathlib import Path
        import sys

        arguments = sys.argv[1:]
        if arguments[0] == 'createdb':
            Path(arguments[2] + '.dbtype').touch()
            if Path(arguments[1]).name == 'query.faa':
                query = Path(arguments[1]).read_text()
                assert query == '>ORF\\nMAEFQ\\n', query
                with open(os.environ['BFF_SEARCH_LOG'], 'a') as handle:
                    handle.write(json.dumps(query) + '\\n')
        elif arguments[0] == 'createindex':
            Path(arguments[1] + '.idx').touch()
        elif arguments[0] == 'convertalis':
            Path(arguments[4]).write_text('ORF\\tNP_000001.1\\t1.0\\t5\\t1.0\\t1.0\\n')
        elif arguments[0] != 'search':
            raise AssertionError(arguments)
    ''')
    executable(tool_bin / "mafft", '''
        from pathlib import Path
        import sys
        content = Path(sys.argv[-1]).read_text()
        assert content == '>seq0\\nMAEFQ\\n>seq1\\nMAEFQ\\n', content
        print(content, end='')
    ''')
    database = tmp_path / "Bio_DBs"
    write_file(database / "refseq_proteins_merged.faa", ">NP_000001.1\nMAEFQ\n")
    assembly = write_file(database / "refseq_assemblies" / "fixture_cds_from_genomic.fna", f">cds [protein_id=NP_000001.1]\n{ORF}\n")
    if compressed_assembly:
        with gzip.open(str(assembly) + ".gz", "wt") as handle:
            handle.write(assembly.read_text())
        assembly.unlink()
    invocation_log = tmp_path / "mmseqs-inputs.jsonl"
    environment = {"PATH": str(tool_bin) + os.pathsep + os.environ["PATH"], "BFF_SEARCH_LOG": str(invocation_log)}
    script = PACKAGE / "core" / "codon_msa_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["-d", database, "-mb", mmseqs, "-t", "1"]
    run_cli(script, ["-i", root, "-o", directory_output, *common], tmp_path, environment)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-f", paths["fasta"], "-o", file_output, *common], tmp_path, environment)
    assert len(invocation_log.read_text().splitlines()) == 4
    for gene in GENES:
        relative = Path(gene) / "CodonMSA" / f"{gene}.codon.msa.fasta"
        assert (file_output / relative).read_text() == (directory_output / relative).read_text()
        assert (file_output / relative).read_text() == ">ORF\nATGGCTGAATTCCAG\n>NP_000001.1\nATGGCTGAATTCCAG\n"
        manifest = Path(gene) / "CodonMSA" / f"{gene}.codon.msa.manifest.tsv"
        assert rows(file_output / manifest) == rows(directory_output / manifest)
        assert [row["status"] for row in rows(file_output / manifest)] == ["PASS", "PASS"]


def write_adabm_parameters(path, alphabet):
    content = "".join(f"h {position} {token} {index / 10}\n" for position in range(5) for index, token in enumerate(alphabet))
    return write_file(path, content)


def test_adabmdca_cli_file_directory_scoring(gene_tree, tmp_path):
    root, files = gene_tree
    parameter_root = tmp_path / "params"
    for gene in GENES:
        write_adabm_parameters(parameter_root / f"{gene}.protein_adabm_params", PROTEIN_ALPHABET)
        write_adabm_parameters(parameter_root / f"{gene}.codon_adabm_params", CODON_ALPHABET)
    script = PACKAGE / "mutation_effects" / "adabmdca_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["--protein-params", parameter_root, "--codon-params", parameter_root, "--skip-train"]
    run_cli(script, ["-f", root, "--msa", root, "-cm", root, "-o", directory_output, *common], tmp_path)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-f", paths["fasta"], "-m", paths["mutations"], "--msa", paths["msa"], "-cm", paths["codon_msa"], "-o", file_output, *common], tmp_path)
    protein = assert_rows_equal(file_output, directory_output, "adabmDCA/{gene}.protein.tsv", ("G4A",))
    codon = assert_rows_equal(file_output, directory_output, "adabmDCA/{gene}.codon.tsv", ("T6C", "C13T"))
    assert float(protein[0]["prediction_protein_independent_adabm"]) == pytest.approx((PROTEIN_ALPHABET.index("T") - PROTEIN_ALPHABET.index("A")) / 10)
    for prediction in codon:
        token = prediction["nt_mutant"]
        if token == "C13T":
            assert prediction["qc_flags"] == "STOP_GAIN"
            assert prediction["prediction_codon_independent_adabm"] == ""
            continue
        wildtype, mutant = ("GCT", "GCC") if token == "T6C" else ("CAG", "TAG")
        expected = (CODON_ALPHABET.index(CODON_TO_CHAR[mutant]) - CODON_ALPHABET.index(CODON_TO_CHAR[wildtype])) / 10
        assert float(prediction["prediction_codon_independent_adabm"]) == pytest.approx(expected, abs=1e-6)


def write_ev_parameters(path, alphabet, focus):
    length, states = len(focus), len(alphabet)
    frequencies = np.tile(np.arange(1, states + 1, dtype="float32"), (length, 1))
    frequencies /= frequencies.sum(axis=1, keepdims=True)
    fields = np.tile(np.arange(states, dtype="float32") / 10, (length, 1))
    pair_count = length * (length - 1) // 2
    content = (
        struct.pack("=5i5f", length, states, 1, 0, 500, 0.2, 0.01, 0.01, 0.0, 1.0)
        + alphabet.encode("ascii") + struct.pack("=f", 1.0) + focus.encode("ascii")
        + np.arange(1, length + 1, dtype="int32").tobytes() + frequencies.tobytes() + fields.tobytes()
        + bytes(pair_count * states * states * 8)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def ev_parameter_root(tmp_path):
    parameters = tmp_path / "params"
    for gene in GENES:
        write_ev_parameters(parameters / f"{gene}.model_params", PROTEIN_ALPHABET, PROTEIN)
        write_ev_parameters(parameters / f"{gene}.codon_model_params", CODON_ALPHABET, "".join(CODON_TO_CHAR[ORF[index:index + 3]] for index in range(0, 15, 3)))
    return parameters


def require_evmutation():
    if importlib.util.find_spec("EVmutation") is None:
        pytest.skip("Real EVmutation dependency required on PYTHONPATH; BFF scoring is not mocked")


def test_evmutation_cli_file_directory_scoring(gene_tree, tmp_path):
    require_evmutation()
    root, files = gene_tree
    parameters = ev_parameter_root(tmp_path)
    script = PACKAGE / "mutation_effects" / "evmutation_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["--model-params", parameters, "--codon-model-params", parameters]
    run_cli(script, ["-f", root, "-o", directory_output, *common], tmp_path)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-f", paths["fasta"], "-m", paths["mutations"], "-o", file_output, *common], tmp_path)
    assert_ev_outputs(file_output, directory_output)


def assert_ev_outputs(file_output, directory_output):
    protein = assert_rows_equal(file_output, directory_output, "EVmutation/{gene}.protein.tsv", ("G4A",))
    codon = assert_rows_equal(file_output, directory_output, "EVmutation/{gene}.codon.tsv", ("T6C", "C13T"))
    assert float(protein[0]["prediction_epistatic"]) == pytest.approx((PROTEIN_ALPHABET.index("T") - PROTEIN_ALPHABET.index("A")) / 10)
    for prediction in codon:
        if prediction["nt_mutant"] == "C13T":
            assert prediction["qc_flags"] == "STOP_GAIN"
            assert prediction["prediction_codon_independent"] == ""
        else:
            assert prediction["prediction_codon_independent"] != ""
            assert float(prediction["prediction_codon_epistatic"]) == pytest.approx(-0.2)


def test_controller_cli_file_directory_nextflow_real_scoring(gene_tree, tmp_path):
    require_evmutation()
    if os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None:
        pytest.skip("Set BFF_TEST_NEXTFLOW=1 with Nextflow installed for actual controller-to-export integration")
    root, files = gene_tree
    runtime = tmp_path / "runtime"
    source = PACKAGE / "mutation_effects"
    runtime.mkdir()
    for source_file in source.glob("*.py"):
        shutil.copy2(source_file, runtime / source_file.name)
    (runtime / "bin").mkdir()
    for source_file in (source / "bin").iterdir():
        if source_file.is_file() and source_file.suffix in {".py", ".nf", ".config"}:
            shutil.copy2(source_file, runtime / "bin" / source_file.name)
    parameters = ev_parameter_root(tmp_path)
    for gene in GENES:
        for suffix in ("model_params", "codon_model_params"):
            destination = root / gene / "models" / f"{gene}.{suffix}"
            destination.parent.mkdir(exist_ok=True)
            shutil.copy2(parameters / f"{gene}.{suffix}", destination)
    hardware = write_file(tmp_path / "hardware.json", json.dumps({"cpus": 2, "memory_gib": 8, "gpus": []}))
    environment = {"NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true", "NXF_ANSI_LOG": "false", "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]}
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["--evmutation-only", "--threads", "1", "--resource-hardware", hardware]
    script = runtime / "mutEffects_controller.py"
    run_cli(script, ["-f", root, "--msa", root, "-cm", root, "--model-params", root, "--codon-model-params", root, "-o", directory_output, *common], tmp_path, environment)
    for paths in detached_files(files, tmp_path):
        gene = paths["fasta"].stem
        run_cli(script, ["-f", paths["fasta"], "-m", paths["mutations"], "--msa", paths["msa"], "-cm", paths["codon_msa"], "--model-params", parameters / f"{gene}.model_params", "--codon-model-params", parameters / f"{gene}.codon_model_params", "-o", file_output, *common], tmp_path, environment)
    assert_ev_outputs(file_output, directory_output)


def test_rare_codon_cli_file_directory_enrichment(gene_tree, tmp_path):
    if importlib.util.find_spec("calc_rare_enrichment") is None:
        pytest.skip("Real cg_cotrans dependency required on PYTHONPATH; BFF enrichment is not mocked")
    from codons import codon_to_aa

    root, files = gene_tree
    reference = "ATG" + "GCT" * 38 + "CAG"
    homolog = "ATG" + "GCC" * 38 + "CAG"
    tokens = ("G97A", "T99C")
    for paths in files.values():
        paths["codon_msa"].write_text(f">ORF\n{reference}\n>homolog\n{homolog}\n")
        paths["mutations"].write_text("mutant\n" + "\n".join(tokens) + "\n")
    usage_by_codon = {codon: 1 / sum(other == amino_acid for other in codon_to_aa.values()) for codon, amino_acid in codon_to_aa.items() if amino_acid != "Stop"}
    usage_by_codon.update({"GCT": 0.05, "GCC": 0.85, "GCA": 0.05, "GCG": 0.05})
    per_sequence = {identifier: usage_by_codon for identifier in ("ORF", "homolog")}
    usage = tmp_path / "codon_usage.p.gz"
    with gzip.open(usage, "wb") as handle:
        pickle.dump({"groups": ["all"], "gene_groups": {gene: 0 for gene in GENES}, "overall_codon_usage": per_sequence, "unweighted_codon_usage": per_sequence, "gene_group_codon_usage": {identifier: {0: usage_by_codon} for identifier in per_sequence}}, handle)
    script = PACKAGE / "rare_codon" / "rare_codon_pipeline.py"
    directory_output, file_output = tmp_path / "directory", tmp_path / "files"
    common = ["-u", usage, "-L", "3"]
    run_cli(script, ["-a", root, "-o", directory_output, *common], tmp_path)
    for paths in detached_files(files, tmp_path):
        run_cli(script, ["-a", paths["codon_msa"], "-m", paths["mutations"], "-o", file_output, *common], tmp_path)
    predictions = assert_rows_equal(file_output, directory_output, "RareCodon/{gene}.rare_codon.tsv", tokens)
    assert all(prediction["codon_position"] == "33" for prediction in predictions)
    assert all(float(prediction["n_rare"]) == 3 for prediction in predictions)
    assert all(float(prediction["f_enriched_wt"]) == 1 for prediction in predictions)
    assert all(0 <= float(prediction["p_enriched"]) <= 1 for prediction in predictions)
