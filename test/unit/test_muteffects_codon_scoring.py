"""Explicit codon missense scoring and preservation of unscorable rows."""

import csv
import importlib
import importlib.util
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from biofeaturefactory.mutation_effects import adabmdca_pipeline as adabm
from biofeaturefactory.mutation_effects.bin.codon_encoding import CODON_ALPHABET, CODON_TO_CHAR


@pytest.fixture(params=["adabm", "evmutation"])
def backend(request, monkeypatch):
    if request.param == "adabm":
        return request.param, adabm
    package = ModuleType("EVmutation")
    model = ModuleType("EVmutation.model")
    model.CouplingsModel = object
    vendor_tools = ModuleType("EVmutation.tools")
    package.model = model
    package.tools = vendor_tools
    monkeypatch.setitem(sys.modules, "EVmutation", package)
    monkeypatch.setitem(sys.modules, "EVmutation.model", model)
    monkeypatch.setitem(sys.modules, "EVmutation.tools", vendor_tools)
    filename = Path(adabm.__file__).with_name("evmutation_pipeline.py")
    specification = importlib.util.spec_from_file_location("bff_codon_scoring_test_evmutation", filename)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return request.param, module


def _score(backend, mutations, orf="ATGGCTGAA", enabled=True, **kwargs):
    name, module = backend
    scorer = module.score_nt_mutations_adabm if name == "adabm" else module.score_nt_mutations
    return scorer(mutations, "NPM1", orf, aa_lookup={}, score_missense_codon=enabled, **kwargs)


def _lookup(name, position=2, codon="ACT"):
    if name == "adabm":
        score = {"indep": 1.5, "epi": 2.5, "pair": 1.0, "concordance": "CONCORDANT", "frequency": 0.2}
    else:
        score = {
            "prediction_codon_independent": 1.5, "prediction_codon_epistatic": 2.5,
            "codon_epistatic_contribution": 1.0, "codon_epistatic_concordance": "CONCORDANT",
            "codon_frequency": 0.2,
        }
    return {(position, codon): score}


def _epistatic_column(name):
    return "prediction_codon_epistatic_adabm" if name == "adabm" else "prediction_codon_epistatic"


def test_explicit_codon_mode_scores_missense_without_protein_model(backend):
    name, module = backend
    protein, codon = _score(backend, ["G4A"], codon_lookup=_lookup(name))
    assert protein == []
    assert len(codon) == 1
    assert codon[0][_epistatic_column(name)] == 2.5
    assert codon[0]["mutation_class"] == "MISSENSE"
    assert "MISSENSE_CODON_LEVEL" in codon[0]["qc_flags"]
    assert "MISSENSE_SCORED" in codon[0]["qc_flags"]
    assert "NO_PROTEIN_MODEL" not in codon[0]["qc_flags"]


def test_default_missense_routing_remains_protein(backend):
    name, module = backend
    protein, codon = _score(backend, ["G4A"], enabled=False, codon_lookup=_lookup(name))
    assert len(protein) == 1
    assert codon == []
    assert "MISSENSE_CODON_LEVEL" not in protein[0]["qc_flags"]


@pytest.mark.parametrize("lookup,reason", [(None, "MISSENSE_UNSCORED"), ({}, "MISSENSE_NOT_IN_CODON_MODEL")])
def test_missing_codon_scores_remain_explicitly_unscored(backend, lookup, reason):
    name, module = backend
    protein, codon = _score(backend, ["G4A"], codon_lookup=lookup)
    assert protein == []
    assert codon[0][_epistatic_column(name)] == ""
    assert reason in codon[0]["qc_flags"]


@pytest.mark.parametrize("mutation,orf,reason", [
    ("not-a-variant", "ATGGCTGAA", "INVALID_MUTATION"),
    ("A99C", "ATGGCTGAA", "OUT_OF_RANGE"),
    ("A4C", "ATGA", "PARTIAL_CODON"),
    ("A4C", "ATGANN", "UNKNOWN_CODON"),
    ("GCT8ACA", "ATGGCTGAA", "OUT_OF_RANGE"),
    ("GC4A", "ATGGCTGAA", "FRAMESHIFT_NOT_REPRESENTABLE_FIXED_L"),
    ("GCTG4G", "ATGGCTGAA", "CODON_LENGTH_CHANGE_UNSCORED"),
])
def test_codon_only_mode_preserves_unscorable_variants(backend, mutation, orf, reason):
    name, module = backend
    protein, codon = _score(backend, [mutation], orf=orf)
    assert protein == []
    assert len(codon) == 1
    assert codon[0]["nt_mutant"] == mutation
    assert codon[0][_epistatic_column(name)] == ""
    assert reason in codon[0]["qc_flags"]


def test_stop_rows_remain_annotation_only_even_when_lookup_has_scores(backend):
    name, module = backend
    protein, codon = _score(backend, ["G7T"], codon_lookup=_lookup(name, 3, "TAA"))
    assert protein == []
    assert codon[0]["mutation_class"] == "STOP_GAIN"
    assert codon[0][_epistatic_column(name)] == ""
    assert "STOP_GAIN" in codon[0]["qc_flags"]


def _codon_models(mutant_codons):
    tokens = {symbol: index for index, symbol in enumerate(CODON_ALPHABET)}
    wild_types = [CODON_TO_CHAR[codon] for codon in ("ATG", "GCT", "GAA")]
    mutant_tokens = [CODON_TO_CHAR[codon] for codon in mutant_codons]
    fields = np.zeros((3, 65))
    fields[1, tokens[mutant_tokens[0]]] = 2
    fields[2, tokens[mutant_tokens[1]]] = 3
    couplings = np.zeros((3, 65, 3, 65))
    couplings[1, tokens[mutant_tokens[0]], 2, tokens[mutant_tokens[1]]] = 5
    couplings[2, tokens[mutant_tokens[1]], 1, tokens[mutant_tokens[0]]] = 5
    frequencies = np.full((3, 65), 1 / 65)
    target_indices = np.array([tokens[symbol] for symbol in wild_types])
    context = adabm.ModelContext(fields, couplings, {0: 1, 1: 2, 2: 3}, target_indices, frequencies, CODON_ALPHABET)
    common = {
        "L": 3, "index_list": [101, 102, 103], "index_map": {101: 0, 102: 1, 103: 2},
        "alphabet_map": tokens, "target_seq": wild_types, "target_seq_mapped": target_indices,
        "h_i": fields, "fi": lambda position, symbol: 1 / 65,
    }
    epistatic = SimpleNamespace(**common, J_ij=couplings.transpose(0, 2, 1, 3))
    independent = SimpleNamespace(**common, J_ij=np.zeros((3, 3, 65, 65)))
    return context, (epistatic, independent)


@pytest.mark.parametrize("mutants,enabled,mutation_class", [
    (("ACA", "TTC"), True, "MISSENSE"),
    (("GCC", "GAG"), False, "SYNONYMOUS"),
])
def test_multicodon_scores_include_joint_coupling_and_correct_codon_encoding(backend, mutants, enabled, mutation_class):
    name, module = backend
    context, models = _codon_models(mutants)
    model_args = {"codon_ctx": context} if name == "adabm" else {"codon_models": models}
    protein, codon = _score(backend, ["GCTGAA4" + "".join(mutants)], enabled=enabled, **model_args)
    assert protein == []
    assert codon[0]["mutation_class"] == mutation_class
    assert codon[0][_epistatic_column(name)] == pytest.approx(10)
    independent_column = "prediction_codon_independent_adabm" if name == "adabm" else "prediction_codon_independent"
    assert codon[0][independent_column] == pytest.approx(5)
    assert "CONCORDANCE_UNDEFINED_MULTICODON" in codon[0]["qc_flags"]
    assert f"{mutation_class}_SCORED" in codon[0]["qc_flags"]


def test_codon_mode_never_builds_protein_params(backend, monkeypatch):
    name, module = backend
    arguments = SimpleNamespace(
        msa="unused-protein.msa", protein_params="unused-protein.params", model_params="unused-protein.params",
        codon_msa=None, codon_params=None, codon_model_params=None,
        skip_codon=False, score_missense_codon=True,
    )
    builder = "_build_adabmdca_params" if name == "adabm" else "_build_model_params"
    monkeypatch.setattr(module, builder, lambda *args, **kwargs: pytest.fail("Protein training was attempted"))
    resolver = module._resolve_per_gene_adabm_params if name == "adabm" else module._resolve_per_gene_models
    assert all(value is None for value in resolver("NPM1", arguments))


def test_intronic_rows_are_kept_in_the_selected_codon_output(backend, monkeypatch, tmp_path):
    name, module = backend
    fasta = tmp_path / "NPM1.fasta"
    fasta.write_text(">ORF\nATGGCTGAA\n")
    mutations = tmp_path / "NPM1_mutations.csv"
    mutations.write_text("mutant\ngd.T5000C\n")
    arguments = SimpleNamespace(validation_log=None, skip_codon=False, score_missense_codon=True, quiet=True)
    if name == "adabm":
        counts = module._process_gene("NPM1", str(fasta), str(mutations), tmp_path, arguments)
        subdirectory = "adabmDCA"
    else:
        counts = module._process_gene("NPM1", str(fasta), str(mutations), None, None, tmp_path, arguments)
        subdirectory = "EVmutation"
    assert counts[:2] == (0, 1)
    with (tmp_path / "NPM1" / subdirectory / "NPM1.codon.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["nt_mutant"] == "gd.T5000C"
    assert "NON_ORF_TOKEN" in rows[0]["qc_flags"]


def test_real_gene_processing_scores_codon_missense_and_skips_protein_loading(backend, monkeypatch, tmp_path):
    name, module = backend
    fasta = tmp_path / "NPM1.fasta"
    fasta.write_text(">ORF\nATGGCTGAA\n")
    mutations = tmp_path / "NPM1_mutations.csv"
    mutations.write_text("mutant\nG4A\n")
    protein_params = tmp_path / "unused-protein.params"
    protein_params.write_text("must not be loaded")
    codon_params = tmp_path / "codon.params"
    codon_params.write_text("fixture model")
    arguments = SimpleNamespace(validation_log=None, skip_codon=False, score_missense_codon=True, quiet=True, codon_focus="ORF")
    if name == "adabm":
        monkeypatch.setattr(module, "_resolve_per_gene_adabm_params", lambda *args: (str(protein_params), str(codon_params), "unused-protein.msa", "encoded-codon.msa"))
        monkeypatch.setattr(module, "_build_protein_lookup_from_params", lambda *args: pytest.fail("Protein model was loaded"))
        monkeypatch.setattr(module, "_build_codon_lookup_from_params", lambda *args: (_lookup(name), None))
        counts = module._process_gene("NPM1", str(fasta), str(mutations), tmp_path, arguments)
        subdirectory = "adabmDCA"
    else:
        loaded = []

        def load(filename):
            loaded.append(filename)
            return SimpleNamespace()

        monkeypatch.setattr(module, "CouplingsModel", load)
        monkeypatch.setattr(module, "build_aa_prediction_lookup", lambda *args: pytest.fail("Protein model was loaded"))
        monkeypatch.setattr(module, "build_codon_lookup", lambda *args: (_lookup(name), None))
        counts = module._process_gene("NPM1", str(fasta), str(mutations), str(protein_params), str(codon_params), tmp_path, arguments)
        assert loaded == [str(codon_params)]
        subdirectory = "EVmutation"
    assert counts[:2] == (0, 1)
    with (tmp_path / "NPM1" / subdirectory / "NPM1.codon.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert float(rows[0][_epistatic_column(name)]) == 2.5
    assert "MISSENSE_CODON_LEVEL" in rows[0]["qc_flags"]


def test_cli_threads_codon_flag_into_real_processing_entrypoint(backend, monkeypatch, tmp_path):
    name, module = backend
    fasta = tmp_path / "NPM1.fasta"
    fasta.write_text(">ORF\nATGGCTGAA\n")
    mutations = tmp_path / "NPM1_mutations.csv"
    mutations.write_text("mutant\nG4A\n")
    captured = []

    def process(*arguments):
        captured.append(arguments[-1])
        return (0, 1) if name == "adabm" else (0, 1, 0, 0)

    monkeypatch.setattr(module, "_process_gene", process)
    monkeypatch.setattr(sys, "argv", [
        "pipeline", "--fasta", str(fasta), "--mutations", str(mutations),
        "--output", str(tmp_path / "output"), "--score-missense-codon",
    ])
    module.main()
    assert len(captured) == 1
    assert captured[0].score_missense_codon is True


def test_codon_missense_flag_rejects_skip_codon(backend, monkeypatch):
    name, module = backend
    with pytest.raises(ValueError, match="skip_codon"):
        _score(backend, ["G4A"], skip_codon=True)
    monkeypatch.setattr(sys, "argv", ["pipeline", "--fasta", "unused", "--score-missense-codon", "--skip-codon"])
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2


def _write_native_codon_params(path, backend_name):
    fields = np.zeros((3, 65), dtype=np.float32)
    fields[1, CODON_ALPHABET.index(CODON_TO_CHAR["ACT"])] = 2.5
    if backend_name == "adabm":
        path.write_text("".join(
            f"h {position} {symbol} {fields[position, symbol_index]}\n"
            for position in range(3) for symbol_index, symbol in enumerate(CODON_ALPHABET)
        ))
        return
    with path.open("wb") as handle:
        np.array([3, 65, 1, 0, 1], dtype=np.int32).tofile(handle)
        np.array([0.8, 0.01, 16.2, 0, 1], dtype=np.float32).tofile(handle)
        np.array(list(CODON_ALPHABET), dtype="S1").tofile(handle)
        np.array([1], dtype=np.float32).tofile(handle)
        np.array([CODON_TO_CHAR[codon] for codon in ("ATG", "GCT", "GAA")], dtype="S1").tofile(handle)
        np.array([1, 2, 3], dtype=np.int32).tofile(handle)
        np.full((3, 65), 1 / 65, dtype=np.float32).tofile(handle)
        fields.tofile(handle)
        for pair in range(3):
            np.full((65, 65), 1 / 65 ** 2, dtype=np.float32).tofile(handle)
        for pair in range(3):
            np.zeros((65, 65), dtype=np.float32).tofile(handle)


@pytest.mark.parametrize("name", ["adabm", "evmutation"])
def test_native_params_cli_writes_actual_codon_scores_and_qc_rows(name, monkeypatch, tmp_path, capsys):
    if name == "adabm":
        module = adabm
    else:
        try:
            module = importlib.import_module("biofeaturefactory.mutation_effects.evmutation_pipeline")
        except ModuleNotFoundError as error:
            if error.name and error.name.startswith("EVmutation"):
                pytest.skip("Native EVmutation model dependency is not installed")
            raise
    fasta = tmp_path / "NPM1.fasta"
    fasta.write_text(">ORF\nATGGCTGAA\n")
    mutations = tmp_path / "NPM1_mutations.csv"
    requested = ["G4A", "A99C", "gd.T5000C", "G7T", "GC4A"]
    mutations.write_text("mutant\n" + "\n".join(requested) + "\n")
    msa = tmp_path / "NPM1.codon.msa.fasta"
    msa.write_text(">ORF\nATGGCTGAA\n>homolog\nATGACTGAA\n")
    params = tmp_path / "native.params"
    _write_native_codon_params(params, name)
    params_option = "--codon-params" if name == "adabm" else "--codon-model-params"
    skip_training = "--skip-train" if name == "adabm" else "--skip-plmc"
    output = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", [
        "pipeline", "--fasta", str(fasta), "--mutations", str(mutations),
        "--codon-msa", str(msa), params_option, str(params), skip_training,
        "--score-missense-codon", "--quiet", "--output", str(output),
    ])
    module.main()
    subdirectory = "adabmDCA" if name == "adabm" else "EVmutation"
    gene_directory = output / "NPM1" / subdirectory
    with (gene_directory / "NPM1.codon.tsv").open() as handle:
        rows = {row["nt_mutant"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with (gene_directory / "NPM1.protein.tsv").open() as handle:
        protein = list(csv.DictReader(handle, delimiter="\t"))
    assert protein == []
    assert set(rows) == set(requested)
    score_column = _epistatic_column(name)
    assert float(rows["G4A"][score_column]) == pytest.approx(2.5)
    assert "MISSENSE_CODON_LEVEL" in rows["G4A"]["qc_flags"]
    for mutation, reason in (
        ("A99C", "OUT_OF_RANGE"), ("gd.T5000C", "NON_ORF_TOKEN"),
        ("G7T", "STOP_GAIN"), ("GC4A", "FRAMESHIFT_NOT_REPRESENTABLE_FIXED_L"),
    ):
        assert rows[mutation][score_column] == ""
        assert reason in rows[mutation]["qc_flags"]
    captured = capsys.readouterr()
    assert "gd.T5000C" in captured.err
