"""Exercise CLI input guards without importing licensed or training backends."""

import __future__
import argparse
import ast
import os
import sys
from pathlib import Path

import pytest

from biofeaturefactory.mutation_effects.bin import resource_planner
from biofeaturefactory.lib.utility import read_fasta, trim_muts


REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINES = {
    "codon_usage": ("codon_usage/codon_usage_pipeline.py", "-f", "main"),
    "rare_codon": ("rare_codon/rare_codon_pipeline.py", "-a", "main"),
    "protein_msa": ("core/msa_generation_pipeline.py", "-i", "main"),
    "codon_msa": ("core/codon_msa_pipeline.py", "-i", "main"),
    "controller": ("mutation_effects/mutEffects_controller.py", "-f", "parse_args"),
    "adabmdca": ("mutation_effects/adabmdca_pipeline.py", "-f", "main"),
    "evmutation": ("mutation_effects/evmutation_pipeline.py", "-f", "main"),
}
MUTATION_PIPELINES = tuple(name for name in PIPELINES if name not in {"protein_msa", "codon_msa"})
EFFECT_PIPELINES = ("controller", "adabmdca", "evmutation")


class BoundaryReached(BaseException):
    def __init__(self, arguments, keywords):
        self.arguments = arguments
        self.keywords = keywords


def stop_before_work(*arguments, **keywords):
    raise BoundaryReached(arguments, keywords)


def load_cli(name):
    relative_path, _, entrypoint = PIPELINES[name]
    source = REPO_ROOT / "biofeaturefactory" / relative_path
    tree = ast.parse(source.read_text())
    functions = {entrypoint, "_resolve_genes", "validate_args", "resolve_focus_fastas", "_find_file_for_gene"}
    body = [
        node for node in tree.body
        if (isinstance(node, ast.FunctionDef) and node.name in functions)
        or (isinstance(node, ast.ImportFrom) and node.module == "biofeaturefactory.lib.utility")
    ]
    namespace = {
        "argparse": argparse,
        "os": os,
        "sys": sys,
        "Path": Path,
        "resource_planner": resource_planner,
    }
    executable = compile(ast.Module(body=body, type_ignores=[]), str(source), "exec",
                         flags=__future__.annotations.compiler_flag)
    exec(executable, namespace)
    for boundary in (
        "process_directory", "process_fasta_with_mutations", "_ensure_codon_usage",
        "generate_msa", "generate_codon_msa_from_focus", "_process_gene",
        "_resolve_per_gene_models",
    ):
        namespace[boundary] = stop_before_work
    return namespace, entrypoint


@pytest.mark.parametrize("name", MUTATION_PIPELINES)
@pytest.mark.parametrize("suffix", (".csv", ".tsv", ".txt"))
def test_explicit_mutation_text_suffixes_preserve_reader_behavior(name, suffix, inputs, monkeypatch):
    mutations = inputs["mutations"].with_suffix(suffix)
    mutations.write_text("mutant\nG4A\nT6C\n")
    assert trim_muts(mutations) == ["G4A", "T6C"]
    result = invoke_cli(name, [PIPELINES[name][1], primary_file(name, inputs), "-m", mutations], inputs, monkeypatch)
    if name == "controller":
        assert result.mutations == mutations
    else:
        assert isinstance(result, BoundaryReached)


@pytest.mark.parametrize("name", EFFECT_PIPELINES)
@pytest.mark.parametrize("suffix", (".a2m", ".a3m", ".fna"))
def test_explicit_fasta_format_msa_suffixes(name, suffix, inputs, monkeypatch):
    msa = inputs["msa"].with_suffix(suffix)
    msa.write_text(">PAM\nMA\n>homolog\nMT\n")
    assert read_fasta(msa) == {"PAM": "MA", "homolog": "MT"}
    result = invoke_cli(name, ["-f", inputs["fasta"], "-m", inputs["mutations"], "--msa", msa], inputs, monkeypatch)
    if name == "controller":
        assert result.msa == msa
    else:
        assert isinstance(result, BoundaryReached)


@pytest.mark.parametrize("name", EFFECT_PIPELINES)
def test_stockholm_is_not_a_supported_fasta_msa(name, inputs, monkeypatch, capsys):
    msa = inputs["msa"].with_suffix(".sto")
    msa.write_text("# STOCKHOLM 1.0\nPAM MA\nhomolog MT\n//\n")
    with pytest.raises(KeyError):
        read_fasta(msa)
    with pytest.raises(SystemExit) as error:
        invoke_cli(name, ["-f", inputs["fasta"], "-m", inputs["mutations"], "--msa", msa], inputs, monkeypatch)
    assert error.value.code == 2
    assert "expected .a2m" in capsys.readouterr().err


@pytest.fixture
def inputs(tmp_path):
    root = tmp_path / "inputs"
    gene = root / "PAM"
    files = {
        "fasta": gene / "fastas" / "PAM.fasta",
        "mutations": gene / "mappings" / "mutations" / "PAM_mutations.csv",
        "msa": gene / "MSA" / "PAM.msa.a2m",
        "codon_msa": gene / "CodonMSA" / "PAM.codon.msa.fasta",
        "params": gene / "MSA" / "PAM.model_params",
    }
    for path in files.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    files["fasta"].write_text(">ORF\nATGGCCTAA\n")
    files["mutations"].write_text("mutant\nG4A\n")
    files["msa"].write_text(">PAM\nMA\n")
    files["codon_msa"].write_text(">ORF\nATGGCC\n")
    files["params"].write_bytes(b"parameters are not loaded by these tests")
    reference = tmp_path / "reference.tsv"
    reference.write_text("not read by the input-boundary tests\n")
    return {**files, "root": root, "gene": gene, "reference": reference,
            "output": tmp_path / "output", "db_root": tmp_path / "database"}


def invoke_cli(name, arguments, inputs, monkeypatch):
    namespace, entrypoint = load_cli(name)
    argv = [name, *map(str, arguments), "-o", str(inputs["output"])]
    if name == "protein_msa":
        argv.extend(["-d", str(inputs["reference"]), "-j", "not-executed"])
    elif name == "codon_msa":
        argv.extend(["-d", str(inputs["db_root"])])
    monkeypatch.setattr(sys, "argv", argv)
    try:
        return namespace[entrypoint]()
    except BoundaryReached as boundary:
        return boundary


def primary_file(name, inputs):
    return inputs["codon_msa"] if name == "rare_codon" else inputs["fasta"]


@pytest.mark.parametrize("name", MUTATION_PIPELINES)
def test_file_mode_requires_explicit_mutations(name, inputs, monkeypatch, capsys):
    with pytest.raises(SystemExit) as error:
        invoke_cli(name, [PIPELINES[name][1], primary_file(name, inputs)], inputs, monkeypatch)
    assert error.value.code == 2
    assert "--mutations is required in file mode" in capsys.readouterr().err
    assert not inputs["output"].exists()


@pytest.mark.parametrize("name", tuple(PIPELINES))
def test_explicit_files_reach_existing_work_boundary(name, inputs, monkeypatch):
    arguments = [PIPELINES[name][1], primary_file(name, inputs)]
    if name in MUTATION_PIPELINES:
        arguments.extend(["-m", inputs["mutations"]])
    result = invoke_cli(name, arguments, inputs, monkeypatch)
    if name == "controller":
        assert result.mutations == inputs["mutations"]
    else:
        assert isinstance(result, BoundaryReached)


@pytest.mark.parametrize("name", MUTATION_PIPELINES)
def test_parent_directory_preserves_automatic_mutations(name, inputs, monkeypatch):
    result = invoke_cli(name, [PIPELINES[name][1], inputs["root"]], inputs, monkeypatch)
    if name == "controller":
        selected = result.mutations
    elif name == "rare_codon":
        selected = result.arguments[0].mutations
    elif name == "adabmdca":
        selected = result.arguments[2]
    elif name == "evmutation":
        selected = result.arguments[1].mutations
    else:
        selected = result.arguments[1]
    expected = inputs["mutations"] if name == "adabmdca" else inputs["root"]
    assert Path(selected) == expected


@pytest.mark.parametrize("name", tuple(PIPELINES))
@pytest.mark.parametrize("subdirectory", ("", "fastas", "MSA", "CodonMSA"))
def test_nested_input_directories_are_rejected(name, subdirectory, inputs, monkeypatch, capsys):
    nested = inputs["gene"] / subdirectory
    with pytest.raises(SystemExit) as error:
        invoke_cli(name, [PIPELINES[name][1], nested], inputs, monkeypatch)
    assert error.value.code == 2
    message = capsys.readouterr().err
    assert "parent <dir>" in message
    assert str(inputs["root"]) in message
    assert not inputs["output"].exists()


@pytest.mark.parametrize("name", MUTATION_PIPELINES)
@pytest.mark.parametrize("mutations_first", (False, True))
@pytest.mark.parametrize("first_mode", ("file", "directory"))
def test_first_supplied_input_controls_mixed_mode_error(
    name, mutations_first, first_mode, inputs, monkeypatch, capsys,
):
    primary = (PIPELINES[name][1], primary_file(name, inputs))
    companion = ("-m", inputs["mutations"])
    first, second = (companion, primary) if mutations_first else (primary, companion)
    if first_mode == "directory":
        first = (first[0], inputs["root"])
    else:
        second = (second[0], inputs["root"])
    with pytest.raises(SystemExit) as error:
        invoke_cli(name, [*first, *second], inputs, monkeypatch)
    assert error.value.code == 2
    assert f"{first[0]} selected {first_mode} mode" in capsys.readouterr().err
    assert not inputs["output"].exists()


@pytest.mark.parametrize("name", EFFECT_PIPELINES)
@pytest.mark.parametrize("msa_flag, msa_kind", (("--msa", "msa"), ("--codon-msa", "codon_msa")))
def test_msa_option_can_select_file_mode_first(name, msa_flag, msa_kind, inputs, monkeypatch, capsys):
    with pytest.raises(SystemExit) as error:
        invoke_cli(name, [msa_flag, inputs[msa_kind], "-f", inputs["root"], "-m", inputs["root"]],
                   inputs, monkeypatch)
    assert error.value.code == 2
    assert f"{msa_flag} selected file mode" in capsys.readouterr().err


@pytest.mark.parametrize("name", EFFECT_PIPELINES)
def test_explicit_alignment_files_preserve_backend_flags(name, inputs, monkeypatch):
    result = invoke_cli(name, ["-f", inputs["fasta"], "-m", inputs["mutations"],
                              "--msa", inputs["msa"], "-cm", inputs["codon_msa"]],
                        inputs, monkeypatch)
    if name == "controller":
        assert result.msa == inputs["msa"]
        assert result.codon_msa == inputs["codon_msa"]
        assert result.run_evmutation is True
        assert result.run_adabmdca is True
    else:
        assert isinstance(result, BoundaryReached)


@pytest.mark.parametrize("params_flag", (
    "--model-params", "--codon-model-params", "--adabmdca-protein-params", "--adabmdca-codon-params",
))
def test_controller_prebuilt_params_are_inputs(params_flag, inputs, monkeypatch, capsys):
    with pytest.raises(SystemExit) as error:
        invoke_cli("controller", [params_flag, inputs["params"], "-f", inputs["root"]], inputs, monkeypatch)
    assert error.value.code == 2
    assert f"{params_flag} selected file mode" in capsys.readouterr().err


@pytest.mark.parametrize("name, params_flag", (
    ("adabmdca", "--protein-params"), ("evmutation", "--model-params"),
))
def test_standalone_parameter_destinations_remain_outputs(name, params_flag, inputs, monkeypatch):
    destination = inputs["output"].parent / "new_params"
    result = invoke_cli(name, ["-f", inputs["fasta"], "-m", inputs["mutations"],
                              "--msa", inputs["msa"], params_flag, destination], inputs, monkeypatch)
    assert isinstance(result, BoundaryReached)
    assert not destination.exists()


@pytest.mark.parametrize("name", ("protein_msa", "codon_msa"))
def test_core_parent_directory_reaches_generation(name, inputs, monkeypatch):
    result = invoke_cli(name, ["-i", inputs["root"]], inputs, monkeypatch)
    assert isinstance(result, BoundaryReached)


def test_protein_msa_alias_cannot_hide_mixed_inputs(inputs, monkeypatch, capsys):
    with pytest.raises(SystemExit) as error:
        invoke_cli("protein_msa", ["-f", inputs["fasta"], "-i", inputs["root"]], inputs, monkeypatch)
    assert error.value.code == 2
    assert "-f selected file mode" in capsys.readouterr().err


def test_codon_msa_file_only_alias_stays_available(inputs, monkeypatch):
    result = invoke_cli("codon_msa", ["-f", inputs["fasta"]], inputs, monkeypatch)
    assert isinstance(result, BoundaryReached)


def test_codon_msa_file_only_alias_rejects_directory(inputs, monkeypatch, capsys):
    with pytest.raises(SystemExit) as error:
        invoke_cli("codon_msa", ["-f", inputs["root"]], inputs, monkeypatch)
    assert error.value.code == 2
    assert "--fasta is not a file" in capsys.readouterr().err


def test_reference_and_usage_files_do_not_select_directory_mode(inputs, monkeypatch):
    result = invoke_cli("rare_codon", ["-u", inputs["reference"], "-rcu", inputs["reference"],
                                       "-a", inputs["root"]], inputs, monkeypatch)
    assert isinstance(result, BoundaryReached)
    assert Path(result.arguments[0].mutations) == inputs["root"]
