"""Standalone adabmDCA defaults and backend-aware training dispatch."""

import inspect
import sys

import pytest

from biofeaturefactory.mutation_effects import adabmdca_pipeline as pipeline


@pytest.fixture
def parse_cli(tmp_path, monkeypatch):
    fasta_path = tmp_path / "GENE.fasta"
    fasta_path.write_text(">ORF\nATGGCTTAA\n")
    mutations_path = tmp_path / "GENE.csv"
    mutations_path.write_text("Mutation\nG4A\n")
    msa_path = tmp_path / "GENE.msa.fasta"
    msa_path.write_text(">GENE\nMA\n>other\nMT\n")
    captured = []

    def capture_process(gene, fasta_file, mutations_file, output_dir, args):
        captured.append(args)

    monkeypatch.setattr(pipeline, "_process_gene", capture_process)

    def parse_arguments(extra_arguments=()):
        monkeypatch.setattr(sys, "argv", [
            "adabmdca_pipeline.py", "--fasta", str(fasta_path),
            "--mutations", str(mutations_path), "--msa", str(msa_path),
            "--output", str(tmp_path / "output"), "--quiet",
            *extra_arguments,
        ])
        pipeline.main()
        return captured[-1]

    return parse_arguments


def test_standalone_defaults_to_pseudolikelihood_with_backend_epoch_resolution(parse_cli):
    args = parse_cli()

    assert args.adabmdca_model == "pseudoDCA"
    assert args.adabmdca_nepochs is None


@pytest.mark.parametrize("model", [None, "pseudoDCA", "bmDCA", "eaDCA", "edDCA"])
@pytest.mark.parametrize("epochs", [None, 17])
def test_cli_dispatch_preserves_model_and_explicit_epochs(
    parse_cli, tmp_path, monkeypatch, model, epochs,
):
    options = []
    if model is not None:
        options.extend(["--adabmdca-model", model])
    if epochs is not None:
        options.extend(["--adabmdca-nepochs", str(epochs)])
    args = parse_cli(options)
    pseudo_calls = []
    boltzmann_calls = []
    monkeypatch.setattr(
        pipeline, "train_adabmdca_pseudolikelihood_in_process",
        lambda **kwargs: pseudo_calls.append(kwargs),
    )
    monkeypatch.setattr(
        pipeline, "train_adabmdca_in_process",
        lambda **kwargs: boltzmann_calls.append(kwargs),
    )
    params_path = str(tmp_path / "GENE.protein_adabm_params")

    result = pipeline._build_adabmdca_params(
        "GENE", args.msa, "GENE", params_path, args,
    )

    selected_model = model or "pseudoDCA"
    expected_epochs = epochs if epochs is not None else (
        500 if selected_model == "pseudoDCA" else 50000
    )
    assert args.adabmdca_nepochs == epochs
    assert result == params_path
    if selected_model == "pseudoDCA":
        assert not boltzmann_calls
        assert len(pseudo_calls) == 1
        training_call = pseudo_calls[0]
        assert "model" not in training_call
        assert "nchains" not in training_call
        assert "nsweeps" not in training_call
    else:
        assert not pseudo_calls
        assert len(boltzmann_calls) == 1
        training_call = boltzmann_calls[0]
        assert training_call["model"] == selected_model
        assert training_call["nchains"] == 10000
        assert training_call["nsweeps"] == 10
    assert training_call["nepochs"] == expected_epochs
    assert training_call["msa_path"] == args.msa
    assert training_call["output_params_path"] == params_path


def test_default_pseudo_dispatch_does_not_forward_mcmc_options(
    parse_cli, tmp_path, monkeypatch,
):
    args = parse_cli([
        "--adabmdca-nchains", "123", "--adabmdca-nsweeps", "7",
        "--adabmdca-device", "cpu", "--adabmdca-dtype", "float64",
        "--adabmdca-seed", "23", "--adabmdca-lr", "0.02",
        "--adabmdca-tol", "0.004", "--adabmdca-patience", "5",
        "--adabmdca-check-every", "8",
    ])
    training_calls = []
    monkeypatch.setattr(
        pipeline, "train_adabmdca_pseudolikelihood_in_process",
        lambda **kwargs: training_calls.append(kwargs),
    )

    pipeline._build_adabmdca_params(
        "GENE", args.msa, "GENE", str(tmp_path / "params"), args,
    )

    assert len(training_calls) == 1
    training_call = training_calls[0]
    assert "nchains" not in training_call
    assert "nsweeps" not in training_call
    assert "target" not in training_call
    assert training_call["nepochs"] == 500
    assert training_call["device_str"] == "cpu"
    assert training_call["dtype_str"] == "float64"
    assert training_call["seed"] == 23
    assert training_call["lr"] == 0.02
    assert training_call["tol"] == 0.004
    assert training_call["patience"] == 5
    assert training_call["check_every"] == 8


def test_boltzmann_specific_function_defaults_remain_unchanged():
    boltzmann_parameters = inspect.signature(pipeline.train_adabmdca_in_process).parameters
    pseudo_parameters = inspect.signature(
        pipeline.train_adabmdca_pseudolikelihood_in_process,
    ).parameters

    assert boltzmann_parameters["model"].default == "bmDCA"
    assert boltzmann_parameters["nepochs"].default == 50000
    assert pseudo_parameters["nepochs"].default == 500
    assert "nchains" not in pseudo_parameters
    assert "nsweeps" not in pseudo_parameters


def test_cli_help_describes_default_epochs_and_mcmc_scope(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["adabmdca_pipeline.py", "--help"])

    with pytest.raises(SystemExit) as exit_info:
        pipeline.main()

    assert exit_info.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "Training algorithm (default: pseudoDCA)" in help_text
    assert "Default: 500 for pseudoDCA, 50000 for bmDCA/eaDCA/edDCA" in help_text
    assert "MCMC chains for Boltzmann models; unused by pseudoDCA" in help_text
    assert "MCMC sweeps for Boltzmann models; unused by pseudoDCA" in help_text
