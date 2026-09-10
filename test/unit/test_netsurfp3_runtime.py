"""Residue coverage, overlap ownership and incremental NetSurfP3 publication."""

import csv
from pathlib import Path

import numpy as np
import pytest
import torch
from nsp3.base.base_predict import BasePredict

from biofeaturefactory.NetSurfP3 import netsurfp3_pipeline as pipeline
from nsp3.embeddings.esm1b import ESM1bEmbedding


def test_loaded_embedding_preserves_last_residue_with_residue_length_mask(monkeypatch):
    class TokenModel(torch.nn.Module):
        def forward(self, tokens, repr_layers=None):
            return {"representations": {33: tokens.unsqueeze(-1).float()}}

    embedding = object.__new__(ESM1bEmbedding)
    torch.nn.Module.__init__(embedding)
    embedding.max_embedding = 1024
    embedding.offset = 200
    embedding.model = TokenModel()
    model = torch.nn.Module()
    model.embedding = embedding
    monkeypatch.setattr(torch, "load", lambda *args, **kwargs: {"state_dict": {}})
    predictor = object.__new__(BasePredict)
    pipeline._patched_base_init(predictor, model, "unused")
    tokens = torch.tensor([[0, 1, 2, 3, 4, 0]])
    assert predictor.model.embedding(tokens, padding_length=4).flatten().tolist() == [1, 2, 3, 4]
    assert ESM1bEmbedding.forward is not pipeline._bounded_esm_forward


@pytest.mark.parametrize("length", [1, 974, 1021])
def test_bounded_embedding_matches_complete_token_input(length):
    embedding = object.__new__(ESM1bEmbedding)
    torch.nn.Module.__init__(embedding)
    embedding.max_embedding = 1024
    embedding.model = lambda tokens, repr_layers=None: {
        "representations": {33: tokens.unsqueeze(-1).float()}}
    tokens = torch.arange(length + 2).unsqueeze(0)
    result = pipeline._bounded_esm_forward(embedding, tokens, padding_length=length)
    assert result.flatten().tolist() == list(range(1, length + 1))


@pytest.mark.parametrize("length", [974, 1021, 1022, 1480, 1846, 1847, 2045, 2046,
                                         2670, 2671, 3069, 3070, 3418, 3494, 3495,
                                         3685, 4093, 4094])
def test_actual_embedding_forward_has_complete_ordered_wrapper_coverage(tmp_path, monkeypatch, length):
    embedding = object.__new__(ESM1bEmbedding)
    torch.nn.Module.__init__(embedding)
    embedding.max_embedding = 1024
    embedding.offset = 200
    calls = []

    def model(tokens, repr_layers=None):
        calls.append(tokens.shape[1])
        return {"representations": {33: tokens.unsqueeze(-1).float()}}

    embedding.model = model
    sequence = ("ACDEFGHIKLMNPQRSTVWY" * ((length // 19) + 1))[:length]
    fasta = tmp_path / "long.fasta"
    fasta.write_text(f">protein\n{sequence}\n")
    expected_chunks = []

    def predict(config, feature, checkpoint, fasta_path):
        entries = pipeline.read_fasta(fasta_path)
        identifiers, sequences, predictions = [], [], []
        for identifier, residues in entries.items():
            expected_chunks.append(identifier)
            tokens = torch.tensor([[0, *map(ord, residues), 0]])
            values = embedding.forward(tokens).numpy()
            identifiers.append([identifier])
            sequences.append([residues])
            predictions.append([values])
        return identifiers, sequences, predictions

    monkeypatch.setattr(pipeline, "load_config", lambda path: {})
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor",
                        lambda model_path, config_path, device:
                        lambda fasta_path: predict({}, "SecondaryFeatures", model_path, fasta_path))
    monkeypatch.setattr(pipeline, "extract_residue_predictions",
                        lambda tensors, seq_idx, pos_idx, residue:
                        {"residue": residue, "value": tensors[0][seq_idx, pos_idx, 0]})
    results = pipeline.run_nsp3_prediction(fasta, "unused", "unused", max_seq_length=10000)
    assert list(results["protein"]) == list(range(1, length + 1))
    assert [row["value"] for row in results["protein"].values()] == list(map(ord, sequence))
    assert len(calls) == len(expected_chunks)
    assert max(calls) <= 1023


def test_overlap_uses_interior_not_first_chunk(tmp_path, monkeypatch):
    fasta = tmp_path / "seam.fasta"
    fasta.write_text(">protein\n" + "M" * 1500 + "\n")

    def predict(config, feature, checkpoint, fasta_path):
        identifiers, sequences, predictions = [], [], []
        for identifier, residues in pipeline.read_fasta(fasta_path).items():
            chunk_number = int(identifier.rsplit("chunk", 1)[1])
            identifiers.append([identifier])
            sequences.append([residues])
            predictions.append([np.full((1, len(residues), 1), chunk_number)])
        return identifiers, sequences, predictions

    monkeypatch.setattr(pipeline, "load_config", lambda path: {})
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor",
                        lambda model_path, config_path, device:
                        lambda fasta_path: predict({}, "SecondaryFeatures", model_path, fasta_path))
    monkeypatch.setattr(pipeline, "extract_residue_predictions",
                        lambda tensors, seq_idx, pos_idx, residue:
                        {"owner": tensors[0][seq_idx, pos_idx, 0]})
    results = pipeline.run_nsp3_prediction(fasta, "unused", "unused")["protein"]
    assert results[980]["owner"] == 0
    assert results[1000]["owner"] == 1
    assert len(results) == 1500


def test_missing_chunk_is_not_published_as_complete(tmp_path, monkeypatch):
    fasta = tmp_path / "missing.fasta"
    fasta.write_text(">protein\nMMM\n")
    monkeypatch.setattr(pipeline, "load_config", lambda path: {})
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor",
                        lambda *args: lambda fasta_path: ([], [], []))
    with pytest.raises(RuntimeError, match="0/3 residues"):
        pipeline.run_nsp3_prediction(fasta, "unused", "unused")


def test_completed_genes_survive_later_prediction_failure(tmp_path, monkeypatch, capsys):
    input_root = tmp_path / "input"
    output_root = tmp_path / "output"
    for gene in ["FIRST", "BROKEN", "LAST"]:
        fasta_dir = input_root / gene / "fastas"
        mutation_dir = input_root / gene / "mappings" / "mutations"
        fasta_dir.mkdir(parents=True)
        mutation_dir.mkdir(parents=True)
        (fasta_dir / f"{gene}.fasta").write_text(">ORF\nATGAAAACCTAA\n")
        (mutation_dir / f"{gene}_mutations.csv").write_text("mutant\nA4G\n")
    visited = []
    runtimes = []

    def predict(fasta, *args, **kwargs):
        runtimes.append(kwargs.get("runtime"))
        sequences = pipeline.read_fasta(fasta)
        gene = next(iter(sequences)).split("-")[0]
        visited.append(gene)
        if gene == "BROKEN":
            assert (output_root / "FIRST/NetSurfP3/FIRST.netsurfp3.summary.tsv").is_file()
            raise RuntimeError("simulated model OOM")
        predictions = {}
        for identifier, sequence in sequences.items():
            tensors = [np.ones((1, len(sequence), width)) / width for width in (8, 3, 2, 1, 1, 1)]
            predictions[identifier] = {
                position + 1: pipeline.extract_residue_predictions(tensors, 0, position, residue)
                for position, residue in enumerate(sequence)
            }
        return predictions

    monkeypatch.setattr(pipeline, "discover_fasta_files", lambda path: {
        gene: str(input_root / gene / "fastas" / f"{gene}.fasta")
        for gene in ["FIRST", "BROKEN", "LAST"]})
    monkeypatch.setattr(pipeline, "run_nsp3_prediction", predict)
    monkeypatch.setattr("sys.argv", ["netsurfp3", "-i", str(input_root), "-o", str(output_root),
                                     "-M", "unused", "-c", "unused"])
    assert pipeline.main() == 1
    assert visited == ["FIRST", "BROKEN", "LAST"]
    assert all(runtime is not None for runtime in runtimes)
    assert len({id(runtime) for runtime in runtimes}) == 1
    for gene in visited:
        directory = output_root / gene / "NetSurfP3"
        assert len(list(directory.glob("*.tsv"))) == 3
        with (directory / f"{gene}.netsurfp3.summary.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        assert ("PREDICTION_FAILED" in rows[0]["qc_flags"]) == (gene == "BROKEN")
    assert "simulated model OOM" in capsys.readouterr().err


def test_publish_failure_does_not_replace_previous_tables(tmp_path, monkeypatch):
    directory = tmp_path / "GENE" / "NetSurfP3"
    directory.mkdir(parents=True)
    paths = [directory / f"GENE.netsurfp3.{suffix}.tsv" for suffix in ("summary", "residues", "local")]
    for path in paths:
        path.write_text("previous output\n")

    def fail_write(*args):
        raise OSError("disk failure")

    monkeypatch.setattr(pipeline, "write_local_tsv", fail_write)
    with pytest.raises(OSError, match="disk failure"):
        pipeline._publish_gene_outputs(tmp_path, "GENE", [], [], [])
    assert all(path.read_text() == "previous output\n" for path in paths)
