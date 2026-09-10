"""Predictor lifetime, inference gradients, adaptive memory limits and OOM recovery."""

from contextlib import nullcontext
from types import SimpleNamespace
import weakref

import numpy as np
import pytest
import torch

from test_netsurfp3_checkpoint import runtime


@pytest.fixture
def cuda_memory(monkeypatch):
    unit = 1024 * 1024
    state = {
        "allocated": 1000 * unit,
        "reserved": 1000 * unit,
        "free": 1000 * unit,
        "total": 2000 * unit,
        "peak": 1100 * unit,
        "releases": 0,
    }
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "synchronize", lambda *args, **kwargs: None)
    monkeypatch.setattr(torch.cuda, "reset_peak_memory_stats", lambda *args, **kwargs: None)
    monkeypatch.setattr(torch.cuda, "memory_allocated", lambda *args, **kwargs: state["allocated"])
    monkeypatch.setattr(torch.cuda, "memory_reserved", lambda *args, **kwargs: state["reserved"])
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda *args, **kwargs: state["peak"])
    monkeypatch.setattr(torch.cuda, "mem_get_info", lambda *args, **kwargs: (state["free"], state["total"]))
    monkeypatch.setattr(torch.cuda, "device", lambda *args, **kwargs: nullcontext())

    def release():
        state["releases"] += 1

    monkeypatch.setattr(torch.cuda, "empty_cache", release)
    return state


def _raw_predictions(pipeline, fasta_path):
    identifiers, sequences, predictions = [], [], []
    for identifier, sequence in pipeline.read_fasta(fasta_path).items():
        identifiers.append([identifier])
        sequences.append([sequence])
        predictions.append([np.array(list(map(ord, sequence)))[None, :, None]])
    return identifiers, sequences, predictions


def _extract_value(tensors, sequence_index, position, residue):
    return {"residue": residue, "value": tensors[0][sequence_index, position, 0]}


def test_factory_loads_checkpoint_before_device_transfer_and_keeps_config(runtime, monkeypatch):
    pipeline, _, _ = runtime
    config = {"arch": {"args": {"embedding_pretrained": "unused-pretrained.pt"}}}
    events = []

    class Model:
        def to(self, device):
            events.append(("transfer", str(device)))
            return self

        def eval(self):
            events.append(("eval",))

    def construct(module, kind, actual_config):
        assert "embedding_pretrained" not in actual_config["arch"]["args"]
        events.append(("construct", kind))
        return Model()

    def checkpoint_predictor(model, checkpoint):
        events.append(("checkpoint", checkpoint))
        return SimpleNamespace(model=model)

    monkeypatch.setattr(pipeline, "load_config", lambda path: config)
    monkeypatch.setattr(pipeline.nsp3_main, "module_arch", object(), raising=False)
    monkeypatch.setattr(pipeline.nsp3_main, "get_instance", construct, raising=False)
    monkeypatch.setattr(pipeline.nsp3_main, "module_pred",
                        SimpleNamespace(SecondaryFeatures=checkpoint_predictor), raising=False)
    predictor = pipeline._load_nsp3_predictor("model.pth", "config.yml", torch.device("cuda:0"))
    assert isinstance(predictor.model, Model)
    assert events == [("construct", "arch"), ("checkpoint", "model.pth"), ("transfer", "cuda:0"), ("eval",)]
    assert config["arch"]["args"]["embedding_pretrained"] == "unused-pretrained.pt"


def test_one_predictor_is_reused_across_genes_and_batches(runtime, tmp_path, monkeypatch):
    pipeline, _, _ = runtime
    factories, calls, references = [], [], []

    class Predictor:
        def __call__(self, fasta_path):
            calls.append(list(pipeline.read_fasta(fasta_path)))
            assert not torch.is_grad_enabled()
            return _raw_predictions(pipeline, fasta_path)

    def factory(*args):
        factories.append(args)
        predictor = Predictor()
        references.append(weakref.ref(predictor))
        return predictor

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", factory)
    monkeypatch.setattr(pipeline, "extract_residue_predictions", _extract_value)
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        for gene in ("FIRST", "SECOND"):
            fasta = tmp_path / f"{gene}.fasta"
            fasta.write_text(f">{gene}_WT\nMACDE\n>{gene}_MUT\nMAGDE\n")
            results = pipeline.run_nsp3_prediction(
                fasta, "model.pth", "config.yml", batch_size=1, runtime=session)
            assert set(results) == {f"{gene}_WT", f"{gene}_MUT"}
        assert session.predictor is references[0]()
    assert len(factories) == 1
    assert len(calls) == 4
    assert session.predictor is None
    assert references[0]() is None
    assert torch.is_grad_enabled()


@pytest.mark.parametrize("raises", [False, True])
def test_prediction_disables_gradients_and_restores_context(runtime, tmp_path, monkeypatch, raises):
    pipeline, _, _ = runtime

    def predict(path):
        assert not torch.is_grad_enabled()
        if raises:
            raise ValueError("predictor failed")
        return "result"

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: predict)
    with torch.enable_grad():
        with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
            if raises:
                with pytest.raises(ValueError, match="predictor failed"):
                    session.predict(tmp_path / "input.fasta", [("gene", "M")])
            else:
                assert session.predict(tmp_path / "input.fasta", [("gene", "M")]) == "result"
            assert torch.is_grad_enabled()
        assert session.predictor is None


def test_cuda_calibration_uses_free_memory_lengths_and_limits(runtime, cuda_memory, monkeypatch, tmp_path):
    pipeline, _, _ = runtime
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: lambda path: "result")
    sequences = [(f"gene_{index}", "M" * 100) for index in range(30)]
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        assert session.choose_batch_size(sequences, 100) == 1
        session.predict(tmp_path / "unused", sequences[:1])
        assert session.last_peak_bytes == cuda_memory["peak"]
        assert session.choose_batch_size(sequences, 100) == 6
        assert session.choose_batch_size(sequences, 2) == 2
        longer = [(identifier, "M" * 202) for identifier, _ in sequences]
        assert session.choose_batch_size(longer, 100) == 1
        cuda_memory["free"] = 400 * 1024 * 1024
        assert session.choose_batch_size(sequences, 100) == 2
        cuda_memory["reserved"] += 600 * 1024 * 1024
        assert session.choose_batch_size(sequences, 100) == 6
        cuda_memory["free"] = 10000 * 1024 * 1024
        assert session.choose_batch_size(sequences, 100) == 25
    assert cuda_memory["releases"] == 1


def test_cpu_batches_obey_cli_and_upstream_batch_limits(runtime, monkeypatch):
    pipeline, _, _ = runtime
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: lambda path: None)
    sequences = [(str(index), "M") for index in range(40)]
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        assert session.choose_batch_size(sequences, 100) == 25
        assert session.choose_batch_size(sequences, 3) == 3
        assert session.choose_batch_size(sequences[:2], 100) == 2


def test_cpu_batch_stops_at_first_sequence_length_change(runtime, monkeypatch):
    pipeline, _, _ = runtime
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: lambda path: None)
    sequences = [("first", "M" * 20), ("second", "A" * 20),
                 ("third", "C" * 21), ("fourth", "D" * 20)]
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        assert session.choose_batch_size(sequences, 100) == 2
        assert session.choose_batch_size(sequences[1:], 100) == 1
        assert session.choose_batch_size(sequences[2:], 100) == 1


def test_cuda_batch_stops_at_first_sequence_length_change(runtime, cuda_memory, monkeypatch, tmp_path):
    pipeline, _, _ = runtime
    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: lambda path: "result")
    sequences = [("first", "M" * 20), ("second", "A" * 20),
                 ("third", "C" * 21), ("fourth", "D" * 20)]
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        session.predict(tmp_path / "unused", sequences[:1])
        assert session.choose_batch_size(sequences, 100) == 2
        assert session.choose_batch_size(sequences[1:], 100) == 1


def test_mixed_lengths_preserve_input_order_without_padding_neighbors(runtime, tmp_path, monkeypatch):
    pipeline, _, _ = runtime
    fasta = tmp_path / "genes.fasta"
    entries = [("first", "M" * 20), ("second", "A" * 20),
               ("third", "C" * 21), ("fourth", "D" * 20)]
    fasta.write_text("".join(f">{identifier}\n{sequence}\n" for identifier, sequence in entries))
    batches = []

    def predict(path):
        sequences = pipeline.read_fasta(path)
        assert len({len(sequence) for sequence in sequences.values()}) == 1
        batches.append(list(sequences))
        return _raw_predictions(pipeline, path)

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: predict)
    monkeypatch.setattr(pipeline, "extract_residue_predictions", _extract_value)
    result = pipeline.run_nsp3_prediction(fasta, "model.pth", "config.yml", batch_size=100)
    assert batches == [["first", "second"], ["third"], ["fourth"]]
    assert list(result) == [identifier for identifier, _ in entries]
    assert [len(result[identifier]) for identifier, _ in entries] == [20, 20, 21, 20]


def test_cuda_oom_halves_same_batch_without_missing_or_duplicate_residues(
        runtime, cuda_memory, tmp_path, monkeypatch):
    pipeline, _, _ = runtime
    fasta = tmp_path / "genes.fasta"
    fasta.write_text("".join(f">gene_{index}\nMACDE\n" for index in range(5)))
    attempts, completed = [], []

    def predict(path):
        identifiers = list(pipeline.read_fasta(path))
        attempts.append(identifiers)
        if len(identifiers) > 2:
            raise torch.cuda.OutOfMemoryError("simulated CUDA out of memory")
        completed.extend(identifiers)
        return _raw_predictions(pipeline, path)

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: predict)
    monkeypatch.setattr(pipeline, "extract_residue_predictions", _extract_value)
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        monkeypatch.setattr(session, "choose_batch_size", lambda sequences, limit: min(len(sequences), limit))
        results = pipeline.run_nsp3_prediction(
            fasta, "model.pth", "config.yml", batch_size=5, runtime=session)
    assert len(attempts[0]) == 5
    assert attempts[1] == attempts[0][:2]
    assert completed == [f"gene_{index}" for index in range(5)]
    assert set(results) == set(completed)
    for predictions in results.values():
        assert list(predictions) == list(range(1, 6))
        assert [row["value"] for row in predictions.values()] == list(map(ord, "MACDE"))
    assert cuda_memory["releases"] >= 2


@pytest.mark.parametrize("is_oom", [False, True])
def test_single_sequence_failure_is_not_retried(runtime, cuda_memory, tmp_path, monkeypatch, is_oom):
    pipeline, _, _ = runtime
    fasta = tmp_path / "gene.fasta"
    fasta.write_text(">gene\nMACDE\n")
    calls = []

    def predict(path):
        calls.append(path)
        if is_oom:
            raise torch.cuda.OutOfMemoryError("simulated CUDA out of memory")
        raise RuntimeError("unrelated model failure")

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: predict)
    with pipeline._NSP3Runtime("model.pth", "config.yml") as session:
        expected = "(?i)memory|GPU|CUDA" if is_oom else "unrelated model failure"
        with pytest.raises(RuntimeError, match=expected):
            pipeline.run_nsp3_prediction(fasta, "model.pth", "config.yml", batch_size=1, runtime=session)
    assert len(calls) == 1
    assert session.predictor is None


def test_owned_runtime_is_closed_when_prediction_raises(runtime, tmp_path, monkeypatch):
    pipeline, _, _ = runtime
    fasta = tmp_path / "gene.fasta"
    fasta.write_text(">gene\nMACDE\n")
    closed = []
    original_close = pipeline._NSP3Runtime.close

    def close(session):
        original_close(session)
        closed.append(session.predictor)

    def fail(path):
        raise RuntimeError("unrelated model failure")

    monkeypatch.setattr(pipeline, "_load_nsp3_predictor", lambda *args: fail)
    monkeypatch.setattr(pipeline._NSP3Runtime, "close", close)
    with pytest.raises(RuntimeError, match="unrelated model failure"):
        pipeline.run_nsp3_prediction(fasta, "model.pth", "config.yml")
    assert closed == [None]
