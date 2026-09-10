"""Checkpoint reconstruction and strict loading without the vendored NSP3 package."""

from argparse import Namespace
import gc
import importlib.util
from pathlib import Path
import sys
import types
import weakref

import esm
import pytest
import torch
from esm.modules import ESM1bLayerNorm


@pytest.fixture
def runtime(monkeypatch):
    class BasePredict(torch.nn.Module):
        pass

    class ESM1bEmbedding(torch.nn.Module):
        pass

    exports = {
        "nsp3": {},
        "nsp3.main": {},
        "nsp3.cli": {"load_config": lambda path: {}},
        "nsp3.base": {},
        "nsp3.base.base_predict": {"BasePredict": BasePredict},
        "nsp3.embeddings": {},
        "nsp3.embeddings.esm1b": {"ESM1bEmbedding": ESM1bEmbedding},
    }
    for name, attributes in exports.items():
        module = types.ModuleType(name)
        module.__path__ = []
        module.__dict__.update(attributes)
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(sys, "path", list(sys.path))
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    filename = Path(__file__).resolve().parents[2] / "biofeaturefactory/NetSurfP3/netsurfp3_pipeline.py"
    spec = importlib.util.spec_from_file_location("nsp3_checkpoint_test_wrapper", filename)
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)

    def build_model(normalization, dtype=torch.float32, root_embedding=False):
        arguments = Namespace(
            arch="roberta_large", layers=2, embed_dim=8, ffn_embed_dim=16,
            attention_heads=2, max_positions=32, token_dropout=False,
            emb_layer_norm_before=normalization,
        )
        embedding = ESM1bEmbedding()
        embedding.model = esm.ProteinBertModel(arguments, esm.Alphabet.from_architecture("roberta_large"))
        embedding.max_embedding = 1024
        if root_embedding:
            return embedding.to(dtype=dtype)
        model = torch.nn.Module()
        model.embedding = embedding
        model.head = torch.nn.Linear(8, 3)
        return model.to(dtype=dtype)

    def load_model(model, checkpoint_path):
        predictor = BasePredict(model, str(checkpoint_path))
        return predictor.model

    return pipeline, build_model, load_model


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
@pytest.mark.parametrize("existing_normalization", [False, True])
def test_checkpoint_normalization_is_loaded_and_used(runtime, tmp_path, dtype, existing_normalization):
    pipeline, build_model, load_model = runtime
    trained = build_model(True, dtype).eval()
    with torch.no_grad():
        trained.embedding.model.emb_layer_norm_before.weight.copy_(torch.linspace(0.5, 1.5, 8))
        trained.embedding.model.emb_layer_norm_before.bias.copy_(torch.linspace(-0.2, 0.3, 8))
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": trained.state_dict()}, checkpoint)
    reconstructed = build_model(existing_normalization, dtype)
    original_layer = reconstructed.embedding.model.emb_layer_norm_before
    loaded = load_model(reconstructed, checkpoint)
    layer = loaded.embedding.model.emb_layer_norm_before
    assert isinstance(layer, ESM1bLayerNorm)
    if existing_normalization:
        assert layer is original_layer
    assert layer.weight.device == loaded.embedding.model.embed_tokens.weight.device
    assert layer.weight.dtype == dtype
    assert not layer.training
    for name, tensor in trained.state_dict().items():
        torch.testing.assert_close(loaded.state_dict()[name], tensor)
    assert loaded.embedding.forward.__func__ is pipeline._bounded_esm_forward
    tokens = torch.tensor([[0, 5, 6, 7, 2]])
    with torch.no_grad():
        expected = trained.embedding.model(tokens, repr_layers=[2])["representations"][2]
        actual = loaded.embedding.model(tokens, repr_layers=[2])["representations"][2]
    torch.testing.assert_close(actual, expected)


def test_checkpoint_without_normalization_does_not_add_it(runtime, tmp_path):
    _, build_model, load_model = runtime
    trained = build_model(False)
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": trained.state_dict()}, checkpoint)
    loaded = load_model(build_model(False), checkpoint)
    assert loaded.embedding.model.emb_layer_norm_before is None
    assert set(loaded.state_dict()) == set(trained.state_dict())


def test_root_embedding_normalization_key_has_no_leading_dot(runtime, tmp_path):
    _, build_model, load_model = runtime
    trained = build_model(True, root_embedding=True)
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": trained.state_dict()}, checkpoint)
    loaded = load_model(build_model(False, root_embedding=True), checkpoint)
    assert isinstance(loaded.model.emb_layer_norm_before, ESM1bLayerNorm)
    torch.testing.assert_close(
        loaded.model.emb_layer_norm_before.weight,
        trained.model.emb_layer_norm_before.weight,
    )


@pytest.mark.parametrize("defect", ["normalization_weight", "normalization_bias", "head", "unexpected", "shape"])
def test_incompatible_checkpoint_is_rejected(runtime, tmp_path, defect):
    _, build_model, load_model = runtime
    state = build_model(True).state_dict()
    if defect == "normalization_weight":
        state.pop("embedding.model.emb_layer_norm_before.weight")
    elif defect == "normalization_bias":
        state.pop("embedding.model.emb_layer_norm_before.bias")
    elif defect == "head":
        state.pop("head.weight")
    elif defect == "unexpected":
        state["unknown.weight"] = torch.ones(1)
    else:
        state["embedding.model.emb_layer_norm_before.weight"] = torch.ones(7)
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": state}, checkpoint)
    with pytest.raises(RuntimeError, match="Missing key|Unexpected key|size mismatch"):
        load_model(build_model(False), checkpoint)


def test_checkpoint_deserializes_on_cpu_when_cuda_is_available(runtime, tmp_path, monkeypatch):
    _, build_model, load_model = runtime
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": build_model(True).state_dict()}, checkpoint)
    original_load = torch.load
    options = []

    def tracked_load(*args, **kwargs):
        options.append(kwargs)
        return original_load(*args, **kwargs)

    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch, "load", tracked_load)
    loaded = load_model(build_model(False), checkpoint)
    assert len(options) == 1
    assert str(options[0]["map_location"]) == "cpu"
    assert options[0]["weights_only"] is False
    assert next(loaded.parameters()).device.type == "cpu"


def test_embedding_override_does_not_require_cyclic_gc(runtime, tmp_path):
    pipeline, build_model, load_model = runtime
    checkpoint = tmp_path / "model.pth"
    torch.save({"state_dict": build_model(True).state_dict()}, checkpoint)
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        loaded = load_model(build_model(False), checkpoint)
        assert loaded.embedding.forward.__func__ is pipeline._bounded_esm_forward
        model_reference = weakref.ref(loaded)
        embedding_reference = weakref.ref(loaded.embedding)
        esm_reference = weakref.ref(loaded.embedding.model)
        weight_reference = weakref.ref(loaded.embedding.model.embed_tokens.weight)
        del loaded
        assert model_reference() is None
        assert embedding_reference() is None
        assert esm_reference() is None
        assert weight_reference() is None
    finally:
        if was_enabled:
            gc.enable()
