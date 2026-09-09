from pathlib import Path
import sys
import types

import pytest

from biofeaturefactory.mutation_effects import adabmdca_pipeline as adabm
from biofeaturefactory.mutation_effects import mutEffects_controller as controller


@pytest.fixture
def controlled_training(tmp_path, monkeypatch):
    torch = pytest.importorskip("torch")

    class Dataset:
        def __init__(self, **kwargs):
            self.data = torch.tensor([[0, 1], [1, 0]])
            self.weights = torch.ones(2)

        def get_num_residues(self):
            return 2

        def get_num_states(self):
            return 2

        def get_effective_size(self):
            return 2

        def get_frequencies(self, **kwargs):
            return torch.full((2, 2), 0.5), None

    saved = []
    exports = {
        "adabmDCA": {}, "adabmDCA.dataset": {"DatasetDCA": Dataset},
        "adabmDCA.fasta": {"get_tokens": list},
        "adabmDCA.io": {"save_params": lambda *args: saved.append(args)},
        "adabmDCA.utils": {
            "get_device": torch.device, "get_dtype": lambda value: torch.float32,
            "init_parameters": lambda **kwargs: {
                "bias": torch.zeros(2, 2), "coupling_matrix": torch.zeros(2, 2, 2, 2),
            },
        },
    }
    for name, attributes in exports.items():
        module = types.ModuleType(name)
        module.__dict__.update(attributes)
        monkeypatch.setitem(sys.modules, name, module)
    monkeypatch.setattr(adabm, "_configure_task_resources", lambda device: None)

    def run(norms, quiet=True):
        observed = []

        def norm(tensor, *args, **kwargs):
            value = norms[len(observed) // 2] / 2
            observed.append(value)
            return torch.tensor(value, dtype=torch.float64)

        monkeypatch.setattr(torch.Tensor, "norm", norm)
        adabm.train_adabmdca_pseudolikelihood_in_process(
            "stub-dataset", "AC", str(tmp_path / "params"), device_str="cpu",
            nepochs=len(norms), tol=0.001, patience=2, check_every=1, quiet=quiet,
        )
        assert len(saved) == 1
        return len(observed) // 2

    return run, saved


@pytest.mark.parametrize("norms,epochs", [
    ([10, 10.001, 10.002, 10.003], 4),
    ([10, 10, 11, 11, 11, 11], 5),
    ([10, 10, 10, 10], 3),
    ([10, 8, 6, 4], 4),
    ([10, 9.999, 9.998, 9.997], 3),
    ([0, 0, 0, 0], 3),
])
def test_actual_pseudolikelihood_loop_plateau_policy(controlled_training, norms, epochs):
    run, saved = controlled_training
    assert run(norms) == epochs


@pytest.mark.parametrize("nonfinite", [float("nan"), float("inf")])
def test_nonfinite_gradient_is_not_saved_as_trained(controlled_training, nonfinite):
    run, saved = controlled_training
    with pytest.raises(RuntimeError, match="[Nn]on.?finite"):
        run([10, nonfinite, nonfinite])
    assert saved == []


def test_zero_gradient_convergence_diagnostic_is_safe(controlled_training, capsys):
    run, saved = controlled_training
    assert run([0, 0, 0, 0], quiet=False) == 3
    assert "epochs run: 3 of 4" in capsys.readouterr().out


@pytest.fixture
def database_args(tmp_path, monkeypatch):
    source = tmp_path / "GENE.fasta"
    source.write_text(">ORF\nATGGCTTAA\n")
    mutations = tmp_path / "GENE.csv"
    mutations.write_text("mutant\nG4A\n")
    database = tmp_path / "Bio DBs"
    database.mkdir()
    monkeypatch.setattr(sys, "argv", [
        "mutEffects_controller", "--fasta", str(source), "--mutations", str(mutations),
        "--output", str(tmp_path / "results"), "--db-root", str(database),
        "--evmutation-only",
    ])
    args = controller.parse_args()
    return args, controller.build_manifest(["GENE"], args)


def test_gzip_only_database_rejected_before_nextflow(database_args):
    args, manifest = database_args
    (args.db_root / "uniref90.fasta.gz").write_bytes(b"gzip-placeholder")
    with pytest.raises(SystemExit, match="uncompressed.*jackhmmer|jackhmmer.*uncompressed"):
        controller.validate_db_coverage(["GENE"], manifest, args)
    command = controller.build_nextflow_cmd(args, "manifest.json")
    assert "--uniref90_db" not in command


def test_plain_database_validation_and_forwarding_agree(database_args):
    args, manifest = database_args
    database = args.db_root / "uniref90.fasta"
    database.write_text(">protein\nMA\n")
    controller.validate_db_coverage(["GENE"], manifest, args)
    command = controller.build_nextflow_cmd(args, "manifest.json")
    assert command[command.index("--uniref90_db") + 1] == str(database.resolve())


def test_empty_plain_database_is_not_silently_overwritten(database_args):
    args, manifest = database_args
    (args.db_root / "uniref90.fasta").touch()
    (args.db_root / "uniref90.fasta.gz").write_bytes(b"gzip-placeholder")
    with pytest.raises(SystemExit, match="Move aside.*first"):
        controller.validate_db_coverage(["GENE"], manifest, args)
    assert (args.db_root / "uniref90.fasta").stat().st_size == 0


def test_missing_database_reports_preflight_error(database_args):
    args, manifest = database_args
    with pytest.raises(SystemExit, match="uniref90"):
        controller.validate_db_coverage(["GENE"], manifest, args)


def test_prebuilt_msa_needs_no_uniref_database(database_args):
    args, manifest = database_args
    manifest["msa"] = ["GENE"]
    controller.validate_db_coverage(["GENE"], manifest, args)
    command = controller.build_nextflow_cmd(args, "manifest.json")
    assert "--uniref90_db" not in command
