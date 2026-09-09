# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""AlphaFold3 burst orchestration and ingest regressions."""

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from biofeaturefactory.alphafold3 import burst
from biofeaturefactory.alphafold3.alphafold3_pipeline import (
    AlphaFold3Pipeline,
    MutationContext,
)
from biofeaturefactory.alphafold3.bin.af3_runner import AF3Input
from biofeaturefactory.alphafold3.bin.binding_metrics import (
    BindingMetrics,
    RnaEditSpan,
    ThresholdConfig,
    compute_delta_metrics,
)
from biofeaturefactory.alphafold3.burst import BurstInput, _failed_flag


class TestFailedFlag:
    def test_encodes_the_reason(self):
        assert _failed_flag("boom") == "FAILED:boom"

    def test_prefix_is_machine_greppable(self):
        # Ingest distinguishes a failed shard from a pending one by this prefix,
        # so it must stay a fixed literal rather than free text.
        assert _failed_flag("anything").startswith("FAILED:")


def _submit_args(tmp_path: Path) -> SimpleNamespace:
    return SimpleNamespace(
        output=str(tmp_path / "out"),
        clear_cache=False,
        postar_db="postar.bed.gz",
        rbp_mapping="rbp.tsv",
        rbp_sequences=None,
        msa_dir="msa",
        slurm_log_dir=None,
        slurm_job_name="af3_test",
        slurm_partition="gpu",
        slurm_time="00:10:00",
        slurm_mem="8G",
        model_dir=str(tmp_path / "models"),
        docker_image="alphafold3:test",
        no_submit=True,
        array_throttle=4,
    )


def _mock_burst_services(monkeypatch) -> None:
    monkeypatch.setattr(burst, "_acquire_burst_lock", lambda _burst_dir: object())
    monkeypatch.setattr(burst, "_ingest_warn_in_flight_slurm", lambda _burst_dir: None)
    monkeypatch.setattr(burst, "POSTAR3Database", lambda _path: object())
    monkeypatch.setattr(burst, "RBPSequenceMapper", lambda **_kwargs: object())


def test_submit_assets_are_immutable_per_generation(tmp_path, monkeypatch, capsys):
    args = _submit_args(tmp_path)
    af3_input = AF3Input(
        name="NPM1 A/B:WT",
        rna_sequence="ACGU",
        protein_sequence="MPEP",
    )
    burst_input = BurstInput(
        gene="NPM1",
        pkey="NPM1-A1G",
        rbp_name="RBP1",
        allele="WT",
        window_idx=0,
        af3_input=af3_input,
    )
    alias_input = BurstInput(
        gene="PAM",
        pkey="PAM-C2T",
        rbp_name="RBP2",
        allele="MUT",
        window_idx=3,
        af3_input=AF3Input(
            name="different metadata name",
            rna_sequence="ACGU",
            protein_sequence="MPEP",
        ),
    )
    assert burst_input.pkey != alias_input.pkey
    assert burst_input.af3_input.name != alias_input.af3_input.name
    assert burst_input.input_hash == alias_input.input_hash
    cache_checks = []
    submitted = []

    _mock_burst_services(monkeypatch)
    monkeypatch.setattr(
        burst,
        "iterate_inputs",
        lambda _args, _db, _mapper, skipped=None: [burst_input, alias_input],
    )

    def cache_complete(cache_dir, output_name):
        cache_checks.append((Path(cache_dir), output_name))
        return False

    def run_sbatch(command, **_kwargs):
        submitted.append(command)
        return subprocess.CompletedProcess(
            command, returncode=0, stdout="Submitted batch job 123\n", stderr="",
        )

    monkeypatch.setattr(burst, "is_cache_complete", cache_complete)
    monkeypatch.setattr(burst.subprocess, "run", run_sbatch)

    assert burst.cmd_submit(args) == 0
    first_stderr = capsys.readouterr().err
    submissions_root = Path(args.output) / ".burst" / "submissions"
    first_submission, = list(submissions_root.iterdir())
    first_manifest = first_submission / "manifest.tsv"
    first_script = first_submission / "run.slurm"
    first_snapshot = (first_manifest.read_bytes(), first_script.read_bytes())

    output_name = burst._cache_output_name(burst_input.input_hash)
    first_json = first_submission / "inputs" / f"{burst_input.input_hash}.json"
    assert first_json.is_file()
    assert json.loads(first_json.read_text())["name"] == output_name
    assert first_submission.joinpath("logs").is_dir()
    assert str(first_script) in first_stderr

    args.no_submit = False
    assert burst.cmd_submit(args) == 0

    submission_dirs = set(submissions_root.iterdir())
    assert len(submission_dirs) == 2
    second_submission, = submission_dirs - {first_submission}
    assert first_snapshot == (first_manifest.read_bytes(), first_script.read_bytes())
    for submission_dir in submission_dirs:
        assert submission_dir.joinpath("inputs").is_dir()
        assert submission_dir.joinpath("manifest.tsv").is_file()
        assert submission_dir.joinpath("run.slurm").is_file()
        assert submission_dir.joinpath("logs").is_dir()
        input_jsons = list(submission_dir.joinpath("inputs").glob("*.json"))
        assert len(input_jsons) == 1
        assert json.loads(input_jsons[0].read_text())["name"] == output_name
        assert len(submission_dir.joinpath("manifest.tsv").read_text().splitlines()) == 3

    assert submitted == [[
        "sbatch", "--array=0-0%4", str((second_submission / "run.slurm").resolve()),
    ]]
    last_submit = json.loads(
        (Path(args.output) / ".burst" / "last_submit.json").read_text()
    )
    assert last_submit["generation_id"] == second_submission.name
    assert last_submit["submission_dir"] == str(second_submission.resolve())
    assert last_submit["manifest_path"] == str(
        (second_submission / "manifest.tsv").resolve()
    )
    assert last_submit["script"] == str((second_submission / "run.slurm").resolve())
    assert last_submit["log_dir"] == str((second_submission / "logs").resolve())

    assert burst._cache_output_name(alias_input.input_hash) == output_name
    assert cache_checks == [
        (Path(args.output) / ".cache" / "af3" / burst_input.input_hash, output_name),
        (Path(args.output) / ".cache" / "af3" / burst_input.input_hash, output_name),
    ]
    manifest_text = (second_submission / "manifest.tsv").read_text()
    assert "output_name" in manifest_text.splitlines()[1].split("\t")
    assert manifest_text.splitlines()[-1].endswith(f"\t{output_name}")
    script_text = (second_submission / "run.slurm").read_text()
    assert sys.executable in script_text
    assert str(
        (Path(burst.__file__).parent / "bin" / "burst_manifest.py").resolve()
    ) in script_text


def test_ingest_preserves_the_span_of_each_first_site_set(tmp_path, monkeypatch):
    output_dir = tmp_path / "out"
    (output_dir / ".cache" / "af3").mkdir(parents=True)
    args = SimpleNamespace(
        output=str(output_dir),
        postar_db="postar.bed.gz",
        rbp_mapping="rbp.tsv",
        rbp_sequences=None,
        msa_dir="msa",
    )
    window_zero_span = RnaEditSpan(
        offset=1, ref_len=1, alt_len=2, wt_len=8, mut_len=9,
    )
    window_one_span = RnaEditSpan(
        offset=5, ref_len=1, alt_len=2, wt_len=8, mut_len=9,
    )

    def make_input(allele, window_idx, span, rna_sequence):
        return BurstInput(
            gene="NPM1",
            pkey="NPM1-A2AG",
            rbp_name="RBP1",
            allele=allele,
            window_idx=window_idx,
            af3_input=AF3Input(
                name=f"NPM1_RBP1_{allele}_w{window_idx}",
                rna_sequence=rna_sequence,
                protein_sequence="MPEP",
            ),
            edit_span=span,
        )

    wt_zero = make_input("WT", 0, window_zero_span, "AAAA")
    mut_zero = make_input("MUT", 0, window_zero_span, "AAAGA")
    wt_one = make_input("WT", 1, window_one_span, "CCCC")
    mut_one = make_input("MUT", 1, window_one_span, "CCCGC")
    inputs = [wt_zero, mut_zero, wt_one, mut_one]
    wt_sites = [object()]
    mut_sites = [object()]
    aggregate = SimpleNamespace(
        contact_frequency_rna={1: 1.0},
        contact_frequency_protein={2: 0.5},
    )
    parsed = {
        burst._cache_output_name(wt_zero.input_hash): (object(), wt_sites, aggregate),
        burst._cache_output_name(mut_zero.input_hash): (object(), [], None),
        burst._cache_output_name(wt_one.input_hash): (object(), [], None),
        burst._cache_output_name(mut_one.input_hash): (object(), mut_sites, aggregate),
    }
    parsed_names = []
    formatted = []
    qc_calls = []

    _mock_burst_services(monkeypatch)
    monkeypatch.setattr(
        burst, "iterate_inputs", lambda _args, _db, _mapper, skipped=None: inputs,
    )

    def parse_cache(_cache_dir, output_name, rbp_name, threshold_config):
        parsed_names.append(output_name)
        return parsed[output_name]

    def format_sites(
        _pkey, _rbp_name, allele, sites,
        _contact_frequency_rna=None, _contact_frequency_protein=None,
        *, edit_span, window_idx,
    ):
        formatted.append((allele, sites, edit_span))
        return [{"allele": allele}]

    delta = SimpleNamespace(wt_metrics=object(), mut_metrics=object(), n_windows=None)
    monkeypatch.setattr(burst, "_parse_cache_entry", parse_cache)
    monkeypatch.setattr(
        burst,
        "_aggregate_across_windows",
        lambda metrics, _rbp_name, _config: next(m for m in metrics if m is not None),
    )
    monkeypatch.setattr(burst, "compute_delta_metrics", lambda **_kwargs: delta)
    monkeypatch.setattr(burst, "compute_window_delta", lambda *_args: delta)
    monkeypatch.setattr(burst, "format_sites_rows", format_sites)
    monkeypatch.setattr(burst, "aggregate_mutation_summary", lambda _deltas: {})
    monkeypatch.setattr(burst, "format_events_rows", lambda _pkey, _deltas: [])
    monkeypatch.setattr(burst, "write_tsv", lambda *_args, **_kwargs: None)

    def qc_flag(deltas):
        qc_calls.append(deltas)
        return "PASS"

    monkeypatch.setattr(burst, "qc_flag_for_deltas", qc_flag)

    assert burst.cmd_ingest(args) == 0
    assert parsed_names == [
        burst._cache_output_name(item.input_hash) for item in inputs
    ]
    assert formatted == [
        ("WT", wt_sites, window_zero_span),
        ("MUT", mut_sites, window_one_span),
    ]
    assert qc_calls == [[delta]]


def test_ingest_rejects_sites_without_their_exact_span(tmp_path, monkeypatch):
    output_dir = tmp_path / "out"
    (output_dir / ".cache" / "af3").mkdir(parents=True)
    args = SimpleNamespace(
        output=str(output_dir),
        postar_db="postar.bed.gz",
        rbp_mapping="rbp.tsv",
        rbp_sequences=None,
        msa_dir="msa",
    )
    wt_input = BurstInput(
        gene="NPM1",
        pkey="NPM1-A1G",
        rbp_name="RBP1",
        allele="WT",
        window_idx=0,
        af3_input=AF3Input("NPM1_WT", "AAAA", "MPEP"),
        edit_span=None,
    )

    _mock_burst_services(monkeypatch)
    monkeypatch.setattr(
        burst, "iterate_inputs", lambda _args, _db, _mapper, skipped=None: [wt_input],
    )
    monkeypatch.setattr(
        burst,
        "_parse_cache_entry",
        lambda _cache_dir, output_name, rbp_name, threshold_config: (
            object(), [object()], None,
        ),
    )

    with pytest.raises(RuntimeError, match="exact edit_span is unavailable"):
        burst.cmd_ingest(args)


def _metrics(rbp_name, *, confident):
    return BindingMetrics(
        rbp_name=rbp_name,
        chain_pair_pae_min=5.0 if confident else 99.0,
        interface_contacts=5,
        interface_plddt_rna=70.0 if confident else 1.0,
        interface_plddt_protein=70.0 if confident else 1.0,
        has_binding=True,
    )


def test_local_finalizer_writes_confident_counts_and_partial_qc():
    pipeline = object.__new__(AlphaFold3Pipeline)
    pipeline.summary_rows = []
    pipeline.events_rows = []
    context = MutationContext(
        pkey="NPM1-A1G",
        gene="NPM1",
        mutation="A1G",
        nt_pos=1,
        ref="A",
        alt="G",
        transcript_seq="AAAA",
        wt_rna_window="AAAA",
        mut_rna_window="GAAA",
        window_center=0,
    )
    deltas = [
        compute_delta_metrics(
            "LOWCONF", _metrics("LOWCONF", confident=False),
            _metrics("LOWCONF", confident=False),
        ),
        compute_delta_metrics("FAILED", None, None),
    ]

    pipeline._finalize_mutation_results(context, deltas)

    summary = pipeline.summary_rows[0]
    assert summary["n_rbps_tested"] == 2
    assert summary["n_rbps_binding_wt"] == 0
    assert summary["n_rbps_binding_mut"] == 0
    assert summary["qc_flags"] == "PARTIAL"
    assert [row["cls"] for row in pipeline.events_rows] == [
        "no_binding", "incomplete",
    ]


def test_burst_ingest_writes_confident_counts_and_partial_qc(
    tmp_path, monkeypatch,
):
    output_dir = tmp_path / "out"
    (output_dir / ".cache" / "af3").mkdir(parents=True)
    args = SimpleNamespace(
        output=str(output_dir),
        postar_db="postar.bed.gz",
        rbp_mapping="rbp.tsv",
        rbp_sequences=None,
        msa_dir="msa",
    )

    def make_input(rbp_name, allele, protein_sequence):
        return BurstInput(
            gene="NPM1",
            pkey="NPM1-A1G",
            rbp_name=rbp_name,
            allele=allele,
            window_idx=0,
            af3_input=AF3Input(
                name=f"NPM1_{rbp_name}_{allele}",
                rna_sequence="AAAA" if allele == "WT" else "GAAA",
                protein_sequence=protein_sequence,
            ),
        )

    inputs = [
        make_input("LOWCONF", "WT", "MPEP"),
        make_input("LOWCONF", "MUT", "MPEP"),
        make_input("FAILED", "WT", "MFAIL"),
        make_input("FAILED", "MUT", "MFAIL"),
    ]
    parsed = {}
    for item in inputs:
        output_name = burst._cache_output_name(item.input_hash)
        metrics = (
            _metrics(item.rbp_name, confident=False)
            if item.rbp_name == "LOWCONF" else None
        )
        parsed[output_name] = (metrics, [], None)
    written = {}

    _mock_burst_services(monkeypatch)
    monkeypatch.setattr(
        burst, "iterate_inputs", lambda _args, _db, _mapper, skipped=None: inputs,
    )
    monkeypatch.setattr(
        burst,
        "_parse_cache_entry",
        lambda _cache_dir, output_name, rbp_name, threshold_config: parsed[
            output_name
        ],
    )
    monkeypatch.setattr(
        burst,
        "write_tsv",
        lambda rows, path, _fieldnames: written.setdefault(Path(path).name, rows),
    )

    assert burst.cmd_ingest(args) == 1

    summary = written["NPM1.tsv"][0]
    assert summary["n_rbps_tested"] == 2
    assert summary["n_rbps_binding_wt"] == 0
    assert summary["n_rbps_binding_mut"] == 0
    assert summary["qc_flags"] == "PARTIAL"
    assert [row["cls"] for row in written["NPM1.events.tsv"]] == [
        "no_binding", "incomplete",
    ]


@pytest.mark.parametrize("vcf_layout", ["single", "flat", "nested"])
def test_iterate_inputs_resolves_exact_gene_vcf(
    vcf_layout, tmp_path, monkeypatch,
):
    fasta_path = tmp_path / "NPM1.fasta"
    mutations_path = tmp_path / "NPM1.csv"
    fasta_path.write_text(">NPM1\nAAAA\n")
    mutations_path.write_text("mutation\nA1G\n")

    vcf_root = tmp_path / "vcfs"
    if vcf_layout == "nested":
        vcf_dir = vcf_root / "NPM1" / "vcf"
    else:
        vcf_dir = vcf_root
    vcf_dir.mkdir(parents=True)
    exact_vcf = vcf_dir / "NPM1.vcf"
    exact_vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\nchr5\t1\t.\tA\tG\n")
    (vcf_dir / "NPM1.annotated.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\nchrDERIVATIVE\t1\t.\tA\tG\n"
    )
    vcf_arg = exact_vcf if vcf_layout == "single" else vcf_root

    seen_chromosomes = []
    monkeypatch.setattr(burst, "read_fasta", lambda _path: {"NPM1": "AAAA"})
    monkeypatch.setattr(
        burst, "trim_muts", lambda _path, _validation_log, _gene: ["A1G"],
    )

    def iterate_mutation(**kwargs):
        seen_chromosomes.append(kwargs["chrom"])
        return iter(())

    monkeypatch.setattr(burst, "_iterate_mutation", iterate_mutation)
    args = SimpleNamespace(
        fasta=str(fasta_path),
        mutations=str(mutations_path),
        vcf=str(vcf_arg),
        chromosome_mapping=None,
        chrom="fallback",
        validation_log=None,
    )

    assert list(burst.iterate_inputs(args, object(), object())) == []
    assert seen_chromosomes == ["chr5"]


def test_parse_cache_entry_uses_output_scoped_completion(tmp_path, monkeypatch):
    cache_checks = []
    parsed_paths = []

    def cache_complete(cache_dir, output_name):
        cache_checks.append((Path(cache_dir), output_name))
        return True

    def parse_samples(result_dir):
        parsed_paths.append(Path(result_dir))
        return []

    monkeypatch.setattr(burst, "is_cache_complete", cache_complete)
    monkeypatch.setattr(burst, "parse_all_samples", parse_samples)
    cache_dir = tmp_path / "cache"

    assert burst._parse_cache_entry(
        cache_dir,
        output_name="NPM1_RBP1_WT",
        rbp_name="RBP1",
        threshold_config=ThresholdConfig(),
    ) == (None, None, None)
    assert cache_checks == [(cache_dir, "NPM1_RBP1_WT")]
    assert parsed_paths == [cache_dir / "NPM1_RBP1_WT"]
