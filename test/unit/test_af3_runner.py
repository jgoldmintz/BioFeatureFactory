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
"""AlphaFold3 local runner resume and native directory batching."""

import hashlib
import json
import subprocess
import threading
from pathlib import Path

import pytest

import biofeaturefactory.alphafold3.bin.af3_runner as af3_runner_module
from biofeaturefactory.alphafold3.bin.af3_runner import (
    AF3Input,
    AF3Runner,
    AF3RunnerConfig,
)


def _input(name: str, msa: str = ">query\nMKT\n") -> AF3Input:
    return AF3Input(
        name=name,
        rna_sequence="ACGU",
        protein_sequence="MKT",
        rna_chain_id="R",
        protein_chain_id="P",
        protein_msa=msa,
    )


def _write_primary_outputs(
    result_dir: Path,
    output_name: str,
    missing: str = "",
) -> None:
    result_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "model": (result_dir / f"{output_name}_model.cif", "data_test\n#\n"),
        "confidences": (
            result_dir / f"{output_name}_confidences.json",
            json.dumps({"atom_plddts": [90.0]}),
        ),
        "summary": (
            result_dir / f"{output_name}_summary_confidences.json",
            json.dumps({"ranking_score": 0.9}),
        ),
    }
    for output_type, (path, content) in outputs.items():
        if output_type != missing:
            path.write_text(content)


def _write_existing_job(
    runner: AF3Runner,
    job_id: str,
    af3_input: AF3Input,
    missing: str = "",
    write_provenance: bool = True,
) -> Path:
    job_dir = runner.output_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / f"{af3_input.name}.json").write_text(
        json.dumps(af3_input.to_json_dict(), indent=2)
    )
    result_dir = job_dir / "output" / af3_input.sanitised_name()
    _write_primary_outputs(result_dir, af3_input.sanitised_name(), missing)
    if write_provenance:
        provenance_path = runner._result_provenance_path(result_dir)
        provenance_path.write_text(json.dumps(
            _expected_provenance(runner, af3_input),
            indent=2,
            sort_keys=True,
        ))
    return result_dir


def _expected_provenance(runner: AF3Runner, af3_input: AF3Input) -> dict:
    model_path = runner._model_dir / "af3.bin"
    canonical_input = json.dumps(
        af3_input.to_json_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return {
        "runtime": {
            "schema_version": runner._PROVENANCE_SCHEMA_VERSION,
            "docker_image_id": "sha256:test-image",
            "model_files": [{
                "name": model_path.name,
                "size": model_path.stat().st_size,
                "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
            }],
            "inference_mode": "norun_data_pipeline",
        },
        "input_sha256": hashlib.sha256(canonical_input.encode()).hexdigest(),
    }


def _mounted_path(command: list[str], container_path: str) -> Path:
    marker = f":{container_path}"
    for argument in command:
        if marker in argument:
            host_path, suffix = argument.split(marker, 1)
            if suffix in ("", ":ro"):
                return Path(host_path)
    raise AssertionError(f"No mount found for {container_path}: {command}")


class FakeAF3:
    def __init__(self, plan=None):
        self.plan = plan
        self.calls = []
        self.image_inspections = 0

    def __call__(self, command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            self.image_inspections += 1
            return subprocess.CompletedProcess(
                command,
                0,
                stdout="sha256:test-image\n",
                stderr="",
            )

        input_dir = _mounted_path(command, "/root/af_input")
        output_dir = _mounted_path(command, "/root/af_output")
        payloads = [
            json.loads(path.read_text())
            for path in sorted(input_dir.glob("*.json"))
        ]
        names = [payload["name"] for payload in payloads]
        precreated_output_names = sorted(
            path.name for path in output_dir.iterdir() if path.is_dir()
        )
        nonempty_precreated_outputs = sorted(
            path.name
            for path in output_dir.iterdir()
            if path.is_dir() and any(path.iterdir())
        )
        self.calls.append({
            "command": list(command),
            "names": names,
            "precreated_output_names": precreated_output_names,
            "nonempty_precreated_outputs": nonempty_precreated_outputs,
            "timeout": kwargs.get("timeout"),
        })

        if self.plan is None:
            completed_names = set(names)
            return_code = 0
            log_text = ""
        else:
            completed_names, return_code, log_text = self.plan(
                len(self.calls) - 1, names
            )

        for payload in payloads:
            if payload["name"] not in completed_names:
                continue
            af3_input = AF3Input(
                name=payload["name"],
                rna_sequence=payload["sequences"][0]["rna"]["sequence"],
                protein_sequence=payload["sequences"][1]["protein"]["sequence"],
            )
            output_name = af3_input.sanitised_name()
            _write_primary_outputs(
                output_dir / output_name,
                output_name,
            )

        log_handle = kwargs.get("stdout")
        if log_handle is not None and log_text:
            log_handle.write(log_text)
            log_handle.flush()
        return subprocess.CompletedProcess(command, return_code)


@pytest.fixture
def runner_factory(tmp_path, monkeypatch):
    runners = []
    monkeypatch.setattr(AF3Runner, '_preflight_runtime', lambda self: None)

    def create(
        *, batch_size=16, jax_cache_dir=None, adopt_legacy_results=False,
        max_gpus=1, runner_root=None
    ):
        if runner_root is None:
            runner_root = tmp_path / f"runner-{len(runners)}"
        runner_root = Path(runner_root)
        model_dir = runner_root / "models"
        model_dir.mkdir(parents=True, exist_ok=True)
        (model_dir / "af3.bin").write_bytes(b"weights")
        runner = AF3Runner(AF3RunnerConfig(
            output_base_dir=str(runner_root / "af3_runs"),
            model_dir=str(model_dir),
            max_gpus=max_gpus,
            batch_size=batch_size,
            adopt_legacy_results=adopt_legacy_results,
            jax_cache_dir=(str(jax_cache_dir) if jax_cache_dir else None),
            timeout_per_job=7,
        ))
        runners.append(runner)
        return runner

    yield create

    for runner in runners:
        runner.shutdown()


class TestResume:
    def test_exact_complete_existing_job_is_resumed_without_inference(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        af3_input = _input("NPM1_variant_HNRNPC_WT")
        expected_result = _write_existing_job(
            runner, "NPM1_variant_HNRNPC_WT", af3_input
        )
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        future = runner.submit_job_async(
            af3_input, job_id="NPM1_variant_HNRNPC_WT"
        )
        job = future.result(timeout=1)

        assert future.done()
        assert job.status == "completed"
        assert job.result_path == expected_result
        assert runner.get_job_status(job.job_id) == "completed"
        assert fake_af3.calls == []
        assert fake_af3.image_inspections == 1

    def test_explicit_legacy_adoption_writes_provenance_without_inference(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(adopt_legacy_results=True)
        af3_input = _input("NPM1_variant_HNRNPC_MUT")
        result_dir = _write_existing_job(
            runner,
            "NPM1_variant_HNRNPC_MUT",
            af3_input,
            write_provenance=False,
        )
        provenance_path = runner._result_provenance_path(result_dir)
        assert not provenance_path.exists()

        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        job = runner.submit_job_async(
            af3_input, job_id="NPM1_variant_HNRNPC_MUT"
        ).result(timeout=1)

        assert job.status == "completed"
        assert job.result_path == result_dir
        assert json.loads(provenance_path.read_text()) == (
            runner._expected_job_provenance(af3_input)
        )
        assert fake_af3.calls == []
        assert fake_af3.image_inspections == 1

    def test_legacy_result_without_adoption_is_refused_without_modification(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        af3_input = _input("NPM1_variant_ELAVL1_MUT")
        result_dir = _write_existing_job(
            runner,
            "NPM1_variant_ELAVL1_MUT",
            af3_input,
            write_provenance=False,
        )
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with pytest.raises(RuntimeError, match="--adopt-legacy-results"):
            runner.submit_job_async(
                af3_input, job_id="NPM1_variant_ELAVL1_MUT"
            )

        assert result_dir.is_dir()
        assert not runner._result_provenance_path(result_dir).exists()
        assert fake_af3.calls == []

    def test_provenance_mismatch_is_quarantined_and_rerun(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        af3_input = _input("PAM_variant_FMR1_WT")
        result_dir = _write_existing_job(
            runner, "PAM_variant_FMR1_WT", af3_input
        )
        provenance_path = runner._result_provenance_path(result_dir)
        provenance = json.loads(provenance_path.read_text())
        provenance["runtime"]["docker_image_id"] = "sha256:old-image"
        provenance_path.write_text(json.dumps(provenance))
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                af3_input, job_id="PAM_variant_FMR1_WT"
            )
        job = future.result(timeout=5)

        job_dir = runner.output_dir / "PAM_variant_FMR1_WT"
        assert len(fake_af3.calls) == 1
        assert job.status == "completed"
        assert list(job_dir.glob("output.stale-*"))
        assert json.loads(
            runner._result_provenance_path(job.result_path).read_text()
        ) == _expected_provenance(runner, af3_input)

    @pytest.mark.parametrize("missing", ["model", "confidences", "summary"])
    def test_partial_existing_job_is_rerun(
        self, runner_factory, monkeypatch, missing
    ):
        runner = runner_factory()
        af3_input = _input("NPM1_variant_ELAVL1_WT")
        old_result = _write_existing_job(
            runner, "NPM1_variant_ELAVL1_WT", af3_input, missing=missing
        )
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                af3_input, job_id="NPM1_variant_ELAVL1_WT"
            )
        job = future.result(timeout=5)

        assert len(fake_af3.calls) == 1
        assert fake_af3.calls[0]["names"] == [af3_input.name]
        assert job.status == "completed"
        assert job.result_path == old_result
        assert runner._is_complete_result(old_result, af3_input.sanitised_name())
        assert list(old_result.parent.glob(f"{old_result.name}.stale-*"))

    def test_sample_outputs_do_not_make_partial_job_resumable(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        af3_input = _input("NPM1_variant_YTHDF1_MUT")
        result_dir = _write_existing_job(
            runner, "NPM1_variant_YTHDF1_MUT", af3_input, missing="model"
        )
        for path in result_dir.iterdir():
            path.unlink()
        sample_name = f"{af3_input.sanitised_name()}_seed-1_sample-0"
        _write_primary_outputs(result_dir / "seed-1_sample-0", sample_name)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                af3_input, job_id="NPM1_variant_YTHDF1_MUT"
            )
        job = future.result(timeout=5)

        assert len(fake_af3.calls) == 1
        assert job.status == "completed"
        assert job.result_path == result_dir
        assert runner._is_complete_result(result_dir, af3_input.sanitised_name())

    def test_non_utf8_confidences_make_existing_job_non_resumable(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        af3_input = _input("NPM1_variant_YTHDF1_WT")
        result_dir = _write_existing_job(
            runner, "NPM1_variant_YTHDF1_WT", af3_input
        )
        confidences = (
            result_dir
            / f"{af3_input.sanitised_name()}_confidences.json"
        )
        confidences.write_bytes(b"\xff")
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                af3_input, job_id="NPM1_variant_YTHDF1_WT"
            )
        job = future.result(timeout=5)

        assert len(fake_af3.calls) == 1
        assert job.status == "completed"
        assert job.result_path == result_dir
        assert list(result_dir.parent.glob(f"{result_dir.name}.stale-*"))

    def test_changed_msa_is_rerun(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory()
        old_input = _input("PAM_variant_DGCR8_WT", msa=">query\nMKT\n>old\nMRT\n")
        new_input = _input("PAM_variant_DGCR8_WT", msa=">query\nMKT\n>new\nMST\n")
        assert old_input.get_hash() != new_input.get_hash()
        _write_existing_job(runner, "PAM_variant_DGCR8_WT", old_input)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                new_input, job_id="PAM_variant_DGCR8_WT"
            )
        job = future.result(timeout=5)

        job_dir = runner.output_dir / "PAM_variant_DGCR8_WT"
        assert len(fake_af3.calls) == 1
        assert job.status == "completed"
        assert json.loads((job_dir / f"{new_input.name}.json").read_text()) == (
            new_input.to_json_dict()
        )
        assert list(job_dir.glob("output.stale-*"))


class TestAF3InputHash:
    def test_name_does_not_change_content_identity(self):
        first = _input("first name")
        second = _input("second name")

        assert first.get_hash() == second.get_hash()

    @pytest.mark.parametrize(
        "changed_input",
        [
            AF3Input("job", "UCGU", "MKT", protein_msa=">query\nMKT\n"),
            AF3Input("job", "ACGU", "MRT", protein_msa=">query\nMKT\n"),
            AF3Input(
                "job", "ACGU", "MKT", rna_chain_id="R",
                protein_msa=">query\nMKT\n",
            ),
            AF3Input(
                "job", "ACGU", "MKT", protein_chain_id="P",
                protein_msa=">query\nMKT\n",
            ),
            AF3Input("job", "ACGU", "MKT", protein_msa=">query\nMKT\n>hit\nMRT\n"),
        ],
    )
    def test_each_prediction_field_changes_content_identity(self, changed_input):
        baseline = AF3Input(
            "job", "ACGU", "MKT", protein_msa=">query\nMKT\n"
        )

        assert changed_input.get_hash() != baseline.get_hash()

    def test_hash_is_stable_full_sha256(self):
        af3_input = AF3Input(
            "ignored name",
            "ACTU",
            "MKT",
            rna_chain_id="R",
            protein_chain_id="P",
            protein_msa=">query\nMKT\n>hit\nMRT\n",
        )

        assert af3_input.get_hash() == (
            "092f23094b746265c9b28a78faa0b5bf82c01d8cd6fe7b1105c442dafa742b24"
        )


class TestBatchExecution:
    def test_duplicate_inflight_job_id_is_rejected(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            future = runner.submit_job_async(
                _input("first_input"), job_id="shared-job-id"
            )
            with pytest.raises(RuntimeError, match="already in flight"):
                runner.submit_job_async(
                    _input("second_input"), job_id="shared-job-id"
                )

        job = future.result(timeout=5)
        assert job.status == "completed"
        assert fake_af3.calls[0]["names"] == ["first_input"]

    def test_same_output_job_id_is_locked_across_runner_instances(
        self, runner_factory, tmp_path, monkeypatch
    ):
        shared_root = tmp_path / "shared-runner-root"
        first_runner = runner_factory(
            batch_size=8,
            runner_root=shared_root,
        )
        second_runner = runner_factory(
            batch_size=8,
            runner_root=shared_root,
        )
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with first_runner.batch_submissions():
            future = first_runner.submit_job_async(
                _input("first_process_input"),
                job_id="shared-filesystem-job-id",
            )
            with pytest.raises(RuntimeError, match="another process"):
                second_runner.submit_job_async(
                    _input("second_process_input"),
                    job_id="shared-filesystem-job-id",
                )

        job = future.result(timeout=5)
        assert job.status == "completed"
        assert fake_af3.calls[0]["names"] == ["first_process_input"]

    def test_cancelled_deferred_job_is_not_dispatched_and_releases_lock(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)

        with runner.batch_submissions():
            cancelled = runner.submit_job_async(
                _input("cancelled_input"), job_id="reusable-job-id"
            )
            retained = runner.submit_job_async(
                _input("retained_input"), job_id="retained-job-id"
            )
            assert cancelled.cancel()

        assert retained.result(timeout=5).status == "completed"
        assert cancelled.cancelled()
        assert fake_af3.calls[0]["names"] == ["retained_input"]

        with runner.batch_submissions():
            replacement = runner.submit_job_async(
                _input("replacement_input"), job_id="reusable-job-id"
            )
        assert replacement.result(timeout=5).status == "completed"

    def test_cancelled_queued_executor_job_releases_lock(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=1)
        fake_af3 = FakeAF3()
        first_run_started = threading.Event()
        release_first_run = threading.Event()
        inference_count = 0

        def blocking_fake(command, **kwargs):
            nonlocal inference_count
            if command[:2] == ["docker", "run"]:
                inference_count += 1
                if inference_count == 1:
                    first_run_started.set()
                    assert release_first_run.wait(timeout=5)
            return fake_af3(command, **kwargs)

        monkeypatch.setattr(
            af3_runner_module.subprocess, "run", blocking_fake
        )
        first = runner.submit_job_async(
            _input("blocking_input"), job_id="blocking-job"
        )
        assert first_run_started.wait(timeout=5)
        cancelled = runner.submit_job_async(
            _input("queued_input"), job_id="queued-job"
        )
        assert cancelled.cancel()

        replacement = runner.submit_job_async(
            _input("replacement_input"), job_id="queued-job"
        )
        release_first_run.set()
        assert first.result(timeout=5).status == "completed"
        assert replacement.result(timeout=5).status == "completed"
        assert cancelled.cancelled()

    def test_jobs_are_balanced_across_available_gpus(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8, max_gpus=2)
        fake_af3 = FakeAF3()
        first_two_runs = threading.Barrier(2)

        def concurrent_fake(command, **kwargs):
            if command[:2] == ["docker", "run"]:
                first_two_runs.wait(timeout=2)
            return fake_af3(command, **kwargs)

        monkeypatch.setattr(
            af3_runner_module.subprocess, "run", concurrent_fake
        )
        inputs = [_input(f"gpu_job_{index}") for index in range(6)]

        with runner.batch_submissions():
            futures = [
                runner.submit_job_async(af3_input, job_id=af3_input.name)
                for af3_input in inputs
            ]
        jobs = [future.result(timeout=5) for future in futures]

        assert sorted(len(call["names"]) for call in fake_af3.calls) == [3, 3]
        gpu_assignments = {
            call["command"][call["command"].index("--gpus") + 1]
            for call in fake_af3.calls
        }
        assert gpu_assignments == {"device=0", "device=1"}
        assert all(job.status == "completed" for job in jobs)

    def test_directory_batches_split_and_mount_persistent_jax_cache(
        self, runner_factory, tmp_path, monkeypatch
    ):
        jax_cache_dir = tmp_path / "persistent-jax-cache"
        runner = runner_factory(batch_size=2, jax_cache_dir=jax_cache_dir)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        inputs = [_input(f"batch_job_{index}") for index in range(5)]

        with runner.batch_submissions():
            futures = [
                runner.submit_job_async(af3_input, job_id=af3_input.name)
                for af3_input in inputs
            ]
        jobs = [future.result(timeout=5) for future in futures]

        assert [len(call["names"]) for call in fake_af3.calls] == [2, 2, 1]
        assert [call["timeout"] for call in fake_af3.calls] == [14, 14, 7]
        for call in fake_af3.calls:
            command = call["command"]
            assert "--input_dir=/root/af_input" in command
            assert not any(arg.startswith("--json_path=") for arg in command)
            assert "--jax_compilation_cache_dir=/root/jax_cache" in command
            assert f"{jax_cache_dir.resolve()}:/root/jax_cache" in command
        assert all(job.status == "completed" for job in jobs)

    def test_result_directories_are_empty_and_host_precreated_before_docker(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        inputs = [_input("Alpha Job"), _input("Beta(Job)")]

        with runner.batch_submissions():
            futures = [
                runner.submit_job_async(af3_input, job_id=f"job-{index}")
                for index, af3_input in enumerate(inputs)
            ]
        jobs = [future.result(timeout=5) for future in futures]

        assert fake_af3.calls[0]["precreated_output_names"] == sorted(
            af3_input.sanitised_name() for af3_input in inputs
        )
        assert fake_af3.calls[0]["nonempty_precreated_outputs"] == []
        assert all(job.status == "completed" for job in jobs)

    def test_each_job_receives_its_exact_isolated_result_path(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)
        fake_af3 = FakeAF3()
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        inputs = [_input("Alpha Job"), _input("Beta(Job)")]

        with runner.batch_submissions():
            futures = [
                runner.submit_job_async(af3_input, job_id=f"job-{index}")
                for index, af3_input in enumerate(inputs)
            ]
        jobs = [future.result(timeout=5) for future in futures]

        expected_paths = [
            runner.output_dir
            / f"job-{index}"
            / "output"
            / af3_input.sanitised_name()
            for index, af3_input in enumerate(inputs)
        ]
        assert [job.result_path for job in jobs] == expected_paths
        assert len({job.result_path for job in jobs}) == len(jobs)
        for index, (job, af3_input) in enumerate(zip(jobs, inputs)):
            output_name = af3_input.sanitised_name()
            other_output_name = inputs[1 - index].sanitised_name()
            assert (job.result_path / f"{output_name}_model.cif").is_file()
            assert not list(job.result_path.glob(f"{other_output_name}_*"))

    def test_partial_batch_preserves_success_and_retries_only_missing_jobs(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)

        def plan(call_index, names):
            if call_index == 0:
                return {names[0]}, 1, "one biological input failed\n"
            return set(names), 0, ""

        fake_af3 = FakeAF3(plan=plan)
        monkeypatch.setattr(af3_runner_module.subprocess, "run", fake_af3)
        inputs = [_input(f"fallback_job_{index}") for index in range(3)]

        with runner.batch_submissions():
            futures = [
                runner.submit_job_async(af3_input, job_id=af3_input.name)
                for af3_input in inputs
            ]
        jobs = [future.result(timeout=5) for future in futures]

        assert fake_af3.calls[0]["names"] == [item.name for item in inputs]
        assert [call["names"] for call in fake_af3.calls[1:]] == [
            [inputs[1].name],
            [inputs[2].name],
        ]
        assert inputs[0].name not in {
            name for call in fake_af3.calls[1:] for name in call["names"]
        }
        assert all(job.status == "completed" for job in jobs)
        assert all(job.result_path.is_dir() for job in jobs)
        assert not list(runner.output_dir.glob("*/output/*.incomplete-*"))

    def test_timeout_removes_named_container_when_cidfile_is_absent(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)
        run_commands = []
        cleanup_commands = []

        def timeout_fake(command, **kwargs):
            if command[:3] == ["docker", "image", "inspect"]:
                return subprocess.CompletedProcess(
                    command,
                    0,
                    stdout="sha256:test-image\n",
                    stderr="",
                )
            if command[:3] == ["docker", "rm", "-f"]:
                cleanup_commands.append(list(command))
                return subprocess.CompletedProcess(
                    command, 0, stdout="removed\n", stderr=""
                )
            run_commands.append(list(command))
            raise subprocess.TimeoutExpired(command, kwargs["timeout"])

        monkeypatch.setattr(af3_runner_module.subprocess, "run", timeout_fake)
        af3_input = _input("timeout_job")

        with runner.batch_submissions():
            future = runner.submit_job_async(
                af3_input, job_id=af3_input.name
            )
        job = future.result(timeout=5)

        container_name = run_commands[0][run_commands[0].index("--name") + 1]
        assert cleanup_commands == [["docker", "rm", "-f", container_name]]
        assert job.status == "timeout"
        assert job.result_path is None
        assert runner._docker_failed is False
        assert not list(runner.output_dir.glob("*/output/*.incomplete-*"))
        assert not list(runner._work_dir.glob("batch-*"))

    def test_daemon_error_does_not_claim_timed_out_container_is_gone(
        self, tmp_path, monkeypatch
    ):
        commands = []

        def unavailable_daemon(command, **kwargs):
            commands.append(list(command))
            return subprocess.CompletedProcess(
                command,
                1,
                stdout="",
                stderr="Cannot connect to the Docker daemon",
            )

        monkeypatch.setattr(
            af3_runner_module.subprocess, "run", unavailable_daemon
        )

        assert not AF3Runner._stop_container(
            tmp_path / "missing.cid", "bff-af3-timeout"
        )
        assert commands == [
            ["docker", "rm", "-f", "bff-af3-timeout"],
            ["docker", "inspect", "bff-af3-timeout"],
        ]

    def test_read_only_incomplete_tree_is_preserved_before_temp_cleanup(
        self, runner_factory, monkeypatch
    ):
        runner = runner_factory(batch_size=8)

        def incomplete_fake(command, **kwargs):
            if command[:3] == ["docker", "image", "inspect"]:
                return subprocess.CompletedProcess(
                    command,
                    0,
                    stdout="sha256:test-image\n",
                    stderr="",
                )

            input_dir = _mounted_path(command, "/root/af_input")
            output_dir = _mounted_path(command, "/root/af_output")
            payload_path = next(input_dir.glob("*.json"))
            payload = json.loads(payload_path.read_text())
            af3_input = AF3Input(
                name=payload["name"],
                rna_sequence=payload["sequences"][0]["rna"]["sequence"],
                protein_sequence=(
                    payload["sequences"][1]["protein"]["sequence"]
                ),
            )
            output_name = af3_input.sanitised_name()
            staged_result = output_dir / output_name
            sample_name = f"{output_name}_seed-1_sample-0"
            sample_dir = staged_result / "seed-1_sample-0"
            _write_primary_outputs(sample_dir, sample_name)
            sample_dir.chmod(0o555)
            staged_result.chmod(0o555)
            return subprocess.CompletedProcess(command, 1)

        monkeypatch.setattr(
            af3_runner_module.subprocess,
            "run",
            incomplete_fake,
        )
        af3_input = _input("read_only_partial_job")
        preserved_paths = []
        try:
            with runner.batch_submissions():
                future = runner.submit_job_async(
                    af3_input,
                    job_id=af3_input.name,
                )
            job = future.result(timeout=5)

            preserved_paths = list(
                (runner.output_dir / af3_input.name / "output").glob(
                    f"{af3_input.sanitised_name()}.incomplete-*"
                )
            )
            assert job.status == "failed"
            assert len(preserved_paths) == 1
            sample_dir = preserved_paths[0] / "seed-1_sample-0"
            assert sample_dir.is_dir()
            assert list(sample_dir.glob("*_model.cif"))
            assert not list(runner._work_dir.glob("batch-*"))
        finally:
            staged_paths = list(
                runner._work_dir.glob(
                    f"batch-*/outputs/{af3_input.sanitised_name()}"
                )
            )
            for output_path in preserved_paths + staged_paths:
                for path in output_path.rglob("*"):
                    if path.is_dir():
                        path.chmod(0o755)
                output_path.chmod(0o755)
