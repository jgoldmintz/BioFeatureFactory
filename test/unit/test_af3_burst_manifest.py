"""Strict manifest and cache behavior for AlphaFold3 burst jobs."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from biofeaturefactory.alphafold3.bin.burst_manifest import (
    COMMENT_HEADER,
    ManifestRow,
    count_pending,
    is_cache_complete,
    is_result_complete,
    main,
    publish_cache,
    begin_cache_stage, cache_coordination_dir,
    read_manifest,
    write_manifest,
)


def _write_result(
    cache_dir: Path,
    output_name: str,
    *,
    model_text: str = "data_model\n",
) -> Path:
    result_dir = cache_dir / output_name
    result_dir.mkdir(parents=True, exist_ok=True)
    (result_dir / f"{output_name}_model.cif").write_text(model_text)
    (result_dir / f"{output_name}_confidences.json").write_text(
        json.dumps({"atom_plddts": [90.0]})
    )
    (result_dir / f"{output_name}_summary_confidences.json").write_text(
        json.dumps({"ranking_score": 0.9})
    )
    return result_dir


def _row(cache_dir: Path, output_name: str, array_idx: int = 0) -> ManifestRow:
    return ManifestRow(
        array_idx=array_idx,
        input_hash=f"hash-{array_idx}",
        pkey="PAM-C1T",
        rbp_name="HNRNPC",
        allele="WT",
        window_idx=0,
        input_json_path=f"/inputs/{array_idx}.json",
        cache_dir=str(cache_dir),
        output_name=output_name,
    )


class TestManifestV2:
    def test_round_trip_includes_output_name_as_final_field(self, tmp_path):
        path = tmp_path / "manifest.tsv"
        row = _row(tmp_path / "cache" / "hash-0", "PAM_HNRNPC_WT")

        write_manifest([row], path)

        lines = path.read_text().splitlines()
        assert lines[0] == COMMENT_HEADER == "# burst_manifest_version=2"
        assert lines[1].split("\t")[-1] == "output_name"
        assert lines[2].split("\t")[-1] == "PAM_HNRNPC_WT"
        assert read_manifest(path) == [row]

    def test_count_pending_uses_each_rows_exact_output_name(self, tmp_path):
        complete_dir = tmp_path / "cache" / "complete"
        partial_dir = tmp_path / "cache" / "partial"
        _write_result(complete_dir, "expected_complete")
        _write_result(partial_dir / "unrelated", "expected_partial")
        manifest_path = tmp_path / "manifest.tsv"
        write_manifest(
            [
                _row(complete_dir, "expected_complete", array_idx=0),
                _row(partial_dir, "expected_partial", array_idx=1),
            ],
            manifest_path,
        )

        assert count_pending(manifest_path) == 1


class TestStrictCachePredicate:
    def test_accepts_exact_nonempty_parseable_sibling_triplet(self, tmp_path):
        cache_dir = tmp_path / "cache" / "hash"
        _write_result(cache_dir, "PAM_HNRNPC_WT")

        assert is_cache_complete(cache_dir, "PAM_HNRNPC_WT")

    def test_rejects_triplet_found_only_in_nested_sample_directory(self, tmp_path):
        cache_dir = tmp_path / "cache" / "hash"
        _write_result(
            cache_dir / "PAM_HNRNPC_WT" / "seed-1_sample-0",
            "PAM_HNRNPC_WT",
        )

        assert not is_cache_complete(cache_dir, "PAM_HNRNPC_WT")

    def test_rejects_wrong_basenames_in_expected_directory(self, tmp_path):
        cache_dir = tmp_path / "cache" / "hash"
        result_dir = cache_dir / "PAM_HNRNPC_WT"
        result_dir.mkdir(parents=True)
        (result_dir / "other_job_model.cif").write_text("data_model\n")
        (result_dir / "other_job_confidences.json").write_text("{}")
        (result_dir / "other_job_summary_confidences.json").write_text("{}")

        assert not is_cache_complete(cache_dir, "PAM_HNRNPC_WT")

    @pytest.mark.parametrize(
        "present",
        [
            (),
            ("model",),
            ("confidences",),
            ("summary",),
            ("model", "confidences"),
            ("model", "summary"),
            ("confidences", "summary"),
        ],
    )
    def test_rejects_every_incomplete_file_combination(self, tmp_path, present):
        output_name = "PAM_HNRNPC_WT"
        cache_dir = tmp_path / "cache" / "hash"
        result_dir = _write_result(cache_dir, output_name)
        paths = {
            "model": result_dir / f"{output_name}_model.cif",
            "confidences": result_dir / f"{output_name}_confidences.json",
            "summary": result_dir / f"{output_name}_summary_confidences.json",
        }
        for output_type, path in paths.items():
            if output_type not in present:
                path.unlink()

        assert not is_cache_complete(cache_dir, output_name)
        assert not is_result_complete(result_dir, output_name)

    @pytest.mark.parametrize(
        ("relative_name", "bad_content"),
        [
            ("PAM_HNRNPC_WT_model.cif", ""),
            ("PAM_HNRNPC_WT_confidences.json", ""),
            ("PAM_HNRNPC_WT_summary_confidences.json", ""),
            ("PAM_HNRNPC_WT_confidences.json", "{"),
            ("PAM_HNRNPC_WT_summary_confidences.json", "["),
        ],
    )
    def test_rejects_empty_or_malformed_primary_output(
        self, tmp_path, relative_name, bad_content
    ):
        cache_dir = tmp_path / "cache" / "hash"
        result_dir = _write_result(cache_dir, "PAM_HNRNPC_WT")
        (result_dir / relative_name).write_text(bad_content)

        assert not is_cache_complete(cache_dir, "PAM_HNRNPC_WT")

    def test_cli_uses_the_same_predicate(self, tmp_path):
        cache_dir = tmp_path / "cache" / "hash"
        assert main(["check", str(cache_dir), "PAM_HNRNPC_WT"]) == 1

        _write_result(cache_dir, "PAM_HNRNPC_WT")

        assert main(["check", str(cache_dir), "PAM_HNRNPC_WT"]) == 0

    def test_rejects_non_utf8_json(self, tmp_path):
        cache_dir = tmp_path / "cache" / "hash"
        result_dir = _write_result(cache_dir, "PAM_HNRNPC_WT")
        (result_dir / "PAM_HNRNPC_WT_confidences.json").write_bytes(b"\xff")

        assert not is_cache_complete(cache_dir, "PAM_HNRNPC_WT")


class TestCachePublication:
    def test_valid_source_is_atomically_published(self, tmp_path):
        cache_root = tmp_path / "cache"
        cache_dir = cache_root / "hash"
        source_dir = begin_cache_stage(cache_dir)
        _write_result(source_dir, "PAM_HNRNPC_WT")

        published = publish_cache(source_dir, cache_dir, "PAM_HNRNPC_WT")

        assert published == cache_dir
        assert not source_dir.exists()
        assert is_cache_complete(cache_dir, "PAM_HNRNPC_WT")
        assert (cache_coordination_dir(cache_root) / ".bff-af3-cache.lock").is_file()

    def test_existing_destination_is_quarantined_not_nested(self, tmp_path):
        cache_root = tmp_path / "cache"
        cache_dir = cache_root / "hash"
        source_dir = begin_cache_stage(cache_dir)
        old_result = _write_result(
            cache_dir, "PAM_HNRNPC_WT", model_text="old\n"
        )
        _write_result(source_dir, "PAM_HNRNPC_WT", model_text="new\n")

        publish_cache(source_dir, cache_dir, "PAM_HNRNPC_WT")

        new_model = cache_dir / "PAM_HNRNPC_WT" / "PAM_HNRNPC_WT_model.cif"
        quarantined = list(cache_root.glob("hash.stale-*"))
        assert new_model.read_text() == "new\n"
        assert len(quarantined) == 1
        assert (
            quarantined[0]
            / old_result.name
            / "PAM_HNRNPC_WT_model.cif"
        ).read_text() == "old\n"
        assert not (cache_dir / source_dir.name).exists()

    def test_partial_source_is_rejected_before_existing_cache_changes(
        self, tmp_path
    ):
        cache_root = tmp_path / "cache"
        source_dir = cache_root / ".staged"
        cache_dir = cache_root / "hash"
        old_result = _write_result(
            cache_dir, "PAM_HNRNPC_WT", model_text="old\n"
        )
        new_result = _write_result(source_dir, "PAM_HNRNPC_WT")
        (new_result / "PAM_HNRNPC_WT_model.cif").unlink()

        with pytest.raises(ValueError, match="source is incomplete"):
            publish_cache(source_dir, cache_dir, "PAM_HNRNPC_WT")

        assert source_dir.is_dir()
        assert (old_result / "PAM_HNRNPC_WT_model.cif").read_text() == "old\n"
        assert not list(cache_root.glob("hash.stale-*"))

    def test_publish_cli_rejects_partial_output(self, tmp_path, capsys):
        source_dir = tmp_path / "cache" / ".staged"
        source_dir.mkdir(parents=True)
        cache_dir = tmp_path / "cache" / "hash"

        exit_code = main([
            "publish",
            str(source_dir),
            str(cache_dir),
            "PAM_HNRNPC_WT",
        ])

        assert exit_code == 2
        assert "source is incomplete" in capsys.readouterr().err


def test_slurm_template_uses_shared_strict_cache_helper():
    template_path = (
        Path(__file__).parents[2]
        / "biofeaturefactory"
        / "alphafold3"
        / "bin"
        / "slurm_array.sh.tmpl"
    )
    template = template_path.read_text()

    assert "__PYTHON_BIN__" in template
    assert "__CACHE_HELPER__" in template
    assert 'check "$CACHE_DIR" "$OUTPUT_NAME"' in template
    assert 'publish "$TMP_OUT" "$CACHE_DIR" "$OUTPUT_NAME"' in template
    assert 'stage "$CACHE_DIR" "$CACHE_GENERATION"' in template
    assert "trap preserve_incomplete EXIT" in template
    assert 'mv -- "$TMP_OUT" "$partial_path"' in template
    assert 'rm -rf -- "$TMP_OUT"' not in template
    assert "find \"$CACHE_DIR\"" not in template
    assert 'mv "$TMP_OUT" "$CACHE_DIR"' not in template


def test_slurm_failure_preserves_staging_without_publishing(tmp_path):
    output_name = "bff_af3_testhash"
    cache_dir = tmp_path / "cache" / "testhash"
    input_json = tmp_path / "input.json"
    input_json.write_text(json.dumps({"name": output_name}))
    manifest = tmp_path / "manifest.tsv"
    row = _row(cache_dir, output_name)
    row.input_json_path = str(input_json)
    write_manifest([row], manifest)

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_docker = fake_bin / "docker"
    fake_docker.write_text("#!/bin/sh\nexit 23\n")
    fake_docker.chmod(0o755)

    template_path = (
        Path(__file__).parents[2]
        / "biofeaturefactory"
        / "alphafold3"
        / "bin"
        / "slurm_array.sh.tmpl"
    )
    cache_helper = template_path.with_name("burst_manifest.py")
    substitutions = {
        "__JOB_NAME__": "af3_test",
        "__PARTITION__": "gpu",
        "__TIME__": "00:10:00",
        "__MEM__": "8G",
        "__LOG_DIR__": str(tmp_path / "logs"),
        "__ARRAY_MAX__": "0",
        "__MANIFEST__": str(manifest),
        "__MODEL_DIR__": str(tmp_path / "models"),
        "__DOCKER_IMAGE__": "alphafold3:test",
        "__PYTHON_BIN__": sys.executable,
        "__CACHE_HELPER__": str(cache_helper),
        "__CACHE_GENERATION__": "0",
    }
    script = template_path.read_text()
    for placeholder, value in substitutions.items():
        script = script.replace(placeholder, value)
    script_path = tmp_path / "run.slurm"
    script_path.write_text(script)

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    env["SLURM_ARRAY_TASK_ID"] = "0"
    env["SLURM_JOB_ID"] = "123"
    env["SLURM_JOB_GPUS"] = "0"
    result = subprocess.run(
        ["bash", str(script_path)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 23
    incomplete = list(
        cache_coordination_dir(cache_dir.parent).glob('*.incomplete-123-0-*')
    )
    assert len(incomplete) == 1
    assert incomplete[0].is_dir()
    assert not cache_dir.exists()
    assert "Preserved incomplete AF3 output" in result.stderr
