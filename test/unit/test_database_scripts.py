"""Database bootstrap destination and exclusion behavior."""

import os
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
BOOTSTRAP = REPO_ROOT / "scripts" / "bootstrap.sh"
BUILD_DB = REPO_ROOT / "scripts" / "build_db.sh"
BUILD_REFSEQ_PATHS = REPO_ROOT / "scripts" / "build_refseq_ftp_paths.awk"


def _clean_env(tmp_path: Path) -> dict[str, str]:
    env = os.environ.copy()
    for name in (
        "AF3_DOWNLOAD_RBP_MSAS",
        "AF3_RBP_MSA_ARCHIVE_URL",
        "AF3_MSA_URL_TEMPLATE",
        "AF3_MSA_VERSION",
        "ARIA_CONN",
        "ARIA_SPLIT",
        "BFF_BIO_DBS",
        "BIO_DBS_DIR",
        "CAPTURE_FILE",
        "CONDA_PREFIX",
        "DB_ROOT",
        "EXTRA_TAXON_GROUPS",
        "MIRBASE_HSA_URL",
        "MIRBASE_REFRESH",
        "NETWORK_LOG",
        "PARALLEL_JOBS",
        "POSTAR3_TXT_URL",
        "PYTHON",
        "SKIP_COCOPUTS",
        "SKIP_IDMAPPING",
        "SKIP_UNIREF90",
        "TAXON_GROUP",
        "UNIREF90_URL",
    ):
        env.pop(name, None)
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    env["HOME"] = str(home)
    return env


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(0o755)


def _bootstrap_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    repo = tmp_path / "checkout" / "BioFeatureFactory"
    scripts = repo / "scripts"
    scripts.mkdir(parents=True)
    shutil.copy2(BOOTSTRAP, scripts / "bootstrap.sh")
    capture = tmp_path / "builder-env.txt"
    _write_executable(
        scripts / "build_db.sh",
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\n%s\n%s\n%s\n' \
  "$DB_ROOT" \
  "${SKIP_UNIREF90:-0}" \
  "${SKIP_IDMAPPING:-0}" \
  "${SKIP_COCOPUTS:-0}" > "$CAPTURE_FILE"
""",
    )
    return repo, scripts / "bootstrap.sh", capture


def _builder_fixture(tmp_path: Path) -> tuple[Path, Path]:
    repo = tmp_path / "nested" / "BioFeatureFactory"
    scripts = repo / "scripts"
    scripts.mkdir(parents=True)
    shutil.copy2(BUILD_DB, scripts / "build_db.sh")
    shutil.copy2(BUILD_REFSEQ_PATHS, scripts / "build_refseq_ftp_paths.awk")
    return repo, scripts / "build_db.sh"


def _install_failing_downloaders(tmp_path: Path, env: dict[str, str]) -> Path:
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    script = """#!/usr/bin/env bash
printf '%s\n' "$PWD" >> "$NETWORK_LOG"
exit 97
"""
    for name in ("aria2c", "curl", "wget"):
        _write_executable(fake_bin / name, script)
    network_log = tmp_path / "network.log"
    env["NETWORK_LOG"] = str(network_log)
    env["PATH"] = f"{fake_bin}:{env['PATH']}"
    return network_log


def _seed_minimal_builder_database(database_root: Path) -> None:
    assembly = "GCF_TEST"
    database_root.mkdir()
    (database_root / "assembly_summary_refseq.txt").write_text(
        "# assembly_accession\tversion_status\trefseq_category\tassembly_level\t"
        "ftp_path\tgroup\n"
        f"{assembly}\tlatest\treference genome\tComplete Genome\t"
        f"https://example.invalid/{assembly}\tfixture_group\n"
    )
    assembly_root = database_root / "refseq_assemblies" / assembly
    assembly_root.mkdir(parents=True)
    for suffix, content in (
        ("protein.faa", ">protein\nM\n"),
        ("cds_from_genomic.fna", ">cds\nATG\n"),
        ("feature_table.txt", "feature\n"),
    ):
        (assembly_root / f"{assembly}_{suffix}").write_text(content)
    (database_root / "mature_hsa.fasta").write_text(">hsa-test\nAUG\n")


def test_bootstrap_delegates_to_repo_local_database_root(tmp_path):
    repo, bootstrap, capture = _bootstrap_fixture(tmp_path)
    caller = tmp_path / "caller"
    caller.mkdir()
    env = _clean_env(tmp_path)
    env["CAPTURE_FILE"] = str(capture)

    result = subprocess.run(
        ["bash", str(bootstrap), "db-phase", "--exclude-htslib"],
        cwd=caller,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    captured = capture.read_text().splitlines()
    assert Path(captured[0]).resolve() == (repo / "Bio_DBs").resolve()
    assert captured[1:] == ["0", "0", "0"]
    assert not (repo / "scripts" / "_downloads").exists()
    assert "DEFER to scripts/build_db.sh" in result.stdout


def test_bootstrap_resolves_relative_override_and_forwards_exclusions(tmp_path):
    _, bootstrap, capture = _bootstrap_fixture(tmp_path)
    caller = tmp_path / "caller with space"
    caller.mkdir()
    env = _clean_env(tmp_path)
    env["CAPTURE_FILE"] = str(capture)
    env["BFF_BIO_DBS"] = str(tmp_path / "environment databases")
    env["DB_ROOT"] = str(tmp_path / "inherited database root")

    result = subprocess.run(
        [
            "bash",
            str(bootstrap),
            "db-phase",
            "--exclude-htslib",
            "--bio-dbs",
            "relative databases",
            "--exclude-uniref90",
            "--exclude-idmapping",
            "--exclude-cocoputs",
        ],
        cwd=caller,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    captured = capture.read_text().splitlines()
    assert Path(captured[0]).resolve() == (caller / "relative databases").resolve()
    assert captured[1:] == ["1", "1", "1"]


def test_build_db_default_does_not_reuse_ancestor_database_root(tmp_path):
    repo, build_db = _builder_fixture(tmp_path)
    ancestor_db = tmp_path / "Bio_DBs"
    ancestor_db.mkdir()
    sentinel = ancestor_db / "sentinel"
    sentinel.write_text("untouched")
    env = _clean_env(tmp_path)
    network_log = _install_failing_downloaders(tmp_path, env)

    result = subprocess.run(
        ["bash", str(build_db)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 97
    assert network_log.read_text().splitlines() == [
        str((repo / "Bio_DBs").resolve())
    ]
    assert (repo / "Bio_DBs").is_dir()
    assert sentinel.read_text() == "untouched"


def test_build_db_relative_override_and_skip_flags_are_effective(tmp_path):
    _, build_db = _builder_fixture(tmp_path)
    caller = tmp_path / "caller with space"
    caller.mkdir()
    database_root = caller / "relative databases"
    _seed_minimal_builder_database(database_root)
    env = _clean_env(tmp_path)
    network_log = _install_failing_downloaders(tmp_path, env)
    env.update(
        {
            "AF3_DOWNLOAD_RBP_MSAS": "0",
            "DB_ROOT": "relative databases",
            "MIRBASE_REFRESH": "0",
            "SKIP_COCOPUTS": "1",
            "SKIP_IDMAPPING": "1",
            "SKIP_UNIREF90": "1",
            "TAXON_GROUP": "fixture_group",
            "EXTRA_TAXON_GROUPS": "fixture_group",
        }
    )

    result = subprocess.run(
        ["bash", str(build_db)],
        cwd=caller,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert not network_log.exists()
    assert f"DB root: {database_root.resolve()}" in result.stdout
    assert "SKIP (SKIP_IDMAPPING=1)" in result.stdout
    assert "SKIP (SKIP_UNIREF90=1)" in result.stdout
    assert not (database_root / "idmapping.dat.gz").exists()
    assert not (database_root / "uniref90.fasta.gz").exists()


def test_bootstrap_rejects_database_phase_with_every_effective_step_excluded(
    tmp_path,
):
    _, bootstrap, capture = _bootstrap_fixture(tmp_path)
    env = _clean_env(tmp_path)
    env["CAPTURE_FILE"] = str(capture)

    result = subprocess.run(
        [
            "bash",
            str(bootstrap),
            "db-phase",
            "--exclude-build-db",
            "--exclude-htslib",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert "every step in the selected phase(s) is excluded" in result.stderr
    assert not capture.exists()
