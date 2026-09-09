"""Package-version gate and curl-only RBP downloads without installing tools."""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
from packaging.specifiers import SpecifierSet


REPO = Path(__file__).resolve().parents[2]
BOOTSTRAP = REPO / "scripts" / "bootstrap.sh"
BUILD_DB = REPO / "scripts" / "build_db.sh"


@pytest.mark.parametrize("version", ["3.9", "3.10", "3.11", "3.12", "3.13", "3.14"])
@pytest.mark.parametrize("phase", ["env-phase", "db-phase"])
def test_bootstrap_gate_matches_package_policy(tmp_path, version, phase):
    source = BOOTSTRAP.read_text()
    preflight = source.split('ROOT_DIR="$SCRIPT_DIR"', 1)[0]
    script = tmp_path / "bootstrap.sh"
    script.write_text(preflight + "\nprintf 'PREFLIGHT_COMPLETE\\n'\n")
    environment_dir = tmp_path / "custom-env"
    (environment_dir / "bin").mkdir(parents=True)
    interpreter = environment_dir / "bin" / "python"
    interpreter.write_text(f"#!/bin/sh\nprintf '{version}\\n'\n")
    interpreter.chmod(0o755)
    environment = os.environ.copy()
    environment["CONDA_PREFIX"] = str(environment_dir)
    result = subprocess.run(["bash", str(script), phase], env=environment, capture_output=True, text=True)
    requirement = re.search(r'^requires-python = "([^"]+)"',
                            (REPO / "pyproject.toml").read_text(), re.MULTILINE).group(1)
    allowed = phase == "db-phase" or version in SpecifierSet(requirement)
    assert result.returncode == (0 if allowed else 1), result.stdout + result.stderr
    assert ("PREFLIGHT_COMPLETE" in result.stdout) == allowed


@pytest.mark.parametrize("download_fails", [False, True])
def test_curl_only_active_rbp_download(tmp_path, download_fails):
    database = tmp_path / "db"
    rbp_dir = database / "AF3" / "RBP_db"
    rbp_dir.mkdir(parents=True)
    (rbp_dir / "rbp_uniprot_ids.txt").write_text("P12345\n")
    binary_dir = tmp_path / "bin"
    binary_dir.mkdir()
    curl = binary_dir / "curl"
    curl.write_text(f'''#!{sys.executable}
import pathlib
import sys
if "-I" in sys.argv:
    print("200")
    raise SystemExit(0)
assert "--fail" in sys.argv
output = pathlib.Path(sys.argv[sys.argv.index("--output") + 1])
output.write_text(">P12345\\nMKT\\n")
raise SystemExit({22 if download_fails else 0})
''')
    curl.chmod(0o755)
    environment = os.environ.copy()
    environment.update({"PATH": f"{binary_dir}:/usr/bin:/bin", "DB_ROOT": str(database),
                        "SKIP_REFSEQ": "1", "SKIP_IDMAPPING": "1", "SKIP_UNIREF90": "1",
                        "SKIP_MIRNA": "1", "SKIP_COCOPUTS": "1", "SKIP_AF3RBP": "0",
                        "AF3_DOWNLOAD_RBP_MSAS": "1", "AF3_RBP_MSA_ARCHIVE_URL": "",
                        "POSTAR3_TXT_URL": "", "AF3_MSA_VERSION": "v6"})
    result = subprocess.run(["/bin/bash", str(BUILD_DB)], env=environment,
                            capture_output=True, text=True, cwd=tmp_path)
    output = rbp_dir / "msa" / "AF-P12345-F1-msa_v6.a3m"
    assert result.returncode == (22 if download_fails else 0), result.stdout + result.stderr
    assert output.exists() is not download_fails
    if not download_fails:
        assert output.read_text() == ">P12345\nMKT\n"
        assert not Path(str(output) + ".part").exists()
