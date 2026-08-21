#!/usr/bin/env bash
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

set -euo pipefail

# BioFeatureFactory dependency bootstrapper
# - Reads implicit setup needs across module READMEs.
# - Automates git/aria2c/curl/tar/make/conda steps where possible.
# - Emits checks/instructions for licensed/manual dependencies.
#
# Usage:
#   ./bootstrap.sh                 # Full install (all steps)
#   ./bootstrap.sh env-only        # Only pip install
#   ./bootstrap.sh git-only        # Only git clones and builds
#   ./bootstrap.sh db-only         # Only database downloads
#   ./bootstrap.sh git-only --exclude-evmutation --exclude-cg-cotrans
#   ./bootstrap.sh db-only --exclude-uniref90
#
# Subcommands:
#   env-only       Only run pip and conda installs (steps 2, 6-6c).
#   git-only       Only run git clones, builds, and conda installs (steps 3-8).
#   db-only        Only run FTP downloads and build_db (steps 9, 12).
#   (none)         Run everything.
#
# Exclude flags (fine-grained control within any mode):
#   --exclude-pip-install     Skip pip install.
#   --exclude-build-tools     Skip apt/yum install of build-essential (gcc, g++, make).
#   --exclude-evmutation      Skip cloning EVmutation and plmc.
#   --exclude-build-plmc      Skip compiling plmc (clone still happens).
#   --exclude-adabmdca        Skip cloning + installing adabmDCApy (GPU Boltzmann backend).
#   --exclude-nextflow        Skip installing Nextflow + OpenJDK.
#   --exclude-cg-cotrans      Skip downloading cg_cotrans.
#   --exclude-netsurfp3       Skip cloning NetSurfP-3.0.
#   --exclude-signalp         Skip the SignalP 6.0 presence check.
#   --exclude-miranda         Skip conda install of miranda.
#   --exclude-spliceai        Skip conda install of spliceai.
#   --exclude-genesplicer     Skip downloading/building GeneSplicer from JHU source.
#   --exclude-clone-af3       Skip cloning AlphaFold3.
#   --exclude-uniref90        Skip UniRef90 FTP download.
#   --exclude-idmapping       Skip UniProt idmapping FTP download.
#   --exclude-build-db        Skip calling build_db.sh.

# --- Defaults: everything on ---
PIP_INSTALL=1
INSTALL_BUILD_TOOLS=1
CLONE_EVMUTATION=1
BUILD_PLMC=1
CLONE_ADABMDCA=1
DOWNLOAD_CG_COTRANS=1
CLONE_NETSURFP3=1
INSTALL_SIGNALP=1
INSTALL_MIRANDA=1
INSTALL_SPLICEAI=1
BUILD_GENESPLICER=1
INSTALL_NEXTFLOW=1
CLONE_AF3=1
DOWNLOAD_UNIREF90=1
DOWNLOAD_IDMAPPING=1
RUN_BUILD_DB=1
INSTALL_EDITABLE_REPOS=1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# apply_plmc_patch <patch_file> <human_label>
#   Applies a patch to plmc with explicit success/already-applied/rejected
#   reporting. Uses --dry-run probes to distinguish "already applied" (reverse
#   patch would succeed) from "broken" (neither forward nor reverse applies).
#   Returns 0 on success or already-applied, 1 on hunk rejection or missing
#   source. The caller decides whether a rejected patch should abort the
#   bootstrap (currently it does, since these patches are critical for q=64
#   codon Potts correctness).
apply_plmc_patch() {
  local patch_file="$1"
  local label="$2"
  local plmc_dir="$BFF_DIR/mutation_effects/plmc"
  if [[ ! -f "$patch_file" ]]; then
    echo "    NOTE plmc patch '$label' not found at $patch_file; skipping"
    return 0
  fi
  if [[ ! -d "$plmc_dir" ]]; then
    echo "    ERROR plmc source dir missing at $plmc_dir; cannot apply '$label'"
    return 1
  fi
  if (cd "$plmc_dir" && patch -p1 --dry-run < "$patch_file" >/dev/null 2>&1); then
    if (cd "$plmc_dir" && patch -p1 < "$patch_file" >/dev/null 2>&1); then
      echo "    APPLIED plmc patch: $label"
      return 0
    fi
    echo "    ERROR plmc patch '$label' dry-run succeeded but apply failed"
    return 1
  fi
  if (cd "$plmc_dir" && patch -p1 -R --dry-run < "$patch_file" >/dev/null 2>&1); then
    echo "    SKIP plmc patch (already applied): $label"
    return 0
  fi
  echo "    WARN plmc patch '$label' did not apply (hunks rejected)"
  echo "         Diagnose with: cd $plmc_dir && patch -p1 < $patch_file"
  return 1
}
# --- Parse subcommand ---
SUBCOMMAND=""
ARGS=()
for arg in "$@"; do
  case "$arg" in
    env-only|git-only|db-only) SUBCOMMAND="$arg" ;;
    *) ARGS+=("$arg") ;;
  esac
done

# Apply subcommand: disable groups not selected
case "$SUBCOMMAND" in
  env-only)
    INSTALL_BUILD_TOOLS=0
    CLONE_EVMUTATION=0; BUILD_PLMC=0; CLONE_ADABMDCA=0; DOWNLOAD_CG_COTRANS=0
    CLONE_NETSURFP3=0; INSTALL_SIGNALP=0; CLONE_AF3=0
    DOWNLOAD_UNIREF90=0; DOWNLOAD_IDMAPPING=0; RUN_BUILD_DB=0
    # INSTALL_NEXTFLOW stays on (env-only still wants nextflow + OpenJDK present)
    ;;
  git-only)
    PIP_INSTALL=0
    DOWNLOAD_UNIREF90=0; DOWNLOAD_IDMAPPING=0; RUN_BUILD_DB=0
    ;;
  db-only)
    PIP_INSTALL=0
    INSTALL_BUILD_TOOLS=0
    INSTALL_EDITABLE_REPOS=0
    CLONE_EVMUTATION=0; BUILD_PLMC=0; CLONE_ADABMDCA=0; DOWNLOAD_CG_COTRANS=0
    CLONE_NETSURFP3=0; INSTALL_SIGNALP=0; INSTALL_MIRANDA=0
    BUILD_GENESPLICER=0; INSTALL_NEXTFLOW=0; CLONE_AF3=0
    ;;
esac

# --- Parse exclude flags ---
set -- "${ARGS[@]+"${ARGS[@]}"}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exclude-pip-install)   PIP_INSTALL=0 ;;
    --exclude-build-tools)   INSTALL_BUILD_TOOLS=0 ;;
    --exclude-evmutation)    CLONE_EVMUTATION=0 ;;
    --exclude-build-plmc)    BUILD_PLMC=0 ;;
    --exclude-adabmdca)      CLONE_ADABMDCA=0 ;;
    --exclude-nextflow)      INSTALL_NEXTFLOW=0 ;;
    --exclude-cg-cotrans)    DOWNLOAD_CG_COTRANS=0 ;;
    --exclude-netsurfp3)     CLONE_NETSURFP3=0 ;;
    --exclude-signalp)       INSTALL_SIGNALP=0 ;;
    --exclude-miranda)       INSTALL_MIRANDA=0 ;;
    --exclude-spliceai)      INSTALL_SPLICEAI=0 ;;
    --exclude-genesplicer)   BUILD_GENESPLICER=0 ;;
    --exclude-clone-af3)     CLONE_AF3=0 ;;
    --exclude-uniref90)      DOWNLOAD_UNIREF90=0 ;;
    --exclude-idmapping)     DOWNLOAD_IDMAPPING=0 ;;
    --exclude-build-db)      RUN_BUILD_DB=0 ;;
    --pip-install|--build-plmc|--download-uniref90|--download-idmapping|--clone-alphafold3)
      ;; # legacy no-ops
    -h|--help)
      sed -n '/^# Usage:/,/^$/p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
  shift
done

# --- Validate contradictions ---
if [[ "$SUBCOMMAND" == "env-only" && "$PIP_INSTALL" -eq 0 ]]; then
  echo "ERROR: env-only with --exclude-pip-install is contradictory." >&2
  exit 1
fi
if [[ "$SUBCOMMAND" == "db-only" && "$DOWNLOAD_UNIREF90" -eq 0 && "$DOWNLOAD_IDMAPPING" -eq 0 && "$RUN_BUILD_DB" -eq 0 ]]; then
  echo "ERROR: db-only with all database steps excluded leaves nothing to do." >&2
  exit 1
fi
if [[ "$SUBCOMMAND" == "git-only" && "$INSTALL_BUILD_TOOLS" -eq 0 && "$CLONE_EVMUTATION" -eq 0 && "$BUILD_PLMC" -eq 0 && "$CLONE_ADABMDCA" -eq 0 && "$DOWNLOAD_CG_COTRANS" -eq 0 && "$CLONE_NETSURFP3" -eq 0 && "$INSTALL_SIGNALP" -eq 0 && "$INSTALL_MIRANDA" -eq 0 && "$BUILD_GENESPLICER" -eq 0 && "$INSTALL_NEXTFLOW" -eq 0 && "$CLONE_AF3" -eq 0 ]]; then
  echo "ERROR: git-only with all git/build steps excluded leaves nothing to do." >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$ROOT_DIR/.." && pwd)"
BFF_DIR="$REPO_ROOT/biofeaturefactory"
cd "$ROOT_DIR"

mkdir -p _downloads

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: $1" >&2
    exit 1
  }
}

# version_at_least <current> <required>
#   Returns 0 if 'current' is >= 'required' under natural version ordering,
#   1 otherwise. Driver: `sort -V`. Tolerates trailing build suffixes
#   (e.g. "11.0.21+9-LTS") by stripping at the first non-[0-9.] character.
version_at_least() {
  local current="$1" required="$2"
  current="${current%%[^0-9.]*}"
  required="${required%%[^0-9.]*}"
  [[ -z "$current" ]] && return 1
  [[ "$(printf '%s\n%s\n' "$required" "$current" | sort -V | head -n1)" == "$required" ]]
}

# java_major <version-string> -> echoes integer major version
#   Java prints either "1.8.0_321" (Java 8) or "11.0.21" (Java 9+). Handle both.
java_major() {
  local v="$1"
  if [[ "$v" == 1.* ]]; then
    echo "$v" | awk -F. '{print $2}'
  else
    echo "$v" | awk -F. '{print $1}'
  fi
}

download_file() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "  EXISTS $out"
    return
  fi
  if command -v aria2c >/dev/null 2>&1; then
    aria2c --continue=true \
      --max-connection-per-server=4 \
      --split=4 \
      --min-split-size=4M \
      --max-tries=20 \
      --retry-wait=10 \
      --timeout=60 \
      --lowest-speed-limit=50K \
      --file-allocation=none \
      -o "$(basename "$out")" \
      -d "$(dirname "$out")" \
      "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --retry 5 --retry-delay 3 -o "$out" "$url"
  else
    echo "ERROR: need aria2c or curl to download $url" >&2
    exit 1
  fi
}

clone_or_update() {
  local repo_url="$1"
  local dest="$2"
  if [[ -d "$dest/.git" ]]; then
    echo "  UPDATE $dest"
    git -C "$dest" fetch --all --tags
    git -C "$dest" pull --ff-only || true
  elif [[ -d "$dest" ]]; then
    echo "  SKIP $dest exists but is not a git repo"
  else
    echo "  CLONE $repo_url -> $dest"
    git clone "$repo_url" "$dest"
  fi
}

# ── Step 1: Core tool checks ─────────────────────────────────────────────
echo "[1/12] Validating core tools..."
require_cmd git
require_cmd tar

# pkg_install <package_list...>
#   Linux system package installer (apt-get / yum / dnf). Routes through
#   sudo when not running as root. BFF is Linux-only end-to-end (the GPU
#   pipelines and most cluster tooling assume Linux); macOS branches were
#   removed deliberately. Prints WARN if no supported package manager is
#   found.
pkg_install() {
  local sudo_cmd=""
  if [[ "$EUID" -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
    sudo_cmd="sudo"
  fi
  if command -v apt-get >/dev/null 2>&1; then
    echo "  APT install $*"
    $sudo_cmd apt-get update -qq
    $sudo_cmd apt-get install -y "$@"
    return $?
  elif command -v yum >/dev/null 2>&1; then
    echo "  YUM install $*"
    $sudo_cmd yum install -y "$@"
    return $?
  elif command -v dnf >/dev/null 2>&1; then
    echo "  DNF install $*"
    $sudo_cmd dnf install -y "$@"
    return $?
  fi
  echo "  WARN no apt-get/yum/dnf found; install manually: $*" >&2
  return 1
}

# ── Step 1b: Build toolchain (gcc, g++, make) ───────────────────────────
# Required for compiling plmc and (potentially) GeneSplicer + native Python
# extensions during pip install. Skipped in env-only / db-only.
# Policy: verify presence; do NOT change/upgrade if already installed.
echo "[1b/12] Build toolchain..."
if [[ "$INSTALL_BUILD_TOOLS" -eq 1 ]]; then
  have_gcc=0; have_gxx=0; have_make=0
  command -v gcc  >/dev/null 2>&1 && have_gcc=1
  command -v g++  >/dev/null 2>&1 && have_gxx=1
  command -v make >/dev/null 2>&1 && have_make=1
  if [[ "$have_gcc" -eq 1 && "$have_gxx" -eq 1 && "$have_make" -eq 1 ]]; then
    gcc_ver="$(gcc -dumpversion 2>/dev/null || echo unknown)"
    make_ver="$(make --version 2>/dev/null | head -n1)"
    echo "  OK gcc $gcc_ver, g++ $(g++ -dumpversion 2>/dev/null || echo unknown), $make_ver — no install needed"
  else
    missing=()
    [[ "$have_gcc"  -eq 0 ]] && missing+=("gcc")
    [[ "$have_gxx"  -eq 0 ]] && missing+=("g++")
    [[ "$have_make" -eq 0 ]] && missing+=("make")
    echo "  MISSING: ${missing[*]}"
    if command -v apt-get >/dev/null 2>&1; then
      pkg_install build-essential
    else
      pkg_install gcc gcc-c++ make
    fi
  fi
fi

# ── Step 2: Pip install ──────────────────────────────────────────────────
echo "[2/12] Core Python requirements..."
if [[ "$PIP_INSTALL" -eq 1 ]]; then
  echo "  PIP install -e .[all] (editable install with all optional deps)"
  python -m pip install -e "${REPO_ROOT}[all]"
  echo "  PIP reinstall pyarrow (numpy ABI fix)"
  python -m pip install pyarrow --force-reinstall --quiet
fi

# ── Step 3: mutation_effects (EVmutation + plmc) + adabmDCApy + cg_cotrans ──
echo "[3/12] mutation_effects module dependencies..."
if [[ "$CLONE_EVMUTATION" -eq 1 ]]; then
  clone_or_update "https://github.com/debbiemarkslab/EVmutation.git" "$BFF_DIR/mutation_effects/EVmutation"
  clone_or_update "https://github.com/debbiemarkslab/plmc.git" "$BFF_DIR/mutation_effects/plmc"
  echo "  PATCH plmc"
  apply_plmc_patch "$SCRIPT_DIR/patches/plmc_proximal_group_lasso.patch" \
                   "proximal-group-lasso"
  apply_plmc_patch "$SCRIPT_DIR/patches/plmc_case_sensitive_custom_alphabet.patch" \
                   "case-sensitive-custom-alphabet (q=64 codon Potts)"
  if [[ "$BUILD_PLMC" -eq 1 ]]; then
    echo "  BUILD plmc"
    (
      cd "$BFF_DIR/mutation_effects/plmc"
      make all-openmp || make all
    )
  fi
fi

# adabmDCApy: PyTorch implementation of adabmDCA. adabmdca_pipeline.py imports
# it IN-PROCESS (`from adabmDCA.training import ...`) to train/score — the
# `adabmDCA` console script is never invoked, so the package must be IMPORTABLE
# (CLI on PATH is irrelevant). Used by Nextflow processes run_protein_adabmdca /
# run_codon_adabmdca. Repo: spqb/adabmDCApy (NOT spqb/adabmDCA — a different package).
ADABMDCA_MIN="0.7.5"
if [[ "$CLONE_ADABMDCA" -eq 1 ]]; then
  ADABM_DIR="$BFF_DIR/mutation_effects/adabmDCApy"

  # Fresh clone every run. Remove any prior state (partial clone, stale checkout)
  # so `git clone` has a clean target.
  rm -rf "$ADABM_DIR"

  (
    cd "$BFF_DIR/mutation_effects"
    git clone https://github.com/spqb/adabmDCApy.git
  )
  # Editable install (pip install -e) + adabmDCA version verification are handled
  # centrally in step [8b], so that env-only (which does not clone here) and
  # git-only (which does not run the core pip install) both install it.
fi

if [[ "$DOWNLOAD_CG_COTRANS" -eq 1 ]]; then
  CG_DIR="$BFF_DIR/rare_codon/cg_cotrans"
  if [[ -d "$CG_DIR" ]]; then
    echo "  EXISTS $CG_DIR"
  else
    echo "  SKIP cg_cotrans (manual download required)"
  fi
fi

# ── Step 4: NetSurfP3 ───────────────────────────────────────────────────
echo "[4/12] NetSurfP3 module dependency..."
if [[ "$CLONE_NETSURFP3" -eq 1 ]]; then
  clone_or_update "https://github.com/Eryk96/NetSurfP-3.0.git" "$BFF_DIR/NetSurfP3/nsp3"
  # Editable install (pip install -e nsp3/nsp3) is handled centrally in step [8b].
fi

# ── Step 5: SignalP 6.0 ─────────────────────────────────────────────────
# CHECK ONLY. SignalP 6.0 is a licensed manual install like the other DTU tools
# (netNglyc, netphos, netMHC) -- see the step 11 checklist. It is deliberately NOT
# pip-installed here, for three measured reasons:
#
#   1. github.com/fteufel/signalp-6.0 is the TRAINING repo, not the prediction
#      distribution. Its package is src/signalp6/ (training_utils, no predict.py),
#      its setup.py declares NO entry_points and NO package_data, so no install of
#      it can ever create the `signalp6` executable that SignalP6Handler shells out
#      to, and it ships no model weights. The licensed DTU tarball is a DIFFERENT
#      package (module `signalp`, entry_points signalp6=signalp:predict,
#      package_data model_weights/*) and is self-sufficient.
#   2. That repo's requirements.txt pins scikit-learn==0.23.1 / scipy==1.5.0 /
#      sklearn==0.0 / torch==1.7.0. `pip install -r` on it FAILS
#      ("Failed to build 'scikit-learn'"), and because this script runs under
#      `set -euo pipefail` that aborted the entire bootstrap here at step 5 of 12 --
#      before miranda, mmseqs2, hmmer, spliceai, genesplicer, nextflow, the AF3
#      clone, the editable installs and build_db.sh.
#   3. SignalP 6.0 requires torch<2 (installation_instructions.md: "SignalP 6.0 is
#      not compatible with PyTorch 2.0+"). nsp3 and AlphaFold3 require modern
#      torch. The pin is mutually exclusive with the rest of this env, so SignalP
#      MUST live in its own environment. netnglyc_pipeline.py only ever invokes it
#      as a subprocess and reads prediction_results.txt, so a separate env is a
#      supported configuration, not a workaround.
echo "[5/12] SignalP 6.0 (netNglyc dependency, licensed manual install)..."
if [[ "$INSTALL_SIGNALP" -eq 1 ]]; then
  if command -v signalp6 >/dev/null 2>&1; then
    echo "  OK signalp6 on PATH: $(command -v signalp6)"
  else
    echo "  MISSING signalp6 not on PATH — see the step 11 checklist for install steps."
    echo "          netNglyc runs will FAIL at SignalP6Handler until it is present."
  fi
fi

# ── Step 6: Miranda ─────────────────────────────────────────────────────
echo "[6/12] Miranda (conda)..."
if [[ "$INSTALL_MIRANDA" -eq 1 ]]; then
  if command -v miranda >/dev/null 2>&1; then
    echo "  OK miranda already on PATH"
  elif command -v conda >/dev/null 2>&1; then
    echo "  CONDA install miranda"
    conda install -y -c bioconda miranda
  elif command -v mamba >/dev/null 2>&1; then
    echo "  MAMBA install miranda (conda not found)"
    mamba install -y -c bioconda miranda
  else
    echo "  WARN conda/mamba not found; install miranda manually: conda install -c bioconda miranda"
  fi
fi

# ── Step 6b: mmseqs2 (EVmutation codon MSA) ─────────────────────────────
echo "[6b/12] mmseqs2..."
if command -v mmseqs >/dev/null 2>&1; then
  echo "  OK mmseqs2 already on PATH"
elif command -v conda >/dev/null 2>&1; then
  echo "  CONDA install mmseqs2"
  conda install -y -c bioconda mmseqs2
elif command -v mamba >/dev/null 2>&1; then
  echo "  MAMBA install mmseqs2 (conda not found)"
  mamba install -y -c bioconda mmseqs2
else
  echo "  WARN conda/mamba not found; install mmseqs2 manually: conda install -c bioconda mmseqs2"
fi

# ── Step 6c: HMMER / jackhmmer (EVmutation protein MSA) ────────────────
echo "[6c/12] HMMER (jackhmmer)..."
if command -v jackhmmer >/dev/null 2>&1; then
  echo "  OK jackhmmer already on PATH"
elif command -v conda >/dev/null 2>&1; then
  echo "  CONDA install hmmer"
  conda install -y -c bioconda hmmer
elif command -v mamba >/dev/null 2>&1; then
  echo "  MAMBA install hmmer (conda not found)"
  mamba install -y -c bioconda hmmer
else
  echo "  WARN conda/mamba not found; install hmmer manually: conda install -c bioconda hmmer"
fi

# ── Step 6d: SpliceAI ───────────────────────────────────────────────────
# Placed here, with the other conda installs, rather than in the steps 10-12
# block where the spliceai CHECK lives: that block is skipped under `env-only`,
# and env-only is documented as "only run pip and conda installs". A conda
# package belongs in the phase that installs conda packages.
echo "[6d/12] SpliceAI (conda)..."
if [[ "$INSTALL_SPLICEAI" -eq 1 ]]; then
  if command -v spliceai >/dev/null 2>&1; then
    echo "  OK spliceai already on PATH"
  elif command -v conda >/dev/null 2>&1; then
    echo "  CONDA install spliceai"
    conda install -y -c bioconda spliceai
  else
    echo "  WARN conda not found; install spliceai manually: conda install -c bioconda spliceai"
  fi
fi

# ── Step 7: GeneSplicer ─────────────────────────────────────────────────
# SOURCE BUILD, not conda. `conda install -c bioconda genesplicer` installs a
# DIFFERENT binary (md5 b2e2384f..., 70848 bytes) from the one JHU ships
# (md5 61a8648a..., 53448 bytes). Measured on this stack: identical input and
# identical human model, the conda build writes NOTHING to stdout and exits 0,
# while the source build reports splice sites. Both are valid executables with
# the same six flags (-f -a -d -e -i -h), so nothing about the invocation can
# recover it -- and genesplicer_ensemble.py cannot tell an empty stdout with
# rc=0 from a sequence that genuinely has no splice sites, so every gene would
# be silently reported as site-free.
echo "[7/12] GeneSplicer binary (source build from JHU)..."
if [[ "$BUILD_GENESPLICER" -eq 1 ]]; then
  GS_DIR="$BFF_DIR/genesplicer"
  GS_TAR="$GS_DIR/GeneSplicer.tar.gz"
  GS_SRC="$GS_DIR/GeneSplicer"
  GS_BIN="$GS_SRC/sources/genesplicer"

  if [[ -x "$GS_BIN" ]]; then
    echo "  EXISTS $GS_BIN"
  else
    mkdir -p "$GS_DIR"
    download_file "ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz" "$GS_TAR"
    tar -xzf "$GS_TAR" -C "$GS_DIR"
    rm -f "$GS_TAR"
    if [[ -d "$GS_SRC/sources" ]]; then
      echo "  BUILD GeneSplicer from source"
      (cd "$GS_SRC/sources" && make)
    else
      echo "  ERROR GeneSplicer sources/ not found after extracting $GS_TAR" >&2
    fi
  fi

  if [[ -x "$GS_BIN" ]]; then
    echo "  OK genesplicer built: $GS_BIN"
    echo "     human model:       $GS_SRC/human"
    echo "     pass to the pipeline with:"
    echo "       --genesplicer-dir $GS_SRC/sources --model-dir $GS_SRC/human"
    if command -v genesplicer >/dev/null 2>&1; then
      echo "  NOTE a 'genesplicer' is also on PATH ($(command -v genesplicer))."
      echo "       If it is the bioconda build it emits no sites; prefer the path above."
    fi
  else
    echo "  WARN genesplicer was not built; the GeneSplicer pipeline cannot run." >&2
  fi
fi

# ── Step 7b: Nextflow + OpenJDK ──────────────────────────────────────────
# Floors:
#   Java 17  — Nextflow 24+ requires it.
#   Nextflow 22.10.0 — minimum where the DSL2 operators used here are stable.
# Policy: verify version; if at/above the floor, skip install (no version
# change). If below, WARN with the manual upgrade command rather than
# silently replacing — system Java affects unrelated tooling.
JAVA_MIN_MAJOR=17
NEXTFLOW_MIN="22.10.0"
echo "[7b/12] Nextflow + OpenJDK..."
if [[ "$INSTALL_NEXTFLOW" -eq 1 ]]; then
  # ---- OpenJDK ----
  # Nextflow 24+ requires Java 17+. Earlier policy was warn-and-continue,
  # but the official Nextflow installer fails immediately when Java is below
  # the floor — leaving the user with a broken nextflow shim. Now: when Java
  # is too old, install OpenJDK 17 alongside (apt/yum don't replace the older
  # JVM) and re-point `update-alternatives` so `/usr/bin/java` resolves to
  # the new one. The older Java remains on disk for anything that needs it.
  needs_java_install=0
  if command -v java >/dev/null 2>&1; then
    java_raw="$(java -version 2>&1 | head -n 1 | awk -F\" '{print $2}')"
    jmaj="$(java_major "$java_raw")"
    if [[ -n "$jmaj" && "$jmaj" -ge "$JAVA_MIN_MAJOR" ]]; then
      echo "  OK java $java_raw (>= $JAVA_MIN_MAJOR) — no install needed"
    else
      echo "  java $java_raw is below the required floor (Java $JAVA_MIN_MAJOR+ for Nextflow 24+)."
      echo "  Installing OpenJDK $JAVA_MIN_MAJOR alongside; switching default java -> 17."
      needs_java_install=1
    fi
  else
    echo "  java not found — installing OpenJDK $JAVA_MIN_MAJOR"
    needs_java_install=1
  fi

  if [[ "$needs_java_install" -eq 1 ]]; then
    if command -v apt-get >/dev/null 2>&1; then
      pkg_install openjdk-17-jdk
    else
      pkg_install java-17-openjdk-devel
    fi

    # Locate the Java 17 binary and switch update-alternatives so subsequent
    # `java` invocations (including `nextflow`'s installer) pick it up.
    arch="$(dpkg --print-architecture 2>/dev/null || uname -m)"
    JAVA_17_BIN=""
    for candidate in \
        "/usr/lib/jvm/java-17-openjdk-${arch}/bin/java" \
        "/usr/lib/jvm/java-17-openjdk-amd64/bin/java" \
        "/usr/lib/jvm/java-17-openjdk-arm64/bin/java" \
        "/usr/lib/jvm/java-1.17.0-openjdk/bin/java" \
        "/usr/lib/jvm/jre-17-openjdk/bin/java"
    do
      [[ -x "$candidate" ]] && JAVA_17_BIN="$candidate" && break
    done
    if [[ -z "$JAVA_17_BIN" ]]; then
      JAVA_17_BIN="$(find /usr/lib/jvm -maxdepth 3 -name java -path '*java-17*' -type f 2>/dev/null | head -1)"
    fi

    if [[ -x "$JAVA_17_BIN" ]]; then
      sudo_cmd=""
      [[ "$EUID" -ne 0 ]] && command -v sudo >/dev/null 2>&1 && sudo_cmd="sudo"
      $sudo_cmd update-alternatives --set java "$JAVA_17_BIN" 2>/dev/null || true
      export JAVA_HOME="$(dirname "$(dirname "$JAVA_17_BIN")")"
      export JAVA_CMD="$JAVA_17_BIN"
      echo "  Java 17 active: $JAVA_17_BIN"
      echo "  JAVA_HOME=$JAVA_HOME (exported for this bootstrap session + nextflow installer)"
      echo "  NOTE: this updates the system default. To persist for your shell, add to ~/.bashrc:"
      echo "         export JAVA_HOME=$JAVA_HOME"
    else
      echo "  ERROR: installed openjdk-17 but could not locate the Java 17 binary." >&2
      echo "         Looked under /usr/lib/jvm/java-17-openjdk-*. Nextflow install will likely fail." >&2
      exit 1
    fi
  fi

  # ---- Nextflow ----
  if command -v nextflow >/dev/null 2>&1; then
    # Nextflow prints "      version 25.10.4 build 11173" — $1=version, $2=25.10.4.
    nf_ver="$(nextflow -version 2>&1 | awk '$1=="version" {print $2; exit}')"
    if [[ -n "$nf_ver" ]] && version_at_least "$nf_ver" "$NEXTFLOW_MIN"; then
      echo "  OK nextflow $nf_ver (>= $NEXTFLOW_MIN) — no install needed"
    else
      echo "  WARN nextflow ${nf_ver:-unknown} is below the required floor ($NEXTFLOW_MIN)."
      echo "       Not auto-upgrading. Upgrade manually:"
      echo "         self-update:  nextflow self-update"
      echo "         conda:        conda install -y -c bioconda 'nextflow>=$NEXTFLOW_MIN'"
    fi
  else
    echo "  nextflow not found — installing"
    if command -v conda >/dev/null 2>&1; then
      echo "  CONDA install nextflow"
      conda install -y -c bioconda "nextflow>=${NEXTFLOW_MIN}"
    elif command -v mamba >/dev/null 2>&1; then
      echo "  MAMBA install nextflow (conda not found)"
      mamba install -y -c bioconda "nextflow>=${NEXTFLOW_MIN}"
    else
      # Fallback: official installer at https://get.nextflow.io writes the
      # launcher to $PWD/nextflow and chmods +rx. It does NOT place it on
      # PATH; that's our job. Steps:
      #   1) Run installer in a clean tmp dir so we don't pollute scripts/.
      #   2) Try destinations in order: ~/.local/bin (if on PATH and
      #      writable), then /usr/local/bin (sudo-or-root), then leave it
      #      in tmp and surface the absolute path.
      #   3) Verify `nextflow` resolves on PATH after install.
      echo "  Installing nextflow via official installer (https://get.nextflow.io)"
      nf_tmp="$(mktemp -d)"
      (
        cd "$nf_tmp"
        # -fsSL: fail on HTTP errors, silent body, show errors, follow redirects.
        # set -o pipefail (already enabled) propagates curl failure through the pipe.
        curl -fsSL https://get.nextflow.io | bash
      )
      nf_bin="$nf_tmp/nextflow"
      if [[ ! -x "$nf_bin" ]]; then
        echo "  ERROR: installer ran but $nf_bin is missing or not executable" >&2
        rm -rf "$nf_tmp"
        exit 1
      fi

      # Pick a destination.
      install_dest=""
      # Candidate 1: ~/.local/bin (user-owned, no sudo). Only if on PATH.
      local_bin="$HOME/.local/bin"
      case ":$PATH:" in
        *":$local_bin:"*)
          mkdir -p "$local_bin"
          if [[ -w "$local_bin" ]]; then
            install_dest="$local_bin/nextflow"
          fi
          ;;
      esac
      # Candidate 2: /usr/local/bin, sudo if non-root.
      if [[ -z "$install_dest" ]]; then
        if [[ -w "/usr/local/bin" ]]; then
          install_dest="/usr/local/bin/nextflow"
          mv "$nf_bin" "$install_dest"
        elif [[ "$EUID" -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
          echo "  Requesting sudo to install to /usr/local/bin"
          if sudo mv "$nf_bin" "/usr/local/bin/nextflow"; then
            install_dest="/usr/local/bin/nextflow"
          fi
        fi
      else
        mv "$nf_bin" "$install_dest"
      fi

      # Candidate 3 (fallback): leave in a stable location and tell the user.
      if [[ -z "$install_dest" ]]; then
        install_dest="$ROOT_DIR/nextflow"
        mv "$nf_bin" "$install_dest"
        echo "  WARN nextflow installed at $install_dest but no PATH-writable"
        echo "       destination was found. Add it to PATH manually, e.g.:"
        echo "         mkdir -p \$HOME/.local/bin && mv \"$install_dest\" \$HOME/.local/bin/"
        echo "         # then ensure \$HOME/.local/bin is on \$PATH"
      else
        echo "  Installed nextflow -> $install_dest"
      fi
      rm -rf "$nf_tmp"

      # Verify
      if command -v nextflow >/dev/null 2>&1; then
        echo "  OK nextflow on PATH: $(command -v nextflow)"
      else
        echo "  WARN nextflow installed at $install_dest but not yet on PATH"
        echo "       (open a new shell, or 'export PATH=\"$(dirname "$install_dest"):\$PATH\"')"
      fi
    fi
  fi
fi

# ── Step 8: AlphaFold3 ──────────────────────────────────────────────────
echo "[8/12] AlphaFold3 upstream (optional clone)..."
if [[ "$CLONE_AF3" -eq 1 ]]; then
  clone_or_update "https://github.com/google-deepmind/alphafold3.git" "$BFF_DIR/alphafold3/alphafold3"
fi

# ── Step 8b: Editable installs of python-based cloned repos ─────────────
# Runs in full / env-only / git-only (not db-only). Decoupled from the core
# `pip install -e .[all]` (PIP_INSTALL) and from the per-repo CLONE_* flags so:
#   - env-only installs whatever python repos are already cloned, and warns for
#     any repo that is absent (nothing was cloned this run);
#   - git-only installs the repos it just cloned, and warns if the core BFF env
#     was never set up.
# EVmutation is intentionally NOT here: it ships no setup.py/pyproject and is
# imported via sys.path from its clone, so it cannot be `pip install -e`'d.
echo "[8b/12] Editable installs of python-based cloned repos..."
if [[ "$INSTALL_EDITABLE_REPOS" -eq 1 ]]; then
  # Each entry: label|relpath-under-BFF_DIR|install-subdir-suffix|setup-marker
  EDITABLE_REPOS=(
    "nsp3|NetSurfP3/nsp3|/nsp3|setup.py"
    "adabmDCApy|mutation_effects/adabmDCApy||pyproject.toml"
  )

  # git-only (or any mode that skipped the core pip install): the BFF package and
  # its dependencies may be missing. Warn but still attempt the editable installs.
  if [[ "$PIP_INSTALL" -eq 0 ]] && ! python -c "import biofeaturefactory" >/dev/null 2>&1; then
    echo "  WARN BioFeatureFactory core env not detected ('import biofeaturefactory' failed)."
    echo "       Run './bootstrap.sh env-only' or the full bootstrap to install python dependencies."
  fi

  for entry in "${EDITABLE_REPOS[@]}"; do
    IFS='|' read -r er_label er_rel er_suffix er_marker <<< "$entry"
    er_repo_dir="$BFF_DIR/$er_rel"
    er_install_dir="$er_repo_dir$er_suffix"
    if [[ -d "$er_repo_dir" && -f "$er_install_dir/$er_marker" ]]; then
      echo "  PIP install -e $er_label ($er_install_dir)"
      python -m pip install -e "$er_install_dir"
    else
      echo "  WARN $er_label repo not present at $er_repo_dir"
      echo "       Run the full bootstrap, or './bootstrap.sh git-only', to install the proper repos."
    fi
  done

  # adabmDCA version + CLI verification (warn-only), after the editable install above.
  if [[ "$CLONE_ADABMDCA" -eq 1 || -d "$BFF_DIR/mutation_effects/adabmDCApy" ]]; then
    installed_ver="$(python -m pip show adabmDCA 2>/dev/null | awk -F': ' '/^Version:/ {print $2}')"
    if [[ -n "$installed_ver" ]]; then
      if version_at_least "$installed_ver" "$ADABMDCA_MIN"; then
        echo "  OK adabmDCA $installed_ver (>= $ADABMDCA_MIN)"
      else
        echo "  WARN installed adabmDCA $installed_ver is below floor ($ADABMDCA_MIN)"
      fi
    else
      echo "  WARN adabmDCA package not detected via pip show — install may have failed"
    fi
    if command -v adabmDCA >/dev/null 2>&1; then
      echo "  OK adabmDCA CLI on PATH: $(command -v adabmDCA)"
    else
      echo "  WARN adabmDCA CLI not on PATH after install (env activation issue?)"
    fi
  fi
fi

# ── Step 9: FTP datasets ────────────────────────────────────────────────
echo "[9/12] Optional FTP datasets..."
if [[ "$DOWNLOAD_UNIREF90" -eq 1 ]]; then
  download_file \
    "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" \
    "$ROOT_DIR/_downloads/uniref90.fasta.gz"
fi
if [[ "$DOWNLOAD_IDMAPPING" -eq 1 ]]; then
  download_file \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz" \
    "$ROOT_DIR/_downloads/idmapping.dat.gz"
fi

# ── Steps 10-12 only run in full or git-only mode ────────────────────────
if [[ -z "$SUBCOMMAND" || "$SUBCOMMAND" == "git-only" ]]; then
  echo "[10/12] SpliceAI/Nextflow checks..."
  if command -v spliceai >/dev/null 2>&1; then
    echo "  OK spliceai on PATH"
  else
    # Reachable when --exclude-spliceai was passed, when neither conda nor mamba
    # is available, or when step 6d's install failed.
    echo "  WARN spliceai not found (installed at step 6d; see --exclude-spliceai)"
  fi
  if command -v nextflow >/dev/null 2>&1; then
    echo "  OK nextflow on PATH"
  else
    echo "  WARN nextflow not found (required by spliceai/README.md)"
  fi

  echo "[11/12] Licensed/manual dependencies checklist..."
  cat <<'EOF'
Manual installs required (cannot be auto-downloaded):
  - cg_cotrans/: Download from https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
    Extract into biofeaturefactory/rare_codon/cg_cotrans/
  - netNglyc/: NetNGlyc 1.0 (DTU academic license)
    Patch NetNGlyc tcsh SIGNALP path to:
      biofeaturefactory/netNglyc/bin/signalp6_adapter
  - SignalP 6.0 (DTU academic license) -- REQUIRED by netNglyc
    Download the "fast" package: https://services.healthtech.dtu.dk/service.php?SignalP-6.0
    Install it in a SEPARATE python env: SignalP 6.0 pins torch<2, while nsp3 and
    AlphaFold3 need modern torch, so it CANNOT share the BFF env. netNglyc invokes
    it as a subprocess and reads prediction_results.txt, so a separate env is fine.
      conda create -n signalp6 python=3.10 && conda activate signalp6
      pip install signalp-6-package/
      pip install "numpy<2"                       # required for 6.0h
      SIGNALP_DIR=$(python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
      cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
    Then put that env's `signalp6` on PATH for netNglyc runs.
    NOTE: github.com/fteufel/signalp-6.0 is the TRAINING repo -- it declares no
    console_scripts and ships no weights, so it cannot satisfy this dependency.
  - netphos/: NetPhos 3.1 + APE (DTU academic license), requires tcsh
  - netMHC/: NetMHCpan 4.1 / NetMHC 4.0 (DTU academic license)
    NetMHC 4.0 data files: https://services.healthtech.dtu.dk/services/NetMHC-4.0/data.tar.gz
    Extract into the netMHC-4.0 installation directory.
  - alphafold3/: NVIDIA GPU stack + Docker + AF3 model weights
EOF

  echo "[12/12] Summary checks..."
  [[ -d "$BFF_DIR/mutation_effects/EVmutation" ]]   && echo "  OK EVmutation clone"
  [[ -d "$BFF_DIR/mutation_effects/plmc" ]]         && echo "  OK plmc clone"
  [[ -d "$BFF_DIR/mutation_effects/adabmDCApy" ]]   && echo "  OK adabmDCApy clone"
  command -v adabmDCA >/dev/null 2>&1                 && echo "  OK adabmDCA CLI on PATH"
  [[ -d "$BFF_DIR/rare_codon/cg_cotrans" ]]         && echo "  OK cg_cotrans"
  [[ -d "$BFF_DIR/NetSurfP3/nsp3" ]]                && echo "  OK NetSurfP3 clone"
  command -v signalp6 >/dev/null 2>&1                 && echo "  OK signalp6 on PATH"
  command -v miranda >/dev/null 2>&1                  && echo "  OK miranda on PATH"
  [[ -x "$BFF_DIR/genesplicer/GeneSplicer/sources/genesplicer" ]] && echo "  OK genesplicer source build"
  [[ -d "$BFF_DIR/alphafold3/alphafold3" ]]         && echo "  OK AlphaFold3 clone"
fi

if [[ "$RUN_BUILD_DB" -eq 1 ]]; then
  DB_SCRIPT="$ROOT_DIR/build_db.sh"
  LEGACY_DB_SCRIPT="$(cd "$ROOT_DIR/.." && pwd)/Bio_DBs/build_db.sh"
  if [[ -x "$DB_SCRIPT" ]]; then
    echo "  RUN scripts/build_db.sh"
    "$DB_SCRIPT"
  elif [[ -x "$LEGACY_DB_SCRIPT" ]]; then
    echo "  RUN Bio_DBs/build_db.sh (legacy path)"
    "$LEGACY_DB_SCRIPT"
  else
    echo "  WARN build_db.sh not found/executable at $DB_SCRIPT or $LEGACY_DB_SCRIPT"
  fi
fi

echo "Done."