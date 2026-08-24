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
#   ./bootstrap.sh                          # Full install (all phases)
#   ./bootstrap.sh env-phase                # Just the python environment
#   ./bootstrap.sh env-phase git-phase      # Environment + clones/builds
#   ./bootstrap.sh git-phase db-phase       # Clones/builds + datasets
#   ./bootstrap.sh git-phase --exclude-evmutation --exclude-cg-cotrans
#   ./bootstrap.sh db-phase --exclude-uniref90
#
# Phases (ADDITIVE -- name any combination; the run is their UNION):
#   env-phase      pip install, the conda installs, Nextflow/OpenJDK, and editable
#                  installs of repos that are ALREADY cloned (steps 2, 6-6d, 7b,
#                  8b-8c). No clones, no source builds, no dataset downloads.
#   git-phase      Clones, source builds, the conda installs, and editable installs
#                  (steps 1b, 3-8c, 10-12).
#   db-phase       FTP dataset downloads, CoCoPUTs codon table, build_db (steps 9, 9b, 12).
#   (none)         Run every phase.
#
#   Each phase may be written bare: `env`, `git`, `db` are accepted for
#   `env-phase`, `git-phase`, `db-phase`. The selectors are deliberately NOT
#   called "-only": naming two of them runs both, so "only" would be a lie.
#   env-phase and git-phase overlap (both want the conda installs, Nextflow and
#   the editable installs); a step wanted by either phase runs once.
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
#   --exclude-spliceai-env    Install spliceai into the ACTIVE env instead of its own.
#                             Reintroduces the pyarrow/TensorFlow Abseil deadlock; see
#                             step 6d.
#   --spliceai-env NAME       Name of the dedicated spliceai env (default: bff-spliceai).
#   --exclude-mmseqs2         Skip conda install of mmseqs2.
#   --exclude-hmmer           Skip conda install of HMMER (jackhmmer).
#   --exclude-editable-repos  Skip `pip install -e` of cloned python repos (nsp3, adabmDCApy).
#   --exclude-genesplicer     Skip downloading/building GeneSplicer from JHU source.
#   --exclude-clone-af3       Skip cloning AlphaFold3.
#   --exclude-uniref90        Skip UniRef90 FTP download.
#   --exclude-cocoputs        Skip the CoCoPUTs codon-usage table build (rare_codon).
#   --bio-dbs DIR             Prepared-database root (default: nearest Bio_DBs above the
#                             repo, else BFF_BIO_DBS, else scripts/_downloads).
#   --exclude-idmapping       Skip UniProt idmapping FTP download.
#   --exclude-build-db        Skip calling build_db.sh.
#
# Repair flags:
#   --fix-python              Let conda move the env's interpreter to $PY_TARGET when it
#                             is outside the supported range. DESTRUCTIVE: the packages
#                             installed under the current interpreter are orphaned.

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
INSTALL_MMSEQS2=1
INSTALL_HMMER=1
BUILD_GENESPLICER=1
INSTALL_NEXTFLOW=1
CLONE_AF3=1
DOWNLOAD_UNIREF90=1
DOWNLOAD_IDMAPPING=1
BUILD_COCOPUTS_CUT=1
SPLICEAI_OWN_ENV=1
SPLICEAI_ENV_NAME="${SPLICEAI_ENV_NAME:-bff-spliceai}"
RUN_BUILD_DB=1
INSTALL_EDITABLE_REPOS=1
FIX_PYTHON=0
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
# PY_VER / conda_install_pinned
#   MEASURED FAILURE this guards against: `conda install -c bioconda <pkg>` with no
#   python constraint lets the solver satisfy a bioconda dependency by REPLACING the
#   interpreter. On 2026-08-21 that upgraded an env from python 3.11 -> 3.14.7,
#   orphaning all 216 packages under lib/python3.11/site-packages (numpy, torch,
#   pandas, pysam, scipy, biopython) and leaving lib/python3.14/site-packages holding
#   only pip. The subsequent `pip install -e .[all]` then could not resolve
#   numpy<2 to a wheel on cp314 and fell back to compiling numpy from source.
#   Pinning python= makes the solver either honour the interpreter or FAIL LOUDLY.
# resolve_py_bin / PY_BIN / PY_VER
#   Do NOT trust a bare `python`. Run without `conda activate`, modern macOS and
#   most Linux distros have no `python` at all (only `python3`), so a bare probe
#   returns empty and every downstream step silently targets the wrong interpreter
#   -- or none. Prefer the ACTIVE conda env's interpreter, then python3, then
#   python, and use that one binary everywhere in this script.
resolve_py_bin() {
  if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/python" ]]; then
    echo "$CONDA_PREFIX/bin/python"; return 0
  fi
  command -v python3 2>/dev/null && return 0
  command -v python  2>/dev/null && return 0
  return 1
}
py_ver_of() { "$1" -c 'import sys; print("%d.%d" % sys.version_info[:2])' 2>/dev/null || true; }

PY_BIN="$(resolve_py_bin || true)"
PY_VER=""
if [[ -n "$PY_BIN" ]]; then PY_VER="$(py_ver_of "$PY_BIN")"; fi

# conda_env_args
#   Scope every conda write to the ACTIVE env by prefix. Without -p, conda targets
#   whatever it considers active, which for a NON-activated shell is `base` -- so an
#   unscoped `conda install python=3.11` can downgrade the base installation instead
#   of the project env. Empty when no env is active; callers must refuse in that case.
conda_env_args() {
  if [[ -n "${CONDA_PREFIX:-}" ]]; then printf -- '-p\n%s\n' "$CONDA_PREFIX"; fi
}

conda_install_pinned() {
  local label="$1"; shift
  local extra=() pin=()
  while IFS= read -r line; do [[ -n "$line" ]] && extra+=("$line"); done < <(conda_env_args)
  if [[ -n "$PY_VER" ]]; then pin=("python=$PY_VER"); fi
  if command -v conda >/dev/null 2>&1; then
    echo "  CONDA install $label (pinned python=${PY_VER:-unpinned}, env=${CONDA_PREFIX:-<none>})"
    conda install -y "${extra[@]+"${extra[@]}"}" -c bioconda "$@" "${pin[@]+"${pin[@]}"}"
  elif command -v mamba >/dev/null 2>&1; then
    echo "  MAMBA install $label (conda not found; pinned python=${PY_VER:-unpinned})"
    mamba install -y "${extra[@]+"${extra[@]}"}" -c bioconda "$@" "${pin[@]+"${pin[@]}"}"
  else
    echo "  WARN conda/mamba not found; install $label manually: conda install -c bioconda $*"
    return 1
  fi
}

# --- Parse phase selectors ---
# Phases are ADDITIVE. The old env-only/git-only/db-only were mutually exclusive
# and each carried its own hand-maintained disable-list, which is how `db-only`
# ended up still conda-installing spliceai/mmseqs2/hmmer: those three were simply
# missing from its list. Membership is now declared once per phase below and the
# enabled set is the UNION of the phases named, so a step cannot be forgotten.
PHASE_ENV=0; PHASE_GIT=0; PHASE_DB=0
PHASE_SELECTORS=""
ARGS=()
for arg in "$@"; do
  case "$arg" in
    # Handled here, before the phase banner is printed, so `--help` output is clean.
    -h|--help) sed -n '/^# Usage:/,/^$/p' "$0"; exit 0 ;;
    env-phase|env) PHASE_ENV=1; PHASE_SELECTORS="$PHASE_SELECTORS $arg" ;;
    git-phase|git) PHASE_GIT=1; PHASE_SELECTORS="$PHASE_SELECTORS $arg" ;;
    db-phase|db)   PHASE_DB=1;  PHASE_SELECTORS="$PHASE_SELECTORS $arg" ;;
    env-only|git-only|db-only)
      echo "ERROR: '$arg' no longer exists. Phases are additive, so '-only' would be a lie." >&2
      echo "       Use '${arg%%-*}-phase' (or bare '${arg%%-*}'), and combine freely:" >&2
      echo "         ./bootstrap.sh env-phase git-phase" >&2
      exit 1
      ;;
    *) ARGS+=("$arg") ;;
  esac
done

# Flag membership per phase. env-phase and git-phase intentionally overlap on the
# conda installs, Nextflow and the editable installs -- a flag named by EITHER
# phase is enabled, and its step still runs exactly once.
PHASE_ENV_FLAGS="PIP_INSTALL INSTALL_MIRANDA INSTALL_SPLICEAI INSTALL_MMSEQS2 INSTALL_HMMER INSTALL_NEXTFLOW INSTALL_EDITABLE_REPOS"
PHASE_GIT_FLAGS="INSTALL_BUILD_TOOLS CLONE_EVMUTATION BUILD_PLMC CLONE_ADABMDCA DOWNLOAD_CG_COTRANS CLONE_NETSURFP3 INSTALL_SIGNALP INSTALL_MIRANDA INSTALL_SPLICEAI INSTALL_MMSEQS2 INSTALL_HMMER BUILD_GENESPLICER INSTALL_NEXTFLOW CLONE_AF3 INSTALL_EDITABLE_REPOS"
PHASE_DB_FLAGS="DOWNLOAD_UNIREF90 DOWNLOAD_IDMAPPING BUILD_COCOPUTS_CUT RUN_BUILD_DB"
ALL_PHASE_FLAGS="PIP_INSTALL INSTALL_BUILD_TOOLS CLONE_EVMUTATION BUILD_PLMC CLONE_ADABMDCA DOWNLOAD_CG_COTRANS CLONE_NETSURFP3 INSTALL_SIGNALP INSTALL_MIRANDA INSTALL_SPLICEAI INSTALL_MMSEQS2 INSTALL_HMMER BUILD_GENESPLICER INSTALL_NEXTFLOW CLONE_AF3 INSTALL_EDITABLE_REPOS DOWNLOAD_UNIREF90 DOWNLOAD_IDMAPPING BUILD_COCOPUTS_CUT RUN_BUILD_DB"

enable_phase_flags() {
  local f
  for f in $1; do printf -v "$f" 1; done
}

if [[ "$PHASE_ENV" -eq 1 || "$PHASE_GIT" -eq 1 || "$PHASE_DB" -eq 1 ]]; then
  for f in $ALL_PHASE_FLAGS; do printf -v "$f" 0; done
  if [[ "$PHASE_ENV" -eq 1 ]]; then enable_phase_flags "$PHASE_ENV_FLAGS"; fi
  if [[ "$PHASE_GIT" -eq 1 ]]; then enable_phase_flags "$PHASE_GIT_FLAGS"; fi
  if [[ "$PHASE_DB"  -eq 1 ]]; then enable_phase_flags "$PHASE_DB_FLAGS";  fi
  echo "Phases selected:$PHASE_SELECTORS"
else
  PHASE_ENV=1; PHASE_GIT=1; PHASE_DB=1
  echo "Phases selected: env git db (none named -> full bootstrap)"
fi

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
    --exclude-spliceai-env)  SPLICEAI_OWN_ENV=0 ;;
    --spliceai-env)          SPLICEAI_ENV_NAME="${2:?--spliceai-env needs a name}"; shift ;;
    --exclude-mmseqs2)       INSTALL_MMSEQS2=0 ;;
    --exclude-hmmer)         INSTALL_HMMER=0 ;;
    --exclude-editable-repos) INSTALL_EDITABLE_REPOS=0 ;;
    --exclude-genesplicer)   BUILD_GENESPLICER=0 ;;
    --exclude-clone-af3)     CLONE_AF3=0 ;;
    --exclude-uniref90)      DOWNLOAD_UNIREF90=0 ;;
    --exclude-cocoputs)      BUILD_COCOPUTS_CUT=0 ;;
    --bio-dbs)               BIO_DBS_DIR="${2:?--bio-dbs needs a directory}"; shift ;;
    --exclude-idmapping)     DOWNLOAD_IDMAPPING=0 ;;
    --exclude-build-db)      RUN_BUILD_DB=0 ;;
    --fix-python)            FIX_PYTHON=1 ;;
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

# --- Validate: something must remain to do ---
# One generic check replaces the three per-subcommand ones. Those had to be
# hand-updated whenever a flag was added, and the git-only variant had already
# drifted out of sync with the flag list it was meant to cover.
_any_enabled=0
for f in $ALL_PHASE_FLAGS; do
  if [[ "${!f}" -eq 1 ]]; then _any_enabled=1; fi
done
if [[ "$_any_enabled" -eq 0 ]]; then
  echo "ERROR: every step in the selected phase(s) is excluded; nothing to do." >&2
  exit 1
fi
if [[ "$PHASE_ENV" -eq 1 && "$PHASE_GIT" -eq 0 && "$PHASE_DB" -eq 0 && "$PIP_INSTALL" -eq 0 ]]; then
  echo "WARN: env-phase with --exclude-pip-install; only the conda and editable-install steps will run." >&2
fi
# ── Preflight: interpreter must be in the supported range ───────────────
# Ceiling 3.12: pyproject pins numpy>=1.20,<2, and numpy 1.26.4 (the last 1.x) ships
#   no wheel past cp312. Above it pip falls back to an sdist and compiles numpy.
# Floor 3.10: adabmDCA requires it (see pyproject [project.optional-dependencies]).
# 3.11 is the version this stack has actually been run on.
# This gate is a REFUSAL, not a repair. Replacing the interpreter of a populated env
# orphans everything installed under the old one, so the downgrade is opt-in via
# --fix-python rather than automatic.
PY_MIN="3.10"; PY_MAX="3.12"; PY_TARGET="3.11"
py_in_range() {
  local v="$1"
  [[ -n "$v" ]] || return 1
  [[ "$(printf '%s\n%s\n' "$PY_MIN" "$v" | sort -V | head -n1)" == "$PY_MIN" ]] || return 1
  [[ "$(printf '%s\n%s\n' "$v" "$PY_MAX" | sort -V | head -n1)" == "$v" ]] || return 1
  return 0
}

if [[ "$PHASE_ENV" -eq 1 || "$PHASE_GIT" -eq 1 ]]; then
  if [[ -z "$PY_BIN" ]]; then
    echo "ERROR: no python interpreter found (looked at \$CONDA_PREFIX/bin/python, python3, python)." >&2
    echo "       Activate the project environment first:  conda activate bff" >&2
    exit 1
  fi

  # --fix-python is attempted BEFORE any failure is reported: a repair that succeeds
  # is not an error, and printing the refusal first made a successful run look broken.
  if ! py_in_range "$PY_VER" && [[ "$FIX_PYTHON" -eq 1 ]]; then
    echo "[0/12] python $PY_VER is outside $PY_MIN-$PY_MAX; --fix-python given, repairing."
    if ! command -v conda >/dev/null 2>&1; then
      echo "ERROR: --fix-python given but conda is not on PATH." >&2
      exit 1
    fi
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
      echo "ERROR: --fix-python needs an ACTIVE conda environment." >&2
      echo "       With none active, conda would install into 'base' and downgrade it." >&2
      echo "       Run:  conda activate bff   then re-run this script." >&2
      exit 1
    fi
    if [[ "$CONDA_PREFIX" == "$(conda info --base 2>/dev/null)" ]]; then
      echo "ERROR: the active environment IS conda 'base' ($CONDA_PREFIX)." >&2
      echo "       Refusing to change the interpreter of base. Use a project env:" >&2
      echo "         conda create -n bff python=$PY_TARGET && conda activate bff" >&2
      exit 1
    fi
    echo "  asking conda for python=$PY_TARGET in $CONDA_PREFIX"
    echo "  NOTE packages installed under python $PY_VER will be orphaned."
    conda install -y -p "$CONDA_PREFIX" "python=$PY_TARGET"
    PY_BIN="$(resolve_py_bin || true)"
    PY_VER=""
    if [[ -n "$PY_BIN" ]]; then PY_VER="$(py_ver_of "$PY_BIN")"; fi
  fi

  if py_in_range "$PY_VER"; then
    echo "[0/12] python $PY_VER at $PY_BIN is within the supported range ($PY_MIN-$PY_MAX)."
  else
    echo "ERROR: python ${PY_VER:-unknown} is outside the supported range ($PY_MIN-$PY_MAX)." >&2
    echo "       interpreter: ${PY_BIN:-<none found>}" >&2
    echo "       numpy<2 (pyproject) has no wheel above cp312; pip would try to COMPILE numpy." >&2
    if [[ "$FIX_PYTHON" -eq 1 ]]; then
      echo "       --fix-python ran but the interpreter is still out of range." >&2
    else
      echo "       Fix it in ONE of these ways, then re-run:" >&2
      echo "         conda install -y -p \"\$CONDA_PREFIX\" python=$PY_TARGET" >&2
      echo "         conda create -n bff python=$PY_TARGET && conda activate bff" >&2
      echo "         ./bootstrap.sh$PHASE_SELECTORS --fix-python   # let this script do it" >&2
    fi
    exit 1
  fi
fi

# SCRIPT_DIR was already resolved above (apply_plmc_patch needs it before the
# arg parser runs); ROOT_DIR is the same directory, so reuse rather than recompute.
ROOT_DIR="$SCRIPT_DIR"

# Prepared-database root. The CoCoPUTs codon-usage table belongs here with every
# other built database, NOT in the repo: it is derived data (~11 MB of upstream
# zips plus the table), it is shared across checkouts, and rare_codon reads it by
# path rather than importing it.
#
# Bio_DBs is NOT a fixed number of levels above this script. ROOT_DIR is the
# `scripts/` directory, so the pre-existing `$ROOT_DIR/../Bio_DBs` at the legacy
# build_db.sh call resolves to <repo>/Bio_DBs, which is not where it lives on any
# machine checked -- this repo sits at <base>/BFF/BioFeatureFactory while Bio_DBs
# sits at <base>/Bio_DBs, three levels up. Searching upward finds it wherever the
# checkout is nested, and --bio-dbs / BFF_BIO_DBS override when it is elsewhere
# entirely. Falls back to scripts/_downloads so a machine with no Bio_DBs still
# bootstraps rather than aborting.
find_bio_dbs() {
  local d="$ROOT_DIR"
  local i
  for i in 1 2 3 4 5; do
    d="$(cd "$d/.." 2>/dev/null && pwd)" || return 1
    [[ -d "$d/Bio_DBs" ]] && { echo "$d/Bio_DBs"; return 0; }
    [[ "$d" == "/" ]] && break
  done
  return 1
}
BIO_DBS_DIR="${BIO_DBS_DIR:-${BFF_BIO_DBS:-}}"
REPO_ROOT="$(cd "$ROOT_DIR/.." && pwd)"
BFF_DIR="$REPO_ROOT/biofeaturefactory"
cd "$ROOT_DIR"

# Only materialise the download cache when a step will actually write to it.
# (GeneSplicer is excluded: step 7 downloads into $GS_DIR, not _downloads.)
if [[ "$DOWNLOAD_UNIREF90" -eq 1 || "$DOWNLOAD_IDMAPPING" -eq 1 ]]; then
  mkdir -p _downloads
fi

# ── Failure collection ──────────────────────────────────────────────────
# `set -e` makes ANY unguarded failure abort the run, so a single step that
# cannot succeed on this host (e.g. the GeneSplicer source build, which needs
# GNU g++ and hard-errors under Apple clang) used to kill every step after it.
# Work steps now record and continue; the run ends with a summary and exits 1 if
# anything failed. PRECONDITIONS stay fatal -- bad arguments, an out-of-range
# interpreter, and the step-1 tool check all run BEFORE any work, and continuing
# past them produces nothing useful.
BOOTSTRAP_FAILURES=()
record_failure() {
  BOOTSTRAP_FAILURES+=("$1")
  echo "  FAILED: $1" >&2
}

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
    return 1
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
# extensions during pip install. Skipped unless the git phase is selected.
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
      pkg_install build-essential || record_failure "step 1b: install build-essential"
    else
      pkg_install gcc gcc-c++ make || record_failure "step 1b: install gcc/g++/make"
    fi
  fi
fi

# ── Step 2: Pip install ──────────────────────────────────────────────────
echo "[2/12] Core Python requirements..."
if [[ "$PIP_INSTALL" -eq 1 ]]; then
  echo "  PIP install -e .[all] (editable install with all optional deps)"
  "$PY_BIN" -m pip install -e "${REPO_ROOT}[all]" \
    || record_failure "step 2: pip install -e .[all]"
  echo "  PIP reinstall pyarrow (numpy ABI fix)"
  "$PY_BIN" -m pip install pyarrow --force-reinstall --quiet \
    || record_failure "step 2: pyarrow reinstall (numpy ABI fix)"
fi

# ── Step 3: mutation_effects (EVmutation + plmc) + adabmDCApy + cg_cotrans ──
echo "[3/12] mutation_effects module dependencies..."
if [[ "$CLONE_EVMUTATION" -eq 1 ]]; then
  clone_or_update "https://github.com/debbiemarkslab/EVmutation.git" "$BFF_DIR/mutation_effects/EVmutation" \
    || record_failure "step 3: clone EVmutation"
  clone_or_update "https://github.com/debbiemarkslab/plmc.git" "$BFF_DIR/mutation_effects/plmc" \
    || record_failure "step 3: clone plmc"
  echo "  PATCH plmc"
  apply_plmc_patch "$SCRIPT_DIR/patches/plmc_proximal_group_lasso.patch" \
                   "proximal-group-lasso" \
    || record_failure "step 3: plmc patch proximal-group-lasso"
  apply_plmc_patch "$SCRIPT_DIR/patches/plmc_case_sensitive_custom_alphabet.patch" \
                   "case-sensitive-custom-alphabet (q=64 codon Potts)" \
    || record_failure "step 3: plmc patch case-sensitive-custom-alphabet (q=64 codon Potts)"
  if [[ "$BUILD_PLMC" -eq 1 ]]; then
    echo "  BUILD plmc"
    # all-openmp needs GNU gcc; Apple clang rejects -fopenmp, hence the fallback
    # to the single-threaded target. Only a failure of BOTH is recorded.
    (
      cd "$BFF_DIR/mutation_effects/plmc"
      make all-openmp || make all
    ) || record_failure "step 3: build plmc"
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
  ) || record_failure "step 3: clone adabmDCApy"
  # Editable install (pip install -e) + adabmDCA version verification are handled
  # centrally in step [8b], so that env-phase (which does not clone here) and
  # git-phase (which does not run the core pip install) both install it.
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
  clone_or_update "https://github.com/Eryk96/NetSurfP-3.0.git" "$BFF_DIR/NetSurfP3/nsp3" \
    || record_failure "step 4: clone NetSurfP-3.0"
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
  else
    conda_install_pinned miranda miranda || true
  fi
fi

# ── Step 6b: mmseqs2 (EVmutation codon MSA) ─────────────────────────────
echo "[6b/12] mmseqs2..."
if [[ "$INSTALL_MMSEQS2" -eq 1 ]]; then
  if command -v mmseqs >/dev/null 2>&1; then
    echo "  OK mmseqs2 already on PATH"
  else
    conda_install_pinned mmseqs2 mmseqs2 || true
  fi
fi

# ── Step 6c: HMMER / jackhmmer (EVmutation protein MSA) ────────────────
echo "[6c/12] HMMER (jackhmmer)..."
if [[ "$INSTALL_HMMER" -eq 1 ]]; then
  if command -v jackhmmer >/dev/null 2>&1; then
    echo "  OK jackhmmer already on PATH"
  else
    conda_install_pinned hmmer hmmer || true
  fi
fi

# ── Step 6d: SpliceAI ───────────────────────────────────────────────────
# Placed here, with the other conda installs, rather than in the steps 10-12
# block where the spliceai CHECK lives: that block is skipped unless the git phase runs,
# and the env phase is documented as pip/conda installs. A conda
# package belongs in the phase that installs conda packages.
echo "[6d/12] SpliceAI (own conda env)..."
# SpliceAI gets its OWN environment. Two independent reasons, both measured:
#
#  1. Abseil collision. pyarrow and tensorflow each vendor their own Abseil, and
#     whichever loads first owns the symbols. The `spliceai` console script
#     imports pandas -> pyarrow BEFORE tensorflow, so pyarrow wins, and tf.data's
#     prefetch CondVar inside Keras model.predict is then waited on by libarrow's
#     absl and signalled by tensorflow's. It never wakes. Measured on a real SMN2
#     run in the shared `bff` env: 6.31 s of CPU in 38 MINUTES of wall clock at
#     0.0% CPU, blocked in
#       PrefetchDatasetOp::GetNextInternal -> absl::CondVar::WaitCommon
#         -> AbslInternalPerThreadSemWait (in libarrow.2500.dylib)
#     Reduced to a 6-line repro, run both ways in that env:
#       import pyarrow    -> model.predict()  TIMEOUT (deadlock)
#       import tensorflow -> model.predict()  0.06 s, and 0.03 s AFTER pyarrow
#     pandas 3.0.5 REQUIRES pyarrow, so pyarrow cannot simply be dropped from a
#     shared env. bin/main.nf works around it by importing tensorflow first; an
#     env without pyarrow removes the hazard rather than sequencing around it.
#
#  2. setuptools. spliceai/__init__.py:2 does `from pkg_resources import
#     get_distribution`, and pkg_resources was REMOVED in setuptools 81. Pinning
#     the SHARED env to setuptools<81 holds all of BFF on a deprecated setuptools
#     line to satisfy one tool. In its own env the pin costs nothing.
#
# --exclude-spliceai-env puts it back in the active env (and back in reach of both
# problems). --spliceai-env NAME changes the env name.
if [[ "$INSTALL_SPLICEAI" -eq 1 ]]; then
  CONDA_BIN=""
  command -v conda >/dev/null 2>&1 && CONDA_BIN=conda
  [[ -z "$CONDA_BIN" ]] && command -v mamba >/dev/null 2>&1 && CONDA_BIN=mamba

  if [[ "$SPLICEAI_OWN_ENV" -eq 1 && -n "$CONDA_BIN" ]]; then
    # `conda env list` is the only reliable existence test: `conda run -n` on a
    # missing env exits non-zero for several unrelated reasons too.
    if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$SPLICEAI_ENV_NAME"; then
      echo "  OK env '$SPLICEAI_ENV_NAME' already exists"
    else
      echo "  CREATE env '$SPLICEAI_ENV_NAME' (python=${PY_VER:-3.11})"
      "$CONDA_BIN" create -y -n "$SPLICEAI_ENV_NAME" "python=${PY_VER:-3.11}" \
        || record_failure "step 6d: conda create -n $SPLICEAI_ENV_NAME"
    fi
    echo "  INSTALL spliceai + setuptools<81 into '$SPLICEAI_ENV_NAME'"
    "$CONDA_BIN" install -y -n "$SPLICEAI_ENV_NAME" -c bioconda spliceai "python=${PY_VER:-3.11}" \
      || record_failure "step 6d: conda install spliceai -n $SPLICEAI_ENV_NAME"
    "$CONDA_BIN" run -n "$SPLICEAI_ENV_NAME" python -m pip install --quiet "setuptools>=77.0,<81" \
      || record_failure "step 6d: setuptools<81 in $SPLICEAI_ENV_NAME"

    # Verify what actually breaks, not just that the CLI exists. `spliceai --help`
    # is argparse only -- it exits before TensorFlow loads, so it passed on the
    # very env that then deadlocked. This runs the real failure mode with a
    # timeout instead.
    echo "  CHECK pkg_resources + Keras predict after pyarrow (the deadlock)"
    "$CONDA_BIN" run -n "$SPLICEAI_ENV_NAME" python - <<'SPLICEAI_CHECK' \
      || record_failure "step 6d: '$SPLICEAI_ENV_NAME' fails the import-order/pkg_resources check"
import faulthandler, sys, threading
faulthandler.dump_traceback_later(180, exit=True)   # a hang is a FAILURE, not a wait
try:
    import pkg_resources                             # spliceai/__init__.py:2
except Exception as exc:
    print(f"  FAIL pkg_resources: {exc}"); sys.exit(1)
try:
    import pyarrow                                   # only if something pulled it in
    print(f"  WARN pyarrow {pyarrow.__version__} present in this env")
except ImportError:
    print("  OK  pyarrow absent (Abseil collision impossible)")
try:
    import numpy as np, tensorflow as tf
    m = tf.keras.Sequential([tf.keras.layers.Dense(2, input_shape=(4,))])
    m.predict(np.zeros((4, 4), dtype="float32"), verbose=0)
    print("  OK  keras predict completed")
except Exception as exc:
    print(f"  FAIL keras predict: {exc}"); sys.exit(1)
faulthandler.cancel_dump_traceback_later()
SPLICEAI_CHECK
    echo "  OK spliceai env ready: $SPLICEAI_ENV_NAME"
  else
    if [[ "$SPLICEAI_OWN_ENV" -eq 1 ]]; then
      echo "  WARN conda/mamba not found; installing into the active env instead"
    else
      echo "  NOTE --exclude-spliceai-env: installing into the active env"
    fi
    if command -v spliceai >/dev/null 2>&1; then
      echo "  OK spliceai already on PATH"
    else
      conda_install_pinned spliceai spliceai || true
    fi
    if command -v spliceai >/dev/null 2>&1; then
      if ! "$PY_BIN" -c 'import pkg_resources' >/dev/null 2>&1; then
        echo "  FIX  installing setuptools<81 (spliceai needs pkg_resources)"
        "$PY_BIN" -m pip install --quiet "setuptools>=77.0,<81" \
          || record_failure "step 6d: pip install setuptools<81 (spliceai needs pkg_resources)"
      fi
      echo "  WARN shared env: spliceai is exposed to the pyarrow/TensorFlow Abseil"
      echo "       deadlock. bin/main.nf imports tensorflow first to avoid it."
    fi
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
    if ! download_file "ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz" "$GS_TAR"; then
      record_failure "step 7: download GeneSplicer.tar.gz"
    elif ! tar -xzf "$GS_TAR" -C "$GS_DIR"; then
      record_failure "step 7: extract GeneSplicer.tar.gz"
    else
      rm -f "$GS_TAR"
      if [[ -d "$GS_SRC/sources" ]]; then
        echo "  BUILD GeneSplicer from source"
        # genesplicer.cpp:58 declares `main` with no return type. GNU g++ accepts
        # it (warning: ISO C++ forbids declaration of 'main' with no type) and
        # emits the object; Apple clang makes it a hard error that neither
        # -std=gnu++98 nor -fpermissive suppresses. Expect this to fail on macOS
        # and succeed on Linux, which is the deployment target.
        (cd "$GS_SRC/sources" && make) || record_failure "step 7: build GeneSplicer (needs GNU g++; fails under Apple clang)"
      else
        record_failure "step 7: GeneSplicer sources/ missing after extraction"
      fi
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
      pkg_install openjdk-17-jdk || record_failure "step 7b: install openjdk-17-jdk"
    else
      pkg_install java-17-openjdk-devel || record_failure "step 7b: install java-17-openjdk-devel"
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
      record_failure "step 7b: locate the Java 17 binary after install"
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
      conda install -y -c bioconda "nextflow>=${NEXTFLOW_MIN}" \
        || record_failure "step 7b: conda install nextflow"
    elif command -v mamba >/dev/null 2>&1; then
      echo "  MAMBA install nextflow (conda not found)"
      mamba install -y -c bioconda "nextflow>=${NEXTFLOW_MIN}" \
        || record_failure "step 7b: mamba install nextflow"
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
        record_failure "step 7b: nextflow installer produced no binary"
        nf_bin=""
      fi

      # Everything below moves $nf_bin onto PATH, so it must be skipped when the
      # installer produced nothing -- otherwise `mv` fails and set -e aborts the run.
      if [[ -n "$nf_bin" ]]; then
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
  clone_or_update "https://github.com/google-deepmind/alphafold3.git" "$BFF_DIR/alphafold3/alphafold3" \
    || record_failure "step 8: clone AlphaFold3"
fi

# ── Step 8b: Editable installs of python-based cloned repos ─────────────
# Runs whenever the env or git phase is selected. Decoupled from the core
# `pip install -e .[all]` (PIP_INSTALL) and from the per-repo CLONE_* flags so:
#   - env-phase installs whatever python repos are already cloned, and warns for
#     any repo that is absent (nothing was cloned this run);
#   - git-phase installs the repos it just cloned, and warns if the core BFF env
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

  # git-phase (or any selection that skipped the core pip install): the BFF package and
  # its dependencies may be missing. Warn but still attempt the editable installs.
  # Probe a real leaf module: `import biofeaturefactory` succeeds even on a
  # broken/absent install, since the top-level package is an empty namespace.
  if [[ "$PIP_INSTALL" -eq 0 ]] && ! "$PY_BIN" -c "import biofeaturefactory.lib.utility" >/dev/null 2>&1; then
    echo "  WARN BioFeatureFactory core env not detected ('import biofeaturefactory.lib.utility' failed)."
    echo "       Run './bootstrap.sh env-phase' or the full bootstrap to install python dependencies."
  fi

  for entry in "${EDITABLE_REPOS[@]}"; do
    IFS='|' read -r er_label er_rel er_suffix er_marker <<< "$entry"
    er_repo_dir="$BFF_DIR/$er_rel"
    er_install_dir="$er_repo_dir$er_suffix"
    if [[ -d "$er_repo_dir" && -f "$er_install_dir/$er_marker" ]]; then
      echo "  PIP install -e $er_label ($er_install_dir)"
      "$PY_BIN" -m pip install -e "$er_install_dir" \
        || record_failure "step 8b: pip install -e $er_label"
    else
      echo "  WARN $er_label repo not present at $er_repo_dir"
      echo "       Run the full bootstrap, or './bootstrap.sh git-phase', to install the proper repos."
    fi
  done

  # adabmDCA version + CLI verification (warn-only), after the editable install above.
  if [[ "$CLONE_ADABMDCA" -eq 1 || -d "$BFF_DIR/mutation_effects/adabmDCApy" ]]; then
    installed_ver="$("$PY_BIN" -m pip show adabmDCA 2>/dev/null | awk -F': ' '/^Version:/ {print $2}')"
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

# ── Step 8c: Verify the installed package layout ────────────────────────
# The package is layered lib -> core -> pipelines (see biofeaturefactory/core/__init__.py):
#   lib/    import-only modules, no CLI
#   core/   input producers, driven by CLI (including from Nextflow), not imported
#           by the pipelines
#   rest/   pipelines that consume core's output
# A bare `import biofeaturefactory` passes even when that layout is broken, because
# the top-level package is an empty namespace. Import each layer explicitly, then
# resolve the console scripts declared in pyproject.toml the same way a shell
# invocation would. Warn-only: a missing OPTIONAL backend (torch, ViennaRNA, an
# unfetched clone) must not abort a bootstrap that has already done real work.
echo "[8c/12] Verifying installed package layout..."
if [[ "$PIP_INSTALL" -eq 1 || "$INSTALL_EDITABLE_REPOS" -eq 1 ]]; then
  "$PY_BIN" - <<'PYCHECK' || echo "  WARN layout verification reported problems (see above)"
import importlib, sys
from importlib.metadata import distributions, entry_points

LIB = ["primitives", "utility", "annotation", "msa", "codon_metrics", "dtu_outputs"]
CORE = ["variant_mapping", "vcf_converter", "msa_generation_pipeline", "codon_msa_pipeline"]

bad = 0
for layer, mods in (("lib", LIB), ("core", CORE)):
    ok = []
    for m in mods:
        name = f"biofeaturefactory.{layer}.{m}"
        try:
            importlib.import_module(name)
            ok.append(m)
        except Exception as e:
            bad += 1
            print(f"  FAIL {name}: {type(e).__name__}: {e}")
    print(f"  OK  {layer}: {len(ok)}/{len(mods)} modules import")

# A source tree that has been built in place carries BOTH an in-tree
# biofeaturefactory.egg-info and the site-packages dist-info, and a stale
# dist-info from an older pyproject will shadow the current one. Name the
# distribution each result came from so a stale shim is not mistaken for a
# broken pipeline.
for dist in distributions():
    if (dist.metadata["Name"] or "").lower() == "biofeaturefactory":
        n = len([e for e in dist.entry_points if e.name.startswith("bff-")])
        print(f"  metadata: {dist._path}  ({n} bff-* scripts)")

try:
    eps = sorted(entry_points(group="console_scripts"), key=lambda e: e.name)
except TypeError:  # Python 3.9 API
    eps = sorted(entry_points().get("console_scripts", []), key=lambda e: e.name)
eps = [e for e in eps if e.name.startswith("bff-")]
if not eps:
    print("  WARN no bff-* console scripts found; is the package installed?")
else:
    good = 0
    for e in eps:
        try:
            e.load()
            good += 1
        except Exception as ex:
            bad += 1
            print(f"  FAIL {e.name} -> {e.value}: {type(ex).__name__}: {ex}")
    print(f"  OK  console scripts: {good}/{len(eps)} resolve")

sys.exit(1 if bad else 0)
PYCHECK
else
  echo "  SKIP (no python env step ran in this mode)"
fi

# ── Step 9: FTP datasets ────────────────────────────────────────────────
echo "[9/12] Optional FTP datasets..."
if [[ "$DOWNLOAD_UNIREF90" -eq 1 ]]; then
  download_file \
    "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz" \
    "$ROOT_DIR/_downloads/uniref90.fasta.gz" \
    || record_failure "step 9: download uniref90.fasta.gz"
fi
if [[ "$DOWNLOAD_IDMAPPING" -eq 1 ]]; then
  download_file \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz" \
    "$ROOT_DIR/_downloads/idmapping.dat.gz" \
    || record_failure "step 9: download idmapping.dat.gz"
fi

# ── Step 9b: CoCoPUTs codon usage table -- MOVED ────────────────────────
# The table is built by scripts/build_db.sh (step 9b there), not here. It is a
# prepared database like every other artifact under Bio_DBs, and building it from
# bootstrap gave that directory two owners with two different ideas of where it
# lives. BUILD_COCOPUTS_CUT/--exclude-cocoputs still control it -- they are
# forwarded to build_db.sh as SKIP_COCOPUTS at step 12. Both flags sit in
# PHASE_DB_FLAGS and both default to 1, so any invocation that used to reach the
# build here reaches it there.

# ── Steps 10-12: summaries over what the git phase clones/builds ─────────
if [[ "$PHASE_GIT" -eq 1 ]]; then
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
  # `|| true` on every line: these are `cond && echo` probes, and a FALSE one is
  # normal (the thing simply is not installed). Without it the last probe failing
  # makes this whole if-block return 1, which set -e turns into an abort.
  summary_probe() {  # <label> <test...>
    local label="$1"; shift
    if "$@" >/dev/null 2>&1; then echo "  OK $label"; else echo "  --  $label (absent)"; fi
  }
  summary_probe "EVmutation clone"     test -d "$BFF_DIR/mutation_effects/EVmutation"
  summary_probe "plmc clone"           test -d "$BFF_DIR/mutation_effects/plmc"
  summary_probe "plmc binary"          test -x "$BFF_DIR/mutation_effects/plmc/bin/plmc"
  summary_probe "adabmDCApy clone"     test -d "$BFF_DIR/mutation_effects/adabmDCApy"
  summary_probe "adabmDCA CLI on PATH" command -v adabmDCA
  summary_probe "cg_cotrans"           test -d "$BFF_DIR/rare_codon/cg_cotrans"
  summary_probe "NetSurfP3 clone"      test -d "$BFF_DIR/NetSurfP3/nsp3"
  summary_probe "signalp6 on PATH"     command -v signalp6
  summary_probe "miranda on PATH"      command -v miranda
  summary_probe "genesplicer build"    test -x "$BFF_DIR/genesplicer/GeneSplicer/sources/genesplicer"
  summary_probe "AlphaFold3 clone"     test -d "$BFF_DIR/alphafold3/alphafold3"
fi

if [[ "$RUN_BUILD_DB" -eq 1 ]]; then
  DB_SCRIPT="$ROOT_DIR/build_db.sh"
  if [[ -x "$DB_SCRIPT" ]]; then
    echo "  RUN scripts/build_db.sh"
    # DB_ROOT and SKIP_COCOPUTS are the two things bootstrap owns about the
    # database build: --bio-dbs picks the root, --exclude-cocoputs turns off the
    # codon-usage table. Everything else build_db.sh decides for itself.
    if [[ -n "$BIO_DBS_DIR" ]]; then export DB_ROOT="$BIO_DBS_DIR"; fi
    if [[ "$BUILD_COCOPUTS_CUT" -eq 0 ]]; then export SKIP_COCOPUTS=1; fi
    "$DB_SCRIPT" || record_failure "step 12: scripts/build_db.sh"
  else
    echo "  WARN build_db.sh not found/executable at $DB_SCRIPT"
  fi
fi

# ── Final report ────────────────────────────────────────────────────────
if [[ "${#BOOTSTRAP_FAILURES[@]}" -gt 0 ]]; then
  echo
  echo "════════════════════════════════════════════════════════════════"
  echo "Bootstrap finished with ${#BOOTSTRAP_FAILURES[@]} failure(s):"
  for _f in "${BOOTSTRAP_FAILURES[@]}"; do
    echo "  ✗ $_f"
  done
  echo "════════════════════════════════════════════════════════════════"
  echo "Every other step completed. Re-run with the matching --exclude-* flags"
  echo "to skip the failures, or fix them and re-run the same phase."
  exit 1
fi

echo "Done. All selected steps succeeded."