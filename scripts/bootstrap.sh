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
#   ./bootstrap.sh pip-only        # Only pip install
#   ./bootstrap.sh git-only        # Only git clones and builds
#   ./bootstrap.sh db-only         # Only database downloads
#   ./bootstrap.sh git-only --exclude-evmutation --exclude-cg-cotrans
#   ./bootstrap.sh db-only --exclude-uniref90
#
# Subcommands:
#   pip-only       Only run pip install (step 2).
#   git-only       Only run git clones, builds, and conda installs (steps 3-8).
#   db-only        Only run FTP downloads and build_db (steps 9, 12).
#   (none)         Run everything.
#
# Exclude flags (fine-grained control within any mode):
#   --exclude-pip-install     Skip pip install.
#   --exclude-evmutation      Skip cloning EVmutation and plmc.
#   --exclude-build-plmc      Skip compiling plmc (clone still happens).
#   --exclude-cg-cotrans      Skip downloading cg_cotrans.
#   --exclude-netsurfp3       Skip cloning NetSurfP-3.0.
#   --exclude-signalp         Skip cloning SignalP 6.0.
#   --exclude-miranda         Skip conda install of miranda.
#   --exclude-genesplicer     Skip downloading/building GeneSplicer.
#   --exclude-clone-af3       Skip cloning AlphaFold3.
#   --exclude-uniref90        Skip UniRef90 FTP download.
#   --exclude-idmapping       Skip UniProt idmapping FTP download.
#   --exclude-build-db        Skip calling build_db.sh.

# --- Defaults: everything on ---
PIP_INSTALL=1
CLONE_EVMUTATION=1
BUILD_PLMC=1
DOWNLOAD_CG_COTRANS=1
CLONE_NETSURFP3=1
INSTALL_SIGNALP=1
INSTALL_MIRANDA=1
BUILD_GENESPLICER=1
CLONE_AF3=1
DOWNLOAD_UNIREF90=1
DOWNLOAD_IDMAPPING=1
RUN_BUILD_DB=1

# --- Parse subcommand ---
SUBCOMMAND=""
ARGS=()
for arg in "$@"; do
  case "$arg" in
    pip-only|git-only|db-only) SUBCOMMAND="$arg" ;;
    *) ARGS+=("$arg") ;;
  esac
done

# Apply subcommand: disable groups not selected
case "$SUBCOMMAND" in
  pip-only)
    CLONE_EVMUTATION=0; BUILD_PLMC=0; DOWNLOAD_CG_COTRANS=0
    CLONE_NETSURFP3=0; INSTALL_SIGNALP=0
    BUILD_GENESPLICER=0; CLONE_AF3=0
    DOWNLOAD_UNIREF90=0; DOWNLOAD_IDMAPPING=0; RUN_BUILD_DB=0
    ;;
  git-only)
    PIP_INSTALL=0
    DOWNLOAD_UNIREF90=0; DOWNLOAD_IDMAPPING=0; RUN_BUILD_DB=0
    ;;
  db-only)
    PIP_INSTALL=0
    CLONE_EVMUTATION=0; BUILD_PLMC=0; DOWNLOAD_CG_COTRANS=0
    CLONE_NETSURFP3=0; INSTALL_SIGNALP=0; INSTALL_MIRANDA=0
    BUILD_GENESPLICER=0; CLONE_AF3=0
    ;;
esac

# --- Parse exclude flags ---
set -- "${ARGS[@]+"${ARGS[@]}"}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exclude-pip-install)   PIP_INSTALL=0 ;;
    --exclude-evmutation)    CLONE_EVMUTATION=0 ;;
    --exclude-build-plmc)    BUILD_PLMC=0 ;;
    --exclude-cg-cotrans)    DOWNLOAD_CG_COTRANS=0 ;;
    --exclude-netsurfp3)     CLONE_NETSURFP3=0 ;;
    --exclude-signalp)       INSTALL_SIGNALP=0 ;;
    --exclude-miranda)       INSTALL_MIRANDA=0 ;;
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
if [[ "$SUBCOMMAND" == "pip-only" && "$PIP_INSTALL" -eq 0 ]]; then
  echo "ERROR: pip-only with --exclude-pip-install is contradictory." >&2
  exit 1
fi
if [[ "$SUBCOMMAND" == "db-only" && "$DOWNLOAD_UNIREF90" -eq 0 && "$DOWNLOAD_IDMAPPING" -eq 0 && "$RUN_BUILD_DB" -eq 0 ]]; then
  echo "ERROR: db-only with all database steps excluded leaves nothing to do." >&2
  exit 1
fi
if [[ "$SUBCOMMAND" == "git-only" && "$CLONE_EVMUTATION" -eq 0 && "$BUILD_PLMC" -eq 0 && "$DOWNLOAD_CG_COTRANS" -eq 0 && "$CLONE_NETSURFP3" -eq 0 && "$INSTALL_SIGNALP" -eq 0 && "$INSTALL_MIRANDA" -eq 0 && "$BUILD_GENESPLICER" -eq 0 && "$CLONE_AF3" -eq 0 ]]; then
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

# ── Step 2: Pip install ──────────────────────────────────────────────────
echo "[2/12] Core Python requirements..."
if [[ "$PIP_INSTALL" -eq 1 ]]; then
  echo "  PIP install -e .[all] (editable install with all optional deps)"
  python -m pip install -e "${REPO_ROOT}[all]"
  echo "  PIP reinstall pyarrow (numpy ABI fix)"
  python -m pip install pyarrow --force-reinstall --quiet
fi

# ── Step 3: EVmutation + plmc + cg_cotrans ───────────────────────────────
echo "[3/12] EVmutation module dependencies..."
if [[ "$CLONE_EVMUTATION" -eq 1 ]]; then
  clone_or_update "https://github.com/debbiemarkslab/EVmutation.git" "$BFF_DIR/EVmutation/EVmutation"
  clone_or_update "https://github.com/debbiemarkslab/plmc.git" "$BFF_DIR/EVmutation/plmc"

  if [[ "$BUILD_PLMC" -eq 1 ]]; then
    echo "  BUILD plmc"
    (
      cd "$BFF_DIR/EVmutation/plmc"
      if [[ "$(uname -s)" == "Darwin" ]]; then
        if command -v brew >/dev/null 2>&1; then
          brew list libomp >/dev/null 2>&1 || brew install libomp
        fi
        make all-mac-openmp || make all-mac || make all
      else
        make all-openmp || make all
      fi
    )
  fi
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
  if [[ "$PIP_INSTALL" -eq 1 ]] && [[ -f "$BFF_DIR/NetSurfP3/nsp3/requirements.txt" ]]; then
    echo "  PIP install NetSurfP3 requirements"
    python -m pip install -r "$BFF_DIR/NetSurfP3/nsp3/requirements.txt"
  fi
fi

# ── Step 5: SignalP 6.0 ─────────────────────────────────────────────────
echo "[5/12] SignalP 6.0 (netNglyc dependency)..."
if [[ "$INSTALL_SIGNALP" -eq 1 ]]; then
  clone_or_update "https://github.com/fteufel/signalp-6.0" "$BFF_DIR/netNglyc/signalp-6.0"
  if [[ "$PIP_INSTALL" -eq 1 ]] && [[ -f "$BFF_DIR/netNglyc/signalp-6.0/requirements.txt" ]]; then
    echo "  PIP install SignalP 6.0 requirements"
    python -m pip install -r "$BFF_DIR/netNglyc/signalp-6.0/requirements.txt"
  fi
fi

# ── Step 6: Miranda ─────────────────────────────────────────────────────
echo "[6/12] Miranda (conda)..."
if [[ "$INSTALL_MIRANDA" -eq 1 ]]; then
  if command -v miranda >/dev/null 2>&1; then
    echo "  OK miranda already on PATH"
  elif command -v mamba >/dev/null 2>&1; then
    echo "  MAMBA install miranda"
    mamba install -y -c bioconda miranda
  elif command -v conda >/dev/null 2>&1; then
    echo "  CONDA install miranda"
    conda install -y -c bioconda miranda
  else
    echo "  WARN conda/mamba not found; install miranda manually: conda install -c bioconda miranda"
  fi
fi

# ── Step 7: GeneSplicer ─────────────────────────────────────────────────
echo "[7/12] GeneSplicer binary..."
if [[ "$BUILD_GENESPLICER" -eq 1 ]]; then
  if command -v genesplicer >/dev/null 2>&1; then
    echo "  OK genesplicer already on PATH"
  elif command -v mamba >/dev/null 2>&1; then
    echo "  MAMBA install genesplicer"
    mamba install -y -c bioconda genesplicer
  elif command -v conda >/dev/null 2>&1; then
    echo "  CONDA install genesplicer"
    conda install -y -c bioconda genesplicer
  else
    echo "  WARN conda/mamba not found; falling back to source build"
    GS_DIR="$BFF_DIR/genesplicer"
    GS_TAR="$GS_DIR/GeneSplicer.tar.gz"
    GS_SRC="$GS_DIR/GeneSplicer"
    if [[ -x "$GS_SRC/sources/genesplicer" ]]; then
      echo "  EXISTS $GS_SRC/sources/genesplicer"
    else
      mkdir -p "$GS_DIR"
      download_file "ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz" "$GS_TAR"
      tar -xzf "$GS_TAR" -C "$GS_DIR"
      rm -f "$GS_TAR"
      if [[ -d "$GS_SRC/sources" ]]; then
        echo "  BUILD GeneSplicer from source"
        (cd "$GS_SRC/sources" && make)
      else
        echo "  WARN GeneSplicer source directory not found after extraction"
      fi
    fi
  fi
fi

# ── Step 8: AlphaFold3 ──────────────────────────────────────────────────
echo "[8/12] AlphaFold3 upstream (optional clone)..."
if [[ "$CLONE_AF3" -eq 1 ]]; then
  clone_or_update "https://github.com/google-deepmind/alphafold3.git" "$BFF_DIR/alphafold3/alphafold3"
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
    echo "  WARN spliceai not found (required by spliceai/README.md)"
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
    SignalP 6.0 is cloned automatically to netNglyc/signalp-6.0/
    Patch NetNGlyc tcsh SIGNALP path to:
      biofeaturefactory/netNglyc/bin/signalp6_adapter
  - netphos/: NetPhos 3.1 + APE (DTU academic license), requires tcsh
  - netMHC/: NetMHCpan 4.1 (DTU academic license)
  - alphafold3/: NVIDIA GPU stack + Docker + AF3 model weights
EOF

  echo "[12/12] Summary checks..."
  [[ -d "$BFF_DIR/EVmutation/EVmutation" ]]     && echo "  OK EVmutation clone"
  [[ -d "$BFF_DIR/EVmutation/plmc" ]]           && echo "  OK plmc clone"
  [[ -d "$BFF_DIR/rare_codon/cg_cotrans" ]]     && echo "  OK cg_cotrans"
  [[ -d "$BFF_DIR/NetSurfP3/nsp3" ]]            && echo "  OK NetSurfP3 clone"
  [[ -d "$BFF_DIR/netNglyc/signalp-6.0" ]]      && echo "  OK SignalP 6.0 clone"
  command -v miranda >/dev/null 2>&1             && echo "  OK miranda on PATH"
  command -v genesplicer >/dev/null 2>&1             && echo "  OK genesplicer on PATH"
  [[ -d "$BFF_DIR/alphafold3/alphafold3" ]]     && echo "  OK AlphaFold3 clone"
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