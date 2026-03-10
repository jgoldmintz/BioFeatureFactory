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
#   ./bootstrap.sh
#   ./bootstrap.sh --exclude-uniref90 --exclude-idmapping
#
# Flags:
#   --exclude-pip-install     Skip Python requirements installation.
#   --exclude-build-plmc      Skip compiling plmc.
#   --exclude-uniref90        Skip UniRef90 FTP download (very large).
#   --exclude-idmapping       Skip UniProt idmapping FTP download.
#   --exclude-clone-af3       Skip cloning DeepMind AlphaFold3 upstream repo.
#   --exclude-signalp         Skip cloning SignalP 6.0 (netNglyc dependency).
#   --exclude-miranda         Skip conda install of miranda.
#   --exclude-genesplicer     Skip downloading/building GeneSplicer.
#   --exclude-build-db        Skip calling Bio_DBs/build_db.sh at the end.
#
# Legacy include flags are accepted as no-ops for compatibility:
#   --pip-install --build-plmc --download-uniref90 --download-idmapping --clone-alphafold3

PIP_INSTALL=1
BUILD_PLMC=1
DOWNLOAD_UNIREF90=1
DOWNLOAD_IDMAPPING=1
CLONE_AF3=1
INSTALL_SIGNALP=1
INSTALL_MIRANDA=1
BUILD_GENESPLICER=1
RUN_BUILD_DB=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exclude-pip-install)   PIP_INSTALL=0 ;;
    --exclude-build-plmc)    BUILD_PLMC=0 ;;
    --exclude-uniref90)      DOWNLOAD_UNIREF90=0 ;;
    --exclude-idmapping)     DOWNLOAD_IDMAPPING=0 ;;
    --exclude-clone-af3)     CLONE_AF3=0 ;;
    --exclude-signalp)       INSTALL_SIGNALP=0 ;;
    --exclude-miranda)       INSTALL_MIRANDA=0 ;;
    --exclude-genesplicer)   BUILD_GENESPLICER=0 ;;
    --exclude-build-db)      RUN_BUILD_DB=0 ;;
    --pip-install|--build-plmc|--download-uniref90|--download-idmapping|--clone-alphafold3)
      # legacy include flags: defaults are already enabled
      ;;
    -h|--help)
      sed -n '1,42p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
  shift
done

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
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

echo "[1/12] Validating core tools..."
require_cmd git
require_cmd tar

echo "[2/12] Core Python requirements..."
REPO_ROOT="$(cd "$ROOT_DIR/.." && pwd)"
if [[ "$PIP_INSTALL" -eq 1 ]]; then
  echo "  PIP install requirements.txt"
  python -m pip install -r "$REPO_ROOT/requirements.txt"
  echo "  PIP install -e . (editable install)"
  python -m pip install -e "$REPO_ROOT"
  echo "  PIP reinstall pyarrow (numpy ABI fix)"
  python -m pip install pyarrow --force-reinstall --quiet
fi

echo "[3/12] EVmutation module dependencies..."
clone_or_update "https://github.com/debbiemarkslab/EVmutation.git" "$ROOT_DIR/EVmutation/EVmutation"
clone_or_update "https://github.com/debbiemarkslab/plmc.git" "$ROOT_DIR/EVmutation/plmc"

if [[ "$BUILD_PLMC" -eq 1 ]]; then
  echo "  BUILD plmc"
  (
    cd "$ROOT_DIR/EVmutation/plmc"
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

echo "[3/12] rare_codon module dependency (cg_cotrans)..."
CG_DIR="$ROOT_DIR/rare_codon/cg_cotrans"
if [[ -d "$CG_DIR" ]]; then
  echo "  EXISTS $CG_DIR"
else
  TMP_TAR="$ROOT_DIR/rare_codon/cg_cotrans.tar.gz"
  download_file "https://shakhnovich.faculty.chemistry.harvard.edu/files/shakhnovich/files/cg_cotrans.tar.gz" "$TMP_TAR"
  mkdir -p "$ROOT_DIR/rare_codon"
  tar -xzf "$TMP_TAR" -C "$ROOT_DIR/rare_codon"
  rm -f "$TMP_TAR"
fi

echo "[4/12] NetSurfP3 module dependency..."
clone_or_update "https://github.com/Eryk96/NetSurfP-3.0.git" "$ROOT_DIR/NetSurfP3/nsp3"
if [[ "$PIP_INSTALL" -eq 1 ]] && [[ -f "$ROOT_DIR/NetSurfP3/nsp3/requirements.txt" ]]; then
  echo "  PIP install NetSurfP3 requirements"
  python -m pip install -r "$ROOT_DIR/NetSurfP3/nsp3/requirements.txt"
fi

echo "[5/12] SignalP 6.0 (netNglyc dependency)..."
if [[ "$INSTALL_SIGNALP" -eq 1 ]]; then
  clone_or_update "https://github.com/fteufel/signalp-6.0" "$ROOT_DIR/netNglyc/signalp-6.0"
  if [[ "$PIP_INSTALL" -eq 1 ]] && [[ -f "$ROOT_DIR/netNglyc/signalp-6.0/requirements.txt" ]]; then
    echo "  PIP install SignalP 6.0 requirements"
    python -m pip install -r "$ROOT_DIR/netNglyc/signalp-6.0/requirements.txt"
  fi
fi

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

echo "[7/12] GeneSplicer binary..."
if [[ "$BUILD_GENESPLICER" -eq 1 ]]; then
  GS_DIR="$ROOT_DIR/genesplicer"
  GS_TAR="$GS_DIR/GeneSplicer.tar.gz"
  GS_SRC="$GS_DIR/GeneSplicer"
  if [[ -x "$GS_SRC/genesplicer" ]]; then
    echo "  EXISTS $GS_SRC/genesplicer"
  else
    mkdir -p "$GS_DIR"
    download_file "ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz" "$GS_TAR"
    tar -xzf "$GS_TAR" -C "$GS_DIR"
    rm -f "$GS_TAR"
    if [[ -d "$GS_SRC" ]]; then
      echo "  BUILD GeneSplicer"
      (cd "$GS_SRC" && make)
    else
      echo "  WARN GeneSplicer source directory not found after extraction"
    fi
  fi
fi

echo "[8/12] AlphaFold3 upstream (optional clone)..."
if [[ "$CLONE_AF3" -eq 1 ]]; then
  clone_or_update "https://github.com/google-deepmind/alphafold3.git" "$ROOT_DIR/alphafold3/alphafold3"
fi

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
  - netNglyc/: NetNGlyc 1.0 (DTU academic license)
    SignalP 6.0 is cloned automatically to netNglyc/signalp-6.0/
    Patch NetNGlyc tcsh SIGNALP path to:
      biofeaturefactory/netNglyc/bin/signalp6_adapter
  - netphos/: NetPhos 3.1 + APE (DTU academic license), requires tcsh
  - netMHC/: NetMHCpan 4.1 (DTU academic license)
  - alphafold3/: NVIDIA GPU stack + Docker + AF3 model weights
EOF

echo "[12/12] Summary checks..."
[[ -d "$ROOT_DIR/EVmutation/EVmutation" ]]     && echo "  OK EVmutation clone"
[[ -d "$ROOT_DIR/EVmutation/plmc" ]]           && echo "  OK plmc clone"
[[ -d "$ROOT_DIR/rare_codon/cg_cotrans" ]]     && echo "  OK cg_cotrans"
[[ -d "$ROOT_DIR/NetSurfP3/nsp3" ]]            && echo "  OK NetSurfP3 clone"
[[ -d "$ROOT_DIR/netNglyc/signalp-6.0" ]]      && echo "  OK SignalP 6.0 clone"
command -v miranda >/dev/null 2>&1             && echo "  OK miranda on PATH"
[[ -x "$ROOT_DIR/genesplicer/GeneSplicer/genesplicer" ]] && echo "  OK GeneSplicer built"
[[ -d "$ROOT_DIR/alphafold3/alphafold3" ]]     && echo "  OK AlphaFold3 clone"

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