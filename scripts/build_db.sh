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

# End-to-end database build for codon-aware and AF3 pipelines.
#
# Default output root: <repo>/Bio_DBs
#
# Outputs:
#   - assembly_summary_refseq.txt
#   - ftp_paths.txt
#   - refseq_assemblies/<ASM>/{*_protein.faa.gz,*_cds_from_genomic.fna.gz,*_feature_table.txt.gz}
#   - refseq_proteins_merged.faa
#   - idmapping.dat.gz
#   - protein_id_to_refseq.tsv
#   - uniref90.fasta.gz (UniRef90 for jackhmmer MSA pipeline)
#   - mature_hsa.fasta (latest miRBase CURRENT)
#   - AF3/RBP_db/* (POSTAR3 indexed if present)
#
# Optional env vars:
#   DB_ROOT=<path>                         (default: <repo>/Bio_DBs)
#   TAXON_GROUP=vertebrate_mammalian      (RefSeq "group" column filter)
#   ARIA_SPLIT=4
#   ARIA_CONN=4
#   PARALLEL_JOBS=8                       (assemblies downloaded concurrently)
#   UNIREF90_URL=...                      (default: UniProt FTP uniref90.fasta.gz)
#   MIRBASE_HSA_URL=...                   (default: http://www.benoslab.pitt.edu/comir/hsa_data/mature_hsa.fa)
##   MIRBASE_REFRESH=1                    (1=always refresh mature.fa from CURRENT, 0=only if missing)
#   POSTAR3_TXT_URL=...                   (optional; download if human-POSTAR3.txt missing)
#   AF3_RBP_MSA_ARCHIVE_URL=...           (optional; tar(.gz/.zst/.xz/.bz2) with msa/*.a3m)
#   AF3_DOWNLOAD_RBP_MSAS=1               (1=download per-ID MSAs from rbp_uniprot_ids.txt)
#   AF3_MSA_VERSION=v6                    (used in AF3_MSA_URL_TEMPLATE fallback)
#   AF3_MSA_URL_TEMPLATE=...              (default: https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-msa_{VERSION}.a3m)

TAXON_GROUP="${TAXON_GROUP:-vertebrate_mammalian}"
EXTRA_TAXON_GROUPS="${EXTRA_TAXON_GROUPS:-vertebrate_other invertebrate}"
ARIA_SPLIT="${ARIA_SPLIT:-4}"
ARIA_CONN="${ARIA_CONN:-4}"
PARALLEL_JOBS="${PARALLEL_JOBS:-8}"

UNIREF90_URL="${UNIREF90_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz}"

MIRBASE_HSA_URL="${MIRBASE_HSA_URL:-http://www.benoslab.pitt.edu/comir/hsa_data/mature_hsa.fa}"
MIRBASE_REFRESH="${MIRBASE_REFRESH:-1}"
POSTAR3_TXT_URL="${POSTAR3_TXT_URL:-}"
AF3_RBP_MSA_ARCHIVE_URL="${AF3_RBP_MSA_ARCHIVE_URL:-}"
AF3_DOWNLOAD_RBP_MSAS="${AF3_DOWNLOAD_RBP_MSAS:-1}"
AF3_MSA_VERSION="${AF3_MSA_VERSION:-v6}"
AF3_MSA_URL_TEMPLATE="${AF3_MSA_URL_TEMPLATE:-https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-msa_{VERSION}.a3m}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DB_ROOT="${DB_ROOT:-$PROJECT_ROOT/../Bio_DBs}"
AWK_SCRIPT="$SCRIPT_DIR/build_refseq_ftp_paths.awk"

mkdir -p "$DB_ROOT"
cd "$DB_ROOT"

download_file() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "  EXISTS $out"
    return 0
  fi
  aria2c --continue=true \
    --max-connection-per-server="${ARIA_CONN}" \
    --split="${ARIA_SPLIT}" \
    --min-split-size=4M \
    --max-tries=20 \
    --retry-wait=10 \
    --timeout=60 \
    --lowest-speed-limit=50K \
    --file-allocation=none \
    -d "$(dirname "$out")" \
    -o "$(basename "$out")" \
    "$url"
}

echo "[1/10] Downloading RefSeq assembly summary..."
if [[ -s assembly_summary_refseq.txt ]]; then
  echo "  EXISTS assembly_summary_refseq.txt"
else
  if command -v aria2c &>/dev/null; then
    aria2c -x 16 -s 16 -k 1M -c \
      -o assembly_summary_refseq.txt \
      https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
  else
    curl -L -o assembly_summary_refseq.txt \
      https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
  fi
fi

echo "[2/10] Building ftp_paths.txt (group=${TAXON_GROUP})..."
awk -F '\t' -v taxon_group="${TAXON_GROUP}" -f "$AWK_SCRIPT" assembly_summary_refseq.txt > ftp_paths.txt
echo "  ftp paths: $(wc -l < ftp_paths.txt)"

echo "[3/10] Downloading RefSeq assembly files..."
mkdir -p refseq_assemblies logs

pids=()
while IFS= read -r FTP_PATH; do
  # throttle number of concurrent assembly jobs
  while [[ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$PARALLEL_JOBS" ]]; do
    sleep 1
  done

  {
    ASM="$(basename "$FTP_PATH")"
    OUT_DIR="refseq_assemblies/$ASM"
    mkdir -p "$OUT_DIR"

    for suf in protein.faa.gz cds_from_genomic.fna.gz feature_table.txt.gz; do
      URL="${FTP_PATH}/${ASM}_${suf}"
      OUT="${OUT_DIR}/${ASM}_${suf}"
      OUT_UNZIPPED="${OUT%.gz}"

      if [[ -s "$OUT" || -s "$OUT_UNZIPPED" ]]; then
        echo "  EXISTS $OUT or $OUT_UNZIPPED"
        continue
      fi

      code="$(curl -s -o /dev/null -w "%{http_code}" -I "$URL")"
      if [[ "$code" != "200" ]]; then
        echo "  SKIP [$code] $URL" | tee -a logs/download_missing.log
        continue
      fi

      if command -v aria2c &>/dev/null; then
        aria2c --continue=true \
          --max-connection-per-server="${ARIA_CONN}" \
          --split="${ARIA_SPLIT}" \
          --min-split-size=4M \
          --max-tries=20 \
          --retry-wait=10 \
          --timeout=60 \
          --lowest-speed-limit=50K \
          --file-allocation=none \
          -d "$OUT_DIR" \
          -o "$(basename "$OUT")" \
          "$URL" || echo "FAIL $URL" >> logs/download_fail.log
      else
        curl -L -o "$OUT" "$URL" || echo "FAIL $URL" >> logs/download_fail.log
      fi
    done
  } &

  pids+=("$!")
done < ftp_paths.txt

download_fail=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    download_fail=1
  fi
done
if [[ "$download_fail" -ne 0 ]]; then
  echo "ERROR: One or more parallel download workers failed." >&2
  exit 1
fi

if [[ -n "$EXTRA_TAXON_GROUPS" ]]; then
  echo "[3b/9] Downloading extra RefSeq assembly groups: ${EXTRA_TAXON_GROUPS}..."
  for extra_grp in $EXTRA_TAXON_GROUPS; do
    echo "  Group: ${extra_grp}"
    extra_ftp_paths="$(mktemp)"
    awk -F '\t' -v taxon_group="${extra_grp}" -f "$AWK_SCRIPT" assembly_summary_refseq.txt > "$extra_ftp_paths"
    echo "  ftp paths for ${extra_grp}: $(wc -l < "$extra_ftp_paths")"

    pids=()
    while IFS= read -r FTP_PATH; do
      while [[ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$PARALLEL_JOBS" ]]; do
        sleep 1
      done

      {
        ASM="$(basename "$FTP_PATH")"
        OUT_DIR="refseq_assemblies/$ASM"
        mkdir -p "$OUT_DIR"

        for suf in protein.faa.gz cds_from_genomic.fna.gz feature_table.txt.gz; do
          URL="${FTP_PATH}/${ASM}_${suf}"
          OUT="${OUT_DIR}/${ASM}_${suf}"
          OUT_UNZIPPED="${OUT%.gz}"

          if [[ -s "$OUT" || -s "$OUT_UNZIPPED" ]]; then
            echo "  EXISTS $OUT or $OUT_UNZIPPED"
            continue
          fi

          code="$(curl -s -o /dev/null -w "%{http_code}" -I "$URL")"
          if [[ "$code" != "200" ]]; then
            echo "  SKIP [$code] $URL" | tee -a logs/download_missing.log
            continue
          fi

          if command -v aria2c &>/dev/null; then
            aria2c --continue=true \
              --max-connection-per-server="${ARIA_CONN}" \
              --split="${ARIA_SPLIT}" \
              --min-split-size=4M \
              --max-tries=20 \
              --retry-wait=10 \
              --timeout=60 \
              --lowest-speed-limit=50K \
              --file-allocation=none \
              -d "$OUT_DIR" \
              -o "$(basename "$OUT")" \
              "$URL" || echo "FAIL $URL" >> logs/download_fail.log
          else
            curl -L -o "$OUT" "$URL" || echo "FAIL $URL" >> logs/download_fail.log
          fi
        done
      } &

      pids+=("$!")
    done < "$extra_ftp_paths"

    download_fail=0
    for pid in "${pids[@]}"; do
      if ! wait "$pid"; then
        download_fail=1
      fi
    done
    if [[ "$download_fail" -ne 0 ]]; then
      echo "ERROR: One or more parallel download workers failed for group ${extra_grp}." >&2
      exit 1
    fi

    rm -f "$extra_ftp_paths"
  done

  # Force-rebuild merged FAA so extra groups are included.
  rm -f refseq_proteins_merged.faa
  echo "  Removed cached merged FAA; will rebuild in step [4/10] to include extra groups"
fi

echo "[4/10] Building merged RefSeq protein FASTA..."
if [[ -s refseq_proteins_merged.faa ]]; then
  echo "  EXISTS refseq_proteins_merged.faa"
else
  : > refseq_proteins_merged.faa
  while IFS= read -r f; do
    if [[ "$f" == *.gz ]]; then
      gzip -dc "$f" >> refseq_proteins_merged.faa
    else
      cat "$f" >> refseq_proteins_merged.faa
    fi
  done < <(find refseq_assemblies -type f \( -name '*_protein.faa.gz' -o -name '*_protein.faa' \) | sort)
fi

echo "[5/10] Downloading UniProt idmapping.dat.gz (if missing)..."
if [[ ! -s idmapping.dat.gz && ! -s idmapping.dat ]]; then
  if command -v aria2c &>/dev/null; then
    aria2c -x 16 -s 16 -k 1M -c \
      -o idmapping.dat.gz \
      https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
  else
    curl -L -o idmapping.dat.gz \
      https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
  fi
else
  echo "  EXISTS idmapping.dat.gz or idmapping.dat"
fi

echo "[6/10] Building protein_id_to_refseq.tsv..."
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

if [[ -s idmapping.dat.gz ]]; then
  IDMAP_READ_CMD=(gzip -dc idmapping.dat.gz)
elif [[ -s idmapping.dat ]]; then
  IDMAP_READ_CMD=(cat idmapping.dat)
else
  echo "ERROR: Neither idmapping.dat.gz nor idmapping.dat is available." >&2
  exit 1
fi

# UniProt -> RefSeq
"${IDMAP_READ_CMD[@]}" \
  | awk -F $'\t' '$2=="RefSeq"{print $1"\t"$3}' \
  | sort -u > "$tmpdir/uniprot_to_refseq.tsv"

# UniProt -> UniParc
"${IDMAP_READ_CMD[@]}" \
  | awk -F $'\t' '$2=="UniParc"{print $1"\t"$3}' \
  | sort -u > "$tmpdir/uniprot_to_uniparc.tsv"

# UniParc -> UniProt
awk -F $'\t' '{print $2"\t"$1}' "$tmpdir/uniprot_to_uniparc.tsv" \
  | sort -u > "$tmpdir/uniparc_to_uniprot.tsv"

# UniParc -> RefSeq via UniProt join
join -t $'\t' -1 2 -2 1 \
  <(sort -t $'\t' -k2,2 "$tmpdir/uniparc_to_uniprot.tsv") \
  <(sort -t $'\t' -k1,1 "$tmpdir/uniprot_to_refseq.tsv") \
  | awk -F $'\t' '{print $2"\t"$3}' \
  | sort -u > "$tmpdir/uniparc_to_refseq.tsv"

# Combined lookup: (UniProt OR UniParc) -> RefSeq
cat "$tmpdir/uniprot_to_refseq.tsv" "$tmpdir/uniparc_to_refseq.tsv" \
  | sort -u > protein_id_to_refseq.tsv

echo "[7/10] Downloading UniRef90 (jackhmmer MSA pipeline)..."
if [[ -s uniref90.fasta.gz || -s uniref90.fasta ]]; then
  echo "  EXISTS uniref90.fasta.gz or uniref90.fasta"
else
  if command -v aria2c &>/dev/null; then
    aria2c -x 16 -s 16 -k 4M -c \
      --max-tries=20 \
      --retry-wait=10 \
      --timeout=120 \
      --lowest-speed-limit=50K \
      --file-allocation=none \
      -o uniref90.fasta.gz \
      "$UNIREF90_URL"
  else
    curl -L -o uniref90.fasta.gz "$UNIREF90_URL"
  fi
fi

echo "[8/10] Downloading human mature miRNA sequences (mature_hsa.fasta)..."
if [[ "$MIRBASE_REFRESH" -eq 1 ]] || [[ ! -s mature_hsa.fasta ]]; then
  if command -v aria2c &>/dev/null; then
    aria2c --continue=true \
      --allow-overwrite=true \
      --max-connection-per-server="${ARIA_CONN}" \
      --split="${ARIA_SPLIT}" \
      --min-split-size=4M \
      --max-tries=20 \
      --retry-wait=10 \
      --timeout=60 \
      --lowest-speed-limit=50K \
      --file-allocation=none \
      -d "." \
      -o "mature_hsa.fasta" \
      "$MIRBASE_HSA_URL"
  else
    curl -L -o mature_hsa.fasta "$MIRBASE_HSA_URL"
  fi
else
  echo "  EXISTS mature_hsa.fasta"
fi
_hsa_count="$(grep -c '^>' mature_hsa.fasta 2>/dev/null || echo 0)"
echo "  WROTE mature_hsa.fasta (${_hsa_count} human mature miRNAs)"

echo "[9/10] AF3 RBP DB setup (POSTAR3 + tabix indexes)..."
RBP_DB_DIR="AF3/RBP_db"
mkdir -p "$RBP_DB_DIR"

if [[ ! -s "$RBP_DB_DIR/human-POSTAR3.txt" && -n "$POSTAR3_TXT_URL" ]]; then
  download_file "$POSTAR3_TXT_URL" "$RBP_DB_DIR/human-POSTAR3.txt"
fi

if [[ ! -d "$RBP_DB_DIR/msa" && -n "$AF3_RBP_MSA_ARCHIVE_URL" ]]; then
  archive_ext="${AF3_RBP_MSA_ARCHIVE_URL##*.}"
  case "$archive_ext" in
    zst) out_archive="$RBP_DB_DIR/rbp_msa.tar.zst" ;;
    xz) out_archive="$RBP_DB_DIR/rbp_msa.tar.xz" ;;
    bz2) out_archive="$RBP_DB_DIR/rbp_msa.tar.bz2" ;;
    gz) out_archive="$RBP_DB_DIR/rbp_msa.tar.gz" ;;
    *) out_archive="$RBP_DB_DIR/rbp_msa.tar" ;;
  esac
  download_file "$AF3_RBP_MSA_ARCHIVE_URL" "$out_archive"
  mkdir -p "$RBP_DB_DIR/msa"
  case "$out_archive" in
    *.tar.zst) tar --zstd -xf "$out_archive" -C "$RBP_DB_DIR" ;;
    *.tar.xz) tar -xJf "$out_archive" -C "$RBP_DB_DIR" ;;
    *.tar.bz2) tar -xjf "$out_archive" -C "$RBP_DB_DIR" ;;
    *.tar.gz) tar -xzf "$out_archive" -C "$RBP_DB_DIR" ;;
    *.tar) tar -xf "$out_archive" -C "$RBP_DB_DIR" ;;
  esac
fi

if [[ "$AF3_DOWNLOAD_RBP_MSAS" -eq 1 && -s "$RBP_DB_DIR/rbp_uniprot_ids.txt" ]]; then
  mkdir -p "$RBP_DB_DIR/msa"
  echo "  Downloading per-RBP AF MSAs from rbp_uniprot_ids.txt..."
  tmp_ids="$(mktemp)"
  # Accept first token/column as ID, skip comments/empty lines.
  awk '
    /^[[:space:]]*#/ {next}
    NF==0 {next}
    {
      gsub(/\r/, "", $1)
      print $1
    }
  ' "$RBP_DB_DIR/rbp_uniprot_ids.txt" | sort -u > "$tmp_ids"

  while IFS= read -r uniprot_id; do
    [[ -z "$uniprot_id" ]] && continue

    out="$RBP_DB_DIR/msa/AF-${uniprot_id}-F1-msa_${AF3_MSA_VERSION}.a3m"
    if [[ -s "$out" ]]; then
      continue
    fi

    # Try configured template first, then common historical versions as fallback.
    try_urls=()
    try_urls+=("${AF3_MSA_URL_TEMPLATE//\{ID\}/$uniprot_id}")
    try_urls[0]="${try_urls[0]//\{VERSION\}/$AF3_MSA_VERSION}"
    try_urls+=("https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-msa_v6.a3m")
    try_urls+=("https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-msa_v5.a3m")
    try_urls+=("https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-msa_v4.a3m")

    downloaded=0
    for url in "${try_urls[@]}"; do
      code="$(curl -s -o /dev/null -w "%{http_code}" -I "$url")"
      if [[ "$code" == "200" ]]; then
        download_file "$url" "$out"
        downloaded=1
        break
      fi
    done

    if [[ "$downloaded" -ne 1 ]]; then
      echo "  WARN missing MSA for ${uniprot_id}" >> logs/af3_msa_missing.log
    fi
  done < "$tmp_ids"
  rm -f "$tmp_ids"
fi

if [[ -s "$RBP_DB_DIR/human-POSTAR3.txt" ]]; then
  if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    if [[ ! -s "$RBP_DB_DIR/human-POSTAR3.bed.gz" ]]; then
      LC_ALL=C sort -k1,1V -k2,2n "$RBP_DB_DIR/human-POSTAR3.txt" \
        | bgzip > "$RBP_DB_DIR/human-POSTAR3.bed.gz"
    fi
    if [[ ! -s "$RBP_DB_DIR/human-POSTAR3.bed.gz.tbi" ]]; then
      tabix -f -p bed "$RBP_DB_DIR/human-POSTAR3.bed.gz"
    fi
  else
    echo "  WARN bgzip/tabix not found; cannot build human-POSTAR3.bed.gz(.tbi)"
  fi
fi

if [[ -d "$RBP_DB_DIR/msa" ]]; then
  msa_count="$(find "$RBP_DB_DIR/msa" -type f -name '*.a3m' | wc -l | tr -d ' ')"
  echo "  RBP MSA files: $msa_count"
else
  echo "  WARN AF3/RBP_db/msa is missing. Set AF3_RBP_MSA_ARCHIVE_URL to auto-populate."
fi

echo "[10/10] Summary"
echo "Done."
echo "  DB root: ${DB_ROOT}"
echo "  RefSeq assemblies dir: ${DB_ROOT}/refseq_assemblies"
echo "  Merged proteins: ${DB_ROOT}/refseq_proteins_merged.faa"
echo "  ID map: ${DB_ROOT}/protein_id_to_refseq.tsv"
echo "  UniRef90: ${DB_ROOT}/uniref90.fasta.gz"
echo "  miRNA FASTA: ${DB_ROOT}/mature_hsa.fasta"
echo "  AF3 RBP DB dir: ${DB_ROOT}/AF3/RBP_db"
