# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# Licensed under the PolyForm Noncommercial License 1.0.0.
# You may use this software for any purpose other than commercial use.
# See LICENSE, or https://polyformproject.org/licenses/noncommercial/1.0.0/
#
# Required Notice: Copyright (C) 2023-2026 Jacob Goldmintz
# (https://github.com/jgoldmintz)

# Usage:
#   awk -F '\t' -f build_refseq_ftp_paths.awk assembly_summary_refseq.txt > ftp_paths.txt
#
# Optional variables:
#   -v taxon_group="vertebrate_mammalian"
#   -v taxon_group="vertebrate_mammalian|vertebrate_other|invertebrate"
#
# taxon_group accepts a single group name or a |-separated list of group names.
# This script is header-aware and resilient to column order changes.

BEGIN {
  OFS = "\t"
  if (taxon_group == "") {
    taxon_group = "vertebrate_mammalian"
  }
}

# Parse the main header row with column names.
/^#.*assembly_accession/ {
  for (i = 1; i <= NF; i++) {
    h = $i
    gsub(/^# */, "", h)
    col[h] = i
  }
  next
}

# Skip other comment lines.
/^#/ { next }

{
  # Require expected columns; if missing, skip safely.
  if (!("version_status" in col) || !("refseq_category" in col) || !("assembly_level" in col) || !("ftp_path" in col) || !("group" in col)) {
    next
  }

  idx_vs  = col["version_status"]
  idx_cat = col["refseq_category"]
  idx_lvl = col["assembly_level"]
  idx_ftp = col["ftp_path"]
  idx_grp = col["group"]

  vs = $idx_vs
  category = $idx_cat
  lvl = $idx_lvl
  ftp = $idx_ftp
  grp = $idx_grp

  n = split(taxon_group, groups, "|")
  matched = 0
  for (k = 1; k <= n; k++) {
    if (grp == groups[k]) { matched = 1; break }
  }
  if (vs == "latest" && (category == "reference genome" || category == "representative genome") && (lvl == "Complete Genome" || lvl == "Chromosome" || lvl == "Scaffold") && ftp != "" && ftp != "na" && matched) {
    print ftp
  }
}
