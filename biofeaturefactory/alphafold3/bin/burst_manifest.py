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

"""
Burst-mode manifest reader/writer.

The manifest is the contract between the burst submit phase (which generates
AF3 input JSONs and emits a SLURM job array) and the burst ingest phase
(which reads completed AF3 outputs from the L1 cache and computes deltas).

Schema (TSV with a comment header line + a column header line):

    # burst_manifest_version=1
    array_idx<TAB>input_hash<TAB>pkey<TAB>rbp_name<TAB>allele<TAB>window_idx<TAB>input_json_path<TAB>cache_dir

A row is considered "complete" (cache hit) when its cache_dir contains all
three primary AF3 output suffixes: _model.cif, _confidences.json,
_summary_confidences.json. The same hit criterion is used by af3_runner's
_check_cache and by the slurm_array.sh.tmpl idempotent skip block.
"""

import os
from dataclasses import dataclass, fields
from pathlib import Path
from typing import List


MANIFEST_VERSION = 1
COMMENT_HEADER = f"# burst_manifest_version={MANIFEST_VERSION}"

REQUIRED_SUFFIXES = (
    "_model.cif",
    "_confidences.json",
    "_summary_confidences.json",
)


@dataclass
class ManifestRow:
    """A single AF3 invocation in the burst manifest.

    Each row corresponds to one (gene, mutation, RBP, allele, window) AF3 job.
    Two rows per (gene, mutation, RBP) at single-window — one WT, one MUT.
    """
    array_idx: int
    input_hash: str
    pkey: str            # gene-mutation, e.g. 'F9-C123T'
    rbp_name: str
    allele: str          # 'WT' or 'MUT'
    window_idx: int      # 0 if single-window
    input_json_path: str
    cache_dir: str

    @classmethod
    def column_names(cls) -> List[str]:
        return [f.name for f in fields(cls)]

    def to_tsv_line(self) -> str:
        return "\t".join(str(getattr(self, name)) for name in self.column_names())

    @classmethod
    def from_tsv_fields(cls, parts: List[str]) -> "ManifestRow":
        if len(parts) != len(cls.column_names()):
            raise ValueError(
                f"Manifest row has {len(parts)} fields, expected {len(cls.column_names())}"
            )
        return cls(
            array_idx=int(parts[0]),
            input_hash=parts[1],
            pkey=parts[2],
            rbp_name=parts[3],
            allele=parts[4],
            window_idx=int(parts[5]),
            input_json_path=parts[6],
            cache_dir=parts[7],
        )


def write_manifest(rows: List[ManifestRow], path: Path) -> None:
    """Atomically write the manifest to ``path``.

    Writes to ``path.with_suffix(path.suffix + ".tmp")`` first, fsyncs, then
    renames. An ingest reader running concurrently sees either the previous
    manifest or the new one — never a partial file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")

    header = "\t".join(ManifestRow.column_names())
    with open(tmp, "w") as f:
        f.write(COMMENT_HEADER + "\n")
        f.write(header + "\n")
        for row in rows:
            f.write(row.to_tsv_line() + "\n")
        f.flush()
        os.fsync(f.fileno())

    os.replace(tmp, path)


def read_manifest(path: Path) -> List[ManifestRow]:
    """Read and validate the manifest at ``path``.

    Raises ValueError if the version line is absent or mismatched.
    """
    path = Path(path)
    rows: List[ManifestRow] = []
    expected_cols = ManifestRow.column_names()

    with open(path, "r") as f:
        first = f.readline().rstrip("\n")
        if not first.startswith("# burst_manifest_version="):
            raise ValueError(f"{path}: missing version header line")
        try:
            version = int(first.split("=", 1)[1].strip())
        except (IndexError, ValueError):
            raise ValueError(f"{path}: malformed version header: {first!r}")
        if version != MANIFEST_VERSION:
            raise ValueError(
                f"{path}: manifest version {version}, this code expects {MANIFEST_VERSION}"
            )

        header = f.readline().rstrip("\n").split("\t")
        if header != expected_cols:
            raise ValueError(
                f"{path}: manifest column header {header} does not match expected {expected_cols}"
            )

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(ManifestRow.from_tsv_fields(line.split("\t")))

    return rows


def is_cache_complete(cache_dir: Path) -> bool:
    """Return True iff cache_dir contains all three required AF3 output files.

    Mirrors the hit criterion in af3_runner._check_cache and slurm_array.sh.tmpl.
    """
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return False
    for suffix in REQUIRED_SUFFIXES:
        if not list(cache_dir.glob(f"**/*{suffix}")):
            return False
    return True


def count_pending(manifest_path: Path, cache_root: Path = None) -> int:
    """Count manifest rows whose cache_dir is not yet complete.

    ``cache_root`` is unused — the manifest already records absolute cache_dir
    paths, so the lookup is direct. Kept as a parameter for API forward-compat
    in case the manifest format moves to relative paths.
    """
    rows = read_manifest(manifest_path)
    return sum(1 for r in rows if not is_cache_complete(Path(r.cache_dir)))
