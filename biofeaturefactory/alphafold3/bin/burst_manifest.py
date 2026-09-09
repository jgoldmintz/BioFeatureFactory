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
AF3 input JSONs and emits a SLURM job array) and each array worker (which reads
one immutable input row and publishes completed output to the L1 cache).

Schema (TSV with a comment header line + a column header line):

    # burst_manifest_version=2
    array_idx<TAB>input_hash<TAB>pkey<TAB>rbp_name<TAB>allele<TAB>window_idx<TAB>input_json_path<TAB>cache_dir<TAB>output_name

A row is complete only when ``cache_dir/output_name`` contains the exact,
nonempty top-ranked AF3 model/confidences/summary sibling triplet and both JSON
files parse. The SLURM script calls this module's CLI so submit, worker, and
ingest share this predicate and the same atomic publisher.
"""

import argparse
import errno
import fcntl
import json
import os
import shutil
import tempfile
import sys
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass, fields
from pathlib import Path
from typing import List, Optional, Tuple


MANIFEST_VERSION = 2
COMMENT_HEADER = f"# burst_manifest_version={MANIFEST_VERSION}"
_CACHE_LOCK_NAME = ".bff-af3-cache.lock"
_STAGE_TICKET = '.bff-af3-stage.json'


def cache_coordination_dir(cache_root: Path) -> Path:
    """Stable sibling of the AF3 cache; clear_cache never deletes this tree."""
    cache_root = Path(cache_root).resolve()
    return cache_root.with_name(f'.{cache_root.name}-coordination')


@contextmanager
def cache_lock(cache_root: Path):
    coordination = cache_coordination_dir(cache_root)
    coordination.mkdir(parents=True, exist_ok=True)
    with open(coordination / _CACHE_LOCK_NAME, 'a+') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        yield coordination


def _read_generation(coordination: Path) -> str:
    path = coordination / 'generation'
    return path.read_text().strip() if path.exists() else '0'


def cache_generation(cache_root: Path) -> str:
    with cache_lock(cache_root) as coordination:
        return _read_generation(coordination)


def clear_cache(cache_root: Path) -> None:
    """Invalidate old workers and clear completed AF3 entries, never staging."""
    cache_root = Path(cache_root)
    with cache_lock(cache_root) as coordination:
        generation_path = coordination / 'generation'
        temporary = generation_path.with_suffix('.tmp')
        temporary.write_text(uuid.uuid4().hex)
        os.replace(temporary, generation_path)
        if cache_root.exists():
            shutil.rmtree(cache_root)
        cache_root.mkdir(parents=True, exist_ok=True)


def begin_cache_stage(cache_dir: Path, expected_generation: Optional[str] = None) -> Path:
    """Allocate an output tree protected from clearing, with a generation ticket."""
    cache_dir = Path(cache_dir).resolve()
    with cache_lock(cache_dir.parent) as coordination:
        generation = _read_generation(coordination)
        if expected_generation is not None and expected_generation != generation:
            raise ValueError('AF3 submission cache generation was cleared; resubmit it')
        stage = Path(tempfile.mkdtemp(prefix=f'.af3_{cache_dir.name}.', dir=coordination))
        (stage / _STAGE_TICKET).write_text(json.dumps({
            'cache_dir': str(cache_dir), 'generation': generation,
        }))
        return stage


@dataclass
class ManifestRow:
    """A single AF3 invocation in the burst manifest.

    Each row corresponds to one (gene, mutation, RBP, allele, window) AF3 job.
    Two rows per (gene, mutation, RBP) at single-window -- one WT, one MUT.
    """
    array_idx: int
    input_hash: str
    pkey: str            # gene-mutation, e.g. 'F9-C123T'
    rbp_name: str
    allele: str          # 'WT' or 'MUT'
    window_idx: int      # 0 if single-window
    input_json_path: str
    cache_dir: str
    output_name: str

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
            output_name=parts[8],
        )


def write_manifest(rows: List[ManifestRow], path: Path) -> None:
    """Atomically write the manifest to ``path``.

    Writes to ``path.with_suffix(path.suffix + ".tmp")`` first, fsyncs, then
    renames. A worker sees either the complete manifest or no manifest -- never
    a partial file.
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


def _result_paths(result_dir: Path, output_name: str) -> Tuple[Path, Path, Path]:
    """Return the exact primary AF3 paths for one cache entry."""
    if (
        not output_name
        or output_name in {".", ".."}
        or "/" in output_name
        or "\\" in output_name
    ):
        raise ValueError(f"Invalid AF3 output name: {output_name!r}")

    result_dir = Path(result_dir)
    return (
        result_dir / f"{output_name}_model.cif",
        result_dir / f"{output_name}_confidences.json",
        result_dir / f"{output_name}_summary_confidences.json",
    )


def is_result_complete(result_dir: Path, output_name: str) -> bool:
    """Return whether one exact AF3 result has valid primary outputs."""
    try:
        model_path, confidences_path, summary_path = _result_paths(
            result_dir, output_name
        )
    except ValueError:
        return False

    for path in (model_path, confidences_path, summary_path):
        try:
            if not path.is_file() or path.stat().st_size == 0:
                return False
        except OSError:
            return False

    try:
        with open(confidences_path, "r") as confidences_handle:
            json.load(confidences_handle)
        with open(summary_path, "r") as summary_handle:
            json.load(summary_handle)
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return True


def is_cache_complete(cache_dir: Path, output_name: str) -> bool:
    """Return whether a content-addressed cache entry is complete."""
    return is_result_complete(Path(cache_dir) / output_name, output_name)


def _quarantine_path(cache_dir: Path) -> Path:
    """Return an unused same-parent path for an existing cache entry."""
    while True:
        suffix = (
            f"stale-{time.strftime('%Y%m%d-%H%M%S')}-"
            f"{uuid.uuid4().hex[:8]}"
        )
        quarantined = cache_dir.with_name(f"{cache_dir.name}.{suffix}")
        if not os.path.lexists(quarantined):
            return quarantined


def publish_cache(
    source_dir: Path,
    cache_dir: Path,
    output_name: str,
) -> Path:
    """Validate and atomically publish one AF3 output root into the cache."""
    source_dir = Path(source_dir)
    cache_dir = Path(cache_dir)
    _result_paths(source_dir / output_name, output_name)

    if source_dir.is_symlink() or not source_dir.is_dir():
        raise ValueError(f"AF3 cache source is not a directory: {source_dir}")

    source_absolute = Path(os.path.abspath(source_dir))
    cache_absolute = Path(os.path.abspath(cache_dir))
    if (
        source_absolute == cache_absolute
        or cache_absolute in source_absolute.parents
        or source_absolute in cache_absolute.parents
    ):
        raise ValueError("AF3 cache source and destination must not overlap")

    cache_root = cache_dir.parent
    with cache_lock(cache_root) as coordination:
        if not is_cache_complete(source_dir, output_name):
            raise ValueError(
                f"AF3 cache source is incomplete for {output_name}: {source_dir}"
            )
        try:
            ticket = json.loads((source_dir / _STAGE_TICKET).read_text())
        except (OSError, ValueError) as error:
            raise ValueError('AF3 source lacks a valid stage generation ticket') from error
        if ticket != {'cache_dir': str(cache_dir.resolve()),
                      'generation': _read_generation(coordination)}:
            raise ValueError('AF3 staging generation was cleared; output preserved, not published')
        cache_root.mkdir(parents=True, exist_ok=True)
        if source_dir.stat().st_dev != cache_root.stat().st_dev:
            raise OSError(
                errno.EXDEV,
                "AF3 cache publication requires source and destination on "
                "the same filesystem",
                str(source_dir),
            )

        quarantined: Optional[Path] = None
        if os.path.lexists(cache_dir):
            quarantined = _quarantine_path(cache_dir)
            os.replace(cache_dir, quarantined)

        try:
            os.replace(source_dir, cache_dir)
        except BaseException:
            if quarantined is not None and not os.path.lexists(cache_dir):
                os.replace(quarantined, cache_dir)
            raise

    return cache_dir


def count_pending(manifest_path: Path, cache_root: Path = None) -> int:
    """Count manifest rows whose cache_dir is not yet complete.

    ``cache_root`` is unused -- the manifest already records absolute cache_dir
    paths, so the lookup is direct. Kept as a parameter for API forward-compat
    in case the manifest format moves to relative paths.
    """
    rows = read_manifest(manifest_path)
    return sum(
        1
        for row in rows
        if not is_cache_complete(Path(row.cache_dir), row.output_name)
    )


def main(argv: Optional[List[str]] = None) -> int:
    """Run the cache predicate or publisher for a SLURM array task."""
    parser = argparse.ArgumentParser(description="AF3 burst cache helper")
    subparsers = parser.add_subparsers(dest="command", required=True)

    complete_parser = subparsers.add_parser(
        "check", help="exit 0 only for a complete exact cache result"
    )
    complete_parser.add_argument("cache_dir")
    complete_parser.add_argument("output_name")

    publish_parser = subparsers.add_parser(
        "publish", help="validate and atomically publish an AF3 output root"
    )
    publish_parser.add_argument("source_dir")
    publish_parser.add_argument("cache_dir")
    publish_parser.add_argument("output_name")
    stage_parser = subparsers.add_parser('stage', help='allocate generation-protected staging')
    stage_parser.add_argument('cache_dir')
    stage_parser.add_argument('generation')

    args = parser.parse_args(argv)
    if args.command == "check":
        return 0 if is_cache_complete(args.cache_dir, args.output_name) else 1

    try:
        if args.command == 'stage':
            print(begin_cache_stage(args.cache_dir, args.generation))
            return 0
        publish_cache(args.source_dir, args.cache_dir, args.output_name)
    except (OSError, ValueError) as error:
        print(f"ERR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
