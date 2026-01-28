#!/usr/bin/env python3
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
Controller for the SpliceAI Nextflow pipeline.

- Mirrors the pipeline's CLI (mutations, reference, mapping paths, etc.).
- Launches `nextflow run main.nf ...` and streams its output.
- Runs the spliceai-tracker automatically unless --disable_tracker is set.
- Use --resume to continue from a previous run (Nextflow's built-in resume).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

HERE = Path(__file__).resolve().parent
NEXTFLOW_SCRIPT = HERE / "bin/main.nf"
TRACKER_SCRIPT = HERE / "bin/spliceai-tracker.py"
WORK_ROOT = HERE / "work"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Controller for the SpliceAI Nextflow pipeline."
    )
    # Pipeline inputs / options
    parser.add_argument("--mutations_path", type=Path,
                        help="Dir or file with <GENE>_mutations.csv (required unless using --input_vcf_path)")
    parser.add_argument("--input_vcf_path", type=Path,
                        help="Existing VCF (file) or directory of VCFs when skipping VCF generation")
    parser.add_argument("--skip_vcf_generation", action="store_true",
                        help="Use pre-built VCFs (requires --input_vcf_path)")
    parser.add_argument("--reference_genome", type=Path, required=True)
    parser.add_argument("--annotation_file", type=Path, required=True)
    parser.add_argument("--transcript_mapping_path", type=Path, required=True)
    parser.add_argument("--chromosome_mapping_path", type=Path, required=True)
    parser.add_argument("--genomic_mapping_path", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, default=Path("."))
    parser.add_argument("--vcf_output_dir", type=Path)
    parser.add_argument("--validation_log", type=Path)
    parser.add_argument("--chromosome_format", choices=["refseq", "simple", "ucsc"], default="refseq")
    parser.add_argument("--splice_threshold", type=float, default=0.0)
    parser.add_argument("--retry_jitter", type=int, default=10)
    parser.add_argument("--clear_vcf_cache", action="store_true")
    parser.add_argument("--validate_mapping", action="store_true")
    parser.add_argument("--maxforks", type=int, default=0,
                        help="Max concurrent run_spliceai tasks (0=Nextflow default)")
    parser.add_argument("--forceAll_isoforms", action="store_true",
                        help="Process all isoforms regardless of count (default: apply filtering for genes >50 isoforms)")
    parser.add_argument("--max_isoforms_per_gene", type=int, default=50,
                        help="Maximum isoforms per gene before hybrid filtering is applied (default: 50)")

    # Controller options
    parser.add_argument("--poll_seconds", type=float, default=15.0,
                        help="Tracker poll interval")
    parser.add_argument("--disable_tracker", action="store_true",
                        help="Skip launching spliceai-tracker.py")
    parser.add_argument("--pipeline", type=Path, default=NEXTFLOW_SCRIPT,
                        help="Nextflow script to run (default: main.nf)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from a previous run using Nextflow's -resume flag")

    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args: argparse.Namespace) -> None:
    if args.skip_vcf_generation and not args.input_vcf_path:
        raise SystemExit("ERROR: --skip_vcf_generation requires --input_vcf_path")
    if args.input_vcf_path and not args.input_vcf_path.exists():
        raise SystemExit(f"ERROR: --input_vcf_path {args.input_vcf_path} not found")
    if args.input_vcf_path and not args.skip_vcf_generation:
        raise SystemExit("ERROR: When providing --input_vcf_path you must also pass --skip_vcf_generation")
    if not args.input_vcf_path and not args.mutations_path:
        raise SystemExit("ERROR: Must provide --mutations_path when not using --input_vcf_path")


def normalize(path: Optional[Path]) -> Optional[str]:
    return str(path.resolve()) if path else None


def build_nextflow_cmd(args: argparse.Namespace) -> List[str]:
    cmd = ["nextflow", "run", str(args.pipeline)]
    if args.resume:
        cmd.append("-resume")

    def add_param(name: str, value):
        if value is not None:
            cmd.extend([f"--{name}", str(value)])

    def add_flag(name: str, enabled: bool):
        if enabled:
            cmd.append(f"--{name}")

    add_param("mutations_path", args.mutations_path)
    if args.input_vcf_path:
        if args.input_vcf_path.is_file():
            add_param("input_vcf_file", args.input_vcf_path)
        else:
            add_param("input_vcf_dir", args.input_vcf_path)

    add_flag("skip_vcf_generation", args.skip_vcf_generation)
    add_param("reference_genome", args.reference_genome)
    add_param("annotation_file", args.annotation_file)
    add_param("transcript_mapping_path", args.transcript_mapping_path)
    add_param("chromosome_mapping_path", args.chromosome_mapping_path)
    add_param("genomic_mapping_path", args.genomic_mapping_path)
    add_param("output_dir", args.output_dir)
    add_param("vcf_output_dir", args.vcf_output_dir)
    add_param("validation_log", args.validation_log)
    add_param("chromosome_format", args.chromosome_format)
    add_param("splice_threshold", args.splice_threshold)
    add_param("retry_jitter", args.retry_jitter)
    add_flag("clear_vcf_cache", args.clear_vcf_cache)
    add_flag("validate_mapping", args.validate_mapping)
    if args.maxforks:
        add_param("maxforks", args.maxforks)
    add_flag("forceAll_isoforms", args.forceAll_isoforms)
    add_param("max_isoforms_per_gene", args.max_isoforms_per_gene)
    return cmd


class TrackerRunner:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.proc: Optional[subprocess.Popen] = None

    def start(self):
        if self.args.disable_tracker:
            return
        mutations = normalize(self.args.mutations_path)
        if not mutations:
            return
        cmd = [
            sys.executable,
            str(TRACKER_SCRIPT),
            "--mutations-path",
            mutations,
            "--work-root",
            str(WORK_ROOT),
            "--poll-seconds",
            str(max(self.args.poll_seconds, 2.0)),
        ]
        self.proc = subprocess.Popen(cmd, cwd=str(HERE))

    def stop(self):
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.terminate()
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.proc.kill()
        self.proc = None


def run_controller(args: argparse.Namespace):
    tracker = TrackerRunner(args)
    cmd = build_nextflow_cmd(args)
    print(f"[controller] Launching: {' '.join(cmd)}", flush=True)
    nf_proc = subprocess.Popen(cmd, cwd=str(HERE))
    tracker.start()
    try:
        exit_code = nf_proc.wait()
    finally:
        tracker.stop()
    sys.exit(exit_code)


if __name__ == "__main__":
    run_controller(parse_args())
