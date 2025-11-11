#!/usr/bin/env python3
"""
Adaptive controller for the SpliceAI Nextflow pipeline.

- Mirrors the pipeline’s CLI (mutations, reference, mapping paths, etc.).
- Launches `nextflow run main.nf …` for you and streams its output.
- Runs the spliceai-tracker automatically unless --disable_tracker is set.
- Watches .nextflow.log for exit-134 events; when <= tail_threshold genes remain
  and a tail gene fails repeatedly (>=3 within rapid_window_minutes or two within
  rapid_gap_minutes), the controller stops the current run and relaunches with
  `-resume --maxforks tail_maxforks` (typically 1).
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
NEXTFLOW_SCRIPT = HERE / "main.nf"
TRACKER_SCRIPT = HERE / "spliceai-tracker.py"
LOG_PATH = HERE / ".nextflow.log"
WORK_ROOT = HERE / "work"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Adaptive wrapper for the SpliceAI Nextflow pipeline."
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
                        help="Base maxforks value (controller supplies per run; 0=Nextflow default)")

    # Controller knobs
    parser.add_argument("--initial_maxforks", type=int, default=3,
                        help="Max concurrent run_spliceai tasks on the first pass")
    parser.add_argument("--tail_maxforks", type=int, default=1,
                        help="Maxforks used after tail-mode restart (usually 1)")
    parser.add_argument("--tail_threshold", type=int, default=6,
                        help="Genes remaining at which we consider the run to be in the tail")
    parser.add_argument("--rapid_window_minutes", type=float, default=10.0,
                        help="Window for detecting >=3 failures of the same gene")
    parser.add_argument("--rapid_gap_minutes", type=float, default=5.0,
                        help="Threshold for two back-to-back failures")
    parser.add_argument("--poll_seconds", type=float, default=15.0,
                        help="Monitor interval for tail detection")
    parser.add_argument("--disable_tracker", action="store_true",
                        help="Skip launching spliceai-tracker.py")
    parser.add_argument("--pipeline", type=Path, default=NEXTFLOW_SCRIPT,
                        help="Nextflow script to run (default: main.nf)")
    parser.add_argument("--resume", action="store_true",
                        help="Start the very first launch with -resume")

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
    if args.initial_maxforks < 1 or args.tail_maxforks < 1:
        raise SystemExit("ERROR: maxforks values must be >= 1")


def normalize(path: Optional[Path]) -> Optional[str]:
    return str(path.resolve()) if path else None


def discover_gene_ids(args: argparse.Namespace) -> List[str]:
    if args.input_vcf_path:
        p = args.input_vcf_path
        if p.is_file():
            return [p.stem.replace(".vcf", "")]
        genes = sorted(v.stem.replace(".vcf", "") for v in p.glob("*.vcf"))
        if not genes:
            raise SystemExit(f"ERROR: No *.vcf files in {p}")
        return genes

    mp = args.mutations_path
    if not mp:
        raise SystemExit("ERROR: --mutations_path required when not using --input_vcf_path")
    if mp.is_file():
        return [mp.stem.replace("_mutations", "")]
    genes = sorted(csv.stem.replace("_mutations", "") for csv in mp.glob("*.csv"))
    if not genes:
        raise SystemExit(f"ERROR: No *_mutations.csv files in {mp}")
    return genes


def build_nextflow_cmd(args: argparse.Namespace, maxforks: int, resume_flag: bool) -> List[str]:
    cmd = ["nextflow", "run", str(args.pipeline)]
    if resume_flag or args.resume:
        cmd.append("-resume")

    def add_param(name: str, value: Optional[Path | str | float | int]):
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
    add_param("maxforks", maxforks if maxforks else args.maxforks)
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


class LogTailer:
    pattern = re.compile(r"^(?P<ts>\w+-\d+\s+\d+:\d+:\d+\.\d+).+run_spliceai \((?P<gene>[^)]+)\).+exit: 134")

    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.fp: Optional[object] = None
        self.pos = 0
        self.year = dt.datetime.now().year

    def _ensure_open(self) -> bool:
        if not self.log_path.exists():
            return False
        if self.fp is None or self.fp.closed:
            self.fp = self.log_path.open("r")
            self.fp.seek(0, os.SEEK_END)
            self.pos = self.fp.tell()
        return True

    def read_new_failures(self) -> List[Tuple[str, dt.datetime]]:
        events: List[Tuple[str, dt.datetime]] = []
        if not self._ensure_open():
            return events
        self.fp.seek(self.pos)
        for line in self.fp:
            self.pos = self.fp.tell()
            m = self.pattern.search(line)
            if m:
                ts = dt.datetime.strptime(f"{self.year} {m.group('ts')}", "%Y %b-%d %H:%M:%S.%f")
                events.append((m.group("gene").strip(), ts))
        return events


def collect_completed_genes() -> List[str]:
    completed = []
    for sub in WORK_ROOT.glob("*/*"):
        exit_file = sub / ".exitcode"
        if not exit_file.exists():
            continue
        try:
            status = exit_file.read_text().strip()
        except OSError:
            continue
        if status != "0":
            continue
        cmd_file = sub / ".command.sh"
        if not cmd_file.exists():
            continue
        text = cmd_file.read_text()
        marker = '-O "'
        if marker not in text:
            continue
        gene = text.split(marker, 1)[1].split(".spliceai.vcf", 1)[0]
        completed.append(gene)
    return completed


def should_restart_tail(
        tail_genes: List[str],
        failure_history: Dict[str, List[dt.datetime]],
        window: dt.timedelta,
        gap: dt.timedelta,
) -> bool:
    now = dt.datetime.now()
    tail_set = set(tail_genes)
    for gene in tail_set:
        history = [ts for ts in failure_history.get(gene, []) if now - ts <= window]
        failure_history[gene] = history
        if len(history) >= 3 and history[-1] - history[-3] <= window:
            return True
        if len(history) >= 2 and history[-1] - history[-2] <= gap:
            return True
    return False


def monitor_and_maybe_restart(
        nf_proc: subprocess.Popen,
        tracker: TrackerRunner,
        args: argparse.Namespace,
        gene_ids: List[str],
        failure_history: Dict[str, List[dt.datetime]],
) -> bool:
    tailer = LogTailer(LOG_PATH)
    window = dt.timedelta(minutes=args.rapid_window_minutes)
    gap = dt.timedelta(minutes=args.rapid_gap_minutes)

    try:
        while True:
            ret = nf_proc.poll()
            if ret is not None:
                tracker.stop()
                return False  # finished normally

            for gene, ts in tailer.read_new_failures():
                failure_history.setdefault(gene, []).append(ts)

            completed = set(collect_completed_genes())
            remaining = [g for g in gene_ids if g not in completed]
            tail_mode = 0 < len(remaining) <= args.tail_threshold

            if tail_mode and should_restart_tail(remaining, failure_history, window, gap):
                print(
                    f"[controller] Rapid exit-134 detected among {remaining}; restarting with maxforks={args.tail_maxforks}",
                    flush=True,
                )
                nf_proc.terminate()
                try:
                    nf_proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    nf_proc.kill()
                tracker.stop()
                return True

            time.sleep(max(args.poll_seconds, 5.0))
    finally:
        tracker.stop()


def run_controller(args: argparse.Namespace):
    gene_ids = discover_gene_ids(args)
    failure_history: Dict[str, List[dt.datetime]] = collections.defaultdict(list)
    current_max = args.initial_maxforks
    resume_flag = args.resume
    tracker = TrackerRunner(args)

    while True:
        cmd = build_nextflow_cmd(args, current_max, resume_flag)
        print(f"[controller] Launching: {' '.join(cmd)}", flush=True)
        nf_proc = subprocess.Popen(cmd, cwd=str(HERE))
        tracker.start()
        restart_needed = monitor_and_maybe_restart(nf_proc, tracker, args, gene_ids, failure_history)
        if restart_needed and current_max != args.tail_maxforks:
            current_max = args.tail_maxforks
            resume_flag = True  # ensure resume on subsequent run
            continue
        exit_code = nf_proc.wait()
        sys.exit(exit_code)


if __name__ == "__main__":
    run_controller(parse_args())