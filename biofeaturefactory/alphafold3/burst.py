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
AlphaFold3 burst-mode driver.

Two phases, two subcommands:

  submit:
    Walk the same gene -> mutation -> window -> RBP iteration as the LOCAL
    pipeline; for each (rna_window, protein) pair, write an AF3 input JSON,
    deduplicate by content hash, append a manifest row. Render the SLURM
    array script from the bundled template, then sbatch it (or skip with
    --no-submit). The orchestrator exits; AF3 inference happens on the
    cluster.

  ingest:
    Re-walk the same iteration deterministically. For each (mutation, RBP),
    look up the WT and MUT input hashes in the L1 cache, parse the AF3
    outputs, compute BindingMetrics + DeltaMetrics, and write the per-gene
    3-tier TSVs (summary / events / sites).

Both subcommands import the analytical layer from the existing AF3 module:
af3_parser, binding_metrics, rbp_database, rbp_sequence_mapper. No
analytical code is duplicated. The only thing this driver does that
alphafold3_pipeline.py doesn't is decouple AF3 invocation from result
ingestion via the L1 cache.

The cache layout, key, and hit criterion match BFF/AF3_CACHE_PLAN.md:
  cache root     {output}/.cache/af3/
  cache key      AF3Input.get_hash() = MD5(rna_seq + '_' + protein_seq)[:12]
  hit criterion  the cache dir contains all three primary AF3 outputs:
                 _model.cif, _confidences.json, _summary_confidences.json
"""

import argparse
import fcntl
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from biofeaturefactory.alphafold3.bin.rbp_database import POSTAR3Database, RBPBindingSite
from biofeaturefactory.alphafold3.bin.rbp_sequence_mapper import RBPSequenceMapper
from biofeaturefactory.alphafold3.bin.af3_runner import AF3Input, create_rna_protein_input
from biofeaturefactory.alphafold3.bin.af3_parser import (
    AggregatedBindingAnalysis, BindingAnalysis,
    aggregate_binding_analyses, analyze_binding, extract_interface_sites,
    parse_all_samples,
)
from biofeaturefactory.alphafold3.bin.binding_metrics import (
    BindingMetrics, DeltaMetrics, RnaEditSpan, ThresholdConfig,
    aggregate_mutation_summary, compute_delta_metrics,
    format_events_rows, format_sites_rows,
)
from biofeaturefactory.alphafold3.bin.burst_manifest import (
    ManifestRow, is_cache_complete, write_manifest,
)

from biofeaturefactory.lib.utility import (
    mint_pkey,
    Variant, _collect_failures_from_logs, extract_gene_from_filename,
    get_mutation_data_bioAccurate, load_mapping, parse_variant, read_fasta,
    splice_seq, subseq, trim_muts, write_tsv,
)


# Superset test for a nucleotide allele; mirrors alphafold3_pipeline._NT_BASES.
# NOT `allele in "ACGTU"`: substring containment accepts any allele whose
# letters happen to sit CONSECUTIVELY in that literal ('AC', 'CG', 'ACGT' all
# pass) while rejecting 'AG', 'GA' and 'ACC'.
_NT_BASES = frozenset("ACGTU")


# ---------------------------------------------------------------------------
# Internal types
# ---------------------------------------------------------------------------

@dataclass
class BurstInput:
    """One AF3 invocation in the burst plan.

    Submit emits these into the manifest; ingest regenerates them
    deterministically and looks up cache_dir by input_hash.
    """
    gene: str
    pkey: str               # 'GENE-mutation', e.g. 'F9-C123T'
    rbp_name: str
    allele: str             # 'WT' or 'MUT'
    window_idx: int         # 0 if single-window
    af3_input: AF3Input
    distance_to_mutation: int = 0   # bp distance to nearest RBP site edge
    # Edit geometry of THIS window. Carried on the input rather than recomputed
    # in ingest because ingest groups by (pkey, rbp) and no longer holds the
    # transcript coordinates the span was derived from.
    edit_span: Optional[RnaEditSpan] = None
    variant_kind: str = 'snv'       # Variant.kind

    @property
    def input_hash(self) -> str:
        return self.af3_input.get_hash()


@dataclass
class SkippedMutation:
    """A mutation the walk rejected before any AF3 input was built.

    Ingest turns each of these into one summary row, so burst and the LOCAL
    pipeline emit the same row set for the same input.
    """
    gene: str
    pkey: str
    qc_flag: str
    # Filled whenever the token parses. The class is a property of the TOKEN and
    # stays true however the token was rejected, so it is reported rather than
    # blanked -- matching AlphaFold3Pipeline._summary_stub.
    variant_kind: str = ''
    length_delta: str = ''


def _failed_flag(message: str) -> str:
    """qc_flag for a mutation the LOCAL pipeline would fail; same 50-char cut."""
    return f'FAILED:{message[:50]}'


def _acquire_burst_lock(burst_dir: Path):
    """Take a non-blocking exclusive lock on `{burst_dir}/.lock`.

    Returns the open file handle; caller must keep it in scope for the lock
    to persist. Raises RuntimeError if another submit/ingest is running.
    """
    burst_dir.mkdir(parents=True, exist_ok=True)
    lock_path = burst_dir / ".lock"
    fh = open(lock_path, 'w')
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        fh.close()
        raise RuntimeError(
            f"Another burst submit/ingest process holds {lock_path}. "
            "Wait for it to finish, or `rm` the lock file if it is stale."
        )
    return fh


def _ingest_warn_in_flight_slurm(burst_dir: Path) -> None:
    """Warn if `last_submit.json` indicates the SLURM job is still running.

    Soft check: uses `sacct` if available; silently no-ops on hosts that
    don't have it. Does NOT block ingest; the user contract is to wait
    for sacct to clear before running ingest.
    """
    last_submit = burst_dir / "last_submit.json"
    if not last_submit.exists():
        return
    try:
        info = json.loads(last_submit.read_text())
    except (json.JSONDecodeError, OSError):
        return
    job_id = info.get("job_id")
    if not job_id:
        return
    try:
        result = subprocess.run(
            ["sacct", "-j", str(job_id), "--format=State", "--noheader", "--parsable2"],
            capture_output=True, text=True, timeout=10,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        # No sacct on this host (e.g., running ingest from a non-SLURM host)
        print(f"NOTE: cannot check SLURM job {job_id} status (sacct unavailable). "
              "Ensure all array tasks have finished before trusting ingest output.",
              file=sys.stderr)
        return
    if result.returncode != 0:
        return
    terminal = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
                "NODE_FAIL", "OUT_OF_MEMORY", "BOOT_FAIL", "PREEMPTED"}
    states = {line.strip() for line in result.stdout.splitlines() if line.strip()}
    pending = states - terminal
    if pending:
        print(f"WARNING: SLURM job {job_id} has tasks still in non-terminal "
              f"state(s) {sorted(pending)}. Ingest will only see results for "
              "completed array tasks; partial completions become PARTIAL/ALL_FAILED "
              "qc_flags. Consider waiting until sacct shows all tasks terminal.",
              file=sys.stderr)


def _parse_vcf_chrom(vcf_path: Path) -> Optional[str]:
    """Extract CHROM from a per-gene VCF. Mirrors alphafold3_pipeline:60-69."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields:
                return fields[0]
    return None


# ---------------------------------------------------------------------------
# Shared input collection (submit + ingest both call this with same args)
# ---------------------------------------------------------------------------

def iterate_inputs(
    args,
    rbp_db: POSTAR3Database,
    seq_mapper: RBPSequenceMapper,
    skipped: Optional[List[SkippedMutation]] = None,
) -> Iterator[BurstInput]:
    """Walk gene -> mutation -> window -> RBP, yielding every AF3 invocation.

    Mutations that produce no AF3 invocation are appended to ``skipped``;
    they are absent from the yield stream but still owed an output row.

    Determinism contract (REQUIRED for submit -> ingest correctness):
      Submit and ingest must produce identical (gene, mutation, window_idx,
      RBP, allele) sequences in the same order. The data sources that must
      not change between submit and ingest are:

        - The transcript FASTA dir / file at args.fasta (ordering relies on
          ``sorted(fasta_input.glob('*.fasta'))``).
        - The POSTAR3 BED file backing rbp_db (tabix returns sorted output;
          changing the BED invalidates per-RBP distances and may add/drop
          RBPs near a mutation).
        - The chromosome-mapping CSVs and per-gene mutation CSVs (mutation
          order is preserved by trim_muts / list(chrom_mapping.keys())).
        - The RBP UniProt mapping at args.rbp_mapping and the MSA dir at
          args.msa_dir (changing the available RBP set invalidates the
          ``seq_mapper.get_rbp_data`` filter at the per-RBP step).

      Per-mutation iteration:
        - Windows are deduplicated by RNA sequence; window_idx is the
          post-dedup position (matches alphafold3_pipeline.py:441).
        - RBPs near the mutation are iterated in ``sorted(...)`` order over
          the dict returned by ``rbp_db.group_by_rbp(...)`` -- explicit sort,
          stable across runs.
        - Within each RBP, two BurstInputs are yielded (WT then MUT).

      What breaks determinism:
        - Re-downloading POSTAR3 between submit and ingest.
        - Editing rbp_uniprot_ids.txt or rebuilding the MSA dir between
          submit and ingest.
        - Adding/removing transcripts or mutations between phases.

      Practical effect of a determinism break: ingest still finds completed
      AF3 outputs by content hash (cache_root / input_hash), so any matching
      (rna_window, protein_seq) pair still produces correct deltas. Newly
      added inputs that submit never enqueued will appear as cache misses
      and become PARTIAL/ALL_FAILED in the qc_flags column. Removed inputs
      from ingest's iteration simply don't appear in the TSVs.

    Mirrors alphafold3_pipeline.py:process_gene + _process_mutation +
    _generate_windows + _submit_rbp_jobs but does not actually run AF3.
    """
    fasta_input = Path(args.fasta)
    mutations_input = Path(args.mutations) if args.mutations else None
    vcf_input = Path(args.vcf) if args.vcf else None
    chrom_map_input = Path(args.chromosome_mapping) if args.chromosome_mapping else None

    if fasta_input.is_dir():
        fasta_files = sorted(fasta_input.glob('*.fasta'))
    else:
        fasta_files = [fasta_input]

    for fasta_file in fasta_files:
        gene_name = extract_gene_from_filename(str(fasta_file))

        # Resolve chromosome name
        chrom = args.chrom
        if vcf_input:
            if vcf_input.is_dir():
                cands = list(vcf_input.glob(f'{gene_name}.vcf'))
                if cands:
                    chrom = _parse_vcf_chrom(cands[0])
            else:
                chrom = _parse_vcf_chrom(vcf_input)

        # Resolve chromosome mapping (provides per-mutation chromosomal positions)
        chrom_mapping: Optional[Dict[str, str]] = None
        if chrom_map_input:
            if chrom_map_input.is_dir():
                cands = list(chrom_map_input.glob(f'*{gene_name}*.csv'))
                if cands:
                    chrom_mapping = load_mapping(str(cands[0]), mapType="chromosome")
            else:
                chrom_mapping = load_mapping(str(chrom_map_input), mapType="chromosome")

        # Resolve mutations source
        mut_path: Optional[str] = None
        if mutations_input:
            if mutations_input.is_dir():
                cands = list(mutations_input.glob(f'*{gene_name}*.csv'))
                if cands:
                    mut_path = str(cands[0])
            else:
                mut_path = str(mutations_input)

        if not mut_path and not chrom_mapping:
            print(f"No mutations source for {gene_name}", file=sys.stderr)
            continue

        # Load transcript
        fasta_data = read_fasta(str(fasta_file))
        if not fasta_data:
            continue
        transcript_seq = None
        for key in ('ORF', 'orf', 'transcript', gene_name):
            if key in fasta_data:
                transcript_seq = fasta_data[key]
                break
        if transcript_seq is None:
            transcript_seq = next(iter(fasta_data.values()))

        # Load mutations list
        if mut_path:
            mutations = trim_muts(mut_path, args.validation_log, gene_name)
        else:
            mutations = list(chrom_mapping.keys())
            if args.validation_log:
                failures = _collect_failures_from_logs(args.validation_log)
                skip = failures.get(gene_name.upper(), set()) if gene_name else set()
                mutations = [m for m in mutations if m not in skip]

        for mutation in mutations:
            yield from _iterate_mutation(
                args=args,
                rbp_db=rbp_db,
                seq_mapper=seq_mapper,
                gene_name=gene_name,
                mutation=mutation,
                transcript_seq=transcript_seq,
                chrom=chrom,
                chrom_mapping=chrom_mapping,
                skipped=skipped,
            )


def _iterate_mutation(
    args,
    rbp_db: POSTAR3Database,
    seq_mapper: RBPSequenceMapper,
    gene_name: str,
    mutation: str,
    transcript_seq: str,
    chrom: Optional[str],
    chrom_mapping: Optional[Dict[str, str]],
    skipped: Optional[List[SkippedMutation]] = None,
) -> Iterator[BurstInput]:
    """Yield AF3 inputs for one mutation across all RBPs and windows.

    A mutation that yields nothing is appended to ``skipped`` carrying the
    qc_flag alphafold3_pipeline.py records for the same token.
    """
    # {GENE}-{sha}, matching alphafold3_pipeline.py so burst rows join its output.
    pkey = mint_pkey(gene_name, mutation)
    kind = ''
    delta_nt: str = ''

    def skip(qc_flag: str) -> None:
        if skipped is not None:
            skipped.append(SkippedMutation(gene=gene_name, pkey=pkey, qc_flag=qc_flag,
                                           variant_kind=kind, length_delta=delta_nt))

    # Parse mutation. parse_variant is length-aware and NEVER raises, so an
    # indel / MNV / delins arrives here as a record instead of dying inside
    # get_mutation_data_bioAccurate's int(token[1:-1]). Non-SNV is the DEFAULT
    # path; there is no flag. Mirrors alphafold3_pipeline._process_mutation.
    variant = parse_variant(mutation, is_nt=True)
    if variant is None:
        # Not a nucleotide token under the shared grammar. Fall back to the
        # ORIGINAL parser, unchanged, so the NA-vs-FAILED split -- including the
        # exact ValueError text -- is preserved to the byte.
        try:
            pos_data = get_mutation_data_bioAccurate(mutation, is_nt=False)  # AF3 burst opt-out of primitive nt-validation; local ref-check + alt-guard below
        except (ValueError, IndexError) as e:
            skip(_failed_flag(str(e)))
            return
        if pos_data[0] is None:
            skip('NA:non_nucleotide_token')
            return
        nt_pos = pos_data[0]
        wt_nt, mut_nt = pos_data[1]
        # aa/Stop tokens are out of scope for this nucleotide-only pipeline: N/A, not FAILED.
        if not (set(wt_nt.upper()) <= _NT_BASES and set(mut_nt.upper()) <= _NT_BASES):
            skip('NA:non_nucleotide_token')
            return
        # Reachable: int() accepts spellings the grammar does not ('A-12T',
        # 'A 12T', 'A1_2T'). Variant rejects pos < 1, so keep the original
        # out-of-range wording for the one case that used to produce it.
        if nt_pos < 1:
            skip(_failed_flag(f"Position {nt_pos} out of range"))
            return
        variant = Variant(pos=nt_pos, ref=wt_nt, alt=mut_nt)

    kind = variant.kind
    delta_nt = variant.length_delta
    nt_pos = variant.pos
    pos_0 = variant.pos0
    ref, alt = variant.ref, variant.alt

    # Bound the END of the REF span, not its start: a multi-base REF can begin
    # inside the transcript and run off the end. Identical to
    # `pos_0 >= len(transcript_seq)` when len(ref) == 1.
    if pos_0 < 0 or pos_0 + len(ref) > len(transcript_seq):
        skip(_failed_flag(f"Position {nt_pos} out of range"))
        return
    # REF guard over the WHOLE span; transcript_seq[pos_0] alone passes any
    # multi-base REF whose first base happens to match.
    observed = transcript_seq[pos_0:pos_0 + len(ref)]
    if observed.upper() != ref.upper():
        skip(_failed_flag(
            f"Reference mismatch at {nt_pos}: expected {ref}, "
            f"found {observed}"
        ))
        return

    # Genomic position (for POSTAR3 lookup + distance computation)
    genomic_pos: Optional[int] = None
    if chrom_mapping and mutation in chrom_mapping:
        entry = chrom_mapping[mutation]
        # parse_variant, not int(entry[1:-1]). The genomic spelling of a non-SNV
        # carries multi-base alleles and int() raises on every one of them;
        # widening the token parser above WITHOUT this one would leave
        # genomic_pos None for every indel and report each as
        # 'no_rbps_in_region' -- a fabricated negative, strictly worse than the
        # FAILED row it produces today.
        genomic_variant = parse_variant(entry, is_nt=True)
        if genomic_variant is not None:
            genomic_pos = genomic_variant.pos
    elif args.tx_start is not None:
        if args.strand == '+':
            genomic_pos = args.tx_start + nt_pos - 1
        else:
            genomic_pos = args.tx_start - nt_pos + 1

    # "We looked and found nothing" and "we never looked" are different
    # findings; reporting the second as the first is a fabricated observation.
    if chrom is None:
        skip('NA:no_chromosome')
        return
    if genomic_pos is None:
        skip('NA:unresolved_genomic_position')
        return

    # POSTAR3 BED intervals are 0-based half-open; genomic_pos is 1-based
    bed_pos = genomic_pos - 1

    # Windows (post-dedup; window_idx is the position in the deduped list).
    # splice_seq honours len(ref); the previous one-base concatenation cannot
    # express an indel. The mutant window covers the SAME BIOLOGICAL SPAN as the
    # WT one, so it is longer/shorter by exactly delta. The window is centred on
    # the MIDPOINT of the REF span (len(ref)//2 == 0 for an SNV).
    # Mirrors alphafold3_pipeline._generate_windows.
    mut_seq = splice_seq(transcript_seq, pos_0, ref, alt, validate=False)
    delta = len(alt) - len(ref)
    centre = min(pos_0 + len(ref) // 2, len(transcript_seq) - 1)
    windows: List[Tuple[str, str, int]] = []
    if not args.multi_window:
        wt_w = subseq(transcript_seq, centre, args.window_size)
        window_start = max(0, centre - args.window_size // 2)
        mut_w = mut_seq[window_start:window_start + len(wt_w) + delta]
        windows = [(wt_w, mut_w, pos_0 - window_start)]
    else:
        offsets = [float(x) for x in args.multi_window_offsets.split(',')]
        seen: set = set()
        for frac in offsets:
            target_center = int(frac * args.window_size)
            window_start = centre - target_center
            window_start = max(0, min(window_start, len(transcript_seq) - args.window_size))
            window_end = window_start + args.window_size
            wt_w = transcript_seq[window_start:window_end]
            mut_w = mut_seq[window_start:window_end + delta]
            if wt_w in seen:
                continue
            seen.add(wt_w)
            windows.append((wt_w, mut_w, pos_0 - window_start))

    if not windows:
        skip('no_rbps_tested')
        return

    # The variant must sit inside its own window for the WT/MUT comparison to
    # cover the edit at all. Unreachable for an SNV.
    bad = [w for w in windows if not (0 <= w[2] and w[2] + len(ref) <= len(w[0]))]
    if bad:
        skip(_failed_flag(
            f"REF span outside window: offset {bad[0][2]}, "
            f"len(ref) {len(ref)}, window {len(bad[0][0])}"))
        return

    # RBPs near the mutation, queried across the WHOLE REF span. With
    # len(ref) == 1 the interval is exactly query_position(chrom, bed_pos, w).
    rbp_sites = rbp_db.query(chrom, bed_pos - args.rbp_window,
                             bed_pos + len(ref) + args.rbp_window)
    if not rbp_sites:
        skip('no_rbps_in_region')
        return
    rbps_to_test = rbp_db.group_by_rbp(rbp_sites)

    n_yielded = 0
    for rbp_name in sorted(rbps_to_test.keys()):
        sites = rbps_to_test[rbp_name]
        rbp_data = seq_mapper.get_rbp_data(rbp_name)
        if not rbp_data:
            continue
        protein_seq = rbp_data.sequence
        protein_msa = rbp_data.msa_content
        # Token cap; mirrors AlphaFold3Pipeline._submit_rbp_jobs. Tested against
        # the LONGER of the two windows: an insertion makes the mutant window
        # longer than the WT one, so sizing on the WT alone would enqueue a MUT
        # job that exceeds the limit.
        if max(len(windows[0][0]), len(windows[0][1])) + len(protein_seq) > 5000:
            continue
        # Per-RBP distance to mutation; same value for both alleles + all
        # windows. Measured from the whole REF span: a deletion that straddles a
        # binding site has neither endpoint inside it.
        distance = min(site.distance_to_span(bed_pos, len(ref)) for site in sites)

        for win_idx, (wt_w, mut_w, offset) in enumerate(windows):
            span = RnaEditSpan(offset, len(ref), len(alt), len(wt_w), len(mut_w))
            suffix = f"_w{win_idx}" if args.multi_window else ""
            wt_input = create_rna_protein_input(
                job_name=f"{pkey}_{rbp_name}_WT{suffix}",
                rna_seq=wt_w, protein_seq=protein_seq, protein_msa=protein_msa,
            )
            mut_input = create_rna_protein_input(
                job_name=f"{pkey}_{rbp_name}_MUT{suffix}",
                rna_seq=mut_w, protein_seq=protein_seq, protein_msa=protein_msa,
            )
            yield BurstInput(
                gene=gene_name, pkey=pkey, rbp_name=rbp_name,
                allele='WT', window_idx=win_idx, af3_input=wt_input,
                distance_to_mutation=distance, edit_span=span, variant_kind=kind,
            )
            yield BurstInput(
                gene=gene_name, pkey=pkey, rbp_name=rbp_name,
                allele='MUT', window_idx=win_idx, af3_input=mut_input,
                distance_to_mutation=distance, edit_span=span, variant_kind=kind,
            )
            n_yielded += 2

    if n_yielded == 0:
        skip('no_rbps_tested')


# ---------------------------------------------------------------------------
# Submit subcommand
# ---------------------------------------------------------------------------

def cmd_submit(args) -> int:
    output_dir = Path(args.output)
    burst_dir = output_dir / ".burst"
    inputs_dir = burst_dir / "inputs"
    cache_root = output_dir / ".cache" / "af3"
    log_dir = Path(args.slurm_log_dir) if args.slurm_log_dir else (burst_dir / "logs")

    burst_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir.mkdir(parents=True, exist_ok=True)
    cache_root.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Hold an exclusive lock for the duration of submit; prevents concurrent
    # submit + submit or submit + ingest. Released when this function returns.
    _lock_fh = _acquire_burst_lock(burst_dir)

    # Optional pre-run cache wipe
    if args.clear_cache:
        ext_cache = output_dir / ".cache"
        if ext_cache.exists():
            shutil.rmtree(ext_cache)
            print(f"Cleared cache at {ext_cache}", file=sys.stderr)
        cache_root.mkdir(parents=True, exist_ok=True)

    # Initialize POSTAR3 + RBP mapper (same as LOCAL pipeline)
    print("Loading POSTAR3 database...", file=sys.stderr)
    rbp_db = POSTAR3Database(args.postar_db)
    print("Loading RBP sequence mapper...", file=sys.stderr)
    seq_mapper = RBPSequenceMapper(
        mapping_file=args.rbp_mapping,
        sequence_fasta=args.rbp_sequences,
        msa_dir=args.msa_dir,
    )

    # Walk inputs, dedupe by hash, write manifest entries for pending work
    manifest_rows: List[ManifestRow] = []
    skipped_muts: List[SkippedMutation] = []
    seen_hashes: set = set()
    skipped_cached = 0
    skipped_dup = 0
    total_inputs = 0

    for bi in iterate_inputs(args, rbp_db, seq_mapper, skipped=skipped_muts):
        total_inputs += 1
        h = bi.input_hash
        if h in seen_hashes:
            skipped_dup += 1
            continue
        seen_hashes.add(h)

        cache_dir = cache_root / h
        if is_cache_complete(cache_dir):
            skipped_cached += 1
            continue

        # Write input JSON (idempotent: skip if already written)
        input_json_path = inputs_dir / f"{h}.json"
        if not input_json_path.exists():
            with open(input_json_path, 'w') as f:
                json.dump(bi.af3_input.to_json_dict(), f, indent=2)

        manifest_rows.append(ManifestRow(
            array_idx=len(manifest_rows),
            input_hash=h,
            pkey=bi.pkey,
            rbp_name=bi.rbp_name,
            allele=bi.allele,
            window_idx=bi.window_idx,
            input_json_path=str(input_json_path.resolve()),
            cache_dir=str(cache_dir.resolve()),
        ))

    n_pending = len(manifest_rows)
    print(
        f"Inputs walked: {total_inputs}; pending: {n_pending}; "
        f"already cached: {skipped_cached}; duplicate-hash skipped: {skipped_dup}",
        file=sys.stderr,
    )
    if skipped_muts:
        print(f"Mutations with no AF3 input: {len(skipped_muts)} "
              f"(ingest emits one summary row for each)", file=sys.stderr)

    manifest_path = burst_dir / "manifest.tsv"
    write_manifest(manifest_rows, manifest_path)
    print(f"Wrote manifest: {manifest_path}", file=sys.stderr)

    if n_pending == 0:
        print("Nothing to submit; all inputs already cached.", file=sys.stderr)
        print(f"Run ingest: python -m biofeaturefactory.alphafold3.burst ingest "
              f"--output {output_dir} ...", file=sys.stderr)
        return 0

    # Render SLURM array script
    template_path = Path(__file__).parent / "bin" / "slurm_array.sh.tmpl"
    if not template_path.exists():
        print(f"ERR: SLURM template not found at {template_path}", file=sys.stderr)
        return 2
    template = template_path.read_text()

    substitutions = {
        "__JOB_NAME__":     args.slurm_job_name,
        "__PARTITION__":    args.slurm_partition,
        "__TIME__":         args.slurm_time,
        "__MEM__":          args.slurm_mem,
        "__LOG_DIR__":      str(log_dir.resolve()),
        "__MANIFEST__":     str(manifest_path.resolve()),
        "__MODEL_DIR__":    str(Path(args.model_dir).resolve()),
        "__DOCKER_IMAGE__": args.docker_image,
    }
    rendered = template
    for k, v in substitutions.items():
        rendered = rendered.replace(k, v)
    script_path = burst_dir / "run.slurm"
    script_path.write_text(rendered)
    script_path.chmod(0o755)
    print(f"Rendered SLURM script: {script_path}", file=sys.stderr)

    if args.no_submit:
        print("--no-submit: skipping sbatch.", file=sys.stderr)
        print(f"To submit manually:", file=sys.stderr)
        print(f"  sbatch --array=0-{n_pending - 1}%{args.array_throttle} {script_path}",
              file=sys.stderr)
        return 0

    sbatch_cmd = [
        "sbatch",
        f"--array=0-{n_pending - 1}%{args.array_throttle}",
        str(script_path),
    ]
    print(f"Submitting: {' '.join(sbatch_cmd)}", file=sys.stderr)
    result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"sbatch failed (rc={result.returncode}):", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return result.returncode

    job_id: Optional[str] = None
    m = re.search(r'Submitted batch job (\d+)', result.stdout)
    if m:
        job_id = m.group(1)
    print(f"Submitted SLURM job {job_id} ({n_pending} tasks, throttle %{args.array_throttle})",
          file=sys.stderr)

    last_submit_path = burst_dir / "last_submit.json"
    with open(last_submit_path, 'w') as f:
        json.dump({
            "job_id": job_id,
            "manifest_path": str(manifest_path.resolve()),
            "submit_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "n_tasks": n_pending,
            "throttle": args.array_throttle,
            "script": str(script_path.resolve()),
        }, f, indent=2)
    print(f"Wrote {last_submit_path}", file=sys.stderr)
    print(f"\nWhen the job completes, run ingest:", file=sys.stderr)
    print(f"  python -m biofeaturefactory.alphafold3.burst ingest --output {output_dir} \\",
          file=sys.stderr)
    print(f"      --postar-db ... --rbp-mapping ... --msa-dir ... --fasta ...", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Ingest subcommand
# ---------------------------------------------------------------------------

def cmd_ingest(args) -> int:
    output_dir = Path(args.output)
    burst_dir = output_dir / ".burst"
    cache_root = output_dir / ".cache" / "af3"

    if not cache_root.exists():
        print(f"No AF3 cache at {cache_root}; nothing to ingest.", file=sys.stderr)
        return 1

    # Hold the burst lock so ingest can't race a concurrent submit.
    _lock_fh = _acquire_burst_lock(burst_dir)
    # Soft-warn if SLURM tasks are still in flight (does not block ingest).
    _ingest_warn_in_flight_slurm(burst_dir)

    print("Loading POSTAR3 database...", file=sys.stderr)
    rbp_db = POSTAR3Database(args.postar_db)
    print("Loading RBP sequence mapper...", file=sys.stderr)
    seq_mapper = RBPSequenceMapper(
        mapping_file=args.rbp_mapping,
        sequence_fasta=args.rbp_sequences,
        msa_dir=args.msa_dir,
    )

    threshold_config = ThresholdConfig()

    # Group BurstInputs by (pkey, rbp_name) so we can pair WT+MUT across
    # windows. inputs_by_winallele maps (window_idx, allele) -> BurstInput.
    grouped: Dict[Tuple[str, str], Dict[Tuple[int, str], BurstInput]] = {}
    skipped_muts: List[SkippedMutation] = []
    n_inputs = 0
    for bi in iterate_inputs(args, rbp_db, seq_mapper, skipped=skipped_muts):
        n_inputs += 1
        grouped.setdefault((bi.pkey, bi.rbp_name), {})[(bi.window_idx, bi.allele)] = bi

    print(f"Walked {n_inputs} expected inputs across {len(grouped)} "
          f"(mutation, RBP) pairs; {len(skipped_muts)} mutation(s) yielded none",
          file=sys.stderr)

    # Per-gene accumulators
    summary_rows: Dict[str, List[dict]] = {}
    events_rows:  Dict[str, List[dict]] = {}
    sites_rows:   Dict[str, List[dict]] = {}
    pkey_to_gene: Dict[str, str] = {}
    pkey_to_deltas: Dict[str, List[DeltaMetrics]] = {}

    for (pkey, rbp_name), inputs_by_winallele in grouped.items():
        any_input = next(iter(inputs_by_winallele.values()))
        gene = any_input.gene
        distance = any_input.distance_to_mutation
        pkey_to_gene[pkey] = gene
        pkey_to_deltas.setdefault(pkey, [])

        win_indices = sorted({wi for (wi, _) in inputs_by_winallele.keys()})

        wt_metrics_per_win: List[Optional[BindingMetrics]] = []
        mut_metrics_per_win: List[Optional[BindingMetrics]] = []
        first_wt_sites: Optional[list] = None
        first_mut_sites: Optional[list] = None
        first_wt_freq_rna: Optional[Dict[int, float]] = None
        first_wt_freq_prot: Optional[Dict[int, float]] = None
        first_mut_freq_rna: Optional[Dict[int, float]] = None
        first_mut_freq_prot: Optional[Dict[int, float]] = None

        for win_idx in win_indices:
            wt_bi = inputs_by_winallele.get((win_idx, 'WT'))
            mut_bi = inputs_by_winallele.get((win_idx, 'MUT'))

            wt_m, wt_sites, wt_agg = _parse_cache_entry(
                cache_root / wt_bi.input_hash if wt_bi else None,
                rbp_name=rbp_name, threshold_config=threshold_config,
            )
            mut_m, mut_sites, mut_agg = _parse_cache_entry(
                cache_root / mut_bi.input_hash if mut_bi else None,
                rbp_name=rbp_name, threshold_config=threshold_config,
            )
            wt_metrics_per_win.append(wt_m)
            mut_metrics_per_win.append(mut_m)

            if first_wt_sites is None and wt_sites is not None:
                first_wt_sites = wt_sites
                first_wt_freq_rna = wt_agg.contact_frequency_rna if wt_agg else None
                first_wt_freq_prot = wt_agg.contact_frequency_protein if wt_agg else None
            if first_mut_sites is None and mut_sites is not None:
                first_mut_sites = mut_sites
                first_mut_freq_rna = mut_agg.contact_frequency_rna if mut_agg else None
                first_mut_freq_prot = mut_agg.contact_frequency_protein if mut_agg else None

        if len(win_indices) == 1:
            wt_metrics = wt_metrics_per_win[0]
            mut_metrics = mut_metrics_per_win[0]
        else:
            wt_metrics = _aggregate_across_windows(
                wt_metrics_per_win, rbp_name, threshold_config,
            )
            mut_metrics = _aggregate_across_windows(
                mut_metrics_per_win, rbp_name, threshold_config,
            )

        delta = compute_delta_metrics(
            rbp_name=rbp_name,
            wt_metrics=wt_metrics, mut_metrics=mut_metrics,
            distance_to_mutation=distance, config=threshold_config,
        )
        if len(win_indices) > 1:
            delta.n_windows = len(win_indices)
        pkey_to_deltas[pkey].append(delta)

        if first_wt_sites:
            sites_rows.setdefault(gene, []).extend(format_sites_rows(
                pkey, rbp_name, 'WT', first_wt_sites,
                first_wt_freq_rna, first_wt_freq_prot,
            ))
        if first_mut_sites:
            sites_rows.setdefault(gene, []).extend(format_sites_rows(
                pkey, rbp_name, 'MUT', first_mut_sites,
                first_mut_freq_rna, first_mut_freq_prot,
            ))

    # Per-mutation summaries + events
    for pkey, deltas in pkey_to_deltas.items():
        gene = pkey_to_gene[pkey]
        summary = aggregate_mutation_summary(deltas)
        summary['pkey'] = pkey
        summary['Gene'] = gene
        has_complete = any(d.wt_metrics is not None and d.mut_metrics is not None for d in deltas)
        has_partial = any(d.wt_metrics is not None or d.mut_metrics is not None for d in deltas)
        if not deltas:
            summary['qc_flags'] = 'no_rbps_tested'
        elif has_complete:
            summary['qc_flags'] = 'PASS'
        elif has_partial:
            summary['qc_flags'] = 'PARTIAL'
        else:
            summary['qc_flags'] = 'ALL_FAILED'
        summary_rows.setdefault(gene, []).append(summary)
        events_rows.setdefault(gene, []).extend(format_events_rows(pkey, deltas))

    for sk in skipped_muts:
        # aggregate_mutation_summary([]) supplies the full column set; write_tsv
        # takes its header from rows[0], so a short row would drop columns.
        summary = aggregate_mutation_summary([])
        summary['pkey'] = sk.pkey
        summary['Gene'] = sk.gene
        summary['qc_flags'] = sk.qc_flag
        summary_rows.setdefault(sk.gene, []).append(summary)

    def _union_fieldnames(rows):
        """Header = union of every row's keys, first-appearance order.

        write_tsv defaults to fieldnames=rows[0].keys() with extrasaction='ignore',
        so a first row narrower than a later one silently truncates the file.
        Mirrors AlphaFold3Pipeline._write_rows in alphafold3_pipeline.py -- the two
        drivers must agree on the schema they emit.
        """
        return list(dict.fromkeys(key for row in rows for key in row))

    # Flush per-gene TSVs
    all_genes = set(summary_rows) | set(events_rows) | set(sites_rows)
    for gene in sorted(all_genes):
        out_dir = output_dir / gene / "AlphaFold3"
        out_dir.mkdir(parents=True, exist_ok=True)
        if summary_rows.get(gene):
            write_tsv(summary_rows[gene], str(out_dir / f"{gene}.tsv"),
                      _union_fieldnames(summary_rows[gene]))
            print(f"Wrote {len(summary_rows[gene])} summary rows to "
                  f"{out_dir / f'{gene}.tsv'}", file=sys.stderr)
        if events_rows.get(gene):
            write_tsv(events_rows[gene], str(out_dir / f"{gene}.events.tsv"),
                      _union_fieldnames(events_rows[gene]))
            print(f"Wrote {len(events_rows[gene])} events rows to "
                  f"{out_dir / f'{gene}.events.tsv'}", file=sys.stderr)
        if sites_rows.get(gene):
            write_tsv(sites_rows[gene], str(out_dir / f"{gene}.sites.tsv"),
                      _union_fieldnames(sites_rows[gene]))
            print(f"Wrote {len(sites_rows[gene])} sites rows to "
                  f"{out_dir / f'{gene}.sites.tsv'}", file=sys.stderr)

    return 0


def _parse_cache_entry(
    cache_dir: Optional[Path],
    rbp_name: str,
    threshold_config: ThresholdConfig,
) -> Tuple[Optional[BindingMetrics], Optional[list], Optional[AggregatedBindingAnalysis]]:
    """Parse one L1 cache directory into (metrics, sites, aggregation).

    Returns (None, None, None) when the cache dir is missing or incomplete.
    Mirrors alphafold3_pipeline.py:_parse_af3_output:551-598.
    """
    if cache_dir is None or not is_cache_complete(cache_dir):
        return None, None, None

    structures = parse_all_samples(str(cache_dir))
    if not structures:
        return None, None, None

    analyses = [analyze_binding(s, rna_chain="R", protein_chain="P") for s in structures]
    sites = extract_interface_sites(structures[0]) if structures else None

    if len(structures) == 1:
        binding = analyses[0]
        if not binding:
            return None, sites, None
        metrics = BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=binding.chain_pair_pae_min,
            interface_contacts=binding.n_contacts,
            interface_plddt_rna=binding.interface_plddt_rna,
            interface_plddt_protein=binding.interface_plddt_protein,
            has_binding=binding.n_contacts >= threshold_config.min_contacts,
        )
        return metrics, sites, None

    agg = aggregate_binding_analyses(analyses)
    if not agg:
        return None, sites, None

    metrics = BindingMetrics(
        rbp_name=rbp_name,
        chain_pair_pae_min=agg.mean.chain_pair_pae_min,
        interface_contacts=agg.mean.n_contacts,
        interface_plddt_rna=agg.mean.interface_plddt_rna,
        interface_plddt_protein=agg.mean.interface_plddt_protein,
        has_binding=agg.mean.n_contacts >= threshold_config.min_contacts,
        n_samples=agg.n_samples,
        std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
        std_interface_contacts=agg.std_n_contacts,
        std_plddt_rna=agg.std_interface_plddt_rna,
        std_plddt_protein=agg.std_interface_plddt_protein,
    )
    return metrics, sites, agg


def _aggregate_across_windows(
    metrics_per_win: List[Optional[BindingMetrics]],
    rbp_name: str,
    threshold_config: ThresholdConfig,
) -> Optional[BindingMetrics]:
    """Aggregate per-window BindingMetrics into one (mean across windows).

    Mirrors alphafold3_pipeline.py:_aggregate_parsed_results:600-640.
    """
    valid = [m for m in metrics_per_win if m is not None]
    if not valid:
        return None
    analyses = [
        BindingAnalysis(
            rna_chain="R", protein_chain="P",
            n_contacts=m.interface_contacts,
            min_contact_distance=0.0,
            mean_contact_distance=0.0,
            interface_plddt_rna=m.interface_plddt_rna,
            interface_plddt_protein=m.interface_plddt_protein,
            chain_pair_pae_min=m.chain_pair_pae_min,
            contact_residues_rna=[], contact_residues_protein=[],
        )
        for m in valid
    ]
    agg = aggregate_binding_analyses(analyses)
    if not agg:
        return None
    return BindingMetrics(
        rbp_name=rbp_name,
        chain_pair_pae_min=agg.mean.chain_pair_pae_min,
        interface_contacts=agg.mean.n_contacts,
        interface_plddt_rna=agg.mean.interface_plddt_rna,
        interface_plddt_protein=agg.mean.interface_plddt_protein,
        has_binding=agg.mean.n_contacts >= threshold_config.min_contacts,
        n_samples=agg.n_samples,
        std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
        std_interface_contacts=agg.std_n_contacts,
        std_plddt_rna=agg.std_interface_plddt_rna,
        std_plddt_protein=agg.std_interface_plddt_protein,
    )


# ---------------------------------------------------------------------------
# Preflight subcommand
# ---------------------------------------------------------------------------

def cmd_preflight(args) -> int:
    """Verify cluster preconditions before submit.

    Hard checks (failure -> exit 1):
      - docker on PATH
      - AF3 Docker image present locally
      - model_dir exists, readable, non-empty
      - output dir writable
      - postar_db readable
      - rbp_mapping readable

    Soft checks (failure -> warning, exit 0):
      - sbatch on PATH (preflight may run on a non-SLURM host)
      - nvidia-smi works (preflight host may not have a GPU)
    """
    failures: List[str] = []
    warnings: List[str] = []

    def report(name: str, ok: bool, msg: str, hard: bool = True) -> None:
        status = "OK" if ok else ("FAIL" if hard else "WARN")
        print(f"  [{status}] {name}: {msg}", file=sys.stderr)
        if not ok:
            (failures if hard else warnings).append(name)

    print("Preflight checks:", file=sys.stderr)

    # docker on PATH
    docker_path = shutil.which("docker")
    report("docker on PATH", docker_path is not None,
           f"{docker_path}" if docker_path else "not found")

    # AF3 image present
    if docker_path:
        result = subprocess.run(
            ["docker", "image", "inspect", args.docker_image],
            capture_output=True, text=True,
        )
        report(f"docker image '{args.docker_image}'", result.returncode == 0,
               "present" if result.returncode == 0
               else f"missing (build with `docker build -t {args.docker_image} -f alphafold3/docker/Dockerfile .`)")
    else:
        report(f"docker image '{args.docker_image}'", False,
               "skipped (docker not on PATH)")

    # model_dir
    model_dir = Path(args.model_dir)
    if not model_dir.exists():
        report("model_dir", False, f"{model_dir} does not exist")
    elif not model_dir.is_dir():
        report("model_dir", False, f"{model_dir} is not a directory")
    else:
        contents = list(model_dir.iterdir())
        if not contents:
            report("model_dir", False, f"{model_dir} is empty")
        else:
            report("model_dir", True, f"{model_dir} ({len(contents)} entries)")

    # output writable
    output_dir = Path(args.output)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        probe = output_dir / ".burst_preflight_probe"
        probe.write_text("ok")
        probe.unlink()
        report("output writable", True, str(output_dir))
    except OSError as e:
        report("output writable", False, f"{output_dir}: {e}")

    # postar_db
    postar = Path(args.postar_db)
    if not postar.exists():
        report("postar_db", False, f"{postar} does not exist")
    else:
        size_mb = postar.stat().st_size / (1024 * 1024)
        tbi = postar.with_suffix(postar.suffix + ".tbi")
        if str(postar).endswith(".gz") and not tbi.exists():
            report("postar_db", True,
                   f"{postar} ({size_mb:.0f} MB) — note: no .tbi index found, "
                   "tabix queries will fall back to slower in-memory mode")
        else:
            report("postar_db", True, f"{postar} ({size_mb:.0f} MB)")

    # rbp_mapping
    rbp_map = Path(args.rbp_mapping)
    if not rbp_map.exists():
        report("rbp_mapping", False, f"{rbp_map} does not exist")
    else:
        n_lines = sum(1 for _ in open(rbp_map))
        report("rbp_mapping", True, f"{rbp_map} ({n_lines} lines)")

    # sbatch (soft)
    sbatch_path = shutil.which("sbatch")
    report("sbatch on PATH", sbatch_path is not None,
           f"{sbatch_path}" if sbatch_path
           else "not found (preflight running on non-SLURM host?)",
           hard=False)

    # nvidia-smi (soft)
    nvidia_path = shutil.which("nvidia-smi")
    if nvidia_path:
        result = subprocess.run([nvidia_path], capture_output=True, text=True, timeout=10)
        report("nvidia-smi", result.returncode == 0,
               f"{nvidia_path}" if result.returncode == 0
               else f"present but exited {result.returncode}",
               hard=False)
    else:
        report("nvidia-smi", False,
               "not found (preflight host may not have GPU; compute nodes still need it)",
               hard=False)

    print(f"\nPreflight summary: {len(failures)} failure(s), {len(warnings)} warning(s)",
          file=sys.stderr)
    if failures:
        print(f"Failed checks: {', '.join(failures)}", file=sys.stderr)
        return 1
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_shared_input_args(p: argparse.ArgumentParser) -> None:
    """Args needed by both submit and ingest (input collection)."""
    p.add_argument('--postar-db', required=True, help='POSTAR3 database file')
    p.add_argument('--rbp-mapping', required=True, help='Gene-UniProt mapping TSV')
    p.add_argument('-rs', '--rbp-sequences', help='Protein sequences FASTA (optional if --msa-dir provided)')
    p.add_argument('-md', '--msa-dir', help='Directory with A3M MSA files (preferred)')
    p.add_argument('-f', '--fasta', required=True, help='Transcript FASTA file or directory')
    p.add_argument('-mu', '--mutations', help='Mutations CSV file or directory')
    p.add_argument('-v', '--vcf', help='Per-gene VCF file or directory (for chromosome resolution)')
    p.add_argument('-cm', '--chromosome-mapping', help='Chromosome mapping CSV file or directory')
    p.add_argument('-ch', '--chrom', help='Chromosome (alternative to --vcf)')
    p.add_argument('-ts', '--tx-start', type=int,
                   help='Transcript start (alternative to --chromosome-mapping)')
    p.add_argument('-s', '--strand', default='+', choices=['+', '-'], help='Strand')
    p.add_argument('-ws', '--window-size', type=int, default=101, help='RNA window size (odd)')
    p.add_argument('-rw', '--rbp-window', type=int, default=50, help='Window for RBP site lookup (+/-bp)')
    p.add_argument('-vl', '--validation-log', help='Validation log for filtering mutations')
    p.add_argument('-mw', '--multi-window', action='store_true', default=False,
                   help='Multi-window mode (multiplies AF3 runs)')
    p.add_argument('-mwo', '--multi-window-offsets', default='0.3,0.5,0.7',
                   help='Mutation positions as fractions of window')


def main() -> int:
    parser = argparse.ArgumentParser(
        prog='biofeaturefactory.alphafold3.burst',
        description='AlphaFold3 burst-mode driver: SLURM submit + cache ingest',
    )
    sub = parser.add_subparsers(dest='subcommand', required=True)

    p_sub = sub.add_parser('submit', help='Build inputs, write manifest, submit SLURM array')
    _add_shared_input_args(p_sub)
    p_sub.add_argument('--output', '-o', required=True, help='Output base directory')
    p_sub.add_argument('--model-dir', required=True, help='Path to AF3 model weights directory')
    p_sub.add_argument('--docker-image', default='alphafold3', help='AF3 Docker image name')
    p_sub.add_argument('-sjn', '--slurm-job-name', default='af3_burst', help='SLURM --job-name')
    p_sub.add_argument('-sp', '--slurm-partition', default='gpu', help='SLURM --partition')
    p_sub.add_argument('-st', '--slurm-time', default='01:00:00', help='SLURM --time per task')
    p_sub.add_argument('-sm', '--slurm-mem', default='64G', help='SLURM --mem per task')
    p_sub.add_argument('-sld', '--slurm-log-dir',
                       help='Directory for SLURM logs (default: {output}/.burst/logs)')
    p_sub.add_argument('-at', '--array-throttle', type=int, default=256,
                       help='Concurrent array task limit (sbatch --array=...%%N)')
    p_sub.add_argument('-ns', '--no-submit', action='store_true', default=False,
                       help='Render manifest + script but skip sbatch')
    p_sub.add_argument('-cc', '--clear-cache', action='store_true', default=False,
                       help='Delete {output}/.cache before submitting (forces re-run)')
    p_sub.set_defaults(func=cmd_submit)

    p_ing = sub.add_parser('ingest', help='Read AF3 cache, compute deltas, write 3-tier TSVs')
    _add_shared_input_args(p_ing)
    p_ing.add_argument('--output', '-o', required=True, help='Output base directory')
    p_ing.set_defaults(func=cmd_ingest)

    p_pre = sub.add_parser('preflight', help='Verify cluster preconditions (Docker, image, model_dir, etc.)')
    p_pre.add_argument('--output', '-o', required=True, help='Output base directory (writability check)')
    p_pre.add_argument('--model-dir', required=True, help='Path to AF3 model weights directory')
    p_pre.add_argument('--postar-db', required=True, help='POSTAR3 database file (readability check)')
    p_pre.add_argument('--rbp-mapping', required=True, help='Gene-UniProt mapping TSV (readability check)')
    p_pre.add_argument('--docker-image', default='alphafold3', help='AF3 Docker image name')
    p_pre.set_defaults(func=cmd_preflight)

    args = parser.parse_args()

    # Per-subcommand input validation; preflight uses its own minimal arg set.
    if args.subcommand in ('submit', 'ingest'):
        if not args.msa_dir and not args.rbp_sequences:
            parser.error("Provide either --msa-dir or --rbp-sequences")
        if not args.mutations and not args.chromosome_mapping:
            parser.error("Provide --mutations or --chromosome-mapping")
        if not args.chrom and not args.vcf:
            parser.error("Provide --chrom or --vcf for POSTAR3 chromosome lookup")
        if args.mutations and not args.chromosome_mapping and args.tx_start is None:
            parser.error("--tx-start is required when using --mutations without --chromosome-mapping")

    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
