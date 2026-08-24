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
AlphaFold3 RNA-RBP Interaction Pipeline

Analyzes how mutations affect RNA-binding protein interactions
using AlphaFold3 structure predictions.

Outputs:
- summary.tsv: One row per mutation with aggregated RBP binding changes
- events.tsv: One row per mutation-RBP pair with detailed metrics
- sites.tsv: Per-position confidence scores at RNA-protein interface
"""

import argparse
import csv
import sys
from pathlib import Path
from concurrent.futures import Future
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field

# Local imports
from biofeaturefactory.alphafold3.bin.rbp_database import POSTAR3Database, RBPBindingSite
from biofeaturefactory.alphafold3.bin.rbp_sequence_mapper import RBPSequenceMapper
from biofeaturefactory.alphafold3.bin.af3_runner import AF3Runner, AF3RunnerConfig, ExecutionMode, create_rna_protein_input
from biofeaturefactory.alphafold3.bin.af3_parser import (
    AF3Parser, AF3Structure, analyze_binding, BindingAnalysis,
    parse_all_samples, aggregate_binding_analyses, AggregatedBindingAnalysis,
    extract_interface_sites
)
from biofeaturefactory.alphafold3.bin.binding_metrics import (
    BindingMetrics, DeltaMetrics, ThresholdConfig, RnaEditSpan,
    compute_delta_metrics, aggregate_mutation_summary,
    format_events_rows, format_sites_rows
)

from biofeaturefactory.lib.utility import (
    derive_mutations_root,
    derive_mapping_root,
    mint_pkey,
    read_fasta, trim_muts, get_mutation_data_bioAccurate,
    extract_gene_from_filename, subseq, load_mapping,
    _collect_failures_from_logs, write_tsv,
    Variant, parse_variant, splice_seq,
    is_intronic_token, INTRONIC_PREFIX
)


# Superset test for a nucleotide allele. NOT `allele in "ACGTU"`: substring
# containment accepts any allele whose letters happen to sit CONSECUTIVELY in
# that literal -- 'AC', 'CG' and 'ACGT' all pass it -- while rejecting 'AG',
# 'GA' and 'ACC'. Whether a multi-base allele was admitted therefore depended on
# the spelling of the constant, not on the alphabet.
_NT_BASES = frozenset("ACGTU")


def parse_vcf_chrom(vcf_path: str) -> Optional[str]:
    """Extract CHROM from a per-gene VCF (all data rows share same chromosome)."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 1:
                return fields[0]
    return None


@dataclass
class MutationContext:
    """Context for a single mutation."""
    pkey: str
    gene: str
    mutation: str
    nt_pos: int  # 1-based nucleotide position of the FIRST REF base
    ref: str     # REF allele; multi-base for an MNV / deletion / delins
    alt: str     # ALT allele; multi-base for an MNV / insertion / delins
    transcript_seq: str
    wt_rna_window: str
    mut_rna_window: str
    window_center: int  # 0-based index of the first REF base inside the window
    variant_kind: str = "snv"  # Variant.kind: snv | mnv | insertion | deletion | delins
    genomic_pos: Optional[int] = None  # Chromosomal position (1-based) for RBP distance
    substrate: str = "transcript"      # transcript | pre_mRNA -- which molecule was windowed

    @property
    def bed_pos(self) -> Optional[int]:
        """genomic_pos in the 0-based half-open frame used by POSTAR3 BED intervals."""
        return None if self.genomic_pos is None else self.genomic_pos - 1

    @property
    def length_delta(self) -> int:
        """len(ALT) - len(REF); 0 for an SNV or MNV."""
        return len(self.alt) - len(self.ref)

    @property
    def edit_span(self) -> RnaEditSpan:
        """Edit geometry of the PRIMARY RNA window (windows[0])."""
        return RnaEditSpan(self.window_center, len(self.ref), len(self.alt),
                           len(self.wt_rna_window), len(self.mut_rna_window))


def _window_edit_span(context: 'MutationContext',
                      window: Optional[Tuple[str, str, int]] = None) -> RnaEditSpan:
    """RnaEditSpan for one RNA window of a mutation.

    `window` is a (wt_win, mut_win, centre) triple from _generate_windows; with
    it omitted the primary window on the context is used. Both drivers derive
    the projection from RnaEditSpan and nowhere else, so the sites tables they
    emit agree with each other and with every other pipeline's alignment.
    """
    if window is None:
        return context.edit_span
    wt_win, mut_win, centre = window
    return RnaEditSpan(centre, len(context.ref), len(context.alt),
                       len(wt_win), len(mut_win))


@dataclass
class _ParsedResult:
    """Internal: parsed AF3 output with metrics, structures, and aggregation."""
    metrics: Optional[BindingMetrics]
    structures: List[AF3Structure]
    aggregation: Optional[AggregatedBindingAnalysis]


class AlphaFold3Pipeline:
    """
    Main pipeline for AF3 RNA-RBP interaction analysis.
    """

    def __init__(
        self,
        postar_db: str,
        rbp_mapping: str,
        output_dir: str,
        rbp_sequences: Optional[str] = None,
        msa_dir: Optional[str] = None,
        execution_mode: str = "local",
        af3_binary: str = "alphafold3",
        docker_image: str = "alphafold3",
        model_dir: Optional[str] = None,
        window_size: int = 101,
        rbp_window: int = 50,
        validation_log: Optional[str] = None,
        multi_window: bool = False,
        multi_window_offsets: Optional[List[float]] = None,
        max_gpus: Optional[int] = None
    ):
        """
        Initialize pipeline.

        Args:
            postar_db: Path to POSTAR3 database file
            rbp_mapping: Path to gene-UniProt mapping TSV
            output_dir: Output directory
            rbp_sequences: Path to protein sequences FASTA (optional if msa_dir provided)
            msa_dir: Directory containing A3M MSA files (preferred over rbp_sequences)
            execution_mode: 'local', 'batch', or 'cloud'
            af3_binary: Path to AF3 executable
            docker_image: Docker image name for AF3
            model_dir: Path to AF3 model weights directory
            window_size: RNA window size around mutation (odd number)
            rbp_window: Window to search for RBP binding sites (+/-bp)
            validation_log: Optional validation log for filtering mutations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.window_size = window_size
        self.rbp_window = rbp_window
        self.validation_log = validation_log
        self.multi_window = multi_window
        self.multi_window_offsets = multi_window_offsets or [0.3, 0.5, 0.7]

        # Initialize components
        print("Loading POSTAR3 database...", file=sys.stderr)
        self.rbp_db = POSTAR3Database(postar_db)

        print("Loading RBP sequence mapper...", file=sys.stderr)
        self.seq_mapper = RBPSequenceMapper(
            mapping_file=rbp_mapping,
            sequence_fasta=rbp_sequences,
            msa_dir=msa_dir
        )

        # AF3 runner config
        af3_config = AF3RunnerConfig(
            af3_binary=af3_binary,
            output_base_dir=str(self.output_dir / "af3_runs"),
            execution_mode=ExecutionMode(execution_mode),
            docker_image=docker_image,
            model_dir=model_dir,
            max_gpus=max_gpus
        )
        self.af3_runner = AF3Runner(af3_config)

        # Thresholds
        self.threshold_config = ThresholdConfig()

        # Results storage
        self.summary_rows: List[dict] = []
        self.events_rows: List[dict] = []
        self.sites_rows: List[dict] = []

    def process_gene(
        self,
        fasta_path: str,
        mutations_path: Optional[str] = None,
        gene_name: Optional[str] = None,
        chrom: Optional[str] = None,
        tx_start: Optional[int] = None,
        strand: str = "+",
        chrom_mapping: Optional[Dict[str, str]] = None,
        premrna_mapping: Optional[Dict[str, str]] = None
    ):
        """
        Process all mutations for a gene.

        Args:
            fasta_path: Path to ORF FASTA
            mutations_path: Path to mutations CSV (optional when chrom_mapping provided)
            gene_name: Gene symbol (extracted from filename if not provided)
            chrom: Chromosome for POSTAR3 lookup
            tx_start: Transcript start for coordinate conversion
            strand: Strand ('+' or '-')
            chrom_mapping: Dict mapping mutation -> chromosome entry (e.g., C123T -> C87504250T)
            premrna_mapping: Dict mapping mutation -> pre-mRNA entry. Required for
                intronic ('gd.') tokens: an intron has no ORF or transcript
                coordinate, so those variants are scored against the pre_mRNA
                record instead. That is also the correct substrate biologically --
                intronic splicing regulatory elements (branch point,
                polypyrimidine tract, ISE/ISS) are RBP sites that exist only in
                the pre-mRNA, and the 5' splice site spans the exon|intron
                junction, which a bare intron record would truncate.
        """
        # Get gene name
        if gene_name is None:
            gene_name = extract_gene_from_filename(fasta_path)

        print(f"\nProcessing {gene_name}...", file=sys.stderr)

        # Load transcript
        fasta_data = read_fasta(fasta_path)
        if not fasta_data:
            print(f"  Error: Empty FASTA file", file=sys.stderr)
            return

        # Get transcript sequence (prefer 'ORF' or 'transcript' header)
        transcript_seq = None
        for key in ['ORF', 'orf', 'transcript', gene_name]:
            if key in fasta_data:
                transcript_seq = fasta_data[key]
                break
        if transcript_seq is None:
            transcript_seq = next(iter(fasta_data.values()))

        # Second substrate for intronic variants. Absent for a FASTA that
        # predates the pre_mRNA record, in which case intronic tokens fall
        # through to the existing NA path rather than being scored wrongly.
        premrna_seq = fasta_data.get('pre_mRNA')

        # Load mutations
        if mutations_path:
            mutations = trim_muts(mutations_path, self.validation_log, gene_name)
        elif chrom_mapping:
            mutations = list(chrom_mapping.keys())
            if self.validation_log:
                failures = _collect_failures_from_logs(self.validation_log)
                skip_set = failures.get(gene_name.upper(), set()) if gene_name else set()
                mutations = [m for m in mutations if m not in skip_set]
        else:
            print(f"  Error: No mutations source", file=sys.stderr)
            return

        print(f"  {len(mutations)} mutations to process", file=sys.stderr)

        # Process each mutation
        for mutation in mutations:
            try:
                seq_for_mut = transcript_seq
                coord_token = None
                substrate = 'transcript'
                if is_intronic_token(mutation):
                    pre_tok = (premrna_mapping or {}).get(mutation)
                    if premrna_seq is None or pre_tok is None:
                        # No substrate or no coordinate -> the existing NA path.
                        # Scoring an intronic token against the transcript would
                        # read its position as an ORF offset and silently window
                        # the wrong region.
                        self._add_na_mutation(gene_name, mutation)
                        continue
                    seq_for_mut = premrna_seq
                    coord_token = pre_tok
                    substrate = 'pre_mRNA'
                self._process_mutation(
                    gene_name=gene_name,
                    mutation=mutation,
                    transcript_seq=seq_for_mut,
                    chrom=chrom,
                    tx_start=tx_start,
                    strand=strand,
                    chrom_mapping=chrom_mapping,
                    coord_token=coord_token,
                    substrate=substrate
                )
            except Exception as e:
                print(f"  Error processing {mutation}: {e}", file=sys.stderr)
                self._add_failed_mutation(gene_name, mutation, str(e))

    def _process_mutation(
        self,
        gene_name: str,
        mutation: str,
        transcript_seq: str,
        chrom: Optional[str],
        tx_start: Optional[int],
        strand: str,
        chrom_mapping: Optional[Dict[str, str]] = None,
        coord_token: Optional[str] = None,
        substrate: str = 'transcript'
    ):
        """Process a single mutation.

        `mutation` is always the IDENTITY of the variant and is what the pkey is
        built from. `coord_token` overrides which token supplies the coordinate
        and alleles -- for an intronic variant the identity is 'gd.T5000C' while
        the coordinate lives in pre-mRNA space. Leaving them the same object for
        the transcript path keeps that path byte-identical.
        """
        # {GENE}-{sha}. `mutation` is a verbatim token from trim_muts or the
        # chromosome-mapping keys (:278-284); neither normalises.
        pkey = mint_pkey(gene_name, mutation)
        token = coord_token if coord_token is not None else mutation

        # Parse mutation. parse_variant is length-aware and NEVER raises, so an
        # indel / MNV / delins arrives here as a record instead of dying inside
        # get_mutation_data_bioAccurate's int(token[1:-1]). Non-SNV is the
        # DEFAULT path -- there is no flag, because whether a given column
        # survives is decided per column from len(ref) != len(alt), a fact of
        # the record rather than a user preference.
        variant = parse_variant(token, is_nt=True)
        if variant is None:
            # Not a nucleotide token under the shared grammar. Fall back to the
            # ORIGINAL parser, unchanged, so the existing NA-vs-FAILED split for
            # aa tokens, Stop tokens and malformed tokens -- including the exact
            # ValueError text those produce -- is preserved to the byte.
            pos_data = get_mutation_data_bioAccurate(token, is_nt=False)  # AF3 opt-out of primitive nt-validation; local ref-check + alt-guard below
            if pos_data[0] is None:
                self._add_na_mutation(gene_name, mutation)
                return
            nt_pos = pos_data[0]  # 1-based
            wt_nt, mut_nt = pos_data[1]
            # aa/Stop tokens are out of scope for this nucleotide-only pipeline: N/A, not FAILED.
            if not (set(wt_nt.upper()) <= _NT_BASES and set(mut_nt.upper()) <= _NT_BASES):
                self._add_na_mutation(gene_name, mutation)
                return
            # Reachable: int() accepts spellings the grammar does not ('A-12T',
            # 'A 12T', 'A1_2T'). Variant rejects pos < 1, so keep the original
            # out-of-range wording for the one case that used to produce it.
            if nt_pos < 1:
                raise ValueError(f"Position {nt_pos} out of range")
            variant = Variant(pos=nt_pos, ref=wt_nt, alt=mut_nt)

        nt_pos = variant.pos      # 1-based, first REF base
        pos_0 = variant.pos0
        ref, alt = variant.ref, variant.alt

        # Validate mutation. Bound the END of the REF span, not its start: a
        # multi-base REF can begin inside the transcript and run off the end.
        # Identical to `pos_0 >= len(transcript_seq)` when len(ref) == 1.
        if pos_0 < 0 or pos_0 + len(ref) > len(transcript_seq):
            raise ValueError(f"Position {nt_pos} out of range")

        # REF guard over the WHOLE span. Checking transcript_seq[pos_0] alone
        # passes any multi-base REF whose first base happens to match.
        observed = transcript_seq[pos_0:pos_0 + len(ref)]
        if observed.upper() != ref.upper():
            raise ValueError(f"Reference mismatch at {nt_pos}: expected {ref}, found {observed}")

        # Generate RNA windows
        windows = self._generate_windows(transcript_seq, pos_0, ref, alt)

        # Primary window (centered or first offset)
        wt_window, mut_window, window_center = windows[0]

        # The variant must sit inside its own window for the WT/MUT comparison
        # to cover the edit at all. Unreachable for an SNV; reachable for a REF
        # span longer than the window, or one pushed out by 5'-end truncation.
        if not (0 <= window_center and window_center + len(ref) <= len(wt_window)):
            raise ValueError(
                f"REF span outside window: offset {window_center}, "
                f"len(ref) {len(ref)}, window {len(wt_window)}")

        context = MutationContext(
            pkey=pkey,
            gene=gene_name,
            mutation=mutation,
            nt_pos=nt_pos,
            ref=ref,
            alt=alt,
            transcript_seq=transcript_seq,
            wt_rna_window=wt_window,
            mut_rna_window=mut_window,
            window_center=window_center,
            variant_kind=variant.kind,
            substrate=substrate
        )

        # Resolve chromosomal position
        genomic_pos = None
        if chrom_mapping and mutation in chrom_mapping:
            entry = chrom_mapping[mutation]
            # parse_variant, not int(entry[1:-1]). The genomic spelling of a
            # non-SNV carries multi-base alleles and int() raises on every one
            # of them. Widening the token parser above WITHOUT this one would
            # leave genomic_pos None for every indel, so each would reach
            # _add_no_rbps_result and be reported as 'no_rbps_in_region' -- a
            # fabricated negative, strictly worse than today's FAILED row.
            genomic_variant = parse_variant(entry, is_nt=True)
            if genomic_variant is not None:
                genomic_pos = genomic_variant.pos
            else:
                print(f"    Warning: Could not parse chromosome mapping for {mutation}: {entry}", file=sys.stderr)
        elif tx_start is not None:
            if strand == '+':
                genomic_pos = tx_start + nt_pos - 1
            else:
                genomic_pos = tx_start - nt_pos + 1

        context.genomic_pos = genomic_pos

        # Query RBPs across the variant's WHOLE REF span, not just its first base
        rbp_sites = self._get_nearby_rbps(chrom, context.bed_pos, len(ref))

        if not rbp_sites:
            # "We looked and found nothing" and "we never looked" are different
            # findings. Reporting the second as the first is a fabricated
            # observation, in the same class as a coalesced 0.0.
            if chrom is None:
                reason = 'NA:no_chromosome'
            elif context.bed_pos is None:
                reason = 'NA:unresolved_genomic_position'
            else:
                reason = 'no_rbps_in_region'
            self._add_no_rbps_result(context, reason)
            return

        # Group by RBP
        rbps_to_test = self.rbp_db.group_by_rbp(rbp_sites)
        print(f"    {pkey}: {len(rbps_to_test)} RBPs to test", file=sys.stderr)

        # Phase 1: Submit all RBP jobs (non-blocking)
        pending_list = []
        skipped_rbps: List[str] = []
        for rbp_name, sites in rbps_to_test.items():
            pending = self._submit_rbp_jobs(
                context, rbp_name, sites,
                windows=windows if len(windows) > 1 else None,
                skipped=skipped_rbps
            )
            if pending:
                pending_list.append(pending)

        # Phase 2: Collect all results (blocks on futures as they complete)
        delta_list = []
        for pending in pending_list:
            delta = self._collect_rbp_results(context, pending)
            if delta:
                delta_list.append(delta)

        # Aggregate results
        self._finalize_mutation_results(context, delta_list, skipped_rbps)

    def _get_nearby_rbps(
        self,
        chrom: Optional[str],
        bed_pos: Optional[int],
        span_len: int = 1
    ) -> List[RBPBindingSite]:
        """Query POSTAR3 for RBPs near the variant's REF span (0-based BED frame).

        The neighbourhood is [bed_pos - rbp_window, bed_pos + span_len + rbp_window).
        With span_len == 1 that interval is exactly what
        query_position(chrom, bed_pos, rbp_window) produces, so SNV behaviour is
        unchanged; a multi-base REF widens it by its own length rather than
        anchoring the whole neighbourhood on the span's 5' base.
        """
        if chrom is None or bed_pos is None:
            return []

        return self.rbp_db.query(chrom,
                                 bed_pos - self.rbp_window,
                                 bed_pos + span_len + self.rbp_window)

    def _generate_windows(
        self,
        transcript_seq: str,
        pos_0: int,
        ref: str,
        alt: str
    ) -> List[Tuple[str, str, int]]:
        """
        Generate RNA windows around the mutation.

        Returns list of (wt_window, mut_window, window_center) tuples, where
        window_center is the 0-based index of the FIRST REF base inside the
        window. Single window when multi_window is disabled.

        Two length-aware properties, both no-ops for an SNV:

        1. The window is centred on the MIDPOINT of the REF span. Anchoring on
           its first base gives a multi-base REF the full 5' flank but only
           (half - len(ref)) of 3' flank, and at len(ref) > half the window
           stops covering the edited span at all. len(ref)//2 == 0 for an SNV.
        2. The mutant window covers the SAME BIOLOGICAL SPAN as the WT one, so
           it is longer/shorter by exactly (len(alt) - len(ref)). Taking a
           second window centred on the mutant sequence instead would slide it
           (alt - ref) bases 3', making WT vs MUT a comparison of two different
           regions of the transcript.
        """
        # splice_seq honours len(ref); the previous one-base concatenation
        # (`seq[:pos_0] + mut_nt + seq[pos_0+1:]`) cannot express an indel.
        # REF was verified against the transcript by the caller.
        mut_seq = splice_seq(transcript_seq, pos_0, ref, alt, validate=False)
        delta = len(alt) - len(ref)
        centre = min(pos_0 + len(ref) // 2, len(transcript_seq) - 1)

        if not self.multi_window:
            wt_window = subseq(transcript_seq, centre, self.window_size)
            window_start = max(0, centre - self.window_size // 2)
            mut_window = mut_seq[window_start:window_start + len(wt_window) + delta]
            return [(wt_window, mut_window, pos_0 - window_start)]

        seen = set()
        windows = []
        for frac in self.multi_window_offsets:
            target_center = int(frac * self.window_size)
            window_start = centre - target_center
            window_start = max(0, min(window_start, len(transcript_seq) - self.window_size))
            window_end = window_start + self.window_size

            wt_win = transcript_seq[window_start:window_end]
            mut_win = mut_seq[window_start:window_end + delta]

            if wt_win not in seen:
                seen.add(wt_win)
                windows.append((wt_win, mut_win, pos_0 - window_start))

        # No `if windows else <fallback>`: self.multi_window_offsets is
        # `offsets or [0.3, 0.5, 0.7]` in __init__, so it is never empty, the
        # first iteration always appends, and the fallback could never run. It
        # also carried its own copy of the single-window arithmetic in the old
        # one-base form, which would have been wrong for every indel.
        return windows

    @dataclass
    class _PendingRBPAnalysis:
        """Tracks submitted async jobs for one RBP."""
        rbp_name: str
        sites: List
        distance: int
        wt_future: Optional[Future] = None
        mut_future: Optional[Future] = None
        window_wt_futures: Optional[List[Future]] = None
        window_mut_futures: Optional[List[Future]] = None
        n_windows: int = 1
        # The window triples these futures were built from. Needed at collection
        # time to project each window's sites back into its own WT frame; the
        # futures alone do not say which window they came from.
        windows: Optional[List[Tuple[str, str, int]]] = None

    def _submit_rbp_jobs(
        self,
        context: MutationContext,
        rbp_name: str,
        sites: List[RBPBindingSite],
        windows: Optional[List[Tuple[str, str, int]]] = None,
        skipped: Optional[List[str]] = None
    ) -> Optional['AlphaFold3Pipeline._PendingRBPAnalysis']:
        """Submit AF3 jobs for one RBP without blocking. Returns pending tracker.

        An RBP that is not submitted is appended to `skipped` by name and
        reason. It is otherwise absent from delta_list AND from n_rbps_tested,
        i.e. indistinguishable from an RBP that was never near the variant.
        """
        def _skip(reason: str) -> None:
            if skipped is not None:
                skipped.append(f"{rbp_name}:{reason}")

        rbp_data = self.seq_mapper.get_rbp_data(rbp_name)
        if not rbp_data:
            print(f"      {rbp_name}: sequence not found", file=sys.stderr)
            _skip('sequence_not_found')
            return None

        protein_seq = rbp_data.sequence
        protein_msa = rbp_data.msa_content

        # Test the cap against the LONGER of the two windows: an insertion makes
        # the mutant window longer than the WT one, so sizing on the WT alone
        # would submit a MUT job that exceeds the limit.
        total_tokens = max(len(context.wt_rna_window),
                           len(context.mut_rna_window)) + len(protein_seq)
        if total_tokens > 5000:
            print(f"      {rbp_name}: token limit exceeded ({total_tokens})", file=sys.stderr)
            _skip(f'token_limit_exceeded_{total_tokens}')
            return None

        # Distance from the whole REF span, not from its 5' base: a deletion
        # that straddles a binding site has neither endpoint inside it.
        distance = (min(site.distance_to_span(context.bed_pos, len(context.ref))
                        for site in sites)
                    if context.bed_pos is not None else 0)

        pending = self._PendingRBPAnalysis(
            rbp_name=rbp_name,
            sites=sites,
            distance=distance,
            windows=windows
        )

        if windows and len(windows) > 1:
            pending.n_windows = len(windows)
            pending.window_wt_futures = []
            pending.window_mut_futures = []
            for i, (wt_win, mut_win, _center) in enumerate(windows):
                wt_in = create_rna_protein_input(
                    job_name=f"{context.pkey}_{rbp_name}_WT_w{i}",
                    rna_seq=wt_win, protein_seq=protein_seq, protein_msa=protein_msa
                )
                mut_in = create_rna_protein_input(
                    job_name=f"{context.pkey}_{rbp_name}_MUT_w{i}",
                    rna_seq=mut_win, protein_seq=protein_seq, protein_msa=protein_msa
                )
                pending.window_wt_futures.append(
                    self.af3_runner.submit_job_async(wt_in, job_id=f"{context.pkey}_{rbp_name}_WT_w{i}")
                )
                pending.window_mut_futures.append(
                    self.af3_runner.submit_job_async(mut_in, job_id=f"{context.pkey}_{rbp_name}_MUT_w{i}")
                )
        else:
            wt_input = create_rna_protein_input(
                job_name=f"{context.pkey}_{rbp_name}_WT",
                rna_seq=context.wt_rna_window, protein_seq=protein_seq, protein_msa=protein_msa
            )
            mut_input = create_rna_protein_input(
                job_name=f"{context.pkey}_{rbp_name}_MUT",
                rna_seq=context.mut_rna_window, protein_seq=protein_seq, protein_msa=protein_msa
            )
            pending.wt_future = self.af3_runner.submit_job_async(
                wt_input, job_id=f"{context.pkey}_{rbp_name}_WT"
            )
            pending.mut_future = self.af3_runner.submit_job_async(
                mut_input, job_id=f"{context.pkey}_{rbp_name}_MUT"
            )

        return pending

    def _collect_rbp_results(
        self,
        context: MutationContext,
        pending: '_PendingRBPAnalysis'
    ) -> Optional[DeltaMetrics]:
        """Block on futures, parse results, compute delta metrics for one RBP."""
        rbp_name = pending.rbp_name
        distance = pending.distance

        if pending.n_windows > 1:
            # Keyed by window index, not appended to a flat list: an incomplete
            # job used to shift every later result's position, and the sites
            # block below needs to know WHICH window a result came from to
            # project it back into that window's own WT frame.
            wt_by_win = {}
            mut_by_win = {}
            for i, (wt_f, mut_f) in enumerate(zip(pending.window_wt_futures,
                                                  pending.window_mut_futures)):
                wt_job = wt_f.result()
                mut_job = mut_f.result()
                if wt_job.status == "completed" and wt_job.result_path:
                    wt_by_win[i] = self._parse_af3_output(wt_job.result_path, rbp_name)
                if mut_job.status == "completed" and mut_job.result_path:
                    mut_by_win[i] = self._parse_af3_output(mut_job.result_path, rbp_name)

            wt_metrics = self._aggregate_parsed_results(list(wt_by_win.values()), rbp_name)
            mut_metrics = self._aggregate_parsed_results(list(mut_by_win.values()), rbp_name)

            for allele, by_win in (('WT', wt_by_win), ('MUT', mut_by_win)):
                for i in sorted(by_win):
                    r = by_win[i]
                    if r and r.structures:
                        sites_data = extract_interface_sites(r.structures[0])
                        freq_rna = r.aggregation.contact_frequency_rna if r.aggregation else None
                        freq_prot = r.aggregation.contact_frequency_protein if r.aggregation else None
                        window = pending.windows[i] if pending.windows else None
                        self.sites_rows.extend(format_sites_rows(
                            context.pkey, rbp_name, allele, sites_data,
                            freq_rna, freq_prot,
                            edit_span=_window_edit_span(context, window)))
                        break

            delta = compute_delta_metrics(
                rbp_name=rbp_name,
                wt_metrics=wt_metrics,
                mut_metrics=mut_metrics,
                distance_to_mutation=distance,
                config=self.threshold_config
            )
            delta.n_windows = pending.n_windows
            return delta

        # Single-window
        wt_job = pending.wt_future.result()
        mut_job = pending.mut_future.result()

        wt_parsed = self._parse_af3_output(wt_job.result_path, rbp_name) if wt_job.status == "completed" else None
        mut_parsed = self._parse_af3_output(mut_job.result_path, rbp_name) if mut_job.status == "completed" else None

        wt_metrics = wt_parsed.metrics if wt_parsed else None
        mut_metrics = mut_parsed.metrics if mut_parsed else None

        edit_span = _window_edit_span(context)
        if wt_parsed and wt_parsed.structures:
            sites_data = extract_interface_sites(wt_parsed.structures[0])
            freq_rna = wt_parsed.aggregation.contact_frequency_rna if wt_parsed.aggregation else None
            freq_prot = wt_parsed.aggregation.contact_frequency_protein if wt_parsed.aggregation else None
            self.sites_rows.extend(format_sites_rows(
                context.pkey, rbp_name, 'WT', sites_data, freq_rna, freq_prot,
                edit_span=edit_span))
        if mut_parsed and mut_parsed.structures:
            sites_data = extract_interface_sites(mut_parsed.structures[0])
            freq_rna = mut_parsed.aggregation.contact_frequency_rna if mut_parsed.aggregation else None
            freq_prot = mut_parsed.aggregation.contact_frequency_protein if mut_parsed.aggregation else None
            self.sites_rows.extend(format_sites_rows(
                context.pkey, rbp_name, 'MUT', sites_data, freq_rna, freq_prot,
                edit_span=edit_span))

        return compute_delta_metrics(
            rbp_name=rbp_name,
            wt_metrics=wt_metrics,
            mut_metrics=mut_metrics,
            distance_to_mutation=distance,
            config=self.threshold_config
        )

    def _parse_af3_output(
        self,
        output_dir: Path,
        rbp_name: str
    ) -> Optional[_ParsedResult]:
        """Parse AF3 output with ensemble sample aggregation."""
        if output_dir is None:
            return None

        structures = parse_all_samples(str(output_dir))
        if not structures:
            return None

        analyses = [analyze_binding(s, rna_chain="R", protein_chain="P") for s in structures]

        if len(structures) == 1:
            binding = analyses[0]
            if not binding:
                return _ParsedResult(metrics=None, structures=structures, aggregation=None)
            metrics = BindingMetrics(
                rbp_name=rbp_name,
                chain_pair_pae_min=binding.chain_pair_pae_min,
                interface_contacts=binding.n_contacts,
                interface_plddt_rna=binding.interface_plddt_rna,
                interface_plddt_protein=binding.interface_plddt_protein,
                has_binding=binding.n_contacts >= self.threshold_config.min_contacts
            )
            return _ParsedResult(metrics=metrics, structures=structures, aggregation=None)

        # Multi-sample aggregation
        agg = aggregate_binding_analyses(analyses)
        if not agg:
            return _ParsedResult(metrics=None, structures=structures, aggregation=None)

        metrics = BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=agg.mean.chain_pair_pae_min,
            interface_contacts=agg.mean.n_contacts,
            interface_plddt_rna=agg.mean.interface_plddt_rna,
            interface_plddt_protein=agg.mean.interface_plddt_protein,
            has_binding=agg.mean.n_contacts >= self.threshold_config.min_contacts,
            n_samples=agg.n_samples,
            std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
            std_interface_contacts=agg.std_n_contacts,
            std_plddt_rna=agg.std_interface_plddt_rna,
            std_plddt_protein=agg.std_interface_plddt_protein
        )
        return _ParsedResult(metrics=metrics, structures=structures, aggregation=agg)

    def _aggregate_parsed_results(
        self,
        results: List[Optional[_ParsedResult]],
        rbp_name: str
    ) -> Optional[BindingMetrics]:
        """Aggregate metrics across multiple windows."""
        all_analyses = []
        for r in results:
            if r and r.metrics:
                # Build a BindingAnalysis-like from each result for aggregation
                all_analyses.append(BindingAnalysis(
                    rna_chain="R", protein_chain="P",
                    n_contacts=r.metrics.interface_contacts,
                    min_contact_distance=0.0,
                    mean_contact_distance=0.0,
                    interface_plddt_rna=r.metrics.interface_plddt_rna,
                    interface_plddt_protein=r.metrics.interface_plddt_protein,
                    chain_pair_pae_min=r.metrics.chain_pair_pae_min,
                    contact_residues_rna=[], contact_residues_protein=[]
                ))

        if not all_analyses:
            return None

        agg = aggregate_binding_analyses(all_analyses)
        if not agg:
            return None

        return BindingMetrics(
            rbp_name=rbp_name,
            chain_pair_pae_min=agg.mean.chain_pair_pae_min,
            interface_contacts=agg.mean.n_contacts,
            interface_plddt_rna=agg.mean.interface_plddt_rna,
            interface_plddt_protein=agg.mean.interface_plddt_protein,
            has_binding=agg.mean.n_contacts >= self.threshold_config.min_contacts,
            n_samples=agg.n_samples,
            std_chain_pair_pae_min=agg.std_chain_pair_pae_min,
            std_interface_contacts=agg.std_n_contacts,
            std_plddt_rna=agg.std_interface_plddt_rna,
            std_plddt_protein=agg.std_interface_plddt_protein
        )

    @staticmethod
    def _variant_columns(context: MutationContext, skipped_rbps: List[str]) -> dict:
        """Variant-class columns carried by every row built from a MutationContext.

        Every AF3 metric other than sites.res_id is a whole-complex scalar, so
        these are the only columns the variant class itself determines. They are
        populated for EVERY class -- an SNV reports kind 'snv', delta 0 and a
        fully aligned window, which are all true statements, not placeholders.
        """
        return {
            'variant_class': context.variant_kind,
            'length_delta_nt': context.length_delta,
            'substrate': context.substrate,
            'rna_window_alignment': context.edit_span.alignment,
            'n_rbps_skipped': len(skipped_rbps),
            'rbps_skipped': ';'.join(skipped_rbps),
        }

    def _finalize_mutation_results(
        self,
        context: MutationContext,
        delta_list: List[DeltaMetrics],
        skipped_rbps: Optional[List[str]] = None
    ):
        """Aggregate and store results for a mutation."""
        # Summary row
        summary = aggregate_mutation_summary(delta_list)
        summary['pkey'] = context.pkey
        summary['Gene'] = context.gene

        # QC flags: check whether AF3 predictions actually succeeded
        has_complete = any(d.wt_metrics is not None and d.mut_metrics is not None for d in delta_list)
        has_partial = any(d.wt_metrics is not None or d.mut_metrics is not None for d in delta_list)
        if not delta_list:
            summary['qc_flags'] = 'no_rbps_tested'
        elif has_complete:
            summary['qc_flags'] = 'PASS'
        elif has_partial:
            summary['qc_flags'] = 'PARTIAL'
        else:
            summary['qc_flags'] = 'ALL_FAILED'

        summary.update(self._variant_columns(context, skipped_rbps or []))
        self.summary_rows.append(summary)

        # Events rows
        events = format_events_rows(context.pkey, delta_list)
        self.events_rows.extend(events)

    def _add_no_rbps_result(self, context: MutationContext,
                            qc_flag: str = 'no_rbps_in_region'):
        """Add result for a mutation with no RBP comparison to report.

        The literal dict is deliberate: aggregate_mutation_summary([]) rounds
        the two delta columns, which would render them '0.0' where this row has
        always written '0'.
        """
        row = {
            'pkey': context.pkey,
            'Gene': context.gene,
            'n_rbps_tested': 0,
            'n_rbps_binding_wt': 0,
            'n_rbps_binding_mut': 0,
            'global_count_gained': 0,
            'global_count_lost': 0,
            'global_count_strengthened': 0,
            'global_count_weakened': 0,
            'global_max_abs_delta_pae': 0,
            'top_event_rbp': '',
            'top_event_class': 'none',
            'top_event_delta_pae': 0,
            'qc_flags': qc_flag
        }
        row.update(self._variant_columns(context, []))
        self.summary_rows.append(row)

    def _summary_stub(self, gene: str, mutation: str, qc_flag: str) -> dict:
        """Summary row for a mutation that never reached a window.

        aggregate_mutation_summary([]) supplies the full column set; write_tsv
        takes its header from rows[0], so a short row would drop columns.

        variant_class is filled whenever the token parses, because the class is
        a property of the TOKEN and stays true however the token was rejected.
        The window columns stay empty: no window was built, and qc_flag already
        names why.
        """
        summary = aggregate_mutation_summary([])
        summary['pkey'] = mint_pkey(gene, mutation)
        summary['Gene'] = gene
        summary['qc_flags'] = qc_flag
        variant = parse_variant(mutation, is_nt=True)
        if variant is not None:
            summary['variant_class'] = variant.kind
            summary['length_delta_nt'] = variant.length_delta
        return summary

    def _add_na_mutation(self, gene: str, mutation: str):
        """Add result for a token this nucleotide-only pipeline cannot score."""
        self.summary_rows.append(
            self._summary_stub(gene, mutation, 'NA:non_nucleotide_token')
        )

    def _add_failed_mutation(self, gene: str, mutation: str, error: str):
        """Add result for failed mutation."""
        self.summary_rows.append(
            self._summary_stub(gene, mutation, f'FAILED:{error[:50]}')
        )

    def _write_rows(self, rows, path):
        """Write a list of row dicts to a TSV file.

        Header is the UNION of every row's keys, in first-appearance order.
        write_tsv defaults to fieldnames=rows[0].keys() with extrasaction='ignore',
        so a first row narrower than a later one silently truncates the whole file
        and neither direction errors. That fires the moment any row carries a key
        the first row lacks -- a suppression reason code, an optional metric.
        _summary_stub works around it for the summary table only; events and sites
        had no protection at all.
        """
        if rows:
            fieldnames = list(dict.fromkeys(key for row in rows for key in row))
            write_tsv(rows, path, fieldnames)
            print(f"Wrote {len(rows)} rows to {path}", file=sys.stderr)

    def flush_gene(self, gene_name: str):
        """Write accumulated rows for one gene to nested per-gene/tool paths and clear them."""
        out_dir = self.output_dir / gene_name / "AlphaFold3"
        out_dir.mkdir(parents=True, exist_ok=True)
        self._write_rows(self.summary_rows, out_dir / f"{gene_name}.tsv")
        self._write_rows(self.events_rows,  out_dir / f"{gene_name}.events.tsv")
        self._write_rows(self.sites_rows,   out_dir / f"{gene_name}.sites.tsv")
        self.summary_rows.clear()
        self.events_rows.clear()
        self.sites_rows.clear()


def main():
    parser = argparse.ArgumentParser(
        description='AlphaFold3 RNA-RBP Interaction Pipeline'
    )

    # Required inputs
    parser.add_argument('-pd', '--postar-db', required=True,
                       help='POSTAR3 database file')
    parser.add_argument('-rm', '--rbp-mapping', required=True,
                       help='Gene-UniProt mapping TSV')
    parser.add_argument('-rs', '--rbp-sequences',
                       help='Protein sequences FASTA (optional if --msa-dir provided)')
    parser.add_argument('-md', '--msa-dir',
                       help='Directory with A3M MSA files (preferred over --rbp-sequences)')

    # Gene inputs (file or directory, auto-detected)
    parser.add_argument('-f', '--fasta',
                        help='Transcript FASTA file or directory of FASTA files')
    parser.add_argument('-mu', '--mutations',
                        help='Mutations CSV file or directory (optional when --chromosome-mapping provided)')

    # VCF-based coordinate resolution (replaces --chrom/--tx-start/--strand)
    parser.add_argument('-v', '--vcf',
                        help='Per-gene VCF file or directory from vcf_converter.py (provides chromosome)')
    parser.add_argument('-cm', '--chromosome-mapping',
                        help='Chromosome mapping CSV file or directory (provides mutations and chromosomal positions)')
    parser.add_argument('-pm', '--premrna-mapping',
                        help='pre-mRNA mapping CSV file or directory (premrna_mapping_<GENE>.csv). Required to score intronic (gd.) variants: they have no ORF or transcript coordinate and are windowed on the pre_mRNA record.')

    # Legacy genomic coordinates (still supported)
    parser.add_argument('-ch', '--chrom', help='Chromosome (alternative to --vcf)')
    parser.add_argument('-ts', '--tx-start', type=int, help='Transcript start position (alternative to --chromosome-mapping)')
    parser.add_argument('-s', '--strand', default='+', help='Strand (+/-)')

    # Output
    parser.add_argument('--output', '-o', required=True,
                       help='Output base directory')

    # Execution
    parser.add_argument('-em', '--execution-mode', default='local',
                       choices=['local'],
                       help='AF3 execution mode (only local is supported here; '
                            'for SLURM use `python -m biofeaturefactory.alphafold3.burst submit`)')
    parser.add_argument('-ab', '--af3-binary', default='alphafold3',
                       help='Path to AF3 executable')
    parser.add_argument('-di', '--docker-image', default='alphafold3',
                       help='Docker image name for AF3')
    parser.add_argument('-mdi', '--model-dir',
                       help='Path to AF3 model weights directory')

    # Parameters
    parser.add_argument('-ws', '--window-size', type=int, default=101,
                       help='RNA window size (odd number)')
    parser.add_argument('-rw', '--rbp-window', type=int, default=50,
                       help='Window for RBP site lookup (+/-bp)')
    parser.add_argument('-vl', '--validation-log',
                       help='Validation log for filtering mutations')

    # Multi-window mode (optional, multiplies AF3 runs)
    parser.add_argument('-mw', '--multi-window', action='store_true', default=False,
                       help='Run multiple windows per mutation (multiplies AF3 runs)')
    parser.add_argument('-mwo', '--multi-window-offsets', type=str, default='0.3,0.5,0.7',
                       help='Mutation position as fraction of window (default: 0.3,0.5,0.7)')
    parser.add_argument('-mg', '--max-gpus', type=int, default=None,
                       help='Max GPUs for parallel AF3 execution (default: auto-detect)')

    args = parser.parse_args()


    # One root supplies these: <root>/<GENE>/mappings/{mutations,chromosome,

    # intron_premRNA}/ sit beside <root>/<GENE>/fastas/. --rbp-mapping is NOT

    # derived -- it is a POSTAR3 RBP table, not a variant_mapping product.

    args.mutations = derive_mutations_root(args.mutations, args.fasta, label="af3")

    args.chromosome_mapping = derive_mapping_root(

        args.chromosome_mapping, args.fasta, "chromosome", label="af3")

    args.premrna_mapping = derive_mapping_root(

        args.premrna_mapping, args.fasta, "premrna", label="af3")

    # Validate that we have either MSA dir or sequences
    if not args.msa_dir and not args.rbp_sequences:
        parser.error("Provide either --msa-dir or --rbp-sequences")

    # Model weights required for local execution
    if args.execution_mode == 'local' and not args.model_dir:
        parser.error("--model-dir is required for local execution mode")

    # Initialize pipeline
    pipeline = AlphaFold3Pipeline(
        postar_db=args.postar_db,
        rbp_mapping=args.rbp_mapping,
        output_dir=args.output,
        rbp_sequences=args.rbp_sequences,
        msa_dir=args.msa_dir,
        execution_mode=args.execution_mode,
        af3_binary=args.af3_binary,
        docker_image=args.docker_image,
        model_dir=args.model_dir,
        window_size=args.window_size,
        rbp_window=args.rbp_window,
        validation_log=args.validation_log,
        multi_window=args.multi_window,
        multi_window_offsets=[float(x) for x in args.multi_window_offsets.split(',')] if args.multi_window_offsets else None,
        max_gpus=args.max_gpus
    )

    # Resolve inputs
    fasta_input = Path(args.fasta) if args.fasta else None
    mutations_input = Path(args.mutations) if args.mutations else None
    vcf_input = Path(args.vcf) if args.vcf else None
    chrom_map_input = Path(args.chromosome_mapping) if args.chromosome_mapping else None
    premrna_map_input = Path(args.premrna_mapping) if args.premrna_mapping else None

    def _load_premrna_map(gene):
        """pre-mRNA mapping for one gene, or None when not supplied."""
        if premrna_map_input is None:
            return None
        if premrna_map_input.is_file():
            return load_mapping(str(premrna_map_input), mapType="pre_mRNA")
        cands = list(premrna_map_input.glob(f'*{gene}*.csv')) if gene else []
        return load_mapping(str(cands[0]), mapType="pre_mRNA") if cands else None

    if not fasta_input:
        parser.error("--fasta is required")
    if not mutations_input and not chrom_map_input:
        parser.error("Provide --mutations or --chromosome-mapping")

    # Chromosome name is always required for POSTAR3 lookup — chromosome-mapping only provides
    # genomic positions, not the chromosome name itself
    if not args.chrom and not vcf_input:
        parser.error("--chrom or --vcf is required to resolve chromosome name for POSTAR3 lookup")

    # When using --mutations without --chromosome-mapping, also require --tx-start
    # to convert transcript positions to genomic coordinates
    if mutations_input and not chrom_map_input:
        if args.tx_start is None:
            parser.error("--tx-start is required when using --mutations without --chromosome-mapping")

    if fasta_input.is_dir():
        # --- Directory mode ---
        for fasta_file in sorted(fasta_input.glob('*.fasta')):
            gene_name = extract_gene_from_filename(str(fasta_file))

            # Find matching mutations file
            mut_path = None
            if mutations_input and mutations_input.is_dir():
                candidates = list(mutations_input.glob(f'*{gene_name}*.csv'))
                if candidates:
                    mut_path = str(candidates[0])

            # Resolve chromosome from VCF
            chrom = args.chrom
            if vcf_input:
                vcf_search = vcf_input if vcf_input.is_file() else None
                if vcf_input.is_dir():
                    vcf_candidates = list(vcf_input.glob(f'{gene_name}.vcf'))
                    if vcf_candidates:
                        vcf_search = vcf_candidates[0]
                if vcf_search:
                    chrom = parse_vcf_chrom(str(vcf_search))

            # Load chromosome mapping
            chrom_mapping = None
            if chrom_map_input:
                if chrom_map_input.is_file():
                    chrom_mapping = load_mapping(str(chrom_map_input), mapType="chromosome")
                elif chrom_map_input.is_dir():
                    map_candidates = list(chrom_map_input.glob(f'*{gene_name}*.csv'))
                    if map_candidates:
                        chrom_mapping = load_mapping(str(map_candidates[0]), mapType="chromosome")

            if not mut_path and not chrom_mapping:
                print(f"No mutations source for {gene_name}", file=sys.stderr)
                continue

            pipeline.process_gene(
                fasta_path=str(fasta_file),
                mutations_path=mut_path,
                chrom=chrom,
                tx_start=args.tx_start,
                strand=args.strand,
                chrom_mapping=chrom_mapping,
                premrna_mapping=_load_premrna_map(gene_name)
            )
            pipeline.flush_gene(gene_name)

    else:
        # --- Single file mode ---
        chrom = args.chrom
        chrom_mapping = None

        if vcf_input:
            chrom = parse_vcf_chrom(str(vcf_input))

        if chrom_map_input:
            if chrom_map_input.is_file():
                chrom_mapping = load_mapping(str(chrom_map_input), mapType="chromosome")
            elif chrom_map_input.is_dir():
                gene_name = extract_gene_from_filename(str(fasta_input))
                map_candidates = list(chrom_map_input.glob(f'*{gene_name}*.csv'))
                if map_candidates:
                    chrom_mapping = load_mapping(str(map_candidates[0]), mapType="chromosome")

        mut_path = str(mutations_input) if mutations_input else None
        gene_name = extract_gene_from_filename(str(fasta_input))

        pipeline.process_gene(
            fasta_path=str(fasta_input),
            mutations_path=mut_path,
            chrom=chrom,
            tx_start=args.tx_start,
            strand=args.strand,
            chrom_mapping=chrom_mapping,
            premrna_mapping=_load_premrna_map(gene_name)
        )
        pipeline.flush_gene(gene_name)

    pipeline.af3_runner.shutdown()


if __name__ == '__main__':
    main()
