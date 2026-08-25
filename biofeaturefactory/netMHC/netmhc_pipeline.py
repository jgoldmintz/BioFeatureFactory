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
NetMHC Pipeline for MHC Binding Prediction

Predicts MHC class I binding peptides for WT and mutant protein sequences and
compares the two. Only netMHC-4.0 is supported, invoked natively.

netMHC scores every peptide window, so a length-L sequence yields exactly
L - l + 1 rows per allele, non-binders included. Zero rows is therefore a tool
failure, not a negative result; the usual cause is a missing data/ directory,
which leaves netMHC exiting 0 with no predictions. The pipeline echoes netMHC's
raw output and sets a qc_flags code when that happens.

Key features:
- MHC class I binding prediction across multiple HLA alleles
- WT vs mutant comparison with delta scores
- Three-tier output: summary, events, sites
- Integration with mutation mapping CSVs
"""

import csv
import hashlib
import json
import math
import os
import subprocess
import tempfile
import shutil
from pathlib import Path
import sys
import platform
from concurrent.futures import ThreadPoolExecutor


# Import utility functions. The removed names (combine_batch_outputs,
# discover_fasta_files, ExtractGeneFromFASTA, parse_predictions_with_mutation_filtering,
# load_wt_sequences, extract_mutation_from_sequence_name) plus `time`, `logging` and
# the typing aliases were imported but never referenced anywhere in this module.
from biofeaturefactory.lib.utility import (
    derive_mutations_root,
    discover_mutation_files,
    mint_pkey,
    token_from_name,
    write_fasta,
    write_tsv,
    load_validation_failures,
    should_skip_mutation,
    load_wt_sequence_map,
    extract_gene_from_filename,
    translate_orf_sequence,
    build_mutant_sequences_for_gene,
    detect_alphabet,
    trim_muts,
    parse_variant,
    protein_consequence,
    infer_edit_span,
    align_wt_to_mut,
    splice_seq,
    is_intronic_token,
)


# write_tsv is called with these lists and extrasaction='ignore', so a row key absent
# here is dropped silently -- add the column here whenever a row dict gains one.
#
# The first 17 event columns and the first 18 summary columns are the pre-existing
# layout, in their original order. Everything non-SNV support needed is APPENDED,
# so a consumer reading the old columns positionally is unaffected.
EVENT_FIELDNAMES = [
    'pkey', 'Gene', 'mutation', 'wt_peptide', 'mut_peptide', 'pos', 'mhc_allele',
    'wt_rank', 'mut_rank', 'delta_rank',
    'wt_affinity', 'mut_affinity', 'delta_affinity',
    'bind_level_wt', 'bind_level_mut',
    'classification', 'classification_code',
    # --- appended for non-SNV support ---
    'variant_class', 'aa_consequence', 'mut_pos', 'align_status', 'qc_flags',
]

SUMMARY_FIELDNAMES = [
    'pkey', 'Gene', 'mutation',
    'n_epitopes_wt', 'n_epitopes_mut',
    'count_gained', 'count_lost', 'count_strengthened', 'count_weakened', 'count_stable',
    'max_abs_delta_rank', 'sum_abs_delta_rank',
    'top_event_type', 'top_event_allele', 'top_event_wt_peptide', 'top_event_mut_peptide',
    'top_event_delta_rank',
    'qc_flags',
    # --- appended for non-SNV support ---
    'variant_class', 'aa_consequence',
    'count_wt_only', 'count_mut_only',
    'n_registers_aligned', 'n_registers_total',
]

SITE_FIELDNAMES = [
    'pkey', 'Gene', 'sequence_type',
    'pos', 'mhc_allele', 'peptide', 'core',
    'affinity', 'rank', 'bind_level', 'identity',
]

# The 20 standard residues. parse_variant's amino-acid grammar accepts any run of
# letters as an allele, so the alphabet has to be checked separately before an
# aa-level allele is spliced into a protein -- see aa_allele_problem.
AA_RESIDUES = frozenset('ACDEFGHIKLMNPQRSTVWY')


def is_linux_host():
    """Return True when running on a Linux kernel."""
    return platform.system().lower() == "linux"


def resolve_native_netmhc_path(user_path=None, tool_version="netMHC"):
    """
    Resolve a usable native NetMHC executable when available.

    Search order:
    1. Explicit --native-netmhc-path value
    2. $NETMHC_PATH
    3. $NETMHC_HOME/<tool_version>
    4. Common install locations (~/netMHC/, /opt/netMHC/, /usr/local/bin/)

    $NETMHCPAN_PATH is deliberately NOT consulted -- see the comment on
    common_roots below.

    Args:
        user_path: User-specified path to NetMHC binary
        tool_version: Which NetMHC tool to use. Only netMHC-4.0 is supported --
            the invocation and the 15-column/'HLA'-header parser are hardcoded to
            it (see --netmhc-tool, choices=['netMHC']). The default is 'netMHC' so
            a programmatic caller that omits it cannot silently resolve netMHCpan,
            whose output this parser reads as zero predictions (F5).
    """
    candidates = []

    def _add(path):
        if path:
            candidates.append(os.path.expanduser(path))

    _add(user_path)
    _add(os.environ.get("NETMHC_PATH"))

    netmhc_home = os.environ.get("NETMHC_HOME")
    if netmhc_home:
        _add(os.path.join(netmhc_home, tool_version))

    home = Path.home()
    # netMHC-4.0 only (see --netmhc-tool): do NOT auto-discover netMHCpan -- its
    # output format is unreadable by this parser and would silently yield zero
    # predictions.
    common_roots = [
        home / "netMHC" / tool_version,
        Path(f"/opt/netMHC/{tool_version}"),
        Path(f"/usr/local/bin/{tool_version}"),
    ]

    for candidate in common_roots:
        _add(str(candidate))

    for path in candidates:
        if not path:
            continue
        # If given a directory, look for the binary inside it
        if os.path.isdir(path):
            for subpath in [
                os.path.join(path, "bin", tool_version),
                os.path.join(path, tool_version),
            ]:
                if os.path.isfile(subpath) and os.access(subpath, os.X_OK):
                    return os.path.abspath(subpath)
        elif os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)

    return None


def build_netmhc_executor(args, parser):
    """
    Resolve the native NetMHC binary and return a callable executor.

    Returns a callable executor and a description string.
    """
    native_netmhc = resolve_native_netmhc_path(
        getattr(args, "native_netmhc_path", None),
        getattr(args, "netmhc_tool", "netMHC")   # F5: netMHC-4.0 is the only supported tool
    )
    if not native_netmhc:
        parser.error(
            "Native NetMHC binary not found. Provide --native-netmhc-path, "
            "or set NETMHC_PATH / NETMHC_HOME environment variable."
        )

    verbose_flag = getattr(args, "verbose", True)
    if verbose_flag:
        print(f"Execution mode: native NetMHC ({native_netmhc})")

    def _runner(fasta_file, output_file, timeout, alleles):
        return _run_native_netmhc(fasta_file, output_file, timeout, native_netmhc, alleles)
    return _runner, f"native:{native_netmhc}"


def _run_native_netmhc(fasta_file, output_file, timeout, netmhc_path, alleles=None):
    """
    Execute NetMHC using native installation.

    Args:
        fasta_file: Input FASTA file path
        output_file: Output file path
        timeout: Command timeout in seconds
        netmhc_path: Path to NetMHC executable
        alleles: List of HLA alleles to predict

    Returns:
        tuple: (success, output_content, error_message)
    """
    if not alleles:
        alleles = ["HLA-A0201"]  # Default allele (netMHC-4.0 format)

    all_outputs = []

    try:
        # NetMHC processes one allele at a time, so run for each allele
        for allele in alleles:
            # Build native NetMHC command
            # Format: netMHC -a HLA-A*02:01 -f input.fasta
            cmd = [netmhc_path, "-a", allele, "-f", fasta_file]

            # netMHC binary expects $NETMHC pointing to the platform dir
            # (the parent of bin/), e.g. netMHC-4.0/Darwin_x86_64/
            env = os.environ.copy()
            netmhc_bin_dir = os.path.dirname(os.path.abspath(netmhc_path))
            env["NETMHC"] = os.path.dirname(netmhc_bin_dir) if os.path.basename(netmhc_bin_dir) == "bin" else netmhc_bin_dir

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                env=env,
            )

            # netMHC-4.0 returns exit code 1 even on success;
            # check stderr and stdout for actual error indicators
            if result.stderr and "error" in result.stderr.lower():
                return False, result.stdout, result.stderr
            if "cannot be found" in result.stdout or (not result.stdout.strip()):
                return False, result.stdout, result.stderr or result.stdout

            all_outputs.append(result.stdout)

        # Combine all outputs
        combined_output = "\n".join(all_outputs)

        # Write to output file
        with open(output_file, 'w') as f:
            f.write(combined_output)

        return True, combined_output, None

    except subprocess.TimeoutExpired:
        return False, "", f"NetMHC command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)


def parse_netmhc_output(output_file):
    """
    Parse NetMHC output file and extract binding predictions.

    NetMHC output format (space-separated):
    pos HLA peptide Core Offset I_pos I_len D_pos D_len iCore Identity 1-log50k(aff) Affinity(nM) %Rank BindLevel

    Returns:
        list: List of prediction dictionaries with keys:
              pos, mhc_allele, peptide, core, affinity, rank, bind_level, identity
    """
    predictions = []

    try:
        with open(output_file, 'r') as f:
            in_prediction_section = False
            current_identity = ""

            for line in f:
                original_line = line
                line = line.strip()

                # Skip empty lines and comment lines
                if not line or line.startswith('#'):
                    continue

                # Skip separator lines (all dashes)
                if line.startswith('---'):
                    continue

                # Detect prediction section header
                if 'pos' in line.lower() and 'peptide' in line.lower() and 'HLA' in line:
                    in_prediction_section = True
                    continue

                # Stop at summary line
                if line.startswith('Protein ') and 'Allele' in line:
                    in_prediction_section = False
                    continue

                if in_prediction_section:
                    # Parse prediction line
                    # Format: "  0  HLA-A0201  TMDKSELVQ  ...  28676.59  43.00"
                    # Or:     "219  HLA-A0201  QLLRDNLTL  ...   167.10   1.50 <= WB"

                    fields = line.split()

                    # Need at least 14 fields for valid prediction
                    if len(fields) < 14:
                        continue

                    try:
                        # Extract identity (sequence name) from field 10
                        identity = fields[10] if len(fields) > 10 else ""

                        # BindLevel is optional and appears as "<= WB" or "<= SB" (2 tokens)
                        bind_level = ""
                        if len(fields) >= 16 and fields[14] == "<=":
                            bind_level = fields[15]  # WB or SB
                        elif len(fields) == 15 and fields[14] not in ["<=", ""]:
                            # Sometimes just "WB" or "SB" without <=
                            bind_level = fields[14]

                        prediction = {
                            'pos': int(fields[0]),
                            'mhc_allele': fields[1],
                            'peptide': fields[2],
                            'core': fields[3],
                            'affinity': float(fields[12]),  # Affinity(nM)
                            'rank': float(fields[13]),      # %Rank
                            'bind_level': bind_level,
                            'identity': identity,
                        }
                        predictions.append(prediction)
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue

    except Exception as e:
        print(f"Error parsing NetMHC output {output_file}: {e}")
        return []

    return predictions


def infer_pos_base(predictions, sequence, sample=25):
    """Determine whether netMHC's `pos` column is 0- or 1-based, by checking it.

    Returns (base, qc_flag). The whole WT<->MUT register join is an arithmetic
    projection of this column onto the mutant protein, so the convention is
    verified against the submitted sequence -- a reported peptide must be exactly
    ``sequence[pos - base : pos - base + len(peptide)]`` -- instead of trusted.
    netMHC-4.0 emits 0-based positions and 0 is preferred when both conventions
    happen to reproduce every sampled peptide (a homopolymer stretch can).

    qc_flag is 'POS_CONVENTION_UNVERIFIED' when NEITHER base reproduces the
    peptides, which means the parser and the tool disagree and the register
    projection in compare_wt_mut_predictions cannot be trusted. It is reported
    on the row rather than silently defaulted away.
    """
    candidates = {0, 1}
    checked = 0
    for pred in predictions:
        peptide = pred.get('peptide') or ''
        pos = pred.get('pos')
        if not peptide or pos is None:
            continue
        for base in list(candidates):
            start = pos - base
            if start < 0 or sequence[start:start + len(peptide)] != peptide:
                candidates.discard(base)
        checked += 1
        if not candidates or checked >= sample:
            break
    if 0 in candidates:
        return 0, ''
    if 1 in candidates:
        return 1, ''
    return 0, 'POS_CONVENTION_UNVERIFIED'


def aa_allele_problem(variant):
    """Name why an aa-level allele cannot be spliced into a protein, or ''.

    parse_variant's amino-acid grammar accepts any run of letters as an allele,
    so tokens that are not residue pairs parse cleanly. format_aa_token renders a
    frameshift as '<X><pos>fs' and documents that form as deliberately NOT
    re-parseable, yet parse_variant(is_nt=False) reads 'L20fs' back as a 1->2
    insertion of the residues 'f' and 's'. Splicing that produces a protein
    carrying the literal text 'fs', which netMHC then scores and this pipeline
    reports as a successful measurement of a variant nobody declared.

    A '*' is allowed through here because it IS meaningful -- the protein ends
    there -- and the caller truncates at it, matching what the nucleotide path
    does in _non_snv_mutant_protein. Both alleles are checked over their WHOLE
    span, not just their first character.

    The frameshift rendering needs its own test rather than an alphabet test,
    because 'fs' upper-cases to Phe-Ser and so passes any alphabet check. Its
    shape is exactly what format_aa_token emits -- a single REF residue and the
    literal 'fs' -- and it is refused case-insensitively. That does refuse a
    genuine 1->2 delins whose ALT really is Phe-Ser, written 'L20FS'; the token
    alone cannot distinguish the two, and a refusal that names itself is
    recoverable where a fabricated mutant protein scored as real is not.
    """
    if len(variant.ref) == 1 and variant.alt.lower() == 'fs':
        return 'AA_FRAMESHIFT_TOKEN_NOT_COMPARABLE'
    for label, allele in (('ref', variant.ref), ('alt', variant.alt)):
        bad = sorted({c for c in allele.upper() if c not in AA_RESIDUES and c != '*'})
        if bad:
            return f"AA_TOKEN_NOT_RESIDUES:{label}_has_{''.join(bad)}"
    return ''


def workdir_stem(gene_name, mutation, limit=80):
    """Filesystem-safe, length-bounded stem for one mutant's temp files.

    The temp FASTA and output paths used to embed the raw token. That is fine
    for an SNV ('A118G') and fails outright for the variant classes this
    pipeline now accepts: a 250 nt insertion token is 256 characters, and the
    resulting name exceeds the 255-byte NAME_MAX of both APFS and ext4, so the
    mutant died with OSError before netMHC was ever invoked. The token is
    therefore bounded and made unique with a digest of the full token.

    A token that is already short and filesystem-safe -- every SNV -- is used
    unchanged, so existing temp paths are byte-identical. These files live in a
    per-gene temp directory deleted in main's `finally`; nothing reads the name.
    """
    safe = ''.join(c if (c.isalnum() or c in '.-_') else '_' for c in mutation)
    stem = f"{gene_name}_{safe}"
    if safe != mutation or len(stem) > limit:
        digest = hashlib.sha1(mutation.encode('utf-8')).hexdigest()[:12]
        stem = f"{stem[:limit]}_{digest}"
    return stem


def resolve_aa_alignment(mutation, wt_aa_seq, mut_aa_seq, wt_nt_seq=None, is_nt=True):
    """Project WT protein indices onto the mutant protein for one mutation.

    Returns (wt_to_mut, variant_class, aa_consequence, note).

    wt_to_mut is a list where entry i is the mutant index of WT residue i, or None
    where the edit removed it -- the shared `align_wt_to_mut` convention. It is
    None (meaning "identity, index for index") only when the two proteins are
    provably the same length and no frameshift is involved, which is exactly the
    SNV/MNV case; the SNV path therefore behaves exactly as it did before.

    The edit is recovered from the two PROTEIN sequences with `infer_edit_span`,
    because netMHC is a protein tool and this is the pair of sequences it actually
    scored. The one thing prefix/suffix trimming cannot recover is a frameshift --
    it would report the minimal replacement and then pair residues that are not
    counterparts -- so the frameshift flag comes from `protein_consequence` on the
    nucleotide token and is passed in explicitly.
    """
    variant = parse_variant(mutation, is_nt=is_nt)
    variant_class = variant.kind if variant is not None else 'unparsed'
    aa_consequence = ''
    notes = []
    frameshift = False

    # REF guard over the WHOLE span. The non-SNV path is already guarded inside
    # _non_snv_mutant_protein, but the SNV path reaches the mutant protein through
    # infer_aamutation_from_nt -> get_mutant_aa, which reads the base the sequence
    # actually carries and never compares it to the REF the token declares. A token
    # claiming A118G against an ORF holding T at 118 therefore yields a real mutant
    # built from T, i.e. a variant nobody asked for, with nothing to say so.
    #
    # The row is FLAGGED, not dropped: the predictions themselves are genuine
    # measurements of the protein that was submitted, and dropping the row would
    # trade one silent error for another. Comparing seq[pos] alone would not do --
    # a multi-base REF whose first base happens to match would pass.
    reference = wt_nt_seq if is_nt else wt_aa_seq
    if variant is not None and reference:
        end = variant.pos0 + len(variant.ref)
        if end > len(reference):
            notes.append(f'REF_SPANS_PAST_REFERENCE:{end}>{len(reference)}')
        else:
            observed = reference[variant.pos0:end].upper()
            if observed != variant.ref.upper():
                notes.append(f'REF_MISMATCH:reference_has_{observed}')

    if variant is None:
        notes.append('variant_token_unparsed')
    elif is_nt:
        if not wt_nt_seq:
            notes.append('aa_consequence_unavailable:no_orf')
        else:
            cons = protein_consequence(variant, wt_nt_seq)
            if cons is None:
                notes.append('aa_consequence_unavailable:variant_outside_orf')
            else:
                aa_consequence = cons['aa_consequence']
                frameshift = aa_consequence == 'frameshift'
    else:
        # An aa-level token carries its own consequence class: `kind` already
        # distinguishes snv/mnv/insertion/deletion/delins at the residue level.
        aa_consequence = variant.kind
        if '*' in variant.alt:
            # The mutant protein was truncated at this stop when it was built,
            # the same way _non_snv_mutant_protein truncates the nucleotide
            # path. Recorded because the length change on this row comes from
            # the truncation, not from the token's own length delta.
            notes.append('AA_ALT_TRUNCATED_AT_STOP')

    if not frameshift and len(wt_aa_seq) == len(mut_aa_seq):
        # align_wt_to_mut would return [0, 1, ... n-1] here; say so by returning
        # None rather than materialising an identity list per mutation.
        # aa_consequence stays EMPTY when it could not be derived (unparsed token,
        # no ORF, variant outside the ORF). Defaulting it to 'snv' would assert a
        # consequence class nothing here established; `notes` carries the reason.
        return None, variant_class, aa_consequence, ';'.join(notes)

    offset, ref_len, alt_len = infer_edit_span(wt_aa_seq, mut_aa_seq,
                                               frameshift=frameshift)
    wt_to_mut = align_wt_to_mut(len(wt_aa_seq), offset, ref_len, alt_len)
    notes.append(f"aa_edit:offset_{offset};ref_{ref_len}aa;alt_{alt_len}aa")
    return wt_to_mut, variant_class, aa_consequence, ';'.join(notes)


def _event_row(gene_name, mutation, variant_class, aa_consequence, allele,
               wt_pred, mut_pred, threshold, align_status, qc_flags):
    """Build one event row. Shared by the aligned, WT-only and MUT-only paths so
    the three cannot drift apart.

    A side with no counterpart leaves its columns EMPTY and names the reason in
    qc_flags. The former code substituted rank=100.0 / affinity=50000.0 for an
    absent side, which is not a measurement: it produced |delta_rank| ~ 99 that
    dominated max_abs_delta_rank and top_event selection, and it was
    indistinguishable from a genuine non-binder that netMHC actually scored.
    classification_code is likewise empty rather than 0 on an unpaired row -- 0
    means "stable", i.e. measured and unchanged.
    """
    wt_rank = wt_pred['rank'] if wt_pred else ''
    mut_rank = mut_pred['rank'] if mut_pred else ''
    wt_affinity = wt_pred['affinity'] if wt_pred else ''
    mut_affinity = mut_pred['affinity'] if mut_pred else ''

    if wt_pred and mut_pred:
        delta_rank = mut_pred['rank'] - wt_pred['rank']
        delta_affinity = mut_pred['affinity'] - wt_pred['affinity']

        # Classify by NetMHC binding band (per the tool's own definitions):
        #   2 = SB (strong binder, %Rank < threshold, default 0.5)
        #   1 = WB (weak binder,   threshold <= %Rank < 2.0)
        #   0 = NB (non-binder,    %Rank >= 2.0)
        # Strength order SB(2) > WB(1) > NB(0). Classification is a band transition;
        # delta_rank/delta_affinity are retained as reported magnitudes only.
        # (The former delta_rank <>+/-5 test was dead: inside the both-bind branch
        # delta_rank is bounded to [-2, 2] and could never reach +/-5.)
        wt_band = 2 if wt_rank < threshold else (1 if wt_rank < 2.0 else 0)
        mut_band = 2 if mut_rank < threshold else (1 if mut_rank < 2.0 else 0)

        if wt_band == 0 and mut_band > 0:
            classification, classification_code = "gained", 2
        elif wt_band > 0 and mut_band == 0:
            classification, classification_code = "lost", -2
        elif mut_band > wt_band:
            classification, classification_code = "strengthened", 1
        elif mut_band < wt_band:
            classification, classification_code = "weakened", -1
        else:
            classification, classification_code = "stable", 0
    else:
        # No pair, no band transition. Naming one here would fabricate a
        # comparison against a sentinel.
        delta_rank = delta_affinity = ''
        classification = classification_code = ''

    return {
        'Gene': gene_name,
        'pkey': mint_pkey(gene_name, mutation),
        'mutation': mutation,
        'wt_peptide': wt_pred['peptide'] if wt_pred else '',
        'mut_peptide': mut_pred['peptide'] if mut_pred else '',
        'pos': wt_pred['pos'] if wt_pred else '',
        'mhc_allele': allele,
        'wt_rank': wt_rank,
        'mut_rank': mut_rank,
        'delta_rank': delta_rank,
        'wt_affinity': wt_affinity,
        'mut_affinity': mut_affinity,
        'delta_affinity': delta_affinity,
        'bind_level_wt': wt_pred['bind_level'] if wt_pred else '',
        'bind_level_mut': mut_pred['bind_level'] if mut_pred else '',
        'classification': classification,
        'classification_code': classification_code,
        'variant_class': variant_class,
        'aa_consequence': aa_consequence,
        'mut_pos': mut_pred['pos'] if mut_pred else '',
        'align_status': align_status,
        'qc_flags': qc_flags,
    }


def compare_wt_mut_predictions(gene_name, mutation, wt_preds, mut_preds, threshold=0.5,
                               wt_to_mut=None, variant_class='snv', aa_consequence='snv',
                               wt_pos_base=0, mut_pos_base=0):
    """
    Compare WT and MUT predictions to classify epitope changes.

    Classification is a NetMHC binding-band transition (see _event_row):
    gained / lost / strengthened / weakened / stable.

    Args:
        gene_name: Gene name
        mutation: Mutation identifier
        wt_preds: List of WT predictions
        mut_preds: List of MUT predictions
        threshold: Rank threshold for strong binder (default 0.5)
        wt_to_mut: WT->MUT residue projection from resolve_aa_alignment, or None
            for the identity case (equal-length proteins: SNV, MNV, synonymous).
        variant_class / aa_consequence: carried onto every row.
        wt_pos_base / mut_pos_base: the 0/1 base of netMHC's `pos` column for each
            sequence, from infer_pos_base.

    Returns:
        List of event dictionaries, one per register in the UNION of the two
        alleles' registers: first every WT register sorted by (allele, wt_pos),
        then the mutant-only registers sorted by (allele, mut_pos). Both blocks
        are ordered by a property of the record, so the output is deterministic.
        (The former set-of-keys iteration left row order dependent on
        per-process string hashing.)
    """
    # Lookup maps keyed on the REGISTER (allele, pos) -- NOT the peptide. A window
    # covering the mutated residue has different WT vs MUT peptide strings by
    # construction; keying on the sequence would give them different keys so they
    # would never be paired at all.
    #
    # The register is NOT the same integer on both sides once the edit changes the
    # protein's length: WT residue i and MUT residue i are then different residues.
    # `wt_to_mut` is the projection that says which mutant register a WT register
    # became, and registers with no counterpart are reported as such instead of
    # being paired with a sentinel.
    wt_map = {(p['mhc_allele'], p['pos']): p for p in wt_preds}
    mut_map = {(p['mhc_allele'], p['pos']): p for p in mut_preds}

    def _project(pos):
        """WT register position -> (MUT register position, reason it has none).

        The reason is returned rather than inferred by the caller: a projection
        can come back empty for two different causes, and naming both of them
        'register_deleted' would report a deletion that did not happen.
        """
        if wt_to_mut is None:
            # Identity in residue space, but still rebased: the two sequences are
            # only guaranteed to share netMHC's pos convention, not to have been
            # measured with the same one. Identical when the bases agree.
            return pos - wt_pos_base + mut_pos_base, ''
        idx = pos - wt_pos_base
        if not (0 <= idx < len(wt_to_mut)):
            # netMHC reported a window start that is not a residue of the WT
            # protein the projection was built from, i.e. the tool and this
            # pipeline disagree about the submitted sequence. Nothing was
            # deleted here.
            return None, f'no_mut_counterpart:wt_pos_{pos}_outside_alignment'
        mapped = wt_to_mut[idx]
        if mapped is None:
            return None, 'no_mut_counterpart:register_deleted'
        return mapped + mut_pos_base, ''

    events = []
    claimed_mut_keys = set()

    for (allele, wt_pos), wt_pred in sorted(wt_map.items()):
        mut_pos, gap_reason = _project(wt_pos)
        mut_pred = mut_map.get((allele, mut_pos)) if mut_pos is not None else None
        if mut_pred is not None:
            claimed_mut_keys.add((allele, mut_pos))
            events.append(_event_row(gene_name, mutation, variant_class, aa_consequence,
                                     allele, wt_pred, mut_pred, threshold, 'aligned', ''))
        elif mut_pos is None:
            events.append(_event_row(gene_name, mutation, variant_class, aa_consequence,
                                     allele, wt_pred, None, threshold, 'wt_only',
                                     gap_reason))
        else:
            events.append(_event_row(gene_name, mutation, variant_class, aa_consequence,
                                     allele, wt_pred, None, threshold, 'wt_only',
                                     f'no_mut_counterpart:no_prediction_at_mut_pos_{mut_pos}'))

    # Registers the mutant protein gained have no WT coordinate and so appear in
    # none of the rows above. Emit them explicitly; without them an insertion or a
    # frameshift is invisible here and the run still reports full alignment,
    # because every WT register does keep an entry.
    for (allele, mut_pos), mut_pred in sorted(mut_map.items()):
        if (allele, mut_pos) in claimed_mut_keys:
            continue
        events.append(_event_row(gene_name, mutation, variant_class, aa_consequence,
                                 allele, None, mut_pred, threshold, 'mut_only',
                                 'no_wt_counterpart:register_inserted'))

    return events


def blank_summary_row(gene_name, mutation, qc_flags, variant_class='', aa_consequence=''):
    """A named summary row for a mutation that produced no epitope comparison.

    Every metric column is EMPTY and qc_flags carries the reason. A rejected or
    failed token gets one of these instead of vanishing: a missing row is
    indistinguishable from "this gene never had that mutation", and a zero-filled
    row is worse still because it reads as a real measurement of no change.
    """
    row = {name: '' for name in SUMMARY_FIELDNAMES}
    row.update({
        'pkey': mint_pkey(gene_name, mutation),
        'Gene': gene_name,
        'mutation': mutation,
        'variant_class': variant_class,
        'aa_consequence': aa_consequence,
        'qc_flags': qc_flags,
    })
    return row


def summarize_epitope_changes(gene_name, mutation, events, threshold=0.5,
                              variant_class='snv', aa_consequence='snv',
                              align_note='', wt_aa_len=None, mut_aa_len=None):
    """
    Generate summary statistics for a mutation's epitope changes.

    Args:
        gene_name: Gene name
        mutation: Mutation identifier
        events: List of epitope event dictionaries
        threshold: Strong-binder %Rank cutoff. Was hardcoded to 0.5 here while
            --threshold reached only the band classifier, so the two disagreed
            whenever the user changed it. n_epitopes_* now means exactly "band 2
            (SB)" -- the same `< threshold` test the classifier uses -- rather than
            the former `<= 0.5`.
        variant_class / aa_consequence / align_note: carried through from
            resolve_aa_alignment.
        wt_aa_len / mut_aa_len: protein lengths, used only to report the length
            change in qc_flags.

    Returns:
        Dictionary with summary statistics
    """
    # Deltas and band transitions exist only where a register has BOTH sides.
    aligned = [e for e in events if e['align_status'] == 'aligned']
    wt_only = [e for e in events if e['align_status'] == 'wt_only']
    mut_only = [e for e in events if e['align_status'] == 'mut_only']

    count_gained = sum(1 for e in aligned if e['classification'] == 'gained')
    count_lost = sum(1 for e in aligned if e['classification'] == 'lost')
    count_strengthened = sum(1 for e in aligned if e['classification'] == 'strengthened')
    count_weakened = sum(1 for e in aligned if e['classification'] == 'weakened')
    count_stable = sum(1 for e in aligned if e['classification'] == 'stable')

    # An epitope is an epitope whether or not it has a counterpart in the other
    # allele, so these count over ALL events with a rank on the relevant side --
    # a WT epitope the edit deleted still existed in the WT.
    n_epitopes_wt = sum(1 for e in events if e['wt_rank'] != '' and e['wt_rank'] < threshold)
    n_epitopes_mut = sum(1 for e in events if e['mut_rank'] != '' and e['mut_rank'] < threshold)

    abs_deltas = [abs(e['delta_rank']) for e in aligned]
    # Ties on |delta_rank| are common -- every register of a synonymous change ties
    # at 0.0 -- and a bare max() then returns whichever tied row happened to come
    # first, which made this column depend on the order the events were built in.
    # (allele, pos) breaks the tie by a property of the record instead.
    top_event = (max(aligned, key=lambda e: (abs(e['delta_rank']), e['mhc_allele'], e['pos']))
                 if aligned else None)

    qc_flags = []
    if align_note:
        qc_flags.append(align_note)
    if wt_aa_len is not None and mut_aa_len is not None and wt_aa_len != mut_aa_len:
        # Count over the UNION of both alleles' registers, not just WT ones.
        # "aligned N/N" over WT registers alone reports an insertion as fully
        # aligned however large it is, because every WT register does keep an
        # entry -- the registers with no counterpart are all on the mutant side.
        qc_flags.append(
            f"length_changed:{mut_aa_len - wt_aa_len:+d}aa;"
            f"aligned_{len(aligned)}/{len(events)};"
            f"wt_only_{len(wt_only)};mut_only_{len(mut_only)}"
        )
    if n_epitopes_wt == 0:
        qc_flags.append("no_wt_epitopes")
    if n_epitopes_mut == 0:
        qc_flags.append("no_mut_epitopes")
    if not events:
        qc_flags.append("no_predictions")
    elif not aligned:
        qc_flags.append("no_aligned_registers")

    return {
        'pkey': mint_pkey(gene_name, mutation),
        'Gene': gene_name,
        'mutation': mutation,
        'n_epitopes_wt': n_epitopes_wt,
        'n_epitopes_mut': n_epitopes_mut,
        'count_gained': count_gained,
        'count_lost': count_lost,
        'count_strengthened': count_strengthened,
        'count_weakened': count_weakened,
        'count_stable': count_stable,
        # Empty, not 0.0, when nothing aligned: 0.0 reads as "measured, no change".
        'max_abs_delta_rank': max(abs_deltas) if abs_deltas else '',
        # fsum, not sum: plain float addition is not associative, so the total
        # depended on the order the events happened to be built in and the column
        # was not reproducible across runs. fsum is the correctly-rounded total.
        'sum_abs_delta_rank': math.fsum(abs_deltas) if abs_deltas else '',
        'top_event_type': top_event['classification'] if top_event else '',
        'top_event_allele': top_event['mhc_allele'] if top_event else '',
        'top_event_wt_peptide': top_event['wt_peptide'] if top_event else '',
        'top_event_mut_peptide': top_event['mut_peptide'] if top_event else '',
        'top_event_delta_rank': top_event['delta_rank'] if top_event else '',
        'qc_flags': ';'.join(qc_flags) if qc_flags else '',
        'variant_class': variant_class,
        'aa_consequence': aa_consequence,
        'count_wt_only': len(wt_only),
        'count_mut_only': len(mut_only),
        'n_registers_aligned': len(aligned),
        'n_registers_total': len(events),
    }


# The three writers used to `return` before calling write_tsv when their row list
# was empty, so a gene with no rows got NO FILE -- indistinguishable downstream from
# a gene the pipeline never reached, and a stale file from a previous run survived
# in its place. write_tsv already emits a header-only file when fieldnames are
# supplied (utility.py), which is what every other BFF pipeline produces.
def write_summary_tsv(summary_rows, output_file):
    """Write summary TSV file (header-only when there are no rows)."""
    write_tsv(summary_rows, output_file, SUMMARY_FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(summary_rows)} summary rows to {output_file}")


def write_events_tsv(events, output_file):
    """Write events TSV file (header-only when there are no rows)."""
    write_tsv(events, output_file, EVENT_FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(events)} events to {output_file}")


def write_sites_tsv(sites, output_file):
    """Write sites TSV file (header-only when there are no rows)."""
    write_tsv(sites, output_file, SITE_FIELDNAMES, extrasaction='ignore')
    print(f"Wrote {len(sites)} site predictions to {output_file}")


def first_line(text, default='unknown'):
    """First non-blank line of a tool's error text, for a one-field reason code.

    `.strip().splitlines()[0]` raises IndexError on whitespace-only text, which
    would turn a netMHC failure into a crash while building the row that reports
    that failure.
    """
    for line in str(text or '').splitlines():
        line = line.strip()
        if line:
            return line
    return default


# netMHC-4.0's default peptide length. Not passed on the command line -- the
# invocation deliberately stays as it was -- so this is only used to state the
# expected window count in a diagnostic message.
NETMHC_PEPTIDE_LENGTH = 9

# How much of the tool's own output to echo when it produced nothing parseable.
NETMHC_OUTPUT_ECHO_LINES = 40


def netmhc_output_problem(predictions, sequence, peptide_length=NETMHC_PEPTIDE_LENGTH):
    """Name why a netMHC run yielded no rows, or '' when it yielded some.

    netMHC scores EVERY window of the submitted sequence and prints a row for
    each one: a non-binder is a row with a blank BindLevel, not an absent row.
    A protein of length L therefore has exactly L - l + 1 rows per allele, and
    ZERO rows is never a measurement -- it means the tool refused the input, or
    the tool's output layout and parse_netmhc_output disagree.

    Both of those were previously reported as "0 predictions" on stdout and then
    counted as a successful mutation, so a run in which netMHC produced nothing
    at all ended with "mutations 8/8 ok". A run that measured nothing is named
    as a failure here instead.
    """
    if predictions:
        return ''
    if len(sequence) < peptide_length:
        return f'SEQUENCE_SHORTER_THAN_PEPTIDE:{len(sequence)}<{peptide_length}'
    expected = len(sequence) - peptide_length + 1
    return f'NETMHC_NO_PARSEABLE_ROWS:expected_{expected}_windows_per_allele'


def echo_netmhc_output(label, raw_output, limit=NETMHC_OUTPUT_ECHO_LINES):
    """Print the tool's own output when nothing could be parsed from it.

    The per-gene workdir is deleted in main's `finally`, so the .out file is gone
    by the time the user sees "0 predictions" -- there was no way to tell a tool
    refusal from a parser mismatch without re-running netMHC by hand. The text is
    already in memory: the executor returns it and both call sites discarded it.
    """
    lines = str(raw_output or '').splitlines()
    print(f"  [netmhc] {label}: no parseable rows; "
          f"tool wrote {len(lines)} lines:", file=sys.stderr)
    for line in lines[:limit]:
        print(f"  | {line}", file=sys.stderr)
    if len(lines) > limit:
        print(f"  | ... {len(lines) - limit} more lines", file=sys.stderr)


def iter_mutation_tokens(mapping_file):
    """Yield every mutation token in a mapping file, in order, deduplicated.

    ACCOUNTING ONLY. This mirrors the single-column / CSV detection that
    build_mutant_sequences_for_gene performs internally so that the pipeline can
    name the tokens that function chose not to emit. It decides nothing: the
    sequences themselves always come from build_mutant_sequences_for_gene.

    The duplication exists because build_mutant_sequences_for_gene returns only
    {header: sequence} and prints its skipped-token list instead of returning it,
    so a caller cannot otherwise distinguish "this token was rejected" from "this
    token was never in the file" -- which contract D requires. utility.py is not
    modified here; the shared fix would be for it to return (sequences, skipped).
    """
    if not mapping_file or not os.path.exists(mapping_file):
        return []

    mutant_keys = ['mutant', 'mutation', 'nt_mutation', 'ntmutant']
    try:
        with open(mapping_file, 'r') as handle:
            lines = handle.readlines()
    except OSError:
        return []

    is_single_column = True
    if lines and ',' in lines[0]:
        first_line_lower = lines[0].lower()
        if any(k in first_line_lower for k in ['mutant', 'mutation', 'aamutant']):
            is_single_column = False

    tokens = []
    if is_single_column:
        for line in lines:
            token = line.strip()
            if not token or token.lower() == 'mutant':
                continue
            tokens.append(token.replace(" ", ""))
    else:
        with open(mapping_file, 'r') as handle:
            for row in csv.DictReader(handle):
                for key in mutant_keys:
                    if key in row and row[key]:
                        tokens.append(row[key].strip().replace(" ", ""))
                        break

    seen = set()
    ordered = []
    for token in tokens:
        if token not in seen:
            seen.add(token)
            ordered.append(token)
    return ordered


def pending_mutation_tokens(mapping_file, gene_name, log_path, failure_map):
    """Tokens that SHOULD have produced a mutant, after the deliberate filters.

    Applies the same two filters build_mutant_sequences_for_gene applies -- the
    validation-log allow-list and should_skip_mutation -- so a token that was
    deliberately excluded is not then reported as a failure.
    """
    tokens = iter_mutation_tokens(mapping_file)
    if not tokens:
        return []

    allowed = None
    if log_path:
        try:
            allowed = {
                entry.split(',')[0].strip().upper()
                for entry in trim_muts(mapping_file, log=log_path, gene_name=gene_name)
                if entry
            }
        except Exception:
            allowed = None

    pending = []
    for token in tokens:
        if allowed and token.upper() not in allowed:
            continue
        if should_skip_mutation(gene_name, token, failure_map):
            continue
        pending.append(token)
    return pending


def explain_missing_token(token, wt_nt_seq, wt_aa_seq, is_nt):
    """Name why a token produced no mutant protein: (reason, variant_class, aa_consequence).

    The checks are the same ones build_mutant_sequences_for_gene /
    _non_snv_mutant_protein apply, re-run here purely to attach a reason code to
    the row. Bounds are tested on the END of the REF span, and the REF guard
    compares the WHOLE span -- a single-base check passes on any multi-base REF
    whose first base happens to match.
    """
    # A gd./ch. token is WELL FORMED and out of SCOPE for a protein tool, not
    # malformed. UNPARSEABLE_TOKEN gives it the same reason as actual garbage;
    # build_mutant_sequences_for_gene already names it correctly on stderr, and
    # this triple is the only record that reaches the TSV.
    if is_intronic_token(token):
        return 'NON_ORF_TOKEN:no_reading_frame_at_protein_level', 'non_orf', ''
    variant = parse_variant(token, is_nt=is_nt)
    if variant is None:
        return 'UNPARSEABLE_TOKEN', 'unparsed', ''

    variant_class = variant.kind
    ref_seq = wt_nt_seq if is_nt else wt_aa_seq
    if not ref_seq:
        return 'NO_REFERENCE_SEQUENCE', variant_class, ''
    if variant.pos0 + len(variant.ref) > len(ref_seq):
        return 'REF_SPANS_PAST_ORF', variant_class, ''
    observed = ref_seq[variant.pos0:variant.pos0 + len(variant.ref)].upper()
    if observed != variant.ref.upper():
        return f'REF_MISMATCH:orf_has_{observed}', variant_class, ''

    aa_consequence = ''
    if is_nt:
        cons = protein_consequence(variant, wt_nt_seq)
        if cons is not None:
            aa_consequence = cons['aa_consequence']
        # infer_aamutation_from_nt -- the SNV path inside
        # build_mutant_sequences_for_gene -- returns None for a stop-gain, because
        # codon_to_aa renders a stop as the 4-character string 'Stop' which must
        # not be spliced in as a residue. Real, pre-existing, and named here
        # rather than dropped.
        if variant.is_snv and aa_consequence == 'stop_gained':
            return 'STOP_GAINED:snv_path_drops_stop_tokens', variant_class, aa_consequence
    elif not variant.is_snv:
        # Same two checks the aa backfill in main() applies, in the same order,
        # so the reason on the row is the reason the mutant was not built.
        problem = aa_allele_problem(variant)
        if problem:
            return problem, variant_class, variant.kind
        # validate=False is safe: the REF was compared over its whole span above.
        spliced = splice_seq(wt_aa_seq, variant.pos0, variant.ref, variant.alt,
                             validate=False).split('*')[0]
        if not spliced:
            return 'AA_MUTANT_EMPTY_BEFORE_FIRST_STOP', variant_class, variant.kind
        return 'AA_NON_SNV_SPLICE_FAILED', variant_class, variant.kind

    return 'MUTANT_PROTEIN_NOT_BUILT', variant_class, aa_consequence


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="NetMHC pipeline for MHC binding prediction with WT/mutant comparison",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # -i/-o, matching every other non-Nextflow pipeline. The positional forms are
    # kept as hidden optional trailing args so existing invocations still parse;
    # the flag wins when both are given.
    parser.add_argument('-i', '--input', dest='input_flag', metavar='INPUT',
                        help='DIRECTORY MODE: variant_mapping output root '
                             '(<root>/<GENE>/fastas/ + <root>/<GENE>/mappings/). '
                             'Also accepts a single WT FASTA or a flat directory of them.')
    parser.add_argument('-o', '--output', dest='output_flag', metavar='OUTPUT',
                        help='Output base directory; writes <output>/<GENE>/NetMHC/')
    parser.add_argument('input', nargs='?', help=argparse.SUPPRESS)
    parser.add_argument('output', nargs='?', help=argparse.SUPPRESS)

    # MHC-specific options
    parser.add_argument('-a', '--alleles', nargs='+',
                       help='HLA alleles to predict (e.g., HLA-A*02:01 HLA-B*07:02). If not specified, uses default set.')
    # Only netMHC-4.0 is supported: the invocation (allele format, `-a`/`-f`) and
    # the parser (15-column layout, 'HLA' header) are hardcoded to it. netMHCpan /
    # netMHCII use different output layouts the parser cannot read and would
    # silently yield zero predictions, so they are not offered.
    parser.add_argument('-nt', '--netmhc-tool', choices=['netMHC'],
                       default='netMHC',
                       help='NetMHC tool to use (only netMHC-4.0 is supported)')
    # Execution backend
    parser.add_argument('-nnp', '--native-netmhc-path',
                       help='Path to native NetMHC executable')

    # Processing options
    parser.add_argument('--mutations', '-m',
                       help='Mutation file or directory of mutation CSVs (single-column NT mutations)')
    parser.add_argument('-l', '--log',
                       help='Validation log file to skip failed mutations')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                       help='Binding rank threshold for strong binders (default: 0.5)')
    parser.add_argument('-bs', '--batch-size', type=int, default=100,
                       help='(deprecated; unused -- netMHC input is no longer batched)')
    parser.add_argument('-ti', '--timeout', type=int, default=600,
                       help='Command timeout in seconds (default: 600)')
    parser.add_argument('--max-workers', type=int, default=4,
                       help='Concurrent netMHC runs per gene (default: 4; set 1 for serial)')

    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Validate arguments
    # Flag wins over the positional; both fold into args.input/args.output so
    # nothing downstream has to know which form the caller used.
    args.input = args.input_flag or args.input
    args.output = args.output_flag or args.output

    if not args.input or not args.output:
        parser.error("input and output arguments are required")

    # One root supplies both: <root>/<GENE>/mappings/ sits beside
    # <root>/<GENE>/fastas/, so the input IS the mapping location in directory
    # mode. Explicit values always win; FILE MODE is unaffected.
    args.mutations = derive_mutations_root(args.mutations, args.input, label="netmhc")
    if not args.mutations:
        parser.error("--mutations is REQUIRED for full-pipeline mode")

    # Build NetMHC executor
    executor, exec_desc = build_netmhc_executor(args, parser)

    # Load validation failures if provided
    failure_map = None
    if args.log:
        failure_map = load_validation_failures(args.log)
        if args.verbose:
            print(f"Loaded validation failures from {args.log}")

    # Load WT sequences - handles both files and directories
    wt_sequences, temp_holder = load_wt_sequence_map(args.input, wt_header='ORF')
    if not wt_sequences:
        print(f"Error: No FASTA sequences found in {args.input}")
        return 1

    if args.verbose:
        print(f"Loaded {len(wt_sequences)} WT sequences")

    # Discover mutation files.
    #
    # Per-gene layout first. The glob below is non-recursive, so at a
    # variant_mapping root it saw only <GENE>/ directories and left
    # mutation_files empty -- every gene then had no tokens to request while
    # <GENE>/mappings/mutations/<GENE>_mutations.csv sat three levels down.
    # discover_mutation_files selects by directory rather than by filename, which
    # is what keeps the six other CSV types in that tree out.
    mutation_files = {}
    if args.mutations:
        mut_path = Path(args.mutations)
        if mut_path.is_file():
            mutation_files[extract_gene_from_filename(mut_path.name)] = str(mut_path)
        elif mut_path.is_dir():
            mutation_files = discover_mutation_files(str(mut_path))
            if not mutation_files:
                for csv_file in mut_path.glob("*.csv"):
                    mutation_files[extract_gene_from_filename(csv_file.name)] = str(csv_file)

    # Run-level accounting. Every gene and every token lands in exactly one bucket,
    # so "processed + skipped" reconciles against the input instead of the run
    # ending with an unexplained shortfall.
    run_summary = {
        "genes_total": len(wt_sequences), "genes_processed": 0, "genes_skipped": 0,
        "mutations_total": 0, "mutations_successful": 0, "mutations_unsuccessful": 0,
        "skipped_genes": [], "genes_without_mutation_file": [], "unsuccessful": [],
    }

    # Process each gene
    for gene_name, wt_seq in wt_sequences.items():
        gene_summary_rows = []
        gene_events = []
        gene_sites = []
        if args.verbose:
            print(f"\nProcessing gene: {gene_name}")

        # Auto-detect the WT alphabet: nucleotide -> translate; protein -> use
        # directly; codon-encoded -> skip. Prevents the unconditional-translate
        # crash (TranslationError) that amino-acid input previously caused.
        try:
            detected = detect_alphabet(wt_seq)
        except ValueError:
            print(f"Warning: {gene_name} WT sequence is empty, skipping")
            run_summary["genes_skipped"] += 1
            run_summary["skipped_genes"].append({"gene": gene_name, "reason": "empty_wt_sequence"})
            continue
        if detected == 'nucleotide':
            wt_nt_seq = wt_seq
            wt_aa_seq = translate_orf_sequence(wt_nt_seq)
            build_input_type = 'nt'
            if not wt_aa_seq:
                print(f"Warning: Could not translate {gene_name}, skipping")
                run_summary["genes_skipped"] += 1
                run_summary["skipped_genes"].append({"gene": gene_name, "reason": "untranslatable_orf"})
                continue
        elif detected == 'protein':
            wt_nt_seq = None
            wt_aa_seq = wt_seq
            build_input_type = 'aa'
        else:  # codon-encoded
            print(f"Warning: {gene_name} WT looks codon-encoded; netMHC needs nt or aa, skipping")
            run_summary["genes_skipped"] += 1
            run_summary["skipped_genes"].append({"gene": gene_name, "reason": "codon_encoded_wt"})
            continue

        # Build mutant sequences. non_snp=True unconditionally: whether a token is
        # a substitution or an indel is a fact of the record, not a user
        # preference, and the token grammar is uniquely decodable so no parser mode
        # has to be selected. A per-run boolean could only permit everything or
        # refuse everything, neither of which is what the data needs.
        mapping_file = mutation_files.get(gene_name)
        if not mapping_file:
            run_summary["genes_without_mutation_file"].append(gene_name)
        # Sequence names are {GENE}-{sha}; pkey_map carries name -> token so the
        # token spelling is recovered from the map, not by parsing the name.
        pkey_map = {}
        mutant_seqs = build_mutant_sequences_for_gene(
            gene_name, wt_nt_seq, wt_aa_seq, mapping_file, args.log, failure_map,
            input_type=build_input_type, non_snp=True, pkey_map=pkey_map
        )

        pending_tokens = pending_mutation_tokens(mapping_file, gene_name, args.log, failure_map)

        # An aa-level non-SNV token cannot be built by build_mutant_sequences_for_gene:
        # its non-SNV path splices at the NUCLEOTIDE level and returns None without an
        # ORF, and its aa fallback goes through get_mutation_data_bioAccurate, which is
        # single-residue only. With a protein FASTA there is no ORF, but the edit is
        # fully expressible on the protein itself, so splice it here rather than
        # rejecting a token that is perfectly well defined.
        if build_input_type == 'aa':
            for token in pending_tokens:
                header = mint_pkey(gene_name, token)
                if header in mutant_seqs:
                    continue
                variant = parse_variant(token, is_nt=False)
                if variant is None or variant.is_snv:
                    continue
                if aa_allele_problem(variant):
                    # An allele that is not a residue string. Splicing it would
                    # write literal non-residue text into the protein handed to
                    # netMHC; explain_missing_token names it on the row instead.
                    continue
                try:
                    spliced = splice_seq(wt_aa_seq, variant.pos0,
                                         variant.ref, variant.alt, validate=True)
                except ValueError:
                    # Out of range or REF disagreement. Left unbuilt on purpose:
                    # the missing-token accounting below names the exact reason.
                    continue
                # Truncate at the first stop, exactly as the nucleotide path does
                # in _non_snv_mutant_protein: translation ends there, and a '*'
                # left embedded in the sequence makes netMHC score residues the
                # protein does not have.
                spliced = spliced.split('*')[0]
                if spliced:
                    mutant_seqs[header] = spliced
                    pkey_map[header] = token

        if args.verbose:
            print(f"  Generated {len(mutant_seqs)} mutant sequences")

        # Create temp directory for this gene
        gene_workdir = tempfile.mkdtemp(prefix=f"netmhc_{gene_name}_")
        gene_failure = None
        # Reason per mutation, filled by the workers. A mutation that produced no
        # comparison gets a NAMED summary row below rather than vanishing.
        mutation_failures = {}

        try:
            # Write WT FASTA
            wt_fasta = os.path.join(gene_workdir, f"{gene_name}_wt.fasta")
            write_fasta(Path(wt_fasta), {f"{gene_name}_WT": wt_aa_seq})

            # Run netMHC on the single-record WT FASTA. netMHC digests the
            # sequence into all k-mer windows itself, so no input batching is
            # needed (audit F8 -- the former split-by-record batching was a no-op).
            wt_output = os.path.join(gene_workdir, f"{gene_name}_wt.out")
            success, wt_raw, error = executor(wt_fasta, wt_output, args.timeout, args.alleles)
            if not success:
                # Was a bare `continue`, which skipped the output writes entirely
                # and left the gene with no files at all.
                print(f"Warning: NetMHC failed for {gene_name} WT: {error}")
                gene_failure = f"WT_NETMHC_FAILED:{first_line(error)}"
                wt_predictions = []
            else:
                wt_predictions = parse_netmhc_output(wt_output)
                if args.verbose:
                    print(f"  WT: {len(wt_predictions)} predictions")
                # netMHC reported an exit the executor accepted, but nothing came
                # back through the parser. That is a failure of the gene, not a WT
                # with no epitopes -- see netmhc_output_problem. Without this the
                # mutants were all compared against an EMPTY WT register map and
                # every one of them returned a full set of rows saying nothing
                # changed.
                wt_problem = netmhc_output_problem(wt_predictions, wt_aa_seq)
                if wt_problem:
                    echo_netmhc_output(f"{gene_name} WT", wt_raw)
                    gene_failure = f"WT_{wt_problem}"

            # The register join is an arithmetic projection of netMHC's `pos`
            # column, so its 0/1 base is verified against the WT protein once per
            # gene rather than assumed.
            wt_pos_base, wt_pos_note = infer_pos_base(wt_predictions, wt_aa_seq)

            # Add to sites output (all WT predictions)
            for pred in wt_predictions:
                gene_sites.append({
                    'Gene': gene_name,
                    'pkey': f"{gene_name}-WT",
                    'sequence_type': 'wt',
                    **pred
                })

            # Process each mutant in parallel. Every netMHC run is a separate OS
            # process (subprocess), so threads give real concurrency (the GIL is
            # released during subprocess.run) while keeping ONE netMHC invocation
            # per mutant -- unambiguous output attribution and per-mutant failure
            # isolation. Genes with thousands of mutations no longer pay a fully
            # serial process-spawn cost.
            def _process_mutant(item):
                mut_header, mut_aa_seq = item
                # Key is the {GENE}-{sha} sequence name (utility.py mint_pkey);
                # gene-name prefix so hyphenated genes (NKX2-1, HLA-A, MT-CO1)
                # parse correctly instead of split('-',1) corrupting the token.
                mutation = token_from_name(mut_header, gene_name, pkey_map)
                try:
                    # Bounded stem, not the raw token: an indel token can be
                    # hundreds of characters and blew past NAME_MAX (see
                    # workdir_stem). Unchanged for every SNV.
                    stem = workdir_stem(gene_name, mutation)
                    mut_fasta = os.path.join(gene_workdir, f"{stem}.fasta")
                    write_fasta(Path(mut_fasta), {mut_header: mut_aa_seq})
                    mut_output = os.path.join(gene_workdir, f"{stem}.out")
                    success, mut_raw, error = executor(mut_fasta, mut_output, args.timeout, args.alleles)
                    if not success:
                        print(f"Warning: NetMHC failed for {gene_name} {mutation}: {error}")
                        return None, mutation, f"MUT_NETMHC_FAILED:{first_line(error)}"
                    mut_predictions = parse_netmhc_output(mut_output)
                    if args.verbose:
                        print(f"  {mutation}: {len(mut_predictions)} predictions")
                    mut_problem = netmhc_output_problem(mut_predictions, mut_aa_seq)
                    if mut_problem:
                        # Same reasoning as the WT branch: zero rows out of a run
                        # netMHC did not report as failed is a defect, and the row
                        # it used to produce was an all-empty comparison counted
                        # as a successful mutation.
                        echo_netmhc_output(f"{gene_name} {mutation}", mut_raw)
                        return None, mutation, f"MUT_{mut_problem}"
                    site_rows = [
                        {'Gene': gene_name, 'pkey': mint_pkey(gene_name, mutation),
                         'sequence_type': 'mut', **pred}
                        for pred in mut_predictions
                    ]
                    # WT residue i and MUT residue i are different residues once the
                    # edit changes the protein's length, so the WT->MUT register
                    # join goes through the projection, not integer equality.
                    wt_to_mut, variant_class, aa_consequence, align_note = resolve_aa_alignment(
                        mutation, wt_aa_seq, mut_aa_seq, wt_nt_seq,
                        is_nt=(build_input_type == 'nt')
                    )
                    mut_pos_base, mut_pos_note = infer_pos_base(mut_predictions, mut_aa_seq)
                    notes = [n for n in (align_note, wt_pos_note, mut_pos_note) if n]
                    events = compare_wt_mut_predictions(
                        gene_name, mutation, wt_predictions, mut_predictions, args.threshold,
                        wt_to_mut=wt_to_mut, variant_class=variant_class,
                        aa_consequence=aa_consequence,
                        wt_pos_base=wt_pos_base, mut_pos_base=mut_pos_base
                    )
                    summary = summarize_epitope_changes(
                        gene_name, mutation, events, args.threshold,
                        variant_class=variant_class, aa_consequence=aa_consequence,
                        align_note=';'.join(notes),
                        wt_aa_len=len(wt_aa_seq), mut_aa_len=len(mut_aa_seq)
                    )
                    return (site_rows, events, summary), mutation, None
                except Exception as e:
                    print(f"Warning: netMHC mutant {gene_name} {mutation} failed: {e}")
                    return None, mutation, f"MUT_EXCEPTION:{type(e).__name__}:{e}"

            # map() preserves input order (deterministic output) and runs concurrently.
            if not gene_failure:
                with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as pool:
                    for result, mutation, reason in pool.map(_process_mutant, list(mutant_seqs.items())):
                        if result is None:
                            mutation_failures[mutation] = reason
                            continue
                        site_rows, events, summary = result
                        gene_sites.extend(site_rows)
                        gene_events.extend(events)
                        gene_summary_rows.append(summary)

        finally:
            # Clean up temp directory
            if os.path.exists(gene_workdir):
                shutil.rmtree(gene_workdir, ignore_errors=True)

        if gene_failure:
            run_summary["genes_skipped"] += 1
            run_summary["skipped_genes"].append({"gene": gene_name, "reason": gene_failure})

        # Every pending token that produced no summary row gets one, with the
        # metric columns empty and the reason named. Silence here is the failure
        # mode this accounting exists to remove.
        scored = {row['mutation'] for row in gene_summary_rows}
        # Denominator and the missing list are both taken over the UNION of the
        # tokens read from the file and the mutants build_mutant_sequences_for_gene
        # actually emitted. Taking either alone leaves a gap: a token the two read
        # differently would either be reported twice or vanish without a row.
        run_summary["mutations_total"] += len(set(pending_tokens) | scored | set(mutation_failures))
        missing = [t for t in pending_tokens if t not in scored]
        missing += [t for t in mutation_failures
                    if t not in scored and t not in pending_tokens]
        for token in missing:
            # Classify the token in every case. variant_class and aa_consequence
            # are properties of the RECORD, not of whether netMHC happened to
            # run: blanking them on a tool failure loses information that is
            # sitting in the token, and leaves the row unusable for the
            # per-class accounting the other rows support.
            fallback_reason, variant_class, aa_consequence = explain_missing_token(
                token, wt_nt_seq, wt_aa_seq, is_nt=(build_input_type == 'nt'))
            reason = gene_failure or mutation_failures.get(token) or fallback_reason
            gene_summary_rows.append(
                blank_summary_row(gene_name, token, reason, variant_class, aa_consequence))
            run_summary["unsuccessful"].append(
                {"gene": gene_name, "mutation": token, "reason": reason})

        gene_out = Path(args.output) / gene_name / "NetMHC"
        gene_out.mkdir(parents=True, exist_ok=True)
        write_summary_tsv(gene_summary_rows, str(gene_out / f"{gene_name}.tsv"))
        write_events_tsv(gene_events, str(gene_out / f"{gene_name}.events.tsv"))
        write_sites_tsv(gene_sites, str(gene_out / f"{gene_name}.sites.tsv"))

    run_summary["genes_processed"] = run_summary["genes_total"] - run_summary["genes_skipped"]
    run_summary["mutations_unsuccessful"] = len(run_summary["unsuccessful"])
    run_summary["mutations_successful"] = (
        run_summary["mutations_total"] - run_summary["mutations_unsuccessful"])
    Path(args.output).mkdir(parents=True, exist_ok=True)
    with open(Path(args.output) / "netmhc.run_summary.json", "w") as handle:
        json.dump(run_summary, handle, indent=2)
    print(f"[SUMMARY] genes {run_summary['genes_processed']}/{run_summary['genes_total']}, "
          f"mutations {run_summary['mutations_successful']}/{run_summary['mutations_total']} ok, "
          f"{run_summary['mutations_unsuccessful']} unsuccessful -> netmhc.run_summary.json")

    if args.verbose:
        print(f"\nPipeline complete!")

    # Cleanup temp holder if used
    if temp_holder:
        temp_holder.cleanup()

    return 0


if __name__ == '__main__':
    sys.exit(main())
