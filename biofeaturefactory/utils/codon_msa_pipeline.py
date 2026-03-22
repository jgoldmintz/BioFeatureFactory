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
BioFeatureFactory: Codon-Aware MSA Generator

Generates codon-aware multiple sequence alignments from protein alignments
and corresponding nucleotide sequences.

Workflow:
  1. Load protein MSA (from MUSCLE/MAFFT/etc.)
  2. Load corresponding nucleotide sequences for each sequence ID
  3. Back-translate protein alignment to codons
  4. Output codon MSA preserving alignment structure

Output format: FASTA with codon sequences (gaps as '---')
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
import json
from pathlib import Path

from biofeaturefactory.utils.utility import (
    read_fasta,
    write_fasta,
    codon_to_aa,
    extract_gene_from_filename,
    compute_neff,
)


def _translate_codon(codon):
    """Translate a codon to amino acid."""
    if codon == '---' or '-' in codon:
        return '-'
    aa = codon_to_aa.get(codon.upper(), 'X')
    return '*' if aa == 'Stop' else aa



def _backtranslate_with_qc(protein_seq, nucleotide_seq, seq_id=None):
    """
    Back-translate sequence and compute observed translation identity.

    Returns:
        tuple:
            - codon_seq (str)
            - observed_seqid (float or None)
            - compared_positions (int)
            - remaining_nt (int)
            - mismatch_count (int)

    Raises:
        ValueError: If nucleotide sequence is too short for the protein residues.
    """
    codon_seq = []
    nt_pos = 0
    compared_positions = 0
    match_count = 0
    mismatch_count = 0

    for i, aa in enumerate(protein_seq):
        if aa in '-.':
            codon_seq.append('---')
            continue

        if nt_pos + 3 > len(nucleotide_seq):
            raise ValueError(
                f"Nucleotide sequence too short for protein at position {i} "
                f"(seq: {seq_id or 'unknown'})"
            )

        codon = nucleotide_seq[nt_pos:nt_pos + 3]
        translated_aa = _translate_codon(codon)
        codon_seq.append(codon)
        nt_pos += 3

        if aa.upper() == 'X':
            continue
        compared_positions += 1
        if translated_aa.upper() == aa.upper():
            match_count += 1
        else:
            mismatch_count += 1

    remaining_nt = len(nucleotide_seq) - nt_pos
    observed_seqid = (match_count / compared_positions) if compared_positions > 0 else None
    return ''.join(codon_seq), observed_seqid, compared_positions, remaining_nt, mismatch_count


def _iter_fasta_records(path):
    """Yield (header, sequence) from FASTA (plain or gz)."""
    p = Path(path)
    if p.suffix != '.gz':
        for header, seq in read_fasta(str(p)).items():
            yield header, seq
        return

    with gzip.open(p, 'rt') as handle:
        header = None
        seq_chunks = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, ''.join(seq_chunks)


def _normalize_protein_id(pid):
    """Normalize protein accession by removing version suffix."""
    if not pid:
        return ''
    return pid.split('.')[0]


def _extract_protein_ids(text):
    """
    Extract protein-like accessions from an arbitrary identifier string.

    Examples matched:
      NP_000509.1, XP_123456.7, WP_012345678.1, YP_009724390.1
    """
    if not text:
        return []
    if text.startswith('UniRef90_'):
        core = text.split('UniRef90_', 1)[1]
        core = core.split('/', 1)[0]
        core = core.split()[0]
        return [core]
    matches = re.findall(r'([A-Z]{2}_[0-9]+(?:\.[0-9]+)?)', text)
    if matches:
        return matches
    token = text.split()[0]
    return [token]


def _extract_gene_symbol_hint(text):
    """
    Extract a likely gene symbol from sequence IDs like 'SMN2'.

    Returns uppercase symbol or None if the ID looks like an accession-based header.
    """
    if not text:
        return None
    token = text.split()[0]
    if token.startswith('UniRef90_'):
        return None
    if '_' in token and re.match(r'^[A-Z]{2}_[0-9]+(?:\.[0-9]+)?$', token):
        return None
    if re.match(r'^[A-Za-z0-9\-\.]+$', token):
        return token.upper()
    return None


def _discover_cds_files(assembly=None, assembly_dir=None):
    """Discover RefSeq CDS FASTA files from an assembly path or assembly directory."""
    cds_files = []
    if assembly:
        ap = Path(assembly)
        if ap.is_file():
            cds_files.append(ap)
        elif ap.is_dir():
            cds_files.extend(sorted(ap.glob('*_cds_from_genomic.fna.gz')))
            cds_files.extend(sorted(ap.glob('*_cds_from_genomic.fna')))
    if assembly_dir:
        ad = Path(assembly_dir)
        if ad.is_dir():
            cds_files.extend(sorted(ad.rglob('*_cds_from_genomic.fna.gz')))
            cds_files.extend(sorted(ad.rglob('*_cds_from_genomic.fna')))
    uniq = []
    seen = set()
    for fp in cds_files:
        s = str(fp.resolve())
        if s not in seen:
            seen.add(s)
            uniq.append(fp)
    return uniq



def _resolve_nt_from_assemblies(seq_ids, assembly=None, assembly_dir=None, id_map=None):
    """
    Resolve nucleotide CDS sequences for protein MSA IDs from RefSeq assemblies.

    Returns:
        tuple:
          - matched: dict seq_id -> [candidate dicts]
                       candidate fields include:
                       {'nucleotide_seq', 'nucleotide_seq_id', 'source_assembly',
                        'matched_protein_id', 'lookup_protein_id', 'raw_extracted_id'}
          - unresolved_reason: dict seq_id -> reason
    """
    id_map = id_map or {}
    candidates_by_seq = {}
    gene_hint_by_seq = {}
    wanted_gene_hints = set()
    target_norm_ids = set()
    lookup_by_seq = {}
    for sid in seq_ids:
        cands = _extract_protein_ids(sid)
        candidates_by_seq[sid] = cands
        gene_hint = _extract_gene_symbol_hint(sid)
        gene_hint_by_seq[sid] = gene_hint
        if gene_hint:
            wanted_gene_hints.add(gene_hint)
        lookups = []
        for c in cands:
            c_norm = _normalize_protein_id(c)
            lookup_norms = id_map.get(c_norm, [c_norm])
            lookups.append((c, lookup_norms))
            for lookup_norm in lookup_norms:
                target_norm_ids.add(lookup_norm)
        lookup_by_seq[sid] = lookups

    cds_files = _discover_cds_files(assembly=assembly, assembly_dir=assembly_dir)
    if not cds_files:
        raise ValueError("No CDS files found in refseq_assemblies/ under --db-root")

    cds_by_norm = {}
    cds_by_gene = {}
    protein_id_re = re.compile(r'\[protein_id=([^\]]+)\]')
    gene_re = re.compile(r'\[gene=([^\]]+)\]')
    for cds_file in cds_files:
        assembly_name = cds_file.name.replace('_cds_from_genomic.fna.gz', '').replace('_cds_from_genomic.fna', '')
        for header, seq in _iter_fasta_records(cds_file):
            m = protein_id_re.search(header)
            if not m:
                continue
            pid = m.group(1).strip()
            pid_norm = _normalize_protein_id(pid)
            g = gene_re.search(header)
            gene_sym = g.group(1).strip().upper() if g else ''
            keep_for_pid = pid_norm in target_norm_ids
            keep_for_gene = bool(gene_sym and gene_sym in wanted_gene_hints)
            if not keep_for_pid and not keep_for_gene:
                continue
            rec = {
                    'nucleotide_seq': seq.upper(),
                    'nucleotide_seq_id': pid,
                    'source_assembly': assembly_name,
            }
            rec_key = (rec['nucleotide_seq_id'], rec['source_assembly'], rec['nucleotide_seq'])
            if keep_for_pid:
                current = cds_by_norm.get(pid_norm, [])
                existing_keys = {
                    (x['nucleotide_seq_id'], x['source_assembly'], x['nucleotide_seq']) for x in current
                }
                if rec_key not in existing_keys:
                    current.append(rec)
                cds_by_norm[pid_norm] = current
            if keep_for_gene:
                gcur = cds_by_gene.get(gene_sym, [])
                gkeys = {
                    (x['nucleotide_seq_id'], x['source_assembly'], x['nucleotide_seq']) for x in gcur
                }
                if rec_key not in gkeys:
                    gcur.append(rec)
                cds_by_gene[gene_sym] = gcur

    matched = {}
    unresolved_reason = {}
    refseq_like_re = re.compile(r'^[A-Z]{2}_[0-9]+(?:\.[0-9]+)?$')
    for sid in seq_ids:
        cands = candidates_by_seq.get(sid, [])
        lookups = lookup_by_seq.get(sid, [])
        if not cands:
            unresolved_reason[sid] = "No protein accession detected in sequence ID"
            continue
        candidates = []
        for raw_cand, lookup_norms in lookups:
            for lookup_norm in lookup_norms:
                for pick in cds_by_norm.get(lookup_norm, []):
                    candidates.append({
                        'nucleotide_seq': pick['nucleotide_seq'],
                        'nucleotide_seq_id': pick['nucleotide_seq_id'],
                        'source_assembly': pick['source_assembly'],
                        'matched_protein_id': pick['nucleotide_seq_id'],
                        'lookup_protein_id': lookup_norm,
                        'raw_extracted_id': raw_cand,
                    })

        if not candidates:
            gene_hint = gene_hint_by_seq.get(sid)
            if gene_hint and gene_hint in cds_by_gene:
                for pick in cds_by_gene.get(gene_hint, []):
                    candidates.append({
                        'nucleotide_seq': pick['nucleotide_seq'],
                        'nucleotide_seq_id': pick['nucleotide_seq_id'],
                        'source_assembly': pick['source_assembly'],
                        'matched_protein_id': pick['nucleotide_seq_id'],
                        'lookup_protein_id': gene_hint,
                        'raw_extracted_id': gene_hint,
                    })

        if not candidates:
            if not any(refseq_like_re.match(c) for c in cands):
                unresolved_reason[sid] = (
                    "Protein ID appears non-RefSeq (e.g., UniProt/UPI from UniRef90); "
                    "no direct match to RefSeq [protein_id=...] in CDS files"
                )
            elif gene_hint_by_seq.get(sid):
                unresolved_reason[sid] = (
                    f"No CDS entries with [gene={gene_hint_by_seq.get(sid)}] found in assembly files"
                )
            else:
                unresolved_reason[sid] = "No matching protein_id found in CDS assembly files"
            continue
        # Keep stable ordering but remove duplicate candidates.
        uniq = []
        seen = set()
        for c in candidates:
            key = (c['lookup_protein_id'], c['nucleotide_seq_id'], c['source_assembly'])
            if key in seen:
                continue
            seen.add(key)
            uniq.append(c)
        matched[sid] = uniq

    return matched, unresolved_reason



def _translate_nt_to_aa(nt_seq):
    """Translate full nucleotide sequence (in-frame from position 0) to amino acids."""
    aas = []
    n_codons = len(nt_seq) // 3
    for i in range(n_codons):
        codon = nt_seq[i * 3:(i + 1) * 3]
        aas.append(_translate_codon(codon))
    return ''.join(aas)


def _find_best_coding_offset(protein_seq_aligned, nt_seq):
    """
    Find the best amino-acid start offset in candidate CDS translation.

    This compares the ungapped protein query against all same-length windows in the
    translated CDS and picks the offset with maximal identity.

    Returns:
      (best_aa_offset, best_window_seqid, compared_positions)
    """
    query = ''.join(aa for aa in protein_seq_aligned if aa not in '-.')
    if not query:
        return 0, 0.0, 0

    candidate_aa = _translate_nt_to_aa(nt_seq)
    qlen = len(query)
    if len(candidate_aa) < qlen:
        return 0, 0.0, 0

    best = None
    for start in range(0, len(candidate_aa) - qlen + 1):
        compared = 0
        matches = 0
        window = candidate_aa[start:start + qlen]
        for qaa, caa in zip(query, window):
            if qaa.upper() == 'X':
                continue
            compared += 1
            if qaa.upper() == caa.upper():
                matches += 1
        seqid = (matches / compared) if compared > 0 else 0.0
        score = (seqid, compared, matches)
        if best is None or score > best['score']:
            best = {
                'score': score,
                'start': start,
                'seqid': seqid,
                'compared': compared,
            }

    if best is None:
        return 0, 0.0, 0
    return best['start'], best['seqid'], best['compared']



def _load_focus_nt_from_fasta(focus_nt_fasta):
    """
    Load focus ORF nucleotide sequence from a FASTA produced by exon_aware_mapping.

    Preference:
      1) Record named exactly ORF (case-insensitive)
      2) First FASTA record

    Returns:
      tuple: (selected_header, sequence_upper)
    """
    seqs = read_fasta(focus_nt_fasta)
    if not seqs:
        raise ValueError(f"No sequences found in focus nt FASTA: {focus_nt_fasta}")

    for header, seq in seqs.items():
        if header.strip().upper() == 'ORF':
            return header, seq.upper()

    first_header = next(iter(seqs))
    return first_header, seqs[first_header].upper()


def _vprint(msg, verbose):
    """Verbose print with immediate flush."""
    if verbose:
        print(msg, flush=True)


def _write_manifest(manifest_rows, manifest_path):
    """Write protein->nucleotide mapping and back-translation status manifest."""
    fieldnames = [
        'protein_seq_id',
        'raw_extracted_id',
        'lookup_protein_id',
        'nucleotide_seq_id',
        'status',
        'reason',
        'source_assembly',
        'matched_protein_id',
        'min_seqid_threshold',
        'observed_seqid',
        'compared_positions',
        'translation_mismatch_count',
        'protein_alignment_length',
        'protein_residue_count',
        'nucleotide_length',
        'nucleotide_codon_count',
        'expected_coding_nt',
        'remaining_nt_after_backtranslation',
    ]
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote manifest to {manifest_path}")


def _discover_protein_files(assembly_dir):
    """Find protein FAA files recursively in assembly_dir."""
    ad = Path(assembly_dir)
    files = list(ad.rglob('*_protein.faa.gz')) + list(ad.rglob('*_protein.faa'))
    uniq = []
    seen = set()
    for fp in files:
        s = str(fp.resolve())
        if s not in seen:
            seen.add(s)
            uniq.append(fp)
    return sorted(uniq)


def _mmseqs_db_exists(db_base):
    """Check if an mmseqs2 database exists at db_base."""
    return Path(f"{db_base}.dbtype").exists()


def _mmseqs_index_exists(db_base):
    """Check if an mmseqs2 on-disk index exists at db_base."""
    return Path(f"{db_base}.idx.index").exists()


def _build_or_reuse_merged_refseq_proteins(source_dir, merged_faa_path, verbose=False):
    """
    Build merged RefSeq protein FAA from all protein files in source_dir,
    or reuse if already present and non-empty.

    Returns merged_faa_path.
    """
    merged = Path(merged_faa_path)
    if merged.exists() and merged.stat().st_size > 0:
        _vprint(f"[mmseqs] Reusing existing merged FAA: {merged_faa_path}", verbose)
        return merged_faa_path

    protein_files = _discover_protein_files(source_dir)
    if not protein_files:
        raise ValueError(f"No protein FAA files found in {source_dir}")

    _vprint(f"[mmseqs] Building merged FAA from {len(protein_files)} protein files...", verbose)
    with open(merged, 'w') as out_f:
        for pf in protein_files:
            if str(pf).endswith('.gz'):
                with gzip.open(pf, 'rt') as f:
                    for line in f:
                        out_f.write(line)
            else:
                with open(pf) as f:
                    for line in f:
                        out_f.write(line)
    _vprint(f"[mmseqs] Wrote merged FAA to {merged_faa_path}", verbose)
    return merged_faa_path


def _default_mmseqs_target_paths(target_db_base, db_root):
    """
    Determine default mmseqs target DB base path and merged FAA path.

    Returns (target_db_base, merged_faa_path).
    """
    if target_db_base:
        merged_faa = str(Path(target_db_base).parent / 'refseq_proteins_merged.faa')
        return target_db_base, merged_faa

    if db_root:
        db_dir = Path(db_root)
        db_base = str(db_dir / 'refseq_proteins_db')
        merged_faa = str(db_dir / 'refseq_proteins_merged.faa')
        return db_base, merged_faa

    raise ValueError("Must provide --db-root or --target-db-base")


def _search_focus_protein_against_refseq(
    focus_protein_seq, focus_seq_id,
    db_root=None,
    mmseqs_binary='mmseqs', min_seqid=0.5, min_qcov=0.5,
    max_hits=500, threads=None, target_db_base=None,
    mmseqs_tmp_dir=None, verbose=False,
):
    """
    Search a single focus protein against the local RefSeq protein database using mmseqs2.

    Returns a list of RefSeq protein accessions (strings) for all hits passing the
    identity/coverage thresholds, up to max_hits.
    """
    target_db_base, merged_faa_path = _default_mmseqs_target_paths(target_db_base, db_root)

    merged = Path(merged_faa_path)
    if not (merged.exists() and merged.stat().st_size > 0):
        raise ValueError(
            f"refseq_proteins_merged.faa not found at {merged_faa_path}. "
            "Run scripts/build_db.sh first."
        )

    work_dir = tempfile.mkdtemp(prefix='codon_msa_mmseqs_', dir=mmseqs_tmp_dir)
    try:
        query_fasta = os.path.join(work_dir, 'query.faa')
        query_db    = os.path.join(work_dir, 'queryDB')
        result_db   = os.path.join(work_dir, 'resultDB')
        result_tsv  = os.path.join(work_dir, 'results.tsv')
        tmp_search  = os.path.join(work_dir, 'tmp')
        os.makedirs(tmp_search, exist_ok=True)

        with open(query_fasta, 'w') as f:
            f.write(f">{focus_seq_id}\n{focus_protein_seq}\n")

        devnull = subprocess.DEVNULL
        out_kw = {} if verbose else {'stdout': devnull, 'stderr': devnull}
        threads_args = ['--threads', str(threads)] if threads else []

        # Build target DB and index if needed.
        if not _mmseqs_db_exists(target_db_base):
            _vprint(f"[mmseqs] Building target DB at {target_db_base}...", verbose)
            subprocess.run(
                [mmseqs_binary, 'createdb', merged_faa_path, target_db_base],
                check=True, **out_kw,
            )
        else:
            _vprint(f"[mmseqs] Reusing existing target DB: {target_db_base}", verbose)

        if not _mmseqs_index_exists(target_db_base):
            _vprint(f"[mmseqs] Building on-disk index at {target_db_base}.idx...", verbose)
            idx_tmp = target_db_base + '_idx_tmp'
            os.makedirs(idx_tmp, exist_ok=True)
            try:
                subprocess.run(
                    [mmseqs_binary, 'createindex', target_db_base, idx_tmp]
                    + threads_args,
                    check=True, **out_kw,
                )
            finally:
                shutil.rmtree(idx_tmp, ignore_errors=True)
        else:
            _vprint(f"[mmseqs] Reusing existing on-disk index: {target_db_base}.idx", verbose)

        # Create query DB.
        subprocess.run(
            [mmseqs_binary, 'createdb', query_fasta, query_db],
            check=True, **out_kw,
        )

        # Run search.
        search_cmd = [
            mmseqs_binary, 'search',
            query_db, target_db_base, result_db, tmp_search,
            '--min-seq-id', str(min_seqid),
            '-c', str(min_qcov),
            '--max-seqs', str(max_hits),
            '--alignment-mode', '3',
        ] + threads_args
        _vprint(f"[mmseqs] Running: {' '.join(search_cmd)}", verbose)
        subprocess.run(search_cmd, check=True, **out_kw)

        # Convert to TSV.
        subprocess.run(
            [
                mmseqs_binary, 'convertalis',
                query_db, target_db_base, result_db, result_tsv,
                '--format-output', 'query,target,fident,alnlen,qcov,tcov',
            ],
            check=True, **out_kw,
        )

        # Parse hits.
        accessions = []
        seen_acc = set()
        if os.path.exists(result_tsv):
            with open(result_tsv) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 2:
                        continue
                    # mmseqs uses the first space-delimited token of the FASTA header as ID.
                    acc = parts[1].split()[0]
                    acc_norm = _normalize_protein_id(acc)
                    if acc_norm not in seen_acc:
                        seen_acc.add(acc_norm)
                        accessions.append(acc)

        _vprint(f"[mmseqs] Found {len(accessions)} unique hits", verbose)
        return accessions

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _align_protein_sequences(sequences, aligner='mafft'):
    """
    Align {seq_id: ungapped_protein_seq} using mafft --auto (or muscle).
    Returns {seq_id: aligned_protein_seq}.
    """
    if not sequences:
        return {}
    if len(sequences) == 1:
        return dict(sequences)

    work_dir = tempfile.mkdtemp(prefix='codon_msa_align_')
    try:
        in_fasta  = os.path.join(work_dir, 'input.faa')
        out_fasta = os.path.join(work_dir, 'aligned.faa')

        # Use numeric safe IDs to avoid special-character issues with aligners.
        id_fwd = {}  # orig_id -> safe_id
        id_rev = {}  # safe_id -> orig_id
        with open(in_fasta, 'w') as f:
            for i, (sid, seq) in enumerate(sequences.items()):
                safe = f"seq{i}"
                id_fwd[sid]  = safe
                id_rev[safe] = sid
                f.write(f">{safe}\n{seq}\n")

        if aligner == 'mafft':
            result = subprocess.run(
                ['mafft', '--auto', '--quiet', in_fasta],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            with open(out_fasta, 'w') as f:
                f.write(result.stdout)
        elif aligner == 'muscle':
            subprocess.run(
                ['muscle', '-align', in_fasta, '-output', out_fasta],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            raise ValueError(f"Unknown aligner: {aligner!r}")

        aligned = {}
        for header, seq in _iter_fasta_records(out_fasta):
            safe_id = header.split()[0]
            orig_id = id_rev.get(safe_id)
            if orig_id:
                aligned[orig_id] = seq.upper()
        return aligned

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def generate_codon_msa_from_focus(
    focus_nt_fasta,
    db_root=None,
    output_dir=None,
    min_seqid=0.5, min_qcov=0.5,
    min_backtranslation_seqid=1.0,
    max_hits=500,
    mmseqs_binary='mmseqs', aligner='mafft',
    threads=None, target_db_base=None,
    mmseqs_tmp_dir=None, verbose=False,
):
    """
    Generate a codon-aware MSA starting from a focus CDS (nt FASTA).

    Workflow:
      1. Load focus CDS from focus_nt_fasta and translate to protein.
      2. Search focus protein against RefSeq proteins with mmseqs2.
      3. Resolve CDS from assemblies for each hit accession.
      4. Translate each CDS to protein.
      5. Align focus protein + all translated hits with mafft/muscle.
      6. Back-translate each aligned row using its corresponding CDS.
      7. Write codon MSA FASTA and manifest TSV.

    Returns:
        tuple: (codon_msa dict, manifest_rows list)
    """
    # 1. Load focus CDS and translate to protein.
    focus_nt_header, focus_nt_seq = _load_focus_nt_from_fasta(focus_nt_fasta)
    _vprint(f"[focus-mode] Focus CDS: {focus_nt_header} ({len(focus_nt_seq)} nt)", verbose)

    focus_seq_id = focus_nt_header
    focus_protein_ungapped = _translate_nt_to_aa(focus_nt_seq)
    if focus_protein_ungapped.endswith('*'):
        focus_protein_ungapped = focus_protein_ungapped[:-1]
    _vprint(f"[focus-mode] Focus protein: {focus_seq_id} ({len(focus_protein_ungapped)} aa)", verbose)

    # 3. Search focus protein against RefSeq.
    _vprint(f"[focus-mode] Searching RefSeq (max_hits={max_hits}, min_seqid={min_seqid}, min_qcov={min_qcov})...", verbose)
    hit_accessions = _search_focus_protein_against_refseq(
        focus_protein_seq=focus_protein_ungapped,
        focus_seq_id=focus_seq_id,
        db_root=db_root,
        mmseqs_binary=mmseqs_binary,
        min_seqid=min_seqid,
        min_qcov=min_qcov,
        max_hits=max_hits,
        threads=threads,
        target_db_base=target_db_base,
        mmseqs_tmp_dir=mmseqs_tmp_dir,
        verbose=verbose,
    )
    _vprint(f"[focus-mode] mmseqs returned {len(hit_accessions)} hits", verbose)
    if not hit_accessions:
        raise ValueError("No RefSeq hits found for focus protein; check database and thresholds")

    # 4. Resolve CDS from assemblies for each hit accession.
    _vprint("[focus-mode] Resolving CDS from assemblies...", verbose)
    assembly_dir = str(Path(db_root) / 'refseq_assemblies') if db_root else None
    assembly_matches, assembly_unresolved = _resolve_nt_from_assemblies(
        seq_ids=hit_accessions,
        assembly_dir=assembly_dir,
    )
    _vprint(
        f"[focus-mode] CDS resolved for {len(assembly_matches)} / {len(hit_accessions)} hits",
        verbose,
    )

    # 5. Translate each resolved CDS to protein (strip stop codon if present).
    hit_proteins  = {}  # acc -> ungapped translated protein
    hit_cds       = {}  # acc -> CDS nt sequence (full, including stop if present)
    hit_assembly  = {}  # acc -> source_assembly name
    hit_nt_id     = {}  # acc -> nucleotide_seq_id from CDS file
    hit_cand_meta = {}  # acc -> candidate dict

    for acc, candidates in assembly_matches.items():
        if not candidates:
            continue
        cand = candidates[0]
        nt_seq = cand['nucleotide_seq']
        translated = _translate_nt_to_aa(nt_seq)
        if translated.endswith('*'):
            translated = translated[:-1]
        if '*' in translated:
            continue  # internal stop codon — corrupted/frameshifted CDS
        hit_proteins[acc]  = translated
        hit_cds[acc]       = nt_seq
        hit_assembly[acc]  = cand.get('source_assembly', '')
        hit_nt_id[acc]     = cand.get('nucleotide_seq_id', '')
        hit_cand_meta[acc] = cand

    _vprint(f"[focus-mode] Translated {len(hit_proteins)} CDS sequences", verbose)

    # 6. Align focus protein + all translated hits.
    seqs_to_align = {focus_seq_id: focus_protein_ungapped}
    seqs_to_align.update(hit_proteins)
    _vprint(f"[focus-mode] Aligning {len(seqs_to_align)} sequences with {aligner}...", verbose)
    aligned_proteins = _align_protein_sequences(seqs_to_align, aligner=aligner)
    if not aligned_proteins:
        raise ValueError("Protein alignment produced no output")

    # 7. Back-translate each row.
    codon_msa     = {}
    manifest_rows = []
    failed        = []

    # Focus sequence first.
    focus_aligned = aligned_proteins.get(focus_seq_id, focus_protein_ungapped)
    focus_residue_count = len(focus_protein_ungapped)
    try:
        f_codon, f_seqid, f_compared, f_remaining, f_mismatches = _backtranslate_with_qc(
            focus_aligned, focus_nt_seq, focus_seq_id
        )
        codon_msa[focus_seq_id] = f_codon
        f_status = 'PASS'
        f_reason = ''
        if f_remaining > 3:
            f_status = 'PASS_WITH_EXTRA_TAIL_NT'
            f_reason = f'{f_remaining} unused nucleotide(s) at sequence tail'
    except ValueError as e:
        failed.append(focus_seq_id)
        f_status, f_reason = 'FAIL_BACKTRANSLATION', str(e)
        f_seqid, f_compared, f_remaining, f_mismatches = None, 0, 0, 0
        print(f"Warning: Focus backtranslation failed: {e}", file=sys.stderr)

    manifest_rows.append({
        'protein_seq_id':                 focus_seq_id,
        'raw_extracted_id':               focus_seq_id,
        'lookup_protein_id':              focus_seq_id,
        'nucleotide_seq_id':              focus_nt_header,
        'status':                         f_status,
        'reason':                         f_reason,
        'source_assembly':                'focus_nt_fasta',
        'matched_protein_id':             focus_nt_header,
        'min_seqid_threshold':            min_backtranslation_seqid,
        'observed_seqid':                 f"{f_seqid:.6f}" if f_seqid is not None else '',
        'compared_positions':             f_compared,
        'translation_mismatch_count':     f_mismatches,
        'protein_alignment_length':       len(focus_aligned),
        'protein_residue_count':          focus_residue_count,
        'nucleotide_length':              len(focus_nt_seq),
        'nucleotide_codon_count':         len(focus_nt_seq) // 3,
        'expected_coding_nt':             focus_residue_count * 3,
        'remaining_nt_after_backtranslation': f_remaining,
    })

    total_hits = len(hit_proteins)
    for idx, acc in enumerate(hit_proteins, start=1):
        _vprint(f"[focus-mode seq {idx}/{total_hits}] {acc}", verbose)

        nt_seq   = hit_cds[acc]
        prot_seq = aligned_proteins.get(acc)

        if prot_seq is None:
            failed.append(acc)
            pr_count = len(hit_proteins[acc])
            manifest_rows.append({
                'protein_seq_id':                 acc,
                'raw_extracted_id':               acc,
                'lookup_protein_id':              acc,
                'nucleotide_seq_id':              hit_nt_id.get(acc, ''),
                'status':                         'FAIL_NO_ALIGNMENT',
                'reason':                         'Sequence absent from alignment output',
                'source_assembly':                hit_assembly.get(acc, ''),
                'matched_protein_id':             hit_nt_id.get(acc, ''),
                'min_seqid_threshold':            min_backtranslation_seqid,
                'observed_seqid':                 '',
                'compared_positions':             '',
                'translation_mismatch_count':     '',
                'protein_alignment_length':       '',
                'protein_residue_count':          pr_count,
                'nucleotide_length':              len(nt_seq),
                'nucleotide_codon_count':         len(nt_seq) // 3,
                'expected_coding_nt':             pr_count * 3,
                'remaining_nt_after_backtranslation': '',
            })
            _vprint(f"[focus-mode seq {idx}/{total_hits}] FAIL_NO_ALIGNMENT", verbose)
            continue

        protein_residue_count = sum(1 for aa in prot_seq if aa not in '-.')
        expected_coding_nt    = protein_residue_count * 3

        # Isoform-aware offset (should be 0 for direct RefSeq CDS, but harmless).
        aa_offset, _, _ = _find_best_coding_offset(prot_seq, nt_seq)
        nt_seq_adj = nt_seq[aa_offset * 3:]

        try:
            codon_seq, observed_seqid, compared_positions, remaining_nt, mismatch_count = (
                _backtranslate_with_qc(prot_seq, nt_seq_adj, acc)
            )
        except ValueError as e:
            failed.append(acc)
            manifest_rows.append({
                'protein_seq_id':                 acc,
                'raw_extracted_id':               acc,
                'lookup_protein_id':              acc,
                'nucleotide_seq_id':              hit_nt_id.get(acc, ''),
                'status':                         'FAIL_BACKTRANSLATION',
                'reason':                         str(e),
                'source_assembly':                hit_assembly.get(acc, ''),
                'matched_protein_id':             hit_nt_id.get(acc, ''),
                'min_seqid_threshold':            min_backtranslation_seqid,
                'observed_seqid':                 '',
                'compared_positions':             '',
                'translation_mismatch_count':     '',
                'protein_alignment_length':       len(prot_seq),
                'protein_residue_count':          protein_residue_count,
                'nucleotide_length':              len(nt_seq),
                'nucleotide_codon_count':         len(nt_seq) // 3,
                'expected_coding_nt':             expected_coding_nt,
                'remaining_nt_after_backtranslation': '',
            })
            _vprint(f"[focus-mode seq {idx}/{total_hits}] FAIL_BACKTRANSLATION ({e})", verbose)
            continue

        status = 'PASS'
        reason = ''
        if observed_seqid is None:
            status = 'FAIL_NO_COMPARABLE_POSITIONS'
            reason = 'No comparable protein positions (all gaps/X)'
        elif observed_seqid < min_backtranslation_seqid:
            status = 'FAIL_LOW_SEQID'
            reason = (
                f'Observed seqid {observed_seqid:.4f} below threshold '
                f'{min_backtranslation_seqid:.4f} '
                f'({mismatch_count}/{compared_positions} mismatches)'
            )
        elif remaining_nt > 3:
            status = 'PASS_WITH_EXTRA_TAIL_NT'
            reason = f'{remaining_nt} unused nucleotide(s) at sequence tail'

        if status.startswith('PASS'):
            codon_msa[acc] = codon_seq
        else:
            failed.append(acc)

        manifest_rows.append({
            'protein_seq_id':                 acc,
            'raw_extracted_id':               acc,
            'lookup_protein_id':              acc,
            'nucleotide_seq_id':              hit_nt_id.get(acc, ''),
            'status':                         status,
            'reason':                         reason,
            'source_assembly':                hit_assembly.get(acc, ''),
            'matched_protein_id':             hit_nt_id.get(acc, ''),
            'min_seqid_threshold':            min_backtranslation_seqid,
            'observed_seqid':                 f"{observed_seqid:.6f}" if observed_seqid is not None else '',
            'compared_positions':             compared_positions,
            'translation_mismatch_count':     mismatch_count,
            'protein_alignment_length':       len(prot_seq),
            'protein_residue_count':          protein_residue_count,
            'nucleotide_length':              len(nt_seq),
            'nucleotide_codon_count':         len(nt_seq) // 3,
            'expected_coding_nt':             expected_coding_nt,
            'remaining_nt_after_backtranslation': remaining_nt,
        })
        _vprint(f"[focus-mode seq {idx}/{total_hits}] {status}", verbose)

    print(f"Successfully processed {len(codon_msa)} sequences ({len(failed)} failed)")

    # 8. Write output.
    gene = extract_gene_from_filename(Path(focus_nt_fasta).stem) or Path(focus_nt_fasta).stem
    out_dir = Path(output_dir) / gene / "CodonMSA" if output_dir else Path(focus_nt_fasta).parent / gene / "CodonMSA"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"{gene}.codon.msa.fasta"
    manifest_path = out_dir / f"{gene}.codon.msa.manifest.tsv"

    write_fasta(output_path, codon_msa)
    print(f"Wrote codon MSA to {output_path}")
    _write_manifest(manifest_rows, str(manifest_path))

    # Compute and write stats JSON (mirrors msa_generation_pipeline.py output)
    focus_id = gene
    query_length = len(focus_nt_seq) // 3  # codon positions
    n_seqs = len(codon_msa)
    neff = compute_neff(codon_msa, identity_threshold=0.8)
    neff_ratio = neff / query_length if query_length > 0 else 0.0

    # Mean identity to focus
    focus_seq = codon_msa.get(focus_id, '')
    identities = []
    for seq_id, seq in codon_msa.items():
        if seq_id != focus_id and len(seq) == len(focus_seq):
            matches = sum(1 for a, b in zip(seq, focus_seq) if a == b and a not in '-.')
            aligned = sum(1 for a, b in zip(seq, focus_seq) if a not in '-.' and b not in '-.')
            if aligned > 0:
                identities.append(matches / aligned)
    mean_identity = float(np.mean(identities)) if identities else 0.0

    stats = {
        'query_id': focus_id,
        'query_length_codons': query_length,
        'n_sequences': n_seqs,
        'n_eff': float(round(neff, 1)),
        'n_eff_ratio': float(round(neff_ratio, 2)),
        'mean_identity': float(round(mean_identity, 3)),
        
    }

    stats_path = out_dir / f"{gene}.codon.msa.stats.json"
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Statistics written to {stats_path}")

    return codon_msa, manifest_rows


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Codon-Aware MSA Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python codon_msa_pipeline.py \\
    --fasta SMN2.orf.fasta \\
    --db-root /path/to/Bio_DBs/ \\
    --output results/ \\
    --max-hits 500 \\
    --search-min-seqid 0.5 \\
    --min-seqid 1.0
  # Writes: results/SMN2/CodonMSA/SMN2.codon.msa.fasta
  #         results/SMN2/CodonMSA/SMN2.codon.msa.manifest.tsv
"""
    )

    parser.add_argument('--fasta', required=True,
                        help='FASTA with focus ORF nt sequence (exon_aware_mapping output)')

    parser.add_argument('--db-root', required=True,
                        help='Bio_DBs root directory produced by scripts/build_db.sh '
                             '(contains refseq_assemblies/ and refseq_proteins_merged.faa)')

    parser.add_argument('--output', '-o', required=True,
                        help='Output base directory (writes {GENE}/CodonMSA/{GENE}.codon.msa.fasta and .manifest.tsv)')
    parser.add_argument('--min-seqid', type=float, default=1.0,
                        help='Minimum observed AA identity for back-translation QC (default: 1.0)')
    parser.add_argument('--search-min-seqid', type=float, default=0.5,
                        help='Minimum sequence identity for mmseqs2 search (default: 0.5)')
    parser.add_argument('--search-min-qcov', type=float, default=0.5,
                        help='Minimum query coverage for mmseqs2 search (default: 0.5)')
    parser.add_argument('--max-hits', type=int, default=500,
                        help='Max RefSeq hits to retrieve from mmseqs2 (default: 500)')
    parser.add_argument('--aligner', default='mafft', choices=['mafft', 'muscle'],
                        help='Protein aligner (default: mafft)')
    parser.add_argument('--mmseqs-binary', default='mmseqs',
                        help='Path to mmseqs2 binary (default: mmseqs)')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of threads for mmseqs2 (default: auto)')
    parser.add_argument('--target-db-base', default=None,
                        help='Base path for mmseqs2 target DB (default: auto from --db-root)')
    parser.add_argument('--mmseqs-tmp-dir', default=None,
                        help='Temp directory for mmseqs2 (default: system temp)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress updates')

    args = parser.parse_args()

    if args.min_seqid < 0 or args.min_seqid > 1:
        parser.error("--min-seqid must be between 0 and 1")

    generate_codon_msa_from_focus(
        focus_nt_fasta=args.fasta,
        db_root=args.db_root,
        output_dir=args.output,
        min_seqid=args.search_min_seqid,
        min_qcov=args.search_min_qcov,
        min_backtranslation_seqid=args.min_seqid,
        max_hits=args.max_hits,
        mmseqs_binary=args.mmseqs_binary,
        aligner=args.aligner,
        threads=args.threads,
        target_db_base=args.target_db_base,
        mmseqs_tmp_dir=args.mmseqs_tmp_dir,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
