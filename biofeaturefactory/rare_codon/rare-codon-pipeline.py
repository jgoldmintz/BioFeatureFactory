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
BioFeatureFactory: Rare Codon Enrichment Pipeline

Wrapper for cg_cotrans library (Copyright William M. Jacobs, GPL v3).
Detects regions enriched/depleted in rare codons using sliding window analysis.

Requires:
  - cg_cotrans library: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
  - Codon-aware MSA (FASTA)
  - Pre-computed codon usage file (.p.gz from calc_codon_usage.py)

Output format: TSV with BFF-standard pkey column
"""

import argparse
import csv
import gzip
import os
import pickle
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
CG_DIR = SCRIPT_DIR / "cg_cotrans"

from biofeaturefactory.utils.utility import (
    read_fasta,
    trim_muts,
    get_mutation_data_bioAccurate,
    extract_gene_from_filename,
    load_validation_failures,
    should_skip_mutation,
)

# Import cg_cotrans library (GPL v3 licensed, Copyright 2017 William M. Jacobs)
# Download from: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis

try:
    from calc_rare_enrichment import (
        read_fasta as rc_read_fasta,
        clean_sequences,
        sorted_gis,
        get_codons,
        align_sequences,
        aa_identity,
        load_null_model,
        msa_rare_codon_analysis_wtalign_nseq,
    )
    from codons import codon_to_aa
except ImportError:
    if not CG_DIR.exists():
        raise FileNotFoundError(f"cg_cotrans not found: {CG_DIR} you must download cg_cotrans from https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis")
    py_files = sorted(CG_DIR.glob("*.py"))
    if not py_files:
        raise FileNotFoundError(f"No *.py files found in {CG_DIR}")
      # copy into pipeline script directory (parent of cg_cotrans)
    import subprocess
    subprocess.run(
        ["cp", "-f", *[str(p) for p in py_files], str(SCRIPT_DIR)],
        check=True,
    )
    import importlib
    importlib.invalidate_caches()

    from calc_rare_enrichment import (
        read_fasta as rc_read_fasta,
        clean_sequences,
        sorted_gis,
        get_codons,
        align_sequences,
        aa_identity,
        load_null_model,
        msa_rare_codon_analysis_wtalign_nseq,
    )
    from codons import codon_to_aa


FIELDNAMES = [
    'pkey',
    'Gene',
    'codon_position',
    'p_enriched',
    'p_depleted',
    'f_enriched_wt',
    'frac_seq_enriched',
    'frac_seq_depleted',
    'n_rare',
    'window_size',
    'qc_flags',
]


def run_rare_codon_analysis(gene, msa_path, usage_path, wt_gi, window_size=15,
                            rare_model='no_norm', rare_threshold=0.1,
                            null_model='genome', max_len_diff=0.2, min_aa_iden=0.5):
    """
    Run rare codon enrichment analysis on an MSA.

    Args:
        gene: Gene name
        msa_path: Path to codon-aware MSA FASTA
        usage_path: Path to pre-computed codon usage file (.p.gz)
        wt_gi: Identifier for the wildtype/focus sequence
        window_size: Sliding window width in codons
        rare_model: 'no_norm' or 'cmax_norm'
        rare_threshold: Threshold for defining rare codons
        null_model: 'genome', 'eq', or 'groups'
        max_len_diff: Max relative sequence length difference vs WT
        min_aa_iden: Min amino acid identity vs WT

    Returns:
        dict: Position -> {p_enriched, p_depleted, f_enriched_wt, etc.}
    """
    # Load and clean MSA; sanitize ambiguous codons before clean_sequences
    seqs = rc_read_fasta(msa_path)
    seqs = {gi: ''.join(c if (c == '---' or c in codon_to_aa) else '---'
                        for c in (seq[i:i+3] for i in range(0, len(seq), 3)))
            for gi, seq in seqs.items()}
    seqs = clean_sequences(seqs)

    if wt_gi not in seqs:
        #raise ValueError(f"WT sequence '{wt_gi}' not found in MSA")
        wt_gi = next(iter(seqs))
    # Filter sequences by length
    wt_len = len(seqs[wt_gi]) - seqs[wt_gi].count('-')
    for gi in list(seqs.keys()):
        if gi == wt_gi:
            continue
        gi_len = len(seqs[gi]) - seqs[gi].count('-')
        wt_gi_overlap = sum(1 for i in range(len(seqs[wt_gi]))
                           if seqs[wt_gi][i] != '-' and seqs[gi][i] != '-')
        if abs(1. - gi_len / wt_len) > max_len_diff or \
           1. - wt_gi_overlap / wt_len > max_len_diff:
            del seqs[gi]

    # Filter by AA identity
    gis = sorted_gis(seqs, wt_gi)
    msa_codons = {gi: get_codons(seq) for gi, seq in seqs.items()}
    aa_perc_id = aa_identity(msa_codons, wt_gi)
    for gi in list(seqs.keys()):
        if gi == wt_gi:
            continue
        val = aa_perc_id.get(gi, 0)
        # np.nan < threshold is False; use "not >=" to treat nan as failing
        if not (val >= min_aa_iden):
            del seqs[gi]

    # Rebuild after filtering
    gis = sorted_gis(seqs, wt_gi)
    msa_codons = {gi: get_codons(seq) for gi, seq in seqs.items()}

    # Remove sequences with no valid sense codons after position 30 (nstart used
    # by gene_avg_codon_probabilities) to prevent division by zero there.
    _NSTART = 30
    msa_codons = {gi: codons for gi, codons in msa_codons.items()
                  if gi == wt_gi or
                  sum(1 for c in codons[_NSTART:]
                      if c != '---' and c in codon_to_aa and codon_to_aa[c] != 'Stop') > 0}
    gis = sorted_gis(msa_codons, wt_gi)

    msa_index_dict = align_sequences(msa_codons, wt_gi)

    # Load null model
    rare_codons, gene_codon_usage, rare_codon_prob, _ = load_null_model(
        msa_codons, gis, usage_path, rare_model, False,
        rare_threshold, null_model, gene
    )

    # Run analysis
    rc_analysis = msa_rare_codon_analysis_wtalign_nseq(
        msa_codons, wt_gi, msa_index_dict, rare_codons, rare_codon_prob,
        L=window_size
    )

    # Extract per-position data
    results = {}
    for pos in rc_analysis['p_nseq_enriched'].keys():
        results[pos] = {
            'p_enriched': rc_analysis['p_nseq_enriched'][pos],
            'p_depleted': rc_analysis['p_nseq_depleted'][pos],
            'f_enriched_wt': rc_analysis['f_enriched'][wt_gi].get(pos),
            'frac_seq_enriched': rc_analysis['fseq_enriched'].get(pos),
            'frac_seq_depleted': rc_analysis['fseq_depleted'].get(pos),
            'n_rare': rc_analysis['n_rare'][wt_gi].get(pos, 0),
        }

    return results, len(seqs)


def process_mutations(mutations_list, gene, orf_sequence, rc_results, window_size, failure_map=None):
    """
    Map mutations to rare codon enrichment statistics.

    Args:
        mutations_list: List of mutation strings
        gene: Gene symbol
        orf_sequence: ORF nucleotide sequence
        rc_results: Dict from run_rare_codon_analysis
        window_size: Window size used in analysis
        failure_map: Optional validation failure map

    Returns:
        list: Result dictionaries
    """
    results = []
    failure_map = failure_map or {}

    for ntposnt in mutations_list:
        if should_skip_mutation(gene, ntposnt, failure_map):
            continue

        qc_flags = []
        pkey = f"{gene}-{ntposnt}"

        # Get mutation position
        pos_data = get_mutation_data_bioAccurate(ntposnt)
        if pos_data[0] is None:
            qc_flags.append('INVALID_MUTATION')
            results.append({
                'pkey': pkey,
                'Gene': gene,
                'codon_position': '',
                'p_enriched': '',
                'p_depleted': '',
                'f_enriched_wt': '',
                'frac_seq_enriched': '',
                'frac_seq_depleted': '',
                'n_rare': '',
                'window_size': window_size,
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
            })
            continue

        nt_pos = pos_data[0]
        codon_pos = (nt_pos - 1) // 3 + 1  # 1-based codon position

        # Look up enrichment data (center position of window)
        rc_data = rc_results.get(codon_pos)

        if rc_data is None:
            # Position may be near edges (not covered by window)
            qc_flags.append('POSITION_NOT_IN_WINDOW')
            results.append({
                'pkey': pkey,
                'Gene': gene,
                'codon_position': codon_pos,
                'p_enriched': '',
                'p_depleted': '',
                'f_enriched_wt': '',
                'frac_seq_enriched': '',
                'frac_seq_depleted': '',
                'n_rare': '',
                'window_size': window_size,
                'qc_flags': ';'.join(qc_flags) if qc_flags else 'PASS',
            })
            continue

        results.append({
            'pkey': pkey,
            'Gene': gene,
            'codon_position': codon_pos,
            'p_enriched': rc_data['p_enriched'],
            'p_depleted': rc_data['p_depleted'],
            'f_enriched_wt': rc_data['f_enriched_wt'],
            'frac_seq_enriched': rc_data['frac_seq_enriched'],
            'frac_seq_depleted': rc_data['frac_seq_depleted'],
            'n_rare': rc_data['n_rare'],
            'window_size': window_size,
            'qc_flags': 'PASS',
        })

    return results


def write_output(results, output_path):
    """Write results to TSV file."""
    if not results:
        print("No results to write")
        return

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {output_path}")


def _resolve_per_gene(path_arg, gene, extensions):
    """Return a matching file for *gene* from a file or directory path.

    If path_arg is a file, return it.  If it is a directory, glob for
    ``{gene}.*`` files whose suffix (case-insensitive) is in *extensions*
    and return the first match, or None.
    """
    if path_arg is None:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        for f in sorted(p.iterdir()):
            if f.stem.upper() == gene.upper() and f.suffix.lower() in extensions:
                return str(f)
        # broader search: gene anywhere in stem
        for f in sorted(p.iterdir()):
            if gene.upper() in f.stem.upper() and f.suffix.lower() in extensions:
                return str(f)
    return None


def _collect_msa_inputs(msa_arg):
    """Resolve one or more codon MSA files from --msa."""
    p = Path(msa_arg)
    if p.is_file():
        return [str(p.resolve())]
    if not p.is_dir():
        return []

    patterns = ("*.msa.fasta", "*.fasta", "*.fa")
    files = []
    for pat in patterns:
        files.extend(sorted(p.glob(pat)))
    # stable de-dup preserving order
    seen = set()
    out = []
    for f in files:
        rf = str(f.resolve())
        if rf not in seen:
            seen.add(rf)
            out.append(rf)
    return out


def _ensure_codon_usage(args):
    """Return a usable codon usage .p.gz path, auto-building if missing."""
    if args.usage and Path(args.usage).is_file():
        return args.usage

    if args.usage and not Path(args.usage).exists():
        print(f"Usage file not found: {args.usage}; auto-building codon usage .p.gz", file=sys.stderr)
    elif not args.usage:
        print("--usage not provided; auto-building codon usage .p.gz")

    msa_inputs = _collect_msa_inputs(args.msa)
    if not msa_inputs:
        raise FileNotFoundError(f"Could not resolve MSA inputs from --msa: {args.msa}")

    usage_workdir = Path(args.output).resolve() / "_rare_codon_usage_cache"
    usage_workdir.mkdir(parents=True, exist_ok=True)

    abund_path = usage_workdir / "auto_abundances.tsv"
    # Empty abundances is valid; calc_codon_usage.py falls back to unit weights.
    abund_path.write_text("")

    calc_script = SCRIPT_DIR / "calc_codon_usage.py"
    if not calc_script.exists():
        raise FileNotFoundError(f"Required script missing: {calc_script}")

    cmd = [sys.executable, str(calc_script), *msa_inputs, str(abund_path)]
    # cg_cotrans implementation may require two passes before codon_usage.p.gz exists.
    try:
        for i in range(2):
            print(f"Building codon usage pass {i + 1}/2...")
            subprocess.run(
                cmd,
                cwd=str(usage_workdir),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        err_tail = (e.stderr or "").strip().splitlines()
        err_msg = err_tail[-1] if err_tail else str(e)
        print(
            "Auto-build via calc_codon_usage.py failed; "
            f"falling back to internal usage builder for sparse/synthetic MSAs. ({err_msg})",
            file=sys.stderr,
        )
        usage_path = usage_workdir / "codon_usage.p.gz"
        _build_usage_pgz_from_msa(msa_inputs, usage_path)
        print(f"Using codon usage: {usage_path}")
        return str(usage_path)

    usage_path = usage_workdir / "codon_usage.p.gz"
    if not usage_path.exists():
        # Defensive fallback in case calc script exits 0 but does not emit codon_usage.p.gz.
        _build_usage_pgz_from_msa(msa_inputs, usage_path)
    print(f"Using codon usage: {usage_path}")
    return str(usage_path)


def _build_usage_pgz_from_msa(msa_inputs, usage_path, pseudocount=1.0):
    """Build a minimal codon_usage.p.gz directly from codon MSA(s).

    This fallback is robust for very small/synthetic alignments where the original
    cg_cotrans builder can hit divide-by-zero for absent amino-acid classes.
    """
    aa_codons = defaultdict(list)
    sense_codons = []
    for codon, aa in codon_to_aa.items():
        if aa != "Stop":
            sense_codons.append(codon)
            aa_codons[aa].append(codon)

    # gi -> codon count table
    gi_counts = {}
    gene_groups = {}

    def _init_counts():
        d = defaultdict(float)
        for c in sense_codons:
            d[c] = float(pseudocount)
        return d

    def _sanitize_seqs(seqs):
        """Replace ambiguous codons (not in codon_to_aa) with '---' so that
        clean_sequences does not KeyError on IUPAC ambiguity codes (e.g. CCN)."""
        out = {}
        for gi, seq in seqs.items():
            sanitized = []
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                sanitized.append(codon if (codon == '---' or codon in codon_to_aa) else '---')
            out[gi] = ''.join(sanitized)
        return out

    for msa_path in msa_inputs:
        gene = extract_gene_from_filename(msa_path)
        gene_groups[gene] = 0
        seqs = clean_sequences(_sanitize_seqs(rc_read_fasta(msa_path)))
        for gi, seq in seqs.items():
            if gi not in gi_counts:
                gi_counts[gi] = _init_counts()
            if len(seq) % 3 != 0:
                raise ValueError(f"MSA sequence length must be multiple of 3 for {gi} in {msa_path}")
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                if codon == "---":
                    continue
                if codon in codon_to_aa and codon_to_aa[codon] != "Stop":
                    gi_counts[gi][codon] += 1.0

    if not gi_counts:
        raise ValueError("No usable GI sequences found while building fallback codon usage")

    relative_usage = {}
    for gi, counts in gi_counts.items():
        rel = {}
        for aa, codons in aa_codons.items():
            denom = sum(counts[c] for c in codons)
            if denom <= 0:
                # Should not happen due to pseudocount, but keep stable.
                denom = float(len(codons))
            for c in codons:
                rel[c] = counts[c] / denom
        relative_usage[gi] = rel

    usage_obj = {
        "groups": ["all"],
        "gene_groups": gene_groups,
        "overall_codon_usage": relative_usage,
        "unweighted_codon_usage": relative_usage,
        "gene_group_codon_usage": {gi: {0: relative_usage[gi]} for gi in relative_usage},
    }

    with gzip.open(usage_path, "wb") as f:
        pickle.dump(usage_obj, f)


def _run_single_gene(gene, fasta_file, msa_file, mut_file, args, wt_gi, output_dir):
    """Run rare codon analysis for a single gene and write per-gene output."""
    fasta = read_fasta(fasta_file)

    if 'ORF' not in fasta:
        print(f"Error: No 'ORF' found in {fasta_file}", file=sys.stderr)
        return

    orf_sequence = fasta['ORF']
    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    mut_list = trim_muts(mut_file, log=args.validation_log, gene_name=gene)

    if not mut_list:
        print(f"Error: No mutations found in {mut_file}", file=sys.stderr)
        return

    print(f"Running rare codon enrichment analysis for {gene}...")
    print(f"  MSA: {msa_file}")
    print(f"  Usage: {args.usage}")
    print(f"  Window size: {args.window_size}")

    try:
        rc_results, n_seqs = run_rare_codon_analysis(
            gene=gene,
            msa_path=msa_file,
            usage_path=args.usage,
            wt_gi=wt_gi,
            window_size=args.window_size,
            rare_model=args.rare_model,
            rare_threshold=args.rare_threshold,
            null_model=args.null_model,
            max_len_diff=args.max_len_diff,
            min_aa_iden=args.min_aa_iden,
        )
        print(f"  MSA sequences used: {n_seqs}")
        print(f"  Positions analyzed: {len(rc_results)}")
    except Exception as e:
        print(f"Error in rare codon analysis for {gene}: {e}", file=sys.stderr)
        return

    results = process_mutations(
        mut_list, gene, orf_sequence, rc_results,
        args.window_size, failure_map
    )

    out_dir = Path(output_dir) / gene / "RareCodon"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_output(results, str(out_dir / f"{gene}.rare_codon.tsv"))

    if results:
        n_enriched = sum(1 for r in results if r.get('p_enriched') and float(r['p_enriched']) < 0.05)
        n_depleted = sum(1 for r in results if r.get('p_depleted') and float(r['p_depleted']) < 0.05)
        print(f"\nSummary for {gene}:")
        print(f"  Total mutations: {len(results)}")
        print(f"  In enriched regions (p<0.05): {n_enriched}")
        print(f"  In depleted regions (p<0.05): {n_depleted}")


def main():
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: Rare Codon Enrichment Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single gene
  python rare-codon-pipeline.py \\
    --msa /path/to/codon_msa.fasta \\
    --usage /path/to/codon_usage.p.gz \\
    --wt-gi "focus_sequence_id" \\
    --fasta /path/to/gene.fasta \\
    --mutations /path/to/mutations.csv \\
    --output results/

  # Directory mode
  python rare-codon-pipeline.py \\
    --fasta fastas/ --msa msas/ --usage codon_usage.p.gz \\
    --mutations mutations/ --output results/

Required preprocessing:
  1. Generate codon-aware MSA (use msa/codon-msa-pipeline.py)
  2. Optional: provide --usage codon_usage.p.gz (auto-generated if omitted)

Copyright notice:
  cg_cotrans library is Copyright (C) 2017 William M. Jacobs (GPL v3)
  Source: https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis
"""
    )

    # MSA and codon usage inputs
    parser.add_argument('--msa', required=True,
                        help='Path to codon-aware MSA FASTA file or directory')
    parser.add_argument('--usage', help='Path to codon usage .p.gz file (optional; auto-built if missing)')
    parser.add_argument('--wt-gi', help='Identifier for WT/focus sequence in MSA (single-gene mode)')

    # Mutation inputs
    parser.add_argument('--fasta', required=True,
                        help='FASTA file with ORF sequence, or directory of FASTA files')
    parser.add_argument('--mutations', required=True,
                        help='Mutations CSV file or directory of CSV files')
    parser.add_argument('--validation-log', help='Validation log for filtering')

    # Analysis parameters
    parser.add_argument('--window-size', '-L', type=int, default=15,
                        help='Sliding window width in codons (default: 15)')
    parser.add_argument('--rare-model', choices=['no_norm', 'cmax_norm'], default='no_norm',
                        help='Rare codon definition model (default: no_norm)')
    parser.add_argument('--rare-threshold', type=float, default=0.1,
                        help='Threshold for codon rarity (default: 0.1)')
    parser.add_argument('--null-model', choices=['genome', 'eq', 'groups'], default='genome',
                        help='Null model for codon usage (default: genome)')
    parser.add_argument('--max-len-diff', type=float, default=0.2,
                        help='Max relative length difference vs WT (default: 0.2)')
    parser.add_argument('--min-aa-iden', type=float, default=0.5,
                        help='Min AA identity vs WT (default: 0.5)')

    # Output
    parser.add_argument('--output', '-o', required=True, help='Output base directory')

    args = parser.parse_args()
    args.usage = _ensure_codon_usage(args)

    fasta_path = Path(args.fasta)
    if fasta_path.is_dir():
        fasta_files = sorted(
            f for ext in ('*.fasta', '*.fa', '*.fna')
            for f in fasta_path.glob(ext)
        )
        if not fasta_files:
            print(f"Error: No FASTA files found in {args.fasta}", file=sys.stderr)
            sys.exit(1)
        for fasta_file in fasta_files:
            gene = extract_gene_from_filename(str(fasta_file))
            msa_file = _resolve_per_gene(args.msa, gene, ('.fasta', '.fa', '.msa.fasta'))
            mut_file = _resolve_per_gene(args.mutations, gene, ('.csv',))
            if not msa_file or not mut_file:
                print(f"  Skipping {gene}: missing MSA or mutations")
                continue
            _run_single_gene(gene, str(fasta_file), msa_file, mut_file,
                             args, args.wt_gi or gene, args.output)
    else:
        gene = extract_gene_from_filename(args.fasta)
        _run_single_gene(gene, args.fasta, args.msa, args.mutations,
                         args, args.wt_gi or gene, args.output)


if __name__ == "__main__":
    main()
