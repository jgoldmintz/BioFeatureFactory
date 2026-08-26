# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# Licensed under the PolyForm Noncommercial License 1.0.0.
# You may use this software for any purpose other than commercial use.
# See LICENSE, or https://polyformproject.org/licenses/noncommercial/1.0.0/
#
# Required Notice: Copyright (C) 2023-2026 Jacob Goldmintz
# (https://github.com/jgoldmintz)

'''
Code adapted from boltzmann_sampling.py example
in the ViennaRNA repo: https://github.com/ViennaRNA/ViennaRNA
'''

import os
import sys
import math
import argparse
import json
import statistics
import concurrent.futures
import multiprocessing
from pathlib import Path
import RNA

os.environ["OMP_NUM_THREADS"] = "1"

from biofeaturefactory.utils.utility import (
    _collect_failures_from_logs,
    read_fasta,
    update_str,
    subseq,
    validate_fasta_content,
    extract_gene_from_filename,
    trim_muts,
    get_mutation_data,
    get_mutation_data_bioAccurate,
    detect_alphabet,
)

R = 1.98717e-3  # kcal/mol/K


SUMMARY_HEADER = "pkey\ttranscript_pos\tref_mfe_G\talt_mfe_G\tddg_mfe_kcalmol\tddg_ensemble_kcalmol\td_meanE_kcalmol\tref_sdE_kcalmol\talt_sdE_kcalmol\tjsd_unpaired_bits\tdelta_central"
POS_HEADER = "pkey\ttranscript_pos\tpos\ttx_pos\tdelta_u\tchange_flag\tdirection\tmfe_change_flag\tmfe_change_dir"


def _write_rnafold_tsv(path, header, rows):
    with open(path, 'w') as f:
        f.write(header + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def _task(pkey, transcript_pos, seq_ref, seq_alt, offset, samples, tau, gene_name):
    res = compare_ref_alt(seq_ref, seq_alt, offset=offset, samples=samples)
    summary = (
        pkey, transcript_pos,
        res["ref_mfe_G"], res["alt_mfe_G"],
        res["ddg_mfe_kcalmol"], res["ddg_ensemble_kcalmol"], res["d_meanE_kcalmol"],
        res["ref_sdE_kcalmol"], res["alt_sdE_kcalmol"], res["JSD_unpaired_bits"],
        res["delta_central"]
    )
    ref_struct = res["ref_mfe_struct"]
    alt_struct = res["alt_mfe_struct"]
    rows = []
    for i, du in enumerate(res["d_unpaired_profile"]):
        pos1 = i + 1
        tx_pos_i = transcript_pos - offset + i   # transcript coordinate of this window base
        rdu = round(du, 3)
        change_flag = 1 if abs(rdu) >= tau else 0
        direction = 1 if rdu > 0 else (-1 if rdu < 0 else 0)
        ref_c = ref_struct[i] if i < len(ref_struct) else '.'
        alt_c = alt_struct[i] if i < len(alt_struct) else '.'
        ref_paired = 1 if ref_c in '()' else 0
        alt_paired = 1 if alt_c in '()' else 0
        mfe_change_flag = 1 if ref_paired != alt_paired else 0
        mfe_change_dir = 0 if not mfe_change_flag else (0 if (ref_paired == 0 and alt_paired == 1) else 1)
        rows.append((pkey, transcript_pos, pos1, tx_pos_i, rdu, change_flag, direction, mfe_change_flag, mfe_change_dir))
    return summary, rows, gene_name

def analyze_seq(seq: str, samples: int = 1000):
    md = RNA.md()
    md.temperature = 37.0
    md.dangles = 2
    md.noLP = True
    md.uniq_ML = 1

    fc = RNA.fold_compound(seq, md)

    mfe_struct, mfe_G = fc.mfe()

    fc.exp_params_rescale(mfe_G)
    _, ens_G = fc.pf()  # ensemble free energy (kcal/mol)

    energies = []
    unpaired = [0] * len(seq)

    for s in fc.pbacktrack(samples):
        E = fc.eval_structure(s)
        energies.append(E)
        for i, c in enumerate(s):
            if c == '.':
                unpaired[i] += 1

    mean_E = statistics.fmean(energies)
    sd_E = statistics.pstdev(energies)
    u_prob = [c / float(samples) for c in unpaired]

    return {
        "mfe_G": mfe_G,
        "ensemble_G": ens_G,
        "mfe_struct": mfe_struct,
        "mean_sampled_E": mean_E,
        "sd_sampled_E": sd_E,
        "unpaired_prob": u_prob,
    }

def jsd_unpaired(p, q, eps=1e-12):
    def H(pi):
        a = pi + eps
        b = 1.0 - pi + eps
        return -(a * math.log2(a) + b * math.log2(b))
    m = [(pi + qi) / 2.0 for pi, qi in zip(p, q)]
    return sum(H(mi) - 0.5 * H(pi) - 0.5 * H(qi) for mi, pi, qi in zip(m, p, q)) / len(p)


def window_around(seq: str, pos: int, l: int) -> tuple[str, int]:
    """Return a length-l window centered on pos and pos's 0-based index within it.

    Near a terminus the shortfall is shifted to the opposite side so the window
    stays length l (F21); if the sequence is shorter than l the whole sequence is
    returned. Unlike subseq (symmetric truncation), this keeps windows the same
    length across variants and reports where the variant actually sits, so the
    caller can index the variant's own base instead of assuming the midpoint.
    """
    assert isinstance(l, int) and l > 0 and (l % 2 == 1), "l must be a positive odd integer"
    assert 0 <= pos < len(seq), "pos out of range"
    half = l // 2
    start = pos - half
    end = pos + half + 1
    if start < 0:            # 5' shortfall -> extend the 3' side
        end -= start
        start = 0
    if end > len(seq):       # 3' shortfall -> extend the 5' side
        start -= (end - len(seq))
        end = len(seq)
    start = max(0, start)    # sequence shorter than l: clamp to whole seq
    return seq[start:end], pos - start

def compare_ref_alt(seq_ref: str, seq_alt: str, offset=None, samples: int = 1000):
    assert len(seq_ref) == len(seq_alt), "Windows must be same length"
    ref = analyze_seq(seq_ref, samples=samples)
    alt = analyze_seq(seq_alt, samples=samples)

    ddg_mfe = alt["mfe_G"] - ref["mfe_G"]
    ddg_ens = alt["ensemble_G"] - ref["ensemble_G"]
    d_meanE = alt["mean_sampled_E"] - ref["mean_sampled_E"]
    jsd_u = jsd_unpaired(ref["unpaired_prob"], alt["unpaired_prob"])

    delta_u = [a - r for r, a in zip(ref["unpaired_prob"], alt["unpaired_prob"])]
    # F21: report the variant's own base. `offset` is its 0-based index in the
    # window (from window_around); fall back to the midpoint for callers that
    # don't pass it. Clamp guards ORFs shorter than the window.
    central_idx = offset if offset is not None else len(delta_u) // 2
    central_idx = max(0, min(central_idx, len(delta_u) - 1))
    central_delta = delta_u[central_idx]

    return {
        "ref_mfe_G": ref["mfe_G"],
        "alt_mfe_G": alt["mfe_G"],
        "ddg_mfe_kcalmol": ddg_mfe,
        "ddg_ensemble_kcalmol": ddg_ens,
        "d_meanE_kcalmol": d_meanE,
        "ref_sdE_kcalmol": ref["sd_sampled_E"],
        "alt_sdE_kcalmol": alt["sd_sampled_E"],
        "JSD_unpaired_bits": jsd_u,
        "d_unpaired_profile": delta_u,
        "delta_central": central_delta,
        "ref_mfe_struct": ref["mfe_struct"],
        "alt_mfe_struct": alt["mfe_struct"],
    }

def _autodetect_workers(n_tasks: int, cap: int = 8) -> int:
    n_cpu = os.cpu_count() or multiprocessing.cpu_count() or 1
    return max(1, min(n_cpu // 2 if n_cpu > 1 else 1, cap, n_tasks))

def main():
    parser = argparse.ArgumentParser(
        description="Compute ddG, JSD, and per-position deltau for variant-centered RNA folding windows using ViennaRNA."
    )
    parser.add_argument("-i", "--input", required=True, help="Input fasta sequence file/dir")
    parser.add_argument("-o", "--output", required=True, help="Output base directory")
    parser.add_argument("-w", "--window", type=int, default=151, help="Window size (odd; truncates near ends)")
    parser.add_argument("--transcript-mapping", help="Path to transcript mapping file/directory")
    parser.add_argument("--log", help="Validation log (file or dir) used to filter failed mutations")
    parser.add_argument("--samples", type=int, default=1000, help="Number of Boltzmann samples per sequence")
    parser.add_argument("--tau", type=float, default=0.05, help="Threshold for change_flag on deltau")
    parser.add_argument("--workers", type=int, default=None, help="Max parallel workers (processes)")
    args = parser.parse_args()

    if args.window < 1 or args.window % 2 == 0:
        print(f"Error: --window must be a positive odd integer (got {args.window})", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
        sys.exit(1)

    input_path = Path(args.input)
    transcpt = Path(args.transcript_mapping) if args.transcript_mapping else None
    valid_exts_fa = ('.fasta', '.fa', '.faa', '.fna')
    valid_exts_map = ('.tsv', '.csv')

    files = [input_path] if input_path.is_file() else [i for i in input_path.iterdir() if i.suffix in valid_exts_fa]
    if not files:
        raise ValueError(f'No valid input files found in {input_path.name}')
    if transcpt is None:
        raise ValueError("--transcript-mapping is required")

    maps = [transcpt] if transcpt.is_file() else [t for t in transcpt.iterdir() if t.suffix in valid_exts_map]
    if not maps:
        raise ValueError(f'No valid transcript maps found in {transcpt.name if hasattr(transcpt, "name") else transcpt}')
    maps_by_gene = {}
    for m in maps:
        g = extract_gene_from_filename(str(m))
        if g:
            maps_by_gene.setdefault(g,[]).append(m)

    work_items = []
    run_summary = {
        "genes_total": len(files), "genes_processed": 0, "genes_skipped": 0,
        "mutations_total": 0, "mutations_successful": 0, "mutations_unsuccessful": 0,
        "skipped_genes": [], "unsuccessful": [],
    }

    for file in files:
        gene_name = extract_gene_from_filename(str(file))
        if not validate_fasta_content(file):
            print(f'fasta structure not valid for {file.name}')
            run_summary["genes_skipped"] += 1
            run_summary["skipped_genes"].append({"gene": gene_name or file.name, "reason": "invalid_fasta"})
            continue
        fasta_dict = read_fasta(str(file))
        if "transcript" in fasta_dict:
            transcript_seq = fasta_dict["transcript"]
        else:
            transcript_seq = next(iter(fasta_dict.values()))

        # F22: reject protein/codon FASTA at the gene level (validate_fasta_content is
        # shared with protein pipelines and cannot be tightened). detect_alphabet raises
        # on an empty/gap-only sequence -> treat as non-nucleotide.
        try:
            is_nt = detect_alphabet(transcript_seq) == "nucleotide"
        except ValueError:
            is_nt = False
        if not is_nt:
            print(f"Skipping {gene_name}: transcript is not a nucleotide sequence", file=sys.stderr)
            run_summary["genes_skipped"] += 1
            run_summary["skipped_genes"].append({"gene": gene_name or file.name, "reason": "gene_not_nucleotide"})
            continue

        if len(maps) == 1:
            transcript_map = maps[0]
        else:
            cands = maps_by_gene.get(gene_name, [])
            if not cands:
                print(f'No transcript map matched {gene_name}', file=sys.stderr)
                run_summary["genes_skipped"] += 1
                run_summary["skipped_genes"].append({"gene": gene_name, "reason": "no_transcript_map"})
                continue
            if len(cands) > 1:
                print(f'WARN multiple maps for {gene_name}: {[str(c) for c in cands]} -- picking first',file=sys.stderr)
            transcript_map = cands[0]

        mutants = trim_muts(transcript_map, args.log, gene_name)

        for item in mutants:
            run_summary["mutations_total"] += 1
            toks = [t.strip() for t in str(item).split(",")]
            if len(toks) != 2:
                print(f"Skipping malformed mutation entry: {item}", file=sys.stderr)
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": None, "mapped_mut": str(item), "reason": "malformed_token"})
                continue
            reference_mut, mapped_mut = toks
            pkey = f"{gene_name}-{reference_mut}"

            # F24: an aa/consequence token ('Stop'/'Sto') has no nt position to fold;
            # bioAccurate returns None for it. A nt nonsense variant (e.g. C175T) never
            # matches and folds normally. Unparseable tokens are logged, not crashed on.
            try:
                transcript_pos, _ = get_mutation_data_bioAccurate(mapped_mut, is_nt=True)
            except (ValueError, IndexError):
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": "unparseable_token"})
                continue
            if transcript_pos is None:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": "aa_level_token"})
                continue

            pos, (ref, alt) = get_mutation_data(mapped_mut)
            if not (0 <= pos < len(transcript_seq)):   # F25: out-of-range coord -> skip this mutation, not the whole gene
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": "position_out_of_range"})
                continue
            if ref and transcript_seq[pos].upper() != ref.upper():   # F26: mapping/transcript disagree -> phantom mutation
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": "wt_ref_mismatch"})
                continue

            transcript_mutseq = update_str(transcript_seq, alt, pos)
            seq_ref, offset = window_around(transcript_seq, pos, args.window)
            seq_alt, _ = window_around(transcript_mutseq, pos, args.window)

            work_items.append((pkey, transcript_pos, seq_ref, seq_alt, offset, args.samples, args.tau, gene_name))

    results_by_gene = {}
    if work_items:
        max_workers = args.workers or _autodetect_workers(len(work_items))
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(_task, *w): (w[0], w[7]) for w in work_items}  # w[0]=pkey, w[7]=gene_name
            for fut in concurrent.futures.as_completed(futures):
                pkey_f, gene_f = futures[fut]
                try:
                    summary, rows, gene_name = fut.result()
                except Exception:
                    run_summary["unsuccessful"].append({
                        "gene": gene_f,
                        "reference_mut": pkey_f.split("-", 1)[1] if "-" in pkey_f else None,
                        "mapped_mut": None, "reason": "fold_error"})
                    continue
                results_by_gene.setdefault(gene_name, {'summary': [], 'positions': []})
                results_by_gene[gene_name]['summary'].append(summary)
                results_by_gene[gene_name]['positions'].extend(rows)

    for gname, data in results_by_gene.items():
        out_dir = Path(args.output) / gname / "RNAfold"
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_rnafold_tsv(out_dir / f"{gname}.rnafold.tsv", SUMMARY_HEADER, data['summary'])
        _write_rnafold_tsv(out_dir / f"{gname}.rnafold.positions.tsv", POS_HEADER, data['positions'])

    # F25: run-level accounting written regardless of how many tasks ran
    run_summary["genes_processed"] = run_summary["genes_total"] - run_summary["genes_skipped"]
    run_summary["mutations_unsuccessful"] = len(run_summary["unsuccessful"])
    run_summary["mutations_successful"] = run_summary["mutations_total"] - run_summary["mutations_unsuccessful"]
    Path(args.output).mkdir(parents=True, exist_ok=True)
    with open(Path(args.output) / "rnafold.run_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)
    print(f"[SUMMARY] genes {run_summary['genes_processed']}/{run_summary['genes_total']}, "
          f"mutations {run_summary['mutations_successful']}/{run_summary['mutations_total']} ok, "
          f"{run_summary['mutations_unsuccessful']} unsuccessful -> rnafold.run_summary.json")

if __name__ == "__main__":
    main()
