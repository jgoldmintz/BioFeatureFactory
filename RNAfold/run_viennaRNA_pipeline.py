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

'''
Code adapted from boltzmann_sampling.py example
in the ViennaRNA repo: https://github.com/ViennaRNA/ViennaRNA
'''

import os
import sys
import math
import argparse
import statistics
import concurrent.futures
import multiprocessing
from pathlib import Path
import RNA

os.environ["OMP_NUM_THREADS"] = "1"

from utils.utility import (
    _collect_failures_from_logs,
    read_fasta,
    load_mapping,
    update_str,
    subseq,
    validate_fasta_content,
    extract_gene_from_filename,
    trim_muts,
    get_mutation_data,
    get_mutation_data_bioAccurate,
)

R = 1.98717e-3  # kcal/mol/K


def _task(pkey, transcript_pos, seq_ref, seq_alt, samples, tau):
    res = compare_ref_alt(seq_ref, seq_alt, samples=samples)
    summary = (
        pkey, transcript_pos,
        res["ddg_mfe_kcalmol"], res["ddg_ensemble_kcalmol"], res["d_meanE_kcalmol"],
        res["ref_sdE_kcalmol"], res["alt_sdE_kcalmol"], res["JSD_unpaired_bits"],
        res["delta_central"]
    )
    ref_struct = res["ref_mfe_struct"]
    alt_struct = res["alt_mfe_struct"]
    rows = []
    for i, du in enumerate(res["d_unpaired_profile"]):
        pos1 = i + 1
        rdu = round(du, 3)
        change_flag = 1 if abs(rdu) >= tau else 0
        direction = 1 if rdu > 0 else (-1 if rdu < 0 else 0)
        ref_c = ref_struct[i] if i < len(ref_struct) else '.'
        alt_c = alt_struct[i] if i < len(alt_struct) else '.'
        ref_paired = 1 if ref_c in '()' else 0
        alt_paired = 1 if alt_c in '()' else 0
        mfe_change_flag = 1 if ref_paired != alt_paired else 0
        mfe_change_dir = 0 if not mfe_change_flag else (0 if (ref_paired == 0 and alt_paired == 1) else 1)
        rows.append((pkey, transcript_pos, pos1, rdu, change_flag, direction, mfe_change_flag, mfe_change_dir))
    return summary, rows

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

def compare_ref_alt(seq_ref: str, seq_alt: str, samples: int = 1000):
    assert len(seq_ref) == len(seq_alt), "Windows must be same length"
    ref = analyze_seq(seq_ref, samples=samples)
    alt = analyze_seq(seq_alt, samples=samples)

    ddg_mfe = alt["mfe_G"] - ref["mfe_G"]
    ddg_ens = alt["ensemble_G"] - ref["ensemble_G"]
    d_meanE = alt["mean_sampled_E"] - ref["mean_sampled_E"]
    jsd_u = jsd_unpaired(ref["unpaired_prob"], alt["unpaired_prob"])

    delta_u = [a - r for r, a in zip(ref["unpaired_prob"], alt["unpaired_prob"])]
    central_idx = len(delta_u) // 2
    central_delta = delta_u[central_idx]

    return {
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
    parser.add_argument("-o", "--output", required=True, help="Output TSV (summary table)")
    parser.add_argument("-w", "--window", type=int, default=151, help="Window size (odd; truncates near ends)")
    parser.add_argument("--transcript-mapping", help="Path to transcript mapping file/directory")
    parser.add_argument("--log", help="Validation log (file or dir) used to filter failed mutations")
    parser.add_argument("--samples", type=int, default=1000, help="Number of Boltzmann samples per sequence")
    parser.add_argument("--tau", type=float, default=0.05, help="Threshold for change_flag on deltau")
    parser.add_argument("--positions-out", default=None, help="Per-position output TSV (default: <output>.positions.tsv)")
    parser.add_argument("--workers", type=int, default=None, help="Max parallel workers (processes)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist", file=sys.stderr)
        sys.exit(1)

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    positions_out = args.positions_out or (args.output + ".positions.tsv")

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

    need_sum_header = not os.path.exists(args.output) or os.path.getsize(args.output) == 0
    need_pos_header = not os.path.exists(positions_out) or os.path.getsize(positions_out) == 0
    if need_sum_header:
        with open(args.output, "a") as f:
            f.write("pkey\ttranscript_pos\tddg_mfe_kcalmol\tddg_ensemble_kcalmol\td_meanE_kcalmol\tref_sdE_kcalmol\talt_sdE_kcalmol\tjsd_unpaired_bits\tdelta_central\n")
    if need_pos_header:
        with open(positions_out, "a") as fpos:
            fpos.write("pkey\ttranscript_pos\tpos\tdelta_u\tchange_flag\tdirection\tmfe_change_flag\tmfe_change_dir\n")

    work_items = []

    for file in files:
        if not validate_fasta_content(file):
            print(f'fasta structure not valid for {file.name}')
            continue
        gene_name = extract_gene_from_filename(str(file))
        fasta_dict = read_fasta(str(file))
        transcript_seq = next(iter(fasta_dict.values()))

        if len(maps) == 1:
            transcript_map = maps[0]
        else:
            cands = maps_by_gene.get(gene_name, [])
            if not cands:
                print(f'No transcript map matched {gene_name}', file=sys.stderr)
                continue
            if len(cands) > 1:
                print(f'WARN multiple maps for {gene_name}: {[str(c) for c in cands]} -- picking first',file=sys.stderr)
            transcript_map = cands[0]

        mutants = trim_muts(transcript_map, args.log, gene_name)

        for item in mutants:
            toks = [t.strip() for t in str(item).split(",")]
            if len(toks) != 2:
                print(f"Skipping malformed mutation entry: {item}", file=sys.stderr)
                continue
            reference_mut, mapped_mut = toks

            transcript_pos, _ = get_mutation_data_bioAccurate(mapped_mut)
            pkey = f"{gene_name}-{reference_mut}"

            pos, [_, alt] = get_mutation_data(mapped_mut)

            transcript_mutseq = update_str(transcript_seq, alt, pos)
            seq_ref = subseq(transcript_seq, pos, args.window)
            seq_alt = subseq(transcript_mutseq, pos, args.window)

            work_items.append((pkey, transcript_pos, seq_ref, seq_alt, args.samples, args.tau))

    n_tasks = len(work_items)
    if n_tasks == 0:
        return
    max_workers = args.workers or _autodetect_workers(n_tasks)


    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_task, *args_) for args_ in work_items]
        with open(args.output, "a") as fs, open(positions_out, "a") as fp:
            for fut in concurrent.futures.as_completed(futures):
                summary, rows = fut.result()
                fs.write("\t".join(map(str, summary)) + "\n")
                for r in rows:
                    fp.write("\t".join(map(str, r)) + "\n")

if __name__ == "__main__":
    main()
