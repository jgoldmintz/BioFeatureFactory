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
import csv
import math
import argparse
import json
import statistics
import concurrent.futures
import multiprocessing
from pathlib import Path
import RNA

os.environ["OMP_NUM_THREADS"] = "1"

from biofeaturefactory.lib.utility import (
    discover_mapping_files,
    discover_fasta_files,
    mint_pkey,
    parse_piece_token,
    split_piece_cell,
    _collect_failures_from_logs,
    read_fasta,
    update_str,
    subseq,
    validate_fasta_content,
    extract_gene_from_filename,
    trim_muts,
    get_mutation_data,
    get_mutation_data_bioAccurate,
    get_variant_data,
    get_variant_data_bioAccurate,
    splice_seq,
    align_wt_to_mut,
    aligned_pairs,
    detect_alphabet,
)

R = 1.98717e-3  # kcal/mol/K


# substrate / substrate_pos are APPENDED, never inserted: a consumer reading by
# column index keeps working on the original 12 columns.
#
# substrate names the molecule that was folded:
#   transcript  spliced mRNA -- exonic variants (unchanged behaviour)
#   intron      the excised intron record, transcript orientation
#   pre_mRNA    the full tx_start..tx_end span, transcript orientation
# An intronic variant is folded TWICE, once against each of the latter two. They
# answer different questions: the intron is the lariat as a free molecule, the
# pre-mRNA is the co-transcriptional context in which splice-site selection
# happens. They diverge most for a variant within half a window of an intron
# boundary, where the intron record has no junction context to offer, and least
# for short introns (SMN2 intron 4 is 159 nt against a 151 nt default window).
#
# transcript_pos stays the TRANSCRIPT coordinate and is empty on non-transcript
# rows; substrate_pos is the coordinate within whatever was folded.
SUMMARY_HEADER = "pkey\ttranscript_pos\tref_mfe_G\talt_mfe_G\tddg_mfe_kcalmol\tddg_ensemble_kcalmol\td_meanE_kcalmol\tref_sdE_kcalmol\talt_sdE_kcalmol\tjsd_unpaired_bits\tdelta_central\tqc_flags\tsubstrate\tsubstrate_pos"

# positions.tsv is emitted in the WT frame with one addition: bases the ALT
# INSERTED have no WT coordinate, so they carry pos/tx_pos empty and are marked
# in align_status. Without those rows a 20 nt insertion is invisible here and
# still reports full alignment, because every WT base does have a counterpart.
#
# align_status: aligned | deleted (WT base with no mutant counterpart) |
#               inserted (mutant base with no WT counterpart)
# change_flag stays strictly 0/1 and is empty on non-aligned rows -- it used to
# carry the string "deleted", which broke int() on the column.
POS_HEADER = ("pkey\ttranscript_pos\tpos\ttx_pos\tmut_pos\talign_status"
              "\tdelta_u\tchange_flag\tdirection\tmfe_change_flag\tmfe_change_dir"
              "\tsubstrate\tsubstrate_pos")


def _write_rnafold_tsv(path, header, rows):
    with open(path, 'w') as f:
        f.write(header + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def _task(pkey, transcript_pos, seq_ref, seq_alt, offset, samples, tau, gene_name,
          ref_len=1, alt_len=1, substrate="transcript", substrate_pos=""):
    res = compare_ref_alt(seq_ref, seq_alt, offset=offset, samples=samples,
                          ref_len=ref_len, alt_len=alt_len)
    profile = res["d_unpaired_profile"]
    n_wt = len(profile)
    n_deleted = sum(1 for d in profile if d is None)
    n_inserted = max(0, alt_len - ref_len)
    if ref_len != alt_len:
        # Count over the UNION of both alleles' bases, not just WT positions.
        # aligned_N/N over WT alone reports an insertion as fully aligned however
        # large it is, because every WT base does keep a counterpart -- the bases
        # with no counterpart are all on the mutant side.
        qc = (f"length_changed:{alt_len - ref_len:+d}nt;"
              f"aligned_{n_wt - n_deleted}/{n_wt + n_inserted};"
              f"deleted_{n_deleted};inserted_{n_inserted}")
    else:
        qc = ""
    summary = (
        pkey, transcript_pos,
        res["ref_mfe_G"], res["alt_mfe_G"],
        res["ddg_mfe_kcalmol"], res["ddg_ensemble_kcalmol"], res["d_meanE_kcalmol"],
        res["ref_sdE_kcalmol"], res["alt_sdE_kcalmol"], res["JSD_unpaired_bits"],
        "" if res["delta_central"] is None else res["delta_central"], qc,
        substrate, substrate_pos
    )
    ref_struct = res["ref_mfe_struct"]
    alt_struct = res["alt_mfe_struct"]
    rows = []
    # Rows are emitted in the WT frame: one per WT window base, indexed by WT
    # coordinate. A base the edit deleted has no counterpart to subtract, so its
    # delta columns are empty and change_flag is 'deleted' rather than 0 -- a 0
    # would read as "measured, unchanged".
    mut_index = align_wt_to_mut(n_wt, offset, ref_len, alt_len)
    claimed_mut = set()
    # transcript_pos is "" on the intron and pre-mRNA substrates -- those bases
    # have no transcript coordinate, which is the whole reason they are folded
    # elsewhere. Without this guard the arithmetic below is "" - int and every
    # intronic fold dies as a bare "fold_error" with no indication of why.
    has_tx = isinstance(transcript_pos, int)
    for i, du in enumerate(profile):
        pos1 = i + 1
        tx_pos_i = (transcript_pos - offset + i) if has_tx else ""
        j = mut_index[i]
        if j is not None:
            claimed_mut.add(j)
        if du is None:
            rows.append((pkey, transcript_pos, pos1, tx_pos_i, "", "deleted",
                         "", "", "", "", "", substrate, substrate_pos))
            continue
        rdu = round(du, 3)
        change_flag = 1 if abs(rdu) >= tau else 0
        direction = 1 if rdu > 0 else (-1 if rdu < 0 else 0)
        ref_c = ref_struct[i] if i < len(ref_struct) else '.'
        # Index the ALT structure by the PROJECTED index, not by i. Using i here
        # compares the WT base at window position i to whatever the indel slid
        # into that slot -- the same misalignment the unpaired-probability zip
        # had, in the MFE-pairing columns.
        alt_c = alt_struct[j] if (j is not None and j < len(alt_struct)) else '.'
        ref_paired = 1 if ref_c in '()' else 0
        alt_paired = 1 if alt_c in '()' else 0
        mfe_change_flag = 1 if ref_paired != alt_paired else 0
        mfe_change_dir = 0 if not mfe_change_flag else (0 if (ref_paired == 0 and alt_paired == 1) else 1)
        rows.append((pkey, transcript_pos, pos1, tx_pos_i, (j + 1) if j is not None else "",
                     "aligned", rdu, change_flag, direction, mfe_change_flag, mfe_change_dir,
                     substrate, substrate_pos))

    # Bases the ALT inserted have no WT coordinate and so appear in none of the
    # rows above. Emit them explicitly, keyed on the mutant window position, with
    # the WT columns empty -- otherwise an insertion is entirely invisible in the
    # per-position output.
    for j in range(len(alt_struct)):
        if j in claimed_mut:
            continue
        alt_c = alt_struct[j]
        rows.append((pkey, transcript_pos, "", "", j + 1, "inserted",
                     "", "", "", "", 1 if alt_c in '()' else 0,
                     substrate, substrate_pos))
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

def compare_ref_alt(seq_ref: str, seq_alt: str, offset=None, samples: int = 1000,
                    ref_len=1, alt_len=1):
    """Compare two folded windows covering the SAME biological span.

    ref_len/alt_len are the variant's allele lengths. When they differ the two
    windows differ in length by exactly that amount -- which is correct, because
    both cover the same region of the transcript. The old equal-length assert is
    therefore replaced by a check that the observed difference is the expected
    one; an unexplained difference is still a bug.

    Every per-position quantity is computed over the WT-frame alignment rather
    than a positional zip, so base i of the WT is compared to the base it
    actually became, and deleted bases contribute nothing instead of being
    silently paired with whatever slid into their index.
    """
    expected = alt_len - ref_len
    assert len(seq_alt) - len(seq_ref) == expected, (
        f"window length difference {len(seq_alt) - len(seq_ref)} does not match "
        f"the allele length difference {expected}")
    ref = analyze_seq(seq_ref, samples=samples)
    alt = analyze_seq(seq_alt, samples=samples)

    ddg_mfe = alt["mfe_G"] - ref["mfe_G"]
    ddg_ens = alt["ensemble_G"] - ref["ensemble_G"]
    d_meanE = alt["mean_sampled_E"] - ref["mean_sampled_E"]

    edit_offset = offset if offset is not None else len(seq_ref) // 2
    pairs = list(aligned_pairs(ref["unpaired_prob"], alt["unpaired_prob"],
                               edit_offset, ref_len, alt_len))
    jsd_u = jsd_unpaired([p for _, p, _ in pairs], [q for _, _, q in pairs])

    # WT-frame profile: one entry per WT window position, None where the edit
    # deleted that base. Consumers index this by WT coordinate.
    paired = {i: (p, q) for i, p, q in pairs}
    delta_u = [
        (paired[i][1] - paired[i][0]) if i in paired else None
        for i in range(len(ref["unpaired_prob"]))
    ]
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

def _prepare_fold(seq, coord_token, window, pkey, transcript_pos, gene_name,
                  samples, tau, substrate):
    """Build one fold work item, or return (None, reason).

    Every substrate goes through this single function. The transcript, intron and
    pre-mRNA paths differ ONLY in which sequence and which coordinate token they
    are handed -- the bounds check, the REF guard, the window centring and the
    same-biological-span rule are identical, and duplicating them per substrate is
    how two of the three come to be subtly wrong.

    coord_token is "REF<pos>ALT" in `seq`'s own coordinate space, 1-based.
    """
    pos, refalt = get_variant_data(coord_token)
    if pos is not None:
        ref, alt = refalt
    else:
        # A decomposed piece: 'GTGAGGTCG1del' has no ALT because the single
        # anchor base of the whole deletion belongs to a different piece.
        # parse_variant refuses it by design (Variant rejects an empty allele),
        # but a pure deletion is perfectly foldable -- splice_seq(seq, pos, ref,
        # "") is exactly the WT window minus those bases.
        piece = parse_piece_token(coord_token)
        if piece is None:
            return None, "unparseable_token"
        pos, ref, alt = piece["pos"] - 1, piece["ref"], piece["alt"]

    # Bound on the END of the REF span, not its start: a multi-base REF can begin
    # inside the sequence and run off the end.
    if not (0 <= pos and pos + len(ref) <= len(seq)):
        return None, "position_out_of_range"
    # Compare the WHOLE REF span; checking only seq[pos] would pass on a
    # multi-base REF whose first base happens to match.
    if ref and seq[pos:pos + len(ref)].upper() != ref.upper():
        return None, "wt_ref_mismatch"

    mutseq = splice_seq(seq, pos, ref, alt, validate=False)

    # Centre on the MIDPOINT of the REF span, not its first base -- see the long
    # note on the transcript path for why anchoring on pos skews a multi-base REF.
    centre = min(pos + len(ref) // 2, len(seq) - 1)
    seq_ref, centre_offset = window_around(seq, centre, window)
    offset = centre_offset - (len(ref) // 2)

    # The mutant window must cover the SAME BIOLOGICAL SPAN, not merely the same
    # number of bases -- a second window_around on the mutant would re-centre it.
    win_start = centre - centre_offset
    win_end = win_start + len(seq_ref)
    seq_alt = mutseq[win_start:win_end + (len(alt) - len(ref))]

    if not (0 <= offset and offset + len(ref) <= len(seq_ref)):
        return None, (f"ref_span_outside_window:offset_{offset}_len_{len(ref)}"
                      f"_of_{len(seq_ref)}")

    return (pkey, transcript_pos, seq_ref, seq_alt, offset, samples, tau,
            gene_name, len(ref), len(alt), substrate, pos + 1), None


def _load_mapping_pairs(path, log, gene_name):
    """Read a mapping CSV into {mutant: coordinate_token}.

    Column-addressed, not positional. This used to require EXACTLY two
    comma-fields per line, which silently rejected every row of any mapping
    carrying a third column: adding the `pkey` column took RNAfold from 30/38
    scored to 20/38, reporting each row as "Skipping malformed mutation entry".

    trim_muts is still the row source because it applies the validation-log
    filtering that a bare csv reader would skip, but the header is consulted to
    locate the mutant and value columns rather than assuming positions. `pkey`
    and `orf` are identity/derived columns, never the coordinate this returns.
    """
    if path is None:
        return {}
    try:
        with open(str(path)) as fh:
            header = [h.strip().lower() for h in fh.readline().strip().split(",")]
        rows = trim_muts(str(path), log, gene_name)
    except Exception:
        return {}
    if not header:
        return {}
    mut_i = header.index("mutant") if "mutant" in header else 0
    val_i = next((i for i, h in enumerate(header)
                  if h not in ("pkey", "mutant", "orf")), len(header) - 1)
    pairs = {}
    for item in rows:
        toks = [t.strip() for t in str(item).split(",")]
        if len(toks) <= max(mut_i, val_i):
            continue
        if toks[mut_i] and toks[val_i]:
            pairs[toks[mut_i]] = toks[val_i]
    return pairs


def _load_intron_premrna(path, gene_name):
    """Read intron_premRNA_mapping_<GENE>.csv into per-substrate work lists.

    Returns (intron_jobs, premrna_jobs), each a list of (mutant, coord_token).

    Columns are mutant, orf, intron, pre-mRNA-Transcript. `orf` and `intron` are
    mutually exclusive per row, so one variant spanning a splice site occupies
    TWO rows under one mutant -- which is why these are LISTS, not dicts: a dict
    keyed on mutant would silently keep only the last row of such a variant.

    Only intronic rows are taken for the intron substrate. Exonic rows are
    already covered on the transcript substrate, and scanning them again here
    would double every exonic result.

    "Covered on the transcript substrate" means HAS AN ORF ADDRESS, not "is
    exonic". The transcript mapping is built from ORF-space tokens only, so a
    gd./ch. variant landing in a UTR or a non-coding exon has no transcript row
    and no amino-acid row -- pre-mRNA is its ONLY nucleotide coordinate. Gating
    the pre-mRNA job on a non-empty `intron` cell therefore skipped exactly those
    variants, and because they still appear in the chromosome and gDNA CSVs they
    did not look dropped. Measured on BRCA1: 4 UTR rows present in the mapping,
    0 reaching premrna_jobs, 0 present in the transcript mapping -- folded on no
    substrate at all. The anti-double-counting intent is preserved: an exonic row
    that DOES carry an orf token is still excluded here.
    """
    if path is None:
        return [], []
    intron_jobs, premrna_jobs = [], []
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                mut = (row.get("mutant") or "").strip()
                orf = (row.get("orf") or "").strip()
                intron = (row.get("intron") or "").strip()
                pre = (row.get("pre-mRNA-Transcript") or "").strip()
                if not mut:
                    continue
                # A bracketed cell is a multi-piece variant; each piece is its own
                # fold against its own record. split_piece_cell is the shared
                # reader in variant_mapping, beside parse_piece_token, which
                # already parses each piece below -- the split half used to be
                # reimplemented here and in miranda.
                for piece in split_piece_cell(intron):
                    intron_jobs.append((mut, piece))
                if intron or not orf:
                    for piece in split_piece_cell(pre):
                        premrna_jobs.append((mut, piece))
    except Exception:
        return [], []
    return intron_jobs, premrna_jobs


def _resolve_map_for_gene(map_root, gene_name, valid_exts, name_contains=None):
    """Pick the mapping file for a gene from a file or a directory tree.

    RECURSIVE. variant_mapping groups its output by gene
    (<out>/<GENE>/mappings/<type>/<file>.csv), matching every other pipeline, so
    a flat iterdir of the root finds nothing -- the files are three levels down.
    rglob handles that and the older flat <out>/mappings/<type>/ layout alike.

    name_contains filters by mapping TYPE, which recursion makes necessary: a
    tree holds transcript, chromosome, gDNA, aa, intron_premRNA and pkey files
    that all carry the same gene, so without it whichever sorts first wins.
    """
    if map_root is None:
        return None
    p = Path(map_root)
    if p.is_file():
        return p
    if not p.is_dir():
        return None
    cands = [m for m in sorted(p.rglob("*"))
             if m.is_file() and m.suffix in valid_exts
             and (not name_contains or name_contains.lower() in m.stem.lower())
             and extract_gene_from_filename(str(m)) == gene_name]
    return cands[0] if cands else None


def _autodetect_workers(n_tasks: int, cap: int = 8) -> int:
    n_cpu = os.cpu_count() or multiprocessing.cpu_count() or 1
    return max(1, min(n_cpu // 2 if n_cpu > 1 else 1, cap, n_tasks))

def main():
    parser = argparse.ArgumentParser(
        description="Compute ddG, JSD, and per-position deltau for variant-centered RNA folding windows using ViennaRNA."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="variant_mapping OUTPUT ROOT (<root>/<GENE>/fastas/), or a single "
                             "FASTA. Given the root, --transcript-mapping and "
                             "--intron-premrna-mapping both default to it.")
    parser.add_argument("-o", "--output", required=True, help="Output base directory")
    parser.add_argument("-w", "--window", type=int, default=151, help="Window size (odd; truncates near ends)")
    parser.add_argument("-tm", "--transcript-mapping",
                        help="FILE MODE ONLY. In directory mode this defaults to --input, since "
                             "variant_mapping writes mappings/transcript/ beside fastas/ under "
                             "the same <GENE>/. Supply it only for a mapping outside that layout.")
    parser.add_argument("-ipm", "--intron-premrna-mapping",
                        help="FILE MODE ONLY. In directory mode this is derived from the same "
                             "root, since mappings/intron_premRNA/ is a sibling of "
                             "mappings/transcript/. Intronic variants are folded against BOTH "
                             "the intron record and the pre_mRNA record.")
    parser.add_argument("-l", "--log", help="Validation log (file or dir) used to filter failed mutations")
    parser.add_argument("-s", "--samples", type=int, default=1000, help="Number of Boltzmann samples per sequence")
    parser.add_argument("-ta", "--tau", type=float, default=0.05, help="Threshold for change_flag on deltau")
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

    if input_path.is_file():
        files = [input_path]
    else:
        # discover_fasta_files first: it understands variant_mapping's per-gene
        # layout (<root>/<GENE>/fastas/<GENE>.fasta), which a flat iterdir of the
        # root cannot see -- that root holds only directories. Measured before
        # this: `-i out/` raised "No valid input files found in out" while
        # out/NPM1/fastas/ and out/PAM/fastas/ each held one.
        _disc = discover_fasta_files(str(input_path))
        files = ([Path(p) for _, p in sorted(_disc.items())] if _disc
                 else [i for i in input_path.iterdir() if i.suffix in valid_exts_fa])
    if not files:
        raise ValueError(f'No valid input files found in {input_path.name}')
    # Same defaulting for the intron/pre-mRNA mapping: variant_mapping writes it as
    # a sibling of mappings/transcript/ under the same <GENE>/, so a gene-layout
    # root supplies both.
    if not args.intron_premrna_mapping:
        _ipm_root = args.transcript_mapping or args.input
        if _ipm_root and discover_mapping_files(str(_ipm_root), "intron_premrna"):
            args.intron_premrna_mapping = _ipm_root
            print(f"[rnafold] --intron-premrna-mapping not given; using {_ipm_root}")

    if transcpt is None:
        # Default it to the input root when that root carries the per-gene layout:
        # variant_mapping writes <root>/<GENE>/mappings/transcript/ beside
        # <root>/<GENE>/fastas/, so the two are the same directory and requiring
        # the flag makes the caller repeat themselves. discover_mapping_files
        # confirms a transcript mapping is actually there before defaulting;
        # otherwise the original error stands.
        if (not input_path.is_file()
                and discover_mapping_files(str(input_path), "transcript")):
            transcpt = input_path
            print(f"[rnafold] --transcript-mapping not given; using {input_path} "
                  f"(per-gene layout detected)")
        else:
            raise ValueError("--transcript-mapping is required")

    # Recursive + type-filtered, for the same reason as _resolve_map_for_gene:
    # the mappings sit at <root>/<GENE>/mappings/transcript/ under the gene-first
    # layout, and a flat scan of the root sees none of them.
    maps = ([transcpt] if transcpt.is_file()
            else [t for t in sorted(transcpt.rglob("*"))
                  if t.is_file() and t.suffix in valid_exts_map
                  and "transcript" in t.stem.lower()])
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
        # Column-addressed, like _load_mapping_pairs. A hardcoded "exactly two
        # comma-fields" rejected every row of a mapping with a third column;
        # adding `pkey` took this loop from 30/38 scored to 20/38, each row
        # reported as "malformed_token". The header is read once per gene.
        try:
            with open(str(transcript_map)) as _hf:
                _hdr = [h.strip().lower() for h in _hf.readline().strip().split(",")]
        except Exception:
            _hdr = []
        _mut_i = _hdr.index("mutant") if "mutant" in _hdr else 0
        _val_i = next((i for i, h in enumerate(_hdr)
                       if h not in ("pkey", "mutant", "orf")), 1 if len(_hdr) < 2 else len(_hdr) - 1)
        _pkey_i = _hdr.index("pkey") if "pkey" in _hdr else None

        for item in mutants:
            run_summary["mutations_total"] += 1
            toks = [t.strip() for t in str(item).split(",")]
            if len(toks) <= max(_mut_i, _val_i) or not toks[_mut_i] or not toks[_val_i]:
                print(f"Skipping malformed mutation entry: {item}", file=sys.stderr)
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": None, "mapped_mut": str(item), "reason": "malformed_token"})
                continue
            reference_mut, mapped_mut = toks[_mut_i], toks[_val_i]
            # Prefer the minted key from the mapping over re-deriving it. The
            # local mint is kept as the fallback for a mapping written before
            # the pkey column existed.
            pkey = (toks[_pkey_i] if _pkey_i is not None and len(toks) > _pkey_i and toks[_pkey_i]
                    else mint_pkey(gene_name, reference_mut))

            # F24: an aa/consequence token ('Stop'/'Sto') has no nt position to fold;
            # bioAccurate returns None for it. A nt nonsense variant (e.g. C175T) never
            # matches and folds normally. Unparseable tokens are logged, not crashed on.
            # get_variant_data_bioAccurate never raises; it returns (None, None)
            # for a stop/aa/off-alphabet token. The legacy bioAccurate call also
            # raised ValueError on any multi-base token, which meant an indel was
            # rejected here before any of the checks below could see it.
            transcript_pos, _ = get_variant_data_bioAccurate(mapped_mut, is_nt=True)
            if transcript_pos is None:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": "aa_level_token"})
                continue

            # Non-SNV tokens are processed by default. There is no flag: the
            # grammar is uniquely decodable so no parser mode has to be selected,
            # and whether a column survives is decided per column from
            # len(ref) != len(alt) -- a fact of the record, not a user preference.
            #
            # All bounds/REF/window logic now lives in _prepare_fold, shared with
            # the intron and pre-mRNA substrates. No separate "deletion consumed
            # the window" guard: it is unreachable (verified across deletion
            # lengths 2..399: 248 caught by the offset check, 0 here).
            item, reason = _prepare_fold(
                transcript_seq, mapped_mut, args.window, pkey, transcript_pos,
                gene_name, args.samples, args.tau, "transcript")
            if item is None:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": reference_mut, "mapped_mut": mapped_mut, "reason": reason})
                continue
            work_items.append(item)

        # ---- intronic variants -------------------------------------------
        # These never appear in the transcript mapping and never can: an intron
        # has no transcript coordinate. They arrive through their own mappings
        # and are folded against BOTH the excised intron and the pre-mRNA.
        ipm_f = _resolve_map_for_gene(args.intron_premrna_mapping, gene_name, valid_exts_map,
                                      name_contains="intron_premRNA")
        intron_jobs, premrna_jobs = _load_intron_premrna(ipm_f, gene_name)

        def _attempt(mut, coord_tok, substrate, seq, seq_key):
            """One fold attempt. Counted individually.

            An intronic variant is TWO attempts, not one, and a variant spanning
            a splice site can be more. mutations_total must count ATTEMPTS:
            counting the token once while allowing several failure entries made
            mutations_successful (= total - unsuccessful) go NEGATIVE -- observed
            as "-2/4 ok" on the first run.
            """
            run_summary["mutations_total"] += 1
            # Bounded key, identical to the one variant_mapping wrote into the
            # mapping's pkey column: both hash the same verbatim token. Minting
            # '{gene}-{mut}' here made the key as long as the variant -- 2,215
            # chars for the 2,209 nt deletion in the test fixture.
            pkey = mint_pkey(gene_name, mut)
            if not coord_tok:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": mut, "mapped_mut": None, "reason": f"{substrate}:no_coordinate"})
                return
            if seq is None:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": mut, "mapped_mut": coord_tok, "reason": f"missing_fasta_record:{seq_key}"})
                return
            item, reason = _prepare_fold(
                seq, coord_tok, args.window, pkey, "", gene_name,
                args.samples, args.tau, substrate)
            if item is None:
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": mut, "mapped_mut": coord_tok, "reason": f"{substrate}:{reason}"})
            else:
                work_items.append(item)

        # Substrate 1: the excised intron. "intron<N>:REF<pos>ALT" carries the
        # record name alongside the coordinate, so a row can never be attributed
        # to the wrong intron.
        for mut, coord in intron_jobs:
            if ":" not in coord:
                run_summary["mutations_total"] += 1
                run_summary["unsuccessful"].append({"gene": gene_name, "reference_mut": mut, "mapped_mut": coord, "reason": "malformed_intron_coord"})
                continue
            label, tok = coord.split(":", 1)
            record = next((k for k in fasta_dict if k.split("|", 1)[0] == label), None)
            _attempt(mut, tok, label, fasta_dict.get(record) if record else None, label)

        # Substrate 2: the same pieces in continuous pre-mRNA context.
        for mut, coord in premrna_jobs:
            _attempt(mut, coord, "pre_mRNA", fasta_dict.get("pre_mRNA"), "pre_mRNA")

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
    # Written PER GENE, under <out>/<GENE>/RNAfold/, so the run leaves nothing
    # loose at the output root -- every other artifact this pipeline emits is
    # already grouped that way. The per-gene file carries only that gene's
    # entries, with the counters recomputed for it, so it is readable on its own
    # rather than being a copy of a run-wide total.
    Path(args.output).mkdir(parents=True, exist_ok=True)

    def _gene_of(entry):
        # Substrate keys are '<GENE>~<SUBSTRATE>'; the summary is per gene.
        g = str(entry.get("gene", "") or "")
        return g.split("~", 1)[0]

    genes_seen = {g for g in (
        [_gene_of(e) for e in run_summary.get("unsuccessful", [])]
        + [_gene_of(e) for e in run_summary.get("skipped_genes", [])]
        + list(results_by_gene.keys())
    ) if g}
    # One timestamp for the whole run, so every gene's summary from a single
    # invocation shares a suffix and sorts together. Matches variant_mapping's
    # validation_<YYYYmmdd_HHMMSS>.log convention.
    from datetime import datetime
    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for gname in sorted(genes_seen):
        unsucc = [e for e in run_summary.get("unsuccessful", []) if _gene_of(e) == gname]
        skipped = [e for e in run_summary.get("skipped_genes", []) if _gene_of(e) == gname]
        per = dict(run_summary)
        per["gene"] = gname
        per["unsuccessful"] = unsucc
        per["skipped_genes"] = skipped
        per["mutations_unsuccessful"] = len(unsucc)
        gene_dir = Path(args.output) / gname / "RNAfold"
        gene_dir.mkdir(parents=True, exist_ok=True)
        summary_path = gene_dir / f"rnafold.run_summary.{run_stamp}.json"
        with open(summary_path, "w") as f:
            json.dump(per, f, indent=2)
        # Name the actual file, not the "<gene>" template. The summary is written
        # PER GENE, so a run over several genes produces several of these and the
        # placeholder told the reader neither which genes nor where to look.
        print(f"[SUMMARY] {gname}: mutations "
              f"{per['mutations_total'] - per['mutations_unsuccessful']}/{per['mutations_total']} ok, "
              f"{per['mutations_unsuccessful']} unsuccessful -> {summary_path}")
    print(f"[SUMMARY] genes {run_summary['genes_processed']}/{run_summary['genes_total']}, "
          f"mutations {run_summary['mutations_successful']}/{run_summary['mutations_total']} ok, "
          f"{run_summary['mutations_unsuccessful']} unsuccessful")

if __name__ == "__main__":
    main()
