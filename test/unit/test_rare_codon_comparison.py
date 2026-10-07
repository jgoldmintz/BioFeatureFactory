"""Allele comparison against real cg_cotrans WT windows and fixed references."""

import csv
import gzip
import importlib
import math
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys

import pytest


pytest.importorskip("calc_rare_enrichment", reason="Real cg_cotrans is required; supply its directory on PYTHONPATH")
rc = importlib.import_module("biofeaturefactory.rare_codon.rare_codon_pipeline")

REPO = Path(__file__).resolve().parents[2]
LEGACY_FIELDS = (
    "pkey", "Gene", "codon_position", "p_enriched", "p_depleted", "f_enriched_wt",
    "frac_seq_enriched", "frac_seq_depleted", "n_rare", "window_size", "qc_flags",
)
DELTA_FIELDS = ("n_rare_mut", "f_enriched_mut", "delta_n_rare", "delta_f_enriched")


def write_text(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


def reference_usage():
    families = {}
    for codon, amino_acid in rc.codon_to_aa.items():
        if amino_acid != "Stop":
            families.setdefault(amino_acid, []).append(codon)
    usage = {codon: 1 / len(family) for family in families.values() for codon in family}
    usage.update({"GCT": 0.70, "GCC": 0.04, "GCA": 0.20, "GCG": 0.06})
    usage.update({"GTT": 0.05, "GTC": 0.25, "GTA": 0.30, "GTG": 0.40})
    return usage


def write_usage(path, usage=None, null_usage=None, genes=("GENE",)):
    usage = reference_usage() if usage is None else usage
    null_usage = usage if null_usage is None else null_usage
    identifiers = ("ORF", "homolog")
    payload = {
        "groups": ["all"],
        "gene_groups": {gene: 0 for gene in genes},
        "overall_codon_usage": {identifier: usage for identifier in identifiers},
        "unweighted_codon_usage": {identifier: null_usage for identifier in identifiers},
        "gene_group_codon_usage": {identifier: {0: null_usage} for identifier in identifiers},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        pickle.dump(payload, handle)
    return path


def write_reference(path, usage):
    return write_text(path, "codon\trelative_usage_within_aa\n" + "".join(
        f"{codon}\t{frequency}\n" for codon, frequency in sorted(usage.items())
    ))


def focus_codons(length=42):
    return ["GCC" if index % 7 == 0 else "GCT" for index in range(length)]


def analyze(tmp_path, codons=None, width=5, usage=None, null_usage=None, homolog_sequence=None, **options):
    codons = focus_codons() if codons is None else codons
    sequence = "".join(codons)
    homolog_sequence = sequence if homolog_sequence is None else homolog_sequence
    msa = write_text(tmp_path / "GENE.codon.msa.fasta", f">ORF\n{sequence}\n>homolog\n{homolog_sequence}\n")
    usage_path = write_usage(tmp_path / "usage.p.gz", usage, null_usage)
    context = {}
    results, count = rc.run_rare_codon_analysis(
        "GENE", str(msa), str(usage_path), "ORF", window_size=width,
        comparison_context=context, **options,
    )
    assert context
    return sequence.replace("-", "").upper(), results, context, count


def compare(tokens, sequence, results, context, width=5):
    return rc.process_mutations(tokens, "GENE", sequence, results, width, comparison_context=context)


def assert_legacy_equal(actual, expected):
    for field in LEGACY_FIELDS:
        if isinstance(expected[field], float) and math.isnan(expected[field]):
            assert math.isnan(actual[field]), field
        else:
            assert actual[field] == expected[field], field


def assert_delta(row, wt_count, mut_count, width):
    assert row["comparison_status"] == "PASS", row
    assert row["n_rare"] == wt_count
    assert row["n_rare_mut"] == mut_count
    assert row["f_enriched_wt"] == pytest.approx(wt_count / width)
    assert row["f_enriched_mut"] == pytest.approx(mut_count / width)
    assert row["delta_n_rare"] == mut_count - wt_count
    assert row["delta_f_enriched"] == pytest.approx((mut_count - wt_count) / width)


def assert_refused(row):
    assert row["comparison_status"] not in ("", "PASS", "NOT_REQUESTED"), row
    assert all(row[field] == "" for field in DELTA_FIELDS), row


@pytest.mark.parametrize(
    "source,token,mutant,expected_delta",
    [
        ("GCT", "T99C", "GCC", 1),
        ("GCC", "C99T", "GCT", -1),
        ("GCT", "T99A", "GCA", 0),
        ("GCC", "C99G", "GCG", 0),
        ("GCT", "C98T", "GTT", 1),
    ],
    ids=["synonymous-gain", "synonymous-loss", "synonymous-common", "synonymous-rare", "missense-gain"],
)
def test_exact_allele_counts_preserve_wt_statistics(tmp_path, source, token, mutant, expected_delta):
    codons = focus_codons()
    codons[32] = source
    sequence, results, context, count = analyze(tmp_path, codons)
    original_results = {position: values.copy() for position, values in results.items()}
    expected_wt = sum(codon in {"GCC", "GCG", "GTT"} for codon in codons[30:35])
    mutated = codons.copy()
    mutated[32] = mutant
    expected_mut = sum(codon in {"GCC", "GCG", "GTT"} for codon in mutated[30:35])
    assert expected_mut - expected_wt == expected_delta
    row = compare([token], sequence, results, context)[0]
    assert_delta(row, expected_wt, expected_mut, 5)
    legacy = rc.process_mutations([token], "GENE", sequence, results, 5)[0]
    assert_legacy_equal(row, legacy)
    assert legacy["comparison_status"] == "NOT_REQUESTED"
    assert all(legacy[field] == "" for field in DELTA_FIELDS)
    assert results == original_results
    assert count == 2


def test_equal_length_mnv_across_codons_uses_whole_allele(tmp_path):
    codons = focus_codons()
    codons[32:34] = ["GCT", "GCT"]
    sequence, results, context, _ = analyze(tmp_path, codons)
    row = compare(["TGCT99CGCC"], sequence, results, context)[0]
    expected_wt = codons[31:36].count("GCC")
    assert_delta(row, expected_wt, expected_wt + 2, 5)
    assert row["codon_position"] == 34


def test_same_position_variants_are_independent_and_do_not_recompute_wt(tmp_path, monkeypatch):
    invocations = []
    original = rc.msa_rare_codon_analysis_wtalign_nseq

    def count_analysis(*arguments, **keywords):
        invocations.append(1)
        return original(*arguments, **keywords)

    monkeypatch.setattr(rc, "msa_rare_codon_analysis_wtalign_nseq", count_analysis)
    codons = focus_codons()
    codons[32] = "GCT"
    sequence, results, context, _ = analyze(tmp_path, codons)
    tokens = ["T99C", "T99A", "T99C"]
    rows = compare(tokens, sequence, results, context)
    assert invocations == [1]
    assert [row["delta_n_rare"] for row in rows] == [1, 0, 1]
    assert [row["pkey"] for row in rows] == [rc.mint_pkey("GENE", token) for token in tokens]
    assert rows[0] == rows[2]
    assert [row["n_rare"] for row in rows] == [rows[0]["n_rare"]] * 3


@pytest.mark.parametrize("width", [1, 3, 4, 15])
def test_exact_centered_full_window_and_edge_refusals(tmp_path, width):
    codons = focus_codons()
    sequence, results, context, _ = analyze(tmp_path, codons, width)
    first_center = width // 2
    last_center = len(codons) - width + width // 2
    assert sorted(results) == list(range(first_center + 1, last_center + 2))
    centers = sorted({0, first_center, 20, last_center, len(codons) - 1})
    for center in centers:
        token = f"{codons[center][2]}{center * 3 + 3}A"
        row = compare([token], sequence, results, context, width)[0]
        if center < first_center or center > last_center:
            assert_refused(row)
            assert "POSITION_NOT_IN_WINDOW" in row["qc_flags"]
        else:
            start = center - width // 2
            window = codons[start:start + width]
            wt_count = window.count("GCC")
            mut_count = wt_count - int(codons[center] == "GCC")
            assert_delta(row, wt_count, mut_count, width)


@pytest.mark.parametrize("length,width", [(5, 5), (5, 6), (1, 1)])
def test_small_msa_obeys_available_full_windows(tmp_path, length, width):
    codons = ["GCT"] * length
    sequence, results, context, _ = analyze(tmp_path, codons, width)
    center = min(width // 2, length - 1)
    row = compare([f"T{center * 3 + 3}C"], sequence, results, context, width)[0]
    if width > length:
        assert_refused(row)
        assert results == {}
    else:
        assert_delta(row, 0, 1, width)


@pytest.mark.parametrize("token", ["A99C", "G999A", "G0A", "bad-token", "GCT97G", "G97GCTG", "G97GT"])
def test_invalid_or_unsupported_alleles_never_emit_deltas(tmp_path, token):
    codons = focus_codons()
    codons[32] = "GCT"
    sequence, results, context, _ = analyze(tmp_path, codons)
    row = compare([token], sequence, results, context)[0]
    assert_refused(row)
    legacy = rc.process_mutations([token], "GENE", sequence, results, 5)[0]
    assert_legacy_equal(row, legacy)


@pytest.mark.parametrize("source,token", [("CAG", "C97T"), ("TAA", "T97C")])
def test_stop_gain_and_loss_are_not_sense_comparisons(tmp_path, source, token):
    codons = focus_codons()
    codons[32] = source
    sequence, results, context, _ = analyze(tmp_path, codons)
    assert_refused(compare([token], sequence, results, context)[0])


def test_terminal_stop_does_not_block_upstream_comparison(tmp_path):
    codons = focus_codons()
    codons[32] = "GCT"
    codons.append("TAA")
    sequence, results, context, _ = analyze(tmp_path, codons)
    row = compare(["T99C"], sequence, results, context)[0]
    assert_delta(row, codons[30:35].count("GCC"), codons[30:35].count("GCC") + 1, 5)


def test_whole_codon_alignment_gaps_do_not_shift_mutation_coordinates(tmp_path):
    codons = focus_codons()
    codons[32] = "GCT"
    sequence, results, context, _ = analyze(tmp_path / "plain", codons)
    plain = compare(["T99C"], sequence, results, context)[0]
    codons[10:10] = ["---", "---"]
    codons[:0] = ["---"]
    gapped_sequence, gapped_results, gapped_context, _ = analyze(tmp_path / "gapped", codons)
    gapped = compare(["T99C"], gapped_sequence, gapped_results, gapped_context)[0]
    assert gapped_sequence == sequence
    assert gapped == plain


def test_homolog_insertion_at_focus_gap_does_not_change_focus_delta(tmp_path):
    codons = focus_codons()
    codons[32] = "GCT"
    codons.insert(31, "---")
    homolog_codons = codons.copy()
    homolog_codons[31] = "GCC"
    sequence, results, context, _ = analyze(tmp_path, codons, homolog_sequence="".join(homolog_codons))
    row = compare(["T99C"], sequence, results, context)[0]
    assert_delta(row, 0, 1, 5)
    assert row["codon_position"] == 33


def test_lowercase_focus_uses_unchanged_raw_nt_coordinates(tmp_path):
    codons = focus_codons()
    codons[32] = "GCT"
    sequence, results, context, _ = analyze(tmp_path / "upper", codons)
    expected = compare(["T99C"], sequence, results, context)[0]
    sequence, results, context, _ = analyze(tmp_path / "lower", [codon.lower() for codon in codons])
    assert compare(["T99C"], sequence, results, context)[0] == expected


@pytest.mark.parametrize("hole", ["NNN", "G--"])
def test_sanitized_hole_invalidates_overlapping_and_later_frames_only(tmp_path, hole):
    codons = focus_codons()
    codons[20] = hole
    sequence, results, context, _ = analyze(tmp_path, codons)
    upstream = compare(["T30C"], sequence, results, context)[0]
    assert upstream["comparison_status"] == "PASS"
    for center in (19, 32):
        position = center * 3 + 3
        if hole == "G--" and center > 20:
            position -= 2
        reference_base = sequence[position - 1]
        alternate = "A" if reference_base != "A" else "C"
        row = compare([f"{reference_base}{position}{alternate}"], sequence, results, context)[0]
        assert_refused(row)


def test_wt_disagreement_with_legacy_counts_refuses_comparison(tmp_path):
    sequence, results, context, _ = analyze(tmp_path)
    corrupted = {position: values.copy() for position, values in results.items()}
    corrupted[33]["n_rare"] += 1
    row = compare(["T99C"], sequence, corrupted, context)[0]
    assert_refused(row)
    assert row["n_rare"] == corrupted[33]["n_rare"]


def test_wt_fraction_disagreement_refuses_comparison(tmp_path):
    sequence, results, context, _ = analyze(tmp_path)
    corrupted = {position: values.copy() for position, values in results.items()}
    corrupted[33]["f_enriched_wt"] += 0.01
    assert_refused(compare(["T99C"], sequence, corrupted, context)[0])


@pytest.mark.parametrize("mode", ["empty-context", "different-focus", "different-width", "different-gene"])
def test_mismatched_or_missing_comparison_context_is_not_reused(tmp_path, mode):
    sequence, results, context, _ = analyze(tmp_path)
    gene, width = "GENE", 5
    if mode == "empty-context":
        context = {}
    elif mode == "different-focus":
        sequence = "AAA" + sequence[3:]
    elif mode == "different-width":
        width = 3
    else:
        gene = "OTHER"
    row = rc.process_mutations(["T99C"], gene, sequence, results, width, comparison_context=context)[0]
    assert_refused(row)


def test_wide_mnv_is_not_misrepresented_by_partial_window_delta(tmp_path):
    sequence, results, context, _ = analyze(tmp_path)
    reference = sequence[90:111]
    alternate = "GCC" * 7
    token = f"{reference}91{alternate}"
    row = compare([token], sequence, results, context)[0]
    assert row["comparison_status"] == "VARIANT_SPANS_WINDOW"
    assert_refused(row)


def test_analysis_with_comparison_context_preserves_original_windows(tmp_path):
    _, results, _, count = analyze(tmp_path)
    legacy_results, legacy_count = rc.run_rare_codon_analysis(
        "GENE", str(tmp_path / "GENE.codon.msa.fasta"), str(tmp_path / "usage.p.gz"),
        "ORF", window_size=5,
    )
    assert results == legacy_results
    assert count == legacy_count


def test_groups_null_model_remains_fixed_during_mutant_comparison(tmp_path):
    codons = focus_codons()
    codons[32] = "GCT"
    sequence, results, context, _ = analyze(tmp_path, codons, null_model="groups")
    rows = compare(["T99C", "T99A"], sequence, results, context)
    wt_count = codons[30:35].count("GCC")
    assert_delta(rows[0], wt_count, wt_count + 1, 5)
    assert_delta(rows[1], wt_count, wt_count, 5)


@pytest.mark.parametrize("rare_model,threshold,expected_delta", [("no_norm", 0.03, 0), ("no_norm", 0.10, 1), ("cmax_norm", 0.05, 0), ("cmax_norm", 0.10, 1)])
def test_custom_rarity_rules_are_reused_unchanged(tmp_path, rare_model, threshold, expected_delta):
    sequence, results, context, _ = analyze(tmp_path, rare_model=rare_model, rare_threshold=threshold)
    row = compare(["T99C"], sequence, results, context)[0]
    assert row["comparison_status"] == "PASS"
    assert row["delta_n_rare"] == expected_delta


@pytest.mark.parametrize("null_model,expected_delta", [("genome", 0), ("eq", 1)])
def test_zero_null_probability_retains_existing_eligibility_rule(tmp_path, null_model, expected_delta):
    null_usage = reference_usage()
    null_usage.update({"GCC": 0.0, "GCG": 0.0, "GCT": 0.8, "GCA": 0.2})
    sequence, results, context, _ = analyze(tmp_path, null_usage=null_usage, null_model=null_model)
    row = compare(["T99C"], sequence, results, context)[0]
    assert row["comparison_status"] == "PASS"
    assert row["delta_n_rare"] == expected_delta


def test_reference_override_changes_definition_without_mutating_source_usage(tmp_path):
    usage = reference_usage()
    overridden = dict(usage, GCT=0.04, GCC=0.70, GCA=0.20, GCG=0.06)
    reference = write_reference(tmp_path / "reference.tsv", overridden)
    sequence, results, context, _ = analyze(tmp_path, reference_usage=str(reference))
    row = compare(["T99C"], sequence, results, context)[0]
    assert row["comparison_status"] == "PASS"
    assert row["delta_n_rare"] == -1
    with gzip.open(tmp_path / "usage.p.gz", "rb") as handle:
        assert pickle.load(handle)["overall_codon_usage"]["ORF"] == usage


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize("reference_flag", [False, True], ids=["usage-cache", "reference-override"])
def test_real_cli_detached_file_and_gene_root_comparisons_match(tmp_path, reference_flag):
    genes = ("NPM1", "PAM")
    root = tmp_path / "input root"
    detached = tmp_path / "detached files"
    detached.mkdir()
    tokens = ("T99C", "T99A", "C98T", "G97GT")
    codons = focus_codons()
    codons[32] = "GCT"
    sequence = "".join(codons)
    for gene in genes:
        msa = write_text(root / gene / "CodonMSA" / f"{gene}.codon.msa.fasta", f">ORF\n{sequence}\n>homolog\n{sequence}\n")
        mutations = write_text(root / gene / "mappings" / "mutations" / f"{gene}_mutations.csv", "mutant\n" + "\n".join(tokens) + "\n")
        write_text(root / gene / "mappings" / "aa" / f"{gene}_aa_mapping.csv", "pkey,mutant\nWRONG,C1T\n")
        shutil.copy2(msa, detached / msa.name)
        shutil.copy2(mutations, detached / mutations.name)
    usage = write_usage(tmp_path / "usage.p.gz", genes=genes)
    common = ["-u", str(usage), "-L", "5"]
    if reference_flag:
        reference = write_reference(tmp_path / "reference.tsv", dict(reference_usage(), GCT=0.04, GCC=0.70, GCA=0.20, GCG=0.06))
        common.extend(["-rcu", str(reference)])
    environment = dict(os.environ)
    environment.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": os.pathsep.join([str(REPO), str(Path(importlib.import_module("calc_rare_enrichment").__file__).parent), os.environ.get("PYTHONPATH", "")]),
        "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    })
    script = REPO / "biofeaturefactory" / "rare_codon" / "rare_codon_pipeline.py"
    runs = [(root, None, tmp_path / "directory output")]
    runs.extend((detached / f"{gene}.codon.msa.fasta", detached / f"{gene}_mutations.csv", tmp_path / "file output") for gene in genes)
    for source, mutations, output in runs:
        arguments = [sys.executable, str(script), "-a", str(source), "-o", str(output), *common]
        if mutations is not None:
            arguments.extend(["-m", str(mutations)])
        completed = subprocess.run(arguments, cwd=tmp_path, env=environment, text=True, capture_output=True, timeout=60)
        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert "Error" not in completed.stderr, completed.stdout + completed.stderr
    for gene in genes:
        suffix = Path(gene) / "RareCodon" / f"{gene}.rare_codon.tsv"
        directory_rows = read_rows(tmp_path / "directory output" / suffix)
        file_rows = read_rows(tmp_path / "file output" / suffix)
        assert directory_rows == file_rows
        assert len(file_rows) == 4
        assert [row["pkey"] for row in file_rows] == [rc.mint_pkey(gene, token) for token in tokens]
        expected_deltas = [-1, -1, 0] if reference_flag else [1, 0, 1]
        assert [float(row["delta_n_rare"]) for row in file_rows[:3]] == expected_deltas
        assert all(row["comparison_status"] == "PASS" for row in file_rows[:3])
        assert_refused(file_rows[-1])
        assert all(row["p_enriched"] and row["p_depleted"] for row in file_rows)
        assert set(LEGACY_FIELDS).issubset(file_rows[0])
