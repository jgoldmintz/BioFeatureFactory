import json
import sys

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller


def write_file(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)
    return path


@pytest.fixture
def gene_input(tmp_path):
    root = tmp_path / "input"
    files = {
        "fasta": write_file(root / "GENE" / "fastas" / "GENE.fasta", ">ORF\nATGGCTTAA\n"),
        "mutations": write_file(root / "GENE" / "mappings" / "mutations" / "GENE_mutations.csv", "mutant\nG4A\n"),
        "msa": write_file(root / "GENE" / "MSA" / "GENE.msa.a2m", ">ORF\nMA\n"),
        "codon_msa": write_file(root / "GENE" / "CodonMSA" / "GENE.codon.msa.fasta", ">ORF\nATGGCTTAA\n"),
    }
    hardware = write_file(tmp_path / "hardware.json", json.dumps({"cpus": 8, "memory_gib": 32, "gpus": []}))
    return root, files, hardware, tmp_path / "results"


def parse_cli(monkeypatch, gene_input, *extra):
    root, files, hardware, output = gene_input
    monkeypatch.setattr(sys, "argv", [
        "mutEffects_controller", "--fasta", str(root), "--output", str(output),
        "--adabmdca-only", "--resource-hardware", str(hardware), *map(str, extra),
    ])
    return controller.parse_args()


@pytest.mark.parametrize("mutations,expected", [
    ("G4A\n", {"protein"}),
    ("T6C\n", {"codon"}),
    ("T7C\n", {"codon"}),
    ("T6C\nT7C\n", {"codon"}),
    ("G4A\nT6C\nT7C\n", {"protein", "codon"}),
])
def test_automatic_gene_routes_plan_only_required_sides(gene_input, monkeypatch, mutations, expected):
    root, files, hardware, output = gene_input
    files["mutations"].write_text("mutant\n" + mutations)
    args = parse_cli(monkeypatch, gene_input)
    manifest = controller.build_manifest(["GENE"], args)
    config, plans = controller.prepare_resource_plans(args, manifest, ["GENE"])
    assert {plan["side"] for plan in plans} == expected
    for side in ("protein", "codon"):
        assert manifest["routing"]["GENE"][side] == (side in expected)


@pytest.mark.parametrize("mutations,unused", [("G4A\n", "codon_msa"), ("T6C\n", "msa")])
def test_unused_alignment_needs_no_database(gene_input, monkeypatch, mutations, unused):
    root, files, hardware, output = gene_input
    files["mutations"].write_text("mutant\n" + mutations)
    files[unused].unlink()
    args = parse_cli(monkeypatch, gene_input)
    manifest = controller.build_manifest(["GENE"], args)
    controller.validate_db_coverage(["GENE"], manifest, args)


def test_oversized_unused_codon_does_not_abort_protein_run(gene_input, monkeypatch, capsys):
    root, files, hardware, output = gene_input
    overrides = write_file(root.parent / "overrides.json", json.dumps({"GENE.codon": {"cpu_memory_gib": 1000000}}))
    args = parse_cli(monkeypatch, gene_input, "--resource-plan-only", "--resource-overrides", overrides)
    capsys.readouterr()
    controller.run_controller(args)
    report = json.loads(capsys.readouterr().out)
    assert [plan["side"] for plan in report["tasks"]] == ["protein"]
    assert not output.exists()


@pytest.mark.parametrize("option,source,mutations,side,warning", [
    ("--msa", "msa", "T6C\nT7C\n", "protein", "this will not produce biologically accurate results"),
    ("-cm", "codon_msa", "G4A\n", "codon", "codon"),
])
def test_explicit_mode_overrides_and_warns(gene_input, monkeypatch, capsys, option, source, mutations, side, warning):
    root, files, hardware, output = gene_input
    files["mutations"].write_text("mutant\n" + mutations)
    args = parse_cli(monkeypatch, gene_input, option, root, "--resource-plan-only")
    capsys.readouterr()
    controller.run_controller(args)
    printed = capsys.readouterr()
    report = json.loads(printed.out)
    assert warning in printed.err
    assert [plan["side"] for plan in report["tasks"]] == [side]
    settings = report["tasks"][0]["settings"]
    assert settings["skip_codon"] if side == "protein" else settings["score_missense_codon"]


def test_both_explicit_sources_enable_both_sides(gene_input, monkeypatch):
    root, files, hardware, output = gene_input
    args = parse_cli(monkeypatch, gene_input, "--msa", root, "-cm", root)
    manifest = controller.build_manifest(["GENE"], args)
    config, plans = controller.prepare_resource_plans(args, manifest, ["GENE"])
    assert {plan["side"] for plan in plans} == {"protein", "codon"}
    assert not next(plan for plan in plans if plan["side"] == "codon")["settings"]["score_missense_codon"]


def test_mode_changes_invalidate_codon_fingerprint(gene_input, monkeypatch):
    root, files, hardware, output = gene_input
    args = parse_cli(monkeypatch, gene_input, "--msa", root, "-cm", root)
    manifest = controller.build_manifest(["GENE"], args)
    config, plans = controller.prepare_resource_plans(args, manifest, ["GENE"])
    combined = next(plan for plan in plans if plan["side"] == "codon")
    args = parse_cli(monkeypatch, gene_input, "-cm", root)
    manifest = controller.build_manifest(["GENE"], args)
    config, plans = controller.prepare_resource_plans(args, manifest, ["GENE"])
    assert plans[0]["fingerprint"] != combined["fingerprint"]


def test_skip_codon_keeps_existing_protein_override(gene_input, monkeypatch):
    root, files, hardware, output = gene_input
    files["mutations"].write_text("mutant\nT6C\n")
    files["msa"].unlink()
    args = parse_cli(monkeypatch, gene_input, "--skip-codon")
    manifest = controller.build_manifest(["GENE"], args)
    assert controller.side_enabled(args, manifest, "GENE", "adabmdca", "protein")
    assert not controller.side_enabled(args, manifest, "GENE", "adabmdca", "codon")
    assert manifest["routing"]["GENE"]["warnings"]
    with pytest.raises(SystemExit, match="protein MSA"):
        controller.validate_db_coverage(["GENE"], manifest, args)
