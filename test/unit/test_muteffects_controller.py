from argparse import Namespace
from pathlib import Path

import pytest

from biofeaturefactory.mutation_effects import mutEffects_controller as controller


def _write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _gene_tree(root, gene, protein_msa=True, codon_msa=True):
    paths = {
        "fasta": root / gene / "fastas" / f"{gene}.fasta",
        "mutations": root / gene / "mappings" / "mutations" / f"{gene}_mutations.csv",
        "msa": root / gene / "MSA" / f"{gene}.msa.a2m",
        "codon_msa": root / gene / "CodonMSA" / f"{gene}.codon.msa.fasta",
    }
    _write(paths["fasta"], ">ORF\nATGGCTTAA\n")
    _write(paths["mutations"], "mutant\nG4A\nT6C\n")
    if protein_msa:
        _write(paths["msa"], ">ORF\nMA\n")
    if codon_msa:
        _write(paths["codon_msa"], ">ORF\nATGGCT\n")
    return paths


def _args(fasta, output, **overrides):
    values = {
        "fasta": fasta,
        "mutations": fasta,
        "output": output,
        "msa": None,
        "codon_msa": None,
        "model_params": None,
        "codon_model_params": None,
        "adabmdca_protein_params": None,
        "adabmdca_codon_params": None,
        "db_root": None,
        "run_evmutation": False,
        "run_adabmdca": True,
        "skip_codon_evmutation": False,
        "skip_codon_adabmdca": False,
    }
    values.update(overrides)
    return Namespace(**values)


def _manifest(args):
    genes = controller._resolve_genes(args.fasta)
    return genes, controller.build_manifest(genes, args)


def test_nested_inputs_are_resolved_when_output_is_elsewhere(tmp_path):
    source = tmp_path / "source"
    expected = {gene: _gene_tree(source, gene) for gene in ("NPM1", "PAM")}
    _write(source / "NPM1" / "mappings" / "aa" / "NPM1_aa_mapping.csv",
           "mutant,aa\nG4A,A2T\n")
    _write(source / "NPM1" / "CodonMSA" / "NPM1_100.codon.msa.fasta",
           ">ORF\nATGGCT\n")
    args = _args(source, tmp_path / "results")

    genes, manifest = _manifest(args)

    assert genes == ["NPM1", "PAM"]
    for gene, paths in expected.items():
        assert manifest["input_files"][gene] == {
            kind: str(path.resolve()) for kind, path in paths.items()
        }
    assert manifest["msa"] == genes
    assert manifest["codon_msa"] == genes
    controller.validate_db_coverage(genes, manifest, args)


def test_flat_inputs_preserve_gene_matching(tmp_path):
    source = tmp_path / "fastas"
    mutations = tmp_path / "mutations"
    protein = tmp_path / "protein"
    codon = tmp_path / "codon"
    expected = {}
    for gene in ("NPM1", "PAM"):
        expected[gene] = {
            "fasta": _write(source / f"{gene}.fasta", ">ORF\nATGGCTTAA\n"),
            "mutations": _write(mutations / f"{gene}_mutations.csv", "mutant\nG4A\n"),
            "msa": _write(protein / f"{gene}.msa.a2m", ">ORF\nMA\n"),
            "codon_msa": _write(codon / f"{gene}.codon.msa.fasta", ">ORF\nATGGCT\n"),
        }
    args = _args(source, tmp_path / "results", mutations=mutations,
                 msa=protein, codon_msa=codon)

    genes, manifest = _manifest(args)

    assert set(genes) == {"NPM1", "PAM"}
    for gene, paths in expected.items():
        assert manifest["input_files"][gene] == {
            kind: str(path.resolve()) for kind, path in paths.items()
        }
    controller.validate_db_coverage(genes, manifest, args)


def test_single_fasta_and_explicit_files_are_preserved(tmp_path):
    fasta = _write(tmp_path / "NPM1.fasta", ">ORF\nATGGCTTAA\n")
    mutations = _write(tmp_path / "selected_variants.csv", "mutant\nG4A\n")
    protein = _write(tmp_path / "selected_alignment.a2m", ">ORF\nMA\n")
    codon = _write(tmp_path / "selected_codons.fasta", ">ORF\nATGGCT\n")
    args = _args(fasta, tmp_path / "results", mutations=mutations,
                 msa=protein, codon_msa=codon)

    genes, manifest = _manifest(args)

    assert genes == ["NPM1"]
    assert manifest["input_files"]["NPM1"] == {
        "fasta": str(fasta.resolve()),
        "mutations": str(mutations.resolve()),
        "msa": str(protein.resolve()),
        "codon_msa": str(codon.resolve()),
    }


@pytest.mark.parametrize("kind,suffix", [
    ("msa", "msa.a2m"),
    ("codon_msa", "codon.msa.fasta"),
])
def test_explicit_msa_directory_overrides_automatic_sources(tmp_path, kind, suffix):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    output = tmp_path / "results"
    _gene_tree(output, "NPM1")
    explicit = tmp_path / "explicit"
    chosen = _write(explicit / f"NPM1.{suffix}", ">ORF\nATGGCT\n")
    args = _args(source, output, **{kind: explicit})

    _, manifest = _manifest(args)

    assert manifest["input_files"]["NPM1"][kind] == str(chosen.resolve())


@pytest.mark.parametrize("kind", ["msa", "codon_msa"])
def test_missing_explicit_msa_does_not_fall_back(tmp_path, kind):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    output = tmp_path / "results"
    _gene_tree(output, "NPM1")
    args = _args(source, output, **{kind: tmp_path / "missing"})

    genes, manifest = _manifest(args)

    assert manifest["input_files"]["NPM1"].get(kind) is None
    assert manifest[kind] == []
    with pytest.raises(SystemExit, match="--db-root is required"):
        controller.validate_db_coverage(genes, manifest, args)


def test_partial_protein_coverage_requires_database(tmp_path):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    _gene_tree(source, "PAM", protein_msa=False)
    args = _args(source, tmp_path / "results")

    genes, manifest = _manifest(args)

    assert manifest["msa"] == ["NPM1"]
    assert manifest["codon_msa"] == ["NPM1", "PAM"]
    with pytest.raises(SystemExit, match="protein MSA for 1 gene.*PAM"):
        controller.validate_db_coverage(genes, manifest, args)

    args.db_root = tmp_path / "Bio_DBs"
    _write(args.db_root / "uniref90.fasta", ">protein\nMA\n")
    controller.validate_db_coverage(genes, manifest, args)


def test_partial_codon_coverage_requires_database(tmp_path):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    _gene_tree(source, "PAM", codon_msa=False)
    args = _args(source, tmp_path / "results")

    genes, manifest = _manifest(args)

    assert manifest["codon_msa"] == ["NPM1"]
    with pytest.raises(SystemExit, match="codon MSA for 1 gene.*PAM"):
        controller.validate_db_coverage(genes, manifest, args)


def test_flat_output_cache_has_priority_over_both_nested_roots(tmp_path):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    output = tmp_path / "results"
    _gene_tree(output, "NPM1")
    protein = _write(output / "MSA" / "NPM1.msa.a2m", ">ORF\nMA\n")
    codon = _write(output / "CodonMSA" / "NPM1.codon.msa.fasta", ">ORF\nATGGCT\n")

    _, manifest = _manifest(_args(source, output))

    assert manifest["input_files"]["NPM1"]["msa"] == str(protein.resolve())
    assert manifest["input_files"]["NPM1"]["codon_msa"] == str(codon.resolve())


def test_nested_output_cache_has_priority_over_source(tmp_path):
    source = tmp_path / "source"
    source_paths = _gene_tree(source, "NPM1")
    output = tmp_path / "results"
    cached = _gene_tree(output, "NPM1")

    _, manifest = _manifest(_args(source, output))

    assert manifest["input_files"]["NPM1"]["fasta"] == str(source_paths["fasta"].resolve())
    assert manifest["input_files"]["NPM1"]["mutations"] == str(source_paths["mutations"].resolve())
    for kind in ("msa", "codon_msa"):
        assert manifest["input_files"]["NPM1"][kind] == str(cached[kind].resolve())


def test_mixed_cache_locations_are_resolved_per_gene(tmp_path):
    source = tmp_path / "source"
    npm1 = _gene_tree(source, "NPM1", protein_msa=False)
    pam = _gene_tree(source, "PAM", codon_msa=False)
    output = tmp_path / "results"
    npm1_protein = _write(output / "MSA" / "NPM1.msa.a2m", ">ORF\nMA\n")
    pam_codon = _write(output / "PAM" / "CodonMSA" / "PAM.codon.msa.fasta",
                       ">ORF\nATGGCT\n")
    args = _args(source, output)

    genes, manifest = _manifest(args)

    expected = {"NPM1": {**npm1, "msa": npm1_protein},
                "PAM": {**pam, "codon_msa": pam_codon}}
    for gene, paths in expected.items():
        assert manifest["input_files"][gene] == {
            kind: str(path.resolve()) for kind, path in paths.items()
        }
    controller.validate_db_coverage(genes, manifest, args)


@pytest.mark.parametrize("flat_subdirs", [False, True])
def test_source_root_flat_msa_fallbacks(tmp_path, flat_subdirs):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1", protein_msa=False, codon_msa=False)
    protein_root = source / "MSA" if flat_subdirs else source
    codon_root = source / "CodonMSA" if flat_subdirs else source
    protein = _write(protein_root / "NPM1.msa.a2m", ">ORF\nMA\n")
    codon = _write(codon_root / "NPM1.codon.msa.fasta", ">ORF\nATGGCT\n")

    _, manifest = _manifest(_args(source, tmp_path / "results"))

    assert manifest["input_files"]["NPM1"]["msa"] == str(protein.resolve())
    assert manifest["input_files"]["NPM1"]["codon_msa"] == str(codon.resolve())


def test_gene_prefixes_do_not_share_mutations_or_msas(tmp_path):
    source = tmp_path / "source"
    smn = _gene_tree(source, "SMN", protein_msa=False, codon_msa=False)
    smn2 = _gene_tree(source, "SMN2")

    genes, manifest = _manifest(_args(source, tmp_path / "results"))

    assert genes == ["SMN", "SMN2"]
    assert manifest["msa"] == ["SMN2"]
    assert manifest["codon_msa"] == ["SMN2"]
    assert manifest["input_files"]["SMN"].get("msa") is None
    assert manifest["input_files"]["SMN"].get("codon_msa") is None
    assert manifest["input_files"]["SMN"]["mutations"] == str(smn["mutations"].resolve())
    assert manifest["input_files"]["SMN2"]["mutations"] == str(smn2["mutations"].resolve())


def test_manifest_paths_are_absolute_for_relative_cli_inputs(tmp_path, monkeypatch):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    monkeypatch.chdir(tmp_path)

    _, manifest = _manifest(_args(Path("source"), Path("results")))

    assert all(Path(path).is_absolute()
               for path in manifest["input_files"]["NPM1"].values())


def test_missing_mutations_are_not_replaced_by_another_mapping(tmp_path):
    source = tmp_path / "source"
    paths = _gene_tree(source, "NPM1")
    paths["mutations"].unlink()
    _write(source / "NPM1" / "mappings" / "aa" / "NPM1_aa_mapping.csv",
           "mutant,aa\nG4A,A2T\n")

    _, manifest = _manifest(_args(source, tmp_path / "results"))

    assert manifest["input_files"]["NPM1"]["mutations"] is None


def test_existing_parameter_and_tsv_inventory_is_retained(tmp_path):
    source = tmp_path / "source"
    _gene_tree(source, "NPM1")
    output = tmp_path / "results"
    _write(output / "model_params" / "NPM1.model_params", "params")
    _write(output / "NPM1" / "EVmutation" / "NPM1.protein.tsv", "pkey\nNPM1-test\n")

    _, manifest = _manifest(_args(source, output))

    assert manifest["model_params"] == ["NPM1"]
    assert manifest["EVmutation"] == ["NPM1"]
    assert manifest["codon_EVmutation"] == []
