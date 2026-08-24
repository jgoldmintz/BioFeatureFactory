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
Controller for the mutation-effects Nextflow pipeline (bin/main.nf).

Drives two parallel DCA backends over each gene:
  EVmutation/plmc  — pseudolikelihood Potts inference (CPU, external plmc binary)
  adabmDCA         — Boltzmann / pseudolikelihood Potts inference, run in-process
                     via the adabmDCApy Python API (GPU); see adabmdca_pipeline.py

Each backend scores a protein side (missense/stop) and a codon side
(synonymous). The full chain per gene:
  1. Protein MSA (jackhmmer -> UniRef90)  — skipped if pre-built
  2. Codon MSA (mmseqs2 -> MAFFT)         — skipped if pre-built
  3. Scoring: EVmutation (plmc) and/or adabmDCA, each protein + codon

Before launching Nextflow, inventories existing artifacts per gene
(MSAs, plmc + adabmDCA params, prior TSVs) and writes a manifest so
Nextflow knows what to generate vs skip for each gene. Backend selection
is via --evmutation-only / --adabmdca-only; codon side via --skip-codon.

Use --resume to continue from a previous Nextflow run.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

from biofeaturefactory.lib.utility import (
    discover_fasta_files,
    extract_gene_from_filename,
)
from biofeaturefactory.mutation_effects.evmutation_pipeline import _find_file_for_gene

HERE = Path(__file__).resolve().parent
NEXTFLOW_SCRIPT = HERE / "bin" / "main.nf"


def _resolve_genes(fasta_path: Path) -> List[str]:
    """
    Resolve gene list from --fasta, accepting either form:
      - Directory: enumerate via discover_fasta_files (recursive glob).
      - Single FASTA file: derive a single gene name from the filename via
        extract_gene_from_filename (same logic the rest of BFF uses).
    """
    if fasta_path.is_dir():
        return list(discover_fasta_files(str(fasta_path)).keys())
    if fasta_path.is_file():
        gene = extract_gene_from_filename(fasta_path.stem)
        return [gene] if gene else []
    return []


def build_manifest(genes: List[str], args: argparse.Namespace) -> Dict[str, List[str]]:
    """
    Inventory existing artifacts per gene.

    Returns dict mapping artifact type to list of genes that have it.
    """
    out_dir = args.output

    # Paths match publishDir structure from main.nf:
    #   {output}/MSA/{GENE}.msa.a2m
    #   {output}/CodonMSA/{GENE}.codon.msa.fasta
    #   {output}/model_params/{GENE}.model_params
    #   {output}/codon_model_params/{GENE}.codon_model_params
    #   {output}/{GENE}/EVmutation/{GENE}.protein.tsv
    #   {output}/{GENE}/EVmutation/{GENE}.codon.tsv
    #
    # --msa / --codon-msa / --model-params / --codon-model-params can override
    # the MSA/params locations; final TSVs always come from output dir.

    artifact_checks = {
        "msa": (args.msa or out_dir / "MSA", ["*.a2m", "*.msa.a2m"]),
        "codon_msa": (args.codon_msa or out_dir / "CodonMSA", ["*.codon.msa.fasta"]),
        "model_params": (args.model_params or out_dir / "model_params", ["*.model_params"]),
        "codon_model_params": (args.codon_model_params or out_dir / "codon_model_params", ["*.codon_model_params"]),
        # adabmDCA params (opt-in via --run-adabmdca; manifest entries always tracked
        # so resume works regardless of which run produced the artifact)
        "adabmdca_protein_params": (
            args.adabmdca_protein_params or out_dir / "adabmdca_protein_params",
            ["*.protein_adabm_params", "*.protein.dat", "*.dat"],
        ),
        "adabmdca_codon_params": (
            args.adabmdca_codon_params or out_dir / "adabmdca_codon_params",
            ["*.codon_adabm_params", "*.codon.dat", "*.dat"],
        ),
    }

    # Final TSVs:
    #   EVmutation:  {output}/{GENE}/EVmutation/{GENE}.{protein,codon}.tsv
    #   adabmDCA:    {output}/{GENE}/adabmDCA/{GENE}.{protein,codon}.tsv
    tsv_checks = {
        "EVmutation":         ("EVmutation", "protein.tsv"),
        "codon_EVmutation":   ("EVmutation", "codon.tsv"),
        "adabmdca_protein":   ("adabmDCA",   "protein.tsv"),
        "adabmdca_codon":     ("adabmDCA",   "codon.tsv"),
    }

    manifest = {key: [] for key in list(artifact_checks) + list(tsv_checks)}

    for gene in genes:
        for artifact, (path, patterns) in artifact_checks.items():
            if Path(str(path)).exists() and _find_file_for_gene(gene, str(path), patterns):
                manifest[artifact].append(gene)

        for tsv_key, (subdir, suffix) in tsv_checks.items():
            if (out_dir / gene / subdir / f"{gene}.{suffix}").exists():
                manifest[tsv_key].append(gene)

    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Controller for the EVmutation Nextflow pipeline."
    )

    # Required
    parser.add_argument("-f", "--fasta", type=Path, required=True,
                        help="ORF FASTA file or directory of per-gene FASTA files")
    parser.add_argument("-m", "--mutations", type=Path, required=True,
                        help="Mutations CSV file or directory of per-gene CSVs")
    parser.add_argument("-pb", "--plmc-binary", type=str,
                        help="Path to plmc executable (required when the EVmutation backend runs; "
                             "omit only when --adabmdca-only is set)")
    parser.add_argument("-dr", "--db-root", type=Path,
                        help="Bio_DBs root directory (contains uniref90.fasta, refseq_assemblies/, etc.). "
                             "Required only when at least one gene needs MSA generation; omit when "
                             "--msa and --codon-msa cover every gene in --fasta.")
    parser.add_argument("--output", "-o", type=Path, default=Path("."),
                        help="Output base directory")

    # Pre-built MSA sources (optional, skips generation for genes that have them)
    parser.add_argument("-ms", "--msa", type=Path,
                        help="Protein MSA file or directory")
    parser.add_argument("-cm", "--codon-msa", type=Path,
                        help="Codon MSA file or directory")

    # Tool binaries
    parser.add_argument("-jb", "--jackhmmer-binary", type=str, default="jackhmmer",
                        help="Path to jackhmmer (default: jackhmmer)")
    parser.add_argument("-ji", "--jackhmmer-iterations", type=int, default=5,
                        help="Number of jackhmmer iterations (default: 5)")
    parser.add_argument("-mb", "--mmseqs-binary", type=str, default="mmseqs",
                        help="Path to mmseqs2 (default: mmseqs)")
    parser.add_argument("-a", "--aligner", type=str, default="mafft",
                        choices=["mafft", "muscle"],
                        help="Protein aligner (default: mafft)")

    # Pre-built model params
    parser.add_argument("-mp", "--model-params", type=Path,
                        help="Protein model params file or directory")
    parser.add_argument("-cmp", "--codon-model-params", type=Path,
                        help="Codon model params file or directory")

    # Backend selection — both run in parallel by default.
    # --evmutation-only and --adabmdca-only are mutually exclusive.
    backend_group = parser.add_mutually_exclusive_group()
    backend_group.add_argument("-eo", "--evmutation-only", action="store_true",
                               help="Run only the EVmutation/plmc backend; skip adabmDCA.")
    backend_group.add_argument("-ao", "--adabmdca-only", action="store_true",
                               help="Run only the adabmDCA backend; skip EVmutation/plmc.")

    # --skip-codon: optional argument selects which backend(s) to skip codon-side for.
    #   --skip-codon                  → both backends skip codon-side
    #   --skip-codon evmutation       → only EVmutation skips codon-side
    #   --skip-codon adabmdca         → only adabmDCA skips codon-side
    parser.add_argument("-sc", "--skip-codon", nargs="?", const="both", default=None,
                        choices=["both", "evmutation", "adabmdca"],
                        help="Skip codon-side scoring. No arg = both backends; "
                             "'evmutation' or 'adabmdca' targets one backend. "
                             "Synonymous + stop variants route to the protein TSV instead.")

    # adabmDCA pre-built params
    parser.add_argument("-app", "--adabmdca-protein-params", type=Path,
                        help="Pre-built adabmDCA protein params file or directory "
                             "(default: <output>/adabmdca_protein_params/{GENE}.protein_adabm_params)")
    parser.add_argument("-acp", "--adabmdca-codon-params", type=Path,
                        help="Pre-built adabmDCA codon params file or directory "
                             "(default: <output>/adabmdca_codon_params/{GENE}.codon_adabm_params)")

    # adabmDCA tunables (consulted only when the adabmDCA backend runs).
    # NOTE: there is no --adabmdca-binary flag — the `adabmDCA` console script
    # is installed by `pip install adabmDCA` (the Python/torch implementation,
    # not a compiled binary) and is expected on PATH.
    parser.add_argument("-am", "--adabmdca-model", default="bmDCA",
                        choices=["bmDCA", "eaDCA", "edDCA", "pseudoDCA"],
                        help="Training routine (default: bmDCA). bmDCA/eaDCA/edDCA are "
                             "Boltzmann-learning variants (high memory); pseudoDCA is "
                             "pseudolikelihood (~2x less peak GPU memory, no MCMC). The "
                             "Boltzmann path emits an OOM-triggered hint suggesting "
                             "pseudoDCA when memory becomes the bottleneck.")
    # None => adabmdca_pipeline.py picks per backend (500 pseudoDCA / 50000 Boltzmann).
    # A single 50000 default silently multiplied the pseudoDCA path by 100x.
    parser.add_argument("-an", "--adabmdca-nepochs", type=int, default=None,
                        help="Max epochs. Default: 500 for pseudoDCA, 50000 otherwise.")
    parser.add_argument("-at", "--adabmdca-tol", type=float, default=1e-3,
                        help="pseudoDCA convergence threshold on ||grad||/||grad||_0 (default: 1e-3; 0 disables)")
    parser.add_argument("-ap", "--adabmdca-patience", type=int, default=3)
    parser.add_argument("-ace", "--adabmdca-check-every", type=int, default=10)
    parser.add_argument("-ata", "--adabmdca-target", type=float, default=0.95,
                        help="Pearson Cij target (default: 0.95)")
    parser.add_argument("-al", "--adabmdca-lr", type=float, default=0.01)
    parser.add_argument("-anc", "--adabmdca-nchains", type=int, default=10000,
                        help="PCD chain count (default: 10000)")
    parser.add_argument("-ans", "--adabmdca-nsweeps", type=int, default=10)
    parser.add_argument("-ad", "--adabmdca-device", default="cuda",
                        help="adabmDCA device (default: cuda)")
    parser.add_argument("-adt", "--adabmdca-dtype", default="float32",
                        choices=["float32", "float64"])
    parser.add_argument("-as", "--adabmdca-seed", type=int, default=0)

    # Options
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("-vl", "--validation-log", type=Path)
    parser.add_argument("-r", "--resume", action="store_true",
                        help="Resume previous Nextflow run")

    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args: argparse.Namespace) -> None:
    # Derived backend flags (computed once, reused everywhere downstream).
    args.run_evmutation = not args.adabmdca_only
    args.run_adabmdca = not args.evmutation_only

    # --skip-codon → per-backend booleans.
    args.skip_codon_evmutation = args.skip_codon in ("both", "evmutation")
    args.skip_codon_adabmdca = args.skip_codon in ("both", "adabmdca")

    genes = _resolve_genes(args.fasta)
    if not genes:
        raise SystemExit(f"ERROR: No FASTA files found in {args.fasta}")

    # plmc is only needed to BUILD params. Providing --model-params (protein) /
    # --codon-model-params (codon) lets evmutation_pipeline.py skip plmc and score
    # directly, so the binary is required only when a side still builds params.
    need_protein_plmc = args.run_evmutation and not args.model_params
    need_codon_plmc = args.run_evmutation and not args.skip_codon_evmutation and not args.codon_model_params
    if (need_protein_plmc or need_codon_plmc) and not args.plmc_binary:
        raise SystemExit("ERROR: --plmc-binary is required when the EVmutation backend builds params "
                         "(provide --model-params / --codon-model-params to skip plmc, or --adabmdca-only)")

    # --db-root is validated at runtime in validate_db_coverage(), once the
    # manifest has been built and we know which genes still need MSA generation.


def validate_db_coverage(genes: List[str], manifest: Dict[str, List[str]],
                         args: argparse.Namespace) -> None:
    """
    Decide whether --db-root is needed based on MSA coverage from the manifest.

    Called after build_manifest. Errors only if at least one gene still needs
    MSA generation (protein or codon) and --db-root wasn't supplied / doesn't
    contain the required DB.
    """
    have_protein = set(manifest.get("msa", []))
    have_codon   = set(manifest.get("codon_msa", []))
    need_protein_gen = [g for g in genes if g not in have_protein]
    need_codon_gen   = [g for g in genes if g not in have_codon]

    if not need_protein_gen and not need_codon_gen:
        return  # all MSAs pre-built; --db-root genuinely unnecessary

    if not args.db_root:
        missing = []
        if need_protein_gen:
            missing.append(f"protein MSA for {len(need_protein_gen)} gene(s) "
                           f"(e.g. {', '.join(need_protein_gen[:3])})")
        if need_codon_gen:
            missing.append(f"codon MSA for {len(need_codon_gen)} gene(s) "
                           f"(e.g. {', '.join(need_codon_gen[:3])})")
        raise SystemExit(
            "ERROR: --db-root is required because "
            + " and ".join(missing)
            + " need to be generated. Provide --db-root, or supply pre-built MSAs "
              "via --msa / --codon-msa that cover every gene."
        )

    if need_protein_gen:
        uniref90    = args.db_root / "uniref90.fasta"
        uniref90_gz = args.db_root / "uniref90.fasta.gz"
        if not uniref90.exists() and not uniref90_gz.exists():
            raise SystemExit(
                f"ERROR: uniref90.fasta(.gz) not found in --db-root ({args.db_root}) "
                f"but {len(need_protein_gen)} gene(s) still need protein MSA generation."
            )


def normalize(path: Optional[Path]) -> Optional[str]:
    return str(path.resolve()) if path else None


def build_nextflow_cmd(args: argparse.Namespace, manifest_path: str) -> List[str]:
    cmd = ["nextflow", "run", str(NEXTFLOW_SCRIPT)]
    if args.resume:
        cmd.append("-resume")

    def add_param(name: str, value):
        if value is not None:
            cmd.extend([f"--{name}", str(value)])

    add_param("fasta", normalize(args.fasta))
    add_param("mutations", normalize(args.mutations))
    if args.plmc_binary:
        add_param("plmc_binary", normalize(Path(args.plmc_binary)))
    if args.db_root:
        add_param("db_root", normalize(args.db_root))
        add_param("uniref90_db", normalize(args.db_root / "uniref90.fasta"))
    add_param("output_dir", normalize(args.output))
    add_param("threads", args.threads)
    add_param("manifest", manifest_path)

    add_param("jackhmmer_binary", args.jackhmmer_binary)
    add_param("jackhmmer_iterations", args.jackhmmer_iterations)
    add_param("mmseqs_binary", args.mmseqs_binary)
    add_param("aligner", args.aligner)

    if args.msa:
        add_param("msa", normalize(args.msa))
    if args.codon_msa:
        add_param("codon_msa", normalize(args.codon_msa))
    if args.model_params:
        add_param("model_params", normalize(args.model_params))
    if args.codon_model_params:
        add_param("codon_model_params", normalize(args.codon_model_params))

    # Backend skip flags — only emit "true" when the user wants the skip on.
    # Never emit "false" (Groovy's `if ("false")` is truthy and would invert the gate).
    if not args.run_evmutation:
        add_param("skip_evmutation", "true")
    if not args.run_adabmdca:
        add_param("skip_adabmdca", "true")

    # Codon-side skip flags — same string-truthy convention.
    if args.skip_codon_evmutation:
        add_param("skip_codon_evmutation", "true")
    if args.skip_codon_adabmdca:
        add_param("skip_codon_adabmdca", "true")

    if args.adabmdca_protein_params:
        add_param("adabmdca_protein_params", normalize(args.adabmdca_protein_params))
    if args.adabmdca_codon_params:
        add_param("adabmdca_codon_params", normalize(args.adabmdca_codon_params))

    # adabmDCA tunables — only forward when the adabmDCA backend runs
    if args.run_adabmdca:
        add_param("adabmdca_model",   args.adabmdca_model)
        if args.adabmdca_nepochs is not None:
            add_param("adabmdca_nepochs", args.adabmdca_nepochs)
        add_param("adabmdca_tol",         args.adabmdca_tol)
        add_param("adabmdca_patience",    args.adabmdca_patience)
        add_param("adabmdca_check_every", args.adabmdca_check_every)
        add_param("adabmdca_target",  args.adabmdca_target)
        add_param("adabmdca_lr",      args.adabmdca_lr)
        add_param("adabmdca_nchains", args.adabmdca_nchains)
        add_param("adabmdca_nsweeps", args.adabmdca_nsweeps)
        add_param("adabmdca_device",  args.adabmdca_device)
        add_param("adabmdca_dtype",   args.adabmdca_dtype)
        add_param("adabmdca_seed",    args.adabmdca_seed)

    if args.validation_log:
        add_param("validation_log", normalize(args.validation_log))

    return cmd


def validate_backend_tools(args: argparse.Namespace, manifest: Dict[str, List[str]],
                            genes: List[str]) -> None:
    """
    Pre-flight check that adabmDCApy is importable BEFORE launching the workflow.

    adabmdca_pipeline.py drives training IN-PROCESS via the adabmDCApy Python API
    (`from adabmDCA.training import ...`), NOT by shelling out to the `adabmDCA`
    console script. The correct preflight therefore resolves the package import,
    not a binary on PATH: a console script can be present while the package fails
    to import, and the in-process path needs the import regardless of PATH.
    Surfacing this here avoids a cryptic ImportError deep inside a Nextflow task.

    Skipped when every gene that still needs inference already has pre-built params.
    """
    if not args.run_adabmdca:
        return

    need_protein_inf = [g for g in genes if g not in manifest.get("adabmdca_protein_params", [])]
    need_codon_inf = [] if args.skip_codon_adabmdca else \
        [g for g in genes if g not in manifest.get("adabmdca_codon_params", [])]

    if not need_protein_inf and not need_codon_inf:
        return  # all adabmDCA params pre-built; adabmDCApy not needed

    sides = []
    if need_protein_inf: sides.append(f"protein ({len(need_protein_inf)} gene(s))")
    if need_codon_inf:   sides.append(f"codon ({len(need_codon_inf)} gene(s))")
    side_str = " and ".join(sides)

    # adabmdca_pipeline.py imports these in-process when training fires. find_spec
    # resolves the module without importing it (no torch/CUDA init side effects).
    missing = [m for m in ("adabmDCA", "torch") if importlib.util.find_spec(m) is None]
    if missing:
        raise SystemExit(
            f"ERROR: adabmDCA training needs {side_str}, but these Python package(s)\n"
            f"       are not importable in the active interpreter: {', '.join(missing)}.\n"
            f"       adabmdca_pipeline.py trains in-process via the adabmDCApy API, so the\n"
            f"       package must import here — a console script on PATH is neither\n"
            f"       sufficient nor required. Install into THIS interpreter:\n"
            f"         pip install adabmDCA torch\n"
            f"       Then verify: python -c 'import adabmDCA, torch'"
        )


def run_controller(args: argparse.Namespace):
    genes = _resolve_genes(args.fasta)
    manifest = build_manifest(genes, args)
    validate_db_coverage(genes, manifest, args)
    validate_backend_tools(args, manifest, genes)

    # Summary
    backends = []
    if args.run_evmutation: backends.append("EVmutation")
    if args.run_adabmdca:   backends.append("adabmDCA")
    flags = [f"backends={'+'.join(backends)}"]
    if args.skip_codon_evmutation and args.skip_codon_adabmdca:
        flags.append("skip-codon=both")
    elif args.skip_codon_evmutation:
        flags.append("skip-codon=evmutation")
    elif args.skip_codon_adabmdca:
        flags.append("skip-codon=adabmdca")
    print(f"[mutEffects-controller] {len(genes)} gene(s) [{', '.join(flags)}]", flush=True)
    for artifact, gene_list in manifest.items():
        if gene_list:
            preview = ", ".join(gene_list[:5])
            print(f"  {artifact}: {len(gene_list)} pre-built ({preview})")

    # "Ready to score" set depends on which backends are enabled and codon gating.
    needed_artifacts: List[str] = []
    if args.run_evmutation:
        needed_artifacts.append("model_params")
        if not args.skip_codon_evmutation:
            needed_artifacts.append("codon_model_params")
    if args.run_adabmdca:
        needed_artifacts.append("adabmdca_protein_params")
        if not args.skip_codon_adabmdca:
            needed_artifacts.append("adabmdca_codon_params")
    if needed_artifacts:
        ready = [g for g in genes if all(g in manifest[a] for a in needed_artifacts)]
        if ready:
            print(f"  ready to score: {len(ready)} gene(s)")

    # Write manifest
    out_dir = args.output.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = str(out_dir / ".evmutation_manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    cmd = build_nextflow_cmd(args, manifest_path)
    print(f"[mutEffects-controller] Launching: {' '.join(cmd)}", flush=True)
    nf_proc = subprocess.Popen(cmd, cwd=str(HERE))
    sys.exit(nf_proc.wait())


if __name__ == "__main__":
    run_controller(parse_args())