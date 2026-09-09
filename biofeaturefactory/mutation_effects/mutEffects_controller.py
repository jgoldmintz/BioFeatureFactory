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
  EVmutation/plmc  -- pseudolikelihood Potts inference (CPU, external plmc binary)
  adabmDCA         -- Boltzmann / pseudolikelihood Potts inference, run in-process
                     via the adabmDCApy Python API (GPU); see adabmdca_pipeline.py

Each backend has a protein side and a codon side. Without explicit MSA flags,
mutation classes select the required sides per gene. The full chain per gene:
  1. Protein MSA (jackhmmer -> UniRef90)  -- skipped if pre-built
  2. Codon MSA (mmseqs2 -> MAFFT)         -- skipped if pre-built
  3. Scoring: EVmutation (plmc) and/or adabmDCA on selected sides as they become ready

Before launching Nextflow, inventories existing artifacts per gene
(MSAs, plmc + adabmDCA params, prior TSVs) and writes a manifest so
Nextflow knows what to generate vs skip for each gene. Backend selection
is via --evmutation-only / --adabmdca-only; codon side via --skip-codon.

Use --resume to continue from a previous Nextflow run.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from biofeaturefactory.lib.utility import (
    derive_mutations_root,
    discover_fasta_files,
    discover_mutation_files,
    extract_gene_from_filename,
    find_gene_file,
)
from biofeaturefactory.mutation_effects.bin import evmutation_cache, plmc_resources, resource_planner
from biofeaturefactory.mutation_effects.bin.adabmdca_task import verify_completion
from biofeaturefactory.mutation_effects.bin.mutation_routing import classify_gene, choose_route

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


def build_manifest(genes: List[str], args: argparse.Namespace) -> Dict[str, Any]:
    """
    Inventory existing artifacts per gene.

    Returns artifact gene lists and resolved per-gene input_files for Nextflow.
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
    fasta_files = (
        discover_fasta_files(str(args.fasta)) if args.fasta.is_dir()
        else {gene: str(args.fasta) for gene in genes}
    )
    mutation_files = discover_mutation_files(args.mutations)
    input_files = {}
    param_files = {}

    for gene in genes:
        param_files[gene] = {}
        mutation_file = (
            args.mutations if args.mutations.is_file()
            else mutation_files.get(gene)
        )
        input_files[gene] = {
            "fasta": normalize(Path(fasta_files[gene])),
            "mutations": normalize(Path(mutation_file)) if mutation_file else None,
        }
        for artifact, (path, patterns) in artifact_checks.items():
            candidates = [path]
            if artifact in ("msa", "codon_msa") and not getattr(args, artifact):
                candidates.append(out_dir)
                if args.fasta.is_dir():
                    candidates.extend([args.fasta, args.fasta / path.name])
            resolved = next((
                match for candidate in candidates
                if (match := find_gene_file(str(candidate), gene, patterns))
            ), None)
            if resolved:
                manifest[artifact].append(gene)
                if artifact in ("msa", "codon_msa"):
                    input_files[gene][artifact] = normalize(Path(resolved))
                else:
                    param_files[gene][artifact] = normalize(Path(resolved))

        for tsv_key, (subdir, suffix) in tsv_checks.items():
            if (out_dir / gene / subdir / f"{gene}.{suffix}").exists():
                manifest[tsv_key].append(gene)

    manifest["input_files"] = input_files
    manifest["param_files"] = param_files
    manifest["routing"] = {}
    manifest["ev_fingerprints"] = {}
    for gene, inputs in input_files.items():
        classes = classify_gene(inputs["fasta"], inputs["mutations"], gene,
                                getattr(args, "validation_log", None)) if inputs["mutations"] else []
        route = choose_route(classes, protein_explicit=bool(args.msa),
                             codon_explicit=bool(args.codon_msa))
        if route["codon"]:
            for backend in ("evmutation", "adabmdca"):
                if getattr(args, f"run_{backend}") and getattr(args, f"skip_codon_{backend}"):
                    route["warnings"].append(
                        f"--skip-codon overrides codon processing for {backend}; "
                        "this will not produce biologically accurate results for synonymous/stop-codon effects."
                    )
        manifest["routing"][gene] = route
        manifest["ev_fingerprints"][gene] = {}
        for side in ("protein", "codon"):
            if not side_enabled(args, manifest, gene, "evmutation", side):
                continue
            fingerprint = evmutation_cache.routing_fingerprint(
                gene, side, inputs["fasta"], inputs["mutations"],
                validation_log=getattr(args, "validation_log", None),
                skip_codon=side == "protein" and not side_enabled(args, manifest, gene, "evmutation", "codon"),
                score_missense_codon=side == "codon" and route["score_missense_codon"],
            )
            manifest["ev_fingerprints"][gene][side] = fingerprint
            artifact = "EVmutation" if side == "protein" else "codon_EVmutation"
            params_artifact = "model_params" if side == "protein" else "codon_model_params"
            model_params = param_files[gene].get(params_artifact)
            alignment = inputs.get("msa" if side == "protein" else "codon_msa")
            folder = out_dir / gene / "EVmutation"
            if gene in manifest[artifact] and not (model_params and alignment and evmutation_cache.verify_completion(
                folder / f"{gene}.{side}.routing.json", fingerprint,
                folder / f"{gene}.{side}.tsv", model_params, msa=alignment,
            )):
                manifest[artifact].remove(gene)
    return manifest


def side_enabled(args, manifest, gene, backend, side):
    if not getattr(args, f"run_{backend}"):
        return False
    route = manifest.get("routing", {}).get(gene, {"protein": True, "codon": True})
    if getattr(args, f"skip_codon_{backend}"):
        return side == "protein" and (route["protein"] or route["codon"])
    return route[side]


def side_pending(args, manifest, gene, backend, side):
    artifact = ("EVmutation" if side == "protein" else "codon_EVmutation") if backend == "evmutation" else f"adabmdca_{side}"
    blocked = any(
        error["gene"] == gene and error["backend"] == backend and error["side"] == side
        for error in manifest.get("resource_errors", [])
    )
    return not blocked and side_enabled(args, manifest, gene, backend, side) and gene not in manifest.get(artifact, [])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Controller for the EVmutation Nextflow pipeline."
    )

    # Required
    parser.add_argument("-f", "--fasta", type=Path, required=True,
                        help="ORF FASTA file or directory of per-gene FASTA files")
    parser.add_argument("-m", "--mutations", type=Path, required=False,
                        help="Mutations CSV file or directory of per-gene CSVs")
    parser.add_argument("-pb", "--plmc-binary", type=str,
                        help="Path to plmc executable (required when the EVmutation backend runs; "
                             "omit only when --adabmdca-only is set)")
    parser.add_argument("-dr", "--db-root", type=Path,
                        help="Bio_DBs root directory (contains uniref90.fasta, refseq_assemblies/, etc.). "
                             "Required only when at least one gene needs MSA generation; omit when "
                             "pre-built MSAs cover every gene in --fasta.")
    parser.add_argument("--output", "-o", type=Path, default=Path("."),
                        help="Output base directory")

    # Pre-built MSA sources (optional, skips generation for genes that have them)
    parser.add_argument("-ms", "--msa", type=Path,
                        help="Protein MSA source; without -cm explicitly selects protein-only mode. "
                             "Without either flag, select per gene and discover MSAs under --output/--fasta.")
    parser.add_argument("-cm", "--codon-msa", type=Path,
                        help="Codon MSA source; without --msa explicitly selects codon-only mode, "
                             "including codon-level missense scoring. Both flags enable both sides.")

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

    # Backend selection -- both run in parallel by default.
    # --evmutation-only and --adabmdca-only are mutually exclusive.
    backend_group = parser.add_mutually_exclusive_group()
    backend_group.add_argument("-eo", "--evmutation-only", action="store_true",
                               help="Run only the EVmutation/plmc backend; skip adabmDCA.")
    backend_group.add_argument("-ao", "--adabmdca-only", action="store_true",
                               help="Run only the adabmDCA backend; skip EVmutation/plmc.")

    # --skip-codon: optional argument selects which backend(s) to skip codon-side for.
    #   --skip-codon                  -> both backends skip codon-side
    #   --skip-codon evmutation       -> only EVmutation skips codon-side
    #   --skip-codon adabmdca         -> only adabmDCA skips codon-side
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
    # NOTE: there is no --adabmdca-binary flag -- the `adabmDCA` console script
    # is installed by `pip install adabmDCA` (the Python/torch implementation,
    # not a compiled binary) and is expected on PATH.
    parser.add_argument("-am", "--adabmdca-model", default="pseudoDCA",
                        choices=["bmDCA", "eaDCA", "edDCA", "pseudoDCA"],
                        help="Training routine (default: pseudoDCA). bmDCA/eaDCA/edDCA are "
                             "Boltzmann-learning variants (high memory); pseudoDCA is "
                             "pseudolikelihood (no MCMC). Automatic device fallback "
                             "preserves the selected model.")
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
                        help="Boltzmann-only PCD chain count (default: 10000; unused by pseudoDCA)")
    parser.add_argument("-ans", "--adabmdca-nsweeps", type=int, default=10,
                        help="Boltzmann-only sweeps per step (default: 10; unused by pseudoDCA)")
    parser.add_argument("-ad", "--adabmdca-device", default="auto",
                        help="auto routes to eligible GPUs or CPU; cpu/cuda/cuda:N force placement")
    parser.add_argument("-adt", "--adabmdca-dtype", default="float32",
                        choices=["float32", "float64"])
    parser.add_argument("-as", "--adabmdca-seed", type=int, default=0)

    # Options
    parser.add_argument("-t", "--threads", type=int,
                        help="Threads per task (default: share usable CPUs across estimated concurrent jobs)")
    parser.add_argument("-vl", "--validation-log", type=Path)
    parser.add_argument("-r", "--resume", action="store_true",
                        help="Resume previous Nextflow run")
    parser.add_argument("--resource-cpus", type=int,
                        help="Maximum CPUs shared by all Nextflow tasks")
    parser.add_argument("--resource-memory-gib", type=float,
                        help="Maximum shared host RAM budget in GiB")
    parser.add_argument("--resource-headroom", type=float, default=0.9,
                        help="Usable fraction of available RAM and total GPU VRAM (default: 0.9)")
    parser.add_argument("--resource-memory-margin", type=float, default=1.15,
                        help="Safety multiplier for uncalibrated memory estimates (default: 1.15)")
    parser.add_argument("--resource-overrides", type=Path,
                        help="JSON per gene.side with measured gpu_memory_gib/cpu_memory_gib/gpu_host_memory_gib/threads")
    parser.add_argument("--resource-hardware", type=Path,
                        help="Explicit allocated hardware JSON instead of local autodetection")
    parser.add_argument("--resource-plan-only", action="store_true",
                        help="Print placement estimates without launching tasks or writing outputs")
    parser.add_argument("--gpu-lease-dir", type=Path,
                        default=Path.home() / ".cache" / "biofeaturefactory" / "gpu-leases",
                        help="Shared GPU coordination directory; use the same directory across local runs")
    parser.add_argument("--gpu-wait-timeout", type=float, default=600,
                        help="Maximum seconds waiting for an eligible GPU held by another task")
    parser.add_argument("--msa-memory-gib", type=float, default=8,
                        help="Host RAM request per MSA generation task; tune for the databases used")
    parser.add_argument("--evmutation-memory-gib", type=float,
                        help="Minimum host RAM request per EVmutation task in GiB; default: automatic model/workspace estimate")

    args = parser.parse_args()
    if args.adabmdca_device not in {"auto", "cpu", "cuda"}:
        if not args.adabmdca_device.startswith("cuda:") or not args.adabmdca_device[5:].isdigit():
            parser.error("--adabmdca-device must be auto, cpu, cuda or cuda:N")
    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be positive")
    if not 0 < args.resource_headroom < 1:
        parser.error("--resource-headroom must be between zero and one")
    for name in ("resource_memory_margin", "gpu_wait_timeout", "msa_memory_gib", "evmutation_memory_gib"):
        if getattr(args, name) is None:
            continue
        try:
            resource_planner.positive(getattr(args, name), name)
        except ValueError as error:
            parser.error(str(error))
    if args.resource_memory_margin < 1:
        parser.error("--resource-memory-margin cannot be less than one")

    # One root supplies both; see lib/utility.derive_mutations_root.

    args.mutations = derive_mutations_root(args.mutations, args.fasta, label="mutEffects")

    if not args.mutations:

        parser.error("--mutations is required (no <GENE>/mappings/mutations/ under "

                     f"--fasta {args.fasta})")

    validate_args(args)
    return args


def validate_args(args: argparse.Namespace) -> None:
    # Derived backend flags (computed once, reused everywhere downstream).
    args.run_evmutation = not args.adabmdca_only
    args.run_adabmdca = not args.evmutation_only

    # --skip-codon -> per-backend booleans.
    args.skip_codon_evmutation = args.skip_codon in ("both", "evmutation")
    args.skip_codon_adabmdca = args.skip_codon in ("both", "adabmdca")

    genes = _resolve_genes(args.fasta)
    if not genes:
        raise SystemExit(f"ERROR: No FASTA files found in {args.fasta}")

    # --db-root is validated at runtime in validate_db_coverage(), once the
    # manifest has been built and we know which genes still need MSA generation.


def validate_db_coverage(genes: List[str], manifest: Dict[str, Any],
                         args: argparse.Namespace) -> None:
    """
    Decide whether --db-root is needed based on MSA coverage from the manifest.

    Called after build_manifest. Errors only if at least one gene still needs
    MSA generation (protein or codon) and --db-root wasn't supplied / doesn't
    contain the required DB.
    """
    have_protein = set(manifest.get("msa", []))
    have_codon   = set(manifest.get("codon_msa", []))
    need_protein_gen = [gene for gene in genes if gene not in have_protein and any(
        side_pending(args, manifest, gene, backend, "protein") for backend in ("evmutation", "adabmdca")
    )]
    need_codon_gen = [gene for gene in genes if gene not in have_codon and any(
        side_pending(args, manifest, gene, backend, "codon") for backend in ("evmutation", "adabmdca")
    )]

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
        if not uniref90.is_file() or uniref90.stat().st_size == 0:
            if uniref90_gz.is_file():
                destination_note = (
                    f"Move aside the empty or non-file destination {shlex.quote(str(uniref90))} first. "
                    if uniref90.exists() else ""
                )
                raise SystemExit(
                    "ERROR: jackhmmer requires an uncompressed, rewindable UniRef90 FASTA. "
                    f"Only gzip data is available in --db-root ({args.db_root}). "
                    f"{destination_note}"
                    "Prepare uniref90.fasta before rerunning (allow space for the expanded database):\n"
                    f"  gzip -dk {shlex.quote(str(uniref90_gz.resolve()))}"
                )
            raise SystemExit(
                f"ERROR: nonempty uniref90.fasta not found in --db-root ({args.db_root}) "
                f"but {len(need_protein_gen)} gene(s) still need protein MSA generation."
            )


def normalize(path: Optional[Path]) -> Optional[str]:
    return str(path.resolve()) if path else None


def build_nextflow_cmd(args: argparse.Namespace, manifest_path: str) -> List[str]:
    cmd = ["nextflow", "run", str(NEXTFLOW_SCRIPT)]
    if getattr(args, "resource_runtime_config", None):
        cmd.extend(["-c", str(args.resource_runtime_config)])
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
        uniref90 = args.db_root / "uniref90.fasta"
        if uniref90.is_file() and uniref90.stat().st_size > 0:
            add_param("uniref90_db", normalize(uniref90))
    add_param("output_dir", normalize(args.output))
    threads = getattr(args, "resource_config", {}).get("threads", args.threads)
    add_param("threads", threads)
    add_param("manifest", manifest_path)
    add_param("resource_errors", getattr(args, "resource_errors_path", None))
    if getattr(args, "resource_config_path", None):
        add_param("resource_config", str(args.resource_config_path))
        add_param("resource_executor", "local")
        add_param("gpu_slots", len(args.resource_config["hardware"]["gpus"]))
        add_param("msa_cpus", threads)
        add_param("msa_memory", f"{args.msa_memory_gib} GB")
        add_param("evmutation_cpus", threads)
        if args.evmutation_memory_gib is not None:
            add_param("evmutation_memory", f"{args.evmutation_memory_gib} GB")

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

    # Backend skip flags -- only emit "true" when the user wants the skip on.
    # Never emit "false" (Groovy's `if ("false")` is truthy and would invert the gate).
    if not args.run_evmutation:
        add_param("skip_evmutation", "true")
    if not args.run_adabmdca:
        add_param("skip_adabmdca", "true")

    # Codon-side skip flags -- same string-truthy convention.
    if args.skip_codon_evmutation:
        add_param("skip_codon_evmutation", "true")
    if args.skip_codon_adabmdca:
        add_param("skip_codon_adabmdca", "true")

    if args.adabmdca_protein_params:
        add_param("adabmdca_protein_params", normalize(args.adabmdca_protein_params))
    if args.adabmdca_codon_params:
        add_param("adabmdca_codon_params", normalize(args.adabmdca_codon_params))

    # adabmDCA tunables -- only forward when the adabmDCA backend runs
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


def validate_backend_tools(args: argparse.Namespace, manifest: Dict[str, Any],
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
    need_plmc = any(
        side_pending(args, manifest, gene, "evmutation", side)
        and gene not in manifest.get("model_params" if side == "protein" else "codon_model_params", [])
        for gene in genes for side in ("protein", "codon")
    )
    if need_plmc and not args.plmc_binary:
        raise SystemExit("ERROR: --plmc-binary is required when a selected EVmutation task builds params "
                         "(provide pre-built params or --adabmdca-only)")
    if not args.run_adabmdca:
        return

    need_protein_inf = [gene for gene in genes if side_pending(args, manifest, gene, "adabmdca", "protein")
                        and gene not in manifest.get("adabmdca_protein_params", [])]
    need_codon_inf = [gene for gene in genes if side_pending(args, manifest, gene, "adabmdca", "codon")
                      and gene not in manifest.get("adabmdca_codon_params", [])]

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
            f"       package must import here -- a console script on PATH is neither\n"
            f"       sufficient nor required. Install into THIS interpreter:\n"
            f"         pip install adabmDCA torch\n"
            f"       Then verify: python -c 'import adabmDCA, torch'"
        )


def prepare_resource_plans(args, manifest, genes):
    config = resource_planner.make_config(args)
    config["routing"] = manifest.get("routing", {})
    config["evmutation"] = {"memory_gib": args.evmutation_memory_gib}
    manifest["resource_errors"] = []
    plans = []
    if args.run_adabmdca:
        manifest["adabmdca_protein"] = []
        manifest["adabmdca_codon"] = []
        for side in ("protein", "codon"):
            if side == "codon" and args.skip_codon_adabmdca:
                continue
            explicit_params = getattr(args, f"adabmdca_{side}_params")
            if not explicit_params:
                manifest[f"adabmdca_{side}_params"] = []
            for gene in genes:
                if not side_enabled(args, manifest, gene, "adabmdca", side):
                    continue
                inputs = manifest["input_files"][gene]
                msa = inputs.get("msa" if side == "protein" else "codon_msa")
                if not msa or not inputs["mutations"]:
                    continue
                try:
                    resolved_params = find_gene_file(str(explicit_params), gene, [f"*.{side}_adabm_params", f"*.{side}.dat", "*.dat"]) if explicit_params else None
                    fingerprint = resource_planner.task_fingerprint(gene, side, inputs["fasta"], msa, inputs["mutations"], config, resolved_params)
                    artifact_dir = args.output / gene / "adabmDCA"
                    params_name = f"{gene}.{side}_adabm_params"
                    tsv_name = f"{gene}.{side}.tsv"
                    artifacts = {
                        tsv_name: artifact_dir / tsv_name,
                        params_name: args.output / f"adabmdca_{side}_params" / params_name,
                    }
                    complete = verify_completion(artifact_dir / f"{gene}.{side}.complete.json", fingerprint, artifacts)
                    if complete:
                        manifest[f"adabmdca_{side}"].append(gene)
                        if gene not in manifest[f"adabmdca_{side}_params"]:
                            manifest[f"adabmdca_{side}_params"].append(gene)
                        plan = {"gene": gene, "side": side, "device": "cached", "cpu_memory_gib": 0, "gpu_memory_gib": 0}
                    else:
                        plan = resource_planner.plan_task(gene, side, inputs["fasta"], msa, inputs["mutations"], config, resolved_params)
                except (ValueError, OSError) as error:
                    manifest["resource_errors"].append(resource_planner.resource_error(gene, side, "adabmdca", error))
                    continue
                plan["verified_complete"] = complete
                plans.append(plan)
    config["evmutation_plans"] = []
    for gene in genes:
        for side in ("protein", "codon"):
            if not side_pending(args, manifest, gene, "evmutation", side):
                continue
            inputs = manifest["input_files"][gene]
            msa = inputs.get("msa" if side == "protein" else "codon_msa")
            if not msa or not inputs["mutations"]:
                continue
            artifact = "model_params" if side == "protein" else "codon_model_params"
            params = manifest["param_files"][gene].get(artifact)
            try:
                plan = plmc_resources.plan_evmutation_task(gene, side, msa, config, params)
            except (ValueError, OSError) as error:
                manifest["resource_errors"].append(resource_planner.resource_error(gene, side, "evmutation", error))
                continue
            config["evmutation_plans"].append(plan)
    if args.threads is None:
        tasks = pending_resource_tasks(args, manifest, genes, config, plans)
        config["cpu_allocation"] = resource_planner.automatic_threads(config["hardware"], tasks)
        config["threads"] = config["cpu_allocation"]["threads"]
        for plan in plans:
            if not plan["verified_complete"]:
                override = config["overrides"].get(f"{plan['gene']}.{plan['side']}", {})
                plan["threads"] = override.get("threads", config["threads"])
    else:
        config["cpu_allocation"] = {"threads": config["threads"], "concurrent_jobs": None}
    for plan in config["evmutation_plans"]:
        plan["threads"] = config["threads"]
    return config, plans


def pending_resource_tasks(args, manifest, genes, config, plans):
    planned = {(plan["gene"], plan["side"]): plan for plan in plans}
    ev_planned = {(plan["gene"], plan["side"]): plan for plan in config["evmutation_plans"]}
    tasks = []
    for gene in genes:
        for side in ("protein", "codon"):
            backends = [backend for backend in ("evmutation", "adabmdca")
                        if side_pending(args, manifest, gene, backend, side)]
            if not backends:
                continue
            inputs = manifest["input_files"][gene]
            if not inputs.get("msa" if side == "protein" else "codon_msa"):
                tasks.append({"id": f"{gene}.{side}.msa", "memory_gib": args.msa_memory_gib,
                              "eligible_gpu_uuids": [], "threads": None})
                continue
            for backend in backends:
                task = {"id": f"{gene}.{side}.{backend}", "eligible_gpu_uuids": [], "threads": None}
                if backend == "adabmdca":
                    plan = planned.get((gene, side))
                    if plan is None or plan["device"] == "cached":
                        continue
                    gpu = plan["device"] == "cuda"
                    task["memory_gib"] = plan["gpu_host_memory_gib" if gpu else "cpu_memory_gib"]
                    if gpu:
                        task["eligible_gpu_uuids"] = plan["eligible_gpu_uuids"]
                    task["threads"] = config["overrides"].get(f"{gene}.{side}", {}).get("threads")
                else:
                    plan = ev_planned.get((gene, side))
                    if plan is None:
                        continue
                    task["memory_gib"] = plan["memory_gib"]
                tasks.append(task)
    return tasks


def write_resource_snapshot(directory, contents, suffix):
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / (hashlib.sha256(contents.encode()).hexdigest() + suffix)
    descriptor, temporary = tempfile.mkstemp(dir=directory)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(contents)
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return destination


def print_final_diagnostics(warnings, errors):
    if not warnings and not errors:
        return
    sys.stdout.flush()
    print(f"\n[mutEffects-controller] Final summary: {len(warnings)} warning(s), {len(errors)} error(s)", file=sys.stderr)
    for warning in warnings:
        print(f"WARNING: {warning['gene']}: {warning['message']}", file=sys.stderr)
    for error in errors:
        context = "/".join(str(error.get(field, "")) for field in ("gene", "side", "backend") if error.get(field))
        print(f"ERROR: {context}: {error['message']}", file=sys.stderr)
    sys.stderr.flush()


def read_resource_errors(filename):
    errors = json.loads(Path(filename).read_text())
    if not isinstance(errors, list) or any(
        not isinstance(error, dict)
        or any(not isinstance(error.get(field), str) for field in ("gene", "side", "backend", "message"))
        or error["side"] not in {"protein", "codon"}
        or error["backend"] not in {"evmutation", "adabmdca"}
        for error in errors
    ):
        raise ValueError("Invalid Nextflow resource-error report")
    return errors


def run_controller(args: argparse.Namespace):
    genes = _resolve_genes(args.fasta)
    manifest = build_manifest(genes, args)
    warnings = [{"gene": gene, "message": warning}
                for gene, route in manifest["routing"].items() for warning in route["warnings"]]
    try:
        config, plans = prepare_resource_plans(args, manifest, genes)
    except (ValueError, OSError, KeyError) as error:
        print_final_diagnostics(warnings, manifest.get("resource_errors", []))
        raise SystemExit(f"ERROR: {error}") from error
    args.resource_config = config
    errors = list(manifest["resource_errors"])
    try:
        validate_db_coverage(genes, manifest, args)
    except SystemExit:
        print_final_diagnostics(warnings, errors)
        raise
    if args.resource_plan_only:
        print(json.dumps({"hardware": config["hardware"], "cpu_allocation": config["cpu_allocation"],
                          "routing": manifest["routing"], "tasks": plans + config["evmutation_plans"],
                          "resource_errors": errors, "warnings": warnings}, indent=2), flush=True)
        print_final_diagnostics(warnings, errors)
        if errors:
            raise SystemExit(1)
        return
    pending = any(side_pending(args, manifest, gene, backend, side)
                  for gene in genes for backend in ("evmutation", "adabmdca") for side in ("protein", "codon"))
    if errors and not pending:
        print("[mutEffects-controller] No schedulable tasks remain.", flush=True)
        print_final_diagnostics(warnings, errors)
        raise SystemExit(1)
    try:
        validate_backend_tools(args, manifest, genes)
    except SystemExit:
        print_final_diagnostics(warnings, errors)
        raise

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
    print(f"[resources] {config['hardware']['cpus']} CPUs, {config['hardware']['memory_gib']:.2f} GiB shared RAM, {len(config['hardware']['gpus'])} visible GPUs")
    allocation = config["cpu_allocation"]
    allocation_basis = f"estimated {allocation['concurrent_jobs']} concurrent jobs" if config["threads_explicit"] is False else "explicit --threads"
    print(f"[resources] {config['threads']} threads per task ({allocation_basis}); per-task overrides take precedence")
    for plan in plans:
        print(f"  {plan['gene']}/{plan['side']}: {plan['device']}; CPU RAM {plan['cpu_memory_gib']:.2f} GiB, GPU VRAM {plan['gpu_memory_gib']:.2f} GiB; verified_complete={plan['verified_complete']}")
    for plan in config["evmutation_plans"]:
        print(f"  {plan['gene']}/{plan['side']}: EVmutation; CPU RAM {plan['memory_gib']:.2f} GiB; threads={plan['threads']}; {plan['estimate_basis']}")
    for artifact, gene_list in manifest.items():
        if artifact in {"input_files", "param_files", "routing", "ev_fingerprints", "resource_errors"}:
            continue
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
    try:
        out_dir = args.output.resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        snapshots = out_dir / ".bff-resources"
        args.resource_config_path = write_resource_snapshot(snapshots, json.dumps(config, indent=2) + "\n", ".json")
        args.resource_runtime_config = write_resource_snapshot(snapshots,
            f"executor.cpus = {config['hardware']['cpus']}\n"
            f"executor.memory = '{config['hardware']['memory_gib']:.6f} GB'\n", ".config"
        )
        manifest_path = str(out_dir / ".evmutation_manifest.json")
        with open(manifest_path, 'w') as handle:
            json.dump(manifest, handle, indent=2)

        immutable_manifest = write_resource_snapshot(snapshots, json.dumps(manifest, indent=2), ".manifest.json")
        descriptor, report_path = tempfile.mkstemp(prefix="run-", suffix=".resource-errors.json", dir=snapshots)
        with os.fdopen(descriptor, "w") as handle:
            handle.write("null\n")
        args.resource_errors_path = Path(report_path)
        cmd = build_nextflow_cmd(args, str(immutable_manifest))
        print(f"[mutEffects-controller] Launching: {' '.join(cmd)}", flush=True)
        nf_proc = subprocess.Popen(cmd, cwd=str(HERE))
        exit_code = nf_proc.wait()
    except OSError as error:
        errors.append({"backend": "Nextflow", "message": str(error)})
        print_final_diagnostics(warnings, errors)
        raise SystemExit(1) from error
    try:
        errors.extend(read_resource_errors(args.resource_errors_path))
    except (ValueError, OSError) as error:
        errors.append({"backend": "Nextflow", "message": f"Cannot read resource-error report: {error}"})
    if exit_code:
        errors.append({"backend": "Nextflow", "message": f"Exited with status {exit_code}; see the Nextflow log for details."})
    print_final_diagnostics(warnings, errors)
    sys.exit(exit_code or (1 if errors else 0))


if __name__ == "__main__":
    run_controller(parse_args())
