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
Controller for the EVmutation Nextflow pipeline.

Orchestrates the full MSA generation + plmc + scoring chain:
  1. Protein MSA (jackhmmer -> UniRef90)  — skipped if pre-built
  2. Codon MSA (mmseqs2 -> MAFFT)         — skipped if pre-built
  3. EVmutation scoring (plmc + mutation scoring)

Before launching Nextflow, inventories existing artifacts per gene
(MSAs, model params) and writes a manifest so Nextflow knows what
to generate vs skip for each gene.

Use --resume to continue from a previous Nextflow run.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

from biofeaturefactory.utils.utility import (
    discover_fasta_files,
    extract_gene_from_filename,
)
from biofeaturefactory.EVmutation.evmutation_pipeline import _find_file_for_gene

HERE = Path(__file__).resolve().parent
NEXTFLOW_SCRIPT = HERE / "bin" / "main.nf"


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
    }

    # Final TSVs at {output}/{GENE}/EVmutation/
    tsv_checks = {
        "EVmutation": "protein.tsv",
        "codon_EVmutation": "codon.tsv",
    }

    manifest = {key: [] for key in list(artifact_checks) + list(tsv_checks)}

    for gene in genes:
        for artifact, (path, patterns) in artifact_checks.items():
            if Path(str(path)).exists() and _find_file_for_gene(gene, str(path), patterns):
                manifest[artifact].append(gene)

        for tsv_key, suffix in tsv_checks.items():
            if (out_dir / gene / "EVmutation" / f"{gene}.{suffix}").exists():
                manifest[tsv_key].append(gene)

    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Controller for the EVmutation Nextflow pipeline."
    )

    # Required
    parser.add_argument("--fasta", type=Path, required=True,
                        help="ORF FASTA file or directory of per-gene FASTA files")
    parser.add_argument("--mutations", type=Path, required=True,
                        help="Mutations CSV file or directory of per-gene CSVs")
    parser.add_argument("--plmc-binary", type=str, required=True,
                        help="Path to plmc executable")
    parser.add_argument("--db-root", type=Path, required=True,
                        help="Bio_DBs root directory (contains uniref90.fasta, refseq_assemblies/, etc.)")
    parser.add_argument("--output", "-o", type=Path, default=Path("."),
                        help="Output base directory")

    # Pre-built MSA sources (optional, skips generation for genes that have them)
    parser.add_argument("--msa", type=Path,
                        help="Protein MSA file or directory")
    parser.add_argument("--codon-msa", type=Path,
                        help="Codon MSA file or directory")

    # Tool binaries
    parser.add_argument("--jackhmmer-binary", type=str, default="jackhmmer",
                        help="Path to jackhmmer (default: jackhmmer)")
    parser.add_argument("--jackhmmer-iterations", type=int, default=5,
                        help="Number of jackhmmer iterations (default: 5)")
    parser.add_argument("--mmseqs-binary", type=str, default="mmseqs",
                        help="Path to mmseqs2 (default: mmseqs)")
    parser.add_argument("--aligner", type=str, default="mafft",
                        choices=["mafft", "muscle"],
                        help="Protein aligner (default: mafft)")

    # Pre-built model params
    parser.add_argument("--model-params", type=Path,
                        help="Protein model params file or directory")
    parser.add_argument("--codon-model-params", type=Path,
                        help="Codon model params file or directory")

    # Options
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--validation-log", type=Path)
    parser.add_argument("--pipeline", type=Path, default=NEXTFLOW_SCRIPT,
                        help="Nextflow script (default: bin/main.nf)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume previous Nextflow run")

    args = parser.parse_args()
    validate_args(args)
    return args


def validate_args(args: argparse.Namespace) -> None:
    genes = list(discover_fasta_files(str(args.fasta)).keys())
    if not genes:
        raise SystemExit(f"ERROR: No FASTA files found in {args.fasta}")

    uniref90 = args.db_root / "uniref90.fasta"
    uniref90_gz = args.db_root / "uniref90.fasta.gz"
    if not uniref90.exists() and not uniref90_gz.exists():
        raise SystemExit(f"ERROR: uniref90.fasta(.gz) not found in --db-root ({args.db_root})")


def normalize(path: Optional[Path]) -> Optional[str]:
    return str(path.resolve()) if path else None


def build_nextflow_cmd(args: argparse.Namespace, manifest_path: str) -> List[str]:
    cmd = ["nextflow", "run", str(args.pipeline)]
    if args.resume:
        cmd.append("-resume")

    def add_param(name: str, value):
        if value is not None:
            cmd.extend([f"--{name}", str(value)])

    add_param("fasta", normalize(args.fasta))
    add_param("mutations", normalize(args.mutations))
    add_param("plmc_binary", normalize(Path(args.plmc_binary)))
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

    if args.validation_log:
        add_param("validation_log", normalize(args.validation_log))

    return cmd


def run_controller(args: argparse.Namespace):
    genes = list(discover_fasta_files(str(args.fasta)).keys())
    manifest = build_manifest(genes, args)

    # Summary
    print(f"[evmutation-controller] {len(genes)} gene(s)", flush=True)
    for artifact, gene_list in manifest.items():
        if gene_list:
            preview = ", ".join(gene_list[:5])
            print(f"  {artifact}: {len(gene_list)} pre-built ({preview})")

    ready = [g for g in genes if g in manifest["model_params"] and g in manifest["codon_model_params"]]
    if ready:
        print(f"  ready to score: {len(ready)} gene(s)")

    # Write manifest
    out_dir = args.output.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = str(out_dir / ".evmutation_manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    cmd = build_nextflow_cmd(args, manifest_path)
    print(f"[evmutation-controller] Launching: {' '.join(cmd)}", flush=True)
    nf_proc = subprocess.Popen(cmd, cwd=str(HERE))
    sys.exit(nf_proc.wait())


if __name__ == "__main__":
    run_controller(parse_args())