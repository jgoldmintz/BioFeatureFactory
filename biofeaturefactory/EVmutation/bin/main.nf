#!/usr/bin/env nextflow
// BioFeatureFactory
// Copyright (C) 2023-2026  Jacob Goldmintz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// BioFeatureFactory — EVmutation full pipeline (MSA generation + plmc + scoring)
// DAG: ORF FASTA -> [protein MSA || codon MSA] -> EVmutation scoring
nextflow.enable.dsl = 2

// ---------------- PARAMETERS ----------------
// Mirrors evmutation_pipeline.py argument names
params.fasta              = null    // file or directory
params.mutations          = null    // file or directory
params.plmc_binary        = null
params.output_dir         = '.'
params.validation_log     = null
params.threads            = 4

// Protein MSA: pre-built file/dir (skips jackhmmer) or generate from UniRef90
params.msa                = null    // file or directory
params.uniref90_db        = null
params.jackhmmer_binary   = 'jackhmmer'
params.jackhmmer_iterations = 5

// Codon MSA: pre-built file/dir (skips mmseqs2/MAFFT) or generate from Bio_DBs
params.codon_msa          = null    // file or directory
params.db_root            = null
params.mmseqs_binary      = 'mmseqs'
params.aligner            = 'mafft'
params.manifest           = null    // JSON manifest from controller

// Required checks
if (!params.fasta)        error "ERROR: --fasta is required"
if (!params.mutations)    error "ERROR: --mutations is required"
if (!params.plmc_binary)  error "ERROR: --plmc_binary is required"

if (!params.msa && !params.uniref90_db)
    error "ERROR: --uniref90_db is required when not providing --msa"
if (!params.codon_msa && !params.db_root)
    error "ERROR: --db_root is required when not providing --codon_msa"

// ---------------- HELPERS ----------------
def resolveInputPath(pathParam) {
    def f = new File(pathParam as String)
    return f.isDirectory() ? f : (f.isFile() ? f : null)
}

def resolveMutationCsv(gene_id) {
    def mut_base = new File(params.mutations as String)
    if (mut_base.isFile()) return mut_base.getAbsolutePath()

    def patterns = ["${gene_id}_mutations.csv", "combined_${gene_id}.csv", "${gene_id}.csv"]
    for (pat in patterns) {
        def f = new File(mut_base, pat)
        if (f.exists()) return f.getAbsolutePath()
    }
    def files = (mut_base.listFiles() ?: []) as List
    def match = files.find { it.name.toLowerCase().contains(gene_id.toLowerCase()) && it.name.endsWith('.csv') }
    return match ? match.getAbsolutePath() : null
}

def resolveMsaFile(msa_param, gene_id, extensions) {
    // msa_param is file or directory; resolve to per-gene MSA file
    def base = new File(msa_param as String)
    if (base.isFile()) return base.getAbsolutePath()

    // Directory: search for gene-matching file
    def files = (base.listFiles() ?: []) as List
    for (ext in extensions) {
        def match = files.find { it.name.toLowerCase().startsWith(gene_id.toLowerCase()) && it.name.endsWith(ext) }
        if (match) return match.getAbsolutePath()
    }
    return null
}

// ---------------- WORKFLOW ----------------
workflow {
    // Load manifest (written by controller)
    def manifest = [:]
    if (params.manifest) {
        manifest = new groovy.json.JsonSlurper().parse(new File(params.manifest as String))
    }

    // Build per-gene FASTA channel
    def fasta_base = new File(params.fasta as String)
    def fasta_ch

    if (fasta_base.isDirectory()) {
        fasta_ch = Channel
            .fromPath("${params.fasta}/*.{fasta,fa,fas}")
            .map { f ->
                def gene = f.baseName.replaceAll(/_(nt|aa)$/, '')
                tuple(gene, f)
            }
    } else {
        fasta_ch = Channel
            .of(tuple(fasta_base.name.replaceAll(/\.(fasta|fa|fas)$/, '').replaceAll(/_(nt|aa)$/, ''), file(params.fasta)))
    }

    // Step 1: Protein MSA — split genes by manifest: pre-built vs needs generation
    def has_protein_msa = (manifest.msa ?: []) as Set
    def fasta_need_protein = fasta_ch.filter { gene_id, f -> !(gene_id in has_protein_msa) }
    def fasta_have_protein = fasta_ch.filter { gene_id, f -> gene_id in has_protein_msa }

    def generated_protein = generate_protein_msa(fasta_need_protein)
        .map { items -> tuple(items[0], items[1]) }
    def prebuilt_protein = fasta_have_protein.map { gene_id, fasta_file ->
        def msa_path = resolveMsaFile(params.msa, gene_id, ['.a2m', '.msa.a2m', '.fasta'])
        tuple(gene_id, file(msa_path))
    }
    def protein_msa_ch = generated_protein.mix(prebuilt_protein)

    // Step 2: Codon MSA — same split
    def has_codon_msa = (manifest.codon_msa ?: []) as Set
    def fasta_need_codon = fasta_ch.filter { gene_id, f -> !(gene_id in has_codon_msa) }
    def fasta_have_codon = fasta_ch.filter { gene_id, f -> gene_id in has_codon_msa }

    def generated_codon = generate_codon_msa(fasta_need_codon)
        .map { items -> tuple(items[0], items[1]) }
    def prebuilt_codon = fasta_have_codon.map { gene_id, fasta_file ->
        def msa_path = resolveMsaFile(params.codon_msa, gene_id, ['.codon.msa.fasta', '.codon.fasta', '.fasta'])
        tuple(gene_id, file(msa_path))
    }
    def codon_msa_ch = generated_codon.mix(prebuilt_codon)

    // Step 3: Protein EVmutation — runs as soon as protein MSA is ready
    def has_protein_tsv = (manifest.EVmutation ?: []) as Set
    def protein_ev_input = fasta_ch
        .filter { gene_id, f -> !(gene_id in has_protein_tsv) }
        .join(protein_msa_ch)
        .map { gene_id, fasta_file, protein_msa ->
            def mut_csv = resolveMutationCsv(gene_id)
            tuple(gene_id, fasta_file, protein_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'))
        }
        .filter { gene_id, fasta_file, protein_msa, mut_csv ->
            if (mut_csv.name == 'NO_MUTATIONS') {
                println "WARNING: No mutation CSV found for ${gene_id}, skipping protein EVmutation"
                return false
            }
            return true
        }
    run_protein_evmutation(protein_ev_input)

    // Step 4: Codon EVmutation — runs as soon as codon MSA is ready (independent of protein)
    def has_codon_tsv = (manifest.codon_EVmutation ?: []) as Set
    def codon_ev_input = fasta_ch
        .filter { gene_id, f -> !(gene_id in has_codon_tsv) }
        .join(codon_msa_ch)
        .map { gene_id, fasta_file, codon_msa ->
            def mut_csv = resolveMutationCsv(gene_id)
            tuple(gene_id, fasta_file, codon_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'))
        }
        .filter { gene_id, fasta_file, codon_msa, mut_csv ->
            if (mut_csv.name == 'NO_MUTATIONS') {
                println "WARNING: No mutation CSV found for ${gene_id}, skipping codon EVmutation"
                return false
            }
            return true
        }
    run_codon_evmutation(codon_ev_input)
}

// ---------------- PROCESSES ----------------
process generate_protein_msa {
    publishDir "${params.output_dir}/MSA", mode: 'copy'
    tag { gene_id }

    input:
    tuple val(gene_id), path(fasta_file)

    output:
    tuple val(gene_id), path("${gene_id}.msa.a2m"), path("${gene_id}.msa.stats.json")

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/../../utils/msa_generation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --database "${params.uniref90_db}" \\
        --jackhmmer-binary "${params.jackhmmer_binary}" \\
        --output . \\
        --threads ${params.threads} \\
        --iterations ${params.jackhmmer_iterations}

    # Pipeline writes to {GENE}/MSA/ — flatten
    if [ -d "${gene_id}/MSA" ]; then
        mv ${gene_id}/MSA/${gene_id}.msa.a2m . 2>/dev/null || true
        mv ${gene_id}/MSA/${gene_id}.msa.stats.json . 2>/dev/null || true
    fi
    """
}

process generate_codon_msa {
    publishDir "${params.output_dir}/CodonMSA", mode: 'copy'
    tag { gene_id }

    input:
    tuple val(gene_id), path(fasta_file)

    output:
    tuple val(gene_id), path("${gene_id}.codon.msa.fasta"), path("${gene_id}.codon.msa.manifest.tsv"), path("${gene_id}.codon.msa.stats.json")

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/../../utils/codon_msa_pipeline.py \\
        --fasta "${fasta_file}" \\
        --db-root "${params.db_root}" \\
        --output . \\
        --mmseqs-binary "${params.mmseqs_binary}" \\
        --aligner ${params.aligner} \\
        --threads ${params.threads}

    # Pipeline writes to {GENE}/CodonMSA/ — flatten
    if [ -d "${gene_id}/CodonMSA" ]; then
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.fasta . 2>/dev/null || true
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.manifest.tsv . 2>/dev/null || true
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.stats.json . 2>/dev/null || true
    fi
    """
}

process run_protein_evmutation {
    publishDir "${params.output_dir}/${gene_id}/EVmutation", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/model_params", mode: 'copy', pattern: '*.model_params'
    tag { "${gene_id} protein" }

    input:
    tuple val(gene_id), path(fasta_file), path(protein_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.protein.tsv"),
          path("${gene_id}.model_params"), optional: true

    script:
    def logArg = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    """
    set -euo pipefail
    python3 ${projectDir}/../evmutation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --msa "${protein_msa}" \\
        --plmc-binary "${params.plmc_binary}" \\
        --output . \\
        ${logArg}

    if [ -d "${gene_id}/EVmutation" ]; then
        mv ${gene_id}/EVmutation/${gene_id}.protein.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.model_params" -exec mv {} . \\; 2>/dev/null || true
    """
}

process run_codon_evmutation {
    publishDir "${params.output_dir}/${gene_id}/EVmutation", mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/codon_model_params", mode: 'copy', pattern: '*.codon_model_params'
    tag { "${gene_id} codon" }

    input:
    tuple val(gene_id), path(fasta_file), path(codon_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.codon.tsv"),
          path("${gene_id}.codon_model_params"), optional: true

    script:
    def logArg = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    """
    set -euo pipefail
    python3 ${projectDir}/../evmutation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --codon-msa "${codon_msa}" \\
        --plmc-binary "${params.plmc_binary}" \\
        --output . \\
        ${logArg}

    if [ -d "${gene_id}/EVmutation" ]; then
        mv ${gene_id}/EVmutation/${gene_id}.codon.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.codon_model_params" -exec mv {} . \\; 2>/dev/null || true
    """
}
