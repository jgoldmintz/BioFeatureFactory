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

// BioFeatureFactory -- EVmutation full pipeline (MSA generation + plmc + scoring)
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

// Pre-built plmc params: when provided, evmutation_pipeline.py skips plmc and
// scores the mutations directly. File or directory of {GENE}.(codon_)model_params.
params.model_params       = null
params.codon_model_params = null

params.manifest           = null    // JSON manifest from controller

// Backend gating -- defaults match the controller (both backends run).
// Skip flags are emitted as the string "true" by the controller only when
// disabling; never as "false" (Groovy's `if ("false")` is truthy).
params.skip_evmutation          = false   // when true: EVmutation/plmc processes skipped
params.skip_adabmdca            = false   // when true: adabmDCA processes skipped
params.skip_codon_evmutation    = false   // when true: codon-side EVmutation skipped
params.skip_codon_adabmdca      = false   // when true: codon-side adabmDCA skipped

// adabmDCA backend tunables (consulted when adabmDCA runs).
// The `adabmDCA` CLI is the Python/torch console script installed by
// `pip install adabmDCA` -- it is expected on PATH and is not configurable here.
params.adabmdca_protein_params     = null          // pre-built protein params: file or directory
params.adabmdca_codon_params       = null          // pre-built codon params:   file or directory
// Alphabet selection is owned by adabmdca_pipeline.py: protein side uses the
// 21-char default ("-ACDEFGHIKLMNPQRSTVWY"); codon side uses the 65-char
// alphabet from bin/codon_encoding.py after encoding.
params.adabmdca_model              = 'pseudoDCA'
// null => adabmdca_pipeline.py picks per backend (500 pseudoDCA / 50000 Boltzmann).
// A fixed 50000 here forced the pseudoDCA path to run 100x its own default.
params.adabmdca_nepochs            = null
params.adabmdca_tol                = 0.001         // pseudoDCA early-stop threshold
params.adabmdca_patience           = 3
params.adabmdca_check_every        = 10
params.adabmdca_target             = 0.95          // Pearson target on Cij
params.adabmdca_lr                 = 0.01
params.adabmdca_nchains            = 10000         // PCD chain count
params.adabmdca_nsweeps            = 10            // sweeps per gradient step
params.adabmdca_device             = 'cuda'
params.adabmdca_dtype              = 'float32'
params.adabmdca_seed               = 0
params.resource_config            = null
params.resource_errors            = null
params.resource_executor          = 'local'
params.gpu_slots                  = 1
params.msa_cpus                   = 4
params.msa_memory                 = '8 GB'
params.evmutation_cpus            = 4
params.evmutation_memory          = null

// Required checks live inside the workflow block (Nextflow 26+ forbids
// top-level statements outside process / workflow / function bodies).
// See workflow { validateRequiredParams(); ... }.

def backendSideEnabled(gene_id, side, backend, manifest) {
    def backendSkipped = backend == 'EVmutation' ? params.skip_evmutation : params.skip_adabmdca
    if (backendSkipped) return false
    def routing = manifest.routing?.get(gene_id)
    def protein = routing == null ? true : routing.protein
    def codon = routing == null ? true : routing.codon
    def codonSkipped = backend == 'EVmutation' ? params.skip_codon_evmutation : params.skip_codon_adabmdca
    if (codonSkipped)
        return side == 'protein' && (protein || codon)
    return side == 'protein' ? protein : codon
}

def backendSidePending(gene_id, side, backend, manifest) {
    def artifact = backend == 'EVmutation' ?
        (side == 'protein' ? 'EVmutation' : 'codon_EVmutation') : "adabmdca_${side}"
    return backendSideEnabled(gene_id, side, backend, manifest) &&
        !(gene_id in (manifest[artifact] ?: [])) &&
        !(manifest.resource_errors ?: []).any {
            it.gene == gene_id && it.side == side && it.backend == backend.toLowerCase()
        }
}

def pendingSide(manifest, backend, side) {
    if (manifest.input_files != null)
        return manifest.input_files.keySet().any { backendSidePending(it, side, backend, manifest) }
    return backendSideEnabled('', side, backend, manifest)
}

def msaNeeded(gene_id, side, manifest) {
    return backendSidePending(gene_id, side, 'EVmutation', manifest) ||
        backendSidePending(gene_id, side, 'adabmDCA', manifest)
}

def resolveGeneMsa(gene_id, side, manifest) {
    def artifact = side == 'protein' ? 'msa' : 'codon_msa'
    def resolved = manifest.input_files?.get(gene_id)?.get(artifact)
    if (resolved) return resolved
    def supplied = side == 'protein' ? params.msa : params.codon_msa
    if (!supplied) return null
    def extensions = side == 'protein' ? ['.a2m', '.msa.a2m', '.fasta'] : ['.codon.msa.fasta', '.codon.fasta', '.fasta']
    return resolveMsaFile(supplied, gene_id, extensions)
}

def resolveEvParams(gene_id, side, manifest) {
    def artifact = side == 'protein' ? 'model_params' : 'codon_model_params'
    def resolved = manifest.param_files?.get(gene_id)?.get(artifact)
    if (resolved) return resolved
    def supplied = side == 'protein' ? params.model_params : params.codon_model_params
    if (supplied) return resolveMsaFile(supplied, gene_id, [".${artifact}"])
    def existing = new File("${params.output_dir}/${artifact}/${gene_id}.${artifact}")
    return existing.isFile() ? existing.getAbsolutePath() : null
}

def pendingPlmc(manifest, side) {
    if (manifest.input_files != null)
        return manifest.input_files.keySet().any {
            backendSidePending(it, side, 'EVmutation', manifest) && !resolveEvParams(it, side, manifest)
        }
    def supplied = side == 'protein' ? params.model_params : params.codon_model_params
    return pendingSide(manifest, 'EVmutation', side) && !supplied
}

def scoreOptions(gene_id, side, backend, manifest) {
    return [
        skip_codon: !backendSideEnabled(gene_id, 'codon', backend, manifest),
        score_missense_codon: manifest.routing?.get(gene_id)?.score_missense_codon ?: false,
        model_params: resolveEvParams(gene_id, side, manifest),
        fingerprint: manifest.ev_fingerprints?.get(gene_id)?.get(side),
    ]
}

def validateRequiredParams(manifest) {
    if (!params.fasta)        error "ERROR: --fasta is required"
    if (!params.mutations)    error "ERROR: --mutations is required"
    if (!params.resource_errors && manifest.resource_errors)
        error "ERROR: Manifest contains resource planning errors; run through mutEffects_controller.py to defer them"
    // plmc is only needed when params have to be built. Providing --model_params
    // (protein) / --codon_model_params (codon) skips plmc for that side, so the
    // binary is required only when a side still has to run inference.
    def needProteinPlmc = pendingPlmc(manifest, 'protein')
    def needCodonPlmc = pendingPlmc(manifest, 'codon')
    if ((needProteinPlmc || needCodonPlmc) && !params.plmc_binary)
        error "ERROR: --plmc_binary is required when the EVmutation backend builds params " +
              "(provide --model_params / --codon_model_params to skip plmc, or --skip_evmutation true)"
    def inputs = manifest.input_files
    def proteinMissing = inputs != null ? inputs.keySet().any {
        msaNeeded(it, 'protein', manifest) && !resolveGeneMsa(it, 'protein', manifest)
    } : (pendingSide(manifest, 'EVmutation', 'protein') || pendingSide(manifest, 'adabmDCA', 'protein')) && !params.msa
    def codonMissing = inputs != null ? inputs.keySet().any {
        msaNeeded(it, 'codon', manifest) && !resolveGeneMsa(it, 'codon', manifest)
    } : (pendingSide(manifest, 'EVmutation', 'codon') || pendingSide(manifest, 'adabmDCA', 'codon')) && !params.codon_msa
    if (proteinMissing && !params.uniref90_db)
        error "ERROR: --uniref90_db is required when not providing --msa"
    if (codonMissing && !params.db_root)
        error "ERROR: --db_root is required when not providing --codon_msa"
    if ((pendingSide(manifest, 'adabmDCA', 'protein') || pendingSide(manifest, 'adabmDCA', 'codon')) && !params.resource_config)
        error "ERROR: --resource_config is required for adabmDCA resource planning"
    if (params.resource_executor != 'local')
        error "ERROR: --resource_executor must be local"
    if ((params.gpu_slots as int) < 0)
        error "ERROR: --gpu_slots cannot be negative"
}

// ---------------- HELPERS ----------------
def resolveInputPath(pathParam) {
    def f = new File(pathParam as String)
    return f.isDirectory() ? f : (f.isFile() ? f : null)
}

def resolveMutationCsv(gene_id, manifest) {
    def resolved = manifest.input_files?.get(gene_id)?.mutations
    if (resolved) return resolved
    def mut_base = new File(params.mutations as String)
    if (mut_base.isFile()) return mut_base.getAbsolutePath()

    // Nextflow 26+ removed `for` loops in DSL scripts; use findResult for short-circuit search.
    def patterns = ["${gene_id}_mutations.csv", "combined_${gene_id}.csv", "${gene_id}.csv"]
    def matched = patterns.findResult { pat ->
        def f = new File(mut_base, pat)
        f.exists() ? f.getAbsolutePath() : null
    }
    if (matched) return matched

    // Fuzzy fallback, but anchored. `contains` matched SMN2.csv for gene SMN
    // (verified under nextflow), so a gene whose name is a prefix of another
    // silently scored against the wrong gene's mutation list. Require the stem
    // to BE the gene, or to start with the gene followed by a separator, so
    // SMN2_variants.csv still resolves for SMN2 and never for SMN.
    def files = (mut_base.listFiles() ?: []) as List
    def g = gene_id.toLowerCase()
    def fuzzy = files.find { f ->
        def n = f.name.toLowerCase()
        if (!n.endsWith('.csv')) return false
        def stem = n.substring(0, n.length() - 4)
        return stem == g || stem.startsWith(g + '.') || stem.startsWith(g + '_') || stem.startsWith(g + '-')
    }
    return fuzzy ? fuzzy.getAbsolutePath() : null
}

def resolveMsaFile(msa_param, gene_id, extensions) {
    // msa_param is file or directory; resolve to per-gene MSA file
    def base = new File(msa_param as String)
    if (base.isFile()) return base.getAbsolutePath()

    // Directory: search for gene-matching file. findResult returns the first non-null match.
    // Anchored on a '.' boundary for the same reason as resolveMutationCsv:
    // `startsWith` matched SMN2.msa.a2m for gene SMN, so SMN would have been
    // scored against SMN2's Potts model with no warning. Files are named
    // {GENE}.msa.a2m / {GENE}.codon.msa.fasta, so requiring "{gene}." is exact.
    def files = (base.listFiles() ?: []) as List
    def g = gene_id.toLowerCase()
    return extensions.findResult { ext ->
        def m = files.find { it.name.toLowerCase().startsWith(g + '.') && it.name.endsWith(ext) }
        m ? m.getAbsolutePath() : null
    }
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def evCacheCommand(gene_id, side, msa_file, score_options) {
    if (!score_options.fingerprint) return ':'
    def model_suffix = side == 'protein' ? 'model_params' : 'codon_model_params'
    def model_path = score_options.model_params ?: "${gene_id}.${model_suffix}"
    def arguments = [
        'python3', "${projectDir}/evmutation_cache.py", 'write',
        '--marker', "${gene_id}.${side}.routing.json",
        '--fingerprint', score_options.fingerprint,
        '--tsv', "${gene_id}.${side}.tsv", '--params', model_path, '--msa', msa_file,
    ]
    return arguments.collect { shellQuote(it) }.join(' ')
}

def adabmdcaCommand(side, fasta_file, msa_file, mutations_csv, plan_file, plan, device, threads) {
    def settings = plan.settings
    def arguments = [
        'python3', "${projectDir}/adabmdca_task.py",
        '--plan', plan_file, '--device', device, '--threads', threads, '--',
        '--fasta', fasta_file, '--mutations', mutations_csv,
        side == 'protein' ? '--msa' : '--codon-msa', msa_file,
        '--output', '.', '--adabmdca-model', settings.model,
        '--adabmdca-tol', settings.tol,
        '--adabmdca-patience', settings.patience,
        '--adabmdca-check-every', settings.check_every,
        '--adabmdca-target', settings.target,
        '--adabmdca-lr', settings.lr,
        '--adabmdca-nchains', settings.nchains,
        '--adabmdca-nsweeps', settings.nsweeps,
        '--adabmdca-dtype', settings.dtype,
        '--adabmdca-seed', settings.seed,
    ]
    if (settings.nepochs != null)
        arguments.addAll(['--adabmdca-nepochs', settings.nepochs])
    if (side == 'protein' && settings.skip_codon)
        arguments.add('--skip-codon')
    if (side == 'codon' && settings.score_missense_codon)
        arguments.add('--score-missense-codon')
    if (params.validation_log)
        arguments.addAll(['--validation-log', file(params.validation_log)])
    return arguments.collect { shellQuote(it) }.join(' ')
}

// ---------------- WORKFLOW ----------------
workflow {
    // Load manifest (written by controller)
    def manifest = [:]
    if (params.manifest) {
        manifest = new groovy.json.JsonSlurperClassic().parse(new File(params.manifest as String))
    }
    def resource_errors = new java.util.concurrent.CopyOnWriteArrayList()
    def resource_error_report = params.resource_errors
    workflow.onComplete {
        def collected = resource_errors.toList().sort { entry -> "${entry.gene}/${entry.side}/${entry.backend}" }
        if (resource_error_report) {
            def report = new File(resource_error_report as String)
            report.parentFile?.mkdirs()
            report.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(collected)) + '\n'
        }
    }
    validateRequiredParams(manifest)

    // Build per-gene FASTA channel
    def fasta_base = new File(params.fasta as String)
    def fasta_ch

    if (manifest.input_files != null) {
        fasta_ch = Channel.fromList(manifest.input_files.collect { gene_id, inputs ->
            tuple(gene_id, file(inputs.fasta, checkIfExists: true))
        })
    } else if (fasta_base.isDirectory()) {
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

    // Step 1: Protein MSA -- split genes by manifest: pre-built vs needs generation
    def protein_fasta = fasta_ch.filter { gene_id, fasta_file -> msaNeeded(gene_id, 'protein', manifest) }
    def fasta_need_protein = protein_fasta.filter { gene_id, fasta_file -> !resolveGeneMsa(gene_id, 'protein', manifest) }
    def fasta_have_protein = protein_fasta.filter { gene_id, fasta_file -> resolveGeneMsa(gene_id, 'protein', manifest) }

    def generated_protein = generate_protein_msa(fasta_need_protein)
        .map { items -> tuple(items[0], items[1]) }
    def prebuilt_protein = fasta_have_protein.map { gene_id, fasta_file ->
        tuple(gene_id, file(resolveGeneMsa(gene_id, 'protein', manifest), checkIfExists: true))
    }
    def protein_msa_ch = generated_protein.mix(prebuilt_protein)

    // Step 2: Codon MSA -- same split
    def codon_fasta = fasta_ch.filter { gene_id, fasta_file -> msaNeeded(gene_id, 'codon', manifest) }
    def fasta_need_codon = codon_fasta.filter { gene_id, fasta_file -> !resolveGeneMsa(gene_id, 'codon', manifest) }
    def fasta_have_codon = codon_fasta.filter { gene_id, fasta_file -> resolveGeneMsa(gene_id, 'codon', manifest) }

    def generated_codon = generate_codon_msa(fasta_need_codon)
        .map { items -> tuple(items[0], items[1]) }
    def prebuilt_codon = fasta_have_codon.map { gene_id, fasta_file ->
        tuple(gene_id, file(resolveGeneMsa(gene_id, 'codon', manifest), checkIfExists: true))
    }
    def codon_msa_ch = generated_codon.mix(prebuilt_codon)

    // Step 3-4: EVmutation backend -- gated on !skip_evmutation.
    if (pendingSide(manifest, 'EVmutation', 'protein') || pendingSide(manifest, 'EVmutation', 'codon')) {
        // Step 3: Protein EVmutation -- runs as soon as protein MSA is ready
        def protein_ev_input = fasta_ch
            .filter { gene_id, fasta_file -> backendSidePending(gene_id, 'protein', 'EVmutation', manifest) }
            .join(protein_msa_ch)
            .map { gene_id, fasta_file, protein_msa ->
                def mut_csv = resolveMutationCsv(gene_id, manifest)
                tuple(gene_id, 'protein', fasta_file, protein_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'), scoreOptions(gene_id, 'protein', 'EVmutation', manifest))
            }
            .filter { gene_id, side, fasta_file, protein_msa, mut_csv, score_options ->
                if (mut_csv.name == 'NO_MUTATIONS') {
                    println "WARNING: No mutation CSV found for ${gene_id}, skipping protein EVmutation"
                    return false
                }
                return true
            }
        def codon_ev_input = fasta_ch
            .filter { gene_id, fasta_file -> backendSidePending(gene_id, 'codon', 'EVmutation', manifest) }
            .join(codon_msa_ch)
            .map { gene_id, fasta_file, codon_msa ->
                def mut_csv = resolveMutationCsv(gene_id, manifest)
                tuple(gene_id, 'codon', fasta_file, codon_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'), scoreOptions(gene_id, 'codon', 'EVmutation', manifest))
            }
            .filter { gene_id, side, fasta_file, codon_msa, mut_csv, score_options ->
                if (mut_csv.name == 'NO_MUTATIONS') {
                    println "WARNING: No mutation CSV found for ${gene_id}, skipping codon EVmutation"
                    return false
                }
                return true
            }
        def ev_config = params.resource_config ? file(params.resource_config, checkIfExists: true) : []
        def ev_planned = plan_evmutation(protein_ev_input.mix(codon_ev_input), ev_config)
            .map { gene_id, side, fasta_file, msa_file, mutations_csv, score_options, plan_file ->
                def plan = new groovy.json.JsonSlurperClassic().parse(plan_file.toFile())
                tuple(gene_id, side, fasta_file, msa_file, mutations_csv, score_options, plan_file, plan)
            }
            .filter { gene_id, side, fasta_file, msa_file, mutations_csv, score_options, plan_file, plan ->
                if (plan.resource_error) {
                    resource_errors.add(plan.resource_error)
                    return false
                }
                return true
            }
            .branch { gene_id, side, fasta_file, msa_file, mutations_csv, score_options, plan_file, plan ->
                protein: side == 'protein'
                codon: side == 'codon'
            }
        run_protein_evmutation(ev_planned.protein)
        run_codon_evmutation(ev_planned.codon)
    }

    if (pendingSide(manifest, 'adabmDCA', 'protein') || pendingSide(manifest, 'adabmDCA', 'codon')) {
        def protein_adabm_input = fasta_ch
            .filter { gene_id, fasta_file -> backendSidePending(gene_id, 'protein', 'adabmDCA', manifest) }
            .join(protein_msa_ch)
            .map { gene_id, fasta_file, msa_file ->
                tuple(gene_id, 'protein', fasta_file, msa_file)
            }
        def codon_adabm_input = fasta_ch
            .filter { gene_id, fasta_file ->
                backendSidePending(gene_id, 'codon', 'adabmDCA', manifest)
            }
            .join(codon_msa_ch)
            .map { gene_id, fasta_file, msa_file ->
                tuple(gene_id, 'codon', fasta_file, msa_file)
            }
        def adabm_input = protein_adabm_input.mix(codon_adabm_input)
            .map { gene_id, side, fasta_file, msa_file ->
                def mut_csv = resolveMutationCsv(gene_id, manifest)
                tuple(gene_id, side, fasta_file, msa_file, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'))
            }
            .filter { gene_id, side, fasta_file, msa_file, mut_csv ->
                if (mut_csv.name == 'NO_MUTATIONS') {
                    println "WARNING: No mutation CSV found for ${gene_id}, skipping ${side} adabmDCA"
                    return false
                }
                return true
            }
        def planned = plan_adabmdca(adabm_input, file(params.resource_config, checkIfExists: true))
            .map { gene_id, side, fasta_file, msa_file, mutations_csv, plan_file ->
                def plan = new groovy.json.JsonSlurperClassic().parse(plan_file.toFile())
                tuple(gene_id, side, fasta_file, msa_file, mutations_csv, plan_file, plan)
            }
            .filter { gene_id, side, fasta_file, msa_file, mutations_csv, plan_file, plan ->
                if (plan.resource_error) {
                    resource_errors.add(plan.resource_error)
                    return false
                }
                if (!(plan.device in ['cpu', 'cuda']))
                    error "ERROR: Invalid planned device for ${gene_id} ${side}: ${plan.device}"
                if (plan.device == 'cuda' && (params.gpu_slots as int) < 1)
                    error "ERROR: GPU plans require at least one allocated GPU slot"
                return true
            }
            .branch { gene_id, side, fasta_file, msa_file, mutations_csv, plan_file, plan ->
                gpu: plan.device == 'cuda'
                cpu: plan.device == 'cpu'
            }
        run_adabmdca_gpu(planned.gpu)
        def cpu_retry = run_adabmdca_gpu.out.gpu_oom
            .map { gene_id, side, fasta_file, msa_file, mutations_csv, plan_file, plan, oom_file ->
                if (!plan.cpu_fallback_allowed)
                    error "ERROR: CPU fallback is unavailable for ${gene_id} ${side}"
                tuple(gene_id, side, fasta_file, msa_file, mutations_csv, plan_file, plan)
            }
        run_adabmdca_cpu(planned.cpu.mix(cpu_retry))
    }
}

// ---------------- PROCESSES ----------------
process generate_protein_msa {
    executor params.resource_executor
    cpus { params.msa_cpus as int }
    memory { params.msa_memory }
    publishDir "${params.output_dir}/MSA", mode: 'copy'
    tag { gene_id }

    input:
    tuple val(gene_id), path(fasta_file)

    output:
    tuple val(gene_id), path("${gene_id}.msa.a2m"), path("${gene_id}.msa.stats.json")

    script:
    """
    set -euo pipefail
    python3 -m biofeaturefactory.core.msa_generation_pipeline \\
        --fasta "${fasta_file}" \\
        --database "${params.uniref90_db}" \\
        --jackhmmer-binary "${params.jackhmmer_binary}" \\
        --output . \\
        --threads ${task.cpus} \\
        --iterations ${params.jackhmmer_iterations}

    # Pipeline writes to {GENE}/MSA/ -- flatten
    if [ -d "${gene_id}/MSA" ]; then
        mv ${gene_id}/MSA/${gene_id}.msa.a2m . 2>/dev/null || true
        mv ${gene_id}/MSA/${gene_id}.msa.stats.json . 2>/dev/null || true
    fi
    """
}

process generate_codon_msa {
    executor params.resource_executor
    cpus { params.msa_cpus as int }
    memory { params.msa_memory }
    publishDir "${params.output_dir}/CodonMSA", mode: 'copy'
    tag { gene_id }

    input:
    tuple val(gene_id), path(fasta_file)

    output:
    tuple val(gene_id), path("${gene_id}.codon.msa.fasta"), path("${gene_id}.codon.msa.manifest.tsv"), path("${gene_id}.codon.msa.stats.json")

    script:
    """
    set -euo pipefail
    python3 -m biofeaturefactory.core.codon_msa_pipeline \\
        --fasta "${fasta_file}" \\
        --db-root "${params.db_root}" \\
        --output . \\
        --mmseqs-binary "${params.mmseqs_binary}" \\
        --aligner ${params.aligner} \\
        --threads ${task.cpus}

    # Pipeline writes to {GENE}/CodonMSA/ -- flatten
    if [ -d "${gene_id}/CodonMSA" ]; then
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.fasta . 2>/dev/null || true
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.manifest.tsv . 2>/dev/null || true
        mv ${gene_id}/CodonMSA/${gene_id}.codon.msa.stats.json . 2>/dev/null || true
    fi
    """
}

process plan_evmutation {
    executor params.resource_executor
    cpus 1
    memory '512 MB'
    tag { "${gene_id} ${side} EVmutation plan" }
    publishDir "${params.output_dir}/resource_plans", mode: 'copy', pattern: '*.evmutation.plan.json'

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv), val(score_options)
    path resource_config

    output:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv), val(score_options),
          path("${gene_id}.${side}.evmutation.plan.json")

    script:
    def arguments = [
        'python3', "${projectDir}/plmc_resources.py", 'plan',
        '--gene', gene_id, '--side', side, '--msa', msa_file,
        '--threads', params.evmutation_cpus, '--output', "${gene_id}.${side}.evmutation.plan.json",
    ]
    if (params.resource_errors)
        arguments.add('--defer-errors')
    if (resource_config)
        arguments.addAll(['--config', resource_config])
    if (params.evmutation_memory)
        arguments.addAll(['--memory-gib', nextflow.util.MemoryUnit.of(params.evmutation_memory).toBytes() / 1073741824])
    if (score_options.model_params)
        arguments.addAll(['--params', score_options.model_params])
    def command = arguments.collect { shellQuote(it) }.join(' ')
    """
    set -euo pipefail
    ${command}
    """
}

process run_protein_evmutation {
    executor params.resource_executor
    cpus { plan.threads as int }
    memory { "${plan.memory_gib} GB" }
    publishDir { "${params.output_dir}/${gene_id}/EVmutation" }, mode: 'copy', pattern: '*.{tsv,routing.json}'
    publishDir "${params.output_dir}/model_params", mode: 'copy', pattern: '*.model_params'
    tag { "${gene_id} protein" }

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(protein_msa), path(mutations_csv),
          val(score_options), path(plan_file), val(plan)

    output:
    tuple val(gene_id), path("${gene_id}.protein.tsv"), emit: scores
    tuple val(gene_id), path("${gene_id}.model_params"), emit: model_params, optional: true
    tuple val(gene_id), path("${gene_id}.protein.routing.json"), emit: cache_marker, optional: true

    script:
    def logArg     = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def skipArg    = score_options.skip_codon ? "--skip-codon" : ""
    def plmcArg    = params.plmc_binary ? "--plmc-binary \"${params.plmc_binary}\"" : ""
    def mparamsArg = score_options.model_params ? "--model-params ${shellQuote(score_options.model_params)}" : ""
    def cacheCommand = evCacheCommand(gene_id, 'protein', protein_msa, score_options)
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    python3 ${projectDir}/../evmutation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --msa "${protein_msa}" \\
        ${mparamsArg} \\
        ${plmcArg} \\
        --output . \\
        ${skipArg} \\
        ${logArg}

    if [ -d "${gene_id}/EVmutation" ]; then
        mv ${gene_id}/EVmutation/${gene_id}.protein.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.model_params" -exec mv {} . \\; 2>/dev/null || true
    ${cacheCommand}
    """
}

process run_codon_evmutation {
    executor params.resource_executor
    cpus { plan.threads as int }
    memory { "${plan.memory_gib} GB" }
    publishDir { "${params.output_dir}/${gene_id}/EVmutation" }, mode: 'copy', pattern: '*.{tsv,routing.json}'
    publishDir "${params.output_dir}/codon_model_params", mode: 'copy', pattern: '*.codon_model_params'
    tag { "${gene_id} codon" }

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(codon_msa), path(mutations_csv),
          val(score_options), path(plan_file), val(plan)

    output:
    tuple val(gene_id), path("${gene_id}.codon.tsv"), emit: scores
    tuple val(gene_id), path("${gene_id}.codon_model_params"), emit: model_params, optional: true
    tuple val(gene_id), path("${gene_id}.codon.routing.json"), emit: cache_marker, optional: true

    script:
    def logArg      = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def plmcArg     = params.plmc_binary ? "--plmc-binary \"${params.plmc_binary}\"" : ""
    def cmparamsArg = score_options.model_params ? "--codon-model-params ${shellQuote(score_options.model_params)}" : ""
    def missenseArg = score_options.score_missense_codon ? "--score-missense-codon" : ""
    def cacheCommand = evCacheCommand(gene_id, 'codon', codon_msa, score_options)
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    python3 ${projectDir}/../evmutation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --codon-msa "${codon_msa}" \\
        ${cmparamsArg} \\
        ${missenseArg} \\
        ${plmcArg} \\
        --output . \\
        ${logArg}

    if [ -d "${gene_id}/EVmutation" ]; then
        mv ${gene_id}/EVmutation/${gene_id}.codon.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.codon_model_params" -exec mv {} . \\; 2>/dev/null || true
    ${cacheCommand}
    """
}

process plan_adabmdca {
    executor params.resource_executor
    cpus 1
    memory '512 MB'
    tag { "${gene_id} ${side} plan" }
    publishDir "${params.output_dir}/resource_plans", mode: 'copy', pattern: '*.plan.json'

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv)
    path resource_config

    output:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv),
          path("${gene_id}.${side}.plan.json")

    script:
    def arguments = [
        'python3', "${projectDir}/resource_planner.py", 'plan',
        '--gene', gene_id, '--side', side, '--msa', msa_file,
        '--fasta', fasta_file, '--mutations', mutations_csv,
        '--config', resource_config, '--output', "${gene_id}.${side}.plan.json",
    ]
    if (params.resource_errors)
        arguments.add('--defer-errors')
    def model_params = side == 'protein' ? params.adabmdca_protein_params : params.adabmdca_codon_params
    if (model_params)
        arguments.addAll(['--params', model_params])
    if (side == 'protein' && params.skip_codon_adabmdca)
        arguments.add('--skip-codon')
    def command = arguments.collect { shellQuote(it) }.join(' ')
    """
    set -euo pipefail
    ${command}
    """
}

process run_adabmdca_gpu {
    executor params.resource_executor
    maxForks Math.max(1, params.gpu_slots.toInteger())
    cpus { plan.threads as int }
    memory { "${plan.gpu_host_memory_gib} GB" }
    tag { "${gene_id} ${side} GPU" }
    publishDir { "${params.output_dir}/${gene_id}/adabmDCA" }, mode: 'copy', pattern: '*.{tsv,complete.json}'
    publishDir { "${params.output_dir}/adabmdca_${side}_params" }, mode: 'copy', pattern: '*_adabm_params'

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv),
          path(plan_file), val(plan)

    output:
    tuple val(gene_id), val(side), path("${gene_id}.${side}.tsv"),
          path("${gene_id}.${side}_adabm_params"), path("${gene_id}.${side}.complete.json"),
          emit: complete, optional: true
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv),
          path(plan_file), val(plan), path("${gene_id}.${side}.gpu_oom.json"),
          emit: gpu_oom, optional: true

    script:
    def command = adabmdcaCommand(side, fasta_file, msa_file, mutations_csv, plan_file, plan, 'cuda', task.cpus)
    """
    set -euo pipefail
    ${command}
    """
}

process run_adabmdca_cpu {
    executor params.resource_executor
    cpus { plan.threads as int }
    memory { "${plan.cpu_memory_gib} GB" }
    tag { "${gene_id} ${side} CPU" }
    publishDir { "${params.output_dir}/${gene_id}/adabmDCA" }, mode: 'copy', pattern: '*.{tsv,complete.json}'
    publishDir { "${params.output_dir}/adabmdca_${side}_params" }, mode: 'copy', pattern: '*_adabm_params'

    input:
    tuple val(gene_id), val(side), path(fasta_file), path(msa_file), path(mutations_csv),
          path(plan_file), val(plan)

    output:
    tuple val(gene_id), val(side), path("${gene_id}.${side}.tsv"),
          path("${gene_id}.${side}_adabm_params"), path("${gene_id}.${side}.complete.json")

    script:
    def command = adabmdcaCommand(side, fasta_file, msa_file, mutations_csv, plan_file, plan, 'cpu', task.cpus)
    """
    set -euo pipefail
    ${command}
    """
}
