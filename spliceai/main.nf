#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------- PARAMETERS ----------------
params.help                       = false

// input sources
params.mutations_path             = null       // REQUIRED if generating VCFs. Dir or file.
params.mutations_dir              = null       // legacy alias

// prebuilt VCFs
params.input_vcf_dir              = null       // dir of *.vcf
params.input_vcf_file             = null       // single *.vcf
params.skip_vcf_generation        = false

// references
params.reference_genome           = null       // FASTA for SpliceAI -R
params.annotation_file            = null       // GTF/GFF for SpliceAI -A

// mappings (dir OR file)
params.transcript_mapping_path    = null       // REQUIRED. Dir or file.
params.transcript_mapping_dir     = null       // legacy alias
params.chromosome_mapping_path    = null       // REQUIRED. Dir or file.
params.chromosome_mapping_dir     = null       // legacy alias
params.genomic_mapping_path       = null       // REQUIRED. Dir or file. (validated, not consumed here)
params.genomic_mapping_dir        = null       // legacy alias

// run control
params.splice_threshold           = 0.0
params.output_dir                 = "."
params.vcf_output_dir             = null
params.validation_log             = null
params.clear_vcf_cache            = false
params.validate_mapping           = false
params.chromosome_format          = 'refseq'
params.retry_jitter               = 10
params.maxforks                   = 0          // 0 = unbounded

// ---------------- HELP ----------------
if (params.help) {
    println """
Usage:
  nextflow run main.nf \\
    --reference_genome /path/genome.fa \\
    --annotation_file /path/annotations.gtf \\
    --transcript_mapping_path /path/to/transcript_maps/ \\
    --chromosome_mapping_path /path/to/chromosome_maps/ \\
    --genomic_mapping_path /path/to/genomic_maps/ \\
    [ --mutations_path MUTS | --input_vcf_dir VCF_DIR | --input_vcf_file VCF ] \\
    [options]

Inputs
  --mutations_path       Dir with <GENE>_mutations.csv files OR a single <GENE>_mutations.csv.
  --input_vcf_dir        Dir with one or more *.vcf. Implies --skip_vcf_generation.
  --input_vcf_file       Single VCF. Implies --skip_vcf_generation.

Mappings (dir or file)
  --transcript_mapping_path   Required. If a directory, per gene pick: exact '<GENE>.csv', then 'combined_<GENE>.csv', else contains '*<GENE>*.csv' (case-insensitive).
  --chromosome_mapping_path   Required. Same resolution rules.
  --genomic_mapping_path      Required. Same resolution rules. Validated but not passed to parser in this version.

Other options
  --vcf_output_dir       Publish intermediate VCFs.
  --output_dir           Publish final *_spliceai_results.tsv. Default: current dir.
  --validation_log       Path to a log file for parser.
  --splice_threshold     Minimum SpliceAI delta score. Default: 0.0.
  --chromosome_format    vcf_converter chromosome naming. Default: refseq.
  --validate_mapping     Adds --validate-mapping to vcf_converter.py.
  --clear_vcf_cache      Adds --clear-cache to vcf_converter.py.
  --maxforks INT         Max concurrent SpliceAI tasks for process run_spliceai. 0 means unbounded.

Notes
  Directory-mode mapping selection is deterministic and errors if no match.
"""
    exit 0
}

// ---------------- VALIDATION / NORMALIZATION ----------------
// legacy aliases
if (!params.transcript_mapping_path && params.transcript_mapping_dir)
    params.transcript_mapping_path = params.transcript_mapping_dir
if (!params.chromosome_mapping_path && params.chromosome_mapping_dir)
    params.chromosome_mapping_path = params.chromosome_mapping_dir
if (!params.genomic_mapping_path && params.genomic_mapping_dir)
    params.genomic_mapping_path = params.genomic_mapping_dir
if (!params.mutations_path && params.mutations_dir)
    params.mutations_path = params.mutations_dir

// required checks
if (!params.reference_genome)          error "ERROR: --reference_genome is required"
if (!params.annotation_file)           error "ERROR: --annotation_file is required"
if (!params.transcript_mapping_path)   error "ERROR: --transcript_mapping_path is required"
if (!params.chromosome_mapping_path)   error "ERROR: --chromosome_mapping_path is required"
if (!params.genomic_mapping_path)      error "ERROR: --genomic_mapping_path is required"

// concurrency guard
def cpu_count = Runtime.runtime.availableProcessors()
if (params.maxforks.toInteger() > cpu_count)
    error "ERROR: --maxforks (${params.maxforks}) cannot exceed available CPUs (${cpu_count})"

// ---------------- WORKFLOW ----------------
workflow {
    println "DEBUG: reference_genome = ${params.reference_genome}"
    println "DEBUG: annotation_file = ${params.annotation_file}"
    println "DEBUG: transcript_mapping_path = ${params.transcript_mapping_path}"
    println "DEBUG: chromosome_mapping_path = ${params.chromosome_mapping_path}"
    println "DEBUG: genomic_mapping_path = ${params.genomic_mapping_path}"
    println "DEBUG: mutations_path = ${params.mutations_path}"
    println "DEBUG: input_vcf_dir = ${params.input_vcf_dir}"
    println "DEBUG: input_vcf_file = ${params.input_vcf_file}"
    println "DEBUG: splice_threshold = ${params.splice_threshold}"
    println "DEBUG: validation_log = ${params.validation_log}"
    println "DEBUG: clear_vcf_cache = ${params.clear_vcf_cache}"
    println "DEBUG: chromosome_format = ${params.chromosome_format}"
    println "DEBUG: maxforks = ${params.maxforks} (CPU count: ${cpu_count})"
    println "DEBUG: skip_vcf_generation = ${params.skip_vcf_generation}"

    // helper closure to resolve mapping CSV for a given gene
    def resolveMap = { String pathParam, String geneId, boolean allowMissing ->
        if (!pathParam)
            return [null, false]

        File base = new File(pathParam)
        if (!base.exists()) {
            if (allowMissing)
                return [null, false]
            throw new RuntimeException("ERROR: mapping path ${pathParam} not found")
        }

        if (!base.isDirectory())
            return [ base.getAbsolutePath(), true ]

        def geneLower = geneId.toLowerCase()
        def files = (base.listFiles() ?: []) as List

        def exact = files.findAll { f ->
            def n = f.name.toLowerCase()
            n.endsWith(".csv") && n == "${geneLower}.csv"
        }.sort { it.name.toLowerCase() }
        if (exact)
            return [ exact[0].getAbsolutePath(), true ]

        def combined = files.findAll { f ->
            def n = f.name.toLowerCase()
            n.endsWith(".csv") && n == "combined_${geneLower}.csv"
        }.sort { it.name.toLowerCase() }
        if (combined)
            return [ combined[0].getAbsolutePath(), true ]

        def contains = files.findAll { f ->
            def n = f.name.toLowerCase()
            n.endsWith(".csv") && n.contains(geneLower)
        }.sort { it.name.toLowerCase() }
        if (contains)
            return [ contains[0].getAbsolutePath(), true ]

        if (allowMissing)
            return [null, false]
        throw new RuntimeException("ERROR: no mapping CSV in ${pathParam} for gene ${geneId} (exact|combined|contains)")
    }

    // choose VCF source
    def vcf_source
    if (params.input_vcf_file) {
        vcf_source = Channel
            .fromPath(params.input_vcf_file)
            .map { v -> tuple(v.baseName.replaceAll(/\.vcf$/, ''), v) }
    }
    else if (params.skip_vcf_generation || params.input_vcf_dir) {
        if (!params.input_vcf_dir)
            error "Must specify --input_vcf_dir when using --skip_vcf_generation without --input_vcf_file"
        vcf_source = Channel
            .fromPath("${params.input_vcf_dir}/*.vcf")
            .map { v -> tuple(v.baseName.replaceAll(/\.vcf$/, ''), v) }
    }
    else {
        if (!params.mutations_path)
            error "Must specify --mutations_path when not skipping VCF generation"
        def mutBase = new File(params.mutations_path as String)
        if (!mutBase.exists())
            error "ERROR: --mutations_path ${params.mutations_path} not found"

        def mutation_files_ch
        if (mutBase.isDirectory()) {
            mutation_files_ch = Channel
                .fromPath("${params.mutations_path}/*.csv")
                .map { csv ->
                    def gene = csv.baseName.replaceAll(/_mutations$/, '')
                    def (c_map_path, c_ok) = resolveMap(params.chromosome_mapping_path, gene, false)
                    if (!c_ok)
                        throw new RuntimeException("ERROR: chromosome mapping not resolved for ${gene}")
                    tuple(gene, csv, c_map_path)
                }
        } else {
            mutation_files_ch = Channel
                .fromPath(params.mutations_path)
                .map { csv ->
                    def gene = csv.baseName.replaceAll(/_mutations$/, '')
                    def (c_map_path, c_ok) = resolveMap(params.chromosome_mapping_path, gene, false)
                    if (!c_ok)
                        throw new RuntimeException("ERROR: chromosome mapping not resolved for ${gene}")
                    tuple(gene, csv, c_map_path)
                }
        }
        vcf_source = generate_vcfs(mutation_files_ch)
    }

    // compress and SpliceAI
    compress_and_index(vcf_source)
    run_spliceai(
        compress_and_index.out,
        file(params.reference_genome),
        file(params.annotation_file)
    )

    // build inputs for parser; enforce transcript+chromosome presence per gene
    def parser_in = run_spliceai.out.map { gene_id, spliceai_vcf ->
        def (t_map, t_ok) = resolveMap(params.transcript_mapping_path, gene_id, false)
        if (!t_ok)
            throw new RuntimeException("ERROR: transcript mapping not resolved for ${gene_id}")

        def (c_map, c_ok) = resolveMap(params.chromosome_mapping_path, gene_id, false)
        if (!c_ok)
            throw new RuntimeException("ERROR: chromosome mapping not resolved for ${gene_id}")

        def (g_map, g_ok) = resolveMap(params.genomic_mapping_path, gene_id, false)
        if (!g_ok)
            throw new RuntimeException("ERROR: genomic mapping not resolved for ${gene_id}")

        tuple(gene_id, spliceai_vcf, file(t_map), file(c_map))
    }

    parse_results(
        parser_in,
        params.splice_threshold ?: 0.0,
        params.validation_log ?: ""
    )
}

// ---------------- PROCESSES ----------------
process generate_vcfs {
    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.vcf', enabled: params.vcf_output_dir != null
    tag { gene_id }

    input:
    tuple val(gene_id), path(mutations_csv), val(chrom_map_path)

    output:
    tuple val(gene_id), path("${gene_id}.vcf")

    script:
    def mapArg      = chrom_map_path ? "--chromosome-mapping-input \"${chrom_map_path}\"" : ""
    def logArg      = params.validation_log          ? "--log \"${params.validation_log}\"" : ""
    def validateArg = params.validate_mapping        ? "--validate-mapping" : ""
    def clearArg    = params.clear_vcf_cache         ? "--clear-cache" : ""
    def extras      = [mapArg, logArg, validateArg, clearArg].findAll{ it }.join(' ')
    def chromFormat = params.chromosome_format ?: 'refseq'
    """
    set -euo pipefail
    python3 ${projectDir}/../dependencies/vcf_converter.py \\
      -m "${mutations_csv}" \\
      -o . \\
      --chromosome-format "${chromFormat}" \\
      -r "${params.reference_genome}" \\
      -a "${params.annotation_file}" \\
      ${extras}
    """
}

process compress_and_index {
    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.vcf.gz*', enabled: params.vcf_output_dir != null
    tag { gene_id }

    input:
    tuple val(gene_id), path(vcf)

    output:
    tuple val(gene_id), path("${gene_id}.vcf.gz"), path("${gene_id}.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bgzip -c "${vcf}" > "${gene_id}.vcf.gz"
    tabix -p vcf "${gene_id}.vcf.gz"
    """
}

process run_spliceai {
    maxForks params.maxforks.toInteger()
    errorStrategy 'retry'
    maxRetries 3
    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.spliceai.vcf*', enabled: params.vcf_output_dir != null
    tag { gene_id }

    input:
    tuple val(gene_id), path(vcf_gz), path(vcf_tbi)
    val reference_genome
    val annotation_file

    output:
    tuple val(gene_id), path("${gene_id}.spliceai.vcf")

    script:
    def stagger = new Random(task.hashCode()).nextInt(15) + 5
    def jitter  = (task.attempt > 1) ? new Random().nextInt(params.retry_jitter ?: 10) : 0
    """
    set -euo pipefail
    echo "[run_spliceai] ${task.name} delay ${stagger}s (attempt ${task.attempt})"
    sleep ${stagger}
    if [[ ${task.attempt} -gt 1 ]]; then
      echo "[run_spliceai] retry jitter ${jitter}s"
      sleep ${jitter}
    fi

    export TF_CPP_MIN_LOG_LEVEL=3
    export TF_NUM_INTEROP_THREADS=1
    export TF_NUM_INTRAOP_THREADS=1
    export OMP_NUM_THREADS=1

    spliceai \\
      -I "${vcf_gz}" \\
      -O "${gene_id}.spliceai.vcf" \\
      -R "${reference_genome}" \\
      -A "${annotation_file}"
    """
}

process parse_results {
    publishDir "${params.output_dir ?: '.'}", mode: 'copy'
    stageInMode 'copy'
    scratch true
    tag { gene_id }

    input:
    tuple val(gene_id),
          path(spliceai_vcf),
          path(transcript_map),
          path(chromosome_map)
    val splice_threshold
    val validation_log

    output:
    tuple val(gene_id), path("${gene_id}_spliceai_results.tsv")

    script:
    """
    set -euo pipefail

    mkdir -p stage
    cp -f "${spliceai_vcf}"    stage/in.vcf
    cp -f "${transcript_map}"  stage/transcript.csv
    cp -f "${chromosome_map}"  stage/chrom.csv

    echo "[parse_results] gene=${gene_id}"
    echo "[parse_results] in.vcf lines: \$(wc -l < stage/in.vcf) || true"
    echo "[parse_results] transcript.csv lines: \$(wc -l < stage/transcript.csv) || true"
    echo "[parse_results] chrom.csv lines: \$(wc -l < stage/chrom.csv) || true"

    ARGS=( --input stage/in.vcf
           --output "${gene_id}_spliceai_results.tsv"
           --transcript-mapping stage/transcript.csv
           --chromosome-mapping stage/chrom.csv
           --threshold ${splice_threshold} )

    if [[ -n "${validation_log}" ]]; then
      ARGS+=( --log "${validation_log}" )
    fi

    python3 ${projectDir}/spliceai-parser.py "\${ARGS[@]}"

    echo "[parse_results] out.tsv lines: \$(wc -l < "${gene_id}_spliceai_results.tsv") || true"
    """
}
