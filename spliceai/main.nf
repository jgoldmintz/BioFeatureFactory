#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------- PARAMETERS ----------------
params.help = false
params.mutations_dir = null
params.input_vcf_dir = null
params.input_vcf_file = null
params.skip_vcf_generation = false
params.splice_threshold = 0.0
params.output_dir = "."
params.reference_genome = null
params.annotation_file = null
params.chromosome_mapping_dir = null
params.transcript_mapping_dir = null
params.genomic_mapping_dir = null
params.vcf_output_dir = null
params.validation_log = null
params.clear_vcf_cache = false
params.validate_mapping = false
params.chromosome_format = 'refseq'
params.retry_jitter = 10
params.maxforks = 0 // default: unbounded

if (params.help) {
    println """
Usage: nextflow run main.nf [options]
"""
    exit 0
}

// ---------------- VALIDATION ----------------
def cpu_count = Runtime.runtime.availableProcessors()
if (params.maxforks.toInteger() > cpu_count)
    error "ERROR: --maxforks (${params.maxforks}) cannot exceed available CPUs (${cpu_count})"

// ---------------- WORKFLOW ----------------
workflow {

    if (!params.reference_genome) error "ERROR: --reference_genome is required"
    if (!params.annotation_file) error "ERROR: --annotation_file is required"
    if (!params.transcript_mapping_dir) error "ERROR: --transcript_mapping_dir is required"
    if (!params.genomic_mapping_dir) error "ERROR: --genomic_mapping_dir is required"

    println "DEBUG: reference_genome = ${params.reference_genome}"
    println "DEBUG: annotation_file = ${params.annotation_file}"
    println "DEBUG: transcript_mapping_dir = ${params.transcript_mapping_dir}"
    println "DEBUG: input_vcf_dir = ${params.input_vcf_dir}"
    println "DEBUG: input_vcf_file = ${params.input_vcf_file}"
    println "DEBUG: genomic_mapping_dir = ${params.genomic_mapping_dir}"
    println "DEBUG: validation_log = ${params.validation_log}"
    println "DEBUG: clear_vcf_cache = ${params.clear_vcf_cache}"
    println "DEBUG: chromosome_format = ${params.chromosome_format}"
    println "DEBUG: maxforks = ${params.maxforks} (CPU count: ${cpu_count})"

    def vcf_source

    if (params.input_vcf_file) {
        vcf_source = Channel
            .fromPath(params.input_vcf_file)
            .map { vcf -> tuple(vcf.baseName.replaceAll(/\.vcf$/, ''), vcf) }
    }
    else if (params.skip_vcf_generation || params.input_vcf_dir) {
        if (!params.input_vcf_dir)
            error "Must specify --input_vcf_dir when using --skip_vcf_generation"
        vcf_source = Channel
            .fromPath("${params.input_vcf_dir}/*.vcf")
            .map { vcf -> tuple(vcf.baseName.replaceAll(/\.vcf$/, ''), vcf) }
    }
    else {
        if (!params.mutations_dir)
            error "Must specify --mutations_dir when not skipping VCF generation"

        mutation_files = Channel
            .fromPath("${params.mutations_dir}/*.csv")
            .map { csv -> tuple(csv.baseName.replaceAll(/_mutations$/, ''), csv) }

        vcf_source = generate_vcfs(mutation_files)
    }

    // enforce sequential chaining
    compress_and_index(vcf_source)
    run_spliceai(compress_and_index.out, file(params.reference_genome), file(params.annotation_file))
    parse_results(
        run_spliceai.out,
        file(params.transcript_mapping_dir),
        file(params.chromosome_mapping_dir ?: ""),
        params.splice_threshold ?: 0.0,
        params.validation_log ?: ""
    )
}

// ---------------- PROCESSES ----------------
process generate_vcfs {
    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.vcf', enabled: params.vcf_output_dir != null

    input:
    tuple val(gene_id), path(mutations_csv)

    output:
    tuple val(gene_id), path("${gene_id}.vcf")

    script:
    gene_id = mutations_csv.baseName.replaceAll(/_mutations$/, '')
    def mappingArg = params.chromosome_mapping_dir ? "--chromosome-mapping-input ${params.chromosome_mapping_dir}" : ""
    def logArg = params.validation_log ? "--log ${params.validation_log}" : ""
    def validateArg = params.validate_mapping ? "--validate-mapping" : ""
    def clearCacheArg = params.clear_vcf_cache ? "--clear-cache" : ""
    def optionalFlags = [mappingArg, logArg, validateArg, clearCacheArg].findAll { it }.join(' ')
    def chromFormat = params.chromosome_format ?: 'refseq'
    def options = optionalFlags ? " ${optionalFlags}" : ""
    """
    python3 ${projectDir}/../dependencies/vcf_converter.py \\
        -m ${mutations_csv} \\
        -o . \\
        --chromosome-format ${chromFormat} \\
        -r ${params.reference_genome} \\
        -a ${params.annotation_file}${options}
    """
}

process compress_and_index {
    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.vcf.gz*', enabled: params.vcf_output_dir != null

    input:
    tuple val(gene_id), path(vcf)

    output:
    tuple val(gene_id), path("${gene_id}.vcf.gz"), path("${gene_id}.vcf.gz.tbi")

    script:
    """
    bgzip -c ${vcf} > ${gene_id}.vcf.gz
    tabix -p vcf ${gene_id}.vcf.gz
    """
}

process run_spliceai {
    maxForks params.maxforks.toInteger()
    errorStrategy 'retry'
    maxRetries 3

    publishDir params.vcf_output_dir ?: '.', mode: 'copy', pattern: '*.spliceai.vcf*', enabled: params.vcf_output_dir != null

    input:
    tuple val(gene_id), path(vcf_gz), path(vcf_tbi)
    val reference_genome
    val annotation_file

    output:
    tuple val(gene_id), path("${gene_id}.spliceai.vcf")

    script:
    def stagger = new Random(task.hashCode()).nextInt(15) + 5
    def jitter = (task.attempt > 1) ? new Random().nextInt(params.retry_jitter ?: 10) : 0

    """
    echo "[run_spliceai] task ${task.name} delaying ${stagger}s (attempt ${task.attempt})"
    sleep ${stagger}

    if [[ ${task.attempt} -gt 1 ]]; then
        echo "[run_spliceai] retrying after jitter delay ${jitter}s"
        sleep ${jitter}
    fi

    export TF_CPP_MIN_LOG_LEVEL=3
    export TF_NUM_INTEROP_THREADS=1
    export TF_NUM_INTRAOP_THREADS=1
    export OMP_NUM_THREADS=1

    spliceai \\
        -I ${vcf_gz} \\
        -O ${gene_id}.spliceai.vcf \\
        -R ${reference_genome} \\
        -A ${annotation_file}
    """
}

process parse_results {
    publishDir "${params.output_dir ?: '.'}", mode: 'copy'

    input:
    tuple val(gene_id), path(spliceai_vcf)
    val transcript_mapping_dir
    val chromosome_mapping_dir
    val splice_threshold
    val validation_log

    output:
    tuple val(gene_id), path("${gene_id}_spliceai_results.tsv")

    script:
    def chromosome_mapping_arg = chromosome_mapping_dir ? "--chromosome-mapping ${chromosome_mapping_dir}/combined_${gene_id}.csv" : ""
    def log_arg = validation_log ? "--log ${validation_log}" : ""
    def optional_flags = [chromosome_mapping_arg, log_arg].findAll { it }
    def extras = optional_flags ? " ${optional_flags.join(' ')}" : ""
    """
    python3 ${projectDir}/spliceai-parser.py \\
        --input ${spliceai_vcf} \\
        --output ${gene_id}_spliceai_results.tsv \\
        --transcript-mapping ${transcript_mapping_dir}/combined_${gene_id}.csv \\
        --threshold ${splice_threshold}${extras}
    """
}
