  #!/usr/bin/env nextflow
// BioFeatureFactory
// Copyright (C) 2023–2026  Jacob Goldmintz
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

  nextflow.enable.dsl = 2

  // ---------------- PARAMETERS ----------------
  params.mutations_path             = null
  params.mutations_dir              = null
  params.input_vcf_dir              = null
  params.input_vcf_file             = null
  params.skip_vcf_generation        = false
  params.reference_genome           = null
  params.annotation_file            = null
  params.transcript_mapping_path    = null
  params.transcript_mapping_dir     = null
  params.chromosome_mapping_path    = null
  params.chromosome_mapping_dir     = null
  params.genomic_mapping_path       = null
  params.genomic_mapping_dir        = null
  params.splice_threshold           = 0.0
  params.output_dir                 = "."
  params.vcf_output_dir             = null
  params.validation_log             = null
  params.clear_vcf_cache            = false
  params.validate_mapping           = false
  params.chromosome_format          = 'refseq'
params.retry_jitter               = 10
params.maxforks                   = 0
params.forceAll_isoforms          = false
params.max_isoforms_per_gene      = 50
// Dedicated conda env for spliceai, created by bootstrap step 6d. Set to '' (or
// pass --spliceai_env '') to run spliceai from whatever env launched Nextflow.
// See the note above the spliceai invocation in run_spliceai for why it has one.
params.spliceai_env               = 'bff-spliceai'

  // Nextflow 26+ forbids statements outside a process / workflow / function body,
  // so the legacy aliases, the required checks and the concurrency guard all live
  // in functions now and are called from the workflow block. As bare top-level
  // `if`s they failed compilation before a single task ran:
  //     Error main.nf:47:3: Statements cannot be mixed with script declarations
  // Nextflow reports only the FIRST offender, so all three blocks had to move --
  // fixing the aliases alone would have surfaced the checks next, then the guard.
  // mutation_effects/bin/main.nf already carries the same fix (validateRequiredParams).
  def applyLegacyAliases() {
      if (!params.transcript_mapping_path && params.transcript_mapping_dir)
          params.transcript_mapping_path = params.transcript_mapping_dir
      if (!params.chromosome_mapping_path && params.chromosome_mapping_dir)
          params.chromosome_mapping_path = params.chromosome_mapping_dir
      if (!params.genomic_mapping_path && params.genomic_mapping_dir)
          params.genomic_mapping_path = params.genomic_mapping_dir
      if (!params.mutations_path && params.mutations_dir)
          params.mutations_path = params.mutations_dir
  }

  def validateRequiredParams() {
      if (!params.reference_genome)        error "ERROR: --reference_genome is required"
      if (!params.annotation_file)         error "ERROR: --annotation_file is required"
      if (!params.transcript_mapping_path) error "ERROR: --transcript_mapping_path is required"
      if (!params.chromosome_mapping_path) error "ERROR: --chromosome_mapping_path is required"
      // --genomic_mapping_path is deliberately NOT required. It was required,
      // resolved per gene, and then dropped: the tuple handed to parse_results
      // carries only the transcript and chromosome maps, and spliceai-parser.py
      // declares no --genomic-mapping argument at all. So the check refused runs
      // over a file the pipeline had no way to read. The param is still declared
      // and still accepted (spliceai_pipeline_controller.py always passes it), it
      // just no longer gates the run. Re-add the check together with a parser
      // argument that consumes it, never on its own.
  }

  def validateConcurrency() {
      def cpu_count = Runtime.runtime.availableProcessors()
      if (params.maxforks.toInteger() > cpu_count)
          error "ERROR: --maxforks (${params.maxforks}) cannot exceed available CPUs (${cpu_count})"
  }

  // Top-level FUNCTION, not a `def resolveMap = { ... }` closure inside the
  // workflow block. Nextflow 26's parser does not resolve a workflow-local
  // closure from inside the nested `.map { }` closures that call it, so the
  // previous form failed compilation at every call site:
  //     Error main.nf:168:48: `resolveMap` is not defined
  // mutation_effects/bin/main.nf uses the same top-level-def pattern for
  // resolveMutationCsv / resolveMsaFile.
  def resolveMap(String pathParam, String geneId, boolean allowMissing) {
          if (!pathParam) return [null, false]

          File base = new File(pathParam)
          if (!base.exists()) {
              if (allowMissing) return [null, false]
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
          if (exact) return [ exact[0].getAbsolutePath(), true ]

          def combined = files.findAll { f ->
              def n = f.name.toLowerCase()
              n.endsWith(".csv") && n == "combined_${geneLower}.csv"
          }.sort { it.name.toLowerCase() }
          if (combined) return [ combined[0].getAbsolutePath(), true ]

          def contains = files.findAll { f ->
              def n = f.name.toLowerCase()
              n.endsWith(".csv") && n.contains(geneLower)
          }.sort { it.name.toLowerCase() }
          if (contains) return [ contains[0].getAbsolutePath(), true ]

          if (allowMissing) return [null, false]
          throw new RuntimeException("ERROR: no mapping CSV in ${pathParam} for gene ${geneId}")
      }

  // ---------------- WORKFLOW ----------------
  workflow {
      // Aliases FIRST: validateRequiredParams checks *_mapping_path, which the
      // legacy *_mapping_dir forms are only mapped onto by applyLegacyAliases.
      applyLegacyAliases()
      validateRequiredParams()
      validateConcurrency()

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
      println "DEBUG: maxforks = ${params.maxforks} (CPU count: ${Runtime.runtime.availableProcessors()})"
    println "DEBUG: skip_vcf_generation = ${params.skip_vcf_generation}"
    println "DEBUG: forceAll_isoforms = ${params.forceAll_isoforms}"
    println "DEBUG: max_isoforms_per_gene = ${params.max_isoforms_per_gene}"


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

      compress_and_index(vcf_source)
    run_spliceai(
        compress_and_index.out,
        file(params.reference_genome),
        file(params.annotation_file),
        params.forceAll_isoforms,
        params.max_isoforms_per_gene
    )

      def parser_in = run_spliceai.out.map { gene_id, spliceai_vcf ->
          def (t_map, t_ok) = resolveMap(params.transcript_mapping_path, gene_id, false)
          if (!t_ok)
              throw new RuntimeException("ERROR: transcript mapping not resolved for ${gene_id}")

          def (c_map, c_ok) = resolveMap(params.chromosome_mapping_path, gene_id, false)
          if (!c_ok)
              throw new RuntimeException("ERROR: chromosome mapping not resolved for ${gene_id}")

          // No genomic-mapping resolution here: nothing downstream consumes it.
          // Resolving it threw on a missing file that would never have been read.
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
      publishDir { params.vcf_output_dir ? params.vcf_output_dir : "${params.output_dir ?: '.'}/${gene_id}/vcf" }, mode: 'copy', pattern: '*/vcf/*.vcf', saveAs: { fn -> file(fn).name }
      tag { gene_id }

      input:
      tuple val(gene_id), path(mutations_csv), val(chrom_map_path)

      output:
      tuple val(gene_id), path("${gene_id}/vcf/${gene_id}.vcf")

      script:
      def mapArg      = chrom_map_path ? "--chromosome-mapping-input \"${chrom_map_path}\"" : ""
      def logArg      = params.validation_log          ? "--log \"${params.validation_log}\"" : ""
      def validateArg = params.validate_mapping        ? "--validate-mapping" : ""
      def clearArg    = params.clear_vcf_cache         ? "--clear-cache" : ""
      def extras      = [mapArg, logArg, validateArg, clearArg].findAll{ it }.join(' ')
      def chromFormat = params.chromosome_format ?: 'refseq'
      """
      set -euo pipefail
      python3 -m biofeaturefactory.core.vcf_converter \\
        -m "${mutations_csv}" \\
        -o . \\
        --chromosome-format "${chromFormat}" \\
        -r "${params.reference_genome}" \\
        -a "${params.annotation_file}" \\
        ${extras}
      """
  }

  process compress_and_index {
      publishDir { params.vcf_output_dir ? params.vcf_output_dir : "${params.output_dir ?: '.'}/${gene_id}/vcf" }, mode: 'copy', pattern: '*.vcf.gz*'
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
    publishDir { params.vcf_output_dir ? params.vcf_output_dir : "${params.output_dir ?: '.'}/${gene_id}/vcf" }, mode: 'copy', pattern: '*.spliceai.vcf*'
    tag { gene_id }

    input:
    tuple val(gene_id), path(vcf_gz), path(vcf_tbi)
    val reference_genome
    val annotation_file
    val forceAll_isoforms
    val max_isoforms_per_gene

    output:
    tuple val(gene_id), path("${gene_id}.spliceai.vcf")

    script:
    """
    set -euo pipefail

    export TF_CPP_MIN_LOG_LEVEL=3
    export TF_NUM_INTEROP_THREADS=1
    export TF_NUM_INTRAOP_THREADS=1
    export OMP_NUM_THREADS=1

    # Check isoform count and apply filtering if needed
    ISOFORM_COUNT=\$(grep -c "^${gene_id}\t" "${annotation_file}" || echo 0)
    echo "[run_spliceai] ${gene_id}: \$ISOFORM_COUNT isoforms detected"

    if [[ "${forceAll_isoforms}" == "false" ]] && [[ \$ISOFORM_COUNT -gt ${max_isoforms_per_gene} ]]; then
      echo "[run_spliceai] ${gene_id}: Applying hybrid filter (\$ISOFORM_COUNT > ${max_isoforms_per_gene})"
      python3 ${projectDir}/filter_annotation.py \\
        --annotation "${annotation_file}" \\
        --gene "${gene_id}" \\
        --max-isoforms ${max_isoforms_per_gene} \\
        --output "${gene_id}_filtered_annotation.txt" \\
        --log-file filter.log
      ANNOTATION_TO_USE="${gene_id}_filtered_annotation.txt"
    else
      echo "[run_spliceai] ${gene_id}: Using full annotation (no filtering)"
      ANNOTATION_TO_USE="${annotation_file}"
    fi

    # NOT the `spliceai` console script. It imports pandas -> pyarrow BEFORE
    # tensorflow, and pyarrow 25.0.1 and tensorflow 2.21.0 each ship their own
    # Abseil. Whichever loads first owns the symbols, and when pyarrow wins,
    # tf.data's prefetch CondVar inside Keras `model.predict` is waited on by
    # libarrow's absl and signalled by tensorflow's -- it never wakes. Measured:
    # a real SMN2 run burned 6.31s of CPU in 38 MINUTES of wall clock at 0.0%
    # CPU, blocked in PrefetchDatasetOp::GetNextInternal -> absl::CondVar::
    # WaitCommon -> AbslInternalPerThreadSemWait (in libarrow.2500.dylib).
    #
    # Importing tensorflow FIRST is the whole fix; pyarrow may then load freely,
    # and it is kept here as a belt-and-braces measure even when --spliceai_env
    # points at an env with no pyarrow at all.
    # Reduced to a 6-line repro, run both ways in this env:
    #     import pyarrow    -> model.predict()  TIMEOUT (deadlock)
    #     import tensorflow -> model.predict()  0.06s, and 0.03s AFTER pyarrow
    # Unrelated to numpy: reproduced identically on numpy 1.26.4 and 2.4.6.
    #
    # main() is spliceai's declared console entry point ('spliceai ->
    # spliceai.__main__:main'), and argv[0] is reset so its --help/error text
    # still names the tool. Restore the plain `spliceai` call only once pyarrow
    # and tensorflow agree on one Abseil.
    # Run in the dedicated spliceai env when one is configured. `conda run` is
    # used rather than `conda activate` because a process script is a
    # non-interactive shell: `conda activate` needs conda's shell function, which
    # only exists after `conda shell.<sh> hook`, and sourcing conda.sh from an
    # unknown install prefix is not portable. `conda run` needs neither.
    #
    # The env is VERIFIED, not assumed: if --spliceai_env names an env that does
    # not exist, this fails here with a clear message instead of silently falling
    # back to the shared env and reintroducing the Abseil deadlock.
    SPLICEAI_ENV="${params.spliceai_env ?: ''}"
    SPLICEAI_RUN=""
    if [[ -n "\$SPLICEAI_ENV" ]]; then
      CONDA_BIN=""
      command -v conda >/dev/null 2>&1 && CONDA_BIN=conda
      [[ -z "\$CONDA_BIN" ]] && command -v mamba >/dev/null 2>&1 && CONDA_BIN=mamba
      if [[ -z "\$CONDA_BIN" ]]; then
        echo "[run_spliceai] ERROR: --spliceai_env '\$SPLICEAI_ENV' requested but no conda/mamba on PATH" >&2
        exit 1
      fi
      if ! \$CONDA_BIN env list | awk '{print \$1}' | grep -qx "\$SPLICEAI_ENV"; then
        echo "[run_spliceai] ERROR: conda env '\$SPLICEAI_ENV' not found. Create it with" >&2
        echo "               scripts/bootstrap.sh (step 6d), or pass --spliceai_env '' to use the current env." >&2
        exit 1
      fi
      SPLICEAI_RUN="\$CONDA_BIN run --no-capture-output -n \$SPLICEAI_ENV"
      echo "[run_spliceai] ${gene_id}: using conda env \$SPLICEAI_ENV"
    else
      echo "[run_spliceai] ${gene_id}: using the launching environment (--spliceai_env is empty)"
    fi

    \$SPLICEAI_RUN python3 -c "import tensorflow, sys; sys.argv[0] = 'spliceai'; from spliceai.__main__ import main; main()" \\
      -I "${vcf_gz}" \\
      -O "${gene_id}.spliceai.vcf" \\
      -R "${reference_genome}" \\
      -A "\$ANNOTATION_TO_USE"
    """
}

  process parse_results {
      publishDir { "${params.output_dir ?: '.'}/${gene_id}/SpliceAI" }, mode: 'copy'
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
      tuple val(gene_id), path("${gene_id}.tsv")

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
             --output "${gene_id}.tsv"
             --transcript-mapping stage/transcript.csv
             --chromosome-mapping stage/chrom.csv
             --threshold ${splice_threshold}
             --reference "${params.reference_genome}" )

      if [[ -n "${validation_log}" ]]; then
        ARGS+=( --log "${validation_log}" )
      fi

      python3 ${projectDir}/spliceai-parser.py "\${ARGS[@]}"

      echo "[parse_results] out.tsv lines: \$(wc -l < "${gene_id}.tsv") || true"
      """
  }
