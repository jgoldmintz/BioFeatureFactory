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

// Pre-built plmc params: when provided, evmutation_pipeline.py skips plmc and
// scores the mutations directly. File or directory of {GENE}.(codon_)model_params.
params.model_params       = null
params.codon_model_params = null

params.manifest           = null    // JSON manifest from controller

// Backend gating — defaults match the controller (both backends run).
// Skip flags are emitted as the string "true" by the controller only when
// disabling; never as "false" (Groovy's `if ("false")` is truthy).
params.skip_evmutation          = false   // when true: EVmutation/plmc processes skipped
params.skip_adabmdca            = false   // when true: adabmDCA processes skipped
params.skip_codon_evmutation    = false   // when true: codon-side EVmutation skipped
params.skip_codon_adabmdca      = false   // when true: codon-side adabmDCA skipped

// adabmDCA backend tunables (consulted when adabmDCA runs).
// The `adabmDCA` CLI is the Python/torch console script installed by
// `pip install adabmDCA` — it is expected on PATH and is not configurable here.
params.adabmdca_protein_params     = null          // pre-built protein params: file or directory
params.adabmdca_codon_params       = null          // pre-built codon params:   file or directory
// Alphabet selection is owned by adabmdca_pipeline.py: protein side uses the
// 21-char default ("-ACDEFGHIKLMNPQRSTVWY"); codon side uses the 65-char
// alphabet from bin/codon_encoding.py after encoding.
params.adabmdca_model              = 'bmDCA'       // training routine: bmDCA / eaDCA / edDCA
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
// ---- resource bounds for the adabmDCA processes -------------------------
// These processes had NO resource directive of any kind: not cpus, memory,
// accelerator, maxForks or label. With one gene per task and nothing declared,
// Nextflow launches as many concurrent adabmDCA tasks as there are genes, each
// one training on the GPU and then materializing a dense (L, L, q, q) coupling
// tensor on the host. adabmdca_pipeline.load_adabmdca_params documents that
// footprint: 14.5 GB for a codon-alphabet PAM, 53 GB for BRCA1. Two of those at
// once is an OOM, not a slowdown, and on a rented box it is a metered one.
//
// maxForks is the directive that matters on a single-GPU machine -- it is what
// turns "both tasks die" into "the second one waits". memory and cpus feed the
// local executor's admission control; accelerator is honoured by the grid/cloud
// executors and ignored (harmlessly) by the local one.
params.adabmdca_maxforks           = 1             // concurrent adabmDCA tasks
params.adabmdca_cpus               = 4
params.adabmdca_memory             = '32 GB'       // raise for codon-alphabet BRCA1
params.adabmdca_accelerator        = 1             // GPUs per task
// Number of GPUs to SPREAD tasks across. adabmdca_device is 'cuda' with no index
// (adabmdca_pipeline.py:149/:364 take one device string), so without this every
// concurrent task binds cuda:0 and they OOM each other -- raising
// adabmdca_maxforks on a multi-GPU box made things worse, not faster. Each task
// instead gets CUDA_VISIBLE_DEVICES set to one card, so 'cuda' resolves to a
// different physical GPU per task while the device string stays index-free.
// Leave at 1 for a single-GPU box; set to the card count and raise
// adabmdca_maxforks to match to actually use them.
params.adabmdca_gpus               = 1
// ---- codon-side capacity tier -------------------------------------------
// adabmDCA holds the whole (L, q, L, q) coupling tensor on ONE device: the
// clone has no DataParallel / DistributedDataParallel / torch.distributed /
// device_map anywhere, so a second GPU adds throughput, never capacity. Its
// sparse modes do not help either -- eaDCA/edDCA mask a still-dense tensor
// (train.py:198 `mask = torch.ones(size=(L, q, L, q), dtype=torch.bool)`,
// graph.py:214 `params["coupling_matrix"] *= mask`), so selecting one COSTS an
// extra L^2*q^2 bytes rather than saving any.
//
// L here is the codon MSA's ALIGNMENT WIDTH in codons, not the protein length,
// and the gap makes that much larger: F9's protein is 461 aa but its codon MSA
// is 2292 nt wide = 764 codon sites, which is (764/461)^2 = 2.75x the memory a
// protein-length estimate predicts. Measured from F9/CodonMSA/F9.codon.msa.fasta.
// With q = 65 (len(CODON_ALPHABET)) that puts F9 itself at ~23 GB -- already at
// the edge of a 24 GB card, so this tier is live for genes you run today, not
// only for BRCA1.
//
// A gene over the bound is routed to the plmc/EVmutation codon path instead,
// which solves the same Potts model in HOST RAM. That is the whole point: 53 GB
// of system memory is an ordinary box, 53 GB of VRAM is an H200.
// 0 = DETECT from nvidia-smi at start-up (see resolveCodonVram). A positive
// value is taken verbatim and skips detection, which is what you want on a
// scheduler where the head node has no GPU but the execution nodes do.
params.codon_dca_vram_gb           = 0
params.codon_dca_vram_headroom     = 0.9           // fraction of card total treated as usable
// Peak VRAM as a multiple of ONE coupling tensor. NOT a guess and NOT a knob the
// caller should have to tune: it is the count of simultaneously-live (L, q, L, q)
// tensors in adabmDCA's own bmDCA loop, read off the installed clone.
//
//   1.00  params["coupling_matrix"]            utils.py:133        float32
//   1.00  params_prev clone, INSIDE the loop   training.py:252     float32
//   1.00  fij_target (data two-point stats)    stats.py:69         float32
//   1.00  pij (model two-point stats)          stats.py:69         float32
//   1.00  grad["coupling_matrix"] = fij - pij  training.py:65      float32
//   0.25  mask (dense even for bmDCA)          train.py:198        bool
//   ----
//   5.25
//
// The previous default here was 2.5, which I picked with nothing behind it. It
// was low by 2.1x in the DANGEROUS direction: the tier would admit genes that
// then OOM partway through training, on metered hardware. F9 codon (L=764) sits
// at 23.0 GiB under 2.5 and 48.3 GiB under 5.25 -- "fits a 24 GB card" becomes
// "does not", which is the correct answer.
//
// Still approximate, and knowingly so: it counts coupling-shaped tensors only.
// It excludes the PCD chains (nchains x L, small), torch allocator
// fragmentation, and the extra transients the sparse models allocate
// (graph.py:138 exp_J for eaDCA/edDCA). Treat it as a floor. Override with
// --codon_dca_mem_factor once measured with torch.cuda.max_memory_allocated().
params.codon_dca_mem_factor        = 5.25
params.codon_dca_max_sites         = 0             // >0 overrides the VRAM calc with a hard L cap

// Required checks live inside the workflow block (Nextflow 26+ forbids
// top-level statements outside process / workflow / function bodies).
// See workflow { validateRequiredParams(); ... }.

def validateRequiredParams() {
    if (!params.fasta)        error "ERROR: --fasta is required"
    if (!params.mutations)    error "ERROR: --mutations is required"
    // plmc is only needed when params have to be built. Providing --model_params
    // (protein) / --codon_model_params (codon) skips plmc for that side, so the
    // binary is required only when a side still has to run inference.
    def needProteinPlmc = !params.skip_evmutation && !params.model_params
    def needCodonPlmc   = !params.skip_evmutation && !params.skip_codon_evmutation && !params.codon_model_params
    if ((needProteinPlmc || needCodonPlmc) && !params.plmc_binary)
        error "ERROR: --plmc_binary is required when the EVmutation backend builds params " +
              "(provide --model_params / --codon_model_params to skip plmc, or --skip_evmutation true)"
    if (!params.msa && !params.uniref90_db)
        error "ERROR: --uniref90_db is required when not providing --msa"
    if (!params.codon_msa && !params.db_root)
        error "ERROR: --db_root is required when not providing --codon_msa"
}

// ---------------- HELPERS ----------------
def resolveInputPath(pathParam) {
    def f = new File(pathParam as String)
    return f.isDirectory() ? f : (f.isFile() ? f : null)
}

def resolveMutationCsv(gene_id) {
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

// Alignment width of a codon MSA, in CODON SITES. Reads only up to the second
// header, so cost is one record rather than the whole alignment. Returns 0 when
// the file cannot be read, which codonDcaFitsGpu treats as "unknown -> let it
// run", so an unreadable MSA fails in the process with a real error instead of
// being silently rerouted.
def codonAlignmentSites(msa_file) {
    // No `while ((line = reader.readLine()) != null)`: Nextflow 26's parser
    // rejects an assignment inside a condition --
    //     Main: 188: Unexpected input: '=' @ line 188, column 26.
    // and it rejects bare `for` loops in this file too, so the whole record is
    // sliced out with findIndexOf and summed with each(). readLines() holds the
    // alignment in memory; these are plain text of a few MB (F9's is 21,240
    // lines) and it runs once per gene, so that is cheaper than being clever.
    def lines = []
    try {
        lines = new File(msa_file.toString()).readLines()
    } catch (Exception e) {
        return 0
    }
    int first = lines.findIndexOf { it.trim().startsWith('>') }
    if (first < 0) return 0
    int second = lines.findIndexOf(first + 1) { it.trim().startsWith('>') }
    if (second < 0) second = lines.size()
    int width = 0
    lines.subList(first + 1, second).each { width += it.trim().length() }
    return ((width / 3) as int)
}

// True when a codon-alphabet adabmDCA run for `sites` positions fits the
// declared VRAM budget. q = 65 = len(CODON_ALPHABET) in codon_encoding.py.
def codonDcaFitsGpu(int sites, double vram_gb) {
    if (sites <= 0) return true
    int hard_cap = params.codon_dca_max_sites.toInteger()
    if (hard_cap > 0) return sites <= hard_cap
    return codonDcaGiB(sites) <= vram_gb
}

def codonDcaGiB(int sites) {
    double q = 65.0d
    return ((sites as double) * (sites as double) * q * q * 4.0d
            * (params.codon_dca_mem_factor as double)) / (1024.0d * 1024.0d * 1024.0d)
}

// Total VRAM of the SMALLEST visible GPU, in GB, or null when it cannot be read.
// Smallest, not largest: a task can land on any card, so the tier has to be safe
// for the worst one. Total rather than free -- free at planning time says nothing
// about free when the task actually runs, and headroom covers the CUDA context.
def detectVramGb() {
    try {
        def proc = ['bash', '-c',
                    'nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null'].execute()
        proc.waitForOrKill(10000)
        def out = proc.text.trim()
        if (!out) return null
        def mibs = out.readLines().collect { it.trim() }.findAll { it ==~ /\d+/ }.collect { it as long }
        if (!mibs) return null
        return ((mibs.min() as double) / 1024.0d)
    } catch (Exception e) {
        return null
    }
}

// RETURNS the budget in GB. It does not write it back to params: Nextflow 26
// makes params effectively write-once, and a second assignment is DISCARDED with
// only a warning --
//     WARN: `params.a` is defined multiple times -- Assignments following the
//           first are ignored
// so `params.codon_dca_vram_gb = <detected>` silently left the value at 0 and the
// tier would have rejected every gene. The budget is resolved once in the
// workflow body and passed to codonDcaFitsGpu as an argument instead.
def resolveCodonVram() {
    if ((params.codon_dca_vram_gb as double) > 0) {
        println "[TIER] VRAM budget ${params.codon_dca_vram_gb} GB (declared; detection skipped)"
        return (params.codon_dca_vram_gb as double)
    }
    // Only needed when the codon adabmDCA tier is actually consulted. A run that
    // skips that side never touches a GPU and must not be blocked by the absence
    // of one -- the normal case on a dev machine.
    if (params.skip_adabmdca || params.skip_codon_adabmdca) {
        println "[TIER] codon adabmDCA skipped; no VRAM budget needed"
        return 0.0d
    }
    if (params.codon_dca_max_sites.toInteger() > 0) {
        println "[TIER] site cap ${params.codon_dca_max_sites} in force; VRAM not consulted"
        return 0.0d
    }
    def detected = detectVramGb()
    if (detected == null) {
        error "ERROR: cannot determine GPU VRAM for the codon adabmDCA tier.\n" +
              "  nvidia-smi is unavailable on this host and there is no default budget:\n" +
              "  a guessed value either wastes the card or OOMs a paid run.\n" +
              "  Pick one:\n" +
              "    --codon_dca_vram_gb <GB>     state it explicitly (use this when the\n" +
              "                                 head node has no GPU but the nodes do)\n" +
              "    --codon_dca_max_sites <L>    bypass the VRAM maths with a hard site cap\n" +
              "    --skip_codon_adabmdca true   run the codon side through plmc/EVmutation only"
    }
    double usable = detected * (params.codon_dca_vram_headroom as double)
    println "[TIER] detected ${String.format('%.1f', detected)} GB on the smallest visible GPU; " +
            "budget ${String.format('%.1f', usable)} GB at ${params.codon_dca_vram_headroom} headroom"
    return usable
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

// ---------------- WORKFLOW ----------------
workflow {
    validateRequiredParams()
    def codonVramGb = resolveCodonVram()

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

    // Step 3-4: EVmutation backend — gated on !skip_evmutation.
    if (!params.skip_evmutation) {
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

        // Step 4: Codon EVmutation. Runs for every gene as before, PLUS it is
        // the destination for genes the codon adabmDCA tier rejects. When
        // --skip_codon_evmutation is set it no longer skips unconditionally: an
        // oversized gene still runs here, because skipping it would leave that
        // gene with no codon-side result at all. Previously the two backends were
        // wholly independent -- a gene too large for the GPU simply produced a
        // failed adabmDCA task, and whether anything covered it was an accident
        // of which skip flags happened to be set.
        if (!params.skip_codon_evmutation || !params.skip_adabmdca) {
            def ev_overflow_only = params.skip_codon_evmutation
            def has_codon_tsv = (manifest.codon_EVmutation ?: []) as Set
            def codon_ev_input = fasta_ch
                .filter { gene_id, f -> !(gene_id in has_codon_tsv) }
                .join(codon_msa_ch)
                .filter { gene_id, fasta_file, codon_msa ->
                    if (!ev_overflow_only) return true
                    int sites = codonAlignmentSites(codon_msa)
                    boolean fits = codonDcaFitsGpu(sites, codonVramGb)
                    if (!fits) {
                        println "[TIER] ${gene_id}: codon L=${sites} needs " +
                                "${String.format('%.1f', codonDcaGiB(sites))} GiB > " +
                                "VRAM budget ${String.format('%.1f', codonVramGb)} GB; running codon " +
                                "EVmutation (plmc, host RAM) despite --skip_codon_evmutation"
                    }
                    return !fits
                }
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
    }

    // Step 5: adabmDCA backend — mirrors EVmutation's two-process structure.
    // Each process calls adabmdca_pipeline.py end-to-end (train + score);
    // Nextflow only stages files and gates by manifest, the way EVmutation does.
    if (!params.skip_adabmdca) {
        // -- Protein adabmDCA --
        def has_adabm_protein_tsv = (manifest.adabmdca_protein ?: []) as Set
        def protein_adabm_input = fasta_ch
            .filter { gene_id, f -> !(gene_id in has_adabm_protein_tsv) }
            .join(protein_msa_ch)
            .map { gene_id, fasta_file, protein_msa ->
                def mut_csv = resolveMutationCsv(gene_id)
                tuple(gene_id, fasta_file, protein_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'))
            }
            .filter { gene_id, fasta_file, protein_msa, mut_csv ->
                if (mut_csv.name == 'NO_MUTATIONS') {
                    println "WARNING: No mutation CSV found for ${gene_id}, skipping protein adabmDCA"
                    return false
                }
                return true
            }
        run_protein_adabmdca(protein_adabm_input)

        // -- Codon adabmDCA --
        if (!params.skip_codon_adabmdca) {
            def has_adabm_codon_tsv = (manifest.adabmdca_codon ?: []) as Set
            def codon_adabm_input = fasta_ch
                .filter { gene_id, f -> !(gene_id in has_adabm_codon_tsv) }
                .join(codon_msa_ch)
                .map { gene_id, fasta_file, codon_msa ->
                    def mut_csv = resolveMutationCsv(gene_id)
                    tuple(gene_id, fasta_file, codon_msa, mut_csv ? file(mut_csv) : file('NO_MUTATIONS'))
                }
                .filter { gene_id, fasta_file, codon_msa, mut_csv ->
                    if (mut_csv.name == 'NO_MUTATIONS') {
                        println "WARNING: No mutation CSV found for ${gene_id}, skipping codon adabmDCA"
                        return false
                    }
                    // Capacity gate. Submitting a gene that cannot fit the device
                    // buys nothing: the task OOMs partway through training, which
                    // on a rented box is a paid failure. Rejected genes are picked
                    // up by codon EVmutation in Step 4.
                    int sites = codonAlignmentSites(codon_msa)
                    if (!codonDcaFitsGpu(sites, codonVramGb)) {
                        println "[TIER] ${gene_id}: codon L=${sites} needs " +
                                "${String.format('%.1f', codonDcaGiB(sites))} GiB > " +
                                "VRAM budget ${String.format('%.1f', codonVramGb)} GB; " +
                                "routed to codon EVmutation (plmc) instead of adabmDCA"
                        return false
                    }
                    return true
                }
            run_codon_adabmdca(codon_adabm_input)
        }
    }
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
    python3 -m biofeaturefactory.core.msa_generation_pipeline \\
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
    python3 -m biofeaturefactory.core.codon_msa_pipeline \\
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
    publishDir { "${params.output_dir}/${gene_id}/EVmutation" }, mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/model_params", mode: 'copy', pattern: '*.model_params'
    tag { "${gene_id} protein" }

    input:
    tuple val(gene_id), path(fasta_file), path(protein_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.protein.tsv"),
          path("${gene_id}.model_params"), optional: true

    script:
    def logArg     = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def skipArg    = params.skip_codon_evmutation ? "--skip-codon" : ""
    def plmcArg    = params.plmc_binary ? "--plmc-binary \"${params.plmc_binary}\"" : ""
    def mparamsArg = params.model_params ? "--model-params \"${params.model_params}\"" : ""
    """
    set -euo pipefail
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
    """
}

process run_codon_evmutation {
    publishDir { "${params.output_dir}/${gene_id}/EVmutation" }, mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/codon_model_params", mode: 'copy', pattern: '*.codon_model_params'
    tag { "${gene_id} codon" }

    input:
    tuple val(gene_id), path(fasta_file), path(codon_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.codon.tsv"),
          path("${gene_id}.codon_model_params"), optional: true

    script:
    def logArg      = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def plmcArg     = params.plmc_binary ? "--plmc-binary \"${params.plmc_binary}\"" : ""
    def cmparamsArg = params.codon_model_params ? "--codon-model-params \"${params.codon_model_params}\"" : ""
    """
    set -euo pipefail
    python3 ${projectDir}/../evmutation_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --codon-msa "${codon_msa}" \\
        ${cmparamsArg} \\
        ${plmcArg} \\
        --output . \\
        ${logArg}

    if [ -d "${gene_id}/EVmutation" ]; then
        mv ${gene_id}/EVmutation/${gene_id}.codon.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.codon_model_params" -exec mv {} . \\; 2>/dev/null || true
    """
}

// ---------------- adabmDCA PROCESSES ----------------
// One process per side, mirroring run_protein_evmutation / run_codon_evmutation.
// Each invokes adabmdca_pipeline.py end-to-end (train + score). Nextflow only
// stages the inputs and publishes the TSV + params files.

process run_protein_adabmdca {
    publishDir { "${params.output_dir}/${gene_id}/adabmDCA" }, mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/adabmdca_protein_params", mode: 'copy', pattern: '*.protein_adabm_params'
    tag { "${gene_id} adabm protein" }
    maxForks    params.adabmdca_maxforks.toInteger()
    cpus        params.adabmdca_cpus.toInteger()
    memory      params.adabmdca_memory
    accelerator params.adabmdca_accelerator.toInteger()

    input:
    tuple val(gene_id), path(fasta_file), path(protein_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.protein.tsv"),
          path("${gene_id}.protein_adabm_params"), optional: true

    script:
    def logArg    = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def skipArg   = params.skip_codon_adabmdca ? "--skip-codon" : ""
    def paramsArg = params.adabmdca_protein_params ? "--protein-params \"${params.adabmdca_protein_params}\"" : ""
    def nepochsArg = params.adabmdca_nepochs ? "--adabmdca-nepochs ${params.adabmdca_nepochs}" : ""
    """
    set -euo pipefail

    # One GPU per task. See params.adabmdca_gpus.
    GPU_COUNT=${params.adabmdca_gpus}
    if [[ "\$GPU_COUNT" -gt 1 ]]; then
      export CUDA_VISIBLE_DEVICES=\$(( (${task.index} - 1) % \$GPU_COUNT ))
      echo "[adabmdca] task ${task.index} -> CUDA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES"
    fi
    python3 ${projectDir}/../adabmdca_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --msa "${protein_msa}" \\
        ${paramsArg} \\
        --output . \\
        --adabmdca-model ${params.adabmdca_model} \\
        ${nepochsArg} \\
        --adabmdca-tol ${params.adabmdca_tol} \\
        --adabmdca-patience ${params.adabmdca_patience} \\
        --adabmdca-check-every ${params.adabmdca_check_every} \\
        --adabmdca-target ${params.adabmdca_target} \\
        --adabmdca-lr ${params.adabmdca_lr} \\
        --adabmdca-nchains ${params.adabmdca_nchains} \\
        --adabmdca-nsweeps ${params.adabmdca_nsweeps} \\
        --adabmdca-device ${params.adabmdca_device} \\
        --adabmdca-dtype ${params.adabmdca_dtype} \\
        --adabmdca-seed ${params.adabmdca_seed} \\
        ${skipArg} \\
        ${logArg}

    if [ -d "${gene_id}/adabmDCA" ]; then
        mv ${gene_id}/adabmDCA/${gene_id}.protein.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.protein_adabm_params" -exec mv {} . \\; 2>/dev/null || true
    """
}

process run_codon_adabmdca {
    publishDir { "${params.output_dir}/${gene_id}/adabmDCA" }, mode: 'copy', pattern: '*.tsv'
    publishDir "${params.output_dir}/adabmdca_codon_params", mode: 'copy', pattern: '*.codon_adabm_params'
    tag { "${gene_id} adabm codon" }
    maxForks    params.adabmdca_maxforks.toInteger()
    cpus        params.adabmdca_cpus.toInteger()
    memory      params.adabmdca_memory
    accelerator params.adabmdca_accelerator.toInteger()

    input:
    tuple val(gene_id), path(fasta_file), path(codon_msa), path(mutations_csv)

    output:
    tuple val(gene_id),
          path("${gene_id}.codon.tsv"),
          path("${gene_id}.codon_adabm_params"), optional: true

    script:
    def logArg    = params.validation_log ? "--validation-log \"${file(params.validation_log)}\"" : ""
    def paramsArg = params.adabmdca_codon_params ? "--codon-params \"${params.adabmdca_codon_params}\"" : ""
    def nepochsArg = params.adabmdca_nepochs ? "--adabmdca-nepochs ${params.adabmdca_nepochs}" : ""
    """
    set -euo pipefail

    # One GPU per task. See params.adabmdca_gpus.
    GPU_COUNT=${params.adabmdca_gpus}
    if [[ "\$GPU_COUNT" -gt 1 ]]; then
      export CUDA_VISIBLE_DEVICES=\$(( (${task.index} - 1) % \$GPU_COUNT ))
      echo "[adabmdca] task ${task.index} -> CUDA_VISIBLE_DEVICES=\$CUDA_VISIBLE_DEVICES"
    fi

    python3 ${projectDir}/../adabmdca_pipeline.py \\
        --fasta "${fasta_file}" \\
        --mutations "${mutations_csv}" \\
        --codon-msa "${codon_msa}" \\
        ${paramsArg} \\
        --output . \\
        --adabmdca-model ${params.adabmdca_model} \\
        ${nepochsArg} \\
        --adabmdca-tol ${params.adabmdca_tol} \\
        --adabmdca-patience ${params.adabmdca_patience} \\
        --adabmdca-check-every ${params.adabmdca_check_every} \\
        --adabmdca-target ${params.adabmdca_target} \\
        --adabmdca-lr ${params.adabmdca_lr} \\
        --adabmdca-nchains ${params.adabmdca_nchains} \\
        --adabmdca-nsweeps ${params.adabmdca_nsweeps} \\
        --adabmdca-device ${params.adabmdca_device} \\
        --adabmdca-dtype ${params.adabmdca_dtype} \\
        --adabmdca-seed ${params.adabmdca_seed} \\
        ${logArg}

    if [ -d "${gene_id}/adabmDCA" ]; then
        mv ${gene_id}/adabmDCA/${gene_id}.codon.tsv . 2>/dev/null || true
    fi
    find . -name "${gene_id}.codon_adabm_params" -exec mv {} . \\; 2>/dev/null || true
    """
}
