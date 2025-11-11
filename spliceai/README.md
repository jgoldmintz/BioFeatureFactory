# SpliceAI Pipeline

This Nextflow-based pipeline provides automated SpliceAI splice site prediction analysis for disease-associated synonymous mutations, featuring complete end-to-end processing from mutation CSV files to parsed results with pkey mapping.

## Key Features

### Comprehensive Processing Pipeline
- **VCF Generation**: Convert snv CSV files to SpliceAI-compatible VCF format
- **Exon-Aware Mapping**: Coordinate system conversion between ORF, transcript, genomic, and chromosome space
- **Automated SpliceAI Execution**: Parallel splice prediction with TensorFlow race condition mitigation
- **Intelligent Result Parsing**: Extract predictions with mutation-specific pkey mapping and validation log filtering
- **Live Progress Tracking**: Optional monitor that reports per-gene SpliceAI progress in real time

### TensorFlow Stability Enhancements
- **Race Condition Prevention**: Environment variables to control TensorFlow threading
- **Controlled Parallelism**: Configurable concurrent process limits
- **Retry Logic**: Exponential backoff with jitter for failed predictions
- **Process Isolation**: Staggered execution to prevent model loading conflicts

### Flexible Input Handling
- **Multiple Input Modes**: Run SpliceAI on a single VCF (`--input_vcf_file`), a directory of VCFs (`--input_vcf_dir`), or auto-build VCFs from mutation specs (`--mutations_path`)
- **Chromosome Format Support**: RefSeq (NC_000007.14) and simple (7) chromosome naming
- **Validation Integration**: Skip mutations flagged in exon-aware validation logs
- **Cache Management**: Intelligent VCF caching for repeated runs

## Required Files and Setup

### Essential Pipeline Files
- **`main.nf`** - Nextflow pipeline with TensorFlow threading controls and multi-input support
- **`../dependencies/vcf_converter.py`** - Robust SNP to VCF converter with RefSeq support (packaged in dependencies but runnable standalone)
- **`spliceai-parser.py`** - Advanced VCF parser with pkey mapping and log filtering
- **`annot_to_spliceai.py`** - GTF to SpliceAI annotation format converter

### Input Requirements

#### Mutation CSV Files
- Provide the path via `--mutations_path` (directory or single file).
- Each CSV must be named `<gene_id>_mutations.csv` (case-insensitive) so the gene ID seeds downstream channels
- Directory mode scans every matching CSV and emits `<gene_id>.vcf` per input
- File mode reuses the single CSV for one gene

#### Mapping Files (Required for pkey generation)
- Supply mapping paths using `--transcript_mapping_path`, `--genomic_mapping_path`, and optional `--chromosome_mapping_path` (legacy aliases: `--*_dir`)
- Paths may point to a directory (per-gene CSVs) or a single CSV reused for every gene
- Directory mode discovery is case-insensitive and matches `*<GENE>*.csv`
- Transcript and genomic mapping are required per gene; missing files trigger an error
- Chromosome mapping is optional. If absent for a gene, the parser warns and continues without chromosome mapping for that gene
- Chromosome mapping example (`combined_{GENE}.csv`):
  ```csv
  mutant,chromosome
  C123T,C87504250T
  G456A,G87504456A
  ```
- Transcript mapping example (`combined_{GENE}.csv`):
  ```csv
  mutant,transcript
  C123T,I1145I
  G456A,N656D
  ```

#### Reference Files
- Specify paths through `--reference_genome` and `--annotation_file`
- Reference genome: Indexed FASTA (e.g., `GRCh38_reference.fna`)
- SpliceAI annotation: Tab-delimited file produced by `annot_to_spliceai.py`
- Optional: provide `--vcf_output_dir` to keep intermediate gzip-compressed VCFs separate from `--output_dir`

### Annotation File Generation

#### Convert GTF to SpliceAI Format
```bash
# Generate annotation with RefSeq chromosome naming (matches VCF format)
python3 annot_to_spliceai.py \
    path/to/annotation.gtf \
    --chromosome-format ncbi \
    -o spliceai_annotation.txt

# Alternative: Simple chromosome naming
python3 annot_to_spliceai.py \
    path/to/annotation.gtf \
    --chromosome-format simple \
    -o spliceai_annotation_simple.txt
```

## Quick Start

### 1. Prerequisites
```bash
# Ensure SpliceAI conda environment is available
conda activate ssnvs_prediction_tool
spliceai -h  # Verify SpliceAI installation

# Ensure Nextflow is installed
nextflow -version
```

### 2. Generate Annotation File
```bash
python3 annot_to_spliceai.py \
    /path/to/reference.gtf \
    --chromosome-format ncbi \
    -o annotation_ncbi.txt
```

### 3. Run the Pipeline
Use the controller as the single entry point. The commands below show typical invocations; replace the sample paths with your own locations.

```bash
python spliceai-pipeline-controller.py \
    --mutations_path mutations/ \
    --reference_genome /path/to/reference.fna \
    --annotation_file annotation_ncbi.txt \
    --transcript_mapping_path mappings/transcript/ \
    --chromosome_mapping_path mappings/chromosome/ \
    --genomic_mapping_path mappings/genomic/ \
    --output_dir results/ \
    --initial_maxforks 3
```

- To reuse existing VCFs (single file or directory) instead of generating them, supply `--input_vcf_path /path/to/vcfs --skip_vcf_generation` and omit `--mutations_path`.
- To generate new VCFs from mutation CSVs, provide `--mutations_path` and omit `--input_vcf_path`; the controller builds VCFs automatically before running SpliceAI.
- The controller launches Nextflow, streams logs, runs the progress tracker, and automatically drops to `--maxforks tail_maxforks` when the tail genes hit rapid TensorFlow exit‑134 errors.

## Processing Modes

### Input Mode Selection
- **`--mutations_path`**: Process mutation CSV files (directory mode or single CSV).
- **`--input_vcf_path`**: Process pre-existing VCFs. Accepts either a single `.vcf` or a directory of `.vcf` files; requires `--skip_vcf_generation`.
- **`--skip_vcf_generation`**: Signals that VCFs already exist and should be reused (requires `--input_vcf_path`).

### VCF Source Selection
- Provide `--mutations_path` when you want the pipeline to generate VCFs. Leave `--skip_vcf_generation` unset in this mode.
- Use `--skip_vcf_generation --input_vcf_path /path/to/vcfs` when VCFs already exist. If the path is a directory, every `*.vcf` is processed; if it is a file, only that gene is run.
- When `--input_vcf_path` is omitted, VCFs are built from `--mutations_path`.

### VCF Generation Options
- **`--chromosome_format`**: Output chromosome format (`refseq`, `simple`, `ucsc`)
- **`--validate_mapping`**: Cross-check chromosome mappings against reference
- **`--clear_vcf_cache`**: Force regeneration of cached VCF files
- **`--vcf_output_dir`**: Publish intermediate VCFs, bgzipped files, and indexes to a dedicated directory

### SpliceAI Execution Control
- **`--initial_maxforks`** *(controller)*: Max concurrent `run_spliceai` tasks for the first pass (default: 3).
- **`--tail_maxforks`** *(controller)*: Maxforks applied after the controller detects rapid tail retries (default: 1).
- **`--tail_threshold`**, **`--rapid_window_minutes`**, **`--rapid_gap_minutes`** *(controller)*: Heuristics that define “tail mode” and the frequency of exit‑134 failures that triggers a restart.
- **`--splice_threshold`** *(pipeline)*: Minimum delta score threshold passed to the parser (default: 0.0).
- **`--retry_jitter`** *(pipeline)*: Maximum retry delay for `run_spliceai` (default: 10).
- **`--maxforks`** *(pipeline)*: Final limit enforced inside Nextflow (the controller sets this to `initial_maxforks` or `tail_maxforks` depending on the phase).
- The progress tracker runs automatically; disable it with `--disable_tracker` if you do not want the live `processed/total` table.

### Validation and Filtering
- **`--validation_log`**: Path to validation log for filtering failed mutations
- **`--chromosome_mapping_path`**: Optional per-gene chromosome mapping (legacy `--chromosome_mapping_dir`)
- **`--transcript_mapping_path`**: Required for mutation ID mapping (legacy `--transcript_mapping_dir`)
- **`--genomic_mapping_path`**: Required for genomic mapping (legacy `--genomic_mapping_dir`)

## TensorFlow Race Condition Mitigation

### Problem
SpliceAI uses TensorFlow, which has cumulative metrics counters that can race when multiple processes load/unload models simultaneously, causing exit 134 (SIGABRT) crashes.

### Solution
The pipeline automatically applies TensorFlow threading controls:
```bash
export TF_NUM_INTEROP_THREADS=1
export TF_NUM_INTRAOP_THREADS=1
export OMP_NUM_THREADS=1
export TF_CPP_MIN_LOG_LEVEL=3
```

### Process Management
- **Staggered Execution**: 5-20 second delays prevent simultaneous model loading
- **Retry Logic**: Exponential backoff with jitter for failed processes
- **Error Recovery**: Automatic retry up to 3 attempts per sample

## Output Format

### SpliceAI Results TSV
```csv
pkey	gene	chrom	pos	ref	alt	allele	block_label	ds_ag	ds_al	ds_dg	ds_dl	dp_ag	dp_al	dp_dg	dp_dl	max_delta_score
ABCB1-C123T	ABCB1	NC_000007.14	87504250	C	T	C	A	0.05	0.12	0.03	0.08	-15	25	-8	12	0.12
ABCB1-G456A	ABCB1	NC_000007.14	87504456	G	A	G	B	0.15	0.02	0.25	0.04	5	-18	10	-22	0.25
```

### Key Fields
- `pkey`: Mutation identifier (`{GENE}-{MUTATION_ID}`)
- `ds_ag/ds_al/ds_dg/ds_dl`: Delta scores for acceptor gain/loss, donor gain/loss
- `dp_ag/dp_al/dp_dg/dp_dl`: Distance to predicted sites
- `max_delta_score`: Highest delta score for the variant
- `block_label`: SpliceAI evaluates variants across transcript isoforms and emits one “transcript block” per isoform in the VCF INFO field. The parser preserves this order, labelling the first four unique score vectors as `A`, `B`, `C`, and `D`. Additional isoforms—or isoforms whose predictions match an earlier block exactly—are labelled `dup`. This highlights isoform-specific splicing differences while collapsing redundant outputs. For example, in `Data/spliceai/AR_spliceai_results.tsv`, variant `AR-G1130A` appears twice: block `A` captures the first AR isoform, whereas block `dup` represents other isoforms with identical scores.

## Advanced Usage

### Custom Annotation Generation
```bash
# Generate annotation for different assemblies
python3 annot_to_spliceai.py reference.gtf --chromosome-format simple -o simple.txt
python3 annot_to_spliceai.py reference.gtf --chromosome-format ncbi -o ncbi.txt
python3 annot_to_spliceai.py reference.gtf --chromosome-format ucsc -o ucsc.txt
```

### VCF Conversion Only
```bash
python3 dependencies/vcf_converter.py \
    -m ABCB1_mutations.csv \
    -o vcf_output/ \
    --chromosome-format refseq \
    -r reference.fna \
    -a annotation.gtf \
    --chromosome-mapping-input chromosome_mapping/ \
    --validate-mapping
```

### Result Parsing Only
```bash
python3 spliceai-parser.py \
    --input ABCB1.spliceai.vcf \
    --output ABCB1_results.tsv \
    --chromosome-mapping chromosome_mapping/combined_ABCB1.csv \
    --transcript-mapping transcript_mapping/combined_ABCB1.csv \
    --threshold 0.1 \
    --log validation.log
```

### Live Progress Tracker
The controller launches `spliceai-tracker.py` automatically so you always see the live `processed/total` table. Disable it with `--disable_tracker` if you prefer a quiet log.

You can also run the tracker manually if you want to monitor a previously started run:

```bash
python3 spliceai-tracker.py \
  --mutations-path /path/to/mutations_dir \
  --work-root /path/to/work \
  --poll-seconds 5 \
  --log-file tracker.log
```

The tracker is observational only; stop it with `Ctrl+C` when running manually. The controller cleans it up automatically when the pipeline finishes or restarts.

## Troubleshooting

### Common Issues

#### Chromosome Naming Mismatch
```bash
# Error: No SpliceAI predictions generated
# Cause: VCF uses "NC_000007.14" but annotation uses "7"
# Solution: Regenerate annotation with matching format
python3 annot_to_spliceai.py annotation.gtf --chromosome-format ncbi -o annotation_fixed.txt
```

#### TensorFlow Crashes (Exit 134)
```bash
# The pipeline automatically applies threading controls
# If crashes persist, check TensorFlow version compatibility
conda list tensorflow
# SpliceAI works best with TensorFlow 1.x versions
# Additional mitigations:
#   - Increase `--retry_jitter` to spread retries further apart
#   - Set `--maxforks 1` to run SpliceAI tasks sequentially when parallel launches keep failing
```

#### Missing Mapping Files
```bash
# Error: No pkeys generated
# Cause: Missing mapping files
# Solution: Regenerate the mapping CSVs with the exon-aware builder
python3 ../dependencies/exon_aware_mapping.py \
    --mutations /path/to/mutations/ \
    --annotation /path/to/annotation.gtf \
    --reference /path/to/reference_genome.fna \
    --out-chromosome-mapping /path/to/chromosome_mappings/ \
    --out-transcript-mapping /path/to/transcript_mappings/ \
    --out-genomic-mapping /path/to/genomic_mappings/ \
    --out-fasta /path/to/output_fastas/
```

### Performance Optimization
```bash
# For large datasets, adjust cache and retry settings
python spliceai-pipeline-controller.py \
    --mutations_path /path/to/mutations \
    --reference_genome /path/to/reference.fna \
    --annotation_file annotation_ncbi.txt \
    --transcript_mapping_path mappings/transcript/ \
    --chromosome_mapping_path mappings/chromosome/ \
    --genomic_mapping_path mappings/genomic/ \
    --output_dir results/ \
    --clear_vcf_cache \
    --retry_jitter 20 \
    --initial_maxforks 4 \
    --tail_maxforks 1
```

## Prerequisites

1. **Nextflow**: Version 20.04+ with DSL2 support
2. **SpliceAI**: Installed in conda environment (`ssnvs_prediction_tool`)
3. **Reference Data**: Indexed reference genome and gene annotation
4. **Mapping Files**: Exon-aware coordinate mappings (generated by `exon_aware_mapping.py`)
5. **System Requirements**: 8GB+ RAM, multi-core CPU recommended

## Citation

If you use this SpliceAI pipeline, please cite:

- **SpliceAI**: Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24.
- **BioFeatureFactory**: Goldmintz J. BioFeatureFactory: Python toolkit for automated biological feature extraction. GitHub repository. https://github.com/jgoldmintz/BioFeatureFactory

## License

This project is licensed under the MIT License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
