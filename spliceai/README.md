# SpliceAI Pipeline

This Nextflow-based pipeline provides automated SpliceAI splice site prediction analysis for disease-associated synonymous mutations, featuring complete end-to-end processing from mutation CSV files to parsed results with pkey mapping.

## Key Features

### Comprehensive Processing Pipeline
- **VCF Generation**: Convert snv CSV files to SpliceAI-compatible VCF format
- **Exon-Aware Mapping**: Coordinate system conversion between ORF, transcript, genomic, and chromosome space
- **Automated SpliceAI Execution**: Parallel splice prediction with TensorFlow race condition mitigation
- **Intelligent Result Parsing**: Extract predictions with mutation-specific pkey mapping and validation log filtering

### TensorFlow Stability Enhancements
- **Race Condition Prevention**: Environment variables to control TensorFlow threading
- **Controlled Parallelism**: Configurable concurrent process limits
- **Retry Logic**: Exponential backoff with jitter for failed predictions
- **Process Isolation**: Staggered execution to prevent model loading conflicts

### Flexible Input Handling
- **Multiple Input Modes**: Single VCF, VCF directory, or mutation CSV directory
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
- Provide the directory via `--mutations_dir`
- Each CSV lists mutation IDs (e.g., C123T) in a single column with a header row
- Naming is flexible; gene symbols are inferred from filenames

#### Mapping Files (Required for pkey generation)
- Supply mapping directories using `--chromosome_mapping_dir`, `--transcript_mapping_dir`, and `--genomic_mapping_dir`
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

### 3. Run Complete Pipeline
The commands below show typical invocations; replace the sample paths with your own locations.
```bash
nextflow run main.nf \
    --mutations_dir mutations/ \
    --reference_genome /path/to/reference.fna \
    --annotation_file annotation_ncbi.txt \
    --transcript_mapping_dir mutations/combined/transcript/ \
    --genomic_mapping_dir mutations/combined/genomic/ \
    --chromosome_mapping_dir mutations/combined/chromosome/ \
    --output_dir results/

nextflow run main.nf \
    --input_vcf_dir vcf_files/ \
    --skip_vcf_generation \
    --reference_genome /path/to/reference.fna \
    --annotation_file annotation_ncbi.txt \
    --transcript_mapping_dir mutations/combined/transcript/ \
    --genomic_mapping_dir mutations/combined/genomic/ \
    --output_dir results/

nextflow run main.nf \
    --input_vcf_file ABCB1.vcf \
    --reference_genome /path/to/reference.fna \
    --annotation_file annotation_ncbi.txt \
    --transcript_mapping_dir mutations/combined/transcript/ \
    --genomic_mapping_dir mutations/combined/genomic/ \
    --output_dir results/
```

## Processing Modes

### Input Mode Selection
- **`--mutations_dir`**: Process mutation CSV files (generates VCFs first)
- **`--input_vcf_dir`**: Process pre-existing VCF directory
- **`--input_vcf_file`**: Process single VCF file
- **`--skip_vcf_generation`**: Skip VCF generation when using `--input_vcf_dir`

### VCF Generation Options
- **`--chromosome_format`**: Output chromosome format (`refseq`, `simple`, `ucsc`)
- **`--validate_mapping`**: Cross-check chromosome mappings against reference
- **`--clear_vcf_cache`**: Force regeneration of cached VCF files

### SpliceAI Execution Control
- **`--splice_threshold`**: Minimum delta score threshold (default: 0.0)
- **`--retry_jitter`**: Maximum retry delay in seconds (default: 10)
- **`--maxforks`**: Limit concurrent SpliceAI tasks (default: 0 = unbounded, capped by CPU count)

### Validation and Filtering
- **`--validation_log`**: Path to validation log for filtering failed mutations
- **`--chromosome_mapping_dir`**: Required for pkey generation
- **`--transcript_mapping_dir`**: Required for mutation ID mapping

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
- `block_label`: Unique identifier for transcript isoforms - _A, B, C, D,_ or _dup_ (where _dup_ indicates consistent output across isoforms).

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
nextflow run main.nf \
    --clear_vcf_cache \
    --retry_jitter 20 \
    [other parameters]
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
