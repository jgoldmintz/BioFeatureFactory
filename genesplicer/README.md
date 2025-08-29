

# GeneSplicer Pipeline

GeneSplicer pipeline for splice site prediction in genomic DNA sequences.

## Overview

This module provides a Python wrapper for GeneSplicer, a computational method for detecting splice sites (donor and acceptor sites) in eukaryotic genomic DNA. The pipeline supports multiple processing modes and includes optimizations for handling large-scale genomic data.

## Features

- **Multiple Processing Modes**:
  - `full`: Complete pipeline from sequence processing through prediction
  - `tool_only`: Run GeneSplicer predictions only
  - `parse_only`: Parse existing GeneSplicer output files
- **Optimized Performance**:
  - ThreadPoolExecutor for parallel processing
  - Intelligent batch processing for large datasets
  - Memory-efficient handling of genomic sequences
- **Robust File Management**:
  - Automatic temporary file cleanup
  - Configurable output directory structure
  - Resume capability for interrupted runs
- **Comprehensive Output**:
  - Splice site predictions with confidence scores
  - Donor and acceptor site coordinates
  - Confidence levels (Low, Medium, High)
  - Position-specific scoring

## Installation

### Prerequisites

1. **GeneSplicer binary**: Download and install GeneSplicer from the [official site](https://ccb.jhu.edu/software/genesplicer/)
1. **Training data**: Obtain species-specific training data (human, Arabidopsis, etc.)
1. **mafft**: Required for mapping genomic positions. 

### Setup

```bash
# Clone the repository
git clone https://github.com/jgoldmintz/BioFeatureFactory.git
cd BioFeatureFactory/genesplicer

# Install Python dependencies
pip install -r requirements.txt

# Set GeneSplicer path
export GENESPLICER_PATH=/path/to/genesplicer/bin
export GENESPLICER_TRAINING=/path/to/training/data
```

## Usage

### Basic Usage

**Note**: Based on the current implementation, the pipeline uses the command-line script rather than a Python class interface. Here's the actual usage:

```bash
# Process FASTA files with GeneSplicer
python run_genesplicer_pipeline.py \
    --input /path/to/fasta/directory/ \
    --output /path/to/output/ \
    --genesplicer_dir /path/to/genesplicer/bin/ \
    --is_mutant
```

### Command Line Interface

```bash
# Process mutant sequences (default full pipeline)
python run_genesplicer_pipeline.py \
    --input /path/to/fasta/directory/ \
    --output /path/to/output/ \
    --genesplicer_dir /path/to/genesplicer/bin/ \
    --is_mutant

# Process wildtype sequences
python run_genesplicer_pipeline.py \
    --input /path/to/fasta/directory/ \
    --output /path/to/output/ \
    --genesplicer_dir /path/to/genesplicer/bin/

# Disable parallel processing
python run_genesplicer_pipeline.py \
    --input /path/to/fasta/directory/ \
    --output /path/to/output/ \
    --genesplicer_dir /path/to/genesplicer/bin/ \
    --no-parallel
```

### Processing Modes

#### Full Pipeline Mode

Processes raw sequences through the complete pipeline:

1. Sequence validation and formatting
1. GeneSplicer execution with species-specific parameters
1. Output parsing and feature extraction
1. Results aggregation and reporting

#### Tool Only Mode

Runs GeneSplicer on pre-processed sequences:

- Assumes sequences are already in the correct format
- Executes GeneSplicer with configured parameters
- Generates raw prediction output

#### Parse Only Mode

Processes existing GeneSplicer output files:

- Parses GeneSplicer prediction files
- Extracts splice site features
- Formats results for downstream analysis


## Output Format

### Splice Site Predictions

The output file contains tab-delimited splice site predictions:

| Gene Name   | Wt nt | Genomic Position | Mut nt | End5  | End3  | Score | confidence | splice_site_type |
|-------------|-------|------------------|--------|-------|-------|-------|------------|------------------|
| CFTR-G1211T | G     | 87857            | T      | 87856 | 87857 | 3.50  | Medium     | acceptor         |
| CD44-C618T  | C     | 50844            | T      | 50843 | 50844 | 4.84  | Medium     | acceptor         |


## Performance Optimization

### Memory Management

- Sequences processed in configurable batches
- Automatic cleanup of temporary files
- Efficient memory usage for large genomic datasets

### Parallel Processing

- ThreadPoolExecutor for concurrent sequence processing
- Configurable worker threads based on system resources
- Automatic load balancing for optimal performance

### Temporary File Management

- Automatic cleanup of intermediate files
- Configurable temp directory location
- Safe handling of interrupted processes

## Error Handling

The pipeline includes comprehensive error handling:

- Timeout protection for long-running predictions
- Graceful handling of malformed sequences
- Detailed error logging with stack traces
- Resume capability for failed batches

## Logging

The pipeline provides detailed logging:

- **Progress tracking**: Real-time updates on file processing
- **Error reporting**: Detailed error messages with sequence identifiers
- **Performance metrics**: Processing time and throughput statistics
- **Debug information**: Temporary file creation and cleanup status

```bash
# Example log output
Processing BRCA1_all_muts_nt.fasta...
  Processing 1219 mutations in parallel...
Processing CFTR_all_muts_nt.fasta...
  Processing 847 mutations in parallel...
```

## Examples

### Process Single Gene
```bash
# Create test directory with one gene
mkdir test_single_gene
cp /path/to/BRCA1_all_muts_nt.fasta test_single_gene/

# Run GeneSplicer on single gene
python run_genesplicer_pipeline.py \
    --input test_single_gene/ \
    --output results/ \
    --genesplicer_dir /usr/local/bin/genesplicer/ \
    --is_mutant
```

### Process with Custom Workers
```bash
# Use 4 parallel workers
python run_genesplicer_pipeline.py \
    --input /path/to/fasta/directory/ \
    --output /path/to/output/ \
    --genesplicer_dir /path/to/genesplicer/bin/ \
    --max-workers 4 \
    --is_mutant
```

## Troubleshooting

### Common Issues

1. **"GeneSplicer not found" error**:
   - Verify `--genesplicer_dir` points to directory containing `genesplicer` executable
   - Check executable permissions: `chmod +x /path/to/genesplicer`

2. **"No training data" error**:
   - Ensure `human/` or species-specific training directory exists in GeneSplicer installation
   - Verify training files are present and readable

3. **Memory/Performance issues**:
   - Use `--no-parallel` for systems with limited memory
   - Reduce `--max-workers` to lower resource usage
   - Process smaller batches of FASTA files

4. **Temporary file errors**:
   - Check write permissions in GeneSplicer directory
   - Ensure sufficient disk space for temporary files

### File Requirements

- **FASTA format**: Sequences must be in standard FASTA format
- **File naming**: Files should follow pattern `{GENE}_all_muts_nt.fasta`
- **Sequence headers**: Use format `>{GENE}-{MUTATION_ID}` for mutations

## Citation

If you use this pipeline, please cite:

- GeneSplicer: Pertea M, Lin X, Salzberg SL. GeneSplicer: a new computational method for splice site prediction. Nucleic Acids Res. 2001
- BioFeatureFactory: Goldmintz J. BioFeatureFactory: Python toolkit for automated biological feature extraction. GitHub repository. https://github.com/jgoldmintz/BioFeatureFactory

## License

This project is licensed under the MIT License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues