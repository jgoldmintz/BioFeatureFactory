# Miranda Pipeline

This package contains a complete pipeline for predicting miRNA binding sites in wild type and mutant sequences using the miranda algorithm. 

## Contents

This self-contained package includes:

- **59 gene datasets**
- **Pre-compiled miranda executable** ready to use
- **Uncompiled miranda source** in `miranda.tar.gz` for custom compilation
- **Complete mutation annotation data** for transcript-level analysis (i.e. combined\__Gene-Name_\_mutations.csv)
- **Optimized pipeline script** with parallel processing and cleanup features
- **Conda environment** specification for easy setup

## Quick Start

### 1. Environment Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate miranda-env
```

### 2. Running the Pipeline

For mutant sequence analysis (most common use case):

```bash
python run-miranda-pipeline.py \
    --mode full-pipeline \
    -i data/fasta-files/mut/transcript/ \
    -o miranda-output/ \
    -m software/miranda-compiled/src/ \
    -d data/mature_hsa.fasta \
    -M
```

For wildtype sequence analysis:

```bash
python run-miranda-pipeline.py \
    --mode full-pipeline \
    -i data/fasta-files/mut/transcript/ \
    -o miranda-output/ \
    -m software/miranda-compiled/src/ \
    -d /path/to/mirna/database.fasta
```

## Software Components

### miranda Executable

- **Pre-compiled version**: Located in `software/miranda-compiled/src/miranda`
- **Source code**: Available in `miranda.tar.gz` for custom compilation

If you need to recompile miranda on your system:

```bash
tar -xzf miranda.tar.gz
cd miranda/
./configure
make
```

### Pipeline Features

The `run-miranda-pipeline.py` script includes several optimizations:

- **Parallel processing**: Uses all available CPU cores by default
- **Automatic cleanup**: Option to delete raw miranda files after parsing (for large datasets >1000 sequences)
- **Resume capability**: Skips already processed sequences
- **Multiple modes**: Full pipeline, miranda-only, or parse-only
- **Error handling**: Timeout protection and comprehensive logging

## Data Structure

### Input Data

- **FASTA files**: `data/fasta-files/mut/transcript/` - 59 genes with mutation sequences
- **Mutation data**: `data/mutations/combined/transcript/` - CSV files with mutation annotations

### Output Data

- **Raw miranda files**: Individual `.out` files for each sequence (optional cleanup)
- **Parsed results**: `filtered_mut_miranda.tsv` or `filtered_wt_miranda.tsv`

## Pipeline Modes

### 1. Full Pipeline (`--mode full-pipeline`)
Runs miranda analysis and parses results into TSV format.

### 2. Miranda Only (`--mode miranda-only`)
Runs only the miranda analysis, skipping parsing.

### 3. Parse Only (`--mode parse-only`)
Parses existing miranda output files into TSV format.

## Command Line Options

### Required Arguments
- `-i, --input`: Input directory containing FASTA files
- `-o, --output`: Output directory for miranda results
- `-m, --miranda_dir`: Path to directory containing miranda executable
- `-d, --mirna_db`: Path to miRNA database FASTA file

### Processing Options
- `-M, --is_mutant`: Process mutant sequences (default: wildtype)
- `--no-parallel`: Disable parallel processing
- `--max-workers N`: Set maximum number of parallel workers

### Advanced Options
- `--mutation_data_dir`: Path to mutation data directory (auto-detected by default)

## Example Workflows

### Process Single Gene
```bash
# Create subset directory with one gene
mkdir single-gene-test
cp data/fasta-files/mut/transcript/BRCA1_all_muts_nt.fasta single-gene-test/

# Run pipeline
python run-miranda-pipeline.py \
    --mode full-pipeline \
    -i single-gene-test/ \
    -o results/ \
    -m software/miranda-compiled/src/ \
    -d data/mature_hsa.fasta \
    -M
```

### Process with Custom Settings
```bash
# Run with 4 workers and keep raw files
python run-miranda-pipeline.py \
    --mode full-pipeline \
    -i data/fasta-files/mut/transcript/ \
    -o miranda-output/ \
    -m software/miranda-compiled/src/ \
    -d data/mature_hsa.fasta \
    -M \
    --max-workers 4
```

## Output Format

The parsed TSV file contains the following columns:
- `Pkey`: Mutation identifier (e.g., "BRCA1-C123T")
- `miRNA`: miRNA identifier (e.g., "hsa-miR-21-5p")
- `Transcript-Position`: Position in transcript where binding occurs
- `Tot Score`: Total miranda alignment score
- `Tot Energy`: Total binding energy (kcal/mol)
- `Max Score`: Maximum alignment score in the interaction
- `Max Energy`: Maximum binding energy (kcal/mol)
- `Strand`: Strand information and position details
- `Len1`: Length of miRNA sequence
- `Len2`: Length of target sequence
- `query-sequence`: miRNA sequence
- `ref-sequence`: Target mRNA sequence with alignment gaps

## Interpretation

Gaps indicate a mismatch between miRNA and target sequence and vice versa. The number of gaps is inversely related to scoring metrics (`Tot score`, `Tot Energy`, `Max Score`, and `Max Energy`). Nucleotide case also reflects binding confidence: lowercase indicates less confidence, while uppercase indicates higher confidence.


## Troubleshooting

### Common Issues

1. **"Cannot open file" error**: Ensure miRNA database path is absolute or relative to current directory
2. **Permission denied**: Check that miranda executable has execute permissions (`chmod +x software/miranda-compiled/src/miranda`)
3. **Out of memory**: Reduce `--max-workers` for systems with limited RAM

