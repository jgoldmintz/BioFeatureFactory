# Docker-based NetNGlyc

This Docker setup provides accurate NetNGlyc predictions with SignalP 6.0 integration for analyzing N-glycosylation sites in disease-associated synonymous mutations.

## Apple Silicon Solution

**Primary Use Case**: This Docker solution is specifically designed to solve NetNGlyc compatibility issues on **Apple Silicon (M1/M2/M3) Macs**. 

- **Problem**: NetNGlyc's 32-bit Linux binaries cannot run natively on ARM64 macOS
- **Solution**: Docker provides a 32-bit Linux environment where NetNGlyc works perfectly
- **Cross-platform**: Also works on Intel Macs and Linux systems for consistency

## Required Files and Naming Conventions

### Essential Files

#### Shared Docker Infrastructure (`BioFeatureFactory/docker/`)
- **`Dockerfile`** - Container definition with 32-bit Linux environment
- **`build-container.sh`** - Automated build script with ARM64 Mac support  
- **`signalp_stub`** - SignalP 6.0 integration script
- **`netNglyc.tar.gz`** - Original NetNGlyc 1.0 distribution (download from DTU)

#### NetNGlyc-Specific Files (`BioFeatureFactory/netnglyc/`)
- **`full-docker-netnglyc-pipeline.py`** - Complete production pipeline

### FASTA File Naming Requirements

#### Wildtype FASTA Files
- **Directory**: Must be in separate directory (e.g. `wt/aaseq/`)
- **Naming**: `{GENE}_aa.fasta` (exactly this format)
- **Content**: Single reference sequence per file
- **Examples**: `ABCB1_aa.fasta`, `BRCA1_aa.fasta`, `CFTR_aa.fasta`

#### Mutant FASTA Files  
- **Directory**: Must be in separate directory (e.g. `mut/aaseq/`)
- **Naming**: `{GENE}_aa.fasta` (same name as wildtype, different directory)
- **Content**: Multiple mutation sequences per file
- **Header format**: `>{GENE}-{MUTATION_ID}`

### Mapping File Requirements

#### Mapping Directory Structure
- **Directory**: Must be in separate mapping directory (e.g., `mutations/combined/aa/`)
- **Naming**: `{GENE}_nt_to_aa_mapping.csv` (exactly this format)
- **One file per gene**: Each gene must have its own mapping file

#### CSV Format Requirements
```csv
mutant,aamutant
T3435C,I1145I
A1967G,N656D
G2677T,A893S
```

## Quick Start

### 1. Obtain Required Software

#### NetNGlyc 1.0
Download `netNglyc.tar.gz` from:
- **Official source**: https://services.healthtech.dtu.dk/software.php
- **License**: Academic use only - requires institutional email registration
- **Place in**: `BioFeatureFactory/docker/` directory (same location as Dockerfile)

#### SignalP 6.0 (Required for signalp_stub)
-  Download and install SignalP 6.0 from:
- **Official source**: https://services.healthtech.dtu.dk/software.php
- **License**: Academic use only - requires institutional email registration

- **Install location**: Standard system path (e.g., `/usr/local/bin/signalp6` or conda environment)
- **Note**: The `signalp_stub` script calls SignalP 6.0 internally for signal peptide predictions

### 2. Build Docker Container
```bash
cd path/to/BioFeatureFactory/docker/
./build-container.sh
```

### 3. Run Complete Pipeline
```bash
# Process both wildtype and mutant sequences
python3 full-docker-netnglyc-pipeline.py \
    --fasta_wt path/to/wt/AminoAcid/Sequences/ \
    --fasta_mut path/to/mut/AminoAcid/Sequences/ \
    --output_wt path/to/wt/Outputfile \
    --output_mut path/to/mut/Outputfile \
    --mapping_dir path/to/mapping/directory/ \
    --workers 4
```

## Directory Structure Requirements

**_Critical:_** Files Must Be In Separate Directories

```
#Example file structure

project_root/
├── FASTA_files/
│   ├── wt/
│   │   └── aaseq/                    # Wildtype sequences ONLY
│   │       ├── ABCB1_aa.fasta       # Single reference sequence
│   │       ├── BRCA1_aa.fasta       # Single reference sequence  
│   │       └── CFTR_aa.fasta        # Single reference sequence
│   └── mut/
│       └── aaseq/                    # Mutant sequences ONLY
│           ├── ABCB1_aa.fasta       # Multiple mutation sequences
│           ├── BRCA1_aa.fasta       # Multiple mutation sequences
│           └── CFTR_aa.fasta        # Multiple mutation sequences
└── mutations/
    └── combined/
        └── aa/                       # Mapping files ONLY
            ├── ABCB1_nt_to_aa_mapping.csv
            ├── BRCA1_nt_to_aa_mapping.csv
            └── CFTR_nt_to_aa_mapping.csv
```

### File Naming Rules

1. **FASTA files**: Must end with `_aa.fasta`
2. **Mapping files**: Must end with `_nt_to_aa_mapping.csv`
3. **Gene names**: Must match exactly between FASTA and mapping files
4. **Case sensitive**: All names are case sensitive

## Configuration Options

### Adjust Worker Count
```bash
# Faster processing (more RAM usage)
python3 full-docker-netnglyc-pipeline.py --mode full-pipeline --workers 8 \
    path/to/wt/AminoAcid/Sequences/GENE_aa.fasta results.tsv

# Slower but less resource intensive
python3 full-docker-netnglyc-pipeline.py --mode full-pipeline --workers 2 \
    path/to/wt/AminoAcid/Sequences/GENE_aa.fasta results.tsv
```

## Usage Examples

### 1. Test Mode (No Other Args Required)
```bash
python3 full-docker-netnglyc-pipeline.py --test
```

### 2. Process Single FASTA File (Automatic Mode Selection)
```bash
python3 full-docker-netnglyc-pipeline.py \
    path/to/wt/AminoAcid/Sequences/GENE_aa.fasta \
    output-netnglyc.out
```

### 3. Full Pipeline - Wildtype Sequences (Requires Mapping)
```bash
python3 full-docker-netnglyc-pipeline.py --mode full-pipeline \
    --mapping-dir path/to/mapping/directory/ \
    path/to/wt/AminoAcid/Sequences/GENE_aa.fasta \
    results.tsv
```

### 4. Full Pipeline - Mutant Sequences (Requires Mapping)
```bash
python3 full-docker-netnglyc-pipeline.py --mode full-pipeline \
    --is-mutant \
    --mapping-dir path/to/mapping/directory/ \
    path/to/mut/AminoAcid/Sequences/GENE_aa.fasta \
    results.tsv
```

### 5. Force Parallel Processing with Custom Workers
```bash
python3 full-docker-netnglyc-pipeline.py --mode full-pipeline \
    --processing-mode parallel --workers 8 \
    path/to/wt/AminoAcid/Sequences/GENE_aa.fasta \
    results.tsv
```

### 6. Parse Existing NetNGlyc Outputs (Mutant Mode - Requires Mapping)
```bash
python3 full-docker-netnglyc-pipeline.py --mode parse \
    --is-mutant \
    --mapping-dir path/to/mapping/directory/ \
    netnglyc-outputs/ \
    parsed_results.tsv
```

### 7. Parse Existing NetNGlyc Outputs (Wildtype Mode - Requires Mapping)
```bash
python3 full-docker-netnglyc-pipeline.py --mode parse \
    --mapping-dir path/to/mapping/directory/ \
    netnglyc-outputs/ \
    parsed_results.tsv
```

## Processing Modes

### Automatic Mode Selection (Default)
- **auto**: Intelligent selection based on sequence count
  - 1-50 sequences → single (one Docker container)
  - 51-500 sequences → parallel (multiple Docker containers)
  - 501+ sequences → batch (sequential batching)

### Manual Mode Override
- **single**: One Docker container (optimal for small datasets)
- **parallel**: Multiple Docker containers (fast for medium datasets)
- **batch**: Sequential batching (required for very large datasets)
- **sequential**: One-by-one processing (fallback option)

## Required Arguments by Mode

| Mode | Input | Output | Additional Requirements |
|------|-------|--------|------------------------|
| `--mode parse` (mutant) | NetNGlyc output directory | TSV file | `--is-mutant`, `--mapping-dir` |
| `--mode parse` (wildtype) | NetNGlyc output directory | TSV file | `--mapping-dir` |
| `--mode full-pipeline` (mutant) | FASTA file/directory | TSV file | `--is-mutant`, `--mapping-dir` |
| `--mode full-pipeline` (wildtype) | FASTA file/directory | TSV file | `--mapping-dir` |

## Prerequisites

1. **Docker Desktop**: Install from https://www.docker.com/products/docker-desktop
2. **NetNGlyc License**: Academic registration required at DTU website  
3. **SignalP 6.0 License**: Academic registration required at DTU website
4. **4-8GB RAM**: For parallel processing (adjustable)
5. **Storage**: ~2GB for container + temporary processing files

## Troubleshooting

### File Naming Issues
```bash
# Check FASTA file naming
ls path/to/wt/AminoAcid/Sequences/    # Should show: GENE_aa.fasta
ls path/to/mut/AminoAcid/Sequences/    # Should show: GENE_aa.fasta

# Check mapping file naming
ls path/to/mapping/directory/    # Should show: GENE_nt_to_aa_mapping.csv
```

### Directory Structure Issues
```bash
# Verify separate directories
# WT and MUT files must be in different directories
# Mapping files must be in their own directory
```

### FASTA Header Format Issues
```bash
# Check mutant FASTA headers
head -5 path/to/mut/AminoAcid/Sequences/ABCB1_aa.fasta
# Should show: >ABCB1-T3435C format
# NOT: >T3435C or >ABCB1_T3435C
```

### Docker not found
```bash
# Install Docker Desktop and ensure it's running
docker --version
```

### Container build fails
```bash
# Check NetNGlyc source exists
ls -la path/to/netNglyc-1.0/
```

## Recent Bug Fixes (September-08-2025)

### Critical 1:Many Mapping Logic Fix

**Issue Resolved**: The NetNGlyc pipeline previously only checked position when matching predictions to mutations, missing the amino acid verification step.

**Problem Example**:
- Position 541: Wildtype K → mutations A1622G (K→E) and A1624C (K→N)
- **Before Fix**: Both mutations incorrectly received the same wildtype K prediction
- **After Fix**: Each mutation only matches predictions with the correct amino acid (E for A1622G, N for A1624C)

**Technical Changes**:
1. **Enhanced Amino Acid Parsing**: Added `get_mutation_data_bioAccurate_aa()` function to parse amino acid mutations
2. **Wildtype Processing**: Implemented 1:many mapping logic - one prediction can generate multiple pkeys for different mutations at the same position
3. **Mutant Processing**: Added amino acid verification using `pred['sequon'][0]` (first letter of sequon)
4. **Position + Amino Acid Matching**: Both wildtype and mutant processing now verify:
   - `pred['position'] == target_position`
   - `pred['sequon'][0] == target_amino_acid`

**Impact**:
- **Eliminates false positive matches**: Mutations only match predictions with correct amino acids
- **Supports synonymous mutations**: Multiple mutations at same position get distinct, accurate results
- **Improves data integrity**: Ensures biological accuracy in glycosylation predictions

**Files Modified**:
- `full-docker-netnglyc-pipeline.py`: Enhanced mapping logic and amino acid verification

This fix ensures that NetNGlyc predictions are correctly matched to mutations based on both genomic position and amino acid context, providing accurate glycosylation analysis for disease-associated synonymous mutations.

## Citation

If you use this pipeline, please cite:

- NetNGlyc: Gupta R, Jung E, Brunak S. Prediction of N-glycosylation sites in human proteins. NetNglyc 1.0 Server. 2004. Available from: https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
- SignalP: Teufel F, Almagro Armenteros JJ, Johansen AR, Gíslason MH, Pihl SI, Tsirigos KD, Winther O, Brunak S, von Heijne G, Nielsen H. SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol. 2022;40(7):1023-1025.
- BioFeatureFactory: Goldmintz J. BioFeatureFactory: Python toolkit for automated biological feature extraction. GitHub repository. https://github.com/jgoldmintz/BioFeatureFactory

## License

This project is licensed under the MIT License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues