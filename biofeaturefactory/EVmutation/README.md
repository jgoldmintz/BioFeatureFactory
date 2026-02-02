# EVmutation/PLMC Pipeline

## Overview

This pipeline computes evolutionary mutation effects using the epistatic model from protein multiple sequence alignments. It integrates **plmc** (Potts model maximum-likelihood inference) with **EVmutation** scoring to predict the fitness effects of amino acid substitutions.

---

## Requirements

| Component | Notes |
|-----------|-------|
| **EVmutation library** | Clone from Marks Lab (see setup below) |
| **plmc binary** | Compile from https://github.com/debbiemarkslab/plmc |
| Python ≥3.8 | With numpy, pandas, numba |

### EVmutation Library Setup

```bash
cd EVmutation/
git clone https://github.com/debbiemarkslab/EVmutation.git EVmutation
```

The pipeline imports directly from the `EVmutation/` subdirectory.

**Citation:** Hopf TA, et al. "Mutation effects predicted from sequence co-variation." *Nature Biotechnology* 35:128-135 (2017).

---

## Theoretical Background

### Potts Model

EVmutation fits a maximum-entropy model to an MSA:

$P(A_1, ..., A_L) = (1/Z) \cdot exp(\sum_i h_i(A_i) + \sum_{i<j} J_{ij}(A_i, A_j))$

Where:
- $h_i(A_i)$ = single-site fields (conservation)
- $J_{ij}(A_i, A_j)$ = pairwise couplings (coevolution)
- $Z$ = partition function

### Mutation Effect Prediction

The predicted effect of mutation $A_i \to A'_i$:

$\Delta E = h_i(A'_i) - h_i(A_i) + \sum_{j \neq i} [J_{ij}(A'_i, A_j) - J_{ij}(A_i, A_j)]$

**Interpretation**:
- $\Delta E < 0$: Mutation predicted deleterious
- $\Delta E > 0$: Mutation predicted beneficial/neutral
- More negative = more deleterious

---

## Workflow

```
        Protein MSA (FASTA)
              │
              v
    ┌─────────────────────┐
    │     plmc binary     │
    │  (Potts inference)  │
    └─────────────────────┘
              │
              v
       Model Parameters
       (.model_params)
              │
              v
    ┌─────────────────────┐
    │   CouplingsModel    │
    │  (load parameters)  │
    └─────────────────────┘
              │
              v
    ┌─────────────────────┐
    │  single_mutant_     │
    │  matrix()           │
    └─────────────────────┘
              │
              v
      EVmutation Scores
```

---

## Input Requirements

### Protein MSA (Required)
- FASTA format
- High-quality alignment (MUSCLE, MAFFT, HHblits, jackhmmer)
- Sufficient depth (> 100 sequences recommended)
- Focus sequence identified

### plmc Binary (Required unless skipping)
- Compiled plmc executable
- Specify via `--plmc-binary PATH`
- Source: https://github.com/debbiemarkslab/plmc

---

## Output Format

### Main Output (`{gene}.evmutation.tsv`)

| Column | Description | Range |
|--------|-------------|-------|
| `pkey` | Unique identifier (`GENE-mutation`) | string |
| `mutant` | Amino acid mutation string | e.g., `G25A` |
| `pos` | Position in protein (1-based) | integer |
| `wt` | Wildtype amino acid | single letter |
| `subs` | Substituted amino acid | single letter |
| `prediction_epistatic` | EVmutation epistatic score | float (typically -10 to +2) |
| `prediction_independent` | Single-site (field-only) score | float |
| `frequency` | Frequency of substitution at position | 0-1 |
| `column_conservation` | Max single-AA frequency at position | 0-1 |
| `qc_flags` | Quality control flags | string |

### Matrix Output (`{gene}.evmutation.matrix.tsv`)

Position × amino acid matrix of epistatic scores:

```
position	A	C	D	E	F	G	...
25	-4.80	-6.07	-3.95	-5.56	-7.28	0.00	...
26	-2.34	-3.12	-1.89	-2.45	-4.67	-1.23	...
```

### Evolutionary Couplings (`{gene}.ec`)

plmc automatically outputs pairwise coupling scores (not processed further by this pipeline).

---

## Usage

### Full Pipeline (Run plmc)
```bash
python evmutation-pipeline.py \
    --msa /path/to/protein.msa.fasta \
    --focus HUMAN_BRCA1 \
    --plmc-binary /path/to/plmc \
    --output BRCA1.evmutation.tsv
```

### Using Pre-computed Parameters
```bash
python evmutation-pipeline.py \
    --model-params /path/to/BRCA1.model_params \
    --output BRCA1.evmutation.tsv
```

### Map to Nucleotide Mutations
```bash
python evmutation-pipeline.py \
    --msa /path/to/protein.msa.fasta \
    --focus HUMAN_BRCA1 \
    --plmc-binary /path/to/plmc \
    --fasta /path/to/BRCA1.fasta \
    --mutations /path/to/mutations.csv \
    --output BRCA1.evmutation.tsv
```

### With Matrix Output
```bash
python evmutation-pipeline.py \
    --msa /path/to/protein.msa.fasta \
    --focus HUMAN_BRCA1 \
    --plmc-binary /path/to/plmc \
    --output BRCA1.evmutation.tsv \
    --output-matrix BRCA1.evmutation.matrix.tsv
```

---

## Parameters

### Required
| Parameter | Description |
|-----------|-------------|
| `--msa` OR `--model-params` | Protein MSA or pre-computed parameters |
| `--focus` | Focus sequence ID (if using MSA) |
| `--output` | Output TSV path |

### plmc Options
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--plmc-binary` | Required | Path to plmc executable |
| `--alphabet` | `-ACDEFGHIKLMNPQRSTVWY` | Amino acid alphabet |
| `--lambda-e` | 16.2 | J_ij regularization strength |
| `--lambda-h` | 0.01 | h_i regularization strength |
| `--skip-plmc` | False | Use existing model params |

### Mutation Mapping (Optional)
| Parameter | Description |
|-----------|-------------|
| `--fasta` | ORF FASTA for NT->AA mapping |
| `--mutations` | Mutations CSV file |
| `--validation-log` | Validation log for filtering |

---

## Score Interpretation

### Epistatic Score Distribution

| Score Range | Interpretation |
|-------------|----------------|
| > 0 | Likely neutral/beneficial |
| 0 to -2 | Mild effect |
| -2 to -5 | Moderate deleterious |
| < -5 | Strongly deleterious |

### Epistatic vs Independent Scores

| Score Type | Captures | Use Case |
|------------|----------|----------|
| `prediction_epistatic` | Conservation + coevolution | Primary prediction |
| `prediction_independent` | Conservation only | Baseline comparison |

**Difference** = epistatic contribution:
- Large difference -> strong coevolution signal
- Small difference -> effect driven by conservation

### Column Conservation

- High conservation (> 0.8): Strong constraint, most mutations deleterious
- Low conservation (< 0.3): Position tolerates variation

---

## plmc Requirements

### Compilation
```bash
git clone https://github.com/debbiemarkslab/plmc.git
cd plmc/src
make
# Binary at: plmc/bin/plmc
```

### System Requirements
- Multi-threaded (uses OpenMP)
- Memory scales with MSA depth and length
- Typical runtime: 10-60 minutes per gene

### plmc Parameters (Advanced)

| Flag | Default | Description |
|------|---------|-------------|
| `-le` | 16.2 | L2 regularization on J_ij |
| `-lh` | 0.01 | L2 regularization on h_i |
| `-m` | 500 | Maximum iterations |
| `-t` | 0.2 | Step size |
| `-g` | - | Enable reweighting |

---

## MSA Quality Guidelines

### Minimum Requirements
- ≥100 effective sequences
- ≤40% gaps per sequence
- ≤60% gaps per column
- Coverage of full protein length

### Optimal MSA
- 1,000-10,000 sequences
- Diverse taxonomic sampling
- Focus sequence high quality
- Generated with HHblits or jackhmmer

### MSA Sources
- **Pfam**: Pre-computed domain alignments
- **UniRef90/50**: Sequence clustering databases
- **Custom**: jackhmmer search against UniRef

---

## QC Flags

| Flag | Meaning | Action |
|------|---------|--------|
| `PASS` | Successful prediction | Use for analysis |
| `NOT_IN_MODEL` | Position not in MSA | Check MSA coverage |
| `INVALID_MUTATION` | Parse error | Check input format |

---

## Limitations

1. **Protein-level only**: Operates on amino acid sequences
2. **MSA quality dependent**: Poor MSA -> unreliable predictions
3. **Single mutations**: Does not handle complex variants
4. **Coverage gaps**: No prediction for positions outside MSA

---

## References

- Hopf TA, et al. (2017) Mutation effects predicted from sequence co-variation. **Nature Biotechnology**, 35:128-135.
- Figliuzzi M, et al. (2016) Coevolutionary landscape inference. **Mol Biol Evol**, 33:268-280.
- Ekeberg M, et al. (2013) Improved contact prediction in proteins. **Phys Rev E**, 87:012707.
- plmc GitHub: https://github.com/debbiemarkslab/plmc
- EVmutation GitHub: https://github.com/debbiemarkslab/EVmutation

---

## Copyright Notice

The `EVmutation` library (`model.py` and `tools.py`) is from the EVmutation package by Thomas A. Hopf (Marks Lab, Harvard).
- Source: https://github.com/debbiemarkslab/EVmutation
- plmc: https://github.com/debbiemarkslab/plmc

---

## License

This pipeline wrapper is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.
