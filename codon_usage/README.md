# Codon Usage Pipeline

## Overview

This pipeline computes codon and codon-pair usage statistics, quantifying translational efficiency through codon optimality metrics.

---

## Metrics Computed

### Single Codon Metrics

| Metric | Full Name | Description | Range |
|--------|-----------|-------------|-------|
| **RSCU** | Relative Synonymous Codon Usage | Observed frequency / expected frequency for synonymous codons | 0 to ~6 |
| **W** | Relative Adaptiveness | Codon frequency relative to most frequent synonymous codon (gene-specific) | 0 to 1 |
| **CAI_W** | CAI Reference W | Reference W value from human highly expressed genes | 0 to 1 |
| **tAI** | tRNA Adaptation Index | tRNA availability weight based on gene copy numbers | 0 to 1 |

### Gene-Level Metrics

| Metric | Full Name | Description | Range |
|--------|-----------|-------------|-------|
| **CAI_gene** | Codon Adaptation Index | Geometric mean of CAI_W values across gene | 0 to 1 |
| **tAI_gene** | tRNA Adaptation Index | Geometric mean of tAI weights across gene | 0 to 1 |

### Codon Pair Metrics

| Metric | Full Name | Description | Range |
|--------|-----------|-------------|-------|
| **RSCPU** | Relative Synonymous Codon Pair Usage | Observed pair frequency / expected | 0 to ~10 |
| **CPS** | Codon Pair Score | ln(observed/expected) for codon pairs | -∞ to +∞ |
| **noln CPS** | Non-log CPS | observed/expected ratio (no logarithm) | 0 to +∞ |
| **W_CP** | Relative Adaptiveness (Pair) | Pair frequency relative to most frequent pair | 0 to 1 |

---

## Biological Background

### Codon Usage Bias

Synonymous codons are not used equally. Highly expressed genes preferentially use "optimal" codons that:
- Match abundant tRNAs
- Enable faster translation elongation
- Reduce ribosome stalling

### Codon Pair Bias

Adjacent codon pairs also show non-random usage patterns:
- Some pairs are favored (CPS > 0)
- Some pairs are avoided (CPS < 0)
- Pair bias affects translation efficiency and co-translational folding

---

## Workflow

1. **Load ORF Sequence**
   Read nucleotide sequence from FASTA file.

2. **Compute Reference Statistics**
   Calculate genome-wide codon and codon-pair frequencies.

3. **Process Each Mutation**
   For each mutation:
   - Identify the mutated codon
   - Extract flanking codons (bicodons)
   - Look up usage metrics

4. **Output Results**
   Write TSV with per-mutation codon usage statistics.

---

## Output Format

### Main Output (`{gene}.codon_usage.tsv`)

| Column | Description | Example |
|--------|-------------|---------|
| `pkey` | Unique identifier (`GENE-mutation`) | `BRCA1-C123T` |
| `Gene` | Gene symbol | `BRCA1` |
| `codon_number` | Codon position in ORF (1-based) | `41` |
| `codon` | Mutated codon sequence | `CTG` |
| `position_in_codon` | Position of SNV within codon (1-3) | `3` |
| `RSCU` | Relative synonymous codon usage | `1.32` |
| `W` | Relative adaptiveness (gene-specific) | `0.87` |
| `CAI_W` | CAI reference W (human highly expressed) | `0.45` |
| `tAI` | tRNA adaptation weight for codon | `0.344` |
| `CAI_gene` | CAI for entire gene | `0.72` |
| `tAI_gene` | tAI for entire gene | `0.58` |
| `bicodon_3prime` | 3' bicodon (this codon + next) | `CTGGAA` |
| `RSCPU_3prime` | RSCPU for 3' bicodon | `0.95` |
| `CPS_3prime` | Codon pair score for 3' bicodon | `-0.12` |
| `noln_CPS_3prime` | Non-log CPS for 3' bicodon | `0.89` |
| `W_CP_3prime` | Relative adaptiveness for 3' pair | `0.76` |
| `bicodon_5prime` | 5' bicodon (previous codon + this) | `AAACTG` |
| `RSCPU_5prime` | RSCPU for 5' bicodon | `1.12` |
| `CPS_5prime` | Codon pair score for 5' bicodon | `0.08` |
| `noln_CPS_5prime` | Non-log CPS for 5' bicodon | `1.08` |
| `W_CP_5prime` | Relative adaptiveness for 5' pair | `0.82` |
| `bicodon_context` | Position context | `middle_codon_both_directions` |
| `qc_flags` | Quality flags | `PASS` |

---

## Usage

### Single File Processing
```bash
python codon-usage-pipeline.py \
    --fasta /path/to/BRCA1.fasta \
    --mutations /path/to/BRCA1_mutations.csv \
    --output BRCA1.codon_usage.tsv
```

### Directory Processing
```bash
python codon-usage-pipeline.py \
    --fasta-dir /path/to/fastas/ \
    --mutations-dir /path/to/mutations/ \
    --output combined.codon_usage.tsv
```

### Mutant FASTA Mode
For FASTA files where mutations are encoded in sequence names:
```bash
python codon-usage-pipeline.py \
    --fasta-dir /path/to/mutant_fastas/ \
    --is-mutant \
    --output results.tsv
```

---

## Metric Calculations

### RSCU (Relative Synonymous Codon Usage)

$RSCU_i = X_i / ((1/n) \sum_{j=1}^{n} X_j)$

Where:
- $X_i$ = observed count of codon $i$
- $n$ = number of synonymous codons for that amino acid

**Interpretation**:
- RSCU = 1.0: Used at expected frequency
- RSCU > 1.0: Overrepresented (preferred)
- RSCU < 1.0: Underrepresented (avoided)

### W (Relative Adaptiveness)

$W_i = X_i / max(X_j)$

Where $\max(X_j)$ is the count of the most frequent synonymous codon.

**Interpretation**:
- W = 1.0: Most optimal codon for that amino acid
- W < 1.0: Less optimal (proportional to usage)

### CAI (Codon Adaptation Index)

$CAI = \exp((1/L) \sum_{i=1}^{L} \ln(W_i))$

Where:
- $L$ = number of codons (excluding stops)
- $W_i$ = reference W value for codon $i$ (from highly expressed genes)

**Interpretation**:
- CAI ≈ 1.0: Gene uses optimal codons (highly expressed gene pattern)
- CAI < 0.5: Gene uses suboptimal codons
- Reference W values derived from human highly expressed genes (Sharp & Li 1987)

### tAI (tRNA Adaptation Index)

$tAI = \exp((1/L) \sum_{i=1}^{L} \ln(w_i))$

Where:
- $L$ = number of codons
- $w_i$ = tRNA adaptation weight for codon $i$ based on tRNA gene copy numbers

**Interpretation**:
- tAI ≈ 1.0: Codons match abundant tRNAs (fast translation)
- Low tAI: Codons require rare tRNAs (slow translation)
- Weights from dos Reis et al. 2004, Tuller et al. 2010

### CPS (Codon Pair Score)

$CPS_{ij} = \ln(f_{ij} / (f_i \cdot f_j))$

Where:
- $f_{ij}$ = observed frequency of codon pair
- $f_i, f_j$ = individual codon frequencies

**Interpretation**:
- CPS > 0: Pair is favored (occurs more than expected)
- CPS = 0: Pair occurs at expected frequency
- CPS < 0: Pair is avoided (occurs less than expected)

---

## Bicodon Context

The `bicodon_context` field indicates what codon pairs are available:

| Context | Description |
|---------|-------------|
| `first_codon_3prime_only` | Mutation in first codon (no 5' neighbor) |
| `last_codon_5prime_only` | Mutation in last codon (no 3' neighbor) |
| `middle_codon_both_directions` | Both 5' and 3' bicodons available |
| `insufficient_sequence` | Cannot extract bicodons (flagged) |

---

## QC Flags

| Flag | Meaning |
|------|---------|
| `PASS` | All metrics computed successfully |
| `NO_BICODON` | Could not extract codon pairs |
| `INVALID_CODON` | Codon contains non-standard bases |

---

## Reference Data

Codon usage statistics are calculated from:
- Human coding sequences (CDS) from RefSeq/Ensembl
- Weighted by expression levels (optional)
- Excludes first/last codons and rare genes

### Default Reference: Human Genome
- ~20,000 protein-coding genes
- ~40 million codons
- Standard genetic code

---

## Limitations

1. **Reference-dependent**: Metrics depend on reference codon usage table
2. **Context-limited**: Does not consider mRNA structure
3. **Single codon focus**: Analyzes one mutation at a time
4. **Species-specific**: Human codon usage may not apply to other organisms

---

## References

- Sharp PM, Li WH (1987) The codon adaptation index. **Nucleic Acids Res**, 15:1281-1295.
- dos Reis M, Savva R, Wernisch L (2004) Solving the riddle of codon usage preferences: a test for translational selection. **Nucleic Acids Res**, 32:5036-5044.
- Tuller T, et al. (2010) An evolutionarily conserved mechanism for controlling the efficiency of protein translation. **Cell**, 141:344-354.
- Plotkin JB, Kudla G (2011) Synonymous but not the same. **Nat Rev Genet**, 12:32-42.
- Coleman JR, et al. (2008) Virus attenuation by genome-scale changes in codon pair bias. **Science**, 320:1784-1787.
- Quax TE, et al. (2015) Codon bias as a means to fine-tune gene expression. **Mol Cell**, 59:149-161.

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.
