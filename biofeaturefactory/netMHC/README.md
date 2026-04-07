# NetMHC Pipeline

High-throughput MHC binding prediction for WT and mutant protein sequences using NetMHCpan 4.1 or other NetMHC tools.

---

## 1. Overview

NetMHC predicts binding of peptides to MHC (Major Histocompatibility Complex) molecules, which is crucial for understanding:
- T-cell epitope prediction
- Cancer neoantigen identification
- Vaccine design
- Immunogenicity assessment

This pipeline:
- Predicts MHC class I and/or class II binding for multiple HLA alleles
- Compares WT vs mutant sequences to identify gained/lost epitopes
- Generates ensemble TSV outputs with delta scores and classifications
- Supports batch processing for large datasets

---

## 2. Requirements

| Component | Notes |
|-----------|-------|
| NetMHCpan 4.1 or NetMHC 4.0 | Download from [DTU Health Tech](https://services.healthtech.dtu.dk/) (academic license). NetMHC 4.0 also requires [data files](https://services.healthtech.dtu.dk/services/NetMHC-4.0/data.tar.gz) extracted into the installation directory. Set path via `--native-netmhc-path` or `NETMHC_PATH`/`NETMHC_HOME` env var. |
| Python >=3.9 | Uses shared utilities from `../utils/`. |
| Mutation CSVs | Single-column NT mutations (header: `mutant`) or multi-column mapping CSVs. |
| WT FASTAs | ORF FASTAs (header `>ORF`) in a directory or single file. |
| Log files (optional) | Validation logs to filter known failures. |

---

## 3. Installation

### NetMHCpan 4.1

1. Download from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) (academic license)
2. Install per NetMHCpan documentation
3. Point to the binary via `--native-netmhc-path` (accepts the binary directly or its parent directory)

### NetMHC 4.0

1. Download from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) (academic license)
2. Download and extract the [data files](https://services.healthtech.dtu.dk/services/NetMHC-4.0/data.tar.gz) into the installation directory
3. Point to the platform binary via `--native-netmhc-path` (e.g., `netMHC-4.0/{os_arch}/bin/netMHC` or `netMHC-4.0/{os_arch}/`)

**Allele format note**: NetMHC 4.0 uses `HLA-A0201` format. NetMHCpan uses `HLA-A*02:01` format. The default allele matches the tool selected via `--netmhc-tool`.

---

## 4. Running

### Full pipeline (WT FASTA $\rightarrow$ Mutant synthesis $\rightarrow$ NetMHC $\rightarrow$ TSV ensemble)

```bash
python netmhc_pipeline.py \
    ../../FASTA_files/newnt3/ \
    results/ \
    -m ../../mutations/combined/ \
    --alleles HLA-A0201 HLA-A0101 HLA-B0702 \
    --log /path/to/validation.log \
    --keep-intermediates
# Writes per gene: results/{GENE}/NetMHC/{GENE}.{tsv,events.tsv,sites.tsv}
```

This flow:
1. Loads ORF nucleotide FASTAs and translates to amino acids
2. Synthesizes mutant AA sequences from mapping CSVs
3. Runs NetMHC on both WT and mutant sequences
4. Generates ensemble TSV with epitope gain/loss analysis

---

## 5. Output Files

### Output Structure

Output is written per gene to:

```
{output}/
  {GENE}/
    NetMHC/
      {GENE}.tsv          -- per-mutation summary
      {GENE}.events.tsv   -- per-peptide, per-allele events
      {GENE}.sites.tsv    -- raw binding predictions
```

### 5.1 Summary TSV (`{GENE}.tsv`)

Per-mutation summary with epitope changes:

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | `{GENE}-{MUTATION}` primary key | string |
| `Gene` | Gene symbol | string |
| `mutation` | Mutation identifier | string |
| `n_epitopes_wt`, `n_epitopes_mut` | Count of predicted epitopes (rank < threshold) | count |
| `count_gained`, `count_lost`, `count_strengthened`, `count_weakened` | Epitope classification counts | count |
| `max_abs_delta_rank`, `sum_abs_delta_rank` | Maximum and sum of rank changes | percentile (0-100) |
| `top_event_type` | Dominant event (gained/lost/strengthened/weakened) | categorical |
| `top_event_allele` | HLA allele for top event | string |
| `top_event_peptide` | Peptide sequence for top event | string |
| `top_event_delta_rank` | Rank change for top event | percentile |
| `qc_flags` | Quality control flags (missing_wt, missing_mut, no_delta) | string |

### 5.2 Events TSV (`{GENE}.events.tsv`)

Per-mutation, per-peptide, per-allele events:

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | `{GENE}-{MUTATION}` primary key | string |
| `Gene` | Gene symbol | string |
| `mutation` | Mutation identifier | string |
| `peptide` | Peptide sequence | string |
| `pos` | Position in protein | residue index |
| `mhc_allele` | HLA allele | string |
| `wt_rank`, `mut_rank`, `delta_rank` | WT rank, MUT rank, and $\Delta$ | percentile |
| `wt_affinity`, `mut_affinity` | Binding affinity (nM) | nanomolar |
| `delta_affinity` | Change in binding affinity ($\text{MUT} - \text{WT}$) | nanomolar |
| `bind_level_wt`, `bind_level_mut` | SB (strong) / WB (weak) / NB (non-binder) | categorical |
| `classification` | gained/lost/strengthened/weakened/stable | categorical |
| `classification_code` | Numeric encoding (gained=2, lost=-2, etc.) | integer |

### 5.3 Sites TSV (`{GENE}.sites.tsv`)

Raw predictions for all peptides:

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | `{GENE}-{MUTATION}` primary key | string |
| `Gene` | Gene symbol | string |
| `sequence_type` | wt or mut | categorical |
| `pos` | Position in protein | residue index |
| `mhc_allele` | HLA allele | string |
| `peptide` | Peptide sequence | string |
| `core` | Core binding region | string |
| `affinity` | Binding affinity | nM |
| `rank` | Percentile rank | percentile |
| `bind_level` | Strong binder (SB), Weak binder (WB), or non-binder | categorical |
| `identity` | Sequence name from NetMHC output | string |

---

## 6. MHC Allele Selection

### Common HLA Alleles

The pipeline supports any HLA alleles recognized by NetMHCpan. Common choices:

**Class I (most frequent in global populations):**
- HLA-A*02:01 (most common, ~50% of populations)
- HLA-A*01:01, HLA-A*03:01, HLA-A*24:02
- HLA-B*07:02, HLA-B*08:01, HLA-B*44:03
- HLA-C*07:02, HLA-C*07:01

**Class II (via NetMHCIIpan):**
- HLA-DRB1*01:01, HLA-DRB1*15:01
- HLA-DQA1*05:01-DQB1*02:01

### Specifying Alleles

```bash
# Single allele
--alleles HLA-A*02:01

# Multiple alleles
--alleles HLA-A*02:01 HLA-A*01:01 HLA-B*07:02
```

If `--alleles` is not specified, the pipeline defaults to `HLA-A0201` (NetMHC 4.0 format). For NetMHCpan, use `HLA-A*02:01` format. Specify multiple alleles explicitly for population coverage.

---

## 7. Interpretation

### Epitope Classifications

| Classification | Meaning | Clinical Relevance |
|---------------|---------|-------------------|
| **Gained** | New strong binder in mutant (WT was non-binder) | Potential neoantigen - may trigger immune response |
| **Lost** | Strong binder in WT becomes non-binder in mutant | Immune escape - tumor may evade detection |
| **Strengthened** | Both bind, but mutant binds stronger ($\Delta$ rank $> 5\%$) | Enhanced immunogenicity |
| **Weakened** | Both bind, but mutant binds weaker ($\Delta$ rank $> 5\%$) | Reduced immunogenicity |
| **Stable** | Both bind with minimal change ($\Delta$ rank $\leq 5\%$) | No significant immune impact |

### Binding Thresholds

NetMHC uses percentile rank to classify binders:
- **Strong Binder (SB)**: rank $\leq 0.5\%$
- **Weak Binder (WB)**: $0.5\% <$ rank $\leq 2\%$
- **Non-Binder**: rank $> 2\%$

---

## 8. Command-Line Options

### MHC-Specific Options
- `--alleles ALLELE [ALLELE ...]`: HLA alleles to predict
- `--netmhc-tool {netMHCpan, netMHC, netMHCII}`: Which NetMHC tool to use
- `--threshold FLOAT`: Rank threshold for strong binders (default: 0.5)

### Execution Backend
- `--native-netmhc-path PATH`: Path to native NetMHC executable

### Processing Options
- `--mutations`, `-m`: Mutation file or directory of mutation CSVs (required for full-pipeline)
- `--log FILE`: Validation log to skip failed mutations
- `--batch-size N`: Split sequences longer than N amino acids into batches (default: 100). Use 0 to disable batching.
- `--timeout SEC`: Command timeout (default: 600)

### Output
- `--keep-intermediates`: Keep temp files for debugging
- `--verbose`: Enable verbose output

---

## 9. Batch Processing

### When Batching is Used

The pipeline automatically batches sequences longer than the `--batch-size` threshold (default: 100 AA):

```bash
# Default: sequences >100 AA are batched
python netmhc_pipeline.py input/ results/

# Larger batches for better performance on shorter proteins
python netmhc_pipeline.py input/ results/ --batch-size 500

# Disable batching entirely (process full sequences)
python netmhc_pipeline.py input/ results/ --batch-size 0
```

### How Batching Works

1. **Sequence splitting**: Large proteins split into overlapping chunks
2. **Independent processing**: Each batch runs through NetMHC separately
3. **Result merging**: Predictions from all batches combined into single output
4. **Position preservation**: Peptide positions maintained relative to full-length protein

### Performance Considerations

| Sequence Length | Recommended Batch Size | Rationale |
|----------------|------------------------|-----------|
| < 200 AA | 0 (disabled) | Small enough to process whole |
| 200-1000 AA | 100-200 | Balance memory and overhead |
| 1000-3000 AA | 300-500 | Reduce batch count |
| > 3000 AA | 500-1000 | Minimize I/O overhead |

**Note**: Very small batch sizes increase overhead from repeated NetMHC startup time.

---

## 10. Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `NetMHC not found` | Install NetMHCpan and ensure binary is executable. Use `--native-netmhc-path` or set `NETMHCPAN_PATH`/`NETMHC_PATH`. |
| `No mapping file found` | Verify `-m` points to a mutation CSV or directory containing `{GENE}*.csv` files |
| `Peptide too short/long` | Verify input sequences are long enough for the selected NetMHC tool |
| `Invalid allele name` | NetMHC 4.0 uses `HLA-A0201` format; NetMHCpan uses `HLA-A*02:01` format. Check allele list with the tool's `-listMHC` flag. |
| Binary not executing | Ensure executable permissions are set on the NetMHC binary |

---

## 11. Implementation Status

### Completed Features

**All core functionality is now operational:**
- Pipeline structure and CLI interface
- WT/mutant sequence synthesis
- Mapping CSV integration
- Validation log filtering
- Native execution support (`_run_native_netmhc`)
- NetMHC output parsing (`parse_netmhc_output`)
- WT vs MUT comparison and epitope classification
- Three-tier TSV output system (summary, events, sites)
- Automatic batch processing for large sequences (>100 AA)

**Analysis Features:**
- Epitope gain/loss detection
- Binding affinity change quantification
- Per-mutation summary statistics
- Detailed event tracking with classification codes
- QC flag generation for quality control

### Known Limitations

1. **Default alleles**: When `--alleles` not specified, uses HLA-A*02:01 only. For population coverage, specify multiple alleles explicitly.
2. **Parallel processing**: Batches processed sequentially; could parallelize across genes/batches for large datasets.
3. **`--netmhc-tool`**: Selects the binary path but does not adjust command-line flags per tool.

---

## 12. Citations

- NetMHCpan 4.1: Reynisson, B. et al. (2020). "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation." *Nucleic Acids Research*, 48(W1), W449-W454.
- Jurtz, V. et al. (2017). "NetMHCpan-4.0: Improved Peptide-MHC Class I Interaction Predictions Integrating Eluted Ligand and Peptide Binding Affinity Data." *The Journal of Immunology*, 199(9), 3360-3368.

---

## 13. Example Workflow

### Step 1: Verify Data Organization

```bash
# 1. Check WT ORF FASTAs
ls ../../FASTA_files/newnt3/
# Expected: ABCB1_nt.fasta, BRCA1_nt.fasta, ...

# 2. Check mapping CSVs
ls ../../mutations/combined/aa/
# Expected: ABCB1_transcript_mapping.csv, BRCA1_transcript_mapping.csv, ...
```

### Step 2: Run Full Pipeline

```bash
# Create output directory
mkdir -p results

# Run with multiple HLA alleles for population coverage
python netmhc_pipeline.py \
    ../../FASTA_files/newnt3/ \
    results/ \
    -m ../../mutations/combined/ \
    --alleles HLA-A0201 HLA-A0101 HLA-B0702 \
    --threshold 0.5 \
    --verbose
```

### Step 3: Analyze Outputs

Three TSV files are generated per gene:

```bash
# Per-mutation summary with epitope counts and top events
head results/BRCA1/NetMHC/BRCA1.tsv

# Detailed epitope changes (gained/lost/strengthened/weakened)
head results/BRCA1/NetMHC/BRCA1.events.tsv

# Raw binding predictions for all peptides
head results/BRCA1/NetMHC/BRCA1.sites.tsv
```

**Interpretation Example:**
- `count_gained=3`: Mutation creates 3 new strong-binding epitopes (potential neoantigens)
- `count_lost=1`: Mutation eliminates 1 WT epitope (potential immune escape)
- `top_event_delta_rank=-15.2`: Strongest change is 15.2% improvement in binding

---

**Note**: This pipeline integrates with the existing BioFeatureFactory ecosystem and follows the same patterns as NetNglyc and NetPhos pipelines.

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.
