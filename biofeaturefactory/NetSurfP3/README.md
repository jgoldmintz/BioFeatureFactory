# NetSurfP-3.0 Pipeline

High-throughput protein structure prediction for WT and mutant sequences using NetSurfP-3.0. Predicts surface accessibility, secondary structure, and disorder regions.

---

## 1. Overview

NetSurfP-3.0 predicts multiple structural features per residue:
- **Surface Accessibility** (RSA/ASA): Solvent exposure
- **Secondary Structure** (SS3/SS8): Helix, strand, coil assignments
- **Disorder**: Intrinsically disordered regions (two disorder scores: `disorder_pf` and `disorder_pt`)

This pipeline:
- Compares WT vs mutant structural changes at the mutation site and globally
- Identifies buried-to-exposed transitions, disorder changes, and secondary structure shifts
- Generates ensemble TSV outputs with delta scores and classifications
- Supports batch processing for large datasets

---

## 2. Training a model

The pipeline requires a trained model checkpoint (`--model`). To train one, use the nsp3 CLI from within the cloned repo:

```bash
cd NetSurfP3/nsp3/nsp3
nsp3 train -c experiments/netsurfp_3/<CONFIG_DIR>/<CONFIG>.yml
```

Checkpoints are saved to the `save_dir` specified in the config YAML (e.g., `saved/nsp3/CNNbLSTM/`). 

### 2.1 Choosing the right `--config` for your trained model

The `--config` flag must match the architecture used to train the `--model` checkpoint. Config YAMLs live at `nsp3/experiments/netsurfp_3/`:

| Config | `arch.type` | Loss | Output heads | Use when |
|--------|------------|------|-------------|----------|
| `CNNbLSTM/CNNbLSTM.yml` | `CNNbLSTM_ESM1b_Complete` | `multi_task_loss` | SS8, SS3, disorder, RSA, phi, psi | Full structural prediction (recommended) |
| `CNNbLSTM_ESM1b_v2/config.yml` | `CNNbLSTM_Extended` | `multi_task_extended` | SS8, disorder, RSA iso/cpx, phi, psi | Extended variant (n_hidden=64, no explicit SS3 head) |
| `CNNbLSTM_ESM1b/CNNbLSTM_ESM1b.yml` | `CNNbLSTM_ESM1b` | `secondary_structure_loss` | SS8, SS3 only | SS-only (no disorder/RSA/angles — incompatible with this pipeline) |
| `ESM1b/ESM1b.yml` | `ESM1b` | `secondary_structure_loss` | SS8, SS3 only | SS-only, different model class |
| `CNNTrans/CNNTrans.yml` | `CNNTrans` | varies | varies | Transformer variant |

This pipeline expects all 6 output tensors (SS8, SS3, disorder, RSA, phi, psi). Use `CNNbLSTM/CNNbLSTM.yml` unless you trained with a different config.

To verify which config a checkpoint was trained with:

```python
import torch
ckpt = torch.load("your_checkpoint.pth", map_location="cpu", weights_only=False)
print(ckpt["config"]["arch"]["type"])  # e.g., "CNNbLSTM_ESM1b_Complete"
```

---

## 3. Running

### Full pipeline (WT FASTA $\rightarrow$ Mutant synthesis $\rightarrow$ NetSurfP $\rightarrow$ TSV ensemble)

```bash
python netsurfp3_pipeline.py \
    ../../FASTA_files/newnt3/ \
    results/ \
    --mutation-dir ../../mutations/combined/aa \
    --model /path/to/trained_checkpoint.pth \
    --config NetSurfP3/nsp3/experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml \
    --log /path/to/validation.log
# Writes per gene: results/{GENE}/NetSurfP3/{GENE}.netsurfp3.{summary,residues,local}.tsv
```

This flow:
1. Loads ORF nucleotide FASTAs and translates to amino acids
2. Synthesizes mutant AA sequences from mutation CSVs
3. Runs NetSurfP-3.0 on both WT and mutant sequences
4. Generates ensemble TSV with structural change analysis

---

## 4. Output Files

### Output Structure

Output is written per gene to:

```
{output}/
  {GENE}/
    NetSurfP3/
      {GENE}.netsurfp3.summary.tsv   -- per-mutation structural summary
      {GENE}.netsurfp3.residues.tsv  -- per-residue predictions
      {GENE}.netsurfp3.local.tsv     -- local context changes ($\pm$5 residues)
```

### 4.1 Summary TSV (`{GENE}.netsurfp3.summary.tsv`)

Per-mutation summary with structural changes:

| Column | Description                                          | Units |
|--------|------------------------------------------------------|-------|
| `pkey` | `{GENE}-{MUTATION}` primary key                      | string |
| `mutation_pos` | Position of the mutation                             | residue index |
| `wt_aa`, `mut_aa` | WT and mutant amino acids                            | single-letter |
| `delta_rsa` | Change in relative surface accessibility             | 0-1 |
| `delta_disorder_pf`, `delta_disorder_pt` | Change in disorder probability (two scoring methods) | 0-1 |
| `ss3_change` | 1 if 3-state secondary structure changed, else 0     | 0/1 |
| `ss8_change` | 1 if 8-state secondary structure changed, else 0     | 0/1 |
| `burial_classification` | Buried/Intermediate/Exposed change (-2 to +2)        | integer |
| `disorder_classification` | Ordered/Disordered change (-2 to +2)                 | integer |
| `local_structural_impact` | Sum of $\lvert\Delta\rvert$ within ± 5 residues      | various |
| `global_mean_delta_rsa` | Average RSA change across all residues               | 0-1 |
| `global_ss_changes` | Count of secondary structure changes                 | count |
| `qc_flags` | Quality control flags                                | string |

### 4.2 Per-Residue TSV (`{GENE}.netsurfp3.residues.tsv`)

Raw predictions for all residues:

| Column | Description | Units                     |
|--------|-------------|---------------------------|
| `pos` | Residue position | residue index             |
| `aa` | Amino acid | single-letter             |
| `rsa` | Relative Solvent Accessibility | 0-1 (0=buried, 1=exposed) |
| `asa` | Absolute Solvent Accessibility | $\text{Å}^2$              |
| `ss3` | 3-state secondary structure | H/E/C (helix/strand/coil) |
| `ss8` | 8-state secondary structure | G/H/I/B/E/S/T/C           |
| `disorder_pf`, `disorder_pt` | Disorder probabilities (two scoring methods) | 0-1                       |
| `phi`, `psi` | Backbone dihedral angles | degrees                   |
| `sequence_type` | wt or mut | categorical               |

### 4.3 Local Changes TSV (`{GENE}.netsurfp3.local.tsv`)

Changes in the $\pm$5 residue window around each mutation (14 columns):

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | Mutation identifier | string |
| `relative_pos` | Position relative to mutation (-5 to +5) | residue index |
| `absolute_pos` | Absolute residue position in protein | residue index |
| `delta_rsa` | Change in RSA at this position | 0-1 |
| `delta_disorder_pf`, `delta_disorder_pt` | Change in disorder at this position | 0-1 |
| `wt_ss3`, `mut_ss3` | WT and MUT 3-state secondary structure | H/E/C |
| `wt_ss8`, `mut_ss8` | WT and MUT 8-state secondary structure | G/H/I/B/E/S/T/C |
| `delta_phi`, `delta_psi` | Change in backbone dihedral angles | degrees |
| `is_mutation_site` | Whether this is the mutation position | 0/1 |

---

## 5. Structural Feature Interpretation

### 5.1 Surface Accessibility (RSA)

| RSA Range | Classification | Biological Relevance |
|-----------|---------------|---------------------|
| 0 - 0.25 | **Buried** | Core residues, hydrophobic interactions |
| 0.25 - 0.50 | **Intermediate** | Partially exposed, potential conformational flexibility |
| 0.50 - 1.0 | **Exposed** | Surface residues, interaction interfaces, epitopes |

**Delta RSA Interpretation:**
- **$\Delta > +0.2$**: Buried $\rightarrow$ Exposed transition (potential new interaction site, epitope exposure)
- **$\Delta < -0.2$**: Exposed $\rightarrow$ Buried transition (loss of interaction site, epitope masking)
- **$|\Delta| > 0.1$**: Significant accessibility change

### 5.2 Secondary Structure (SS3)

| Code | Structure | Description |
|------|-----------|-------------|
| **H** | Helix | $\alpha$-helix (regular H-bonding) |
| **E** | Strand | $\beta$-strand (extended conformation) |
| **C** | Coil | Loop/turn (irregular structure) |

**Important Transitions:**
- **H $\rightarrow$ C or E $\rightarrow$ C**: Loss of regular structure (potential destabilization)
- **C $\rightarrow$ H or C $\rightarrow$ E**: Gain of regular structure (potential stabilization)
- **H $\leftrightarrow$ E**: Major structural rearrangement

### 5.3 Disorder

| Disorder Score | Classification | Biological Context |
|---------------|---------------|-------------------|
| < 0.3 | **Ordered** | Structured, stable fold |
| 0.3 - 0.7 | **Intermediate** | Context-dependent folding |
| > 0.7 | **Disordered** | Intrinsically disordered region (IDR) |

**Clinical Relevance:**
- Disorder $\rightarrow$ Order: May stabilize protein, affect binding
- Order $\rightarrow$ Disorder: May destabilize, create flexible linkers

---

## 6. Change Classifications

### 6.1 Burial Classification

```
                    WT
           Buried  Intermediate  Exposed
Mutant
Buried       0        -1           -2
Intermediate +1        0           -1
Exposed      +2       +1            0
```

- **+2 (Buried $\rightarrow$ Exposed)**: Potentially exposes epitopes, creates new interfaces
- **-2 (Exposed $\rightarrow$ Buried)**: May mask epitopes, disrupt existing interactions
- **0 (No change)**: Minimal impact on surface properties

### 6.2 Disorder Classification

```
               WT
        Ordered  Intermediate  Disordered
Mutant
Ordered    0        -1            -2
Intermediate +1      0            -1
Disordered   +2     +1             0
```

- **+2 (Ordered $\rightarrow$ Disordered)**: May increase flexibility, affect binding dynamics
- **-2 (Disordered $\rightarrow$ Ordered)**: May rigidify structure, alter conformational ensemble

---

## 7. Command-Line Options

### Positional Arguments
- `input`: Input WT FASTA file or directory of FASTA files (nucleotide or amino acid)
- `output`: Output base directory

### Required
- `--mutation-dir DIR`: Mutation CSV directory
- `--model PATH`: Trained nsp3 model checkpoint from `nsp3 train` (see section 2)
- `--config PATH`: nsp3 experiment config YAML matching the model architecture (see section 2.1)

### Processing Options
- `--log FILE`: Validation log to skip failed mutations
- `--batch-size N`: Number of sequences to process per NSP3 batch (default: 100)
- `--max-seq-length N`: Maximum sequence length before chunking (default: 1500)
- `--input-type {nt, aa}`: Whether input FASTAs contain nucleotide or amino acid sequences (default: `nt`)
- `--verbose`: Enable verbose output

---

## 8. Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `NetSurfP not found` | Confirm `nsp3/` is cloned correctly and PyTorch/ESM weights are available |
| `PyTorch/ESM not found` | Install PyTorch and provide valid `--model` and `--config` paths |
| Execution timeout | Reduce `--max-seq-length` or `--batch-size` (NetSurfP can be slow for long sequences) |
| `No mapping file found` | Verify `--mutation-dir` contains `{GENE}*.csv` files |
| `Out of memory` | Reduce `--batch-size` (default: 100 sequences per batch) |

---

## 9. Use Cases

### 9.1 Immunogenicity Assessment

Mutations can alter:
1. **Surface accessibility** $\rightarrow$ Epitope exposure/masking
2. **Secondary structure** $\rightarrow$ Conformational epitopes
3. **Disorder regions** $\rightarrow$ Protease susceptibility

```bash
python netsurfp3_pipeline.py \
    wt_fastas/ \
    results/ \
    --mutation-dir mappings/ \
    --model /path/to/trained_checkpoint.pth \
    --config NetSurfP3/nsp3/experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml
```

### 9.2 Protein Stability Analysis

Identify mutations causing:
1. **Burial changes** $\rightarrow$ Hydrophobic core disruption
2. **Disorder gain** $\rightarrow$ Destabilization
3. **Secondary structure loss** $\rightarrow$ Unfolding

```bash
python netsurfp3_pipeline.py \
    wt_fastas/ \
    results/ \
    --mutation-dir mappings/ \
    --model /path/to/trained_checkpoint.pth \
    --config NetSurfP3/nsp3/experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml
```

### 9.3 Interaction Interface Mapping

Detect changes in:
1. **Local structural context** $\rightarrow$ Altered residue geometry near interface
2. **RSA transitions** $\rightarrow$ Burial changes at interaction surfaces
3. **Surface patches** $\rightarrow$ Modified solvent-accessible regions

---

## 10. Example Workflow

```bash
# 1. Prepare WT ORF FASTAs
ls ../../FASTA_files/newnt3/
# ABCB1_nt.fasta, BRCA1_nt.fasta, TP53_nt.fasta, ...

# 2. Ensure mutation CSVs exist
ls ../../mutations/combined/aa/
# ABCB1_transcript_mapping.csv, BRCA1_transcript_mapping.csv, ...

# 3. Run full pipeline
python netsurfp3_pipeline.py \
    ../../FASTA_files/newnt3/ \
    results/ \
    --mutation-dir ../../mutations/combined/aa \
    --model /path/to/trained_checkpoint.pth \
    --config NetSurfP3/nsp3/experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml
# Writes per gene: results/{GENE}/NetSurfP3/{GENE}.netsurfp3.{summary,residues,local}.tsv

# 4. Analyze outputs (example for BRCA1)
# results/BRCA1/NetSurfP3/BRCA1.netsurfp3.summary.tsv   - Per-mutation structural summary
# results/BRCA1/NetSurfP3/BRCA1.netsurfp3.residues.tsv  - Per-residue predictions
# results/BRCA1/NetSurfP3/BRCA1.netsurfp3.local.tsv     - Local context changes ($\pm$5 residues)

# 5. Filter for significant structural changes
# Look for:
# - delta_rsa > 0.2 (major burial changes)
# - delta_disorder_pf > 0.2 (disorder transitions)
# - ss3_change = 1 (secondary structure changed)
```

---

## 11. Integration with Other Pipelines

NetSurfP-3.0 predictions complement other BioFeatureFactory analyses:

| Pipeline | NetSurfP Synergy |
|----------|-----------------|
| **NetMHC** | RSA > 0.5 + predicted epitope = accessible epitope |
| **NetNglyc** | Surface glycosylation sites (RSA > 0.4) are functionally relevant |
| **NetPhos** | Exposed phosphorylation sites (RSA > 0.3) are kinase-accessible |
| **SpliceAI** | Splice-altering mutations $\rightarrow$ RSA changes $\rightarrow$ functional impact |

---

**Note**: This pipeline is designed to integrate seamlessly with the existing BioFeatureFactory ecosystem and follows the same patterns as NetNglyc, NetPhos, and NetMHC pipelines.

---

## 12. Citations

- NetSurfP-3.0: Hoie, M. H. et al. (2022). "NetSurfP-3.0: accurate and fast prediction of protein structural features by protein language models and deep learning." *Nucleic Acids Research*, 50(W1), W510-W515.
- Klausen, M. S. et al. (2019). "NetSurfP-2.0: Improved prediction of protein structural features by integrated deep learning." *Proteins*, 87(6), 520-527.

---

## License

This pipeline wrapper is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

The NetSurfP-3.0 library (`nsp3/`) is third-party software with its own license.
