# AlphaFold3 RNA-RBP Interaction Pipeline

## Overview

The AlphaFold3 pipeline performs **structure-based comparative analysis** of RNA-binding protein (RBP) interactions across wild-type (WT) and mutant transcripts.
For each mutation, the pipeline queries POSTAR3/ENCODE eCLIP binding data to identify RBPs with overlapping peaks, then runs AF3 to predict RNA-protein complex structures for both alleles.
This enables quantitative assessment of how mutations perturb RBP binding.

AlphaFold3 (Abramson *et al.*, 2024) provides the structural prediction engine; the Python layer structures inputs, extracts binding metrics, and computes $\Delta$-based features suitable for biological interpretation and downstream modeling.

---

## Capabilities

- **Unified WT $\leftrightarrow$ MUT execution** -- WT predictions cached and reused across mutations in the same gene.
- **RBP discovery** -- Automatic query of POSTAR3/ENCODE eCLIP binding sites within configurable window (±50 bp).
- **Multi-mode execution** -- Local GPU, SLURM batch, or GCP cloud submission.
- **$\Delta$-based comparison** -- Per-RBP delta metrics quantify mutation-driven perturbation.
- **Event classification** -- Gained, lost, strengthened, weakened binding states.
- **Distance-weighted impact modeling** -- Effects scaled by proximity to SNV using inverse distance weighting (Shepard, 1968).
- **Ensemble aggregation** -- Parses all AF3 seed/sample outputs (mean ± std across samples) for robust confidence estimates.
- **Interface sites table** -- Per-residue pLDDT and contact data at the RNA-protein interface, with per-residue contact frequency from ensemble.
- **Multi-window mode** -- Optional placement of the mutation at multiple fractional offsets within the RNA window, aggregating across windows to reduce positional bias.
- **Validation-aware filtering** -- Exclusion of failed variants via `--validation-log`.

---

## Dependencies

- **AlphaFold3 Docker image** (local mode) or batch/cloud infrastructure
- **NVIDIA GPU** with CUDA drivers and [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
- **Docker**
- **Python $\geq$ 3.8**
- **pysam** (for tabix queries)
- **biofeaturefactory.utils.utility** helper functions

---

## AlphaFold3 Setup

Follow the official [AlphaFold3 installation guide](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md) through the **"Obtaining Model Parameters"** step. Stop before **"Obtaining Genetic Databases"** -- this pipeline does not use AlphaFold3's genetic database pipeline. Instead, RBP binding site data is provided by the POSTAR3 database (`--postar-db`) and protein sequences/MSAs are supplied directly (`--msa-dir` or `--rbp-sequences`).

After completing the installation guide steps (NVIDIA drivers, Docker, Container Toolkit, model weights), build the Docker image:

```bash
cd alphafold3
docker build -t alphafold3 -f docker/Dockerfile .
```

Pass the model weights directory to this pipeline via `--model-dir`. At runtime, the pipeline mounts it into the container at `/root/models`.

---

## Inputs

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--postar-db` | Path to tabix-indexed POSTAR3 BED file |
| `--rbp-mapping` | Gene-UniProt mapping TSV |
| `--fasta` | Transcript FASTA file or directory of FASTA files |
| `--output` / `-o` | Output directory |

One of:
- `--rbp-sequences` -- Protein sequences FASTA
- `--msa-dir` -- Directory with A3M MSA files (preferred)

One of:
- `--mutations` -- Mutations CSV file or directory
- `--chromosome-mapping` -- Chromosome mapping CSV (provides mutations and chromosomal positions)

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--execution-mode` | `local` | Execution mode: `local`, `batch`, or `cloud` |
| `--af3-binary` | `alphafold3` | Path to AF3 executable |
| `--docker-image` | `alphafold3` | Docker image name for AF3 |
| `--model-dir` | -- | Path to AF3 model weights directory (required for local mode) |
| `--window-size` | `101` | RNA window size around mutation (nt, odd number) |
| `--rbp-window` | `50` | Window for RBP site lookup (±bp) |
| `--validation-log` | -- | Path to validation log for filtering failed mutations |
| `--vcf` | -- | Per-gene VCF file or directory (provides chromosome) |
| `--chromosome-mapping` | -- | Chromosome mapping CSV file or directory |
| `--chrom` | -- | Chromosome (alternative to `--vcf`) |
| `--tx-start` | -- | Transcript start position (alternative to `--chromosome-mapping`) |
| `--strand` | `+` | Strand (`+`/`-`) |
| `--multi-window` | off | Run multiple windows per mutation (multiplies AF3 runs) |
| `--multi-window-offsets` | `0.3,0.5,0.7` | Mutation position as fraction of window |
| `--max-gpus` | auto-detect | Max GPUs for parallel AF3 execution |

---

## Outputs

### Output Structure

Output is written per gene to:

```
{output}/
  {GENE}/
    AlphaFold3/
      {GENE}.tsv          -- per-mutation summary
      {GENE}.events.tsv   -- per-RBP binding comparison
      {GENE}.sites.tsv    -- per-residue interface data
```

### 1. Summary Table (`{GENE}.tsv`)

Each row represents a single mutation with aggregated RBP binding analysis.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Unique variant identifier in the format `GENE-mutation`. | -- |
| **Gene** | Gene symbol. | string |
| **n_rbps_tested** | Number of RBPs evaluated for this mutation. | count |
| **n_rbps_binding_wt** | RBPs with confident binding in WT structure. | count |
| **n_rbps_binding_mut** | RBPs with confident binding in MUT structure. | count |
| **global_count_gained** | Number of RBPs with newly gained binding ($S_{\text{wt}} < \tau$, $S_{\text{mut}} \geq \tau$). | count |
| **global_count_lost** | Number of RBPs with lost binding ($S_{\text{wt}} \geq \tau$, $S_{\text{mut}} < \tau$). | count |
| **global_count_strengthened** | RBPs with strengthened binding (both bind, MUT has lower PAE). | count |
| **global_count_weakened** | RBPs with weakened binding (both bind, MUT has higher PAE). | count |
| **global_max_abs_delta_pae** | Maximum $\lvert\Delta_{\text{PAE}}\rvert$ across all RBPs. | Å |
| **top_event_rbp** | RBP with highest priority score. | string |
| **top_event_class** | Event classification: `gained`, `lost`, `strengthened`, `weakened`, `none`. | enum |
| **top_event_delta_pae** | $\Delta_{\text{PAE}}$ for top-ranked RBP. | Å |
| **qc_flags** | Semicolon-separated quality flags. | -- |

---

### 2. Events Table (`{GENE}.events.tsv`)

Each row represents a single RBP binding comparison for one mutation.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | -- |
| **rbp_name** | RBP gene symbol (e.g., HNRNPA1, SRSF1). | string |
| **wt_chain_pair_pae_min** | Minimum PAE between RNA-protein chains (WT). | Å |
| **mut_chain_pair_pae_min** | Minimum PAE between RNA-protein chains (MUT). | Å |
| **delta_chain_pair_pae_min** | $\Delta_{PAE} = PAE_{mut} - PAE_{wt}$ | Å |
| **wt_interface_contacts** | Number of cross-chain contacts $< 8 Å$ (WT). | count |
| **mut_interface_contacts** | Number of cross-chain contacts $< 8 Å$ (MUT). | count |
| **delta_interface_contacts** | Change in interface contacts. | count |
| **cls** | Event classification: `gained`, `lost`, `strengthened`, `weakened`, `unchanged`, `no_binding`. | enum |
| **priority** | Ranking score: $\lvert\Delta_{\text{PAE}}\rvert \cdot \frac{1}{1+d/k}$ | numeric |
| **n_samples_wt** | Number of AF3 samples parsed for WT (ensemble mode). | count |
| **n_samples_mut** | Number of AF3 samples parsed for MUT (ensemble mode). | count |
| **std_pae_wt** | Standard deviation of chain-pair PAE across WT samples. | Å |
| **std_pae_mut** | Standard deviation of chain-pair PAE across MUT samples. | Å |
| **n_windows** | Number of windows used (multi-window mode only). | count |

---

### 3. Sites Table (`{GENE}.sites.tsv`)

Each row represents a single interface residue (RNA or protein) within the AF3 prediction.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | -- |
| **rbp_name** | RBP identifier. | string |
| **allele** | Allele label (`WT` or `MUT`). | string |
| **chain** | Chain identifier (e.g., `A` for RNA, `B` for protein). | string |
| **res_id** | Residue sequence number. | int |
| **res_name** | Residue name (3-letter code). | string |
| **plddt** | Per-residue pLDDT confidence score. | 0-100 |
| **is_contact** | Binary flag: 1 if residue contacts the other chain ($< 8 Å$). | 0/1 |
| **min_contact_distance** | Minimum distance to nearest atom in the other chain. | Å |
| **contact_frequency** | Fraction of ensemble samples where this residue is a contact (ensemble mode). | 0-1 |

---

## Key Quantitative Features

### $\Delta$ chain_pair_pae_min (Primary Binding Proxy)

$$\Delta_{PAE} = PAE_{mut} - PAE_{wt}$$

**What is PAE?** Predicted Aligned Error (PAE) is an AlphaFold confidence metric that estimates the positional error (in Angstroms) between two residues after optimal superposition. For inter-chain predictions, PAE quantifies how confidently AF3 predicts the relative positioning of the RNA and protein chains.

**Interpretation:**
- **Low PAE (< $10 Å$)**: High confidence that the two chains interact in the predicted orientation
- **High PAE (> $20 Å$)**: Low confidence; chains may not interact or their relative position is uncertain

**Delta interpretation:**
- **Positive $\Delta$** $\rightarrow$ Weaker binding in mutant (PAE increased, meaning higher uncertainty)
- **Negative $\Delta$** $\rightarrow$ Stronger binding in mutant (PAE decreased, meaning lower uncertainty)
- **$\Delta \approx 0$** $\rightarrow$ No significant change in predicted binding confidence

---

### Priority Score (Distance-Weighted Perturbation)

$$P = |\Delta_{PAE}| \times \frac{1}{1 + d/k}$$

**Variables:**
- $|\Delta_{PAE}|$ = absolute magnitude of PAE change (Å)
- $d$ = genomic distance (bp) between the SNV position and the RBP binding site center
- $k$ = decay constant (default: 50 bp)

**Rationale:** Mutations closer to an RBP binding site are more likely to directly affect binding. The hyperbolic decay term $\frac{1}{1 + d/k}$ down-weights distal effects:

| Distance (d) | Decay Factor | Interpretation |
|--------------|--------------|----------------|
| 0 bp | 1.0 | Full weight (mutation at binding site) |
| 50 bp | 0.5 | Half weight |
| 100 bp | 0.33 | One-third weight |
| 200 bp | 0.2 | One-fifth weight |

This prioritizes RBPs whose binding sites overlap or are proximal to the mutation.

---

### Interface Contacts

$$N_{contacts} = \sum_{i \in RNA} \sum_{j \in protein} \mathbf{1}[d_{ij} < 8 \text{ Å}]$$

**Variables:**
- $i$ = index over RNA heavy atoms
- $j$ = index over protein heavy atoms
- $d_{ij}$ = Euclidean distance between atoms $i$ and $j$
- $\mathbf{1}[\cdot]$ = indicator function (returns 1 if condition is true, 0 otherwise)

**Interpretation:** Counts atom pairs where an RNA atom is within $8 Å$ of a protein atom. Higher contact counts suggest more extensive binding interfaces. The $8 Å$ threshold captures van der Waals contacts, hydrogen bonds, and electrostatic interactions.

---

### Confident Binding Classification

A binding event is classified as **confident** when ALL of the following conditions are met:

1. **Sufficient contacts**: $N_{contacts} \geq 3$
   - At least 3 RNA-protein atom pairs within $8 Å$

2. **Low inter-chain uncertainty**: $PAE_{min} \leq 10 Å$
   - The minimum PAE between any RNA-protein residue pair is below $10 Å$

3. **Adequate local structure quality**: $pLDDT_{RNA} \geq 50$ OR $pLDDT_{protein} \geq 50$
   - At least one chain at the interface has confident per-residue structure (pLDDT scale: 0-100, where >70 is high confidence)

**Combined criterion:**

$$\text{confident} = (N_{contacts} \geq 3) \land (PAE_{min} \leq 10) \land (pLDDT_{RNA} \geq 50 \lor pLDDT_{protein} \geq 50)$$

This three-part filter excludes predictions where AF3 has low confidence in either the individual chain structures or their relative arrangement.

---

## Event Classification

| Class | Condition |
|-------|-----------|
| **gained** | WT has no confident binding, MUT has confident binding |
| **lost** | WT has confident binding, MUT has none |
| **strengthened** | Both bind confidently, $\Delta_{\text{PAE}} < -2.0$ |
| **weakened** | Both bind confidently, $\Delta_{\text{PAE}} > +2.0$ |
| **unchanged** | Both bind confidently, $\lvert\Delta_{\text{PAE}}\rvert \leq 2.0$ |
| **no_binding** | Neither WT nor MUT has confident binding |
| **incomplete** | WT or MUT prediction missing -- no comparison possible |

---

## Example Commands

### Local Execution (GPU Required)

```bash
python alphafold3_pipeline.py \
    --fasta transcripts/SMN2.fasta \
    --chromosome-mapping mappings/SMN2_chromosome_mapping.csv \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping human_uniprot_genes.tsv \
    --rbp-sequences rbp_sequences.fasta \
    --model-dir /path/to/af3_weights \
    --output af3_results/ \
    --execution-mode local
```

### With Multi-Window Mode

```bash
python alphafold3_pipeline.py \
    --fasta transcripts/SMN2.fasta \
    --chromosome-mapping mappings/SMN2_chromosome_mapping.csv \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping human_uniprot_genes.tsv \
    --msa-dir msa_files/ \
    --model-dir /path/to/af3_weights \
    --output af3_results/ \
    --multi-window \
    --multi-window-offsets 0.3,0.5,0.7
```

### Directory Mode (Multiple Genes)

```bash
python alphafold3_pipeline.py \
    --fasta transcripts/ \
    --chromosome-mapping mappings/ \
    --vcf vcf_files/ \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping human_uniprot_genes.tsv \
    --msa-dir msa_files/ \
    --model-dir /path/to/af3_weights \
    --output af3_results/
```

### Batch Execution (SLURM)

```bash
python alphafold3_pipeline.py \
    --fasta transcripts/SMN2.fasta \
    --mutations mutations/SMN2.csv \
    --chrom chr5 --tx-start 70220768 \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping human_uniprot_genes.tsv \
    --rbp-sequences rbp_sequences.fasta \
    --output af3_results/ \
    --execution-mode batch
```

Generates SLURM scripts in output directory. Submit with `sbatch af3_results/<job_id>/submit.slurm`.

---

## Ensemble Aggregation

AF3 produces multiple samples per seed (typically 4 samples per seed). The pipeline parses all `seed-N_sample-N/` subdirectories from each AF3 run and aggregates metrics across samples.

**Aggregated fields:**
- **Mean** of chain-pair PAE, interface contacts, pLDDT (RNA and protein)
- **Standard deviation** (population stdev) of each metric across samples
- **Contact frequency** -- per-residue fraction of samples where that residue appears as an interface contact (reported in the sites table)

When only a single sample is available (e.g., top-ranked model only), the pipeline falls back to single-model behavior with no standard deviation fields.

The ensemble columns in `events.tsv` (`n_samples_wt`, `n_samples_mut`, `std_pae_wt`, `std_pae_mut`) are empty when running in single-model mode.

---

## Multi-Window Mode

By default, the mutation is centered in the RNA window. The `--multi-window` flag generates multiple windows where the mutation is placed at different fractional offsets within the window (default: 0.3, 0.5, 0.7). This reduces positional bias -- AF3 predictions can be sensitive to where the mutation falls relative to window boundaries.

**Behavior:**
- Each offset produces a distinct WT/MUT window pair
- Duplicate windows (possible near transcript ends) are deduplicated
- AF3 is run separately for each window (jobs suffixed `_w0`, `_w1`, etc.)
- Metrics are aggregated across windows (mean ± std), reported identically to ensemble aggregation
- The `n_windows` column in `events.tsv` records how many windows were used

**Compute cost:** With 3 offsets and N RBPs, multi-window mode runs 6xN AF3 jobs per mutation (3 windows x 2 alleles x N RBPs). The flag is off by default.

---

## Interpretation

- **Localized perturbations** (large $|\Delta_{\text{PAE}}|$ at high priority) $\rightarrow$ direct disruption of RBP binding interface.
- **Multiple gained/lost events** $\rightarrow$ mutation affects RBP binding landscape broadly.
- **Low global $\Delta$** with unchanged classifications $\rightarrow$ negligible regulatory change.
- **High $\Delta$ with lost classification** $\rightarrow$ potential functional impact on post-transcriptional regulation.

---

## Thresholds and Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_contacts` | 3 | Minimum interface contacts for confident binding |
| `max_pae_binding` | 10.0 | Maximum PAE for confident binding (Å) |
| `min_plddt_interface` | 50.0 | Minimum pLDDT at interface |
| `delta_pae_significant` | 2.0 | $\Delta$ threshold for strengthened/weakened classification (Å) |
| `delta_contacts_significant` | 2 | Contact change threshold for classification |
| `contact_threshold` | 8.0 | Distance threshold for interface contacts (Å) |

---

## QC Flags

| Flag | Condition | Suggested Action |
|------|-----------|------------------|
| `PASS` | All RBPs produced complete results | Safe for analysis |
| `PARTIAL` | Some RBPs succeeded, some failed | Interpret with caution |
| `ALL_FAILED` | All RBP predictions failed | Exclude variant or debug AF3 |
| `no_rbps_tested` | No RBPs were evaluated (empty result) | Check upstream RBP lookup |
| `no_rbps_in_region` | No RBP peaks found within query window | Exclude variant |
| `FAILED:{error}` | Individual RBP prediction raised an exception | See stderr for details |

---

## Sign Conventions

| Metric | Sign | Biological Interpretation |
|--------|------|---------------------------|
| `delta_chain_pair_pae_min` | **Positive** | Weaker binding (higher uncertainty) |
|  | **Negative** | Stronger binding (lower uncertainty) |
| `delta_interface_contacts` | **Positive** | More contacts in MUT |
|  | **Negative** | Fewer contacts in MUT |

---

## Caching

- **WT predictions**: Cached by `{gene}-{rbp}` hash, reused across mutations.
- **RBP sequences**: Loaded once from local files at startup.
- **Cache location**: `{output_dir}/cache/`

---

## File Structure

```
alphafold3/
|-- alphafold3_pipeline.py      # Main entry point
|-- README.md
`-- bin/
    |-- af3_runner.py           # AF3 execution (local/batch/cloud)
    |-- af3_parser.py           # Parse mmCIF + JSON outputs
    |-- rbp_database.py         # POSTAR3 tabix query interface
    |-- rbp_sequence_mapper.py  # RBP name -> UniProt -> sequence
    `-- binding_metrics.py      # Delta computation and classification
```

## Module Reference

| File | Purpose |
|------|---------|
| `alphafold3_pipeline.py` | Main entry point |
| `bin/af3_runner.py` | AF3 execution (local/batch/cloud) |
| `bin/af3_parser.py` | Parse mmCIF + JSON outputs |
| `bin/rbp_database.py` | POSTAR3 tabix query interface |
| `bin/rbp_sequence_mapper.py` | RBP name $\rightarrow$ UniProt $\rightarrow$ sequence |
| `bin/binding_metrics.py` | Delta computation and classification |

---

## References

- Abramson J. *et al.* (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. **Nature**, 630:493-500.
- Jumper J. *et al.* (2021) Highly accurate protein structure prediction with AlphaFold. **Nature**, 596:583-589.
- Shepard D. (1968) A two-dimensional interpolation function for irregularly-spaced data. **Proceedings of the 23rd ACM National Conference**, pp. 517-524.
- Zhao W. *et al.* (2022) POSTAR3: an updated platform for exploring post-transcriptional regulation coordinated by RNA-binding proteins. **Nucleic Acids Research**, 50:D483-D492.
- Van Nostrand E.L. *et al.* (2020) A large-scale binding and functional map of human RNA-binding proteins. **Nature**, 583:711-719.

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
