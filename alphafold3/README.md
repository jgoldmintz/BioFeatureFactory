# AlphaFold3 RNA-RBP Interaction Pipeline

## Overview

The AlphaFold3 pipeline performs **structure-based comparative analysis** of RNA-binding protein (RBP) interactions across wild-type (WT) and mutant transcripts.
For each mutation, the pipeline queries POSTAR3/ENCODE eCLIP binding data to identify RBPs with overlapping peaks, then runs AF3 to predict RNA-protein complex structures for both alleles.
This enables quantitative assessment of how mutations perturb RBP binding.

AlphaFold3 (Abramson *et al.*, 2024) provides the structural prediction engine; the Python layer structures inputs, extracts binding metrics, and computes Δ-based features suitable for biological interpretation and downstream modeling.

---

## Capabilities

- **Unified WT<->MUT execution** — WT predictions cached and reused across mutations in the same gene.
- **RBP discovery** — Automatic query of POSTAR3/ENCODE eCLIP binding sites within configurable window (±50 bp).
- **Multi-mode execution** — Local GPU, SLURM batch, or GCP cloud submission.
- **Δ-based comparison** — Per-RBP delta metrics quantify mutation-driven perturbation.
- **Event classification** — Gained, lost, strengthened, weakened binding states.
- **Distance-weighted impact modeling** — Effects scaled by proximity to SNV using inverse distance weighting (Shepard, 1968).
- **Validation-aware filtering** — Exclusion of failed variants via `--validation-log`.

---

## Dependencies

- **AlphaFold3** installation (local mode) or batch/cloud infrastructure
- **Python ≥ 3.8**
- **pysam** (for tabix queries)
- **utility.py** helper functions: `get_mutation_data`, `load_validation_failures`, `should_skip_mutation`

---

## Inputs

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--mutations` | TSV file with mutation specifications |
| `--postar3-db` | Path to tabix-indexed POSTAR3 BED file |
| `--uniprot-mapping` | TSV mapping gene names to UniProt IDs |
| `--sequence-fasta` or `--sequence-dir` | Protein sequences (FASTA file or directory) |
| `--output-dir` | Output directory for results |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--execution-mode` | `local` | Execution mode: `local`, `batch`, or `cloud` |
| `--window-size` | `50` | bp radius around mutation for RBP query |
| `--rna-window` | `101` | RNA sequence length for AF3 input (nt) |
| `--validation-log` | — | Path to validation log for filtering failed mutations |
| `--gcp-project` | — | GCP project ID (cloud mode only) |

### Mutation File Format

```
pkey           gene    chrom   pos        ref  alt  transcript_seq
TP53-c.123A>G  TP53    chr17   7676154    A    G    ATGCGA...
```

Required columns: `pkey`, `gene`, `chrom`, `pos`, `ref`, `alt`, `transcript_seq`

---

## Outputs

### 1. Summary Table (`summary.tsv`)

Each row represents a single mutation with aggregated RBP binding analysis.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Unique variant identifier in the format `GENE-mutation`. | — |
| **n_rbps_tested** | Number of RBPs evaluated for this mutation. | count |
| **n_rbps_binding_wt** | RBPs with confident binding in WT structure. | count |
| **n_rbps_binding_mut** | RBPs with confident binding in MUT structure. | count |
| **global_count_gained** | Number of RBPs with newly gained binding ($S_{\text{wt}} < \tau$, $S_{\text{mut}} \geq \tau$). | count |
| **global_count_lost** | Number of RBPs with lost binding ($S_{\text{wt}} \geq \tau$, $S_{\text{mut}} < \tau$). | count |
| **global_count_strengthened** | RBPs with strengthened binding (both bind, MUT has lower PAE). | count |
| **global_count_weakened** | RBPs with weakened binding (both bind, MUT has higher PAE). | count |
| **global_max_abs_delta_pae** | Maximum $\|\Delta_{\text{PAE}}\|$ across all RBPs. | Å |
| **top_event_rbp** | RBP with highest priority score. | string |
| **top_event_class** | Event classification: `gained`, `lost`, `strengthened`, `weakened`, `none`. | enum |
| **top_event_delta_pae** | $\Delta_{\text{PAE}}$ for top-ranked RBP. | Å |
| **qc_flags** | Semicolon-separated quality flags. | — |

---

### 2. Events Table (`events.tsv`)

Each row represents a single RBP binding comparison for one mutation.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | — |
| **rbp_name** | RBP gene symbol (e.g., HNRNPA1, SRSF1). | string |
| **wt_chain_pair_pae_min** | Minimum PAE between RNA-protein chains (WT). | Å |
| **mut_chain_pair_pae_min** | Minimum PAE between RNA-protein chains (MUT). | Å |
| **delta_chain_pair_pae_min** | $\Delta_{PAE} = PAE_{mut} - PAE_{wt}$ | Å |
| **wt_interface_contacts** | Number of cross-chain contacts < 8 Å (WT). | count |
| **mut_interface_contacts** | Number of cross-chain contacts < 8 Å (MUT). | count |
| **delta_interface_contacts** | Change in interface contacts. | count |
| **cls** | Event classification: `gained`, `lost`, `strengthened`, `weakened`, `unchanged`, `no_binding`. | enum |
| **priority** | Ranking score: $\|\Delta_{\text{PAE}}\| \cdot e^{-d/k}$. | numeric |

---

### 3. Sites Table (`sites.tsv`)

Each row represents a single RNA position within the AF3 prediction window.

| Column | Description | Units |
|--------|-------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | — |
| **allele** | Allele label (`WT` or `MUT`). | string |
| **rbp_name** | RBP identifier. | string |
| **rna_pos** | Position within RNA window (1-indexed). | nt |
| **plddt** | Per-residue confidence score. | 0–100 |
| **in_contact** | Binary flag: 1 if residue contacts protein (< 8 Å). | 0/1 |
| **contact_distance** | Distance to nearest protein atom. | Å |

---

## Key Quantitative Features

### Δchain_pair_pae_min (Primary Binding Proxy)

$$\Delta_{PAE} = PAE_{mut} - PAE_{wt}$$

**What is PAE?** Predicted Aligned Error (PAE) is an AlphaFold confidence metric that estimates the positional error (in Angstroms) between two residues after optimal superposition. For inter-chain predictions, PAE quantifies how confidently AF3 predicts the relative positioning of the RNA and protein chains.

**Interpretation:**
- **Low PAE (< 10 Å)**: High confidence that the two chains interact in the predicted orientation
- **High PAE (> 20 Å)**: Low confidence; chains may not interact or their relative position is uncertain

**Delta interpretation:**
- **Positive Δ** -> Weaker binding in mutant (PAE increased, meaning higher uncertainty)
- **Negative Δ** -> Stronger binding in mutant (PAE decreased, meaning lower uncertainty)
- **Δ ≈ 0** -> No significant change in predicted binding confidence

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

**Interpretation:** Counts atom pairs where an RNA atom is within 8 Å of a protein atom. Higher contact counts suggest more extensive binding interfaces. The 8 Å threshold captures van der Waals contacts, hydrogen bonds, and electrostatic interactions.

---

### Confident Binding Classification

A binding event is classified as **confident** when ALL of the following conditions are met:

1. **Sufficient contacts**: $N_{contacts} \geq 3$
   - At least 3 RNA-protein atom pairs within 8 Å

2. **Low inter-chain uncertainty**: $PAE_{min} \leq 10$ Å
   - The minimum PAE between any RNA-protein residue pair is below 10 Å

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
| **unchanged** | Both bind confidently, $\|\Delta_{\text{PAE}}\| \leq 2.0$ |
| **no_binding** | Neither WT nor MUT has confident binding |

---

## Example Commands

### Local Execution (GPU Required)

```bash
python alphafold3_pipeline.py \
    --mutations mutations.tsv \
    --postar3-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --uniprot-mapping human_uniprot_genes.tsv \
    --sequence-fasta rbp_sequences.fasta \
    --output-dir af3_results/ \
    --execution-mode local
```

### Batch Execution (SLURM)

```bash
python alphafold3_pipeline.py \
    --mutations mutations.tsv \
    --postar3-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --uniprot-mapping human_uniprot_genes.tsv \
    --sequence-dir rbp_sequences/ \
    --output-dir af3_results/ \
    --execution-mode batch
```

Generates SLURM scripts in output directory. Submit with `sbatch af3_results/<job_id>/submit.slurm`.

### Cloud Execution (GCP Batch)

```bash
python alphafold3_pipeline.py \
    --mutations mutations.tsv \
    --postar3-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --uniprot-mapping human_uniprot_genes.tsv \
    --sequence-fasta rbp_sequences.fasta \
    --output-dir af3_results/ \
    --execution-mode cloud \
    --gcp-project my-project
```

Generates GCP Batch JSON configs.

---

## Interpretation

- **Localized perturbations** (large $|\Delta_{\text{PAE}}|$ at high priority) -> direct disruption of RBP binding interface.
- **Multiple gained/lost events** -> mutation affects RBP binding landscape broadly.
- **Low global Δ** with unchanged classifications -> negligible regulatory change.
- **High Δ with lost classification** -> potential functional impact on post-transcriptional regulation.

---

## Thresholds and Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_contacts` | 3 | Minimum interface contacts for confident binding |
| `max_pae_binding` | 10.0 | Maximum PAE for confident binding (Å) |
| `min_plddt_interface` | 50.0 | Minimum pLDDT at interface |
| `delta_pae_significant` | 2.0 | Δ threshold for strengthened/weakened classification (Å) |
| `delta_contacts_significant` | 2 | Contact change threshold for classification |
| `contact_threshold` | 8.0 | Distance threshold for interface contacts (Å) |

---

## QC Flags

| Flag | Condition | Suggested Action |
|------|-----------|------------------|
| `no_rbps_in_region` | No RBP peaks within query window | Exclude variant |
| `rbp_sequence_missing` | Could not find protein sequence | Check mapping file |
| `af3_failed_wt` | AF3 prediction failed on WT | Exclude variant |
| `af3_failed_mut` | AF3 prediction failed on MUT | Exclude variant |
| `low_confidence` | pLDDT < 50 for all interface residues | Use cautiously |
| `no_interface_contacts` | No cross-chain contacts detected | May indicate no binding |
| `token_limit_exceeded` | RNA + protein > 5120 tokens | Reduce window size |
| `PASS` | No issues detected | Safe for analysis |

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
├── alphafold3_pipeline.py      # Main entry point
├── README.md
└── bin/
    ├── af3_runner.py           # AF3 execution (local/batch/cloud)
    ├── af3_parser.py           # Parse mmCIF + JSON outputs
    ├── rbp_database.py         # POSTAR3 tabix query interface
    ├── rbp_sequence_mapper.py  # RBP name -> UniProt -> sequence
    └── binding_metrics.py      # Delta computation and classification
```

## Module Reference

| File | Purpose |
|------|---------|
| `alphafold3_pipeline.py` | Main entry point |
| `bin/af3_runner.py` | AF3 execution (local/batch/cloud) |
| `bin/af3_parser.py` | Parse mmCIF + JSON outputs |
| `bin/rbp_database.py` | POSTAR3 tabix query interface |
| `bin/rbp_sequence_mapper.py` | RBP name -> UniProt -> sequence |
| `bin/binding_metrics.py` | Delta computation and classification |

---

## References

- Abramson J. *et al.* (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. **Nature**, 630:493–500.
- Jumper J. *et al.* (2021) Highly accurate protein structure prediction with AlphaFold. **Nature**, 596:583–589.
- Shepard D. (1968) A two-dimensional interpolation function for irregularly-spaced data. **Proceedings of the 23rd ACM National Conference**, pp. 517–524.
- Zhao W. *et al.* (2022) POSTAR3: an updated platform for exploring post-transcriptional regulation coordinated by RNA-binding proteins. **Nucleic Acids Research**, 50:D483–D492.
- Van Nostrand E.L. *et al.* (2020) A large-scale binding and functional map of human RNA-binding proteins. **Nature**, 583:711–719.

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
