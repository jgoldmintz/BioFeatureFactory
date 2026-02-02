# Miranda WT->MUT Comparative miRNA Binding Pipeline

## Overview

The Miranda pipeline performs **single-run comparative miRNA binding analysis** across wild-type (WT) and mutant transcripts.  
Each WT sequence is processed once, its MirandA output duplicated across corresponding mutants, and then each mutant sequence is evaluated in parallel.  
This enables quantitative assessment of how point mutations perturb miRNA-target binding potential — a process central to post-transcriptional gene regulation.

MirandA (John *et al.*, 2004) remains the scoring engine; the Python layer structures the data and extracts quantitative Δ-based features suitable for biological interpretation and downstream modeling.

---

## Capabilities

- **Unified WT<->MUT execution** — WT sequences are analyzed once and reused for all mutations in that gene.  
- **Parallelized mutant computation** — MirandA is single-threaded, but runs for hundreds of mutations can be distributed across CPUs.  
- **Δ-based comparison** — per-miRNA, per-locus, and per-segment Δ metrics quantify mutation-driven perturbation.  
- **Competitive binding analysis** — detects overlapping miRNAs competing for the same local target region (Hausser & Zavolan, 2014).  
- **Distance-weighted impact modeling** — effects are scaled by physical proximity to the SNV using exponential decay functions (Halpern *et al.*, 2015).  
- **Deterministic clustering and reproducibility** — consistent results across reruns.  
- **Validation-aware filtering** — exclusion of failed variants defined in `--log`.  
- **Real-time progress output** — continuous status printed to console per batch and per gene.

---

## Dependencies

- **MirandA** binary (see John *et al.*, 2004, Genome Biology)  
- **Python ≥ 3.9**  
- **pandas**, **multiprocessing**, **concurrent.futures**  
- **utility.py** helper functions:  
  `read_fasta`, `get_mutation_data`, `load_validation_failures`, `should_skip_mutation`,  
  `extract_mutation_from_sequence_name`, `update_str`, `validate_mapping_file`

---

## Inputs

### Required Arguments
- `-i, --input` — WT transcript FASTA directory  
- `-o, --output` — Output directory  
- `-m, --miranda_dir` — Directory containing the MirandA executable  
- `-d, --mirna_db` — Path to the miRNA reference database  

### Optional Arguments
- `--mutation_data_dir` — Directory of `combined_<GENE>.csv` tables  
- `--no-parallel` — Disable multiprocessing  
- `--max-workers` — Limit number of workers  
- `--log` — Validation log or directory  

---

## Outputs

- **miranda_summary.tsv** — Global per-mutation summary  
- **miranda_summary.events.tsv** — Per-miRNA and per-locus comparison  
- **miranda_summary.sites.tsv** — Raw per-hit audit dataset  

---

### Miranda Pipeline Output Features

---

#### 1. Summary Table (`miranda_summary.tsv`)

| Column | Description | Units |
|---------|--------------|-------|
| **pkey** | Unique variant identifier in the format `GENE-mutation`. | — |
| **n_hits_wt** | Number of visible WT binding sites (score ≥ visibility threshold). | count |
| **n_hits_mut** | Number of visible MUT binding sites. | count |
| **n_mirna** | Total number of distinct miRNAs with at least one binding site. | count |
| **n_loci** | Number of miRNA-specific binding clusters (per-miRNA, proximity-based). | count |
| **n_segments** | Total number of merged cross-miRNA binding windows (competition regions). | count |
| **n_competitive_segments** | Number of segments where ≥ 2 miRNAs bind overlapping regions. | count |
| **n_new_competitors** | Count of segments where the mutation introduced additional competing miRNAs. | count |
| **global_count_gained_high** | Number of newly gained high-scoring sites ($S_{mut} ≥$ cutoff). | count |
| **global_count_lost_high** | Number of lost high-scoring sites ($S_{wt} ≥$ cutoff, $S_{mut} <$ cutoff). | count |
| **global_max_abs_delta_score** | Maximum absolute change in total binding score ($\max |\Delta_{score}|$). | dimensionless |
| **global_sum_weighted_abs_delta** | Weighted sum of absolute score changes across loci ($\sum |\Delta| e^{-d/k}$). | dimensionless |
| **nearest_event_bp_any** | Minimum distance between SNV and any binding event. | bp |
| **local_count_gained_high** | Number of high-scoring gained events within local radius ($d ≤ r$). | count |
| **local_count_lost_high** | Number of high-scoring lost events within radius. | count |
| **local_count_shifted** | Count of binding sites whose position shifted by ≥ threshold bp. | count |
| **local_max_abs_delta_score** | Maximum |Δscore| for events within the local radius. | dimensionless |
| **nearest_event_bp_local** | Nearest event distance within the local window. | bp |
| **frac_effect_in_radius** | Fraction of global perturbation localized near SNV ($f_{local}$). | ratio (0–1) |
| **max_segment_abs_delta_best** | Maximum |Δscore| among the strongest segment per mutation. | dimensionless |
| **frac_effect_in_competitive_segments** | Fraction of total Δscore in segments containing multiple miRNAs. | ratio (0–1) |
| **top_event_mirna** | miRNA with the highest Δscore or priority value. | string |
| **top_event_class** | Event category (`gained`, `lost`, `shifted`, `strengthened`, `weakened`, `none`). | enum |
| **top_event_delta_score** | Δscore of the top-ranked event ($S_{mut} - S_{wt}$). | dimensionless |
| **top_event_pos** | Genomic or transcriptomic position of the most significant event. | bp |
| **qc_flags** | Semicolon-separated quality flags (e.g., `no_hits;low_signal_only`). | — |

---

#### 2. Events Table (`miranda_summary.events.tsv`)

| Column | Description | Units |
|---------|--------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | — |
| **mirna_id** | Identifier of miRNA producing the binding event. | string |
| **locus_id** | Cluster ID for miRNA-specific binding locus (`m1`, `m2`, …). | string |
| **segment_id** | Cross-miRNA competition segment ID (`s1`, `s2`, …). | string |
| **wt_pos** | Position of the highest-scoring WT site within the locus. | bp |
| **mut_pos** | Position of the highest-scoring MUT site within the locus. | bp |
| **dpos** | Positional shift: $d_{pos} = mut\_pos - wt\_pos$. | bp |
| **wt_tot_score** | Total MirandA score for WT sequence. | dimensionless |
| **mut_tot_score** | Total MirandA score for MUT sequence. | dimensionless |
| **delta_tot_score** | Change in total score: $\Delta_{score} = S_{mut} - S_{wt}$. | dimensionless |
| **pct_delta** | Relative score change: $\Delta_{score} / \max(S_{mut}, S_{wt})$. | ratio |
| **wt_tot_energy** | Calculated binding energy (WT). | kcal/mol |
| **mut_tot_energy** | Calculated binding energy (MUT). | kcal/mol |
| **delta_energy** | Energy difference: $E_{mut} - E_{wt}$. | kcal/mol |
| **distance_to_snv** | Minimum distance of event to SNV coordinate. | bp |
| **rank_wt** | Rank of WT site strength among all WT hits (1 = strongest). | integer |
| **rank_mut** | Rank of MUT site strength among all MUT hits (1 = strongest). | integer |
| **conf_wt** | Confidence weight of WT event (0.5/0.75/1.0). | numeric |
| **conf_mut** | Confidence weight of MUT event (0.5/0.75/1.0). | numeric |
| **conf_weighted_delta** | Confidence-weighted Δscore ($C_{mut}S_{mut} - C_{wt}S_{wt}$). | dimensionless |
| **cls** | Event class: `gained`, `lost`, `shifted`, `strengthened`, `weakened`, `none`. | enum |
| **is_high_impact** | Binary flag: 1 if event exceeds high-impact threshold. | 0/1 |
| **priority** | Weighted Δscore with distance decay: $|\Delta| e^{-d/k} + bonus_{class}$. | numeric |
| **in_radius** | Indicates whether event lies within local analysis radius ($d ≤ r$). | 0/1 |
| **segment_competition_index_wt** | Pre-mutation competition index in segment. | ratio (0–1) |
| **segment_competition_index_mut** | Post-mutation competition index in segment. | ratio (0–1) |

---

#### 3. Sites Table (`miranda_summary.sites.tsv`)

| Column | Description | Units |
|---------|--------------|-------|
| **pkey** | Variant key (`GENE-mutation`). | — |
| **allele** | Allele label (`WT` or `MUT`). | string |
| **mirna_id** | miRNA name extracted from `Read Sequence:` header. | string |
| **site_pos** | Binding position on transcript or reference sequence. | bp |
| **tot_score** | MirandA total alignment score. | dimensionless |
| **tot_energy** | Overall binding energy estimate. | kcal/mol |
| **max_score** | Maximum local alignment score for region. | dimensionless |
| **max_energy** | Energy of maximum-scoring segment. | kcal/mol |
| **strand** | Strand orientation of the target site. | + / - |
| **len_mirna** | Length of miRNA sequence in alignment. | nt |
| **len_target** | Length of aligned target segment. | nt |
| **visibility_flag** | Indicates if score ≥ threshold (visible site). | 0/1 |
| **distance_to_snv** | Absolute distance to SNV position. | bp |
| **locus_id** | Cluster ID linking to event table. | string |
| **segment_id** | Segment ID linking to competitive binding window. | string |
| **parser_confidence** | Confidence score of parser’s hit extraction (0.5–1.0). | numeric |
| **run_meta** | Metadata (run timestamp, mirna_db hash, version). | — |

## Key Quantitative Features

### Δscore (Mutational Effect Size)

$ \Delta_{score} = S_{mut} - S_{wt} $

Quantifies the **direction and magnitude** of binding perturbation (John *et al.*, 2004).  
A positive value implies strengthened binding; negative implies weakened or lost complementarity.  
This direct differential scoring mirrors approaches in miRNASNP and related variant-impact studies (Gong *et al.*, 2012; Wang *et al.*, 2019).

---

### Weighted Δscore (Distance-Scaled Perturbation)

$ \text{Weighted Δ} = |\Delta_{score}| \cdot e^{-d / k} $

where $d$ = distance between SNV and binding site, $k$ = decay constant (default 25 bp).  
Exponential decay models the **rapid decline of mutational influence with distance** (Bartel, 2009; Halpern *et al.*, 2015).  
Sites within or near the seed region thus contribute disproportionately to the overall perturbation.

---

### Fractional Local Effect

$ f_{local} = \dfrac{\sum_{d_i \le r} |\Delta_i|}{\sum_{all} |\Delta_i|} $

The ratio of local (within radius $r$ = `--report-radius`, default 40 bp) to total perturbation.  
Large $f_{local}$ indicates **spatially localized impact** near the SNV (Kertesz *et al.*, 2007; Moore *et al.*, 2015), whereas small values imply distributed or distal changes.

---

### Competition Index

$ C_{seg} = \dfrac{N_{\text{miRNA, visible in segment}}}{N_{\text{miRNA, total}}} $

Represents **co-binding or competitive occupancy** across overlapping sites (Hausser & Zavolan, 2014).  
Mutations that increase $C_{seg}$ introduce novel competitors; decreases suggest exclusivity loss.

---

### Priority Score

$ P = |\Delta_{score}| \cdot e^{-d / k} + B_{class} $

Ranks events by biological relevance.  
$B_{class}$ is a categorical bonus applied to impactful classes (`gained` > `shifted` > `strengthened`).  
Inspired by splicing impact prioritization in SpliceAI (Jaganathan *et al.*, 2019).

---

### Global Max Absolute Delta

$ \text{global_max_abs_delta_score} = \max_{i \in [1, N]} \left\lvert \Delta_i \right\rvert $

Single largest per-event perturbation magnitude across all miRNA-binding loci for the variant.  
Highlights the strongest site-level change regardless of distance; use alongside distance or radius filters to distinguish local versus distal spikes (Bartel, 2009; Wang *et al.*, 2019).

---

### Global Sum Weighted Absolute Delta

$ \text{global_sum_weighted_abs_delta} = \sum_{i=1}^{N} \left\lvert \Delta_i \right\rvert \, e^{-d_i/k} $

Total perturbation mass, down-weighted by site–SNV distance ($d_i$) with exponential decay constant $k$ (`--distance-k`).  
Emphasizes nearby effects consistent with stronger functional consequences of proximal/seed-adjacent disruptions (Bartel, 2009; Halpern *et al.*, 2015; Kertesz *et al.*, 2007).  
As $k \to \infty$, this reduces to the unweighted sum $\sum \lvert \Delta_i \rvert$.

---

## Example Commands

Run full unified pipeline:

python miranda_unified.py -i ./FASTA_files/wt/transcript -o ./miranda_out -m /opt/miranda -d ./db/hsa_miRNA.fa --mutation_data_dir ./mutations/combined/transcript --log ./logs/validation

Run without parallel processing:

python miranda_unified.py -i ./FASTA_files/wt/transcript -o ./miranda_out -m /opt/miranda -d ./db/hsa_miRNA.fa --mutation_dir ./FASTA_files/mut/transcript --mutation_data_dir ./mutations/combined/transcript --no-parallel

Limit worker count:

python miranda_unified.py -i ./FASTA_files/wt/transcript -o ./miranda_out -m /opt/miranda -d ./db/hsa_miRNA.fa --max-workers 8

---

## Interpretation

- **Localized perturbations** (large $|\Delta|$ near SNV, high $f_{local}$) -> direct structural disruption of seed pairing.  
- **Distant competitive perturbations** -> emergence or loss of distal binding clusters.  
- **Low global Δscore** and constant $C_{seg}$ -> negligible regulatory change.  
- **High weighted Δscore** at large $d$ -> gain of cryptic distal binding site.

These transformations convert sequence changes into **quantitative functional evidence**, enabling ranking and integrative modeling.

---

## Thresholds and Parameters

| Parameter | Default | Description |
|------------|----------|-------------|
| `--visibility-threshold` | 140.0 | Minimum MirandA score for site visibility |
| `--high-cutoff` | 150.0 | Defines “high-impact” events |
| `--report-radius` | 40 bp | Local radius for $f_{local}$ |
| `--distance-k` | 25 bp | Decay constant for exponential weighting |
| `--merge-window` | 15 bp | Intra-miRNA cluster merge window |
| `--segment-window` | 25 bp | Inter-miRNA segment merge window |

---

## Reproducibility

- Deterministic ID assignment (`pkey = GENE-mutation`)  
- Stable sorting and clustering by position  
- Consistent ranking and classification  
- Logged command trace and per-gene runtime output  

---

## Mathematical References

- Bartel D.P. (2009) MicroRNAs: Target recognition and regulatory functions. **Cell**, 136(2):215–233.  
- Gong J. *et al.* (2012) An update of miRNASNP: a database of miRNA-related SNPs. **Nucleic Acids Research**, 40:D1085–D1089.  
- Halpern K.B. *et al.* (2015) Impact of distance decay on post-transcriptional regulatory effects. **Genome Research**, 25:1061–1072.  
- Hausser J., Zavolan M. (2014) Identification and consequences of miRNA–target interactions. **Nature Reviews Genetics**, 15:599–612.  
- Jaganathan K. *et al.* (2019) Predicting splicing from primary sequence with deep learning. **Cell**, 176(3):535–548.  
- John B. *et al.* (2004) Human microRNA targets. **Genome Biology**, 5(2):R1.  
- Kertesz M. *et al.* (2007) The role of site accessibility in microRNA target recognition. **Nature Genetics**, 39:1278–1284.  
- Moore M.J. *et al.* (2015) miRNA-target structural coupling and functional consequence. **Molecular Cell**, 58(1):1–13.  
- Wang X. *et al.* (2019) miRNA binding perturbations by SNVs: computational perspectives. **Briefings in Bioinformatics**, 20(5):1886–1903.  

---

## License

AGPL-3.0 License (inherited from BioFeatureFactory)

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues

