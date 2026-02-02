# GeneSplicer Pipeline: WT<->ALT Ensemble Delta Caller

## Overview

---

Run GeneSplicer on full genomic context (or configurable context), compare WT vs ALT per SNV



## Capabilities
- Full-context scan (`--pipeline full`) with absolute genomic positions
- Optional local marking with `--report-radius`
- Proximity clustering, WT<->ALT pairing, event calling, priority scoring
- Deterministic, append-safe summary; auditable detail tables

## Dependencies
- GeneSplicer binary at `--genesplicer-dir/genesplicer`; model directory at `../human` relative to the binary
- Python 3.9+, `pandas`, `numpy`

## Inputs
- `--input`: directory of per-gene genomic FASTA files  
  Preference: sequence named `genomic`; else first record in file.
- `--mapping-dir`: directory of per-gene mutation CSVs  
  Mutation column auto-detected: one of `genomic | mutant | mutation | nt_mutation | ntmutant`.  
  Tokens like `G1211T` parsed via `get_mutation_data()`.
- `--log` (optional): validation logs used by `should_skip_mutation()`.

## Outputs

---

### File Set
- `--splice-summary` -> `<basename>.tsv`
- `--splice-summary.events.tsv` -> `<basename>.events.tsv`
- `--splice-summary.sites.tsv` -> `<basename>.sites.tsv`

**Conventions**
- Coordinates are absolute in `genome`/`extended` scan modes; window‐relative only in `window` mode.
- Distances are base pairs (bp) from the SNV.
- Scores are raw GeneSplicer scores (dimensionless).
- Confidence mapping: `low=0.5`, `med/medium=0.75`, `high=1.0`.

---

### 1. Summary Table (`--splice-summary`)
Each row represents a single SNV comparison (`pkey = GENE-mutation`). Global view + optional local triage slice.

| Column                      | Description                                                                                | Units                           |
|-----------------------------|--------------------------------------------------------------------------------------------|---------------------------------|
| `pkey`                      | Unique variant key `GENE-mutation`                                                         | —                               |
| `n_sites_wt`                | Count of **visible** WT sites (score ≥ `--visibility-threshold`)                           | count                           |
| `n_sites_mut`               | Count of **visible** MUT sites                                                             | count                           |
| `n_clusters`                | Number of donor/acceptor clusters evaluated                                                | count                           |
| `global_count_gained_high`  | Gained clusters where max(wt,mut) ≥ `--high-cutoff`                                        | count                           |
| `global_count_lost_high`    | Lost clusters where max(wt,mut) ≥ `--high-cutoff`                                          | count                           |
| `global_count_shifted`      | Clusters with positional shift (                                                           | `dpos`                          | ≥ `--shift-bp`) | count |
| `global_max_abs_Δscore`     | Max absolute score change across all clusters                                              | dimensionless                   |
| `global_sum_weighted_abs_Δ` | Σ\|Δscore\| · exp(−distance/`--distance-k`) across clusters                                | dimensionless                   |
| `nearest_event_bp_any`      | Minimum distance to SNV among all events                                                   | bp                              |
| `local_count_gained_high`   | As above, restricted to `distance ≤ --report-radius`                                       | count                           |
| `local_count_lost_high`     | As above, restricted local                                                                 | count                           |
| `local_count_shifted`       | As above, restricted local                                                                 | count                           |
| `local_max_abs_Δscore`      | Max abs Δscore within local radius                                                         | dimensionless                   |
| `nearest_event_bp_local`    | Nearest event within local radius                                                          | bp                              |
| `frac_effect_in_radius`     | (Σ\|Δ\| in radius) / (Σ\|Δ\| global)                                                       | 0–1                             |
| `top_event_type`            | Event class of highest‐priority cluster (`gained/lost/shifted/strengthened/weakened/none`) | enum                            |
| `top_event_Δscore`          | Δscore of the highest‐priority event                                                       | dimensionless                   |
| `top_event_pos`             | Representative position of top event (MUT if present else WT)                              | bp (absolute) or index (window) |
| `dominant_boundary`         | Boundary with largest Σ\|Δ\| (`donor` or `acceptor`)                                       | enum                            |
| `qc_flags`                  | Semicolon‐separated flags: `no_sites;far_event>2kb;low_signal_only`                        | —                               |

**Interpretation**
- Local impact: `local_max_abs_Δscore` and local gained/lost counts.
- Distant cryptic promotion: `global_count_gained_high > 0` with large `nearest_event_bp_any`.
- Clean negatives: `n_clusters=0` or `global_max_abs_Δscore≈0` and `frac_effect_in_radius≈0`.

---

### 2. Events Table (`--splice-summary.events.tsv`)
Each row represents a donor/acceptor **cluster** (merged nearby sites) and compares top WT vs top MUT within that cluster.

| Column                | Description                                                                               | Units/Values                                     |
|-----------------------|-------------------------------------------------------------------------------------------|--------------------------------------------------|
| `pkey`                | Variant key `GENE-mutation`                                                               | —                                                |
| `type`                | Boundary type                                                                             | `donor` / `acceptor`                             |
| `cluster_id`          | Stable ID within `(pkey, type)` (e.g., `d1`, `a2`)                                        | string                                           |
| `wt_pos`              | Top WT site position in cluster                                                           | bp (absolute) or index (window)                  |
| `mut_pos`             | Top MUT site position in cluster                                                          | bp (absolute) or index (window)                  |
| `dpos`                | Positional shift: `mut_pos − wt_pos`                                                      | bp (or indices)                                  |
| `wt_score`            | Top WT score                                                                              | dimensionless                                    |
| `mut_score`           | Top MUT score                                                                             | dimensionless                                    |
| `dscore`              | Score delta: `mut_score − wt_score`                                                       | dimensionless                                    |
| `pct_delta`           | Relative delta: $r_i = \dfrac{\text{dscore}}{\max( \\|\text{WTscore}\\|, \varepsilon)}$   | dimensionless                                    |
| `distance_to_snv`     | Min distance of (WT or MUT) site in the cluster to the SNV                                | bp                                               |
| `rank_wt`             | Rank of WT site among all WT sites of this `type` (1=strongest)                           | integer                                          |
| `rank_mut`            | Rank of MUT site among all MUT sites of this `type`                                       | integer                                          |
| `conf_wt`             | Confidence weight of WT (`0.5/0.75/1.0`)                                                  | numeric                                          |
| `conf_mut`            | Confidence weight of MUT (`0.5/0.75/1.0`)                                                 | numeric                                          |
| `conf_weighted_delta` | `conf_mut·mut_score − conf_wt·wt_score`                                                   | dimensionless                                    |
| `cls`                 | Event class                                                                               | `gained/lost/shifted/strengthened/weakened/none` |
| `is_high_impact`      | High Δ or high gained/lost (policy thresholds)                                            | 0/1                                              |
| `priority`            | Sorting key: $\\|\text{dscore}\\|*\text{e}^{(−distance/\text{"--distance-k"})}$ + bonuses | numeric                                          |
| `in_radius`           | Inside local triage radius (`distance ≤ --report-radius`)                                 | 0/1                                              |

**Event taxonomy**
- `gained`: WT absent, MUT visible (≥ visibility threshold).
- `lost`: MUT absent, WT visible.
- `shifted`: both visible and `|dpos| ≥ --shift-bp`.
- `strengthened / weakened`: same/near position with `|dscore| ≥ 1.0`.
- `none`: otherwise.

---

### 3. Sites Table (`--splice-summary.sites.tsv`)
Per‐allele per‐site audit rows used to build clusters and events.

| Column            | Description                                         | Units/Values                    |
|-------------------|-----------------------------------------------------|---------------------------------|
| `pkey`            | Variant key `GENE-mutation`                         | —                               |
| `allele`          | Allele label                                        | `WT` / `MUT`                    |
| `type`            | Boundary type                                       | `donor` / `acceptor`            |
| `site_pos`        | Site position (End5 for donors, End3 for acceptors) | bp (absolute) or index (window) |
| `score`           | GeneSplicer score                                   | dimensionless                   |
| `confidence`      | Confidence weight (`low/med/high` -> `0.5/0.75/1.0`) | numeric                         |
| `rank`            | Rank within allele×type (1=strongest)               | integer                         |
| `distance_to_snv` | `\|site_pos − snv_pos\|`                            | bp                              |
| `visible_flag`    | 1 if `score ≥ --visibility-threshold` else 0        | 0/1                             |
| `cluster_id`      | Cluster tag that links to `events`                  | string                          |
| `in_radius`       | Inside local triage radius                          | 0/1                             |

**Consistency checks**
- `n_sites_wt` = count of `visible_flag=1` where `allele=WT` for that `pkey`.
- `n_sites_mut` = count of `visible_flag=1` where `allele=MUT` for that `pkey`.
- `n_clusters` ≥ number of unique `(type, cluster_id)` rows in `events` for that `pkey`.


## Quick heuristics
- Prioritize SNVs with `global_count_gained_high + global_count_lost_high > 0` or large `global_max_abs_Δscore`.
- Far‐field cryptic activation: high‐impact `gained` with large `distance_to_snv`.
- Local disruption: high local counts or `local_max_abs_Δscore` near the SNV.

Priority = `|dscore| * exp(−distance/k)` with bonuses for high-score gained/lost and large shifts.

## Clustering and Pairing
- Single-linkage clustering by position within `(pkey, type)` using `--cluster-radius` (default 3 bp)
- In each cluster, compare top WT vs top MUT to compute deltas and class

---

## Modes and Pipeline
- `--pipeline full`  
  Scan entire genomic sequence for WT and ALT; absolute coordinates; global reporting; local triage flag via `--report-radius`.
- `--pipeline custom` with `--scan-mode`:
  - `genome` — same as full
  - `extended` — slice `[pos − L, pos + L]` with `--context-bp L`
  - `window` — centered window with `--window` (odd). Use only when speed is critical.

`--report-radius` marks `in_radius`. It does not suppress global rows.

## CLI Examples
Full genomic scan with local triage:  

```Bash
python genesplicer_full.py \
  -i FASTA/genomic \
  -m mutations/combined \
  -g /opt/genesplicer \
  -o out/splice_summary.tsv \
  --pipeline full --report-radius 150
```
Extended context (±2000 bp), global + local:  

```Bash
python genesplicer_full.py \
  -i FASTA/genomic \
  -m mutations/combined \
  -g /opt/genesplicer \
  -o out/splice_summary.ext.tsv \
  --pipeline custom \
  --scan-mode extended \ 
  --context-bp 2000 \
  --report-radius 250
```
Window mode (speed), local radius 150. _**Use with caution**_: results may not be biologically accurate.

_A local triage slice (`--report-radius`) is a review aid, not a filter._
```bash
python genesplicer_full.py \
  -i FASTA/genomic \
  -m mutations/combined \
  -g /opt/genesplicer \
  -o out/splice_summary.win.tsv \
  --pipeline custom \
  --scan-mode window \
  --window 151 \
  --report-radius 150
```
## Thresholds and Policy
- `--visibility-threshold` (default 1.0)
- `--high-cutoff` (default 5.0)
- `--med-cutoff` (default 3.0)
- `--cluster-radius` in bp (default 3)
- `--shift-bp` in bp (default 3)
- `--distance-k` in bp (default 75)

## Performance
- Cost dominated by GeneSplicer scan length: `genome` ≥ `extended` ≫ `window`.
- Cache WT once per gene in `genome` mode; run ALT per SNV only.
- Memory scales with site count; use per-gene writes if millions of rows are expected.

## File Discovery
- FASTA and mapping files discovered via `discover_fasta_files()` and `discover_mapping_files()` returning `{gene: Path}`.
- Gene key derived from filenames (logic in `utility.py`).

## Reproducibility
- Deterministic sorting: summary by `pkey`; events by `pkey, type, cluster_id`; sites by `pkey, allele, type, site_pos, score`.
- Headers written once. Summary appends; events/sites rewritten each run.
- Thresholds are implied by columns; persist CLI in job metadata if needed.

## Assumptions
- Coordinates are absolute in `genome/extended`; window-relative only in `window`.
- Confidence mapped to numeric weights: `{low: 0.5, med: 0.75, high: 1.0}`.
- Model directory expected at `../human` relative to the GeneSplicer binary.
- `should_skip_mutation()` uses `--log` if provided.


# Outputs: Statistical Field Details

## Notation
For a given `pkey`, index clusters by $i = 1..C$.  
$s_{\text{wt},i}$ = top WT score, $s_{\text{mut},i}$ = top MUT score, $\Delta_i = s_{\text{mut},i} - s_{\text{wt},i}$ (Pertea et al., 2001).  
$d_i$ = distance (bp) from the SNV to the nearer of the WT or MUT site in cluster $i$.  
$w_i = e^{- d_i / k}$ with $k = \text{--distance-k}$.  
$H = \text{--high-cutoff}$ (default 5.0), $V = \text{--visibility-threshold}$ (default 1.0), $S = \text{--shift-bp}$ (default 3), $R = \text{--report-radius}$ (if set).

---

## Summary Table fields

### global_sum_weighted_abs_Δ
Weighted aggregate magnitude emphasizing proximity (Tobler, 1970; Cressie, 1993).
- **Definition**: $S_{\text{weighted}} = \sum_i w_i \cdot |\Delta_i|$, with $w_i = e^{- d_i / k}$.
- **Interpretation**: exponential distance decay. Half-weight distance $d_{1/2} = k \ln 2$. Example: $k=75 \Rightarrow d_{1/2} \approx 52$ bp.

### priority (per event)
Sorting/ranking score for events and `top_event_*` (Keeney & Raiffa, 1993).
- **Base**: $|\Delta_i| \cdot e^{- d_i / k}$.
- **Bonuses**: $+2$ if $\text{class} \in \{\text{gained},\text{lost}\}$ and $\max(s_{\text{wt},i}, s_{\text{mut},i}) \ge H$; $+1$ if $\text{class}=\text{shifted}$ and $|dpos| \ge S$.
- **Definition**: $\text{priority}_i = |\Delta_i| e^{- d_i / k} + 2 B_i + 1 L_i$, where $B_i=1$ for high gained/lost, $L_i=1$ for large shift; else $0$.
- **Interpretation**: large, nearby effects rank highest; boosts favor confident gain/loss and sizable shifts.

### pct_delta (per event)
Scale-free effect size normalized by WT strength (Borenstein et al., 2009).
- **Definition**: $r_i = \dfrac{\Delta_i}{\max(|s_{\text{wt},i}|, \varepsilon)}$.
- **Interpretation**: compare across genes with different baselines. Heuristic bands: $\sim 0.1$ small, $\sim 0.3$ moderate, $\ge 0.5$ large.

### conf_weighted_delta (per event)
Confidence-weighted change using certainty weights (Borenstein et al., 2009).
- **Mapping**: $\text{low}\rightarrow 0.5$, $\text{med}\rightarrow 0.75$, $\text{high}\rightarrow 1.0$.
- **Definition**: $\text{CW}\Delta_i = c_{\text{mut},i}\, s_{\text{mut},i} - c_{\text{wt},i}\, s_{\text{wt},i}$.
- **Interpretation**: downweights low-confidence calls; robust alternative to plain $\Delta$.

### nearest_event_bp_any, nearest_event_bp_local
- **Definitions**: 
  - `nearest_any` $= \min_i d_i$
  - `nearest_local` $= \min_{i: d_i \le R} d_i$
- **Interpretation**: proximity diagnostics; small values indicate local impact, large values imply distal cryptic activation (Jaganathan et al., 2019).

### frac_effect_in_radius
Fraction of effect mass within the local radius.
- **Definition**: $\Phi = \dfrac{\sum_{d_i \le R} |\Delta_i|}{\sum_{\text{all } i} |\Delta_i|}$.
- **Interpretation**: localization metric. Near $1$ -> impact concentrated near SNV; near $0$ -> mostly distal.

### global_max_abs_Δscore, local_max_abs_Δscore
- **Definitions**: 
  - `max_global` = $\max_i |\Delta_i|$
  - `max local` = $\max_{i: d_i \le R} |\Delta_i|$
- **Interpretation**: 
  - peak magnitude indicators.

---

## Events Table fields

### Visibility and class guards
Site is **visible** if $\text{score} \ge V$. Cluster class depends on visibility of WT/MUT top sites.
- `gained`: WT absent/invisible; MUT visible.
- `lost`: MUT absent/invisible; WT visible.
- `shifted`: both visible and $|dpos| \ge S$.
- `strengthened / weakened`: same/near position with $|\Delta| \ge 1.0$.
- `none`: otherwise.

### rank_wt, rank_mut
- **Definition**: rank within `(allele, type)` by descending score; ties resolved deterministically.
- **Interpretation**: relative strength among all donors or all acceptors per allele.

### pct_delta, conf_weighted_delta, priority
- Same formulas as above, computed per cluster.

---

## Sites Table fields

### distance_to_snv
- **Definition**: `distance = |site_pos - snv_pos|`
- **Interpretation**: feeds distance kernel $w = e^{- d / k}$ used by weighted metrics and `priority` (Tobler, 1970; Cressie, 1993).

### visible_flag
- **Definition**: indicator $1$ if $\text{score} \ge V$, else $0$.
- **Interpretation**: governs inclusion in class logic (e.g., gained/lost) and counts.

---

## Design levers

### Distance kernel k (`--distance-k`)
Controls decay in `global_sum_weighted_abs_Δ` and `priority`. Larger $k$ spreads weight to distal events; smaller $k$ concentrates weight locally (Tobler, 1970; Cressie, 1993).

### Shift threshold S (`--shift-bp`)
Higher $S$ demands larger positional moves for `shifted`; lower $S$ captures micro-shifts.

### Visibility V (`--visibility-threshold`) and high cutoff H (`--high-cutoff`)
$V$ prunes low-signal sites; raising it reduces false positives but may miss weak true sites.  
$H$ gates “high-confidence” semantics in counts and `priority` bonuses.

---

## Quick interpretation patterns
- Local disruption: high `max_local` and small `nearest_local`; $\Phi$ high.
- Cryptic activation: `gained` with large $d_i$, high `priority` if $s_{\text{mut},i} \ge H$ (Jaganathan et al., 2019).
- Global silence: `n_clusters=0` or `max_global \approx 0` and $S_{\text{weighted}} \approx 0$.

## Citation

If you use this pipeline, please cite:

- BioFeatureFactory: Goldmintz J. BioFeatureFactory: Python toolkit for automated biological feature extraction. GitHub repository. https://github.com/jgoldmintz/BioFeatureFactory
- Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). *Introduction to Meta-Analysis.* Chichester, UK: John Wiley & Sons. https://doi.org/10.1002/9780470743386  
  Explains precision and reliability weighting, forming the rationale for confidence-weighted deltas.
- Cressie, N. (1993, rev. 2015). *Statistics for Spatial Data.* Revised edition. Hoboken, NJ: John Wiley & Sons. https://doi.org/10.1002/9781119115151  
  Provides the mathematical framework for exponential distance weighting and spatial decay kernels.
- Jaganathan, K., Kyriazopoulou Panagiotopoulou, S., McRae, J. F., Darbandi, S. F., Knowles, D., Li, Y. I., Kosmicki, J. A., Arbelaez, J., Cui, W., Schwartz, G. B., Chow, E. D., Kanterakis, E., Gao, H., Kia, A., Batzoglou, S., Sanders, S. J., & Farh, K. K.-H. (2019). *Predicting Splicing from Primary Sequence with Deep Learning.* *Cell*, 176(3), 535–548.e24. https://doi.org/10.1016/j.cell.2018.12.015  
  Demonstrates long-range splice-context modeling supporting the pipeline’s distance-aware design.
- Keeney, R. L., & Raiffa, H. (1993). *Decisions with Multiple Objectives: Preferences and Value Tradeoffs.* New York, NY: John Wiley & Sons. https://doi.org/10.1002/9780470611876.ch15  
  Forms the theoretical basis for constructing the weighted-utility “priority” metric.
- Pertea, M., Lin, X., & Salzberg, S. L. (2001). *GeneSplicer: A New Computational Method for Splice Site Prediction.* *Nucleic Acids Research*, 29(5), 1185–1190. https://doi.org/10.1093/nar/29.5.1185  
  Describes the core algorithm underlying splice-site scoring in this pipeline.
- Tobler, W. R. (1970). *A Computer Movie Simulating Urban Growth in the Detroit Region.* *Economic Geography*, 46(Suppl.), 234–240. https://doi.org/10.2307/143141  
  Introduces the distance-decay principle, later generalized as a modeling concept for spatial dependency.



## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
