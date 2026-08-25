# GeneSplicer Pipeline

Runs GeneSplicer on WT and ALT genomic sequence, clusters the predicted donor/acceptor sites, pairs them across alleles, and emits per-variant delta calls.

## Requirements

| Component | Notes |
|-----------|-------|
| GeneSplicer | `conda install -c bioconda genesplicer` (binary on `PATH`), or pass its directory with `-g`. |
| Model directory | One containing `config_file`. Resolved from `-g` or the binary location; conda ships it as `share/genesplicer-*/human`. Override with `--model-dir`. |
| Python >= 3.9 | With `pandas`, `numpy`. Uses `biofeaturefactory/lib/`. |
| Input FASTAs | Genomic FASTAs. Sequence named `genomic` preferred, else the first record. |
| Mapping CSVs | Per-gene genomic mappings. Mutation column auto-detected from `genomic`, `mutant`, `mutation`, `nt_mutation`, `ntmutant`. |

## Usage

```bash
# Directory mode: variant_mapping output root
python genesplicer_ensemble.py -i out/ -o results/

# File mode with explicit mapping and a local build
python genesplicer_ensemble.py \
    -i FASTA/genomic -o results/ \
    -g /opt/genesplicer -m mutations/gDNA/ \
    -rr 150 -l validation.log
```

In directory mode the gene is taken from the **directory** name, not the filename; a flat directory
of FASTAs also works. `-m` is derived from `-i` in directory mode, since `mappings/gDNA/` and
`fastas/` are siblings under the same `<GENE>/`.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i, --input` | Yes | -- | variant_mapping output root, or a single genomic FASTA |
| `-o, --output` | Yes | -- | Output base directory |
| `-m, --variant-mapping-root`, `--mapping-dir` | No | derived | File mode only; genomic mapping CSV |
| `-g, --genesplicer-dir` | No | PATH | Directory containing the GeneSplicer binary |
| `--model-dir` | No | derived | Absolute path to a model directory containing `config_file` |
| `-w, --window` | No | `151` | Window / reporting size |
| `-rr, --report-radius` | No | `--window` | Radius (bp) counted as local; sets the `in_radius` column |
| `-dk, --distance-k` | No | `75` | Distance-decay constant (bp) |
| `-vt, --visibility-threshold` | No | `1.0` | Minimum score for a site to be visible |
| `-hc, --high-cutoff` | No | `5.0` | Score cutoff for high-confidence gained/lost |
| `-sb, --shift-bp` | No | `3` | Minimum bp difference to classify a site as shifted |
| `-cr, --cluster-radius` | No | `3` | Single-linkage radius (bp) for grouping sites into clusters |
| `--max-cluster-span` | No | `0` | Max total span (bp) of one cluster. `0` = unbounded and inert |
| `-rmb, --redirect-max-bp` | No | `100` | Max bp between a lost and a gained site to link them as one redirected site |
| `-rmbp, --register-max-bp` | No | `12` | Linked pairs within this distance are labelled `shifted` rather than `redirected` |
| `-rsmb, --redirect-snv-max-bp` | No | `80` | Max distance from the variant span to the lost site for a redirect to count as variant-caused |
| `-wo, --workers` | No | half cores, max 8 | Parallel workers |
| `-l, --log` | No | -- | Validation log file or directory |

`--report-radius` only marks rows; it does not filter them.

## Output

```
{output}/{GENE}/GeneSplicer/
    {GENE}.tsv          -- one row per variant
    {GENE}.events.tsv   -- per-cluster donor/acceptor events
    {GENE}.sites.tsv    -- per-allele, per-site audit rows
```

Coordinates are absolute genomic positions. Distances are bp from the variant. Scores are raw
GeneSplicer scores. Confidence maps `low=0.5`, `med=0.75`, `high=1.0`.

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `n_sites_wt`, `n_sites_mut` | Visible sites per allele (score >= `-vt`) |
| `n_clusters` | Donor/acceptor clusters evaluated |
| `global_count_gained_high`, `global_count_lost_high` | Gained / lost clusters where `max(wt, mut) >= -hc` |
| `global_count_shifted` | Clusters with `|dpos| >= -sb` |
| `global_max_abs_deltascore` | Largest absolute score change |
| `global_sum_weighted_abs_delta` | Distance-weighted sum of absolute score changes |
| `nearest_event_bp_any` | Distance to the nearest event (bp) |
| `local_count_gained_high`, `local_count_lost_high`, `local_count_shifted` | As above, within `-rr` |
| `local_max_abs_deltascore` | Largest absolute change within `-rr` |
| `nearest_event_bp_local` | Nearest event within `-rr` (bp) |
| `frac_effect_in_radius` | Fraction of total change inside `-rr` (0-1) |
| `top_event_type` | Class of the highest-priority cluster |
| `top_event_deltascore` | Score change of that cluster |
| `top_event_pos` | Representative position (MUT if present, else WT) |
| `dominant_boundary` | `donor` or `acceptor`, whichever carries the larger summed absolute change |
| `qc_flags` | `no_sites`, `far_event>2kb`, `low_signal_only` |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `type` | `donor` or `acceptor` |
| `cluster_id` | Stable ID within `(pkey, type)` (e.g. `d1`, `a2`) |
| `wt_pos`, `mut_pos` | Top site position per allele in the cluster |
| `dpos` | `mut_pos - wt_pos` |
| `wt_score`, `mut_score`, `dscore` | Top scores and MUT - WT |
| `pct_delta` | `dscore / max(|wt_score|, eps)` |
| `distance_to_snv` | Nearest cluster site to the variant (bp) |
| `rank_wt`, `rank_mut` | Rank among that allele's sites of this type (1 = strongest) |
| `conf_wt`, `conf_mut` | 0.5 / 0.75 / 1.0 |
| `conf_weighted_delta` | `conf_mut * mut_score - conf_wt * wt_score` |
| `cls` | gained / lost / shifted / strengthened / weakened / none |
| `is_high_impact` | 0/1 |
| `priority` | Distance-decayed absolute delta plus class bonuses |
| `in_radius` | 0/1 |

Event classes: **gained** = WT absent or invisible, MUT visible. **lost** = the reverse.
**shifted** = both visible and `|dpos| >= -sb`. **strengthened**/**weakened** = same or near
position with `|dscore| >= 1.0`. **none** = otherwise.

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `allele` | WT or MUT |
| `type` | `donor` or `acceptor` |
| `site_pos` | Site position (End5 for donors, End3 for acceptors) |
| `score` | GeneSplicer score |
| `confidence` | 0.5 / 0.75 / 1.0 |
| `rank` | Rank within `(allele, type)` (1 = strongest) |
| `distance_to_snv` | `|site_pos - variant_pos|` (bp) |
| `visible_flag` | 1 if `score >= -vt` |
| `cluster_id` | Links to the events table |
| `in_radius` | 0/1 |

Consistency checks: `n_sites_wt` equals the count of `visible_flag=1, allele=WT` rows for that
`pkey`; likewise for MUT. `n_clusters` is at least the number of unique `(type, cluster_id)` pairs
in `events` for that `pkey`.

## Notes

- Sorting is deterministic: summary by `pkey`; events by `pkey, type, cluster_id`; sites by `pkey, allele, type, site_pos, score`.
- Headers are written once. The summary appends; events and sites are rewritten each run.
- Thresholds are implied by the columns, not stored in the output. Persist the CLI in job metadata if you need it.

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.

Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
