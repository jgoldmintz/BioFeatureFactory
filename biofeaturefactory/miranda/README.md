# Miranda Comparative miRNA Binding Pipeline

Runs miRanda on WT and mutant transcripts and reports per-miRNA, per-locus, and per-segment delta metrics for miRNA-target binding. WT is scanned once per gene and reused across that gene's mutations.

## Requirements

| Component | Notes |
|-----------|-------|
| miRanda | `conda install -c bioconda miranda` (binary on `PATH`), or pass a directory with `-m`. |
| Python >= 3.9 | With `pandas`. Uses `biofeaturefactory/lib/`. |
| miRNA database | FASTA of mature miRNAs, passed with `-d`. |
| Input FASTAs | Transcript FASTAs. |
| Mapping CSVs | Transcript mappings; intronic variants also need the pre-mRNA mappings. |

## Usage

```bash
# Directory mode: variant_mapping output root
python miranda_ensemble.py -i out/ -o miranda_out/ -d db/hsa_miRNA.fa

# File mode with explicit mappings and a local miRanda build
python miranda_ensemble.py \
    -i FASTA_files/transcript/ -o miranda_out/ \
    -m /opt/miranda -d db/hsa_miRNA.fa \
    -M mappings/transcript/ -ipm mappings/intron_premRNA/ \
    -l logs/validation --max-workers 8
```

In directory mode the gene is taken from the **directory** name, not the filename; a flat directory
of FASTAs also works. `-M` and `-ipm` are derived from `-i` in directory mode, since `mappings/`
and `fastas/` are siblings under the same `<GENE>/`. Intronic variants are scanned against both
the intron record and the pre-mRNA record.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i, --input` | Yes | -- | variant_mapping output root, or a single WT FASTA |
| `-o, --output` | Yes | -- | Output base directory |
| `-d, --mirna_db` | Yes | -- | miRNA database FASTA |
| `-m, --miranda_dir` | No | PATH | Directory containing the miRanda executable |
| `-M, --variant-mapping-root`, `--mapping-dir` | No | derived | File mode only; transcript mapping CSV/TSV |
| `-ipm, --intron-premrna-mapping` | No | derived | File mode only; `intron_premRNA_mapping_<GENE>.csv` |
| `-si, --strict-introns` | No | off | Pass miRanda `-strict` for the intron and pre-mRNA substrates only; WT and MUT always use the same setting |
| `-np, --no-parallel` | No | off | Disable multiprocessing |
| `--max-workers` | No | auto | Cap on parallel workers |
| `-l, --log` | No | -- | Validation log file or directory |
| `-wh, --wt-header` | No | `transcript` | Preferred WT FASTA header |

## Output

```
{output}/{GENE}/Miranda/
    {GENE}.tsv          -- per-mutation summary
    {GENE}.events.tsv   -- per-miRNA, per-locus comparison
    {GENE}.sites.tsv    -- raw per-hit audit rows
```

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `n_hits_wt`, `n_hits_mut` | Visible binding sites per allele (score >= visibility threshold) |
| `n_mirna` | Distinct miRNAs with at least one site |
| `n_loci` | Per-miRNA proximity-based binding clusters |
| `n_segments` | Merged cross-miRNA binding windows |
| `n_competitive_segments` | Segments with >= 2 miRNAs on overlapping regions |
| `n_new_competitors` | Segments where the mutation added competing miRNAs |
| `global_count_gained_high`, `global_count_lost_high` | High-scoring sites gained / lost |
| `global_count_shifted` | Sites shifted by at least the shift threshold |
| `global_max_abs_delta_score` | Largest absolute score change |
| `global_sum_weighted_abs_delta` | Distance-weighted sum of absolute score changes |
| `nearest_event_bp_any` | Distance from the variant to the nearest event (bp) |
| `local_count_gained_high`, `local_count_lost_high`, `local_count_shifted` | As above, within the report radius |
| `local_max_abs_delta_score` | Largest absolute change within the radius |
| `nearest_event_bp_local` | Nearest event within the radius (bp) |
| `frac_effect_in_radius` | Fraction of total perturbation inside the radius (0-1) |
| `max_segment_abs_delta_best` | Largest absolute change in the strongest segment |
| `frac_effect_in_competitive_segments` | Fraction of total change in multi-miRNA segments (0-1) |
| `top_event_mirna` | miRNA of the top-ranked event |
| `top_event_class` | gained / lost / shifted / strengthened / weakened / none |
| `top_event_delta_score` | Score change of the top-ranked event |
| `top_event_pos` | Position of the top-ranked event (bp) |
| `qc_flags` | Semicolon-separated (e.g. `no_hits;low_signal_only`) |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `mirna_id` | miRNA identifier |
| `locus_id` | Per-miRNA cluster ID (`m1`, `m2`, ...) |
| `segment_id` | Cross-miRNA segment ID (`s1`, `s2`, ...) |
| `wt_pos`, `mut_pos` | Highest-scoring site position per allele (bp) |
| `dpos` | `mut_pos - wt_pos` (bp) |
| `wt_tot_score`, `mut_tot_score`, `delta_tot_score` | miRanda scores and MUT - WT |
| `pct_delta` | `delta / max(mut, wt)` |
| `wt_tot_energy`, `mut_tot_energy`, `delta_energy` | Binding energies (kcal/mol) and MUT - WT |
| `distance_to_snv` | Distance from the variant to the event (bp) |
| `rank_wt`, `rank_mut` | Rank among that allele's hits (1 = strongest) |
| `conf_wt`, `conf_mut` | 0.5 below the visibility threshold, 1.0 at or above |
| `conf_weighted_delta` | `conf_mut * mut_score - conf_wt * wt_score` |
| `cls` | gained / lost / shifted / strengthened / weakened / none |
| `is_high_impact` | 0/1 |
| `priority` | Ranking score: distance-decayed absolute delta plus a class bonus |
| `in_radius` | 0/1 |

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `allele` | WT or MUT |
| `mirna_id` | miRNA name from the `Read Sequence:` header |
| `site_pos` | Binding position (bp) |
| `tot_score`, `tot_energy` | miRanda total score and energy (kcal/mol) |
| `max_score`, `max_energy` | Best local alignment score and its energy |
| `strand` | + / - |
| `len_mirna`, `len_target` | Aligned lengths (nt) |
| `visibility_flag` | 1 if score >= visibility threshold |
| `distance_to_snv` | Distance to the variant (bp) |
| `locus_id`, `segment_id` | Links to the events table |
| `parser_confidence` | Parser extraction confidence (0.5-1.0) |
| `run_meta` | `{gene}-wt` or `{gene}-mut` |

## Thresholds

Defined as constants at the top of `miranda_ensemble.py`; not CLI arguments.

| Constant | Value | Purpose |
|----------|-------|---------|
| `VISIBILITY_THRESHOLD` | 140.0 | Minimum score for a site to count as visible |
| `HIGH_CUTOFF` | 150.0 | Defines a high-impact event |
| `REPORT_RADIUS` | 40 bp | Local radius for `frac_effect_in_radius` |
| `DISTANCE_K` | 25.0 bp | Decay constant for distance weighting |
| `MERGE_WINDOW_NT` | 15 bp | Intra-miRNA cluster merge window |
| `SEGMENT_WINDOW_NT` | 25 bp | Inter-miRNA segment merge window |
| `SHIFT_NT` | 4 bp | Positional shift threshold |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.

Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
