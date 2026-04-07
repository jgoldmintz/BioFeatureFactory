 # NetPhos Pipeline

  High-throughput kinase-site prediction for WT and mutant protein FASTAs using NetPhos 3.1/APE.

  ## 1. Prerequisites

  | Component | Notes |
  |-----------|-------|
  | NetPhos 3.1 + APE | Download from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetPhos-3.1/) (academic license). Requires `tcsh`. Set path via `--native-ape-path` or `NETPHOS_APE_PATH`/`NETPHOS_HOME` env var. |
  | Python $\geq$ 3.9 | With dependencies in `../utils/`. |
  | Mutation mapping CSVs | `{GENE}`-specific CSVs used for pkey generation (`mutations/combined/aa/` style). |
  | FASTA inputs | WT and/or mutant AA FASTAs, one gene per file. Extensions: `.fasta`, `.fa`, `.fas`, `.fna`. |

  ## 2. Required Inputs

  | Path / Flag | Description |
  |-------------|-------------|
  | `--mapping-dir` | Directory containing mutation CSVs named like `*{GENE}*.csv` (case-insensitive). |
  | `INPUT` arg | FASTA file/dir. |
  | `OUTPUT` arg | Output directory for TSV results. |
  | Optional logs | `--log` validation log to skip failed mutations in mutant mode. |

  ## 3. Quick Start Recipes

  ```bash
  # Full pipeline, WT sequences (auto-detect native or Docker)
  python3 netphos_pipeline.py \
      --mapping-dir mutations/combined/aa/ \
      wt/aaseq/ \
      results/
  # Writes per gene: results/{GENE}/NetPhos/{GENE}.{tsv,events.tsv,sites.tsv}
```
```bash
  # Mutant run with validation filtering
  python3 netphos_pipeline.py \
      --log validation.log \
      --mapping-dir mutations/combined/aa/ \
      mut/aaseq/ \
      results/
  # Writes per gene: results/{GENE}/NetPhos/{GENE}.{tsv,events.tsv,sites.tsv}
```
  ## 4. Execution Backend Controls

  | Option | Purpose |
  |--------|---------|
  | `--native-ape-path /opt/netphos/ape-1.0/ape` | Explicit binary path. |
  | Env vars | `NETPHOS_APE_PATH`, `NETPHOS_HOME`, and `NETNGLYC_HOME` influence detection (searched in that order). |

  ## 5. Processing Controls

  - `--batch-size N` and `--timeout SEC` -- tune long FASTA runs. Batch strategy is chosen automatically by sequence count: 1 seq $\rightarrow$ single run; 2-10 $\rightarrow$ batch of 10; 11-100 $\rightarrow$ batch of 25; >100 $\rightarrow$ batch of 50.
  - `--threshold FLOAT` -- filter by score (default: 0.5). Sites with score below threshold are excluded.
  - `--yes-only` -- sets `--threshold` to 0.0 (disables score filtering entirely); does not filter for YES-answered predictions.
  - `--no-cache` / `--clear-cache` -- manage `~/.netphos_cache`.
  - `--wt-header HEADER` -- FASTA header used to identify the WT sequence (default: `ORF`).
  - `--keep-intermediates` -- retain intermediate files for debugging.
  - `--verbose` -- enable verbose output.

  ## 6. Outputs

  The pipeline produces three TSV files per gene.

  Output structure:
  ```
  {output}/
    {GENE}/
      NetPhos/
        {GENE}.tsv          -- per-mutation summary
        {GENE}.events.tsv   -- per-site, per-kinase events
        {GENE}.sites.tsv    -- raw WT and MUT predictions
  ```

  ### 6.1 Summary TSV (`{GENE}.tsv`)

  Per-mutation summary with phosphorylation site changes.

  | Column | Description | Units |
  |--------|-------------|-------|
  | `pkey` | `{GENE}-{MUTATION}` primary key | string |
  | `Gene` | Gene symbol | string |
  | `n_sites_wt`, `n_sites_mut` | Count of sites with score $\geq$ `--threshold` | count |
  | `count_gained`, `count_lost`, `count_strengthened`, `count_weakened`, `count_stable` | Classification tallies | count |
  | `max_abs_delta`, `sum_abs_delta` | Maximum and sum of absolute score changes | probability (0-1) |
  | `n_kinases_affected` | Number of kinases with any event | count |
  | `top_event_type` | Dominant event (`gained`, `lost`, etc.) | categorical |
  | `top_event_delta` | $\Delta$ score for dominant event | probability (0-1) |
  | `top_event_position` | Amino-acid index of dominant event | residue index |
  | `top_event_kinase` | Kinase for dominant event | string |
  | `top_event_classification_code` | Numeric encoding of dominant event | integer |
  | `qc_flags` | Quality control flags | string |

  ### 6.2 Events TSV (`{GENE}.events.tsv`)

  Per-mutation, per-site, per-kinase events.

  | Column | Description | Units |
  |--------|-------------|-------|
  | `pkey` | `{GENE}-{MUTATION}` primary key | string |
  | `Gene` | Gene symbol | string |
  | `position` | 1-based amino-acid position | residue index |
  | `amino_acid_wt`, `amino_acid_mut` | WT and MUT residues | single-letter |
  | `kinase` | Kinase motif classifier | string |
  | `wt_score`, `mut_score` | NetPhos phosphorylation propensity | probability (0-1) |
  | `delta` | $\Delta$ score ($\text{MUT} - \text{WT}$) | probability (0-1) |
  | `wt_answer`, `mut_answer` | NetPhos YES/NO flag | categorical |
  | `classification` | gained/lost/strengthened/weakened/stable/subthreshold | categorical |
  | `classification_code` | Numeric encoding (gained=2, lost=-2, strengthened=1, weakened=-1, stable=0, subthreshold=-3) | integer |

  ### 6.3 Sites TSV (`{GENE}.sites.tsv`)

  Raw WT and MUT predictions for all residues.

  | Column | Description | Units / Format |
  |--------|-------------|----------------|
  | `pkey` | `{GENE}-{MUTATION}` primary key | string |
  | `Gene` | Gene symbol | string |
  | `allele` | Sequence label (wt or mutant ID) | string |
  | `seq_name` | Full sequence name from NetPhos output | string |
  | `position` | 1-based amino-acid position | amino-acid index |
  | `amino_acid` | Residue at position (S/T/Y) | single-letter code |
  | `context` | Flanking peptide window reported by NetPhos | string |
  | `score` | NetPhos phosphorylation propensity | probability (0-1) |
  | `kinase` | Kinase motif classifier (e.g., CKII, cdc2, unsp) | string |
  | `answer` | NetPhos YES/NO confidence flag | categorical |

  Classification codes and $|\Delta|$ threshold: sites are **gained** if score crosses `--threshold` only in MUT, **lost** if the reverse, and **strengthened/weakened** when both are above threshold but $|\Delta| > 0.05$.

  ## 7. Troubleshooting

  | Symptom | Resolution |
  |---------|------------|
  | tcsh: Command not found | Install tcsh (required by APE). |
  | "No mapping file found" | Confirm --mapping-dir contains files named like *{GENE}*.csv. |
  | APE binary not found | Ensure the binary is executable and provide --native-ape-path or set NETPHOS_APE_PATH. |
  | Cache confusion | Run python3 netphos_pipeline.py --clear-cache. |

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.