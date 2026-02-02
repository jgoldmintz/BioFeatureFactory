 # NetPhos Pipeline

  High-throughput kinase-site prediction for WT and mutant protein FASTAs using NetPhos 3.1/APE.

  ## 1. Prerequisites

  | Component | Notes |
  |-----------|-------|
  | NetPhos 3.1 + APE | Install NetPhos, ensure `tcsh` exists, and point to the binary via `--native-ape-path` or `NETPHOS_APE_PATH/NETPHOS_HOME`. |
  | Python ≥3.9 | With dependencies in `../utils/`. |
  | Mutation mapping CSVs | `{GENE}`-specific CSVs used for pkey generation (`mutations/combined/aa/` style). |
  | FASTA inputs | WT and/or mutant AA FASTAs, one gene per file. Extensions: `.fasta`, `.fa`, `.fas`, `.fna`. |

  ## 2. Required Inputs

  | Path / Flag | Description |
  |-------------|-------------|
  | `--mapping-dir` | Directory containing mutation CSVs named like `*{GENE}*.csv` (case-insensitive). Required for every mode except `netphos-only`. |
  | `INPUT` arg | FASTA file/dir (full or NetPhos-only modes) or NetPhos output file/dir (parse-only). |
  | `OUTPUT` arg | TSV (parse/full modes) or `.out` (NetPhos-only). |
  | Optional logs | `--log` validation log to skip failed mutations in mutant mode. |

  ## 3. Quick Start Recipes

  ```bash
  # Full pipeline, WT sequences (auto-detect native or Docker)
  python3 netphos-pipeline.py --mode full-pipeline \
      --mapping-dir mutations/combined/aa/ \
      wt/aaseq/ \
      results.tsv
```
```bash
  # Mutant run with validation filtering
  python3 netphos-pipeline.py --mode full-pipeline \
      --is-mutant \
      --log validation.log \
      --mapping-dir mutations/combined/aa/ \
      mut/aaseq/ \
      mutant_results.tsv
```
```bash
  # Parse existing NetPhos outputs
  python3 netphos-pipeline.py --mode parse-only \
      --mapping-dir mutations/combined/aa/ \
      netphos_outputs/ \
      parsed.tsv
```
```bash
  # Run NetPhos only (no parsing)
  python3 netphos-pipeline.py --mode netphos-only \
      wt/ABCB1.fa \
      ABCB1-netphos.out
```
  ## 4. Execution Backend Controls

  | Option | Purpose |
  |--------|---------|
  | --native-ape-path /opt/netphos/ape-1.0/ape | Explicit binary path. |
  | Env vars | NETPHOS_APE_PATH, NETPHOS_HOME influence detection. |

  ## 5. Processing Controls

  - --mode {full-pipeline, netphos-only, parse-only}
  - --batch-size N and --timeout SEC to tune long FASTA runs.
  - --yes-only or --threshold to filter by score.
  - --no-cache / --clear-cache manage ~/.netphos_cache.
  - --is-mutant toggles mutant vs WT parsing logic (affects amino-acid matching, validation filtering).

  ## 6. Outputs

  - Full/parse modes produce TSVs with pkey, gene, site, kinase, score, answer.
  - NetPhos-only mode writes raw .out files per FASTA; rerun parse-only later.
  - Every matched mutation receives a pkey = {GENE}-{MUT_ID} entry; duplicates per position are tracked (1:many mapping).
  
  ### NetPhos Output Table

  | Column | Description | Units / Format |
  |--------|-------------|----------------|
  | pkey | Unique mutation identifier composed as {GENE}-{MUTATION_ID} | string |
  | Gene | NetPhos-reported sequence header (typically <GENE>_aa) | string |
  | pos | 1-based amino-acid position of the predicted phosphorylation site | amino-acid index |
  | amino_acid | Residue at pos in the evaluated sequence (S/T/Y) | single-letter code |
  | context | Flanking peptide window reported by NetPhos | string |
  | score | NetPhos phosphorylation propensity | probability (0–1 float) |
  | kinase | Kinase motif classifier (e.g., CKII, cdc2, unsp) | string |
  | answer | NetPhos’s YES/NO confidence flag | categorical |

  ## 7. Troubleshooting

  | Symptom | Resolution |
  |---------|------------|
  | tcsh: Command not found | Install tcsh (required by APE). |
  | "No mapping file found" | Confirm --mapping-dir contains files named like *{GENE}*.csv. |
  | APE binary not found | Ensure the binary is executable and provide --native-ape-path or set NETPHOS_APE_PATH. |
  | Cache confusion | Run python3 netphos-pipeline.py --clear-cache. |

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../LICENSE) file in the root BioFeatureFactory directory for details.