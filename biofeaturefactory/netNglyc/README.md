# NetNGlyc Pipeline

N-linked glycosylation site prediction for WT and mutant protein sequences using NetNGlyc 1.0, with host-side SignalP 6.0 for signal-peptide context.

## Requirements

| Component | Notes |
|-----------|-------|
| NetNGlyc 1.0 | [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/), academic license. Point at it with `-nnb/--native-netnglyc-bin`, or set `NETNGLYC_PATH` / `NETNGLYC_HOME`. |
| SignalP 6.0 (fast) | [DTU Health Tech](https://services.healthtech.dtu.dk/services/SignalP-6.0/), academic license. Requires `numpy<2`. Point at it with `-snp/--signalp6-bin`, or leave it on `PATH`. Cache defaults to `~/.signalp6_cache`; `-cd/--cache-dir` relocates it. |
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. |
| Input FASTAs | ORF FASTAs (header `>ORF`). WT amino-acid sequences are synthesized automatically. |
| Mapping CSVs | Per-gene NT -> AA mappings (`mutant`, `aamutant`) from `variant_mapping`. |

### SignalP 6 adapter

NetNGlyc's tcsh wrapper expects a SignalP v3/v4 binary at `$SIGNALP`. Point it at the shim instead:

```tcsh
setenv SIGNALP /path/to/BioFeatureFactory/biofeaturefactory/netNglyc/bin/signalp6_adapter
```

The pipeline runs SignalP 6 itself, exports `SIGNALP6_RESULTS_DIR`, and the adapter re-emits those
results in the legacy 14-column format. Without it NetNGlyc still runs, but with no signal-peptide context.

## Usage

```bash
# Directory mode: variant_mapping output root
python netnglyc_pipeline.py -i out/ -o results/

# Flat FASTA directory + explicit mapping directory
python netnglyc_pipeline.py -i FASTA_files/nt/ -o results/ -md mutations/aa/ -l validation.log
```

In directory mode `-i` is the `variant_mapping` output root (`<root>/<GENE>/fastas/` and
`<root>/<GENE>/mappings/`); the gene is taken from the directory name. `input` and `output`
are also accepted positionally.

## Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | -- | variant_mapping output root, flat FASTA directory, or single FASTA |
| `-o, --output` | -- | Output base directory |
| `-md, --mapping-dir` | -- | Mutation mapping CSV directory (required for parsing modes) |
| `-l, --log` | -- | Validation log file or directory; skips failed mutations |
| `-th, --threshold` | `0.5` | Minimum glycosylation potential |
| `-w, --workers` | `4` | Parallel workers |
| `-bt, --batch-timeout` | `5000` | NetNGlyc execution timeout (seconds) |
| `-cd, --cache-dir` | -- | Cache directory for SignalP/NetNGlyc results |
| `-cc, --clear-cache` | off | Clear all cached results and exit |
| `-nnb, --native-netnglyc-bin` | -- | NetNGlyc binary, or the install directory containing it |
| `-snp, --signalp6-bin` | `signalp6` on PATH | signalp6 executable, or its install directory |
| `-v, --verbose` | off | Verbose output |

## Output

```
{output}/{GENE}/NetNglyc/
    {GENE}.tsv          -- per-mutation summary
    {GENE}.events.tsv   -- per-position classified deltas
    {GENE}.sites.tsv    -- raw WT and MUT predictions
```

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `n_sites_wt`, `n_sites_mut` | Sites at or above `--threshold` |
| `count_gained`, `count_lost`, `count_strengthened`, `count_weakened`, `count_stable` | Classification tallies |
| `max_abs_delta`, `sum_abs_delta` | Max and sum of absolute potential changes |
| `top_event_type` | Dominant event label |
| `top_event_classification_code` | Numeric encoding (gained=2, lost=-2, ...) |
| `top_event_delta` | Potential change for the dominant event |
| `top_event_position` | Residue index of the dominant event |
| `wt_signalp_has_signal`, `mut_signalp_has_signal` | 1 if a signal peptide is predicted |
| `wt_signalp_probability`, `mut_signalp_probability` | SignalP cleavage probability |
| `wt_signalp_cleavage`, `mut_signalp_cleavage` | SignalP cleavage site |
| `frac_effect_post_cleavage` | Fraction of total absolute delta downstream of the cleavage site |
| `qc_flags` | `missing_wt`, `missing_mut`, `no_delta`, `no_signalp` |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `classification` | gained / lost / strengthened / weakened / stable / subthreshold |
| `classification_code` | gained=2, lost=-2, strengthened=1, weakened=-1, stable=0, subthreshold=-1 |
| `wt_potential`, `mut_potential`, `delta` | WT score, MUT score, MUT - WT |
| `wt_sequon`, `mut_sequon` | Motif strings (e.g. `NKSE`) |
| `wt_above_threshold`, `mut_above_threshold`, `post_cleavage` | 0/1 |
| `position` | Residue index of the motif |

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene` | Identifiers |
| `allele` | Sequence label (WT or mutant ID) |
| `seq_name` | Sequence name from NetNGlyc output |
| `position` | Residue index of the sequon |
| `sequon` | N-X-S/T sequon string |
| `potential` | NetNGlyc glycosylation potential |
| `jury_agreement` | Raw `(votes/total)` string |
| `jury_agreement_score` | Parsed votes/total |
| `n_glyc_result`, `n_glyc_result_code` | NetNGlyc symbol and code (`+++`=3 ... `---`=-3) |
| `signalp_has_signal`, `signalp_probability`, `signalp_cleavage` | SignalP prediction |
| `above_threshold` | 0/1 |

Classification: **gained** if the potential crosses `--threshold` only in MUT, **lost** if only in WT,
**strengthened**/**weakened** when both are above threshold and `|delta| > 0.05`.

## Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `ProcessPoolExecutor ... SC_SEM_NSEMS_MAX` on macOS | Run with `-w 1`; macOS limits POSIX semaphores |
| WT/MUT reported missing after a completed run | Confirm `netnglyc_outputs_*` with `wt/` and `mut/` exists; FASTA inputs are not parse targets |
| SignalP columns blank | `-cd` must contain `*_sp6_output/prediction_results.txt`, or omit it to use `~/.signalp6_cache` |
| Mutants flagged `missing_mut` | NetNGlyc emitted "No sites predicted". Lower `-th` if weaker signals are expected |
| NetNGlyc binary not found | Pass `-nnb`, or set `NETNGLYC_PATH` / `NETNGLYC_HOME` |
| SignalP integration inert | `$SIGNALP` in the netNglyc tcsh script must point at `signalp6_adapter`, and it must be executable |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
