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

### Install SignalP 6 fast

Download and extract the licensed fast package from DTU. Install it in a separate
conda environment: its `torch<2` requirement must stay isolated from BFF's newer
PyTorch dependencies. Bootstrap only checks for SignalP; it does not install it.
Set `SIGNALP_PACKAGE` to the extracted `signalp-6-package` directory containing
`setup.py` and `models/`:

```bash
conda create -n signalp6 python=3.10 -y
conda activate signalp6
SIGNALP_PACKAGE="/path/to/signalp6_fast/signalp-6-package"
python -m pip install "numpy<2" "$SIGNALP_PACKAGE"
SIGNALP_DIR=$(python -c 'import os, signalp; print(os.path.dirname(signalp.__file__))')
mkdir -p "$SIGNALP_DIR/model_weights"
cp "$SIGNALP_PACKAGE/models/distilled_model_signalp6.pt" "$SIGNALP_DIR/model_weights/"
signalp6 --version
```

After installation, reactivate your **BFF environment** and run the NetNGlyc
pipeline there. For example, use `conda activate bff` if your BFF environment is
named `bff`; that name is not required.

The pipeline does **not** activate the SignalP conda environment. It invokes the
installed `signalp6` executable as a subprocess; the executable's shebang selects
that environment's Python and dependencies. Leave the executable and weights in
the SignalP environment; do not move them into NetNGlyc's `bin/` directory.

The resolver searches sibling environments named `signalp6` or `signalp6_fast`,
as well as `PATH` and configured locations. For a different name or location,
pass `--signalp6-bin /path/to/envs/your-signalp-env/bin/signalp6` to the NetNGlyc
pipeline. This flag selects an existing executable, does not install anything,
and is **not a bootstrap flag**. The `signalp6_adapter` below is a separate shim,
not the SignalP installation.

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

SignalP execution failures, timeouts, missing rows and malformed predictions are
failures, not negative signal-peptide calls. A valid `OTHER` prediction remains a
measured negative. Explicit `--signalp6-bin` paths are forwarded through all
workers. Successful sequences in a partial run are retained; failed batches or
workers contribute to the final nonzero exit status and unscored-allele QC.
Fresh and cached SignalP tables share a header-based parser for five-column
eukaryotic and nine-column output layouts. Positive calls require a valid cleavage
position; `OTHER` calls may leave it empty. SignalP cache hits reparse the saved
table rather than trusting previously parsed JSON.
Cached NetNGlyc results lacking SignalP summaries or positive-call cleavage sites
are not reused when SignalP is enabled; legacy exports need regeneration to
reflect these checks.

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
