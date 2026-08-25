# NetPhos Pipeline

Kinase phosphorylation-site prediction for WT and mutant protein sequences using NetPhos 3.1 / APE.

## Requirements

| Component | Notes |
|-----------|-------|
| NetPhos 3.1 + APE | [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetPhos-3.1/), academic license. Requires `tcsh`. Point at it with `-nap/--native-ape-path`, or set `NETPHOS_APE_PATH`, `NETPHOS_HOME`, or `NETNGLYC_HOME` (searched in that order). |
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. |
| Input FASTAs | WT ORF FASTAs (`.fasta`, `.fa`, `.fas`, `.fna`). |
| Mapping CSVs | Per-gene mutation CSVs matched as `*{GENE}*.csv`, case-insensitive. |

## Usage

```bash
# Directory mode: variant_mapping output root
python netphos_pipeline.py -i out/ -o results/

# Flat directory or single FASTA + explicit mapping directory
python netphos_pipeline.py -i wt/aaseq/ -o results/ -md mutations/aa/ -l validation.log
```

In directory mode `-i` is the `variant_mapping` output root (`<root>/<GENE>/fastas/` and
`<root>/<GENE>/mappings/`); the gene is taken from the directory name. A flat directory of
FASTAs or a single FASTA also works, in which case supply `-md`.
`input` and `output` are also accepted positionally.

## Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | -- | variant_mapping output root, flat FASTA directory, or single FASTA |
| `-o, --output` | -- | Output base directory |
| `-md, --mapping-dir` | -- | Mutation mapping directory or single CSV |
| `-l, --log` | -- | Validation log file or directory; skips failed mutations |
| `-t, --threshold` | `0.5` | Minimum score for a site to count |
| `-yo, --yes-only` | off | Redefine a site as a NetPhos `YES` answer instead of `score >= --threshold`. It does NOT zero the threshold |
| `-wh, --wt-header` | `ORF` | FASTA header identifying the WT sequence |
| `-bs, --batch-size` | auto | Sequences per APE run. Given explicitly, it bypasses the tiering below |
| `-ti, --timeout` | `300` | Command timeout in seconds |
| `-nc, --no-cache` | off | Disable result caching (`~/.netphos_cache`) |
| `-cc, --clear-cache` | off | Clear the cache and exit |
| `-nap, --native-ape-path` | -- | Explicit path to the APE binary |
| `-v, --verbose` | off | Verbose output |

### Automatic batching

With `-bs` omitted, batch size is chosen from the sequence count:

| Sequences | Behaviour |
|-----------|-----------|
| 1 | Single APE run |
| 2-10 | Single run first; on failure, retry at batch size **1** |
| 11-100 | Batch size 25 |
| > 100 | Batch size 50 |

The 2-10 fallback is 1, not 10, deliberately: on the native APE install a single
sequence predicts, while 2, 3 and 6 all segfault every kinase and emit nothing.
One sequence per invocation is the only configuration observed to work. APE itself
caps at roughly 2000 sequences per run.

Each successful batch is content-cached, so a rerun re-runs only the failed
batches. Mutants in a failed batch are written to `<output>.dropped` and printed --
a partial run is never reported as a success.

## Output

```
{output}/{GENE}/NetPhos/
    {GENE}.tsv          -- per-mutation summary
    {GENE}.events.tsv   -- per-site, per-kinase events
    {GENE}.sites.tsv    -- raw WT and MUT predictions
```

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `Gene` | Gene symbol |
| `n_sites_wt`, `n_sites_mut` | Sites at or above `--threshold` |
| `count_gained`, `count_lost`, `count_strengthened`, `count_weakened`, `count_stable` | Classification tallies |
| `max_abs_delta`, `sum_abs_delta` | Max and sum of absolute score changes |
| `n_kinases_affected` | Kinases with any event |
| `top_event_type` | Dominant event label |
| `top_event_delta` | Score change for the dominant event |
| `top_event_position` | Residue index of the dominant event |
| `top_event_kinase` | Kinase for the dominant event |
| `top_event_classification_code` | Numeric encoding of the dominant event |
| `qc_flags` | Quality control flags |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene` | Identifiers |
| `position` | 1-based residue position |
| `amino_acid_wt`, `amino_acid_mut` | WT and MUT residues |
| `kinase` | Kinase motif classifier |
| `wt_score`, `mut_score`, `delta` | NetPhos scores and MUT - WT |
| `wt_answer`, `mut_answer` | NetPhos YES/NO flag |
| `classification` | gained / lost / strengthened / weakened / stable / subthreshold |
| `classification_code` | gained=2, lost=-2, strengthened=1, weakened=-1, stable=0, subthreshold=-3 |

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene` | Identifiers |
| `allele` | Sequence label (wt or mutant ID) |
| `seq_name` | Sequence name from NetPhos output |
| `position` | 1-based residue position |
| `amino_acid` | Residue (S/T/Y) |
| `context` | Flanking peptide window |
| `score` | NetPhos score |
| `kinase` | Kinase motif classifier (CKII, cdc2, unsp, ...) |
| `answer` | NetPhos YES/NO flag |

Classification: a site is **gained** if it crosses `--threshold` only in MUT, **lost** if only in WT,
and **strengthened**/**weakened** when both are above threshold and `|delta| > 0.05`.

## Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `tcsh: Command not found` | Install tcsh; APE requires it |
| `No mapping file found` | `-md` must contain files named `*{GENE}*.csv` |
| APE binary not found | Make it executable; pass `-nap` or set `NETPHOS_APE_PATH` |
| Stale results | `python3 netphos_pipeline.py -cc` |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
