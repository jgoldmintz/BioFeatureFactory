# NetMHC Pipeline

MHC class I binding prediction for WT and mutant protein sequences, comparing WT vs mutant to call gained/lost epitopes. Only netMHC-4.0 is supported.

## Requirements

| Component | Notes |
|-----------|-------|
| netMHC 4.0 | [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHC-4.0/), academic license. Point at it with `-nnp/--native-netmhc-path` (binary or its parent directory), or set `NETMHC_PATH` / `NETMHC_HOME`. |
| netMHC 4.0 `data.tar.gz` | **Separate download**, ~63 MB compressed / ~210 MB extracted. The distribution tarball ships `<platform>/data` as a symlink to `../data`, which does not exist until this is extracted into the install root. |
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. |
| Input FASTAs | ORF FASTAs (header `>ORF`); translated to amino acids automatically. |
| Mutation CSVs | Single-column NT mutations (header `mutant`) or multi-column mapping CSVs. |

### Verifying the install

Without `data/`, netMHC still exits 0 and prints its banner but emits no prediction rows. Check
against the shipped test file before running the pipeline:

```bash
./netMHC -a HLA-A0201 -f test/test.fsa | awk 'NF>=14 && $1 ~ /^[0-9]+$/' | wc -l
```

`test.fsa` is 245 residues, so a working 9-mer install prints **237** rows (245 - 9 + 1). Zero rows
means `data/` is missing. Also confirm the banner's `-tdir` scratch directory exists and is writable.

netMHC 4.0 uses `HLA-A0201` allele syntax (not the NetMHCpan `HLA-A*02:01` form). `-listMHC` prints
the accepted names.

## Usage

```bash
# Directory mode: variant_mapping output root
python netmhc_pipeline.py -i out/ -o results/ -a HLA-A0201 HLA-A0101 HLA-B0702

# Flat FASTA directory + explicit mutation directory
python netmhc_pipeline.py -i FASTA_files/nt/ -o results/ -m mutations/ -l validation.log
```

In directory mode `-i` is the `variant_mapping` output root (`<root>/<GENE>/fastas/` and
`<root>/<GENE>/mappings/`); the gene is taken from the directory name. `input` and `output`
are also accepted positionally.

## Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | -- | variant_mapping output root, flat FASTA directory, or single FASTA |
| `-o, --output` | -- | Output base directory |
| `-m, --mutations` | -- | Mutation file or directory of mutation CSVs |
| `-l, --log` | -- | Validation log; skips failed mutations |
| `-a, --alleles` | `HLA-A0201` | HLA alleles, space-separated |
| `-t, --threshold` | `0.5` | Rank threshold for strong binders |
| `-nt, --netmhc-tool` | `netMHC` | Only `netMHC` (4.0) is supported |
| `-nnp, --native-netmhc-path` | -- | Path to the netMHC executable |
| `-ti, --timeout` | `600` | Command timeout in seconds |
| `--max-workers` | `4` | Concurrent netMHC runs per gene; `1` for serial |
| `-bs, --batch-size` | `100` | Deprecated and unused; netMHC input is no longer batched |
| `-v, --verbose` | off | Verbose output |

## Output

```
{output}/{GENE}/NetMHC/
    {GENE}.tsv          -- per-mutation summary
    {GENE}.events.tsv   -- per-peptide, per-allele events
    {GENE}.sites.tsv    -- raw binding predictions
```

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene`, `mutation` | Identifiers |
| `n_epitopes_wt`, `n_epitopes_mut` | Peptides with rank below `--threshold` |
| `count_gained`, `count_lost`, `count_strengthened`, `count_weakened` | Classification tallies |
| `max_abs_delta_rank`, `sum_abs_delta_rank` | Max and sum of rank changes (percentile) |
| `top_event_type` | Dominant event label |
| `top_event_allele`, `top_event_peptide` | Allele and peptide for the dominant event |
| `top_event_delta_rank` | Rank change for the dominant event |
| `qc_flags` | `missing_wt`, `missing_mut`, `no_delta`, and `WT_*` / `MUT_*` tool-failure codes |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene`, `mutation` | Identifiers |
| `peptide` | Peptide sequence |
| `pos` | Position in the protein |
| `mhc_allele` | HLA allele |
| `wt_rank`, `mut_rank`, `delta_rank` | Percentile ranks and MUT - WT |
| `wt_affinity`, `mut_affinity`, `delta_affinity` | Binding affinity (nM) and MUT - WT |
| `bind_level_wt`, `bind_level_mut` | SB / WB / NB |
| `classification` | gained / lost / strengthened / weakened / stable |
| `classification_code` | gained=2, lost=-2, ... |

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene` | Identifiers |
| `sequence_type` | wt or mut |
| `pos` | Position in the protein |
| `mhc_allele` | HLA allele |
| `peptide` | Peptide sequence |
| `core` | Core binding region |
| `affinity` | Binding affinity (nM) |
| `rank` | Percentile rank |
| `bind_level` | SB (rank <= 0.5), WB (rank <= 2), otherwise blank |
| `identity` | Sequence name from netMHC output |

netMHC scores every window, so a length-L sequence yields exactly `L - l + 1` rows per allele
including non-binders. Zero rows is a tool failure, not a negative result; the pipeline echoes
netMHC's raw output and sets a `qc_flags` code when that happens.

## Troubleshooting

| Symptom | Resolution |
|---------|------------|
| Runs "succeed" with 0 predictions | `data/` is missing; extract `data.tar.gz` into the install root and re-run the 237-row check above |
| `Unable to open(r) .../data/version` | Same cause: the `data` symlink is dangling |
| `NetMHC not found` | Pass `-nnp`, or set `NETMHC_PATH` / `NETMHC_HOME`; check the executable bit |
| `No mapping file found` | `-m` must point at a CSV or a directory of `{GENE}*.csv` |
| `Invalid allele name` | Use netMHC 4.0 syntax (`HLA-A0201`); check `-listMHC` |
| Peptide too short | Sequences shorter than the peptide length produce no windows |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
