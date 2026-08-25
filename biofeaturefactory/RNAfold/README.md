# RNAfold Variant-Centered Pipeline

Quantifies the effect of a variant on local RNA secondary structure using the ViennaRNA Python API. A window (default 151 nt) centered on the variant is folded for both alleles and compared at MFE, ensemble, and per-base accessibility level.

## Requirements

| Component | Notes |
|-----------|-------|
| ViennaRNA | Python bindings (`import RNA`). `conda install -c bioconda viennarna`. |
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. |
| Input FASTAs | Transcript FASTAs. |
| Mapping CSVs | Transcript mappings; intronic variants also need the pre-mRNA mappings. |

## Usage

```bash
# Directory mode: variant_mapping output root
python run_viennaRNA_pipeline.py -i out/ -o results/

# File mode
python run_viennaRNA_pipeline.py \
    -i BRCA1.fasta -o results/ \
    -tm mappings/transcript/ -ipm mappings/intron_premRNA/ \
    -l validation.log -s 1000 -ta 0.05 -w 151 --workers 8
```

Given the `variant_mapping` output root, `-tm` and `-ipm` both default to it, since
`mappings/transcript/` and `mappings/intron_premRNA/` are siblings of `fastas/` under the same
`<GENE>/`. Supply them only for mappings outside that layout. Intronic variants are folded against
both the intron record and the pre-mRNA record.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-i, --input` | Yes | -- | variant_mapping output root, or a single FASTA |
| `-o, --output` | Yes | -- | Output base directory |
| `-w, --window` | No | `151` | Window size in nt (odd; truncated near sequence ends) |
| `-tm, --transcript-mapping` | No | `-i` | File mode only; transcript mapping file or directory |
| `-ipm, --intron-premrna-mapping` | No | derived | File mode only; intron/pre-mRNA mapping file or directory |
| `-l, --log` | No | -- | Validation log (file or directory) |
| `-s, --samples` | No | `1000` | Boltzmann samples per sequence |
| `-ta, --tau` | No | `0.05` | Threshold on `|delta_u|` for `change_flag` |
| `--workers` | No | auto | Max parallel worker processes |

## Output

```
{output}/{GENE}/RNAfold/
    {GENE}.rnafold.tsv           -- per-mutation summary
    {GENE}.rnafold.positions.tsv -- per-position delta_u table
```

### `{GENE}.rnafold.tsv`

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | `GENE-mutation` | -- |
| `transcript_pos` | Transcript coordinate of the variant (1-based) | nt |
| `ddg_mfe_kcalmol` | MFE difference, Alt - Ref | kcal/mol |
| `ddg_ensemble_kcalmol` | Ensemble free energy difference, Alt - Ref | kcal/mol |
| `d_meanE_kcalmol` | Difference in mean sampled energy | kcal/mol |
| `ref_sdE_kcalmol`, `alt_sdE_kcalmol` | SD of sampled energies per allele | kcal/mol |
| `jsd_unpaired_bits` | Jensen-Shannon divergence between per-base unpaired probability profiles | bits (0-1) |
| `delta_central` | Change in unpaired probability at the window midpoint | 0-1 |

### `{GENE}.rnafold.positions.tsv`

| Column | Description | Values |
|--------|-------------|--------|
| `pkey` | `GENE-mutation` | -- |
| `transcript_pos` | Transcript coordinate of the variant | nt |
| `pos` | Index within the window (1 = start). The variant sits at `(window+1)/2` | integer |
| `delta_u` | Per-position change in unpaired probability, Alt - Ref | float |
| `change_flag` | 1 if `|delta_u| >= tau` | 0/1 |
| `direction` | Sign of `delta_u` | -1/0/1 |
| `mfe_change_flag` | 1 if base-pairing status differs between Ref and Alt MFE structures | 0/1 |
| `mfe_change_dir` | 0 = unpaired -> paired, 1 = paired -> unpaired | 0/1 |

## Notes

- Folding temperature is 37 C; all energies in kcal/mol.
- Boltzmann sampling is stochastic, so minor run-to-run variation is expected.
- Workers run with `OMP_NUM_THREADS=1` to prevent oversubscription.

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
