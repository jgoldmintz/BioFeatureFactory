# Rare Codon Enrichment Pipeline

Sliding-window test for rare-codon enrichment or depletion at each mutation position, computed across a codon-aware multiple sequence alignment. Wraps the `cg_cotrans` library.

## Requirements

| Component | Notes |
|-----------|-------|
| `cg_cotrans` | Download from [Shakhnovich Lab](https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis). The pipeline imports directly from the `cg_cotrans/` subdirectory. |
| Python >= 3.8 | With `numpy`, `scipy` |
| Codon-aware MSA | Triplet FASTA, gaps as `---`. From `core/codon_msa_pipeline.py`. |
| Mutations CSV | Nucleotide mutation tokens |
| Reference codon usage TSV | Strongly recommended; see `-rcu` below |

Generate the codon MSA first:

```bash
python ../core/codon_msa_pipeline.py -i out/ -o codon_msas/
```

## Usage

```bash
# Directory mode
python rare_codon_pipeline.py -a codon_msas/ -m mutations/ \
    -rcu /path/to/human_GRCh38_codon_usage.tsv -o results/

# Single gene
python rare_codon_pipeline.py \
    -a BRCA1.codon.msa.fasta -wg "Human_BRCA1" -m BRCA1_mutations.csv \
    -u codon_usage.p.gz -rcu human_GRCh38_codon_usage.tsv \
    -L 21 -rt 0.15 --min-aa-iden 0.6 -o results/
```

`-wg` defaults to the gene name in directory mode.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-a, --msa` | Yes | -- | Codon-aware MSA FASTA file or directory |
| `-o, --output` | Yes | -- | Output base directory |
| `-m, --mutations` | No | -- | Mutations CSV file or directory |
| `-u, --usage` | No | auto-built | Codon usage `.p.gz` file |
| `-wg, --wt-gi` | No | gene name | Focus/WT sequence identifier in the MSA |
| `-rcu, --reference-codon-usage` | No | -- | Genome-wide codon usage TSV defining which codons are rare. Columns: `codon`, `amino_acid`, `count`, `relative_usage_within_aa` |
| `-L, --window-size` | No | `15` | Sliding window width in codons |
| `-rm, --rare-model` | No | `no_norm` | Rare codon definition (`no_norm`, `cmax_norm`) |
| `-rt, --rare-threshold` | No | `0.1` | Frequency threshold below which a codon is rare |
| `-nm, --null-model` | No | `genome` | Null model (`genome`, `eq`, `groups`) |
| `--max-len-diff` | No | `0.2` | Max relative length difference vs the focus sequence |
| `--min-aa-iden` | No | `0.5` | Min amino acid identity vs the focus sequence |
| `-vl, --validation-log` | No | -- | Validation log for filtering |

**`-rcu` matters.** Without it, "rare" is derived from the single gene under analysis, which makes
the null model self-referential -- the reference distribution and the test data are the same
sequence. See the README beside the table in `<Bio_DBs>/cocoputs/`.

## Output

```
{output}/{GENE}/RareCodon/{GENE}.rare_codon.tsv
```

| Column | Description | Range |
|--------|-------------|-------|
| `pkey` | `GENE-mutation` | string |
| `Gene` | Gene symbol | string |
| `codon_position` | Codon position in the ORF (1-based) | integer |
| `p_enriched` | P-value for rare-codon enrichment in the window | 0-1 |
| `p_depleted` | P-value for rare-codon depletion in the window | 0-1 |
| `f_enriched_wt` | Fraction of the WT window that is rare codons | 0-1 |
| `frac_seq_enriched` | Fraction of MSA sequences enriched at this position | 0-1 |
| `frac_seq_depleted` | Fraction of MSA sequences depleted at this position | 0-1 |
| `n_rare` | Rare codons in the window | integer |
| `window_size` | Window width in codons | integer |
| `qc_flags` | See below | string |

## QC flags

| Flag | Meaning |
|------|---------|
| `POSITION_NOT_IN_WINDOW` | Position near a sequence edge; window is truncated |
| `FRAMESHIFT:downstream_codons_also_change` | Frameshift; every downstream codon changes, so a single-window test does not describe the effect |
| `INVALID_MUTATION` | Token could not be parsed |

## Sequence filtering

Before the test, MSA sequences are filtered against the focus sequence by relative length
(`--max-len-diff`), amino acid identity (`--min-aa-iden`), and gap content.

## Copyright

`cg_cotrans` is Copyright (C) 2017 William M. Jacobs, GPL v3. This pipeline wraps it without
modification. Cite: Jacobs WM, Shakhnovich EI, *PNAS* 114:11434-11439 (2017).

## License

This wrapper is AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
The underlying `cg_cotrans` library remains GPL v3.
