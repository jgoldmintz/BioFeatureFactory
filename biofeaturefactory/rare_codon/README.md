# Rare Codon Enrichment Pipeline

Sliding-window test for rare-codon enrichment or depletion at each mutation position, computed across a codon-aware multiple sequence alignment. Also compares the WT and mutant focus sequence's rare-codon count and fraction in that same window. Wraps the `cg_cotrans` library without requiring a mutant MSA.

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

When `-wg` is omitted, the focus record is selected as `ORF`, then the gene name,
then the first MSA record. An explicit `-wg` must match a record.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-a, --msa` | Yes | -- | Codon-aware MSA FASTA file or directory |
| `-o, --output` | Yes | -- | Output base directory |
| `-m, --mutations` | No | -- | Mutations CSV file or directory |
| `-u, --usage` | No | auto-built | Codon usage `.p.gz` file |
| `-wg, --wt-gi` | No | auto-selected | Focus/WT sequence identifier in the MSA |
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
| `pkey` | Original mutation identity, `GENE-<12-character SHA1>` | string |
| `Gene` | Gene symbol | string |
| `codon_position` | Codon position in the ORF (1-based) | integer |
| `p_enriched` | WT-MSA ensemble p-value for rare-codon enrichment | 0-1 |
| `p_depleted` | WT-MSA ensemble p-value for rare-codon depletion | 0-1 |
| `f_enriched_wt` | Fraction of the WT window that is rare codons | 0-1 |
| `frac_seq_enriched` | Fraction of MSA sequences enriched at this position | 0-1 |
| `frac_seq_depleted` | Fraction of MSA sequences depleted at this position | 0-1 |
| `n_rare` | Rare codons in the WT window | integer |
| `window_size` | Window width in codons | integer |
| `qc_flags` | WT annotation QC; see below | string |
| `n_rare_mut` | Rare codons in the mutant window | integer |
| `f_enriched_mut` | Fraction of the mutant window that is rare codons | 0-1 |
| `delta_n_rare` | Mutant count minus WT `n_rare` | integer |
| `delta_f_enriched` | Mutant rare-codon fraction minus WT fraction | -1 to 1 |
| `comparison_status` | `PASS` for a scored comparison; otherwise the reason it was not scored | string |

## WT–mutant comparison

The CLI computes comparisons automatically, in both file and directory mode.
The MSA and its WT enrichment/null model are computed once per gene. Each variant
is applied independently to the **raw, ungapped focus record**, using the shared
variant parser and REF-checked sequence replacement. The fixed WT rarity rule is
reused; no homolog sequence is changed, no mutant MSA is built, and codon usage is
not re-estimated from mutants.

The four comparison columns are populated for sense-preserving SNVs and MNVs
wholly contained in the existing centered window. This includes synonymous and
missense substitutions. The window starts at `codon_position - window_size // 2`
(1-based) and contains exactly `window_size` codons. One row still describes one
centered window, not every sliding window affected by a variant.

WT and mutant counts use the same counting convention as the WT analysis,
including its exclusion of amino-acid families with zero rare-codon probability
under the selected null model. A zero delta is a measured unchanged rare-codon
count/fraction, **not** evidence of no biological effect. The existing ensemble
p-values remain WT context; neither mutant p-values nor mutation-effect
significance are inferred from them.

Indels, frameshifts and stop-codon changes retain their WT annotations, but their
comparison values are **blank, not zero**. The same applies to invalid REF spans,
missing full windows, or a sanitized/partially gapped focus that no longer maps
that window to the raw ORF frame. A normal terminal stop does not prevent scoring
an upstream sense-preserving substitution. Before accepting a comparison, the
raw WT window count and fraction must agree with the legacy WT results.

Check **both** `qc_flags` and `comparison_status`: WT annotations may pass while
the allele comparison is unsupported. Programmatic callers that omit
`comparison_context` from `process_mutations` retain annotation-only behavior
(`NOT_REQUESTED`). The CLI always supplies the context from the same analysis.

## QC flags

| Flag | Meaning |
|------|---------|
| `POSITION_NOT_IN_WINDOW` | No full centered WT window is available; no truncated window is scored |
| `FRAMESHIFT:downstream_codons_also_change` | Frameshift; every downstream codon changes, so a single-window test does not describe the effect |
| `INVALID_MUTATION` | Token could not be parsed |

Comparison refusal statuses are separate from these WT flags:

| `comparison_status` | Meaning |
|---------------------|---------|
| `UNSUPPORTED_INDEL`, `UNSUPPORTED_FRAMESHIFT` | Length-changing comparison is not implemented |
| `UNSUPPORTED_STOP_CHANGE` | A touched codon is a stop in either allele |
| `UNSUPPORTED_CODON` | A touched codon is ambiguous or incomplete |
| `REF_MISMATCH`, `VARIANT_OUT_OF_RANGE`, `INVALID_MUTATION` | Variant cannot be applied to the raw focus |
| `POSITION_NOT_IN_WINDOW` | No full centered window is available |
| `UNSUPPORTED_FOCUS_FRAME` | Retained MSA codons do not map to the raw window's nucleotide coordinates |
| `VARIANT_SPANS_WINDOW` | The substitution extends beyond the reported window |
| `WT_WINDOW_MISMATCH` | Raw WT count/fraction disagree with the WT analysis |
| `CONTEXT_UNAVAILABLE`, `GENE_MISMATCH`, `FOCUS_SEQUENCE_MISMATCH`, `WINDOW_SIZE_MISMATCH` | Missing or incompatible programmatic comparison context |
| `NOT_REQUESTED` | Programmatic annotation-only call; not emitted by the normal CLI |

## Sequence filtering

Before the test, MSA sequences are filtered against the focus sequence by relative length
(`--max-len-diff`), amino acid identity (`--min-aa-iden`), and gap content.

## Copyright

`cg_cotrans` is Copyright (C) 2017 William M. Jacobs, GPL v3. This pipeline wraps it without
modification. Cite: Jacobs WM, Shakhnovich EI, *PNAS* 114:11434-11439 (2017).

## License

This wrapper is AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
The underlying `cg_cotrans` library remains GPL v3.
