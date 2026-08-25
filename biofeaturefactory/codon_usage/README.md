# Codon Usage Pipeline

Per-mutation codon and codon-pair usage metrics (RSCU, W, CAI, tAI, RSCPU, CPS) for the mutated codon and its 5'/3' bicodons.

## Requirements

| Component | Notes |
|-----------|-------|
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. No external tool binary required. |
| ORF FASTA | Nucleotide ORF per gene. |
| Mutations CSV | Nucleotide mutation tokens. Optional; omit to score the ORF only. |

## Usage

```bash
# Directory mode (one subdirectory per gene)
python codon_usage_pipeline.py -f /path/to/fastas/ -m /path/to/mutations/ -o results/

# Single file
python codon_usage_pipeline.py -f BRCA1.fasta -m BRCA1_mutations.csv -o results/
```

## Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `-f, --fasta` | Yes | FASTA file or directory of FASTA files |
| `-m, --mutations` | No | Mutations CSV file or directory of CSV files |
| `-vl, --validation-log` | No | Validation log for filtering failed mutations |
| `-o, --output` | Yes | Output base directory |

## Output

```
{output}/{GENE}/CodonUsage/{GENE}.codon_usage.tsv
```

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `Gene` | Gene symbol |
| `codon_number` | Codon position in ORF (1-based) |
| `codon` | Mutated codon |
| `position_in_codon` | Position of the change within the codon (1-3) |
| `RSCU` | Relative synonymous codon usage |
| `W` | Relative adaptiveness (gene-specific) |
| `CAI_W` | CAI reference W |
| `tAI` | tRNA adaptation weight for the codon |
| `CAI_gene`, `tAI_gene` | Gene-level geometric means |
| `bicodon_3prime` | This codon + next |
| `RSCPU_3prime`, `CPS_3prime`, `noln_CPS_3prime`, `W_CP_3prime` | 3' bicodon metrics |
| `bicodon_5prime` | Previous codon + this codon |
| `RSCPU_5prime`, `CPS_5prime`, `noln_CPS_5prime`, `W_CP_5prime` | 5' bicodon metrics |
| `bicodon_context` | `first_codon_3prime_only`, `last_codon_5prime_only`, `middle_codon_both_directions`, `insufficient_sequence` |
| `qc_flags` | See below |

Bicodon metrics are `null` when the pair is absent from the reference distribution; the row is flagged rather than filled with a placeholder.

## QC flags

| Flag | Meaning |
|------|---------|
| `PASS` | All metrics computed |
| `NO_BICODON` | Could not extract a codon pair |
| `NO_WT_CODON` | WT codon could not be resolved |
| `BICODON_NOT_IN_REFERENCE:<which>` | Named bicodon(s) absent from the reference distribution; those columns are null |
| `CODON_INSERTED` / `CODON_DELETED` | In-frame indel changed the codon count |
| `INDEL_NOT_CODON_ALIGNED:<n>` | Indel length not a multiple of 3 |
| `REF_SPANS_PAST_ORF` | Reference allele extends beyond the ORF |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
