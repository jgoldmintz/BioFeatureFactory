# Mutation Effects Pipelines (EVmutation / adabmDCA)

Scores variants as the change in a Potts model Hamiltonian fitted to a multiple sequence alignment. Two backends and two levels:

- **Protein level** -- amino acid MSA. Missense and stop variants.
- **Codon level** -- codon-aware MSA. Synonymous variants, which have no protein-level score.
- **EVmutation backend** -- model fitted by `plmc`.
- **adabmDCA backend** -- model fitted by `adabmDCA train`.

Sign convention for every score: positive = tolerated/favoured, negative = deleterious.

## Requirements

| Component | Notes |
|-----------|-------|
| EVmutation library | Cloned from [Marks Lab](https://github.com/debbiemarkslab/EVmutation) |
| `plmc` | Compiled from https://github.com/debbiemarkslab/plmc. Needed when the EVmutation backend runs |
| adabmDCA | Cloned; GPU strongly recommended. Needed when the adabmDCA backend runs |
| Python >= 3.8 | With `numpy`, `pandas`, `numba` |
| Protein MSA | A2M or FASTA; jackhmmer or HHblits against UniRef90 |
| Codon MSA | Triplet FASTA from `core/codon_msa_pipeline.py`; gaps as `---` |
| ORF FASTA | Nucleotide CDS; record ID `ORF`, else the first record |
| Nextflow | Only for the controller. `curl -s https://get.nextflow.io \| bash` |
| jackhmmer / mmseqs2 / mafft | Only when MSAs must be generated |

## Controller (multi-gene, Nextflow)

`mutEffects_controller.py` validates dependencies, then orchestrates MSA generation and both
scoring backends with per-gene parallelism. Protein MSA (jackhmmer) and codon MSA (mmseqs2 + MAFFT)
run in parallel per gene; scoring waits for both.

```bash
# Full run, both backends, generating MSAs
python mutEffects_controller.py \
    -f fastas/ -m mutations/ \
    -dr /path/to/Bio_DBs/ \
    -pb /path/to/plmc \
    -o results/

# Pre-built MSAs, EVmutation only
python mutEffects_controller.py \
    -f fastas/ -m mutations/ \
    -ms prebuilt_protein_msas/ -cm prebuilt_codon_msas/ \
    -pb /path/to/plmc -eo -o results/

# adabmDCA only, pseudolikelihood (lower peak GPU memory)
python mutEffects_controller.py \
    -f fastas/ -m mutations/ \
    -ms prebuilt_protein_msas/ -cm prebuilt_codon_msas/ \
    -ao -am pseudoDCA -o results/
```

### Controller arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-f, --fasta` | required | ORF FASTA file or directory of per-gene FASTAs |
| `-m, --mutations` | -- | Mutations CSV file or directory of per-gene CSVs |
| `-o, --output` | -- | Output base directory |
| `-pb, --plmc-binary` | -- | Path to `plmc`. Required whenever the EVmutation backend runs |
| `-dr, --db-root` | -- | Bio_DBs root (`uniref90.fasta`, `refseq_assemblies/`, ...). Required only when a gene needs MSA generation |
| `-ms, --msa` | -- | Protein MSA file or directory |
| `-cm, --codon-msa` | -- | Codon MSA file or directory |
| `-mp, --model-params` | -- | Pre-built protein plmc params file or directory |
| `-cmp, --codon-model-params` | -- | Pre-built codon plmc params file or directory |
| `-app, --adabmdca-protein-params` | `<output>/adabmdca_protein_params/` | Pre-built adabmDCA protein params |
| `-acp, --adabmdca-codon-params` | `<output>/adabmdca_codon_params/` | Pre-built adabmDCA codon params |
| `-eo, --evmutation-only` | off | Run only the EVmutation/plmc backend |
| `-ao, --adabmdca-only` | off | Run only the adabmDCA backend |
| `-sc, --skip-codon` | -- | Skip codon scoring. No argument = both backends; `evmutation` or `adabmdca` targets one. Synonymous and stop variants route to the protein TSV |
| `-jb, --jackhmmer-binary` | `jackhmmer` | jackhmmer path |
| `-ji, --jackhmmer-iterations` | `5` | jackhmmer iterations |
| `-mb, --mmseqs-binary` | `mmseqs` | mmseqs2 path |
| `-a, --aligner` | `mafft` | Protein aligner |
| `-am, --adabmdca-model` | `bmDCA` | `bmDCA`/`eaDCA`/`edDCA` are Boltzmann-learning (high memory); `pseudoDCA` is pseudolikelihood (~2x less peak GPU memory, no MCMC) |
| `-an, --adabmdca-nepochs` | 500 (pseudoDCA) / 50000 | Max epochs |
| `-at, --adabmdca-tol` | `0.001` | pseudoDCA convergence threshold on `\|\|grad\|\|/\|\|grad\|\|_0`; 0 disables |
| `-ap, --adabmdca-patience` | `3` | Consecutive passing checks required to stop |
| `-ace, --adabmdca-check-every` | `10` | Epochs between convergence checks |
| `-ata, --adabmdca-target` | `0.95` | Pearson Cij target |
| `-al, --adabmdca-lr` | `0.01` | Learning rate |
| `-anc, --adabmdca-nchains` | `10000` | PCD chain count |
| `-ans, --adabmdca-nsweeps` | `10` | Sweeps per step |
| `-ad, --adabmdca-device` | `cuda` | Device |
| `-adt, --adabmdca-dtype` | `float32` | Dtype |
| `-as, --adabmdca-seed` | `0` | Seed |
| `-t, --threads` | `4` | Threads |
| `-vl, --validation-log` | -- | Validation log for mutation filtering |
| `-r, --resume` | off | Resume a previous Nextflow run |

The Boltzmann path emits an OOM hint suggesting `pseudoDCA` when memory becomes the bottleneck.

## Single backend, no Nextflow

### `evmutation_pipeline.py`

```bash
# Pre-built params
python evmutation_pipeline.py -f SMN2.fasta -m SMN2.csv \
    --model-params SMN2.model_params -cmp SMN2.codon_model_params -o results/

# Build params from MSAs
python evmutation_pipeline.py -f SMN2.fasta -m SMN2.csv \
    --msa SMN2.msa.a2m -cm SMN2.codon.msa.fasta \
    -pb /usr/local/bin/plmc -o results/

# Multi-gene directory mode
python evmutation_pipeline.py -f fastas/ -m mutations/ \
    --msa msas/ -cm codon_msas/ -pb /usr/local/bin/plmc -o results/
```

| Flag | Default | Description |
|------|---------|-------------|
| `-f, --fasta` | required | ORF FASTA file or directory |
| `-m, --mutations` | -- | Mutations CSV file or directory |
| `--model-params` | -- | Protein params file, or directory containing `{GENE}.model_params` |
| `-cmp, --codon-model-params` | -- | Codon params file, or directory containing `{GENE}.codon_model_params` |
| `--msa` | -- | Protein MSA file or directory; triggers plmc when params are absent |
| `--focus` | gene name | Focus sequence ID in the protein MSA |
| `-cm, --codon-msa` | -- | Codon MSA file or directory; encoded then run through plmc |
| `-cf, --codon-focus` | `ORF` | Focus sequence ID in the codon MSA |
| `-g, --gene` | -- | Gene name override (single-gene mode) |
| `-pb, --plmc-binary` | -- | plmc path; required when running plmc |
| `-a, --alphabet` | `-ACDEFGHIKLMNPQRSTVWY` | Protein alphabet |
| `-le, --lambda-e` | `16.2` | L2 regularisation on the pairwise terms |
| `-lh, --lambda-h` | `0.01` | L2 regularisation on the site fields |
| `-sp, --skip-plmc` | off | Skip plmc; codon MSA encoding still runs |
| `-sc, --skip-codon` | off | Route synonymous and stop variants to the protein TSV |
| `-vl, --validation-log` | -- | Validation log for mutation filtering |
| `-o, --output` | `.` | Output directory |
| `-q, --quiet` | off | Suppress verbose output |

### `adabmdca_pipeline.py`

Same input contract; `-pp/--protein-params` and `-cp/--codon-params` replace the plmc params flags,
`-st/--skip-train` replaces `-sp/--skip-plmc`, `-pa/--protein-alphabet` replaces `-a/--alphabet`,
and the `-am`/`-an`/`-at`/`-ap`/`-ace`/`-ata`/`-al`/`-anc`/`-ans`/`-ad`/`-adt`/`-as` training
options are as listed in the controller table.

### Parameter resolution

For each params flag: a **file** is used directly; a **directory** is searched for
`{GENE}.<suffix>`; if omitted and the matching MSA flag is given, the params file is written beside
the MSA and the trainer runs if it does not yet exist; if both are omitted that model is skipped.
This lets MSAs and params live in different directories.

## Output

```
{output}/{GENE}/EVmutation/{GENE}.protein.tsv
{output}/{GENE}/EVmutation/{GENE}.codon.tsv
{output}/{GENE}/adabmDCA/{GENE}.protein.tsv
{output}/{GENE}/adabmDCA/{GENE}.codon.tsv
```

### EVmutation `{GENE}.protein.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `nt_mutant` | Source nucleotide token |
| `codon_position` | Codon position in the ORF (1-based) |
| `wt_codon`, `mut_codon` | WT and mutant codons |
| `mutant` | Amino acid substitution string (e.g. `V123A`) |
| `pos`, `wt`, `subs` | AA position (1-based), WT residue, substituted residue |
| `mutation_class` | `MISSENSE`, `SYNONYMOUS`, `STOP_GAIN`, `STOP_LOSS`, `INFRAME_INS`, `INFRAME_DEL`, `INFRAME_DELINS`, `FRAMESHIFT`, `UNKNOWN` |
| `prediction_epistatic` | Full Potts model score (primary at high Neff) |
| `prediction_independent` | Site-field-only score (primary at low Neff) |
| `epistatic_contribution` | `epistatic - independent` |
| `site_entropy` | Shannon entropy at the position (bits) |
| `mean_epistatic_at_pos`, `std_epistatic_at_pos` | Across all substitutions at that position |
| `z_score_epistatic` | Z-score relative to those substitutions |
| `frequency` | Observed frequency of the substitution in the MSA |
| `column_conservation` | Max single-residue frequency at the position |
| `qc_flags` | See below |

### EVmutation `{GENE}.codon.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `nt_mutant`, `codon_position`, `wt_codon`, `mut_codon` | As above |
| `mutation_class` | `SYNONYMOUS`, `STOP_GAIN`, `STOP_LOSS` |
| `prediction_codon_independent` | Site-field-only score; primary score for synonymous variants |
| `prediction_codon_epistatic` | Full Potts model score |
| `codon_epistatic_contribution` | `epistatic - independent` |
| `codon_epistatic_concordance` | `CONCORDANT`, `DISCORDANT`, `NEUTRAL` |
| `codon_frequency` | Observed frequency of the mutant codon at that position |
| `qc_flags` | See below |

Concordance uses a per-position noise floor of `0.5 * std(contributions at that position)`:
above it and matching in sign is `CONCORDANT`, above it and opposite is `DISCORDANT`, otherwise
`NEUTRAL`.

### adabmDCA tables

Same layout with backend-suffixed score columns.

`{GENE}.protein.tsv`: `pkey`, `nt_mutant`, `codon_position`, `wt_codon`, `mut_codon`, `mutant`,
`pos`, `wt`, `subs`, `mutation_class`, `prediction_protein_independent_adabm`,
`prediction_protein_epistatic_adabm`, `protein_pairwise_contribution_adabm`,
`protein_concordance_adabm`, `frequency_adabm`, `qc_flags`.

`{GENE}.codon.tsv`: `pkey`, `nt_mutant`, `codon_position`, `wt_codon`, `mut_codon`,
`mutation_class`, `prediction_codon_independent_adabm`, `prediction_codon_epistatic_adabm`,
`codon_pairwise_contribution_adabm`, `codon_concordance_adabm`, `codon_frequency_adabm`, `qc_flags`.

## QC flags

| Flag | Meaning |
|------|---------|
| `SYNONYMOUS_SCORED` | Synonymous variant scored from the codon model |
| `SYNONYMOUS_UNSCORED` | No codon model supplied; codon annotations only |
| `SYNONYMOUS_NOT_IN_CODON_MODEL` | Codon model loaded but the position or codon is absent |
| `SYNONYMOUS_PROTEIN_LEVEL` | Routed to the protein TSV by `--skip-codon` |
| `NOT_IN_MODEL` | Position not in the model index |
| `NO_MODEL` / `NO_PROTEIN_MODEL` | Required params were not supplied (adabmDCA) |
| `MODEL_WT_MISMATCH` | Model's WT token disagrees with the ORF |
| `NON_ORF_TOKEN:no_residue_or_codon_site_in_potts_model` | Intronic token; a Potts model is indexed by residue or codon site, so there is nothing to evaluate. Score these with RNAfold, miranda, genesplicer, or AlphaFold3 |
| `INSERTION_NOT_REPRESENTABLE_FIXED_L` | Insertion cannot be expressed in a fixed-length model |
| `DELETION_NO_GAP_SYMBOL_IN_MODEL` | Deletion needs a gap token the alphabet lacks |
| `FRAMESHIFT_NOT_REPRESENTABLE_FIXED_L` | Frameshift cannot be expressed in a fixed-length model |
| `NO_CODON_CHANGED` | Token resolves to no codon change |
| `UNKNOWN_CODON` | Codon not in the standard table |
| `INVALID_MUTATION` | Token could not be parsed |
| `OUT_OF_RANGE` | Position outside the ORF |
| `PARTIAL_CODON` | Variant falls in an incomplete codon |
| `REF_MISMATCH` | Reference nucleotide disagrees with the ORF |
| `Z_SCORE_UNDEFINED_MULTISITE` | Z-score undefined for a multi-site token |
| `CONCORDANCE_UNDEFINED_MULTICODON` / `CONCORDANCE_UNDEFINED_MULTISITE` | Concordance undefined for multi-codon / multi-site tokens |

Intronic tokens still get a row, with every metric column empty -- dropping them would leave a hole
for anyone joining pipelines on `pkey`, indistinguishable from a variant that was never submitted.

## MSA depth

The epistatic score needs enough sequences to estimate pairwise statistics; the independent score
does not. Below the threshold, L2 regularisation drives the pairwise terms toward zero and the
epistatic score converges on the independent one.

| Model | Alphabet | Epistatic score reliable at |
|-------|----------|----------------------------|
| Protein | 20 residues + gap | Neff >> 6,500 |
| Codon | 64 codons + gap | Neff >= 640 x L |

The codon threshold is higher because a 64-character alphabet has far more free pairwise
parameters. `prediction_codon_independent` stays meaningful at any depth above ~10 sequences.

Codon coupling reflects co-variation in **codon choice** across the MSA -- largely phylogenetic
GC/AT usage bias -- not structural contact. Do not read `codon_epistatic_contribution` as a physical
interaction between CDS positions.

Generate a codon MSA with:

```bash
python ../core/codon_msa_pipeline.py -i out/ -o codon_msas/
```

## plmc settings

| Flag | Default | Description |
|------|---------|-------------|
| `-le` | 16.2 | L2 on the pairwise terms |
| `-lh` | 0.01 | L2 on the site fields |
| `-m` | 500 | Maximum iterations |
| `-t` | 0.2 | Step size |
| `-g` | -- | Sequence reweighting |

Runtime: roughly 10-60 minutes per gene for protein models, 30-120 for codon models.

## Copyright

`EVmutation` (`model.py`, `tools.py`) is from the EVmutation package by Thomas A. Hopf, Marks Lab,
Harvard -- https://github.com/debbiemarkslab/EVmutation. `plmc`:
https://github.com/debbiemarkslab/plmc. `adabmDCA` is a separate upstream project. These are
vendored clones; configure them from the caller rather than editing them.

## License

These wrappers are AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
