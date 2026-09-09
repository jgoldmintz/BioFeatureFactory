# Mutation Effects Pipelines (EVmutation / adabmDCA)

Scores variants as the change in a Potts model Hamiltonian fitted to a multiple sequence alignment. Two backends and two levels:

- **Protein level** -- amino acid MSA. Primary route for missense and other protein-altering variants.
- **Codon level** -- codon-aware MSA. Synonymous variants and stop annotations; explicit codon-only mode also scores missense variants.
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
scoring backends with per-gene parallelism. Ready alignments start scoring without
waiting for missing alignments on the other side. Nextflow generates missing
protein MSAs with the core MSA pipeline and missing codon MSAs with the core codon
MSA pipeline when the shared CPU/RAM budget permits.

With neither `--msa` nor `-cm`, both paths are available and selection is per gene:

| Mutation classes | Selected paths |
|------------------|----------------|
| Missense plus synonymous/stop-codon variants | Protein and codon |
| Synonymous and/or stop-codon variants only | Codon only |
| Missense only | Protein only |

The controller first searches for existing `MSA/` and `CodonMSA/` alignments under
the output and input roots. Only selected, missing alignments need databases.
Other protein-altering variants retain the protein path; unclassifiable variants
conservatively retain both paths and their QC handling. Genes without remaining
mutations after validation filtering do not launch model or alignment work.

When protein MSA generation is required, `--db-root` must contain a nonempty
`uniref90.fasta`. Jackhmmer needs a rewindable, uncompressed database; gzip-only
input is rejected before Nextflow starts, with a `gzip -dk` preparation command.
Allow space for the expanded database and move aside an existing empty destination
first. The controller does not silently unpack or overwrite database files.
Pre-built protein MSAs do not require this database.

Explicit `--msa` without `-cm` selects protein-only mode, with a warning for
synonymous/stop variants that this will not produce biologically accurate results
for their codon-level effects. Explicit `-cm` without `--msa` selects codon-only
mode, including codon-level missense scores; a warning notes the substantially
higher computation cost relative to amino-acid scoring. These scores need not be
numerically identical. Both flags explicitly enable both paths. `--skip-codon`
remains a per-backend override to protein processing, with a warning when it
overrides a selected codon path. Backend selection itself is unchanged.

Stop-gain/loss rows retain the existing annotation-only scoring behavior; selecting
the codon path does not add numerical stop-effect predictions.

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
| `-ms, --msa` | -- | Protein MSA source; used alone explicitly selects protein-only processing |
| `-cm, --codon-msa` | -- | Codon MSA source; used alone explicitly selects codon-only processing |
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
| `-am, --adabmdca-model` | `pseudoDCA` | Pseudolikelihood without MCMC; explicitly select `bmDCA`/`eaDCA`/`edDCA` for Boltzmann learning |
| `-an, --adabmdca-nepochs` | `500` | Maximum pseudoDCA epochs; an explicitly selected Boltzmann model defaults to 50000; explicit epoch counts override either default |
| `-at, --adabmdca-tol` | `0.001` | pseudoDCA convergence threshold on `\|\|grad\|\|/\|\|grad\|\|_0`; 0 disables |
| `-ap, --adabmdca-patience` | `3` | Consecutive passing checks required to stop |
| `-ace, --adabmdca-check-every` | `10` | Epochs between convergence checks |
| `-ata, --adabmdca-target` | `0.95` | Pearson Cij target |
| `-al, --adabmdca-lr` | `0.01` | Learning rate |
| `-anc, --adabmdca-nchains` | `10000` | Boltzmann-only PCD chain count; unused by pseudoDCA |
| `-ans, --adabmdca-nsweeps` | `10` | Boltzmann-only sweeps per step; unused by pseudoDCA |
| `-ad, --adabmdca-device` | `auto` | Eligible GPU, otherwise CPU; `cpu`, `cuda`, or `cuda:N` forces placement |
| `-adt, --adabmdca-dtype` | `float32` | Dtype |
| `-as, --adabmdca-seed` | `0` | Seed |
| `-t, --threads` | automatic | Share usable CPUs across estimated concurrent jobs; an explicit value sets the per-task default |
| `-vl, --validation-log` | -- | Validation log for mutation filtering |
| `-r, --resume` | off | Resume a previous Nextflow run |

### Resource-aware local execution

The controller uses Nextflow's existing local executor. No additional scheduler,
GPU instance, or scheduler service is required. Both protein and codon tasks use
one GPU process queue, with its concurrency bounded by the number of allocated
GPUs. CPU tasks run in a separate queue. Nextflow admits both queues against a
shared host RAM and CPU budget, including declared MSA and EVmutation requests.

Without `--threads`, the controller estimates concurrency from pending work,
host RAM requests and distinct eligible GPUs, then divides the usable CPU budget
across those jobs. A machine exposing 21 CPUs has a default budget of 20; two
fitting concurrent tasks receive 10 threads each. Five fitting tasks receive four
each. Verified completed tasks are excluded; missing alignments count as generation
tasks, not simultaneously runnable downstream scoring. Explicit per-task thread
overrides take precedence and consume part of the shared budget.
This is a deterministic startup estimate, not optimal packing or live resizing:
threads do not increase as tasks finish, and newly generated alignments are planned
with the same default share. Nextflow still enforces CPU/RAM admission at runtime.

The default `pseudoDCA` model uses a 500-epoch cap and can stop earlier through its
convergence checks. Reaching the cap does not prove convergence. Its epoch count
is not interchangeable with Boltzmann-learning epochs or plmc optimizer iterations.

Each GPU task locks one eligible physical GPU UUID, sets `CUDA_VISIBLE_DEVICES`,
and holds the lock until its backend exits. Tasks select a free eligible card,
not a round-robin index. Eligibility uses total VRAM with headroom; a temporarily
busy card remains eligible, and free VRAM is checked under its lock before launch.
Locks coordinate BFF runs using the same lease directory;
they do not reserve cloud instances or prevent unrelated programs from using CUDA.
Mixed-capacity cards are checked individually. Tasks waiting for a particular
large card can occupy a GPU-process slot; this is not an optimal backfilling scheduler.

```bash
python biofeaturefactory/mutation_effects/mutEffects_controller.py \
    --fasta ~/out --adabmdca-only --output ~/results \
    --adabmdca-nepochs 1 --adabmdca-nchains 64 --resource-plan-only
```

Remove `--resource-plan-only` to execute. Planning-only prints hardware and task
estimates without writing outputs or launching inference. Missing alignments must
be generated before their resource requirements can be estimated; those tasks are
planned inside Nextflow after MSA generation. This is a one-epoch smoke test, not
a recommendation for obtaining converged models.

| Resource flag | Default | Purpose |
|---------------|---------|---------|
| `--resource-cpus` | detected minus one | Shared CPU ceiling; task libraries receive matching thread limits |
| `--resource-memory-gib` | detected available RAM × headroom | Shared host RAM ceiling, not a per-task allowance |
| `--resource-headroom` | `0.9` | Fraction of available host RAM and total GPU VRAM usable for planning |
| `--resource-memory-margin` | `1.15` | Multiplier on estimated peak memory |
| `--resource-overrides` | none | JSON keyed by `GENE.protein` / `GENE.codon` containing measured memory requests |
| `--resource-hardware` | detected | Explicit allocation JSON: `cpus`, `memory_gib`, `gpus` (each has `uuid` and `memory_gib`) |
| `--gpu-lease-dir` | `~/.cache/biofeaturefactory/gpu-leases` | Use the same directory for local BFF runs sharing GPUs |
| `--gpu-wait-timeout` | `600` seconds | Bounded wait for an eligible GPU |
| `--msa-memory-gib` | `8` | Per-MSA-generation host memory request; adjust to database/workload size |
| `--evmutation-memory-gib` | automatic | Optional minimum per-EVmutation RAM request; cannot lower the model estimate |

adabmDCA estimates use alignment width, alphabet, dtype, model, chain count and
sequence count; they include GPU tasks' host-memory requirements for loading and scoring.
The defaults are **uncalibrated conservative envelopes**, not measured guarantees.
Per-task overrides accept `gpu_memory_gib`, `gpu_host_memory_gib`,
`cpu_memory_gib`, and `threads`. Obtain overrides from completed isolated runs,
not from memory observed before an OOM. Nextflow requests are admission accounting,
not OS-enforced memory limits; unrelated processes and underestimated requests can
still exhaust the host. Independent Nextflow runs do not share a host-RAM ledger.

EVmutation requests are calculated separately for each gene's protein or codon
model, not from its mutation count. The plmc estimate follows its gap-reduced
alphabet (`q=20` protein, `q=64` codon) and focus-selected length:
`P = L*q + L*(L-1)*q*q/2`. It reserves nineteen
parameter-sized arrays for the default L-BFGS optimizer and persistent marginals,
assuming double-precision native arithmetic, plus thread-local workspace, alignment
storage and runtime allowances. Thread-local workspace is conservatively reserved
against the usable CPU ceiling, so automatic thread sharing cannot invalidate the
RAM estimate. Admission uses the larger of training memory and dense EVmutation
scoring memory, including the independent-model copy. Supplied native plmc v2
parameters use their header dimensions and require only the scoring estimate.

Prebuilt alignments are planned before launching Nextflow. Missing alignments are
planned immediately after generation, without holding up other ready models.
Nextflow can overlap models only while their combined RAM and CPU requests fit;
splitting a mutation list does not reduce the memory needed to train its model.
A request exceeding the entire RAM budget, or dimensions overflowing plmc's native
integer indexing, is skipped before training instead of starting an unsafe job. These
estimates do not change plmc precision, optimizer settings or scoring results.

Per-gene resource-planning errors skip only the affected backend and side; other
schedulable jobs continue. This also applies to planning after MSA generation.
The controller collects these errors and routing warnings in a final summary,
then exits nonzero if any requested task was blocked. Skipped tasks are not marked
complete, so a later run can retry them with suitable resources while reusing
verified successful results. `--resource-plan-only` includes `resource_errors` and
`warnings` in its JSON report without writing files or launching jobs. Invalid
global configuration, missing required tools/databases, and execution failures
retain their existing fatal behavior; collected diagnostics are still summarized
when the controller can finish normally. Direct Nextflow invocation retains
fail-fast resource planning.

In `auto` mode, tasks too large for any eligible GPU are routed to CPU with the
same adabmDCA model, dtype and training settings. A recognized CUDA OOM can produce
one fresh CPU task with its own RAM request, if that request fits. Explicit CUDA
placement does not silently fall back. Tasks exceeding the entire usable host
budget are skipped and reported rather than waiting indefinitely. CPU execution can be much
slower and is not guaranteed to produce bit-identical floating-point results.

Plans are published in `resource_plans/`. Successful task artifacts have a
`GENE.side.complete.json` manifest with input/model fingerprints, output hashes,
device and observed resource telemetry. Attempts use isolated directories so an
OOM checkpoint cannot be mistaken for a complete model. Existing TSVs without a
matching completion manifest are recomputed; explicitly supplied prebuilt model
parameters remain supported. Nextflow's `--resume` additionally reuses valid work
directories. Resource configurations and input manifests passed to Nextflow are
immutable snapshots under `.bff-resources/`.

EVmutation scores use `GENE.side.routing.json` records to check their inputs,
effective scoring mode, and TSV/model/MSA hashes before reuse. Older EVmutation
tables without these records are rescored using existing parameters when present;
switching between protein-only, codon-only and automatic routing cannot reuse a
table with the wrong mutation distribution.

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

Both standalone backends accept `--score-missense-codon` for exclusive codon-level
scoring, including missense variants. The controller forwards this for forced
codon-only processing when protein-altering or unclassified variants are present.
It cannot be combined with `--skip-codon`. Existing default standalone routing is
unchanged. Codon-mode missense rows carry `MISSENSE_CODON_LEVEL` in `qc_flags`.

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
