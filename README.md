# BioFeatureFactory

Modular bioinformatics framework for automated feature extraction, coordinate mapping, and predictive modeling using gene- and variant-level inputs.

Give it a list of mutations. It resolves them into every coordinate space a
prediction tool might want, then runs twelve feature pipelines against the
result and returns tables that all join on a single key.

## Table of Contents
1. [How It Works, End to End](#how-it-works-end-to-end)
2. [Core Stages](#core-stages)
3. [Architecture Overview](#architecture-overview)
4. [Installation](#installation)
5. [Pipeline Quick-Start Matrix](#pipeline-quick-start-matrix)
6. [Pipeline Summaries](#pipeline-summaries)
   - [NetNGlyc](#netnglyc-pipeline)
   - [NetPhos](#netphos-pipeline)
   - [NetMHC](#netmhc-pipeline)
   - [NetSurfP3](#netsurfp3-pipeline)
   - [Miranda](#miranda-pipeline)
   - [GeneSplicer](#genesplicer-pipeline)
   - [SpliceAI](#spliceai-pipeline)
   - [RNAfold](#rnafold-pipeline)
   - [AlphaFold3](#alphafold3-pipeline)
   - [Mutation Effects: EVmutation / adabmDCA](#mutation-effects-pipeline)
   - [Codon Usage](#codon-usage-pipeline)
   - [Rare Codon](#rare-codon-pipeline)
7. [Shared Infrastructure & Pathing](#shared-infrastructure--pathing)
8. [Extending BioFeatureFactory](#extending-biofeaturefactory)
9. [Citation](#citation)
10. [License](#license)
11. [Support](#support)

---

## How It Works, End to End

Mutations in, feature tables out. Every stage after the first reads from one
directory.

### 1. Map

`core/variant_mapping.py` takes a mutation list, an indexed reference genome and
an annotation, and writes one directory per gene:

```
<out>/<GENE>/fastas/<GENE>.fasta
<out>/<GENE>/mappings/chromosome/      absolute genomic coordinates
<out>/<GENE>/mappings/gDNA/            relative to the gene's genomic slice
<out>/<GENE>/mappings/transcript/      relative to the spliced transcript
<out>/<GENE>/mappings/aa/              amino-acid change notation
<out>/<GENE>/mappings/intron_premRNA/  intronic and pre-mRNA coordinates
<out>/<GENE>/mappings/pkey/            token <-> primary key inversion table
<out>/<GENE>/mappings/mutations/       <GENE>_mutations.csv
```

```bash
python3 biofeaturefactory/core/variant_mapping.py \
  -m mutations/ -r reference.fa -a annotations.gtf -o mapped/ -v
```

This is the only place in the framework where a coordinate is converted. Every
pipeline downstream reads the CSV in which that conversion is already recorded.

### 2. Align (only for the pipelines that need an MSA)

The two MSA generators add to the same tree:

```bash
python3 biofeaturefactory/core/msa_generation_pipeline.py -f mapped/ -d uniref90.fasta -j jackhmmer -o mapped/
python3 biofeaturefactory/core/codon_msa_pipeline.py      -f mapped/ -d Bio_DBs/ -o mapped/
```

```
<out>/<GENE>/MSA/<GENE>.msa.a2m                protein alignment
<out>/<GENE>/CodonMSA/<GENE>.codon.msa.fasta   codon-aware alignment
```

### 3. Extract

Point any pipeline at the same root with a single flag. Mapping directories,
mutation lists and MSAs are derived from it:

```bash
python3 biofeaturefactory/RNAfold/run_viennaRNA_pipeline.py       -i mapped/ -o results/
python3 biofeaturefactory/codon_usage/codon_usage_pipeline.py     -f mapped/ -o results/
python3 biofeaturefactory/rare_codon/rare_codon_pipeline.py       -a mapped/ -o results/
python3 biofeaturefactory/miranda/miranda_ensemble.py             -i mapped/ -o results/
python3 biofeaturefactory/spliceai/spliceai_pipeline_controller.py -mp mapped/ -o results/
```

Supply mapping paths explicitly only when your files live outside this layout;
those flags are marked **FILE MODE ONLY** in each `--help`. In directory mode
they are optional.

### 4. Join

Every row every pipeline writes carries a `pkey` -- `{GENE}-{sha1(token)[:12]}` --
minted once, during mapping. Twelve tools' outputs concatenate on that one
column with no coordinate reconciliation, which is what makes the results usable
as a feature matrix.

---

## Core Stages

The input-producing CLIs under `biofeaturefactory/core/`.

| Stage | Entry point | Produces |
|---|---|---|
| **Variant mapping** | `core/variant_mapping.py` | Per-gene FASTAs and chromosome / gDNA / transcript / AA / intron-pre-mRNA / pkey / mutations CSVs |
| **VCF conversion** | `core/vcf_converter.py` | Sorted per-gene VCFs; cross-validates FASTA headers against the reference and normalises chromosome naming |
| **Protein MSA** | `core/msa_generation_pipeline.py` | `<GENE>/MSA/<GENE>.msa.a2m` via jackhmmer, with N_eff and gap filtering |
| **Codon MSA** | `core/codon_msa_pipeline.py` | `<GENE>/CodonMSA/<GENE>.codon.msa.fasta` via mmseqs2 + back-translation |

### Core as the extension point

`core/` is deliberately the only stage that knows about genome coordinates. A
mutation is converted between ORF, transcript, gDNA and chromosome space exactly
once, inside `core/variant_mapping.py`; no feature pipeline performs that
conversion, and none of them import `core` at all. What crosses the boundary is
files -- plain CSV and FASTA in a fixed per-gene layout -- so the interface
between the mapping stage and everything downstream is inspectable, diffable and
language-agnostic.

That is what makes the framework expandable rather than merely modular. Adding a
new analysis means writing one consumer against artifacts that already exist: it
requires no change to `core/`, no change to `lib/`, and no change to any of the
twelve pipelines already present. The most recent addition is the worked example
-- the adabmDCA mutation-effect backend imports eighteen names from
`biofeaturefactory.lib.utility` and nothing else from the repository, reading the
codon alignment and mutation list that `core/` had already written.

The same property runs in the other direction. Because every pipeline emits the
same primary key for a variant -- minted once, during mapping, and never
recomputed downstream -- an external tool that produces a `pkey` column joins the
existing feature tables with no integration work and no coordinate
reconciliation. A tool that only needs to read can ignore the Python package
entirely and walk the per-gene directories.

See [Extending BioFeatureFactory](#extending-biofeaturefactory) for the concrete
recipe.

### variant_mapping key options

| Flag | Required | Description |
|------|----------|-------------|
| `-m` / `--mutations` | Yes | Mutation CSV file or directory of `<GENE>_mutations.csv` |
| `-r` / `--reference` | Yes | Indexed reference genome FASTA |
| `-a` / `--annotation` | Yes | GTF, GFF3, or the project's 7-column custom format |
| `-o` / `--output` | No | Output directory (defaults to the working directory) |
| `-Ec` / `--exclude-chromosome` | No | Skip the chromosome mapping CSVs |
| `-Eg` / `--exclude-genomic` | No | Skip the gDNA mapping CSVs |
| `-Et` / `--exclude-transcript` | No | Skip the transcript mapping CSVs |
| `-EA` / `--exclude-aa` | No | Skip the amino-acid mapping CSVs |
| `-Ei` / `--exclude-intron` | No | Skip intron mapping CSVs and per-intron FASTA records |
| `-Ep` / `--exclude-premrna` | No | Skip the pre-mRNA mapping CSVs |
| `--orf` | No | Known ORF FASTA (file or directory). If omitted, ORF is inferred from transcript |
| `-fc` / `--force-cds` | No | Force a transcript: one accession for all genes, or a CSV with `gene,transcript_id` |
| `-v` / `--verbose` | No | Detailed ORF/mutation validation messages plus a log file |

The validation log records out-of-range positions, base mismatches and ORF
derivation details. Downstream helpers (`trim_muts`, `load_validation_failures`)
read it to exclude mutations that failed to map; pass `--log` where a pipeline
accepts it.

---

## Architecture Overview

- **Three layers, one direction.** `lib/` (import-only) $\rightarrow$ `core/` (input
  producers) $\rightarrow$ pipelines (consumers). Documented in `core/__init__.py`. No
  module under `lib/` defines a CLI, and no pipeline imports `core` -- core is
  invoked by command line, and its output files are the interface.
- **One coordinate implementation.** Exon-aware ORF $\leftrightarrow$ transcript $\leftrightarrow$ gDNA
  $\leftrightarrow$ chromosome mapping lives only in `core/variant_mapping.py`.
- **Unified mutation processing.** Insertions, deletions, delins, MNVs, stop
  gains and knockouts are handled by the shared variant grammar in `lib/`, so a
  pipeline gets non-SNV support without implementing it.
- **Directory mode everywhere.** Each pipeline discovers genes from the
  variant_mapping output root and derives the rest; single files and flat
  directories still work unchanged.
- **Multiple processing modes** let you run full pipelines, single tools, or
  parse-only passes.
- **Parallel execution** tuned per tool (`ProcessPoolExecutor`,
  `ThreadPoolExecutor`), with Nextflow orchestration for SpliceAI and the
  mutation-effects backends.
- **Comprehensive logging** so validation failures, warnings and errors are
  traceable per gene and per mutation.

---

## Installation

### Full setup

`scripts/bootstrap.sh` handles the Python environment, the external clones and
source builds, and the prepared databases:

```bash
./scripts/bootstrap.sh                     # everything
./scripts/bootstrap.sh env-phase           # python environment only
./scripts/bootstrap.sh env-phase git-phase # environment + clones/builds
./scripts/bootstrap.sh db-phase            # datasets and prepared DBs
```

Phases are additive -- name any combination and the run is their union. Run
`./scripts/bootstrap.sh --help` for the full list of `--exclude-*` flags.

Databases are built into `Bio_DBs/` by `scripts/build_db.sh`, which bootstrap
calls in `db-phase`. SpliceAI is installed into its own conda environment
(`bff-spliceai`) because it cannot share a process with pyarrow; opt out with
`--exclude-spliceai-env`.

### Python package only

```bash
pip install -e .            # core dependencies
pip install -e ".[all]"     # every optional dependency
```

Core dependencies: `biopython`, `pandas`, `numpy>=1.26.4,<2.5`, `requests`,
`pysam`, `scipy`, `setuptools<81`. Requires Python 3.10-3.13.

Pipeline-specific extras, installable individually or via `[all]`:

- **RNAfold**: `pip install -e ".[rnafold]"` or `conda install -c bioconda viennarna`
- **EVmutation**: `pip install -e ".[evmutation]"` (adds `numba`; plmc is built by bootstrap)
- **adabmDCA**: `pip install -e ".[adabmdca]"` (adds `adabmDCA`, `torch`)
- **NetSurfP3**: `pip install -e ".[netsurfp3]"` (bootstrap clones NetSurfP-3.0 itself)
- **Rare Codon**: `pip install -e ".[rare-codon]"` (adds `networkx`, `prody`)
- **Miranda**: `conda install -c bioconda miranda` (conda only)
- **NetNGlyc with SignalP**: download SignalP 6.0 fast from DTU, then `pip install /path/to/signalp-6-package`

Tests: `pytest test/`.

---

## Pipeline Quick-Start Matrix

| Pipeline | Entry Point | Primary Inputs* | Notable Outputs |
|----------|-------------|-----------------|-----------------|
| **NetNGlyc** | `biofeaturefactory/netNglyc/netnglyc_pipeline.py` | variant_mapping root, SignalP 6.0 | Glycosylation calls with SignalP-aware summaries |
| **NetPhos** | `biofeaturefactory/netphos/netphos_pipeline.py` | variant_mapping root | Kinase-specific phosphorylation predictions |
| **NetMHC** | `biofeaturefactory/netMHC/netmhc_pipeline.py` | variant_mapping root, HLA alleles | MHC binding predictions, epitope gain/loss |
| **NetSurfP3** | `biofeaturefactory/NetSurfP3/netsurfp3_pipeline.py` | variant_mapping root, model + config | RSA, secondary structure, disorder predictions |
| **Miranda** | `biofeaturefactory/miranda/miranda_ensemble.py` | variant_mapping root, miRanda binary, miRNA DB | $\Delta$-based miRNA binding summaries |
| **GeneSplicer** | `biofeaturefactory/genesplicer/genesplicer_ensemble.py` | variant_mapping root, GeneSplicer binary | Donor/acceptor delta summaries |
| **SpliceAI** | `biofeaturefactory/spliceai/spliceai_pipeline_controller.py` | variant_mapping root or VCFs, reference genome, annotation | SpliceAI VCFs, parsed consequence tables |
| **RNAfold** | `biofeaturefactory/RNAfold/run_viennaRNA_pipeline.py` | variant_mapping root | $\Delta\Delta G$ summaries, accessibility deltas |
| **AlphaFold3** | `biofeaturefactory/alphafold3/alphafold3_pipeline.py` | variant_mapping root, RBP binding data (POSTAR3/eCLIP) | RNA-RBP structure predictions, PAE metrics |
| **Mutation Effects** | `biofeaturefactory/mutation_effects/mutEffects_controller.py` | variant_mapping root, protein and/or codon MSA | Epistatic mutation-effect scores (plmc and/or adabmDCA) |
| **Codon Usage** | `biofeaturefactory/codon_usage/codon_usage_pipeline.py` | variant_mapping root | RSCU, CAI, tAI, codon pair scores |
| **Rare Codon** | `biofeaturefactory/rare_codon/rare_codon_pipeline.py` | variant_mapping root (CodonMSA), CoCoPUTs codon table | Rare codon enrichment p-values |

*_Every pipeline also accepts single files or flat directories; naming the
variant_mapping root simply removes the need to pass the mapping paths
separately._

Required arguments have short forms throughout (`--transcript-mapping` $\rightarrow$ `-tm`,
`--reference-codon-usage` $\rightarrow$ `-rcu`); long forms are unchanged.

---

## Pipeline Summaries

### NetNGlyc Pipeline

N-linked glycosylation site prediction with modern SignalP 6.0 integration.

- Requires [NetNGlyc 1.0](https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/) and [SignalP 6.0 fast](https://services.healthtech.dtu.dk/services/SignalP-6.0/) from DTU (academic license).
- SignalP 6.0 integration via `signalp6_adapter` bridges modern SignalP output to NetNGlyc's legacy format.
- Auto-detects installations via `--native-netnglyc-bin` or `NETNGLYC_PATH`/`NETNGLYC_HOME` environment variables.
- Per-mutation glyco calls plus SignalP confidence summaries suitable for modeling.

### NetPhos Pipeline

Kinase-specific phosphorylation site prediction using the NetPhos/APE system.

- Requires [NetPhos 3.1](https://services.healthtech.dtu.dk/services/NetPhos-3.1/) from DTU (academic license).
- Auto-detects APE/NetPhos installations via `--native-ape-path` or `NETPHOS_APE_PATH`/`NETPHOS_HOME` environment variables.
- Produces kinase-aware probability tables, with batch-level caching for repeated runs.
- Reports alignment context and residue projection for non-SNV variants that shift downstream positions.

### NetMHC Pipeline

MHC class I and II binding prediction for WT vs mutant sequences.

- Predicts peptide-MHC binding across multiple HLA alleles.
- Identifies gained/lost epitopes and binding strength changes.
- Supports netMHCpan, netMHC, and netMHCII tools from DTU (academic license).
- Explains skipped tokens rather than dropping them silently.

### NetSurfP3 Pipeline

Protein structure feature prediction using NetSurfP-3.0.

- Predicts relative surface accessibility (RSA), secondary structure, and disorder.
- Compares WT vs mutant structural changes at mutation sites.
- Requires the NetSurfP-3.0 library (cloned by bootstrap) plus a model and config.

### Miranda Pipeline

WT $\leftrightarrow$ MUT miRNA binding analysis with $\Delta$-based metrics.

- Processes each WT transcript once, reuses its miRanda hits for every mutant, then evaluates mutants in parallel.
- Emits summary/events/site tables with $\Delta$-score, competitive binding, and distance-weighted impact metrics.
- Validation-aware filtering plus live progress output for large mutation sets.

### GeneSplicer Pipeline

WT $\leftrightarrow$ ALT ensemble delta caller for splice donor/acceptor sites.

- Runs GeneSplicer on full genomic context (or windowed mode) per gene, compares WT vs ALT, and clusters events.
- Generates summary, events, and sites tables with positional shifts, confidence scoring, and QC flags.
- Writes a per-gene `run_summary.json` beside each gene's TSVs; deterministic batching keeps re-runs reproducible.

### SpliceAI Pipeline

Deep learning-based splice site prediction orchestrated by an adaptive controller.

- `spliceai_pipeline_controller.py` wraps `nextflow run main.nf`, so one command handles VCF generation, SpliceAI inference, and result parsing.
- Runs inside the dedicated `bff-spliceai` environment created by bootstrap; the workflow verifies that environment before submitting work.
- Automatically launches `spliceai-tracker.py` to report per-gene processed/total mutation counts in real time; disable with `--disable_tracker`.
- Supports mutation-driven runs (auto-build per-gene VCFs) and reuse of pre-built VCFs via `--input_vcf_path --skip_vcf_generation`.
- Isoform management: genes exceeding `--max_isoforms_per_gene` (default 50) are filtered via hybrid stratified sampling; override with `--forceAll_isoforms`.
- Parsed output resolves indels against the reference, including left-alignment for tandem-repeat representations.

### RNAfold Pipeline

Secondary structure impact analysis via ViennaRNA.

- Compares reference vs. alternate 151-nt (configurable) windows around each variant.
- Folds exonic, intronic and pre-mRNA contexts; the intron/pre-mRNA mapping is derived from the input root.
- Computes $\Delta\Delta G$, minimum free energy shifts, and Boltzmann ensemble statistics.
- Calculates Jensen-Shannon divergence on base accessibility with configurable $\tau$ thresholding.
- Writes `<gene>/RNAfold/rnafold.run_summary.<timestamp>.json` per invocation.

### AlphaFold3 Pipeline

RNA-RBP interaction structure prediction using AlphaFold3 or Boltz-1.

- Queries POSTAR3/ENCODE eCLIP data for RBPs near mutations.
- Predicts RNA-protein complex structures for WT and mutant alleles.
- Computes PAE-based binding metrics and interface contacts.
- `burst.py` submits and ingests SLURM array jobs against a manifest cache for large batches.

### Mutation Effects Pipeline

Evolutionary mutation effect prediction using epistatic Potts models, with two
backends behind one controller.

- **EVmutation / plmc** -- CPU pseudo-likelihood fit on a protein MSA; requires the EVmutation library and the compiled `plmc` binary (both handled by bootstrap).
- **adabmDCA** -- GPU Boltzmann-machine fit, available on a 21-letter protein alphabet or a 64-codon alphabet.
- `mutEffects_controller.py` selects the backend and orchestrates via Nextflow. Codon-level jobs are sized against detected VRAM; genes whose alignment exceeds the budget route to plmc instead of failing on an out-of-memory error.
- Both backends score the same variants and are interpreted the same way, on different scales.

### Codon Usage Pipeline

Codon optimality metrics for translational efficiency analysis.

- Computes RSCU, CAI, tAI for single codons.
- Calculates codon pair scores (CPS, RSCPU) for bicodons.
- Quantifies changes in translational efficiency.

### Rare Codon Pipeline

Rare codon enrichment detection using sliding window analysis.

- Identifies regions enriched/depleted in rare codons.
- Uses evolutionary conservation from codon-aware MSAs (`<GENE>/CodonMSA/`).
- Scores against the CoCoPUTs human codon-usage table built by `scripts/build_db.sh` rather than the gene's own composition.
- Requires the cg_cotrans library from the Shakhnovich Lab (downloaded by bootstrap).

---

## Shared Infrastructure & Pathing

`biofeaturefactory/lib/` is import-only: no CLI, no side effects, safe to import
from anything.

| Module | Provides |
|---|---|
| `lib/primitives.py` | FASTA I/O, `mint_pkey`, pkey loaders, token helpers, alphabet detection |
| `lib/utility.py` | Variant grammar (`parse_variant`, `canonical_token`, `splice_seq`, `protein_consequence`), mutant sequence construction, mapping loaders, discovery and derivation helpers |
| `lib/annotation.py` | GTF / GFF3 / custom annotation parsing and genome location lookup |
| `lib/msa.py` | jackhmmer execution, Stockholm $\rightarrow$ A2M, sequence weighting, N_eff |
| `lib/codon_metrics.py` | Codon counts, CAI, tAI, bicodon extraction |
| `lib/dtu_outputs.py` | Shared output parsing for NetNGlyc, NetPhos and NetMHC |

Discovery and derivation, used by every pipeline in directory mode:

```python
discover_fasta_files(root)                        # {GENE: fasta path}
discover_mapping_files(root, "transcript")        # {GENE: mapping path}
discover_msa_files(root, kind="protein"|"codon")  # {GENE: msa path}
derive_mapping_root(explicit, input_root, type)   # explicit wins; None if nothing found
derive_mutations_root(explicit, input_root)       # explicit wins; None if nothing found
```

All net-tool pipelines (NetNGlyc, NetPhos, NetMHC, NetSurfP3) require native
binary installations under an academic license.

---

## Extending BioFeatureFactory

The argument for why this works is in [Core as the extension
point](#core-as-the-extension-point). In practice, a new tool needs three
things:

1. **A root.** Take one input path and call `discover_fasta_files` (or
   `discover_msa_files`) to enumerate genes; call `derive_mapping_root` /
   `derive_mutations_root` for the mappings you need.
2. **The variant grammar.** Import what you need from
   `biofeaturefactory.lib.utility` -- `parse_variant`, `splice_seq`,
   `protein_consequence`, `mutation_class` -- rather than parsing tokens
   yourself. Non-SNV handling comes with it.
3. **The key.** Emit `pkey` as your first column, taken from the pkey mapping or
   minted with `mint_pkey(gene, token)`. That is what lets your output sit
   beside the other twelve.

Tools written in other languages do not need to import anything: the per-gene
directory layout is the contract, and it is plain CSV and FASTA.

---

## Citation

### BioFeatureFactory

For citation details, see [CITING.md](CITING.md) and [CITATION.cff](CITATION.cff)

### Pipeline-Specific Citations

- **Miranda** -- Enright AJ, et al. *MicroRNA targets in Drosophila.* Genome Biol. 2003;5(1):R1.
- **NetNGlyc** -- Gupta R, et al. *Prediction of N-glycosylation sites in human proteins.* DTU Health Tech. 2004.
- **NetPhos** -- Blom N, et al. *Sequence and structure-based prediction of eukaryotic protein phosphorylation sites.* J Mol Biol. 1999;294(5):1351-1362.
- **NetMHC** -- Reynisson B, et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0.* Nucleic Acids Res. 2020;48(W1):W449-W454.
- **NetSurfP3** -- Hoie MH, et al. *NetSurfP-3.0: accurate and fast prediction of protein structural features by protein language models and deep learning.* Nucleic Acids Res. 2022;50(W1):W510-W515.
- **GeneSplicer** -- Pertea M, et al. *GeneSplicer: a new computational method for splice site prediction.* Nucleic Acids Res. 2001;29(5):1185-1190.
- **SpliceAI** -- Jaganathan K, et al. *Predicting Splicing from Primary Sequence with Deep Learning.* Cell. 2019;176(3):535-548.
- **RNAfold** -- Lorenz R, et al. *ViennaRNA Package 2.0.* Algorithms Mol Biol. 2011;6(1):26.
- **AlphaFold3** -- Abramson J, et al. *Accurate structure prediction of biomolecular interactions with AlphaFold 3.* Nature. 2024;630:493-500.
- **EVmutation** -- Hopf TA, et al. *Mutation effects predicted from sequence co-variation.* Nat Biotechnol. 2017;35:128-135.
- **adabmDCA** -- Rosset L, et al. *adabmDCA 2.0: a flexible but easy-to-use package for Direct Coupling Analysis.* NAR Genomics and Bioinformatics. 2025. <!-- verify volume/article no. -->
- **Rare Codon** -- Jacobs WM, Shakhnovich EI. *Evidence of evolutionary selection for cotranslational folding.* PNAS. 2017;114:11434-11439.

---

## License

BioFeatureFactory is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0).
See the [LICENSE](LICENSE) file for details.

---

## Support

For issues and questions, open a ticket at [GitHub Issues](https://github.com/jgoldmintz/BioFeatureFactory/issues).
