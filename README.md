# BioFeatureFactory

Modular bioinformatics framework for automated feature extraction, coordinate mapping, and predictive modeling using gene- and variant-level inputs.

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [Installation](#installation)
3. [Pipeline Quick-Start Matrix](#pipeline-quick-start-matrix)
4. [Pipeline Summaries](#pipeline-summaries)
   - [NetNGlyc](#netnglyc-pipeline)
   - [NetPhos](#netphos-pipeline)
   - [NetMHC](#netmhc-pipeline)
   - [NetSurfP3](#netsurfp3-pipeline)
   - [Miranda](#miranda-pipeline)
   - [GeneSplicer](#genesplicer-pipeline)
   - [SpliceAI](#spliceai-pipeline)
   - [RNAfold](#rnafold-pipeline)
   - [AlphaFold3](#alphafold3-pipeline)
   - [EVmutation](#evmutation-pipeline)
   - [Codon Usage](#codon-usage-pipeline)
   - [Rare Codon](#rare-codon-pipeline)
5. [Shared Infrastructure & Pathing](#shared-infrastructure--pathing)
6. [End-to-End Data Preparation Workflow](#end-to-end-data-preparation-workflow)
7. [Citation](#citation)
8. [License](#license)
9. [Support](#support)

---

## Architecture Overview

- **Unified mutation processing** keeps mutant sequence analysis consistent across pipelines.  
- **Shared utilities** in `utils/utility.py` provide discovery helpers, mapping loaders, and mutation filters for every pipeline.  
- **Multiple processing modes** let you run full pipelines, single tools, or parse-only passes.  
- **Parallel execution** leverages multiprocessing/threading tuned per tool (`ProcessPoolExecutor`, `ThreadPoolExecutor`).  
- **Cross-platform support** via Docker images and helper scripts for Linux and macOS (Apple Silicon).  
- **Comprehensive logging** ensures validation, warnings, and errors are traceable.  

---

## Pipeline Quick-Start Matrix

| Pipeline | Entry Point | Primary Inputs* | Notable Outputs |
|----------|-------------|-----------------|-----------------|
| **NetNGlyc** | `netNglyc/netnglyc-pipeline.py` | WT & mutant protein FASTA, mapping CSV dir, SignalP 6.0 | Glycosylation calls with SignalP-aware summaries |
| **NetPhos** | `netphos/netphos-pipeline.py` | WT/mutant protein FASTA, mapping CSV dir | Kinase-specific phosphorylation predictions |
| **NetMHC** | `netMHC/netmhc-pipeline.py` | WT/mutant protein FASTA, mapping CSV dir, HLA alleles | MHC binding predictions, epitope gain/loss |
| **NetSurfP3** | `NetSurfP3/netsurfp3-pipeline.py` | WT/mutant protein FASTA, mapping CSV dir | RSA, secondary structure, disorder predictions |
| **Miranda** | `miranda/miranda-ensemble.py` | WT transcript FASTA, MirandA binary, miRNA DB | Δ-based miRNA binding summaries |
| **GeneSplicer** | `genesplicer/genesplicer_ensemble.py` | Genomic FASTA, mutation CSV dir, GeneSplicer binary | Donor/acceptor delta summaries |
| **SpliceAI** | `spliceai/spliceai-pipeline-controller.py` | Mutation CSVs or VCFs, reference genome | SpliceAI VCFs, parsed consequence tables |
| **RNAfold** | `RNAfold/run_viennaRNA_pipeline.py` | Transcript FASTA, transcript-mapping CSV dir | ΔΔG summaries, accessibility deltas |
| **AlphaFold3** | `alphafold3/alphafold3_pipeline.py` | Mutation CSV, RBP binding data (POSTAR3/eCLIP) | RNA-RBP structure predictions, PAE metrics |
| **EVmutation** | `EVmutation/evmutation-pipeline.py` | Protein MSA, plmc binary | Epistatic mutation effect scores |
| **Codon Usage** | `codon_usage/codon-usage-pipeline.py` | ORF FASTA, mutation CSV | RSCU, CAI, tAI, codon pair scores |
| **Rare Codon** | `rare_codon/rare-codon-pipeline.py` | Codon MSA, codon usage file | Rare codon enrichment p-values |

Each script auto-discovers dependencies via `utils/utility.py` if the repository layout is preserved.

*_Either single files or directories (for bulk processing) are accepted inputs_

---

## Pipeline Summaries

### NetNGlyc Pipeline

N-linked glycosylation site prediction with modern SignalP 6.0 integration.

- Dockerized NetNGlyc 1.0/SignalP 6.0 stack that runs on Intel + Apple Silicon.  
- On Linux hosts the pipeline auto-detects native `netNglyc` installations (`--native-netnglyc-bin`, `--force-native`) and falls back to Docker otherwise.  
- Intelligent FASTA/mapping discovery (any filename/extension) for WT + mutant protein sets.  
- Per-mutation glyco calls plus SignalP confidence summaries suitable for modeling.  

### NetPhos Pipeline

Kinase-specific phosphorylation site prediction using the NetPhos/APE system.

- On Linux hosts the pipeline auto-detects native APE/NetPhos installations (`--native-ape-path`, `--force-native`), while macOS/Windows reuse the shared NetNGlyc container.  
- Flexible FASTA/mapping discovery identical to NetNGlyc; supports WT vs mutant comparisons.  
- Produces kinase-aware probability tables with caching for repeated runs.  

### Miranda Pipeline

WT<->MUT miRNA binding analysis with Δ-based metrics.

- Processes each WT transcript once, reuses its MirandA hits for every mutant, then evaluates mutants in parallel.  
- Emits summary/events/site tables with Δ-score, competitive binding, and distance-weighted impact metrics.  
- Validation-aware filtering plus live progress output for large mutation sets.  

### GeneSplicer Pipeline

WT<->ALT ensemble delta caller for splice donor/acceptor sites.

- Runs GeneSplicer on full genomic context (or windowed mode) per gene, compares WT vs ALT, and clusters events.  
- Generates summary, events, and sites tables with positional shifts, confidence scoring, and QC flags.  
- Deterministic batching and append-safe outputs enable reproducible re-runs.  

### SpliceAI Pipeline

Deep learning-based splice site prediction orchestrated by an adaptive controller.

- `spliceai-pipeline-controller.py` wraps `nextflow run main.nf`, so one command handles VCF generation, SpliceAI inference, and result parsing.  
- Automatically launches the live tracker and monitors `.nextflow.log` for rapid `exit 134` events; when ≤6 genes thrash, it restarts with `--maxforks tail_maxforks` (default 1) to finish deterministically.  
- Supports both mutation-driven runs (auto-build per-gene VCFs) and reuse of pre-built VCFs via `--input_vcf_path --skip_vcf_generation`.  
- Outputs raw SpliceAI VCFs plus exon-aware parsed consequence tables for downstream modeling.  

### RNAfold Pipeline

Secondary structure impact analysis via ViennaRNA.

- Compares reference vs. alternate 151-nt (configurable) windows around each variant.
- Computes ΔΔG, minimum free energy shifts, and Boltzmann ensemble statistics.
- Calculates Jensen-Shannon divergence on base accessibility with configurable τ thresholding.
- Outputs per-mutation ΔΔG, MFE change flags, divergence metrics, and ensemble statistics.

### NetMHC Pipeline

MHC class I and II binding prediction for WT vs mutant sequences.

- Predicts peptide-MHC binding across multiple HLA alleles.
- Identifies gained/lost epitopes and binding strength changes.
- Supports netMHCpan, netMHC, and netMHCII tools via Docker.

### NetSurfP3 Pipeline

Protein structure feature prediction using NetSurfP-3.0.

- Predicts relative surface accessibility (RSA), secondary structure, and disorder.
- Compares WT vs mutant structural changes at mutation sites.
- Requires cloning NetSurfP-3.0 library from GitHub.

### AlphaFold3 Pipeline

RNA-RBP interaction structure prediction using AlphaFold3 or Boltz-1.

- Queries POSTAR3/ENCODE eCLIP data for RBPs near mutations.
- Predicts RNA-protein complex structures for WT and mutant alleles.
- Computes PAE-based binding metrics and interface contacts.

### EVmutation Pipeline

Evolutionary mutation effect prediction using epistatic models.

- Fits Potts model to protein MSA via plmc.
- Predicts fitness effects of amino acid substitutions.
- Requires cloning EVmutation library from GitHub.

### Codon Usage Pipeline

Codon optimality metrics for translational efficiency analysis.

- Computes RSCU, CAI, tAI for single codons.
- Calculates codon pair scores (CPS, RSCPU) for bicodons.
- Quantifies changes in translational efficiency.

### Rare Codon Pipeline

Rare codon enrichment detection using sliding window analysis.

- Identifies regions enriched/depleted in rare codons.
- Uses evolutionary conservation from codon-aware MSAs.
- Requires downloading cg_cotrans library from Shakhnovich Lab.

---

## Installation

From the repository root, install the shared utilities as an editable package:

```bash
pip install -e .
```

This registers the `utils` package so all pipelines can import shared helpers without path manipulation. Re-run after pulling updates that modify `pyproject.toml`.

---

## Shared Infrastructure & Pathing

- `utils/utility.py` serves as the core helper module for NetNGlyc, NetPhos, Miranda, GeneSplicer, SpliceAI preprocessors, and RNAfold.
- The `utils` package is installed via `pip install -e .` (see [Installation](#installation)).
- Helpers include directory discovery (`get_input_dir`, `get_output_dir`), mutation CSV loaders, exon-aware mapping filters, FASTA retrieval utilities, and logging wrappers.
- `docker/` contains Dockerfiles and compose snippets for NetNGlyc and NetPhos.

---

## End-to-End Data Preparation Workflow

For the vast majority of the pipelines in this repository properly mapped mutations are required for biologically accurate predictions. Therefore running the `exon_aware_mapping.py` first is highly recommended. 

1. **Generate exon-aware assets**

```bash
   python3 utils/exon_aware_mapping.py \
     --mutations /path/to/mutations/ \
     --annotation /path/to/annotations.gtf \
     --reference /path/to/reference_genome.fa \
     --out-fasta /path/to/output_fastas/ \
     --out-chromosome-mapping /path/to/chromosome_mappings/ \
     --out-genomic-mapping /path/to/genomic_mappings/ \
     --out-transcript-mapping /path/to/transcript_mappings/ \
     --out-aa-mapping /path/to/aa_mappings/ \
     --orf /path/to/orf_fastas/ \
     --force-cds transcript_overrides.csv \
     --verbose
```

   **Key options:**
   - `--orf` — Directory or file containing known ORF sequences. If omitted, ORF is inferred from transcript.
   - `--force-cds` — Force specific transcript isoforms. Accepts either:
     - Single accession (e.g., `NM_022162.3`) applied to all genes
     - CSV file with `gene,transcript_id` columns for per-gene overrides
   - `--out-aa-mapping` — Optional output directory for amino-acid mapping CSVs.
   - `--verbose` — Print detailed ORF/mutation validation messages and write log file.

   **Outputs:**
   - Per-gene ORF, transcript, and genomic FASTA bundles.
   - Chromosome, genomic-slice, transcript, and amino-acid mapping CSVs (`*{GENE}*.csv`).
   - Validation log with mutation issues (out-of-range positions, base mismatches).  

2. **Consume validation logs automatically**

   - Downstream helpers (e.g., `trim_muts`) consult the exon-aware log directory to exclude problematic mutations.  
   - Pass `--log /path/to/log_dir` where applicable.  

3. **Optional: reference checks and VCF conversion**

   - `utils/vcf_converter.py` cross-validates FASTA headers against the reference genome, normalizes chromosome naming, and creates sorted per-gene VCFs for SpliceAI.  
   - Run when VCF input is required or genome mismatches are suspected.  

4. **Stage outputs for pipelines**

   - Place FASTA bundles and `*{GENE}*.csv` mapping files into directories used by each pipeline.  
   - NetNGlyc and NetPhos need ORF FASTAs and mapping CSVs.  
   - RNAfold also requires transcript mappings when run with `--transcript-mapping`.  

5. **Launch pipelines with consistent inputs**

   - All pipelines use shared exon-aware assets, avoiding redundant preprocessing.  

---

## Citation

####

### BioFeatureFactory
For citation details, see [CITING.md](CITING.md) and [CITATION.cff](CITATION.cff)

### Pipeline-Specific Citations

- **Miranda** — Enright AJ, et al. *MicroRNA targets in Drosophila.* Genome Biol. 2003;5(1):R1.
- **NetNGlyc** — Gupta R, et al. *Prediction of N-glycosylation sites in human proteins.* DTU Health Tech. 2004.
- **NetPhos** — Blom N, et al. *Sequence and structure-based prediction of eukaryotic protein phosphorylation sites.* J Mol Biol. 1999;294(5):1351-1362.
- **NetMHC** — Reynisson B, et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0.* Nucleic Acids Res. 2020;48(W1):W449-W454.
- **NetSurfP3** — Klausen MS, et al. *NetSurfP-2.0: Improved prediction of protein structural features.* Proteins. 2019;87(6):520-527.
- **GeneSplicer** — Pertea M, et al. *GeneSplicer: a new computational method for splice site prediction.* Nucleic Acids Res. 2001;29(5):1185-1190.
- **SpliceAI** — Jaganathan K, et al. *Predicting Splicing from Primary Sequence with Deep Learning.* Cell. 2019;176(3):535-548.
- **RNAfold** — Lorenz R, et al. *ViennaRNA Package 2.0.* Algorithms Mol Biol. 2011;6(1):26.
- **AlphaFold3** — Abramson J, et al. *Accurate structure prediction of biomolecular interactions with AlphaFold 3.* Nature. 2024;630:493-500.
- **EVmutation** — Hopf TA, et al. *Mutation effects predicted from sequence co-variation.* Nat Biotechnol. 2017;35:128-135.
- **Rare Codon** — Jacobs WM, Shakhnovich EI. *Evidence of evolutionary selection for cotranslational folding.* PNAS. 2017;114:11434-11439.  

---

## License

BioFeatureFactory is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0).
See the [LICENSE](LICENSE) file for details.

---

## Support

For issues and questions, open a ticket at [GitHub Issues](https://github.com/jgoldmintz/BioFeatureFactory/issues).
