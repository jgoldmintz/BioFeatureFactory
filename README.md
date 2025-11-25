# BioFeatureFactory

Modular bioinformatics framework for automated feature extraction, coordinate mapping, and predictive modeling using gene- and variant-level inputs.

## Table of Contents
1. [Architecture Overview](#architecture-overview)
2. [Pipeline Quick-Start Matrix](#pipeline-quick-start-matrix)
3. [Pipeline Summaries](#pipeline-summaries)
   - [NetNGlyc](#netnglyc-pipeline)
   - [NetPhos](#netphos-pipeline)
   - [Miranda](#miranda-pipeline)
   - [GeneSplicer](#genesplicer-pipeline)
   - [SpliceAI](#spliceai-pipeline)
   - [RNAfold](#rnafold-pipeline)
4. [Shared Infrastructure & Pathing](#shared-infrastructure--pathing)
5. [End-to-End Data Preparation Workflow](#end-to-end-data-preparation-workflow)
6. [Upcoming Pipelines](#upcoming-pipelines)
7. [Citation](#citation)
8. [License](#license)
9. [Support](#support)

---

## Architecture Overview

- **Unified mutation processing** keeps mutant sequence analysis consistent across pipelines.  
- **Shared utilities** in `dependencies/utility.py` provide discovery helpers, mapping loaders, and mutation filters for every pipeline.  
- **Multiple processing modes** let you run full pipelines, single tools, or parse-only passes.  
- **Parallel execution** leverages multiprocessing/threading tuned per tool (`ProcessPoolExecutor`, `ThreadPoolExecutor`).  
- **Cross-platform support** via Docker images and helper scripts for Linux and macOS (Apple Silicon).  
- **Comprehensive logging** ensures validation, warnings, and errors are traceable.  

---

## Pipeline Quick-Start Matrix

| Pipeline | Entry Point | Primary Inputs*                                                                   | Notable Outputs                                                                       |
|----------|-------------|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| **NetNGlyc** | `netNglyc/full-docker-netnglyc-pipeline.py` | WT & mutant protein FASTA dirs, flexible mapping CSV dir, SignalP 6.0 install      | Glycosylation calls with SignalP-aware summaries                                      |
| **NetPhos** | `netphos/netphos-pipeline.py` | WT/mutant protein FASTA dirs plus mapping CSV dir (Dockerized NetPhos/APE)       | Kinase-specific phosphorylation predictions with cached summaries                     |
| **Miranda** | `miranda/run_miranda_pipeline.py` | WT transcript FASTA, MirandA binary, miRNA DB, optional mutation metadata/logs     | Δ-based miRNA binding summaries, events, and per-site audits                          |
| **GeneSplicer** | `genesplicer/genesplicer_ensemble.py` | Genomic FASTA slices, mutation CSV dir, GeneSplicer binary & models                | Donor/acceptor delta summaries, event tables, detailed site audits                    |
| **SpliceAI** | `spliceai/spliceai-pipeline-controller.py` | Mutation CSVs or pre-built VCFs, reference genome, SpliceAI annotation             | SpliceAI VCFs, parsed consequence tables, adaptive restart + tracking logs            |
| **RNAfold** | `RNAfold/run_viennaRNA_pipeline.py` | Transcript FASTA, transcript-mapping CSV dir, reference FASTA                     | $ΔΔG$ summaries, Jensen–Shannon divergence metrics, per-position accessibility deltas |

Each script auto-discovers dependencies via `dependencies/utility.py` if the repository layout is preserved.

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

WT↔MUT miRNA binding analysis with Δ-based metrics.

- Processes each WT transcript once, reuses its MirandA hits for every mutant, then evaluates mutants in parallel.  
- Emits summary/events/site tables with Δ-score, competitive binding, and distance-weighted impact metrics.  
- Validation-aware filtering plus live progress output for large mutation sets.  

### GeneSplicer Pipeline

WT↔ALT ensemble delta caller for splice donor/acceptor sites.

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
- Computes $ΔΔG$, minimum free energy shifts, and Boltzmann ensemble statistics.  
- Calculates Jensen–Shannon divergence on base accessibility ($Δu$) with configurable τ thresholding.  
- Outputs:
  - **Summary**: per-mutation $ΔΔG$, MFE change flags, divergence metrics, and ensemble statistics.  
  - **Ensemble Summary**: per-position $Δu$ values highlighting nucleotides crossing $τ$. 
- Designed for transcript-level mappings; invoke with `--transcript-mapping` to feed exon-aware outputs directly.  

---

## Shared Infrastructure & Pathing

- `dependencies/utility.py` serves as the core helper module for NetNGlyc, NetPhos, Miranda, GeneSplicer, SpliceAI preprocessors, and RNAfold.  
- Each pipeline appends `../dependencies` to `sys.path` at runtime.  
- Helpers include directory discovery (`get_input_dir`, `get_output_dir`), mutation CSV loaders, exon-aware mapping filters, FASTA retrieval utilities, and logging wrappers.  
- `docker/` contains Dockerfiles and compose snippets for NetNGlyc and NetPhos. 

---

## End-to-End Data Preparation Workflow

For the vast majority of the pipelines in this repository properly mapped mutations are required for biologically accurate predictions. Therefore running the `exon_aware_mapping.py` first is highly recommended. 

1. **Generate exon-aware assets**

```bash
   python3 dependencies/exon_aware_mapping.py \
     --mutations /path/to/mutations/ \
     --annotation /path/to/annotations.gtf \
     --reference /path/to/reference_genome.fa \
     --out-fasta /path/to/output_fastas/ \
     --out-chromosome-mapping /path/to/chromosome_mappings/ \
     --out-genomic-mapping /path/to/genomic_mappings/ \
     --out-transcript-mapping /path/to/transcript_mappings/
```

   - Produces per-gene ORF, transcript, and genomic FASTA bundles.  
   - Emits chromosome, genomic-slice, and transcript mapping CSVs (`*{GENE}*.csv`).  
   - Logs validation issues (missing exons, mismatched coordinates).  

2. **Consume validation logs automatically**

   - Downstream helpers (e.g., `trim_muts`) consult the exon-aware log directory to exclude problematic mutations.  
   - Pass `--log /path/to/log_dir` where applicable.  

3. **Optional: reference checks and VCF conversion**

   - `dependencies/vcf_converter.py` cross-validates FASTA headers against the reference genome, normalizes chromosome naming, and creates sorted per-gene VCFs for SpliceAI.  
   - Run when VCF input is required or genome mismatches are suspected.  

4. **Stage outputs for pipelines**

   - Place FASTA bundles and `*{GENE}*.csv` mapping files into directories used by each pipeline.  
   - NetNGlyc and NetPhos need ORF FASTAs and mapping CSVs.  
   - RNAfold also requires transcript mappings when run with `--transcript-mapping`.  

5. **Launch pipelines with consistent inputs**

   - All pipelines use shared exon-aware assets, avoiding redundant preprocessing.  

---

## Upcoming Pipelines

- netSufP3  
- EVmutation  
- netMHC  

---

## Citation

####

### BioFeatureFactory
For citation details, see [CITING.md](CITING.md) and [CITATION.cff](CITATION.cff)

### Pipeline-Specific Citations

- **Miranda** — Enright AJ, John B, Gaul U, Tuschl T, Sander C, Marks DS. *MicroRNA targets in Drosophila.* Genome Biol. 2003;5(1):R1.  
- **NetNGlyc** — Gupta R, Jung E, Brunak S. *Prediction of N-glycosylation sites in human proteins.* NetNglyc 1.0 Server. 2004. <https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0>  
  Teufel F, Almagro Armenteros JJ, Johansen AR, et al. *SignalP 6.0 predicts all five types of signal peptides using protein language models.* Nat Biotechnol. 2022;40(7):1023–1025.  
- **NetPhos** — Blom N, Gammeltoft S, Brunak S. *Sequence and structure-based prediction of eukaryotic protein phosphorylation sites.* J Mol Biol. 1999;294(5):1351–1362.  
- **GeneSplicer** — Pertea M, Lin X, Salzberg SL. *GeneSplicer: a new computational method for splice site prediction.* Nucleic Acids Res. 2001;29(5):1185–1190.  
- **SpliceAI** — Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. *Predicting Splicing from Primary Sequence with Deep Learning.* Cell. 2019;176(3):535–548.e24.  
- **RNAfold** — Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. *ViennaRNA Package 2.0.* Algorithms Mol Biol. 2011;6(1):26.  

---

## License

BioFeatureFactory is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0).
See the [LICENSE](LICENSE) file for details.

---

## Support

For issues and questions, open a ticket at [GitHub Issues](https://github.com/jgoldmintz/BioFeatureFactory/issues).
