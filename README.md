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
- **Cross-platform support** via Docker images and helper scripts for Linux, macOS (Apple Silicon), and Windows via WSL2.  
- **Comprehensive logging** ensures validation, warnings, and errors are traceable.  

---

## Pipeline Quick-Start Matrix

| Pipeline | Entry Point | Primary Inputs*                                                                   | Notable Outputs                                                                       |
|----------|-------------|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| **NetNGlyc** | `netNglyc/run_netnglyc_pipeline.py` | Gene FASTA bundle, `combined_<GENE>.csv`, SignalP 6.0 model                       | Glycosylation site predictions with SignalP-enhanced summaries                        |
| **NetPhos** | `netphos/run_netphos_pipeline.py` | Gene FASTA bundle, `combined_<GENE>.csv`                                          | Kinase-specific phosphorylation predictions with confidence scores                    |
| **Miranda** | `miranda/run_miranda_pipeline.py` | miRNA FASTA, transcript FASTA, `combined_<GENE>.csv`                              | Candidate miRNA binding sites with timeout-protected runs                             |
| **GeneSplicer** | `genesplicer/run_genesplicer_pipeline.py` | Genomic FASTA slices, `combined_<GENE>.csv`                                       | Splice donor/acceptor predictions per mutation window                                 |
| **SpliceAI** | `spliceai/run_spliceai_nextflow.sh` | Mutation CSV(s), reference genome, annotation GTF                                 | SpliceAI VCFs, parsed consequence tables, exon-aware logs                             |
| **RNAfold** | `RNAfold/run_viennaRNA_pipeline.py` | Transcript FASTA, `combined_<GENE>.csv` (`--transcript-mapping`), reference FASTA | $ΔΔG$ summaries, Jensen–Shannon divergence metrics, per-position accessibility deltas |

Each script auto-discovers dependencies via `dependencies/utility.py` if the repository layout is preserved.

*_Either single files or directories (for bulk processing) are accepted inputs_

---

## Pipeline Summaries

### NetNGlyc Pipeline

N-linked glycosylation site prediction with optional SignalP 6.0 integration.

- Dockerized NetNGlyc 1.0 wrapped with modern SignalP predictions.  
- Single-mutation processing for precise context windows.  
- Batch orchestration with automatic mode selection for full runs vs. reprocessing.  
- Aggregated results include SignalP confidence alongside glycosylation calls.  

### NetPhos Pipeline

Phosphorylation site prediction for serine, threonine, and tyrosine residues.

- Dockerized NetPhos 3.1 bundled with the APE scoring system.  
- Per-mutation windowing prevents overscoring due to nearby variants.  
- Kinase-specific probability outputs and consolidated CSV summaries.  
- Shared filtering logic mirrors NetNGlyc mutation handling.  

### Miranda Pipeline

miRNA target site prediction suite.

- Wraps the `miranda` executable with timeout protection and retry logic.  
- Offers raw output retention or cleanup depending on storage constraints.  
- Produces per-mutation binding site calls ready for downstream filtering.  

### GeneSplicer Pipeline

Splice site prediction using GeneSplicer.

- ThreadPoolExecutor-managed batching keeps CPU utilization steady.  
- Temporary file handling cleans up intermediate outputs automatically.  
- Emits donor/acceptor score tables aligned with exon-aware coordinates.  

### SpliceAI Pipeline

Deep learning-based splice site prediction with Nextflow automation.

- Converts mutation CSVs into per-gene VCFs and launches SpliceAI.  
- Handles TensorFlow race conditions for reproducible parallel runs.  
- Supports RefSeq-to-Ensembl chromosome conversions via exon-aware mappings.  
- Provides both raw SpliceAI scores and parsed annotations.  

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

Prepare mutation-driven analyses once, then reuse generated assets across pipelines.

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
   - Emits chromosome, genomic-slice, and transcript mapping CSVs (`combined_<GENE>.csv` variants).  
   - Logs validation issues (missing exons, mismatched coordinates).  

2. **Consume validation logs automatically**

   - Downstream helpers (e.g., `trim_muts`) consult the exon-aware log directory to exclude problematic mutations.  
   - Pass `--log /path/to/log_dir` where applicable.  

3. **Optional: reference checks and VCF conversion**

   - `dependencies/vcf_converter.py` cross-validates FASTA headers against the reference genome, normalizes chromosome naming, and creates sorted per-gene VCFs for SpliceAI.  
   - Run when VCF input is required or genome mismatches are suspected.  

4. **Stage outputs for pipelines**

   - Place FASTA bundles and `combined_<GENE>.csv` files into directories used by each pipeline.  
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

If you use BioFeatureFactory or any of its pipelines, please cite:

### BioFeatureFactory
Goldmintz J. *BioFeatureFactory: Python toolkit for automated biological feature extraction.* GitHub repository. <https://github.com/jgoldmintz/BioFeatureFactory>

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

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Support

For issues and questions, open a ticket at [GitHub Issues](https://github.com/jgoldmintz/BioFeatureFactory/issues).
