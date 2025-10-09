# BioFeatureFactory
Python toolkit for automated feature extraction using gene sequences and SNP data
## Architecture

### Shared Infrastructure
- **docker/**: Shared Docker environment for both NetNGlyc and NetPhos pipelines
- **dependencies/**: Common utility functions and shared processing logic (see below)

#### Exon-Aware Mapping Workflow (`dependencies/exon_aware_mapping.py`)
`exon_aware_mapping.py` prepares  exon-aware assets by reading in per-gene mutation CSV files together with an annotation file and reference genome to generate:
- Three-in-one FASTA bundles (`ORF`, `transcript`, `genomic`) for each gene
- Coordinate mapping CSVs in chromosome, genomic-slice, and transcript space (`combined_{GENE}.csv`)
- Optional validation logs summarising length/mismatch issues
- Timestamped verbose validation reports (when `--verbose` is enabled) capturing per-mutation failures for later review

Typical usage:
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
Run it whenever mutation sets change to keep downstream pipelines in sync; the emitted CSVs supply the mutation→coordinate linkage used for pkey generation and filtering.

### Key Features
- **Unified mutation processing**: Single-mutation logic for accurate mutant sequence analysis
- **Shared utilities**: Common functions for both NetPhos and NetNGlyc pipelines  
- **Multiple processing modes**: full-pipeline, tool-only, parse-only
- **Parallel processing**: Intelligent worker management with ProcessPoolExecutor
- **Apple Silicon support**: Docker-based solutions for ARM64 compatibility
- **Comprehensive logging**: Error handling and processing validation

## Available Pipelines

### **NetNGlyc Pipeline** (`netnglyc/`)
N-linked glycosylation site prediction with SignalP 6.0 integration.
- Docker-based NetNGlyc 1.0 with modern SignalP 6.0 predictions
- Simplified single-file-per-gene processing architecture
- Batch combination and intelligent processing mode selection
- Single-mutation processing for accurate mutant analysis

### **NetPhos Pipeline** (`netphos/`)
Phosphorylation site prediction for serine, threonine, and tyrosine residues.
- Docker-based NetPhos 3.1 with APE system integration
- Single-mutation processing eliminates over-matching
- Kinase-specific predictions with confidence scoring
- Unified mutation filtering logic

### **Miranda Pipeline** (`miranda/`)
miRNA target site prediction software suite.
- Timeout protection for long-running analyses
- Raw file cleanup options for storage management

### **GeneSplicer Pipeline** (`genesplicer/`)
Splice site prediction with ThreadPoolExecutor optimization.
- Temporary file management with automatic cleanup

### **SpliceAI Pipeline** (`spliceai/`)
Deep learning-based splice site prediction with Nextflow automation.
- End-to-end processing from mutation CSV files to parsed results
- TensorFlow race condition mitigation for stable parallel execution
- Exon-aware coordinate mapping with pkey generation
- RefSeq chromosome format support with annotation conversion

## Upcoming:

- netSufP3
- EVmutation
- netMHC
- RNAfold

## Citation

If you use BioFeatureFactory or any of its pipelines, please cite:

### BioFeatureFactory
- Goldmintz J. BioFeatureFactory: Python toolkit for automated biological feature extraction. GitHub repository. https://github.com/jgoldmintz/BioFeatureFactory

### Pipeline-Specific Citations

#### Miranda Pipeline
- Enright AJ, John B, Gaul U, Tuschl T, Sander C, Marks DS. MicroRNA targets in Drosophila. Genome Biol. 2003;5(1):R1.

#### NetNGlyc Pipeline
- Gupta R, Jung E, Brunak S. Prediction of N-glycosylation sites in human proteins. NetNglyc 1.0 Server. 2004. Available from: https://services.healthtech.dtu.dk/service.php?NetNGlyc-1.0
- Teufel F, Almagro Armenteros JJ, Johansen AR, Gíslason MH, Pihl SI, Tsirigos KD, Winther O, Brunak S, von Heijne G, Nielsen H. SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol. 2022;40(7):1023-1025.

#### NetPhos Pipeline
- Blom N, Gammeltoft S, Brunak S. Sequence and structure-based prediction of eukaryotic protein phosphorylation sites. J Mol Biol. 1999;294(5):1351-62.

#### GeneSplicer Pipeline
- Pertea M, Lin X, Salzberg SL. GeneSplicer: a new computational method for splice site prediction. Nucleic Acids Res. 2001;29(5):1185-90.

#### SpliceAI Pipeline
- Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
