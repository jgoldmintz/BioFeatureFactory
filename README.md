# BioFeatureFactory
Python toolkit for automated feature extraction using gene sequences and SNP data

## Repository Key Features:

- **Multiple processing modes** (full-pipeline, tool-only, parse-only)
- **Parallel processing** with intelligent worker management
- **Comprehensive logging** and error handling
- **Resume capability** and caching for interrupted runs
- **Automatic cleanup** for large datasets
- **Mutation position filtering** using mapping files

## Currently Available:

### **miranda pipeline**
miRNA target site prediction software suite.
- Timeout protection for long-running analyses
- Raw file cleanup options for storage management

### **netNglyc**
Docker-based NetNGlyc 1.0 pipeline with SignalP 6.0 integration.
- Apple Silicon compatibility via Docker
- SignalP 6.0 integration for modern signal peptide predictions
- Intelligent processing mode auto-detection based on sequence count

### **netPhos**
Docker-based NetPhos 3.1 pipeline for phosphorylation site prediction.
- Apple Silicon compatibility via Docker APE system integration
- Per-file processing optimization with intelligent caching
- Serine/threonine/tyrosine phosphorylation predictions with kinase specificity

### **GeneSplicer**
GeneSplicer pipeline for splice site prediction.
- ThreadPoolExecutor optimization
- Temporary file management with automatic cleanup

## Upcoming:

- netSufP3
- EVmutation
- spliceAI
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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For issues and questions:

- GitHub Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
