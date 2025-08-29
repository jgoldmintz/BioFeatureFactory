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

### **GeneSplicer**
GeneSplicer pipeline for splice site prediction.
- ThreadPoolExecutor optimization
- Temporary file management with automatic cleanup

## Upcoming:

- netSufP3
- netPhos
- EVmutation
- spliceAI
- netMHC
