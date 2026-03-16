# Utils

Shared utilities and standalone pipelines for coordinate mapping, VCF conversion, MSA generation, and common bioinformatics operations.

---

## Table of Contents

1. [exon_aware_mapping.py](#exon_aware_mappingpy)
2. [vcf_converter.py](#vcf_converterpy)
3. [codon_msa_pipeline.py](#codon_msa_pipelinepy)
4. [msa_generation_pipeline.py](#msa_generation_pipelinepy)
5. [utility.py](#utilitypy)

---

## exon_aware_mapping.py

Generates exon-aware coordinate mappings and FASTA sequences for genes. Produces ORF (coding sequence), transcript (spliced mRNA), and genomic (full tx_start..tx_end) sequences along with multi-coordinate-system mapping CSVs.

### Usage

```bash
python3 exon_aware_mapping.py \
  --mutations /path/to/mutations/ \
  --annotation /path/to/annotations.gtf \
  --reference /path/to/reference_genome.fa \
  --output /path/to/output_dir/ \
  --orf /path/to/orf_fastas/ \
  --force-cds transcript_overrides.csv \
  --verbose
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--mutations` / `-m` | Yes | Mutation CSV file or directory of CSVs |
| `--reference` / `-r` | Yes | Indexed reference genome FASTA |
| `--annotation` / `-a` | Yes | Gene annotation file (GTF, GFF3, or custom) |
| `--output` / `-o` | No | Output directory (defaults to current working directory) |
| `--exclude-chromosome` / `-Ec` | No | Skip chromosome mapping CSV output |
| `--exclude-genomic` / `-Eg` | No | Skip gDNA mapping CSV output (relative to gene slice) |
| `--exclude-transcript` / `-Et` | No | Skip transcript coordinate mapping output |
| `--exclude-aa` / `-EA` | No | Skip amino-acid mapping output |
| `--orf` | No | ORF FASTA file or directory. If omitted, ORF is inferred from transcript |
| `--force-cds` | No | Force specific transcript. Single accession (applied to all genes) or CSV with `gene,transcript_id` columns |
| `--verbose` | No | Detailed ORF/mutation validation messages and log file |

### Outputs

- **FASTAs**: Per-gene ORF, transcript, and genomic sequences.
- **Chromosome mapping CSVs**: `mutant,chromosome_coordinate` with absolute genomic positions.
- **Genomic mapping CSVs**: Positions relative to gene genomic slice.
- **Transcript mapping CSVs**: Positions relative to spliced transcript.
- **AA mapping CSVs**: Amino-acid change notation.
- **Validation log**: Timestamped file listing out-of-range positions, base mismatches, and ORF derivation details.

### Key Functions

| Function | Purpose |
|----------|---------|
| `build_transcript_seq_and_map` | Builds transcript sequence, exon structure, CDS boundaries, and coordinate mappings from annotation + reference |
| `map_orf_mutations_to_transcript_and_genome` | Projects ORF-relative mutations into transcript, absolute genomic, relative genomic, and AA coordinates |
| `derive_orf_from_transcript` | Auto-detects ORF by finding the longest in-frame ATG-to-STOP region matching mutation positions |
| `resolve_orf_from_sources` | Chooses between supplied ORF or derived ORF, runs validation |
| `validate_mutations_against_orf` | Checks mutation positions and reference bases against ORF sequence |

---

## vcf_converter.py

Converts gene-specific SNP CSV files into standard VCF format with RefSeq chromosome naming. Includes exon-aware mapping support, reference sequence verification, and a JSON-based caching system to avoid redundant regeneration.

### Usage

```bash
python3 vcf_converter.py \
  --mutation /path/to/mutations/ \
  --outdir /path/to/vcfs/ \
  --reference /path/to/reference_genome.fa \
  --annotation /path/to/annotations.gtf \
  --chromosome-format refseq \
  --chromosome-mapping-input /path/to/chromosome_mappings/ \
  --log /path/to/validation_logs/ \
  --verify-sequences /path/to/fastas/
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--mutation` | Yes | Mutation CSV or directory |
| `--outdir` | Yes | Output directory for VCF files |
| `--reference` | Yes | Reference genome FASTA |
| `--annotation` | Yes | Gene annotation file |
| `--chromosome-format` | No | `simple`, `refseq` (default), or `ucsc` |
| `--chromosome-mapping-input` | No | Pre-computed chromosome mapping CSV or directory |
| `--validate-mapping` | No | Cross-check chromosome mappings against reference FASTA |
| `--log` | No | Validation log or directory, used to skip known-bad mutations |
| `--clear-cache` | No | Force regeneration by clearing cached VCF metadata |
| `--verify-sequences` | No | FASTA file or directory to verify against reference |

### Outputs

- **VCF files**: One per gene, sorted by position, with standard VCF 4.2 headers.
- **Cache file**: `.vcf_converter_cache.json` in output directory for incremental runs.
- **Summary statistics**: Printed to stdout.

### Key Functions

| Function | Purpose |
|----------|---------|
| `process_single_file` | Converts one gene's mutations CSV to VCF using mapping lookup or annotation-driven computation |
| `chromosome_to_refseq` | Converts chromosome number to RefSeq accession (e.g., `1` -> `NC_000001.11`) |
| `get_reference_nucleotide` | Extracts reference base at a given genomic position, handling multiple chromosome naming conventions |
| `verify_sequences_against_reference` | Validates FASTA sequences against reference genome |

---

## codon_msa_pipeline.py

Builds codon-aware MSAs by running mmseqs2 searches against RefSeq assemblies, aligning proteins with mafft/muscle, and back-translating to codon-level alignments. Gaps in the protein alignment are expanded to codon-triplet gaps (`---`), preserving alignment structure for evolutionary analysis (e.g., rare codon enrichment).

### Usage

```bash
python3 codon_msa_pipeline.py \
  --fasta /path/to/gene.fasta \
  --db-root /path/to/Bio_DBs/ \
  --output /path/to/output_dir/ \
  --aligner mafft \
  --verbose
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--fasta` | Yes | ORF nucleotide FASTA (file or directory) |
| `--db-root` | Yes | Root directory containing RefSeq databases |
| `--output` / `-o` | Yes | Output directory |
| `--min-seqid` | No | Minimum sequence identity for final filtering (default: 1.0) |
| `--search-min-seqid` | No | Minimum sequence identity for mmseqs2 search (default: 0.5) |
| `--search-min-qcov` | No | Minimum query coverage for mmseqs2 search (default: 0.5) |
| `--max-hits` | No | Maximum number of hits to retain (default: 500) |
| `--aligner` | No | Protein aligner: `mafft` (default) or `muscle` |
| `--mmseqs-binary` | No | Path to mmseqs2 binary (default: `mmseqs`) |
| `--threads` | No | CPU threads |
| `--target-db-base` | No | Pre-built mmseqs2 target database |
| `--mmseqs-tmp-dir` | No | Temporary directory for mmseqs2 |
| `--verbose` | No | Verbose output |

### Outputs

- **Codon MSA FASTA**: `{GENE}/CodonMSA/{GENE}.codon.msa.fasta` — aligned nucleotide sequences with `---` for gaps.
- **Manifest TSV**: `.manifest.tsv` with processing metadata.

---

## msa_generation_pipeline.py

Generates protein MSAs for EVmutation/PLMC analysis using iterative jackhmmer search against UniRef90. Outputs A2M format with quality assessment metrics.

Accepts either a protein or nucleotide FASTA as `--query`. Nucleotide sequences are automatically detected by alphabet composition and translated to protein before the jackhmmer search.

### Usage

```bash
python3 msa_generation_pipeline.py \
  --query /path/to/query.fasta \
  --database /path/to/uniref90.fasta \
  --jackhmmer-binary /path/to/jackhmmer \
  --output /path/to/output.a2m \
  --iterations 5 \
  --threads 4
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--query` / `-q` | Yes | Query FASTA -- protein or nucleotide (auto-detected; nucleotide is translated before search) |
| `--database` / `-d` | Yes | Sequence database (e.g., UniRef90) |
| `--jackhmmer-binary` / `-j` | Yes | Path to jackhmmer executable |
| `--output` / `-o` | Yes | Output MSA file (A2M format) |
| `--iterations` / `-N` | No | jackhmmer iterations (default: 5) |
| `--evalue-inclusion` / `-E` | No | E-value inclusion threshold (default: 1e-3) |
| `--bitscore-threshold` | No | Minimum bits per residue (default: 0.5) |
| `--min-neff-ratio` | No | Minimum N_eff / L ratio (default: 10) |
| `--max-seq-gaps` | No | Max gap fraction per sequence (default: 0.4) |
| `--max-col-gaps` | No | Max gap fraction per column (default: 0.6) |
| `--identity-threshold` | No | Clustering threshold for N_eff (default: 0.8) |
| `--threads` / `-t` | No | CPU threads (default: 4) |
| `--keep-intermediate` | No | Retain intermediate Stockholm file |

### Outputs

- **A2M MSA**: Uppercase = match states, lowercase = insertions, `-` = deletions, `.` = query insertions.
- **Stats JSON**: `<output>.stats.json` with N_eff, coverage, mean identity, and quality pass/fail.

---

## utility.py

Central utility library used by all BioFeatureFactory pipelines. Not a standalone script.

### Function Categories

- **Annotation Parsing** — GTF/GFF3/custom annotation loading, chromosome name normalization, exon/CDS extraction
- **FASTA I/O** — Reading and writing FASTA files with configurable header parsing and duplicate handling
- **Mutation Handling** — Parsing mutation strings, loading mutation CSVs, validation log filtering, mutation-to-amino-acid inference
- **Mapping & Discovery** — Loading mapping CSVs, auto-discovering gene-specific FASTA and mapping files in directories
- **WT/Mutant FASTA Synthesis** — Building mutant protein sequences from WT nucleotide + mapping, batch FASTA generation
- **Prediction Filtering** — Filtering NetPhos/NetNGlyc predictions by mutation position, threshold, and failure maps
- **Codon Usage** — RSCU, CAI, tAI, codon pair scores (CPS, RSCPU), relative adaptiveness (W, W_CP)
- **MSA Analysis** — Alphabet detection, protein query preparation, Stockholm-to-A2M conversion, gap filtering, N_eff computation
- **Batch Processing** — Splitting FASTAs into batches, combining NetNGlyc/NetPhos batch outputs
- **Output Writing** — TSV writing, A2M writing, mutation classification
