# Utils

Shared utilities and standalone pipelines for coordinate mapping, VCF conversion, MSA generation, and common bioinformatics operations.

---

## Table of Contents

1. [exon_aware_mapping.py](#exon_aware_mappingpy)
2. [vcf_converter.py](#vcf_converterpy)
3. [codon-msa-pipeline.py](#codon-msa-pipelinepy)
4. [msa-generation-pipeline.py](#msa-generation-pipelinepy)
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
  --out-fasta /path/to/output_fastas/ \
  --out-chromosome-mapping /path/to/chromosome_mappings/ \
  --out-genomic-mapping /path/to/genomic_mappings/ \
  --out-transcript-mapping /path/to/transcript_mappings/ \
  --out-aa-mapping /path/to/aa_mappings/ \
  --orf /path/to/orf_fastas/ \
  --force-cds transcript_overrides.csv \
  --verbose
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--mutations` | Yes | Mutation CSV file or directory of CSVs |
| `--reference` | Yes | Indexed reference genome FASTA |
| `--annotation` | Yes | Gene annotation file (GTF, GFF3, or custom) |
| `--out-fasta` | Yes | Output directory for per-gene FASTAs |
| `--out-chromosome-mapping` | Yes | Output directory for absolute genomic coordinate mappings |
| `--out-genomic-mapping` | No | Output directory for gDNA mapping CSVs (relative to gene slice) |
| `--out-transcript-mapping` | No | Output directory for transcript coordinate mappings |
| `--out-aa-mapping` | No | Output directory for amino-acid mappings |
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

## codon-msa-pipeline.py

Back-translates a protein MSA to a codon-aware MSA using corresponding nucleotide sequences. Gaps in the protein alignment are expanded to codon-triplet gaps (`---`), preserving alignment structure for evolutionary analysis (e.g., rare codon enrichment).

### Usage

```bash
python3 codon-msa-pipeline.py \
  --protein-msa /path/to/protein.msa.fasta \
  --nucleotide-dir /path/to/nt_sequences/ \
  --output /path/to/codon.msa.fasta \
  --focus HUMAN_QUERY \
  --validate
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--protein-msa` | Yes | Protein MSA FASTA |
| `--nucleotide-dir` | No* | Directory of nucleotide FASTAs |
| `--nucleotide-file` | No* | Single FASTA with all nucleotide sequences |
| `--focus` | No | Focus sequence ID (placed first in output) |
| `--output` | No | Output path (default: `<input>_codon.msa.fasta`) |
| `--validate` | No | Validate codon MSA translates back to protein MSA |

*One of `--nucleotide-dir` or `--nucleotide-file` is required.

### Outputs

- **Codon MSA FASTA**: Aligned nucleotide sequences with `---` for gaps, maintaining column correspondence with the protein MSA.
- **Validation report**: Printed to stdout when `--validate` is used.

### Key Functions

| Function | Purpose |
|----------|---------|
| `backtranslate_protein_to_codons` | Maps aligned protein positions to unaligned nucleotide codons, inserting `---` gaps |
| `validate_codon_msa` | Confirms codon MSA translates back to the original protein MSA |

---

## msa-generation-pipeline.py

Generates protein MSAs for EVmutation/PLMC analysis using iterative jackhmmer search against UniRef90. Outputs A2M format with quality assessment metrics.

### Usage

```bash
python3 msa-generation-pipeline.py \
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
| `--query` | Yes | Query protein FASTA |
| `--database` | Yes | Sequence database (e.g., UniRef90) |
| `--jackhmmer-binary` | Yes | Path to jackhmmer executable |
| `--output` | Yes | Output MSA file (A2M format) |
| `--iterations` | No | jackhmmer iterations (default: 5) |
| `--evalue-inclusion` | No | E-value inclusion threshold (default: 1e-3) |
| `--bitscore-threshold` | No | Minimum bits per residue (default: 0.5) |
| `--min-neff-ratio` | No | Minimum N_eff / L ratio (default: 10) |
| `--max-seq-gaps` | No | Max gap fraction per sequence (default: 0.4) |
| `--max-col-gaps` | No | Max gap fraction per column (default: 0.6) |
| `--identity-threshold` | No | Clustering threshold for N_eff (default: 0.8) |
| `--threads` | No | CPU threads (default: 4) |
| `--keep-intermediate` | No | Retain intermediate Stockholm file |

### Outputs

- **A2M MSA**: Uppercase = match states, lowercase = insertions, `-` = deletions, `.` = query insertions.
- **Stats JSON**: `<output>.stats.json` with N_eff, coverage, mean identity, and quality pass/fail.

### Key Functions

| Function | Purpose |
|----------|---------|
| `run_jackhmmer` | Executes jackhmmer iterative search |
| `stockholm_to_a2m` | Converts Stockholm MSA to A2M format with match/insert state annotation |
| `filter_msa_by_gaps` | Removes gappy sequences then gappy columns |
| `compute_neff` | Calculates effective sequence count via clustering-based weighting |

---

## utility.py

Central utility library used by all BioFeatureFactory pipelines. Not a standalone script.

### Annotation Parsing

| Function | Purpose |
|----------|---------|
| `get_genome_loc(genename, annotation_file, assembly, transcript_id)` | Returns `{chrom, strand, tx_start, tx_end, exons, transcript_id}` from GTF/GFF3/custom annotation. Auto-detects format. |
| `_normalize_chrom_name(raw_name, assembly)` | Normalizes chromosome labels across naming conventions (`NC_000001`, `chr1`, `1` -> `1`) |

### FASTA I/O

| Function | Purpose |
|----------|---------|
| `read_fasta(inf, aformat, duplicate)` | Loads FASTA with format modes: `NCBI` (full header), `FIRST` (first word), `WORD` (second word). `duplicate` controls collision handling. |
| `write_fasta(path, name_to_seq)` | Writes sequences with 60-character line wrapping |

### Mutation Handling

| Function | Purpose |
|----------|---------|
| `get_mutation_data(ntposnt)` | Parses mutation string (e.g., `G123A`) to 0-based position and nucleotides |
| `get_mutation_data_bioAccurate(ntposnt)` | Same but returns 1-based position, skips stop codons |
| `trim_muts(ntPosnt, log, gene_name)` | Loads mutations from CSV, optionally filters by validation logs |
| `extract_mutation_from_sequence_name(seq_name)` | Extracts `(gene, mutation_id)` from sequence names like `ZFP36-C330T` |

### Mapping & Discovery

| Function | Purpose |
|----------|---------|
| `load_mapping(mapping_file, mapType)` | Loads 2-column CSV mapping. `mapType` selects column: `'mutant'` or `'mapped'` |
| `discover_mapping_files(mapping_dir)` | Scans directory for mapping CSVs, returns `{GENE: path}` |
| `discover_fasta_files(fasta_dir)` | Scans directory for FASTA files, returns `{GENE: path}` |
| `extract_gene_from_filename(filename)` | Extracts gene name from filenames using pattern matching |

### WT/Mutant FASTA Synthesis

| Function | Purpose |
|----------|---------|
| `build_mutant_sequences_for_gene` | Generates mutant protein sequences from WT nucleotide sequence + mapping file |
| `synthesize_gene_fastas` | Creates WT and mutant FASTAs for all genes in batch |
| `translate_orf_sequence(nt_sequence)` | Translates nucleotide ORF to amino acid sequence |
| `infer_aamutation_from_nt(mutant_id, nt_sequence)` | Infers amino-acid change from nucleotide mutation |

### Codon Usage

| Function | Purpose |
|----------|---------|
| `get_codon_counts(seq)` | Computes RSCU, W, RSCPU, CPS for codons and codon pairs |
| `compute_cai(seq, w_values)` | Codon Adaptation Index (geometric mean of codon weights) |
| `compute_tai(seq, tai_weights)` | tRNA Adaptation Index |
| `extract_codon_with_bicodons(ntposnt, seq)` | Extracts codon and flanking bicodons at a mutation site |

### MSA Analysis

| Function | Purpose |
|----------|---------|
| `stockholm_to_a2m` | Converts Stockholm MSA to A2M format |
| `filter_msa_by_gaps` | Filters gappy sequences and columns |
| `compute_neff` | Effective sequence count via clustering |
| `validate_msa_quality` | Quality metrics: N_eff ratio, coverage, mean identity |

### Batch Processing

| Function | Purpose |
|----------|---------|
| `split_fasta_into_batches(fasta_file, batch_size, temp_dir)` | Splits FASTA into chunks for parallel processing |
| `combine_batch_outputs(batch_output_files, final_output_file, format_type)` | Combines batch results (netnglyc/netphos formats) |

### Reference Genome

| Function | Purpose |
|----------|---------|
| `download_reference_genome(build, cache_dir)` | Downloads and decompresses GRCh37/GRCh38 from NCBI |

### Other

| Function | Purpose |
|----------|---------|
| `convert_position(seq1, seq2, position1, space)` | Projects positions between gapped sequence alignments |
| `retry_request(func, positional_arguments, keyword_arguments, lim, wait)` | Retries HTTP requests with configurable limit and backoff |
| `load_validation_failures(log_path)` | Loads validation failure maps for mutation filtering |