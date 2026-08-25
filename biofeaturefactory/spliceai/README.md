# SpliceAI Pipeline

Nextflow-orchestrated SpliceAI splice-site prediction: mutation CSVs to VCF, SpliceAI execution, then parsing to per-gene TSV with pkey mapping.

## Requirements

| Component | Notes |
|-----------|-------|
| Nextflow | 20.04+ with DSL2. `curl -s https://get.nextflow.io \| bash` |
| SpliceAI | In the active conda environment (`spliceai -h`). Needs `numpy<2` to match conda-forge TensorFlow. |
| Reference genome | Indexed FASTA (e.g. `GRCh38_reference.fna`) |
| Annotation | Tab-delimited SpliceAI annotation from `annot_to_spliceai.py` |
| Mapping CSVs | From `core/variant_mapping.py` |
| Resources | 8 GB+ RAM, multi-core CPU |

### Component files

| File | Purpose |
|------|---------|
| `spliceai_pipeline_controller.py` | Entry point; validates inputs, launches Nextflow, runs the tracker |
| `bin/main.nf` | Nextflow workflow |
| `bin/filter_annotation.py` | Isoform filter for high-isoform genes |
| `bin/spliceai-parser.py` | VCF parser with pkey mapping and log filtering |
| `annot_to_spliceai.py` | GTF to SpliceAI annotation converter |
| `../core/vcf_converter.py` | Mutation token to VCF converter |

### Generating the annotation file

The annotation's chromosome naming must match the VCF's.

```bash
python3 annot_to_spliceai.py reference.gtf --chromosome-format ncbi -o annotation_ncbi.txt
# also: --chromosome-format simple | ucsc
```

## Usage

```bash
# Build VCFs from mutation CSVs, then run
python spliceai_pipeline_controller.py \
    -mp out/ \
    -rg /path/to/reference.fna \
    -af annotation_ncbi.txt \
    -od results/

# Reuse existing VCFs
python spliceai_pipeline_controller.py \
    -ivp vcfs/ -svg \
    -rg /path/to/reference.fna \
    -af annotation_ncbi.txt \
    -od results/
```

Mutation CSVs must be named `<gene_id>_mutations.csv` (case-insensitive); the gene ID seeds the
downstream channels. Given the `variant_mapping` output root, the transcript, chromosome, and
genomic mapping paths are all derived from it -- supply them only for mappings outside that layout.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-rg, --reference_genome` | Yes | -- | Indexed reference FASTA |
| `-af, --annotation_file` | Yes | -- | SpliceAI annotation file |
| `-mp, --mutations_path` | * | -- | Directory or file with `<GENE>_mutations.csv` |
| `-ivp, --input_vcf_path` | * | -- | Existing VCF or directory of VCFs |
| `-svg, --skip_vcf_generation` | No | off | Use pre-built VCFs; requires `-ivp` |
| `-od, --output_dir` | No | -- | Output base directory |
| `-vod, --vcf_output_dir` | No | -- | Publish intermediate VCFs, bgzip files, and indexes here |
| `-tmp, --transcript_mapping_path` | No | derived | File mode only |
| `-cmp, --chromosome_mapping_path` | No | derived | File mode only |
| `-gmp, --genomic_mapping_path` | No | derived | File mode only; `bin/main.nf` no longer requires it |
| `-vl, --validation_log` | No | -- | Validation log for filtering failed mutations |
| `-cf, --chromosome_format` | No | `refseq` | `refseq`, `simple`, or `ucsc` |
| `-st, --splice_threshold` | No | `0.0` | Minimum delta score passed to the parser |
| `-vm, --validate_mapping` | No | off | Cross-check chromosome mappings against the reference |
| `-cvc, --clear_vcf_cache` | No | off | Force VCF regeneration |
| `-ma, --maxforks` | No | `0` | Max concurrent `run_spliceai` tasks; 0 = Nextflow default |
| `-fi, --forceAll_isoforms` | No | off | Process every isoform, skipping the filter |
| `-mipg, --max_isoforms_per_gene` | No | `50` | Isoform count above which filtering applies |
| `-ps, --poll_seconds` | No | `15.0` | Tracker poll interval |
| `-dt, --disable_tracker` | No | off | Skip launching `spliceai-tracker.py` |
| `-pi, --pipeline` | No | `bin/main.nf` | Nextflow script to run |
| `-re, --resume` | No | off | Resume a previous run via Nextflow `-resume` |

\* One of `-mp` or `-ivp` is required.

Genes above `-mipg` isoforms are reduced by hybrid stratified sampling: the 10 longest isoforms
(capturing canonical/MANE transcripts) plus 40 sampled from the remainder, seeded deterministically
on the gene name. Without it, genes such as BRCA1 (368 isoforms) can run for days.

## Output

```
{output_dir}/{GENE}/SpliceAI/{GENE}.tsv
```

One row per variant per isoform block.

| Column | Description | Units |
|--------|-------------|-------|
| `pkey` | `GENE-MUTATION_ID` from the transcript mapping | string |
| `gene` | Gene symbol | string |
| `chrom` | Chromosome, in the annotation's format | string |
| `pos` | Genomic position (1-based) | integer |
| `ref`, `alt` | Reference and alternate alleles | string |
| `allele` | Allele reported by SpliceAI | string |
| `block_label` | Isoform block: first four unique score vectors are `A`-`D`, the rest `dup` | A/B/C/D/dup |
| `ds_ag`, `ds_al` | Delta score, acceptor gain / loss | 0-1 |
| `ds_dg`, `ds_dl` | Delta score, donor gain / loss | 0-1 |
| `dp_ag`, `dp_al`, `dp_dg`, `dp_dl` | Distance to each predicted site; negative = upstream, 0 = none predicted | nt |
| `max_delta_score` | Maximum of the four delta scores | 0-1 |

SpliceAI evaluates each variant across transcript isoforms and emits one block per isoform. The
parser preserves that order; isoforms whose predictions match an earlier block exactly are labelled
`dup`.

## Limitations

SpliceAI declines to score a variant when both `ref` and `alt` are longer than one base (true
delins and MNVs). It emits `ALT|GENE|.|.|.|.|.|.|.|.` before the model runs, so those rows carry
dots rather than scores. Insertions and deletions written in minimal VCF form (one allele of length
1) score normally -- `core/vcf_converter.py` left-normalizes and trims to that form.

## Running components standalone

```bash
# Parse an existing SpliceAI VCF
python3 bin/spliceai-parser.py \
    --input ABCB1.spliceai.vcf --output ABCB1.tsv \
    --chromosome-mapping chromosome_mapping/combined_ABCB1.csv \
    --transcript-mapping transcript_mapping/combined_ABCB1.csv \
    --threshold 0.1 --log validation.log

# Monitor a run started elsewhere
python3 bin/spliceai-tracker.py \
    --mutations-path /path/to/mutations_dir \
    --work-root /path/to/work --poll-seconds 5 --log-file tracker.log
```

The controller launches and cleans up the tracker itself; run it manually only to observe a
pre-existing run, and stop it with Ctrl+C.

## Troubleshooting

| Symptom | Resolution |
|---------|------------|
| No predictions generated | Chromosome naming mismatch: the VCF uses `NC_000007.14` but the annotation uses `7`. Regenerate with a matching `--chromosome-format` |
| `ImportError: numpy.core.umath failed to import` | pip numpy is ABI-incompatible with conda-forge TensorFlow: `pip install "numpy<2"`. Set `-ma 1` if crashes persist |
| No pkeys generated | Mapping files missing; regenerate with `core/variant_mapping.py` |
| A gene runs for hours with no progress | High isoform count. Filtering is on by default; check the work directory's `.command.log` for `[FILTER]` lines |
| Rows are dots, not scores | The variant is a delins or MNV; see Limitations |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.

Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
