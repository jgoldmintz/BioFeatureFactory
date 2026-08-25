# AlphaFold3 RNA-RBP Interaction Pipeline

For each variant, finds RBPs with POSTAR3/ENCODE eCLIP peaks near the site, runs AlphaFold3 on the WT and MUT RNA-protein complexes, and reports per-RBP binding deltas.

## Requirements

| Component | Notes |
|-----------|-------|
| AlphaFold3 Docker image | Built locally; see setup below |
| NVIDIA GPU | With CUDA drivers and the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) |
| Docker | Required for local mode |
| Python >= 3.9 | With `pysam` for tabix queries. Uses `biofeaturefactory/lib/`. |
| POSTAR3 database | Tabix-indexed BED, passed with `-pd` |
| RBP sequences or MSAs | `-rs` FASTA, or `-md` A3M directory (preferred) |

### AlphaFold3 setup

Follow the [AF3 installation guide](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)
through **Obtaining Model Parameters**. Stop before **Obtaining Genetic Databases** -- this pipeline
does not use AF3's genetic database pipeline; RBP binding sites come from POSTAR3 (`-pd`) and
protein sequences or MSAs are supplied directly (`-md` / `-rs`).

```bash
cd alphafold3
docker build -t alphafold3 -f docker/Dockerfile .
```

Pass the weights directory with `-mdi/--model-dir`; it is mounted into the container at `/root/models`.

## Usage

```bash
# Directory mode: variant_mapping output root supplies mutations, mappings, and VCFs
python alphafold3_pipeline.py \
    -f out/ -o af3_results/ \
    -pd RBP_db/human-POSTAR3.sorted.bed.gz \
    -rm human_uniprot_genes.tsv \
    -md msa_files/ -mdi /path/to/af3_weights

# File mode, single gene
python alphafold3_pipeline.py \
    -f transcripts/SMN2.fasta \
    -cm mappings/SMN2_chromosome_mapping.csv \
    -pd RBP_db/human-POSTAR3.sorted.bed.gz \
    -rm human_uniprot_genes.tsv -rs rbp_sequences.fasta \
    -mdi /path/to/af3_weights -o af3_results/
```

In directory mode `-f` is the `variant_mapping` output root and supplies `-mu`, `-cm`, `-pm`, and
`-v` from `<root>/<GENE>/`.

## Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-pd, --postar-db` | Yes | -- | Tabix-indexed POSTAR3 BED |
| `-rm, --rbp-mapping` | Yes | -- | Gene-UniProt mapping TSV |
| `-o, --output` | Yes | -- | Output base directory |
| `-f, --fasta` | No | -- | variant_mapping output root, single transcript FASTA, or flat directory |
| `-rs, --rbp-sequences` | * | -- | Protein sequences FASTA |
| `-md, --msa-dir` | * | -- | Directory of A3M MSA files (preferred over `-rs`) |
| `-mu, --mutations` | No | derived | File mode only; mutations CSV |
| `-cm, --chromosome-mapping` | No | derived | File mode only; chromosome mapping CSV |
| `-pm, --premrna-mapping` | No | derived | File mode only. Required to score intronic (`gd.`) variants, which have no ORF or transcript coordinate and are windowed on the pre-mRNA record |
| `-v, --vcf` | No | derived | File mode only; per-gene VCF from `vcf_converter.py` (provides the chromosome) |
| `-ch, --chrom` | No | -- | Chromosome; alternative to `-v` |
| `-ts, --tx-start` | No | -- | Transcript start position; alternative to `-cm` |
| `-s, --strand` | No | `+` | `+` or `-` |
| `-em, --execution-mode` | No | `local` | Only `local` is supported here; use `burst.py` for SLURM |
| `-ab, --af3-binary` | No | `alphafold3` | Path to the AF3 executable |
| `-di, --docker-image` | No | `alphafold3` | Docker image name |
| `-mdi, --model-dir` | No | -- | AF3 model weights directory (required in local mode) |
| `-ws, --window-size` | No | `101` | RNA window size in nt (odd) |
| `-rw, --rbp-window` | No | `50` | Window for RBP site lookup (+/- bp) |
| `-mw, --multi-window` | No | off | Run multiple windows per mutation |
| `-mwo, --multi-window-offsets` | No | `0.3,0.5,0.7` | Mutation position as a fraction of the window |
| `-mg, --max-gpus` | No | auto | Max GPUs for parallel AF3 execution |
| `-vl, --validation-log` | No | -- | Validation log for filtering mutations |

\* One of `-rs` or `-md` is required. Mutations come from `-mu` or `-cm`.

## Output

```
{output}/{GENE}/AlphaFold3/
    {GENE}.tsv          -- per-mutation summary
    {GENE}.events.tsv   -- per-RBP binding comparison
    {GENE}.sites.tsv    -- per-residue interface data
```

### `{GENE}.tsv`

| Column | Description |
|--------|-------------|
| `pkey`, `Gene` | Identifiers |
| `n_rbps_tested` | RBPs evaluated for this mutation |
| `n_rbps_binding_wt`, `n_rbps_binding_mut` | RBPs with confident binding per allele |
| `global_count_gained`, `global_count_lost` | RBPs whose binding appeared / disappeared |
| `global_count_strengthened`, `global_count_weakened` | Both bind; MUT PAE lower / higher |
| `global_max_abs_delta_pae` | Largest absolute PAE change (angstrom) |
| `top_event_rbp` | RBP with the highest priority score |
| `top_event_class` | gained / lost / strengthened / weakened / none |
| `top_event_delta_pae` | PAE change for that RBP (angstrom) |
| `qc_flags` | See below |

### `{GENE}.events.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `rbp_name` | RBP gene symbol |
| `wt_chain_pair_pae_min`, `mut_chain_pair_pae_min` | Minimum inter-chain PAE per allele (angstrom) |
| `delta_chain_pair_pae_min` | MUT - WT (angstrom) |
| `wt_interface_contacts`, `mut_interface_contacts` | Cross-chain atom pairs within 8 angstrom |
| `delta_interface_contacts` | MUT - WT |
| `cls` | gained / lost / strengthened / weakened / unchanged / no_binding / incomplete |
| `priority` | Ranking score: absolute delta PAE, hyperbolically decayed by distance to the variant |
| `n_samples_wt`, `n_samples_mut` | AF3 samples parsed per allele (ensemble mode) |
| `std_pae_wt`, `std_pae_mut` | SD of chain-pair PAE across samples (angstrom) |
| `n_windows` | Windows used (multi-window mode only) |

### `{GENE}.sites.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `rbp_name` | RBP identifier |
| `allele` | WT or MUT |
| `chain` | Chain identifier (RNA / protein) |
| `res_id`, `res_name` | Residue number and 3-letter name |
| `plddt` | Per-residue pLDDT (0-100) |
| `is_contact` | 1 if within 8 angstrom of the other chain |
| `min_contact_distance` | Nearest atom in the other chain (angstrom) |
| `contact_frequency` | Fraction of ensemble samples where this residue is a contact (0-1) |

### Sign conventions

| Metric | Positive | Negative |
|--------|----------|----------|
| `delta_chain_pair_pae_min` | Weaker binding in MUT | Stronger binding in MUT |
| `delta_interface_contacts` | More contacts in MUT | Fewer contacts in MUT |

## Classification thresholds

Not CLI arguments. The first five are fields of `ThresholdConfig` in
`bin/binding_metrics.py`; `contact_threshold` is a parameter of the contact
extraction in `bin/af3_parser.py`.

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `min_contacts` | 3 | Minimum interface contacts for confident binding |
| `max_pae_binding` | 10.0 | Maximum inter-chain PAE for confident binding (angstrom) |
| `min_plddt_interface` | 50.0 | Minimum pLDDT at the interface |
| `delta_pae_significant` | 2.0 | Threshold for strengthened / weakened (angstrom) |
| `delta_contacts_significant` | 2 | Contact-count threshold for classification |
| `contact_threshold` | 8.0 | Distance defining an interface contact (angstrom) |

Binding counts as **confident** when all three hold: at least `min_contacts` contacts,
inter-chain PAE at or below `max_pae_binding`, and at least one chain at the interface with
pLDDT at or above `min_plddt_interface`.

## QC flags

| Flag | Condition |
|------|-----------|
| `PASS` | All RBPs produced complete results |
| `PARTIAL` | Some RBPs succeeded, some failed |
| `ALL_FAILED` | Every RBP prediction failed |
| `no_rbps_tested` | No RBPs were evaluated |
| `no_rbps_in_region` | No RBP peaks within the query window |
| `FAILED:{error}` | An individual RBP prediction raised an exception |

## Ensemble and multi-window aggregation

AF3 emits several samples per seed. The pipeline parses every `seed-N_sample-N/` subdirectory and
reports the mean and population SD of chain-pair PAE, interface contacts, and pLDDT, plus per-residue
contact frequency. With a single sample it falls back to single-model behaviour and leaves the
`n_samples_*` / `std_pae_*` columns empty.

`-mw/--multi-window` places the mutation at each offset in `-mwo` instead of centering it, which
reduces sensitivity to window boundaries. Each offset produces a distinct WT/MUT pair (deduplicated
near transcript ends), AF3 runs per window with jobs suffixed `_w0`, `_w1`, ..., and metrics are
aggregated across windows. Cost is `3 offsets x 2 alleles x N RBPs` AF3 jobs per mutation, so it is
off by default.

## SLURM execution

`burst.py` is a separate two-phase driver; `alphafold3_pipeline.py` is local/Docker only.
Submit generates input JSONs, dedupes by `AF3Input.get_hash()`, writes a manifest, renders a SLURM
array script, and `sbatch`es it, then exits. Ingest walks the cache, computes deltas, and writes the
per-gene TSVs.

```bash
python -m biofeaturefactory.alphafold3.burst submit \
    --fasta transcripts/ --chromosome-mapping mappings/ --vcf vcf_files/ \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping rbp_uniprot_ids.txt --msa-dir msa_files/ \
    --output af3_results/ --model-dir /path/to/af3_weights \
    --slurm-partition gpu --slurm-time 01:00:00 --slurm-mem 64G \
    --array-throttle 256

# once sacct shows all tasks COMPLETED/FAILED
python -m biofeaturefactory.alphafold3.burst ingest \
    --output af3_results/ --fasta transcripts/ \
    --chromosome-mapping mappings/ --vcf vcf_files/ \
    --postar-db RBP_db/human-POSTAR3.sorted.bed.gz \
    --rbp-mapping rbp_uniprot_ids.txt --msa-dir msa_files/
```

`--execution-mode batch` and `cloud` were removed from `alphafold3_pipeline.py`; the stub generated
broken artifacts and never ingested results.

## Caching

WT predictions are cached by `{gene}-{rbp}` hash and reused across that gene's mutations. RBP
sequences are loaded once at startup. Cache location: `{output}/cache/`.

## Module reference

| File | Purpose |
|------|---------|
| `alphafold3_pipeline.py` | Local entry point |
| `burst.py` | SLURM submit / ingest driver |
| `bin/af3_runner.py` | AF3 execution |
| `bin/af3_parser.py` | Parses mmCIF and JSON outputs |
| `bin/rbp_database.py` | POSTAR3 tabix query interface |
| `bin/rbp_sequence_mapper.py` | RBP name -> UniProt -> sequence |
| `bin/binding_metrics.py` | Delta computation and classification |

## License

AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.

Issues: https://github.com/jgoldmintz/BioFeatureFactory/issues
