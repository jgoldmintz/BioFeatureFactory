# scripts/

Setup and database build scripts for BioFeatureFactory. Run these once before
using the pipelines.

---

## Files

### `bootstrap.sh`

Installs all BioFeatureFactory pipeline dependencies into the current Python
environment and clones/builds required third-party tools.

```bash
cd scripts/
./bootstrap.sh              # Full install (all steps)
./bootstrap.sh pip-only     # Only pip install
./bootstrap.sh git-only     # Only git clones and builds
./bootstrap.sh db-only      # Only database downloads
```

**Subcommands:**

| Subcommand | Scope |
|------------|-------|
| _(none)_ | Everything (steps 1-12) |
| `pip-only` | Only pip install (step 2) |
| `git-only` | Only git clones, builds, and conda installs (steps 3-8) |
| `db-only` | Only FTP downloads and build\_db (steps 9, 12) |

**What it does (in order):**

| Step | Action |
|------|--------|
| 1 | Validates `git` and `tar` are on PATH |
| 2 | `pip install -e ".[all]"` (editable install with all optional deps) + pyarrow ABI fix |
| 3 | Clones EVmutation and plmc; compiles plmc. Downloads cg\_cotrans tarball into `rare_codon/cg_cotrans/` |
| 4 | Clones NetSurfP-3.0 into `NetSurfP3/nsp3/`; pip-installs its requirements |
| 5 | Clones SignalP 6.0 into `netNglyc/signalp-6.0/`; pip-installs its requirements |
| 6 | Installs miranda via conda/mamba (`conda install -c bioconda miranda`) |
| 7 | Downloads and builds GeneSplicer from source into `genesplicer/GeneSplicer/` |
| 8 | Clones AlphaFold3 upstream repo into `alphafold3/alphafold3/` |
| 9 | Downloads optional large FTP datasets (UniRef90, UniProt idmapping) |
| 10 | Checks for SpliceAI and Nextflow on PATH |
| 11 | Prints manual install checklist |
| 12 | Summary checks; calls `build_db.sh` if present and executable |

**Exclude flags** (combine with any subcommand or full install):

```bash
./bootstrap.sh \
  --exclude-pip-install \
  --exclude-evmutation \
  --exclude-build-plmc \
  --exclude-cg-cotrans \
  --exclude-netsurfp3 \
  --exclude-signalp \
  --exclude-miranda \
  --exclude-genesplicer \
  --exclude-clone-af3 \
  --exclude-uniref90 \
  --exclude-idmapping \
  --exclude-build-db
```

Contradictory combinations (e.g., `pip-only --exclude-pip-install`) produce an error.

**Manual installs** (cannot be automated -- require academic licenses):

- `netNglyc/` -- NetNGlyc 1.0 (DTU). After installing, patch the tcsh SIGNALP
  path to point to `netNglyc/bin/signalp6_adapter`.
- `netphos/` -- NetPhos 3.1 + APE (DTU), requires tcsh.
- `netMHC/` -- NetMHCpan 4.1 (DTU).
- `alphafold3/` -- NVIDIA GPU stack + Docker + AF3 model weights (Google DeepMind).

---

### `build_db.sh`

Downloads and builds the reference database used by the codon MSA, EVmutation,
miranda, and AlphaFold3 pipelines. Output root defaults to `../Bio_DBs/`
(sibling of the repo root); override with `DB_ROOT=<path>`.

```bash
./build_db.sh
# or
DB_ROOT=/data/Bio_DBs ./build_db.sh
```

**What it produces:**

```
Bio_DBs/
  assembly_summary_refseq.txt          -- NCBI RefSeq assembly index
  ftp_paths.txt                        -- filtered FTP paths for download
  refseq_assemblies/<ASM>/             -- per-assembly protein FAA + CDS FNA + feature table
  refseq_proteins_merged.faa           -- merged protein FASTA (mmseqs2 target DB source)
  idmapping.dat.gz                     -- UniProt idmapping
  protein_id_to_refseq.tsv            -- UniProt/UniParc -> RefSeq lookup table
  uniref90.fasta.gz                    -- UniRef90 (jackhmmer MSA pipeline)
  mature_hsa.fasta                     -- human mature miRNA sequences (miRBase)
  AF3/RBP_db/                          -- POSTAR3 RBP binding sites + AlphaFold MSAs
```

**Environment variable overrides:**

| Variable | Default | Description |
|----------|---------|-------------|
| `DB_ROOT` | `../Bio_DBs` | Output directory |
| `TAXON_GROUP` | `vertebrate_mammalian` | Primary RefSeq taxon group |
| `EXTRA_TAXON_GROUPS` | `vertebrate_other invertebrate` | Additional groups |
| `PARALLEL_JOBS` | `8` | Concurrent assembly downloads |
| `ARIA_SPLIT` | `4` | aria2c split count |
| `ARIA_CONN` | `4` | aria2c connections per server |
| `UNIREF90_URL` | UniProt FTP | Override UniRef90 download URL |
| `MIRBASE_HSA_URL` | `http://www.benoslab.pitt.edu/comir/hsa_data/mature_hsa.fa` | URL for human mature miRNA FASTA |
| `MIRBASE_REFRESH` | `1` | Always re-fetch mature miRNA FASTA |
| `POSTAR3_TXT_URL` | _(empty)_ | URL to human-POSTAR3.txt if needed |
| `AF3_RBP_MSA_ARCHIVE_URL` | _(empty)_ | URL to pre-built RBP MSA archive |
| `AF3_DOWNLOAD_RBP_MSAS` | `1` | Download per-RBP AF MSAs from rbp\_uniprot\_ids.txt |
| `AF3_MSA_VERSION` | `v6` | AF MSA version string |
| `AF3_MSA_URL_TEMPLATE` | `https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-msa_{VERSION}.a3m` | URL template for per-RBP AlphaFold MSA downloads |

Requires `aria2c` (preferred) or `curl`, plus `bgzip`/`tabix` for POSTAR3 indexing.

---

### `build_refseq_ftp_paths.awk`

AWK helper called internally by `build_db.sh`. Filters the NCBI RefSeq assembly
summary by taxon group, assembly level, and version status, and emits one FTP
path per passing assembly.

Not intended to be called directly, but can be used standalone:

```bash
awk -F '\t' \
  -v taxon_group="vertebrate_mammalian" \
  -f build_refseq_ftp_paths.awk \
  assembly_summary_refseq.txt > ftp_paths.txt
```

Accepts a `|`-separated list of groups:

```bash
awk -F '\t' \
  -v taxon_group="vertebrate_mammalian|vertebrate_other" \
  -f build_refseq_ftp_paths.awk \
  assembly_summary_refseq.txt > ftp_paths.txt
```

---

## Typical Setup Order

```bash
# 1. Install dependencies and third-party tools
./scripts/bootstrap.sh --exclude-uniref90 --exclude-idmapping

# 2. Build reference databases (large download -- may take hours)
./scripts/build_db.sh
```