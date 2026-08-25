# scripts/

Setup and database build scripts. Run these once before using the pipelines.

## `bootstrap.sh`

Installs pipeline dependencies into the current Python environment and clones/builds
third-party tools.

```bash
cd scripts/
./bootstrap.sh                            # everything
./bootstrap.sh env-phase                  # python environment only
./bootstrap.sh env-phase git-phase        # environment + clones/builds
./bootstrap.sh git-phase db-phase         # clones/builds + datasets
./bootstrap.sh db-phase --exclude-uniref90
```

### Phases

Phases are **additive** -- naming two runs their union, and an overlapping step runs once.
They are deliberately not called `-only`. `env`, `git`, and `db` are accepted as bare aliases.

| Phase | Scope |
|-------|-------|
| `env-phase` | pip install, conda installs, Nextflow/OpenJDK, and editable installs of repos already cloned (steps 2, 6-6d, 7b, 8b-8c). No clones, no source builds, no downloads |
| `git-phase` | Clones, source builds, conda installs, editable installs (steps 1b, 3-8c, 10-12) |
| `db-phase` | FTP downloads and `build_db.sh` (steps 9, 9b, 12) |
| _(none)_ | Every phase |

### Steps

| Step | Action |
|------|--------|
| 1 | Core tool checks (`git`, `tar`) |
| 1b | Build toolchain: gcc, g++, make |
| 1c | tcsh -- the DTU tool launchers need it |
| 1d | Perl `Array::Base` -- netNglyc's neural-net driver |
| 2 | `pip install -e ".[all]"` |
| 3 | EVmutation + plmc (clone and compile), adabmDCApy, cg_cotrans |
| 4 | NetSurfP-3.0 clone into `NetSurfP3/nsp3/` |
| 5 | SignalP 6.0 presence check |
| 6 | miranda via conda |
| 6b | mmseqs2 -- codon MSA |
| 6c | HMMER / jackhmmer -- protein MSA |
| 6c2 | HTSlib -- bgzip + tabix |
| 6d | SpliceAI, into its own conda env (default `bff-spliceai`) |
| 7 | GeneSplicer, built from JHU source |
| 7b | Nextflow + OpenJDK |
| 8 | AlphaFold3 clone into `alphafold3/alphafold3/` |
| 8b | `pip install -e` of cloned python repos (nsp3, adabmDCApy) |
| 8c | Verify the installed package layout |
| 9 | FTP datasets (UniRef90, UniProt idmapping) |
| 9b | No-op -- the CoCoPUTs table moved to `build_db.sh` |
| 10 | Nextflow check |
| 11 | Licensed/manual dependency checklist |
| 12 | Summary probes; calls `build_db.sh` if present and executable |

SpliceAI gets a dedicated conda env because sharing one reintroduces the
pyarrow/TensorFlow Abseil deadlock. It is therefore not on the active env's PATH.

### Flags

| Flag | Effect |
|------|--------|
| `--exclude-pip-install` | Skip pip install |
| `--exclude-build-tools` | Skip build-essential (gcc, g++, make) |
| `--exclude-tcsh` | Skip tcsh |
| `--exclude-perl-modules` | Skip `Array::Base` |
| `--exclude-evmutation` | Skip cloning EVmutation and plmc |
| `--exclude-build-plmc` | Skip compiling plmc; the clone still happens |
| `--exclude-adabmdca` | Skip cloning and installing adabmDCApy |
| `--exclude-cg-cotrans` | Skip downloading cg_cotrans |
| `--exclude-netsurfp3` | Skip cloning NetSurfP-3.0 |
| `--exclude-signalp` | Skip the SignalP 6.0 presence check |
| `--exclude-miranda` | Skip conda install of miranda |
| `--exclude-spliceai` | Skip conda install of spliceai |
| `--exclude-spliceai-env` | Install spliceai into the ACTIVE env instead of its own. Reintroduces the pyarrow/TensorFlow deadlock |
| `--spliceai-env NAME` | Name of the dedicated spliceai env (default `bff-spliceai`) |
| `--exclude-mmseqs2` | Skip conda install of mmseqs2 |
| `--exclude-hmmer` | Skip conda install of HMMER |
| `--exclude-htslib` | Skip conda install of HTSlib |
| `--exclude-nextflow` | Skip Nextflow + OpenJDK |
| `--exclude-editable-repos` | Skip `pip install -e` of cloned python repos |
| `--exclude-genesplicer` | Skip building GeneSplicer |
| `--exclude-clone-af3` | Skip cloning AlphaFold3 |
| `--exclude-uniref90` | Skip the UniRef90 download |
| `--exclude-idmapping` | Skip the UniProt idmapping download |
| `--exclude-cocoputs` | Skip the CoCoPUTs codon-usage table; forwarded to `build_db.sh` as `SKIP_COCOPUTS` |
| `--exclude-build-db` | Skip calling `build_db.sh` |
| `--bio-dbs DIR` | Prepared-database root. Default: nearest `Bio_DBs` above the repo, else `BFF_BIO_DBS`, else `scripts/_downloads` |
| `--fix-python` | Let conda move the env's interpreter into the supported range. **Destructive** -- packages installed under the current interpreter are orphaned |

Contradictory combinations produce an error.

### Manual installs

Academic licenses; cannot be automated.

- **NetNGlyc 1.0** (DTU), requires tcsh. Patch its tcsh `SIGNALP` path to
  `biofeaturefactory/netNglyc/bin/signalp6_adapter`.
- **SignalP 6.0** (DTU), the "fast" package. Install in a **separate** env -- SignalP 6.0 pins
  `torch<2` while nsp3 and AlphaFold3 need modern torch. netNglyc invokes it as a subprocess, so
  a separate env is fine:
  ```bash
  conda create -n signalp6 python=3.10 && conda activate signalp6
  pip install signalp-6-package/
  SIGNALP_DIR=$(python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))")
  cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
  ```
  `signalp6` and `signalp6_fast` are both recognised env names; the executable is `signalp6`
  either way. `netnglyc_pipeline.py` resolves
  `$(conda info --base)/envs/{signalp6,signalp6_fast}/bin/signalp6` automatically, so the env
  does not need to be on PATH. For other layouts pass `-snp/--signalp6-bin`, or set
  `SIGNALP6_PATH` / `SIGNALP6_HOME`.
- **NetPhos 3.1 + APE** (DTU), requires tcsh.
- **netMHC 4.0** (DTU), requires tcsh. Also download
  [`data.tar.gz`](https://services.healthtech.dtu.dk/services/NetMHC-4.0/data.tar.gz) and extract
  it into the install directory -- the distribution ships `<platform>/data` as a symlink to a
  directory that does not exist until you do. Without it netMHC exits 0 and emits no predictions.
- **cg_cotrans** -- [Shakhnovich Lab](https://shakhnovich.faculty.chemistry.harvard.edu/software/coarse-grained-co-translational-folding-analysis),
  extract into `biofeaturefactory/rare_codon/cg_cotrans/`.
- **AlphaFold3** -- NVIDIA GPU stack, Docker, and AF3 model weights.

## `build_db.sh`

Builds the reference databases used by the codon MSA, mutation effects, miranda, rare_codon, and
AlphaFold3 pipelines.

```bash
./build_db.sh
DB_ROOT=/data/Bio_DBs ./build_db.sh
```

`DB_ROOT` defaults to the nearest `Bio_DBs` directory found by walking up to five levels above
`scripts/`, falling back to `<project_root>/Bio_DBs`.

### Output

```
Bio_DBs/
  assembly_summary_refseq.txt              -- NCBI RefSeq assembly index
  ftp_paths.txt                            -- filtered FTP paths
  refseq_assemblies/<ASM>/                 -- per-assembly protein FAA + CDS FNA + feature table
  refseq_proteins_merged.faa               -- merged protein FASTA (mmseqs2 target DB source)
  idmapping.dat.gz                         -- UniProt idmapping
  protein_id_to_refseq.tsv                 -- UniProt/UniParc -> RefSeq lookup
  uniref90.fasta.gz                        -- UniRef90 (jackhmmer protein MSA)
  mature_hsa.fasta                         -- human mature miRNAs (miRBase)
  AF3/RBP_db/                              -- POSTAR3 binding sites + AlphaFold MSAs
  cocoputs/human_GRCh38_codon_usage.tsv    -- genome-wide codon usage (rare_codon null model)
```

The CoCoPUTs table is built at step 9b by `build_cocoputs_cut.py` (~11 MB of upstream zips, a few
seconds, idempotent). It matters: without it `rare_codon_pipeline.py` derives "rare" from the
single gene under analysis, making the reference distribution and the test sequence the same data.
Measured on SMN2, that fallback produced 7 false rare codons and missed 2 real ones out of a true
set of 4 -- it called GCC rare at 9.5% when GCC is 39.9% of all human alanine codons. Pass it with
`rare_codon_pipeline.py -rcu $DB_ROOT/cocoputs/human_GRCh38_codon_usage.tsv`.

### Environment overrides

| Variable | Default | Description |
|----------|---------|-------------|
| `DB_ROOT` | nearest `Bio_DBs`, else `<project_root>/Bio_DBs` | Output directory |
| `TAXON_GROUP` | `vertebrate_mammalian` | Primary RefSeq taxon group |
| `EXTRA_TAXON_GROUPS` | `vertebrate_other invertebrate` | Additional groups |
| `PARALLEL_JOBS` | `8` | Concurrent assembly downloads |
| `ARIA_SPLIT` | `4` | aria2c split count |
| `ARIA_CONN` | `4` | aria2c connections per server |
| `UNIREF90_URL` | UniProt FTP | UniRef90 download URL |
| `MIRBASE_HSA_URL` | benoslab mirror | Human mature miRNA FASTA |
| `MIRBASE_REFRESH` | `1` | Always re-fetch the miRNA FASTA |
| `POSTAR3_TXT_URL` | _(empty)_ | URL to `human-POSTAR3.txt` if needed |
| `AF3_RBP_MSA_ARCHIVE_URL` | _(empty)_ | Pre-built RBP MSA archive |
| `AF3_DOWNLOAD_RBP_MSAS` | `1` | Download per-RBP AF MSAs from `rbp_uniprot_ids.txt` |
| `AF3_MSA_VERSION` | `v6` | AF MSA version string |
| `AF3_MSA_URL_TEMPLATE` | `https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-msa_{VERSION}.a3m` | Per-RBP MSA URL template |
| `SKIP_COCOPUTS` | `0` | Skip the codon-usage table |

Requires `aria2c` (preferred) or `curl`, plus `bgzip`/`tabix` for POSTAR3 indexing.

## `build_cocoputs_cut.py`

Builds `cocoputs/human_GRCh38_codon_usage.tsv` from upstream CoCoPUTs data. Called by
`build_db.sh` step 9b; runnable standalone with `--out-dir`.

## `build_refseq_ftp_paths.awk`

AWK helper called by `build_db.sh`. Filters the RefSeq assembly summary by taxon group, assembly
level, and version status, emitting one FTP path per passing assembly.

```bash
awk -F '\t' -v taxon_group="vertebrate_mammalian|vertebrate_other" \
  -f build_refseq_ftp_paths.awk assembly_summary_refseq.txt > ftp_paths.txt
```

## Typical setup order

```bash
./scripts/bootstrap.sh --exclude-uniref90 --exclude-idmapping
./scripts/build_db.sh
```
