# NetSurfP-3.0 Pipeline

Per-residue surface accessibility, secondary structure, disorder, and backbone angles for WT and mutant protein sequences, compared at the mutation site and globally.

## Requirements

| Component | Notes |
|-----------|-------|
| `nsp3` | Cloned under `NetSurfP3/nsp3/`. Requires PyTorch and ESM-1b weights. |
| Trained checkpoint | `-M/--model`; see below. |
| Config YAML | `-c/--config`; must match the checkpoint architecture. |
| Python >= 3.9 | Uses `biofeaturefactory/lib/`. |
| Input FASTAs | Nucleotide or amino acid. |
| Mutation CSVs | NT or AA tokens, matching `--input-type`. |

### Training a checkpoint

```bash
cd NetSurfP3/nsp3/nsp3
nsp3 train -c experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml
```

Checkpoints land in the config's `save_dir` (e.g. `saved/nsp3/CNNbLSTM/`).

### Choosing `--config`

This pipeline needs all six output heads (SS8, SS3, disorder, RSA, phi, psi). Use
`CNNbLSTM/CNNbLSTM.yml` unless you trained with something else.

| Config | `arch.type` | Output heads | Usable here |
|--------|-------------|--------------|-------------|
| `CNNbLSTM/CNNbLSTM.yml` | `CNNbLSTM_ESM1b_Complete` | SS8, SS3, disorder, RSA, phi, psi | Yes (recommended) |
| `CNNbLSTM_ESM1b_v2/config.yml` | `CNNbLSTM_Extended` | SS8, disorder, RSA iso/cpx, phi, psi | Partial (no explicit SS3 head) |
| `CNNbLSTM_ESM1b/CNNbLSTM_ESM1b.yml` | `CNNbLSTM_ESM1b` | SS8, SS3 | No |
| `ESM1b/ESM1b.yml` | `ESM1b` | SS8, SS3 | No |
| `CNNTrans/CNNTrans.yml` | `CNNTrans` | varies | Depends on training |

To read a checkpoint's architecture:

```python
import torch
ckpt = torch.load("checkpoint.pth", map_location="cpu", weights_only=False)
print(ckpt["config"]["arch"]["type"])
```

## Usage

```bash
# Directory mode: variant_mapping output root
python netsurfp3_pipeline.py -i out/ -o results/ \
    -M /path/to/checkpoint.pth \
    -c NetSurfP3/nsp3/experiments/netsurfp_3/CNNbLSTM/CNNbLSTM.yml

# Flat FASTA directory + explicit mutation directory
python netsurfp3_pipeline.py -i FASTA_files/nt/ -o results/ \
    -m mutations/aa/ -it aa \
    -M /path/to/checkpoint.pth -c .../CNNbLSTM.yml -l validation.log
```

In directory mode `-i` is the `variant_mapping` output root (`<root>/<GENE>/fastas/` and
`<root>/<GENE>/mappings/`); the gene is taken from the directory name. `input` and `output`
are also accepted positionally.

## Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | -- | variant_mapping output root, flat FASTA directory, or single FASTA |
| `-o, --output` | -- | Output base directory |
| `-M, --model` | required | Trained nsp3 checkpoint |
| `-c, --config` | required | nsp3 config YAML matching the checkpoint |
| `-m, --mutation-dir` | required in practice | Mutation file or directory |
| `-it, --input-type` | auto | `nt` or `aa`. Omitted: auto-detected from the WT sequence via `detect_alphabet` |
| `-l, --log` | -- | Validation log; skips failed mutations |
| `-bs, --batch-size` | `100` | Sequences per NSP3 batch |
| `--max-seq-length` | `1500` | Length above which sequences are chunked |
| `-v, --verbose` | off | Verbose output |

With `-it nt` the mutation CSVs hold NT tokens (`A1002T`, and non-SNV forms such as `ACAA1002A`,
`T28TGGT`); with `-it aa` they hold AA tokens (`M334V`, `KE100K`). Non-SNV tokens are processed by default.

## Output

```
{output}/{GENE}/NetSurfP3/
    {GENE}.netsurfp3.summary.tsv   -- per-mutation summary
    {GENE}.netsurfp3.residues.tsv  -- per-residue predictions
    {GENE}.netsurfp3.local.tsv     -- +/-5 residue window around each mutation
```

### `{GENE}.netsurfp3.summary.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `mutation_pos` | Residue index of the mutation |
| `wt_aa`, `mut_aa` | WT and mutant residues |
| `delta_rsa` | Change in relative surface accessibility |
| `delta_disorder_pf`, `delta_disorder_pt` | Change in disorder probability (two scores) |
| `ss3_change`, `ss8_change` | 1 if the secondary structure assignment changed |
| `burial_classification` | -2 to +2 (see below) |
| `disorder_classification` | -2 to +2 (see below) |
| `local_structural_impact` | Sum of absolute deltas within +/-5 residues |
| `global_mean_delta_rsa` | Mean RSA change across all residues |
| `global_ss_changes` | Count of secondary-structure changes |
| `qc_flags` | Quality control flags |

### `{GENE}.netsurfp3.residues.tsv`

| Column | Description |
|--------|-------------|
| `pos`, `aa` | Residue index and single-letter code |
| `rsa` | Relative solvent accessibility (0 buried - 1 exposed) |
| `asa` | Absolute solvent accessibility (angstrom^2) |
| `ss3` | H / E / C |
| `ss8` | G / H / I / B / E / S / T / C |
| `disorder_pf`, `disorder_pt` | Disorder probabilities |
| `phi`, `psi` | Backbone dihedral angles (degrees) |
| `sequence_type` | wt or mut |

### `{GENE}.netsurfp3.local.tsv`

| Column | Description |
|--------|-------------|
| `pkey` | `GENE-mutation` |
| `relative_pos` | Offset from the mutation (-5 to +5) |
| `absolute_pos` | Residue index in the protein |
| `delta_rsa` | RSA change at this position |
| `delta_disorder_pf`, `delta_disorder_pt` | Disorder change at this position |
| `wt_ss3`, `mut_ss3`, `wt_ss8`, `mut_ss8` | Secondary structure assignments |
| `delta_phi`, `delta_psi` | Backbone angle changes (degrees) |
| `is_mutation_site` | 0/1 |

### Classification encodings

Both `burial_classification` and `disorder_classification` are `tier(mut) - tier(wt)`, on a
three-tier scale, giving values -2 to +2.

- Burial tiers from RSA: buried (< 0.25), intermediate (0.25 <= RSA < 0.50), exposed (>= 0.50).
- Disorder tiers: ordered (< 0.3), intermediate (0.3 <= d < 0.7), disordered (>= 0.7).

## Troubleshooting

| Symptom | Resolution |
|---------|------------|
| `NetSurfP not found` | Confirm `nsp3/` is cloned and PyTorch / ESM weights are present |
| Model or config errors | `-c` must match the checkpoint's `arch.type` (see the table above) |
| Execution timeout | Lower `--max-seq-length` or `-bs` |
| `No mapping file found` | `-m` must contain `{GENE}*.csv` |
| Out of memory | Lower `-bs` |

## License

This wrapper is AGPL-3.0 - see [LICENSE](../../LICENSE) in the repository root.
The NetSurfP-3.0 library (`nsp3/`) is third-party software under its own license.
