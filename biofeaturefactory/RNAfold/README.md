# RNAfold Variant-Centered Pipeline

## Overview
This pipeline quantifies the effect of **single-nucleotide variants (SNVs)** on **local RNA secondary structure** using the **ViennaRNA Python API**.  
Each analysis window (default = 151 nt) is centered on the variant to capture context-dependent folding changes.  
Both reference and alternate sequences are folded, sampled, and compared at multiple ensemble levels.

---

## Usage

```bash
python run_viennaRNA_pipeline.py \
    -i /path/to/fastas/ \
    -o results/ \
    --transcript-mapping /path/to/mappings/ \
    --log /path/to/validation.log \
    --samples 1000 \
    --tau 0.05 \
    --window 151 \
    --workers 8
# Writes per gene: results/{GENE}/RNAfold/{GENE}.rnafold.{tsv,positions.tsv}
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-i, --input` | Yes | -- | Input FASTA sequence file or directory |
| `-o, --output` | Yes | -- | Output base directory |
| `-w, --window` | No | `151` | Window size in nt (odd; truncates near ends) |
| `--transcript-mapping` | No | -- | Path to transcript mapping file or directory |
| `--log` | No | -- | Validation log (file or directory) used to filter failed mutations |
| `--samples` | No | `1000` | Number of Boltzmann samples per sequence |
| `--tau` | No | `0.05` | Threshold for `change_flag` on delta\_u |
| `--workers` | No | auto | Max parallel workers (processes) |

---

## Workflow Summary

1. **Window Extraction**  
   A sequence window (odd length, default = 151 nt) is extracted around the SNV using the reference transcript and corresponding alternate base.  
   Transcript positions are obtained from `get_mutation_data_bioAccurate` for biological alignment.

2. **Structure Computation**  
   - Builds a `fold_compound` object using ViennaRNA thermodynamic parameters.  
   - Computes the **Minimum Free Energy (MFE)** structure.  
   - Rescales Boltzmann parameters and computes the **partition function** to derive ensemble free energy.

3. **Boltzmann Sampling**  
   - 1,000 random structures (default) are generated via `fc.pbacktrack(samples)` from the Boltzmann ensemble.   
   - Structure energies and per-base pairing probabilities are recorded.  
   - Sampling provides an estimate of ensemble variability beyond a single deterministic structure.

4. **Variant Comparison**  
   - Ensemble and MFE energies compared between reference and alternate windows.  
   - **$\Delta\Delta G$ metrics** quantify energetic shifts.  
   - **Jensen-Shannon divergence (JSD)** captures overall structural rearrangements.  
   - **Per-position $\Delta u$** measures the change in accessibility at every base in the window.

5. **Parallel Execution**  
   - Uses `ProcessPoolExecutor` with automatic worker detection (`_autodetect_workers`).  
   - Each process runs with `OMP_NUM_THREADS=1` for isolation and predictable performance.

---

## Output Tables

### Output Structure

Output is written per gene to:

```
{output}/
  {GENE}/
    RNAfold/
      {GENE}.rnafold.tsv           -- per-mutation summary
      {GENE}.rnafold.positions.tsv -- per-position delta_u table
```

### 1. Summary Table (`{GENE}.rnafold.tsv`)
Each row represents a single SNV-centered comparison (`pkey` = gene + mutation).

| Column | Description                                                                                                                                         | Units |
|--------|-----------------------------------------------------------------------------------------------------------------------------------------------------|-------|
| `pkey` | Unique variant identifier in the form `GENE-mutation`                                                                                               | -- |
| `transcript_pos` | Transcript coordinate of the mutation (1-based)                                                                                                     | nt |
| `ddg_mfe_kcalmol` | Difference in minimum free energy (Alt - Ref). Indicates change in local thermodynamic stability of the MFE structure.                              | kcal/mol |
| `ddg_ensemble_kcalmol` | Difference in ensemble free energy (Alt - Ref). Reflects ensemble-averaged equilibrium stability over all structures.                      | kcal/mol |
| `d_meanE_kcalmol` | Difference in mean sampled energy from the Boltzmann ensemble. Represents average energetic shift across sampled conformations.                     | kcal/mol |
| `ref_sdE_kcalmol`, `alt_sdE_kcalmol` | Standard deviations of sampled energies for reference and alternate. Higher SD = more structural diversity (less rigidity).                         | kcal/mol |
| `jsd_unpaired_bits` | Jensen-Shannon divergence between per-base unpaired probability distributions. Quantifies structural dissimilarity (0 = identical, >0 = divergent). | bits |
| `delta_central` | Change in unpaired probability ($\Delta u$) at the SNV-centered position (midpoint of the window). Useful for pinpointing direct local effects.             | unitless (0-1) |

---

### 2. Per-Position Table (`{GENE}.rnafold.positions.tsv`)
Each row represents a single nucleotide position within the analyzed window.

| Column | Description                                                                                                                                | Values |
|--------|--------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `pkey` | Variant key corresponding to the parent SNV                                                                                                | -- |
| `transcript_pos` | Transcript coordinate of the variant (matches summary table)                                                                               | nt |
| `pos` | Position index within the window (1 = start, window_length = end). Position 76 corresponds to the SNV.                                     | integer |
| `delta_u` | Per-position change in unpaired probability (Alt - Ref). Positive values indicate increased accessibility in the alternate structure.      | float |
| `change_flag` | 1 if $\lvert\Delta u\rvert \geq \tau$ ($\tau = 0.05$ default) else 0. Marks positions with meaningful shifts.                              | 0/1 |
| `direction` | Sign of $\Delta u$: 1 = increased unpaired probability, -1 = decreased, 0 = unchanged. Indicates direction of accessibility change.        | -1/0/1 |
| `mfe_change_flag` | 1 if base-pairing status differs between Ref/Alt MFE structures (paired $\leftrightarrow$ unpaired).                                       | 0/1 |
| `mfe_change_dir` | Encodes the type of base-pairing change: 0 = unpaired $\rightarrow$ paired, 1 = paired $\rightarrow$ unpaired. Indicates structural opening or closing at that base. | 0/1 |

---

## Interpretation

- **$\Delta\Delta G$ (MFE / Ensemble)**
  - Positive $\Delta\Delta G$ $\rightarrow$ Alt structure is less stable (higher energy).
  - Negative $\Delta\Delta G$ $\rightarrow$ Alt structure is more stable (lower energy).
  - Ensemble $\Delta\Delta G$ captures thermodynamic stability averaged across all possible folds, while MFE $\Delta\Delta G$ isolates the most stable structure.

- **JSD_unpaired_bits**  
  - Measures how much the overall accessibility landscape changes.  
  - Values < 0.02 indicate near-identical folds; 0.05-0.1 indicate local rearrangements; >0.1 suggests significant structural change.

- **$\Delta u$ (Per-base accessibility shift)**
  - Indicates how much each nucleotide's unpaired probability changes.
  - Large $\Delta u$ around the variant may signify disruption of local base pairing.
  - Nonlocal $\Delta u$ patterns (far from the SNV) suggest propagated conformational effects.

- **change_flag and direction**  
  - Provide a thresholded binary signal for modeling or visualization.  
  - Useful for identifying structurally "sensitive" regions across variants.

- **mfe_change_flag and mfe_change_dir**  
  - Capture discrete MFE pairing state transitions, complementing $\Delta u$'s continuous probabilities.
  - Allow classification of "opening" vs "closing" effects in the predicted secondary structure.

---
## JSD Function Reference

The `jsd_unpaired` function in this pipeline measures the **difference between reference and alternate unpaired probability profiles** across the window.

**Computation implemented in code:**

- At each position *i*, the unpaired probabilities are:
  - $p_i =$ unpaired probability in the reference sequence  
  - $q_i =$ unpaired probability in the alternate sequence  
  - $m_i = (p_i + q_i)/2$

The function computes the Jensen-Shannon divergence (JSD) using $log_2$:

 $JSD(P \parallel Q) = \frac{1}{L} \sum_{i=1}^{L} \left[ H(m_i) - \frac{1}{2}(H(p_i) + H(q_i)) \right]$

where, $H(x) = -x\log_2(x) - (1-x)\log_2(1-x)$ is the Bernoulli entropy for the unpaired probability.

- **Averaging:**  
  The script takes a simple arithmetic mean across all positions (no weighting).  

- **Stability:**  
  A small epsilon ($10^{-12}$) is added to probabilities to avoid $\log(0)$.

- **Units:**  
  Because the logarithm is base 2, JSD values are in **bits**.  
  Range: $0 \le JSD \le 1$.  
  Higher values indicate larger structural differences between the two ensembles.

**Interpretation:**
- Low JSD (< 0.02): Nearly identical ensembles.  
- Moderate JSD (0.05-0.10): Localized folding rearrangements.  
- High JSD (> 0.10): Significant global restructuring of the window's folding ensemble.

**Implementation source:**
- Derived directly from *Elements of Information Theory* (Cover & Thomas, 2006).  
- Ensemble probabilities obtained via the ViennaRNA partition function (McCaskill, 1990; Lorenz et al., 2011).  


---
## $\Delta u$ and $\tau$

The symbol $u$ denotes the *unpaired probability* of a nucleotide -- the marginal probability that position $i$ remains unpaired across all Boltzmann-weighted secondary structures:

$u_i = 1 - \sum_{j \neq i} p_{ij}$

where $p_{ij}$ is the base-pairing probability between positions $i$ and $j$.

For each position, the pipeline computes:

$\Delta u_i = u_i^{(\text{alt})} - u_i^{(\text{ref})}$

- **Positive $\Delta u$** $\rightarrow$ base becomes more accessible (region opens).
- **Negative $\Delta u$** $\rightarrow$ base becomes more paired (region closes).
- $\Delta u \approx 0$ $\rightarrow$ no local structural change.

The parameter **$\tau$ (tau)** is a *magnitude threshold* applied to $\lvert\Delta u\rvert$ to define significant changes.  
By default, $\tau = 0.05$ ($\geq 5\%$ change in unpaired probability).  
This threshold filters out small fluctuations due to ensemble sampling noise, allowing `change_flag` and `direction` to represent discrete, biologically meaningful **structural perturbations**.


---
## Notes
- Default temperature is 37 °C. All energies are in kcal/mol.  
- Boltzmann sampling is stochastic; minor variation across runs is expected.  
- `OMP_NUM_THREADS=1` prevents oversubscription during multiprocessing.
- Each gene writes to its own nested output directory: `{output}/{GENE}/RNAfold/`.

---

## License

This project is licensed under the AGPL-3.0 License - see the [LICENSE](../../LICENSE) file in the root BioFeatureFactory directory for details.
