#!/usr/bin/env python3
# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
BioFeatureFactory: adabmDCA mutation scoring (parallel backend to EVmutation/plmc).

Inference runs IN-PROCESS via the adabmDCApy Python API (`from adabmDCA.training
import ...`), NOT the `adabmDCA` console script — either Boltzmann learning
(bmDCA/eaDCA/edDCA) or pseudolikelihood maximization (pseudoDCA), PyTorch/GPU
native. This module then loads the resulting text-format params and produces
per-mutation TSVs with the same routing logic as evmutation_pipeline.py:

  protein TSV  missense (and synonymous + stop when --skip-codon is set)
  codon TSV    synonymous + stop variants only

Math per mutation (i = position, a = WT token, b = mutant token):
  ΔH_independent = h_i(b) - h_i(a)
  ΔH_epistatic   = ΔH_independent + Σ_{j≠i} [J_ij(b, s_j) - J_ij(a, s_j)]
                   s_j is the WT (focus) token at column j
  pairwise       = ΔH_epistatic - ΔH_independent

Zero-sum gauge applied before scoring:
  h_i(a)    -= mean_a h_i(a)
  J_ij(a,b) -= row_mean_b + col_mean_a - global_mean (per (i,j))

Concordance (per-position relative threshold, matches EVmutation convention):
  threshold = 0.5 * std(pairwise) at this position
  |pairwise| <= threshold or std == 0 → NEUTRAL
  sign(epistatic) == sign(independent) → CONCORDANT
  otherwise → DISCORDANT
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Codon encoding utilities (bin/codon_encoding.py)
_BIN_DIR = Path(__file__).resolve().parent / "bin"
sys.path.insert(0, str(_BIN_DIR))
try:
    from codon_encoding import CHAR_TO_CODON, CODON_ALPHABET, encode_codon_msa
    _CODON_ENCODING_AVAILABLE = True
except ImportError:
    _CODON_ENCODING_AVAILABLE = False
    CODON_ALPHABET = None
    CHAR_TO_CODON = {}

from biofeaturefactory.lib.utility import (
    codon_to_aa,
    discover_fasta_files,
    extract_gene_from_filename,
    load_validation_failures,
    mutation_class,
    read_fasta,
    should_skip_mutation,
    trim_muts,
    write_tsv,
)


# ---------------------------------------------------------------------------
# Output schemas (distinct from EVmutation)
# ---------------------------------------------------------------------------

PROTEIN_FIELDNAMES_ADABM = [
    "pkey",
    "nt_mutant",
    "codon_position",
    "wt_codon",
    "mut_codon",
    "mutant",
    "pos",
    "wt",
    "subs",
    "mutation_class",
    "prediction_protein_independent_adabm",
    "prediction_protein_epistatic_adabm",
    "protein_pairwise_contribution_adabm",
    "protein_concordance_adabm",
    "frequency_adabm",
    "qc_flags",
]

CODON_FIELDNAMES_ADABM = [
    "pkey",
    "nt_mutant",
    "codon_position",
    "wt_codon",
    "mut_codon",
    "mutation_class",
    "prediction_codon_independent_adabm",
    "prediction_codon_epistatic_adabm",
    "codon_pairwise_contribution_adabm",
    "codon_concordance_adabm",
    "codon_frequency_adabm",
    "qc_flags",
]


# ---------------------------------------------------------------------------
# adabmDCA training (in-process, direct adabmDCApy API)
# ---------------------------------------------------------------------------
# adabmDCApy is a Python package, not a compiled binary. Subprocessing into
# its CLI (`adabmDCA train`) just adds startup overhead and re-introduces the
# exit-code-swallowing bug in adabmDCA.cli (which discards the inner
# train.py's return code, masking CUDA OOMs and other crashes as silent
# successes). We call the training building blocks directly here — same
# sequence train.py:main() runs, minus the argparse + subprocess plumbing.
# Exceptions propagate naturally.

def train_adabmdca_in_process(
    msa_path: str,
    alphabet: str,
    output_params_path: str,
    model: str = "bmDCA",
    nepochs: int = 50000,
    target: float = 0.95,
    lr: float = 0.01,
    nchains: int = 10000,
    nsweeps: int = 10,
    device_str: str = "cuda",
    dtype_str: str = "float32",
    seed: int = 0,
    quiet: bool = False,
) -> str:
    """
    Train an adabmDCA model on an MSA and write the params text file.

    Returns the path to the written params file. Mirrors the call sequence
    in adabmDCApy's `scripts/train.py:main()` but invoked in-process so
    exceptions propagate cleanly and there's no subprocess startup cost.
    """
    # Lazy imports: the score-only path (load params + score mutations)
    # doesn't need adabmDCApy. Importing at function-call time means a host
    # without adabmDCApy can still run scoring against pre-built params.
    import torch
    import numpy as np
    from adabmDCA.dataset import DatasetDCA
    from adabmDCA.fasta import get_tokens
    from adabmDCA.io import save_params
    from adabmDCA.utils import init_chains, init_parameters, get_device, get_dtype
    from adabmDCA.sampling import get_sampler
    from adabmDCA.checkpoint import Checkpoint
    from adabmDCA.training import train_graph, train_eaDCA, train_edDCA

    ROUTINES = {"bmDCA": train_graph, "eaDCA": train_eaDCA, "edDCA": train_edDCA}
    if model not in ROUTINES:
        raise ValueError(f"unknown adabmDCA model {model!r}; expected one of {list(ROUTINES)}")

    device = get_device(device_str)
    dtype = get_dtype(dtype_str)

    if not quiet:
        alpha_label = "codon (65-char)" if len(alphabet) == 65 else alphabet
        print(f"  Training adabmDCA in-process: model={model}, alphabet={alpha_label}")
        print(f"    device={device_str}, dtype={dtype_str}, nepochs={nepochs}, "
              f"target={target}, nchains={nchains}, lr={lr}")

    # 1. Dataset (loads MSA, filters/dedups, computes sequence weights).
    dataset = DatasetDCA(
        path_data=msa_path,
        path_weights=None,
        alphabet=alphabet,
        clustering_th=0.8,
        no_reweighting=False,
        device=device,
        dtype=dtype,
        message=False,
        filter_sequences=True,
        remove_duplicates=True,
    )
    tokens = get_tokens(alphabet)
    torch.manual_seed(seed)
    dataset.shuffle()

    L = dataset.get_num_residues()
    q = dataset.get_num_states()
    M_eff = dataset.get_effective_size()
    pseudocount = 1.0 / M_eff
    if not quiet:
        print(f"    L={L}, q={q}, M_eff={M_eff:.1f}, pseudocount={pseudocount:.6f}")

    # 2. Target frequencies + initial params + mask.
    fi_target, fij_target = dataset.get_frequencies(pseudocount=pseudocount)
    params = init_parameters(fi=fi_target)
    if model in ("bmDCA", "edDCA"):
        mask = torch.ones(size=(L, q, L, q), dtype=torch.bool, device=device)
        mask[torch.arange(L), :, torch.arange(L), :] = 0
    else:  # eaDCA
        mask = torch.zeros(size=(L, q, L, q), device=device, dtype=torch.bool)

    # 3. PCD chains + sampler.
    chains = init_chains(num_chains=nchains, L=L, q=q, fi=fi_target, device=device, dtype=dtype)
    log_weights = torch.zeros(size=(nchains,), device=device, dtype=dtype)
    sampler = torch.jit.script(get_sampler("gibbs"))

    # 4. Checkpoint object — periodic + final save. file_paths["params"] is
    # the canonical write target for the trained model.
    out_dir = os.path.dirname(output_params_path) or "."
    base = os.path.basename(output_params_path)
    label_stem = base[:-len("_params.dat")] if base.endswith("_params.dat") \
                 else os.path.splitext(base)[0]
    file_paths = {
        "log":    os.path.join(out_dir, f"{label_stem}.log"),
        "params": output_params_path,
        "chains": os.path.join(out_dir, f"{label_stem}_chains.fasta"),
    }
    # Every key Checkpoint.__init__ reads must be present (verified against
    # adabmDCA/checkpoint.py at versions 0.7.x): nepochs, label, model, data,
    # alphabet, sampler, nchains, nsweeps, lr, pseudocount, dtype, target,
    # seed, and (eaDCA only) gsteps + factivate.
    args_dict = {
        "label":       label_stem,
        "model":       model,
        "data":        msa_path,
        "alphabet":    alphabet,
        "sampler":     "gibbs",
        "lr":          lr,
        "nsweeps":     nsweeps,
        "nchains":     nchains,
        "nepochs":     nepochs,
        "target":      target,
        "device":      device_str,
        "dtype":       dtype_str,
        "seed":        seed,
        "pseudocount": pseudocount,
        "gsteps":      10,     # eaDCA-only; harmless for bmDCA/edDCA
        "factivate":   0.001,  # eaDCA-only
    }
    checkpoint = Checkpoint(
        file_paths=file_paths,
        tokens=tokens,
        args=args_dict,
        use_wandb=False,
    )

    # 5. Training. The eaDCA/edDCA-specific kwargs (factivate, gsteps, drate,
    # target_density) are part of every training_routine's signature in
    # adabmDCApy; the bmDCA routine just ignores them.
    #
    # OOM catch: Boltzmann learning's peak memory is ~4-5 × L²q². For codon
    # at q=65 the constant is large enough that medium-sized genes can OOM
    # on a 40 GB GPU. Catch and surface the pseudolikelihood hint so the user
    # sees an actionable suggestion in Nextflow's .command.err, not just an
    # opaque RuntimeError trace from inside TorchScript.
    training_routine = ROUTINES[model]
    try:
        _, _, _, history = training_routine(
            sampler=sampler,
            fij_target=fij_target,
            fi_target=fi_target,
            fi_test=None,
            fij_test=None,
            params=params,
            mask=mask,
            chains=chains,
            log_weights=log_weights,
            target_pearson=target,
            pseudo_count=pseudocount,
            nsweeps=nsweeps,
            max_epochs=nepochs,
            lr=lr,
            factivate=0.001,
            gsteps=10,
            drate=0.01,
            target_density=0.02,
            checkpoint=checkpoint,
            l2_reg=0.0,
        )
    except RuntimeError as e:
        if "out of memory" not in str(e).lower():
            raise  # not an OOM — propagate the original error untouched
        peak_low  = 4 * L * L * q * q * 4 / 1e9
        peak_high = 5 * L * L * q * q * 4 / 1e9
        pl_low    = 2 * L * L * q * q * 4 / 1e9
        pl_high   = 3 * L * L * q * q * 4 / 1e9
        raise RuntimeError(
            f"adabmDCA Boltzmann training OOM'd at L={L}, q={q}.\n"
            f"  Peak memory for Boltzmann is ~4–5 × L²q² × 4 bytes (float32) = "
            f"~{peak_low:.1f}–{peak_high:.1f} GiB.\n"
            f"  Pseudolikelihood mode roughly halves this (~{pl_low:.1f}–{pl_high:.1f} GiB)\n"
            f"  by dropping the MCMC chains + fij_chains buffers. Re-run with:\n"
            f"      --adabmdca-model pseudoDCA\n"
            f"  For genes much larger than ~1200 codons at q=65, even pseudolikelihood\n"
            f"  won't fit on a single GPU — use the EVmutation/plmc codon backend\n"
            f"  (pseudolikelihood on CPU, scales via system RAM).\n"
            f"  Original error: {e}"
        ) from e

    if not quiet:
        pearson = history.get("Pearson", [float("nan")])[-1]
        density = history.get("Density", [float("nan")])[-1]
        epochs  = history.get("Epochs",  [-1])[-1]
        print(f"    Training done: epochs={epochs}, Pearson={pearson:.4f}, Density={density:.4f}")

    # 6. Ensure params landed on disk. Checkpoint writes at training end, but
    # save explicitly as a safety net (and so the caller can find the file).
    if not os.path.isfile(output_params_path):
        # Save mask: UPPER-TRIANGULAR i<j only, same contract as the
        # pseudolikelihood save below — save_params writes every mask-True
        # entry and load_adabmdca_params symmetrizes via `J + J.transpose`,
        # so writing the both-triangle training mask would double every
        # coupling on reload.
        save_mask = torch.zeros_like(mask)
        _iu = torch.triu_indices(L, L, offset=1, device=device)  # (2, L*(L-1)/2), i<j
        save_mask[_iu[0], :, _iu[1], :] = True
        save_params(output_params_path, params, tokens, save_mask)

    return output_params_path


# ---------------------------------------------------------------------------
# adabmDCA training via pseudolikelihood maximization
# ---------------------------------------------------------------------------
# Alternative to Boltzmann learning. Drops the MCMC sampling + fij_chains
# buffer, reducing peak GPU memory from ~4-5 × L²q² to ~2-3 × L²q² at the
# same L. The math is per-position conditional likelihood maximization
# (Ekeberg et al. 2013; the approach plmc uses on CPU). GPU-accelerated here
# so it stays fast for the medium-gene regime that doesn't fit Boltzmann.
#
# Note on scaling: the L² scaling on the dense J tensor is unchanged. This
# is a constant-factor improvement (~2×), not a scaling-class change. For
# genes much larger than ~1200 codons at q=65 the dense J still won't fit
# on a single GPU; in that regime fall back to the EVmutation/plmc codon
# backend which holds J in system RAM.

def train_adabmdca_pseudolikelihood_in_process(
    msa_path: str,
    alphabet: str,
    output_params_path: str,
    nepochs: int = 500,
    lr: float = 0.1,
    lambda_h: float = 0.01,
    lambda_J: float = 0.01,
    chunk_size: int = 32,
    device_str: str = "cuda",
    dtype_str: str = "float32",
    seed: int = 0,
    quiet: bool = False,
) -> str:
    """
    Train an adabmDCA Potts model by pseudolikelihood maximization.

    Parameters mirror Ekeberg 2013 / Hopf 2017 / plmc conventions:
      lambda_h, lambda_J: L2 regularization strengths on h and J.
      nepochs: number of gradient steps (plmc default: 500).
      lr: learning rate for plain SGD on h and J.
      chunk_size: positions processed per j-loop chunk in the J gradient
                  accumulator. Larger = faster but more peak intermediate
                  memory (4 × L × q × chunk × q × 4 bytes float32).

    Writes the trained params to output_params_path in adabmDCApy text format.
    Returns that path.
    """
    import torch
    from adabmDCA.dataset import DatasetDCA
    from adabmDCA.fasta import get_tokens
    from adabmDCA.io import save_params
    from adabmDCA.utils import init_parameters, get_device, get_dtype

    device = get_device(device_str)
    dtype = get_dtype(dtype_str)

    if not quiet:
        alpha_label = "codon (65-char)" if len(alphabet) == 65 else alphabet
        print(f"  Training adabmDCA in-process via PSEUDOLIKELIHOOD: alphabet={alpha_label}")
        print(f"    device={device_str}, dtype={dtype_str}, nepochs={nepochs}, "
              f"lr={lr}, lambda_h={lambda_h}, lambda_J={lambda_J}")

    # 1. Dataset (same as Boltzmann path)
    dataset = DatasetDCA(
        path_data=msa_path,
        path_weights=None,
        alphabet=alphabet,
        clustering_th=0.8,
        no_reweighting=False,
        device=device,
        dtype=dtype,
        message=False,
        filter_sequences=True,
        remove_duplicates=True,
    )
    tokens = get_tokens(alphabet)
    torch.manual_seed(seed)

    L = dataset.get_num_residues()
    q = dataset.get_num_states()
    M_eff = dataset.get_effective_size()
    if not quiet:
        print(f"    L={L}, q={q}, M_eff={M_eff:.1f}")

    # 2. Get the encoded sequences + weights for direct gradient computation.
    # dataset.data is shape (M, L) of token indices; dataset.weights is (M,).
    X = dataset.data.long().to(device)
    weights = dataset.weights.to(device=device, dtype=dtype)
    weights_sum = weights.sum()

    # 3. Initialize h and J. Init h from single-frequencies (same as Boltzmann),
    # J at zero (clean initialization for pseudolikelihood).
    pseudocount = 1.0 / M_eff
    fi_target, _ = dataset.get_frequencies(pseudocount=pseudocount)
    init = init_parameters(fi=fi_target)
    h = init["bias"].to(device=device, dtype=dtype)
    J = init["coupling_matrix"].to(device=device, dtype=dtype)

    # Save mask: UPPER-TRIANGULAR i<j only. save_params writes every mask-True
    # entry, and load_adabmdca_params reconstructs J via `J + J.transpose`
    # (symmetrize, see loader docstring). Writing BOTH triangles here would
    # double every coupling on reload; J is already symmetric (Ekeberg grad
    # symmetrization below), so the i<j triangle fully determines it.
    # (L_idx is also reused by the training loop below.)
    L_idx = torch.arange(L, device=device)
    mask = torch.zeros((L, q, L, q), dtype=torch.bool, device=device)
    _iu = torch.triu_indices(L, L, offset=1, device=device)  # (2, L*(L-1)/2), i<j
    mask[_iu[0], :, _iu[1], :] = True

    # 4. Pre-compute the one-hot encoding of data sequences (used in gradient).
    # Memory: B × L × q × 4 bytes. For F9 codon (B=531, L=764, q=65): ~106 MB.
    B = X.shape[0]
    X_oh = torch.nn.functional.one_hot(X, num_classes=q).to(dtype)  # (B, L, q)

    # Heads-up if the user passed a Boltzmann-sized nepochs (default 50000) into
    # pseudolikelihood — plmc's default is 500. Doesn't fail, just warns.
    if not quiet and nepochs > 5000:
        print(f"    NOTE nepochs={nepochs} is large for pseudolikelihood "
              f"(plmc default: 500). Consider --adabmdca-nepochs 500.")

    # 5. Training loop. Plain SGD on negative-PL + L2.
    #
    # MEMORY-CRITICAL CHANGES vs. the initial version:
    #   (a) torch.no_grad() wraps the loop — autograd is irrelevant for
    #       hand-computed gradients but was building a graph during the j-sum
    #       that retained ~10 GiB of per-iteration intermediates.
    #   (b) grad_J and grad_h allocated ONCE before the loop; zeroed in-place
    #       per epoch. Avoids allocator churn that the cached pool can't fully
    #       reuse for the same-sized re-allocation.
    #   (c) Parameter updates use in-place ops (mul_, add_) so we never hold
    #       two copies of h or J simultaneously.
    #   (d) E built and accumulated in-place (add_/sub_) instead of rebinding
    #       to a new tensor each j-step.
    grad_J = torch.zeros_like(J)
    grad_h = torch.zeros_like(h)
    # Pre-allocate the (B, L, q) E buffer once and reuse via copy_().
    E = torch.empty(B, L, q, device=device, dtype=dtype)
    try:
        with torch.no_grad():
            for epoch in range(nepochs):
                # ----- Forward: compute conditional fields E[s, i, a] -----
                # E[s, i, a] = h[i, a] + Σ_{j != i} J[i, a, j, X[s, j]]
                # Compute the full j-sum, then subtract the diagonal j=i term.

                # Initialize E with the broadcast h (in-place into the reused buffer).
                E.copy_(h.unsqueeze(0).expand(B, -1, -1))

                for j in range(L):
                    # J[:, :, j, X[:, j]] → (L, q, B) view → permute to (B, L, q),
                    # add to E in-place.
                    E.add_(J[:, :, j, :][:, :, X[:, j]].permute(2, 0, 1))

                # Subtract the i==i diagonal contribution we summed above.
                # diag_contrib[s, i, a] = J[i, a, i, X[s, i]]
                J_diag = J[L_idx, :, L_idx, :]                          # (L, q, q)
                diag_contrib = J_diag.gather(
                    dim=2,
                    index=X.transpose(0, 1).unsqueeze(1).expand(-1, q, -1),  # (L, q, B)
                ).permute(2, 0, 1)                                       # (B, L, q)
                E.sub_(diag_contrib)
                del diag_contrib

                # ----- Conditional log-likelihood + softmax probability -----
                logZ = torch.logsumexp(E, dim=-1, keepdim=True)         # (B, L, 1)
                log_p = E - logZ                                         # (B, L, q)
                p = torch.exp(log_p)                                     # (B, L, q)

                # Per-sequence summed log-likelihood: Σ_i log p(X[s, i] | s_{-i})
                log_p_at_data = log_p.gather(-1, X.unsqueeze(-1)).squeeze(-1)  # (B, L)
                PL = (weights.unsqueeze(-1) * log_p_at_data).sum() / weights_sum

                # ----- Gradient of PL (data − model) weighted by sequence weights -----
                diff = X_oh - p                                          # (B, L, q)
                # weights/weights_sum, not weights: PL above is normalized by
                # weights_sum, so the gradient must be too — otherwise lr and
                # the L2 term below are effectively scaled by M_eff.
                wdiff = diff * (weights / weights_sum).view(-1, 1, 1)    # weighted (B, L, q)
                # Free intermediates we no longer need so PyTorch's allocator can reuse.
                del log_p, p, log_p_at_data, diff, logZ

                # ∂PL/∂h[i, a] = Σ_s w_s · (X_oh[s, i, a] - p[s, i, a]) / weights_sum
                torch.sum(wdiff, dim=0, out=grad_h)                      # in-place into grad_h

                # ∂PL/∂J[i, a, j, b] — symmetrized as in Ekeberg 2013 eq. 9.
                # Term1[i, a, j, b] = Σ_s w_s · diff[s, i, a] · X_oh[s, j, b] / weights_sum
                # Symmetrized: grad_J_sym[i,a,j,b] = 0.5·(Term1[i,a,j,b] + Term1[j,b,i,a]).
                #
                # We build the symmetrized gradient DIRECTLY via per-chunk dual
                # accumulation so we never need to materialize the transpose of
                # the full (L, q, L, q) grad_J (which would double peak memory).
                # For each j-chunk:
                #   - 0.5 · Term1[i, a, j_chunk, b]  → grad_J[:, :, j_chunk, :]
                #   - 0.5 · Term1[j_chunk, b, i, a]  → grad_J[j_chunk, :, :, :]
                #     (achieved by permuting the chunk axes: (L, q, chunk, q) →
                #      (chunk, q, L, q) and broadcasting into the j_chunk slice
                #      of dim 0.)
                grad_J.zero_()
                for j_start in range(0, L, chunk_size):
                    j_end = min(j_start + chunk_size, L)
                    X_oh_chunk = X_oh[:, j_start:j_end, :]               # (B, chunk, q)
                    term1_chunk = torch.einsum(
                        "sia,sjb->iajb", wdiff, X_oh_chunk
                    )                                                     # (L, q, chunk, q)
                    # Contribution of 0.5 · Term1 to the (∗, ∗, j_chunk, ∗) slot:
                    grad_J[:, :, j_start:j_end, :].add_(term1_chunk, alpha=0.5)
                    # Contribution of 0.5 · Term1^T to the (j_chunk, ∗, ∗, ∗) slot:
                    grad_J[j_start:j_end, :, :, :].add_(
                        term1_chunk.permute(2, 3, 0, 1), alpha=0.5
                    )
                    del term1_chunk
                grad_J[L_idx, :, L_idx, :] = 0.0

                # ----- L2 regularization + SGD step (all in-place) -----
                # max (PL - 0.5 λ ||·||²)  ⇒  param ← param·(1 - lr·λ) + lr·grad
                h.mul_(1.0 - lr * lambda_h).add_(grad_h, alpha=lr)
                J.mul_(1.0 - lr * lambda_J).add_(grad_J, alpha=lr)
                J[L_idx, :, L_idx, :] = 0.0  # keep diagonal at zero after the update

                if not quiet and (epoch % 25 == 0 or epoch == nepochs - 1):
                    grad_norm = (grad_h.norm() + grad_J.norm()).item()
                    print(f"    epoch {epoch:4d}: PL = {PL.item():.4f}, "
                          f"||grad|| = {grad_norm:.4e}")
                del wdiff
    except RuntimeError as e:
        # Same OOM hint pattern as the Boltzmann path. If pseudolikelihood
        # itself OOMs, the user is over the ceiling for any single-GPU dense
        # Potts approach and needs to fall back to EVmutation/plmc codon.
        if "out of memory" not in str(e).lower():
            raise
        pl_low  = 2 * L * L * q * q * 4 / 1e9
        pl_high = 3 * L * L * q * q * 4 / 1e9
        raise RuntimeError(
            f"adabmDCA pseudolikelihood training OOM'd at L={L}, q={q}.\n"
            f"  Peak memory ~{pl_low:.1f}–{pl_high:.1f} GiB. GPU couldn't supply.\n"
            f"  Pseudolikelihood is already the more memory-efficient option;\n"
            f"  this gene exceeds the single-GPU dense-Potts ceiling. Use the\n"
            f"  EVmutation/plmc codon backend instead — it runs the same\n"
            f"  pseudolikelihood math on CPU with system-RAM-scale budgets.\n"
            f"  Original error: {e}"
        ) from e

    # 6. Save the trained params to disk in adabmDCApy text format so the
    # downstream scoring step can load them the same way as the Boltzmann path.
    final_params = {"bias": h, "coupling_matrix": J}
    save_params(output_params_path, final_params, tokens, mask)
    return output_params_path


# ---------------------------------------------------------------------------
# adabmDCApy text-format params loader (numpy-only, no torch dependency)
# ---------------------------------------------------------------------------

def load_adabmdca_params(path: str, alphabet: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse adabmDCApy text-format params.

    Format:
      'J idx0 idx1 aa0 aa1 value'   only upper-triangular i<j entries saved
      'h idx aa value'

    Returns:
      h: (L, q) float64 array
      J: (L, q, L, q) float64 array (symmetrized)
    """
    token_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    q = len(alphabet)

    h_records: List[Tuple[int, int, float]] = []
    j_records: List[Tuple[int, int, int, int, float]] = []

    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "h":
                # h idx aa value
                idx = int(parts[1])
                aa = parts[2]
                val = float(parts[3])
                h_records.append((idx, token_to_idx[aa], val))
            elif tag == "J":
                # J idx0 idx1 aa0 aa1 value
                idx0 = int(parts[1])
                idx1 = int(parts[2])
                aa0 = parts[3]
                aa1 = parts[4]
                val = float(parts[5])
                j_records.append((idx0, idx1, token_to_idx[aa0], token_to_idx[aa1], val))

    if not h_records:
        raise RuntimeError(f"No 'h' entries found in {path}")

    L = max(rec[0] for rec in h_records) + 1
    if j_records:
        L = max(L, max(rec[0] for rec in j_records) + 1, max(rec[1] for rec in j_records) + 1)

    h = np.zeros((L, q), dtype=np.float64)
    for idx, ai, val in h_records:
        h[idx, ai] = val

    J = np.zeros((L, L, q, q), dtype=np.float64)
    for i, j, ai, bj, val in j_records:
        J[i, j, ai, bj] = val
    # Upper-triangular only is saved; symmetrize so J[j,i,b,a] = J[i,j,a,b]
    J = J + J.transpose(1, 0, 3, 2)

    # Reorganize to (L, q, L, q) — matches EVmutation's J_ij[i, a, j, b] layout
    J = J.transpose(0, 2, 1, 3)

    return h, J


# ---------------------------------------------------------------------------
# Zero-sum gauge transformation
# ---------------------------------------------------------------------------

def apply_zero_sum_gauge(h: np.ndarray, J: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply zero-sum gauge in place semantics (returns new arrays).

    Convention (matches plmc / EVmutation / adabmDCA.dca.set_zerosum_gauge):
      h_i(a)    -= mean_a h_i(a)
      J_ij(a,b) -= row_mean_b + col_mean_a - global_mean   per (i,j)

    Resulting h satisfies sum_a h_i(a) ~= 0; coupling rows and columns sum to
    zero per (i,j) block.
    """
    h_g = h - h.mean(axis=1, keepdims=True)

    # J shape: (L, q, L, q) -> dim 1 = a (left), dim 3 = b (right)
    row_mean = J.mean(axis=1, keepdims=True)
    col_mean = J.mean(axis=3, keepdims=True)
    global_mean = J.mean(axis=(1, 3), keepdims=True)
    J_g = J - row_mean - col_mean + global_mean

    return h_g, J_g


# ---------------------------------------------------------------------------
# Position lookup builders
# ---------------------------------------------------------------------------

def _load_focus_sequence(msa_path: str, focus_id: Optional[str]) -> Tuple[str, str]:
    """
    Load the focus sequence from an MSA.

    Prefers focus_id when provided. Falls back to first record otherwise.
    Returns (focus_id_used, focus_seq_with_gaps).
    """
    seqs = read_fasta(msa_path)
    if not seqs:
        raise RuntimeError(f"No sequences found in MSA: {msa_path}")
    if focus_id:
        for sid, seq in seqs.items():
            # exact or prefix match (jackhmmer mangles names)
            if sid == focus_id or sid.startswith(focus_id + "/") or sid.split("/")[0] == focus_id:
                return sid, seq
    first_id, first_seq = next(iter(seqs.items()))
    return first_id, first_seq


def _column_frequencies(msa_path: str, alphabet: str, encoded: bool) -> np.ndarray:
    """
    Per-column non-gap frequency of each token.

    Gaps are excluded from the denominator, so each column's non-gap
    frequencies sum to 1 and the gap token itself reports 0.0. An all-gap
    column reports 0.0 for every token.

    Returns: (L, q) array where L = column count and q = alphabet size.
    For encoded codon MSAs, callers pass encoded=True (already single-char).
    For protein MSAs, encoded=False and the sequences are read as-is.
    """
    seqs = read_fasta(msa_path)
    if not seqs:
        raise RuntimeError(f"No sequences in MSA: {msa_path}")

    seq_list = list(seqs.values())
    L = len(seq_list[0])
    if any(len(s) != L for s in seq_list):
        raise RuntimeError(f"MSA sequences are not all the same length in {msa_path}")

    token_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    gap_idx = token_to_idx.get("-")
    q = len(alphabet)
    counts = np.zeros((L, q), dtype=np.float64)

    for s in seq_list:
        for col, ch in enumerate(s):
            ch_up = ch if encoded else ch.upper()
            idx = token_to_idx.get(ch_up)
            if idx is not None:
                counts[col, idx] += 1.0

    totals = counts.sum(axis=1, keepdims=True)
    if gap_idx is not None:
        totals = totals - counts[:, [gap_idx]]
        counts[:, gap_idx] = 0.0
    totals[totals == 0] = 1.0
    return counts / totals


def build_lookup(
    h: np.ndarray,
    J: np.ndarray,
    focus_seq: str,
    alphabet: str,
    column_freqs: np.ndarray,
) -> Tuple[Dict[Tuple[int, str], Dict], Dict[int, int]]:
    """
    Build the per-(seq_pos, mut_token) lookup table.

    seq_pos is the 1-based focus-sequence position (CDS codon position for the
    codon side, AA position for the protein side). Alignment columns where the
    focus has a gap are skipped — they don't correspond to a sequence position.

    Returns:
      lookup: dict keyed by (seq_pos, mut_token_char) → score dict
      col_to_seq: dict mapping alignment column → seq_pos (1-based)
    """
    L_align = h.shape[0]
    if len(focus_seq) != L_align:
        raise RuntimeError(
            f"Focus sequence length {len(focus_seq)} does not match params L {L_align}"
        )

    token_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    gap_idx = token_to_idx.get("-")

    # Map alignment column → CDS/AA position (1-based, skipping gap columns in focus)
    col_to_seq: Dict[int, int] = {}
    seq_pos = 0
    for col, ch in enumerate(focus_seq):
        idx = token_to_idx.get(ch)
        if idx is None or idx == gap_idx:
            continue
        seq_pos += 1
        col_to_seq[col] = seq_pos

    # Vectorize WT token indices over all alignment columns (gap → 0; we'll skip via col_to_seq)
    wt_idx_per_col = np.array(
        [token_to_idx.get(ch, gap_idx if gap_idx is not None else 0) for ch in focus_seq],
        dtype=np.int64,
    )

    # First pass: compute scores for every (col, b) where col is a non-gap focus column
    raw: Dict[Tuple[int, int, int], Dict] = {}
    pos_contributions: Dict[int, List[float]] = {}

    for col_i, sp in col_to_seq.items():
        wt_idx = wt_idx_per_col[col_i]
        h_i = h[col_i]                            # (q,)
        J_i = J[col_i]                            # (q, L, q): J[col_i, a, j, b]

        # Σ_{j≠col_i} J[col_i, a, j, s_j]  for each a
        # J_i[:, j, wt_idx_per_col[j]]  → shape (q, L)
        col_indices = np.arange(L_align)
        # gather J[col_i, :, j, wt[j]] over all j
        gathered = J_i[:, col_indices, wt_idx_per_col]  # (q, L)
        gathered[:, col_i] = 0.0                         # zero out j == i
        sum_over_j = gathered.sum(axis=1)               # (q,)

        # Scores for every alternative b
        for b_idx in range(h.shape[1]):
            if b_idx == wt_idx:
                continue
            indep = h_i[b_idx] - h_i[wt_idx]
            epi   = indep + (sum_over_j[b_idx] - sum_over_j[wt_idx])
            pair  = epi - indep
            raw[(sp, col_i, b_idx)] = {
                "indep": float(indep),
                "epi":   float(epi),
                "pair":  float(pair),
            }
            pos_contributions.setdefault(sp, []).append(float(pair))

    # Per-position std of pairwise contributions
    pos_std = {sp: (float(np.std(c)) if len(c) > 1 else 0.0) for sp, c in pos_contributions.items()}

    # Second pass: classify concordance and pack the lookup
    lookup: Dict[Tuple[int, str], Dict] = {}
    idx_to_token = list(alphabet)

    for (sp, col_i, b_idx), s in raw.items():
        std = pos_std.get(sp, 0.0)
        threshold = 0.5 * std
        if std == 0.0 or abs(s["pair"]) <= threshold:
            concordance = "NEUTRAL"
        elif (s["epi"] >= 0) == (s["indep"] >= 0):
            concordance = "CONCORDANT"
        else:
            concordance = "DISCORDANT"
        b_char = idx_to_token[b_idx]
        freq = float(column_freqs[col_i, b_idx])
        lookup[(sp, b_char)] = {
            "indep": s["indep"],
            "epi":   s["epi"],
            "pair":  s["pair"],
            "concordance": concordance,
            "frequency":   freq,
        }

    return lookup, col_to_seq


# ---------------------------------------------------------------------------
# Mutation scoring (mirrors evmutation_pipeline.score_nt_mutations)
# ---------------------------------------------------------------------------

def aa_symbol(aa: str) -> str:
    return "*" if aa == "Stop" else aa


def score_nt_mutations_adabm(
    nt_mutations: List[str],
    gene: str,
    orf_seq: str,
    aa_lookup: Optional[Dict[Tuple[int, str], Dict]] = None,
    codon_lookup: Optional[Dict[Tuple[int, str], Dict]] = None,
    failure_map: Optional[Dict] = None,
    skip_codon: bool = False,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Map and score NT mutations through the adabmDCA backend.

    aa_lookup keys: (aa_pos_1based, mut_aa_symbol) → {indep, epi, pair, concordance, frequency}
    codon_lookup keys: (aa_pos_1based, mut_codon_3letter) → {indep, epi, pair, concordance, frequency}

    Routing:
      missense / unknown_codon  → protein TSV (scored via aa_lookup if present)
      synonymous (skip_codon=False) → codon TSV
      synonymous (skip_codon=True)  → protein TSV with score 0 (trivial: M→M)
      stop_gain/loss (skip_codon=False) → codon TSV (annotation only)
      stop_gain/loss (skip_codon=True)  → protein TSV with empty scores
    """
    failure_map = failure_map or {}
    nt_re = re.compile(r"^([ACGT])(\d+)([ACGT])$")
    protein_rows: List[Dict] = []
    codon_rows: List[Dict] = []

    for nt_mut in nt_mutations:
        if should_skip_mutation(gene, nt_mut, failure_map):
            continue

        pkey = f"{gene}-{nt_mut}"
        qc_flags: List[str] = []

        m = nt_re.match(nt_mut)
        if not m:
            prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "INVALID_MUTATION"})
            protein_rows.append(prow)
            continue

        ref_nt, pos_str, alt_nt = m.groups()
        nt_pos = int(pos_str)
        idx = nt_pos - 1
        if idx < 0 or idx >= len(orf_seq):
            prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "OUT_OF_RANGE"})
            protein_rows.append(prow)
            continue

        if orf_seq[idx] != ref_nt:
            qc_flags.append("REF_MISMATCH")

        codon_start = (idx // 3) * 3
        if codon_start + 3 > len(orf_seq):
            prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": "PARTIAL_CODON"})
            protein_rows.append(prow)
            continue

        wt_codon = orf_seq[codon_start:codon_start + 3]
        mut_codon_list = list(wt_codon)
        mut_codon_list[idx % 3] = alt_nt
        mut_codon = "".join(mut_codon_list)

        wt_aa  = codon_to_aa.get(wt_codon, "X")
        mut_aa = codon_to_aa.get(mut_codon, "X")
        aa_pos = (codon_start // 3) + 1
        aa_mutant = f"{aa_symbol(wt_aa)}{aa_pos}{aa_symbol(mut_aa)}"
        mclass = mutation_class(wt_aa, mut_aa)

        shared = {
            "codon_position": aa_pos, "wt_codon": wt_codon, "mut_codon": mut_codon,
            "mutant": aa_mutant, "pos": aa_pos, "wt": aa_symbol(wt_aa),
            "mutation_class": mclass,
        }

        if mclass == "SYNONYMOUS":
            if skip_codon:
                # Route to protein TSV: synonymous AA score is trivially 0.
                prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
                prow.update({"pkey": pkey, "nt_mutant": nt_mut})
                prow.update(shared)
                prow["subs"] = aa_symbol(mut_aa)
                prow.update({
                    "prediction_protein_independent_adabm": 0.0,
                    "prediction_protein_epistatic_adabm":   0.0,
                    "protein_pairwise_contribution_adabm":  0.0,
                    "protein_concordance_adabm":            "NEUTRAL",
                })
                qc_flags.append("SYNONYMOUS_PROTEIN_LEVEL")
                prow["qc_flags"] = ";".join(qc_flags) if qc_flags else "SYNONYMOUS_PROTEIN_LEVEL"
                protein_rows.append(prow)
            else:
                crow = {f: "" for f in CODON_FIELDNAMES_ADABM}
                crow.update({"pkey": pkey, "nt_mutant": nt_mut})
                crow.update(shared)
                if codon_lookup is not None:
                    scored = codon_lookup.get((aa_pos, mut_codon))
                    if scored is None:
                        qc_flags.append("SYNONYMOUS_NOT_IN_CODON_MODEL")
                    else:
                        crow.update({
                            "prediction_codon_independent_adabm": scored["indep"],
                            "prediction_codon_epistatic_adabm":   scored["epi"],
                            "codon_pairwise_contribution_adabm":  scored["pair"],
                            "codon_concordance_adabm":            scored["concordance"],
                            "codon_frequency_adabm":              scored["frequency"],
                        })
                        qc_flags.append("SYNONYMOUS_SCORED")
                else:
                    qc_flags.append("SYNONYMOUS_UNSCORED")
                crow["qc_flags"] = ";".join(qc_flags)
                codon_rows.append(crow)

        elif mclass in {"STOP_GAIN", "STOP_LOSS"}:
            if skip_codon:
                prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
                prow.update({"pkey": pkey, "nt_mutant": nt_mut})
                prow.update(shared)
                prow["subs"] = aa_symbol(mut_aa)
                qc_flags.append(mclass)
                prow["qc_flags"] = ";".join(qc_flags)
                protein_rows.append(prow)
            else:
                crow = {f: "" for f in CODON_FIELDNAMES_ADABM}
                crow.update({"pkey": pkey, "nt_mutant": nt_mut})
                crow.update(shared)
                qc_flags.append(mclass)
                crow["qc_flags"] = ";".join(qc_flags)
                codon_rows.append(crow)

        else:
            prow = {f: "" for f in PROTEIN_FIELDNAMES_ADABM}
            prow.update({"pkey": pkey, "nt_mutant": nt_mut})
            prow.update(shared)
            prow["subs"] = aa_symbol(mut_aa)
            if mclass == "UNKNOWN":
                qc_flags.append("UNKNOWN_CODON")
            elif aa_lookup is None:
                qc_flags.append("NO_PROTEIN_MODEL")
            else:
                scored = aa_lookup.get((aa_pos, aa_symbol(mut_aa)))
                if scored is None:
                    qc_flags.append("NOT_IN_MODEL")
                else:
                    prow.update({
                        "prediction_protein_independent_adabm": scored["indep"],
                        "prediction_protein_epistatic_adabm":   scored["epi"],
                        "protein_pairwise_contribution_adabm":  scored["pair"],
                        "protein_concordance_adabm":            scored["concordance"],
                        "frequency_adabm":                      scored["frequency"],
                    })
                    qc_flags.append("PASS")
            prow["qc_flags"] = ";".join(qc_flags) if qc_flags else "PASS"
            protein_rows.append(prow)

    return protein_rows, codon_rows


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def _read_orf_sequence(fasta_path: str) -> Tuple[str, str]:
    seqs = read_fasta(fasta_path)
    if not seqs:
        raise RuntimeError(f"No sequences found in FASTA: {fasta_path}")
    for sid, seq in seqs.items():
        if sid.upper() == "ORF":
            return sid, seq.upper()
    first_id, first_seq = next(iter(seqs.items()))
    return first_id, first_seq.upper()


def _build_codon_lookup_from_params(
    codon_params_path: str,
    codon_msa_encoded_path: str,
    codon_focus_id: str,
) -> Dict[Tuple[int, str], Dict]:
    """
    Build a codon (codon_pos, mut_codon_3letter) → score lookup.

    codon_msa_encoded_path must be the SINGLE-CHARACTER encoded MSA produced
    by codon_encoding.encode_codon_msa — same input the params were trained on.
    """
    if not _CODON_ENCODING_AVAILABLE or CODON_ALPHABET is None:
        raise RuntimeError("codon_encoding module not available; cannot score codon side")

    h, J = load_adabmdca_params(codon_params_path, CODON_ALPHABET)
    h, J = apply_zero_sum_gauge(h, J)

    _, focus_seq = _load_focus_sequence(codon_msa_encoded_path, codon_focus_id)
    column_freqs = _column_frequencies(codon_msa_encoded_path, CODON_ALPHABET, encoded=True)

    char_lookup, _ = build_lookup(h, J, focus_seq, CODON_ALPHABET, column_freqs)

    # Re-key from (seq_pos, mut_char) → (seq_pos, mut_codon_triplet) for downstream consumption
    out: Dict[Tuple[int, str], Dict] = {}
    for (sp, mut_char), score in char_lookup.items():
        triplet = CHAR_TO_CODON.get(mut_char)
        if triplet is None:
            continue
        out[(sp, triplet)] = score
    return out


def _build_protein_lookup_from_params(
    protein_params_path: str,
    protein_msa_path: str,
    protein_focus_id: str,
    alphabet: str = "-ACDEFGHIKLMNPQRSTVWY",
) -> Dict[Tuple[int, str], Dict]:
    """
    Build a protein (aa_pos, mut_aa_symbol) → score lookup from the protein MSA.
    """
    h, J = load_adabmdca_params(protein_params_path, alphabet)
    h, J = apply_zero_sum_gauge(h, J)

    _, focus_seq = _load_focus_sequence(protein_msa_path, protein_focus_id)
    column_freqs = _column_frequencies(protein_msa_path, alphabet, encoded=False)

    return build_lookup(h, J, focus_seq, alphabet, column_freqs)[0]


def _process_gene(
    gene: str,
    fasta_file: str,
    mutations_file: str,
    output_dir: Path,
    args: argparse.Namespace,
) -> Tuple[int, int]:
    _, orf_seq = _read_orf_sequence(fasta_file)
    nt_mutations = trim_muts(mutations_file, log=args.validation_log, gene_name=gene)
    if not nt_mutations:
        print("  (no mutations) -> skipping")
        return 0, 0
    print(f"  Loaded {len(nt_mutations)} NT mutations")

    # Resolve params (training if needed) + the MSA paths the scorer reads.
    protein_params_path, codon_params_path, protein_msa_path, codon_msa_path = \
        _resolve_per_gene_adabm_params(gene, args)

    aa_lookup = None
    if protein_params_path and protein_msa_path:
        if not args.quiet:
            print(f"  Loading protein adabmDCA params: {protein_params_path}")
        aa_lookup = _build_protein_lookup_from_params(
            protein_params_path, protein_msa_path, args.focus or gene, args.protein_alphabet
        )

    codon_lookup = None
    if codon_params_path and codon_msa_path and not args.skip_codon:
        if not args.quiet:
            print(f"  Loading codon adabmDCA params: {codon_params_path}")
        codon_lookup = _build_codon_lookup_from_params(
            codon_params_path, codon_msa_path, args.codon_focus or "ORF",
        )

    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    protein_rows, codon_rows = score_nt_mutations_adabm(
        nt_mutations, gene, orf_seq,
        aa_lookup=aa_lookup, codon_lookup=codon_lookup,
        failure_map=failure_map, skip_codon=args.skip_codon,
    )

    gene_dir = output_dir / gene / "adabmDCA"
    gene_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(protein_rows, str(gene_dir / f"{gene}.protein.tsv"),
              PROTEIN_FIELDNAMES_ADABM, extrasaction="ignore")
    print(f"  Wrote {len(protein_rows)} rows -> {gene_dir / (gene + '.protein.tsv')}")
    if not args.skip_codon:
        write_tsv(codon_rows, str(gene_dir / f"{gene}.codon.tsv"),
                  CODON_FIELDNAMES_ADABM, extrasaction="ignore")
        print(f"  Wrote {len(codon_rows)} rows -> {gene_dir / (gene + '.codon.tsv')}")

    return len(protein_rows), len(codon_rows)


# ---------------------------------------------------------------------------
# Discovery helpers (mirror evmutation_pipeline._find_file_for_gene)
# ---------------------------------------------------------------------------

def _find_file_for_gene(gene: str, path_arg: Optional[str], glob_patterns: List[str]) -> Optional[str]:
    if not path_arg:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        for pattern in glob_patterns:
            for f in sorted(p.glob(pattern)):
                if extract_gene_from_filename(f.name).upper() == gene.upper():
                    return str(f)
    return None


def _resolve_params_file(gene: str, path_arg: Optional[str], suffix: str) -> Optional[str]:
    """
    Resolve a params file from a direct path or directory.

    File: use directly.
    Directory: looks for {gene}{suffix} inside it.
    Mirrors evmutation_pipeline._resolve_model_params.
    """
    if not path_arg:
        return None
    p = Path(path_arg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        candidate = p / f"{gene}{suffix}"
        if candidate.exists():
            return str(candidate)
    return None


def _build_adabmdca_params(gene: str, msa_file: str, focus: str, output_path: str,
                            args: argparse.Namespace, codon: bool = False) -> str:
    """
    Encode (codon only) and run adabmDCA train to build params.

    Mirrors evmutation_pipeline._build_model_params. The codon side is encoded
    via codon_encoding.encode_codon_msa before being fed to adabmDCA train with
    the 65-char codon alphabet. Returns output_path (may not exist if
    --skip-train was set and no pre-existing params are found).
    """
    if codon:
        if not _CODON_ENCODING_AVAILABLE:
            raise RuntimeError("codon_encoding module not found in bin/")
        encoded = msa_file + ".encoded.fasta"
        if not args.quiet:
            print(f"  Encoding codon MSA: {msa_file} -> {encoded}")
        n_seqs = encode_codon_msa(msa_file, encoded)
        if not args.quiet:
            print(f"    Encoded {n_seqs} sequences")
        effective_msa = encoded
        alphabet = CODON_ALPHABET
    else:
        effective_msa = msa_file
        alphabet = args.protein_alphabet

    if os.path.exists(output_path):
        if not args.quiet:
            label = "codon " if codon else ""
            print(f"  Using existing {label}adabmDCA params: {output_path}")
        return output_path

    if args.skip_train:
        return output_path  # caller checks existence

    # Dispatch on the chosen training algorithm. bmDCA/eaDCA/edDCA are
    # Boltzmann-learning variants and go through train_adabmdca_in_process.
    # pseudoDCA uses pseudolikelihood maximization — different algorithm,
    # different memory footprint (no MCMC chains).
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    if args.adabmdca_model == "pseudoDCA":
        train_adabmdca_pseudolikelihood_in_process(
            msa_path=effective_msa,
            alphabet=alphabet,
            output_params_path=output_path,
            nepochs=args.adabmdca_nepochs,
            lr=args.adabmdca_lr,
            device_str=args.adabmdca_device,
            dtype_str=args.adabmdca_dtype,
            seed=args.adabmdca_seed,
            quiet=args.quiet,
        )
    else:
        train_adabmdca_in_process(
            msa_path=effective_msa,
            alphabet=alphabet,
            output_params_path=output_path,
            model=args.adabmdca_model,
            nepochs=args.adabmdca_nepochs,
            target=args.adabmdca_target,
            lr=args.adabmdca_lr,
            nchains=args.adabmdca_nchains,
            nsweeps=args.adabmdca_nsweeps,
            device_str=args.adabmdca_device,
            dtype_str=args.adabmdca_dtype,
            seed=args.adabmdca_seed,
            quiet=args.quiet,
        )
    return output_path


def _resolve_per_gene_adabm_params(gene: str, args: argparse.Namespace) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Resolve (and optionally build) adabmDCA params for one gene.

    Mirrors evmutation_pipeline._resolve_per_gene_models, but for adabmDCA.
    Returns (protein_params, codon_params, protein_msa_file, codon_msa_file).
    MSA paths are returned alongside the params because the scoring step
    needs them for column-frequency + focus-sequence extraction.
    """
    # -- Protein --
    protein_params: Optional[str] = None
    protein_msa: Optional[str] = None
    if args.msa:
        protein_msa = _find_file_for_gene(
            gene, args.msa, ["*.a2m", "*.msa.fasta", "*.msa.a2m", "*.fasta", "*.fa", "*.fas"]
        )
        if protein_msa:
            focus = args.focus or gene
            if args.protein_params:
                target = _resolve_params_file(gene, args.protein_params, ".protein_adabm_params") \
                         or str(Path(args.protein_params) / f"{gene}.protein_adabm_params")
            else:
                target = str(Path(protein_msa).parent / f"{gene}.protein_adabm_params")
            _build_adabmdca_params(gene, protein_msa, focus, target, args, codon=False)
            protein_params = target if os.path.exists(target) else None
        elif not args.quiet:
            print(f"  Warning: no protein MSA found for {gene} in {args.msa}")
    elif args.protein_params:
        protein_params = _resolve_params_file(gene, args.protein_params, ".protein_adabm_params")

    # -- Codon --
    codon_params: Optional[str] = None
    codon_msa_encoded: Optional[str] = None  # encoded (single-char) form for scoring
    if args.codon_msa and not args.skip_codon:
        codon_msa_raw = _find_file_for_gene(
            gene, args.codon_msa,
            ["*.codon.msa.fasta", "*.codon.fasta", "*.codon.fa",
             "*.codon.msa.fa", "*.fasta", "*.fa", "*.fas"]
        )
        if codon_msa_raw:
            codon_focus = args.codon_focus or "ORF"
            if args.codon_params:
                target_codon = _resolve_params_file(gene, args.codon_params, ".codon_adabm_params") \
                               or str(Path(args.codon_params) / f"{gene}.codon_adabm_params")
            else:
                target_codon = str(Path(codon_msa_raw).parent / f"{gene}.codon_adabm_params")
            _build_adabmdca_params(gene, codon_msa_raw, codon_focus, target_codon, args, codon=True)
            codon_params = target_codon if os.path.exists(target_codon) else None
            # `_build_adabmdca_params(codon=True)` writes the encoded MSA at
            # `<raw>.encoded.fasta`. Scoring needs that, not the raw triplets.
            encoded_candidate = codon_msa_raw + ".encoded.fasta"
            if os.path.exists(encoded_candidate):
                codon_msa_encoded = encoded_candidate
            elif _CODON_ENCODING_AVAILABLE:
                # Encode on demand (e.g., --skip-train with prebuilt params and no cached encoding)
                encode_codon_msa(codon_msa_raw, encoded_candidate)
                codon_msa_encoded = encoded_candidate
        elif not args.quiet:
            print(f"  Warning: no codon MSA found for {gene} in {args.codon_msa}")
    elif args.codon_params and not args.skip_codon:
        codon_params = _resolve_params_file(gene, args.codon_params, ".codon_adabm_params")

    return protein_params, codon_params, protein_msa, codon_msa_encoded


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="BioFeatureFactory: adabmDCA mutation pipeline (train + score)",
    )
    parser.add_argument("--fasta", required=True,
                        help="ORF FASTA file (single gene) or directory")
    parser.add_argument("--mutations", required=True,
                        help="Mutations CSV file or directory")
    # MSAs are the primary inputs: training is invoked here if params don't exist.
    parser.add_argument("--msa",
                        help="Protein MSA file or directory; triggers `adabmDCA train` if params absent")
    parser.add_argument("--codon-msa",
                        help="Codon MSA file or directory (triplets); encoded + trained if params absent")
    # Pre-built params (cache / override) — skip train if these resolve.
    parser.add_argument("--protein-params",
                        help="Pre-built protein params file or directory")
    parser.add_argument("--codon-params",
                        help="Pre-built codon params file or directory")
    parser.add_argument("--protein-alphabet", default="-ACDEFGHIKLMNPQRSTVWY",
                        help="Protein alphabet (default: standard 20 AA + gap)")
    parser.add_argument("--focus",
                        help="Focus sequence ID in protein MSA (default: gene name)")
    parser.add_argument("--codon-focus",
                        help="Focus sequence ID in codon MSA (default: ORF)")
    parser.add_argument("--gene", help="Gene name override (single-gene mode)")
    parser.add_argument("--skip-codon", action="store_true",
                        help="Route synonymous + stop to the protein TSV; do not produce codon TSV")
    parser.add_argument("--skip-train", action="store_true",
                        help="Skip adabmDCA train invocation (encoding still runs for codon MSAs)")
    # adabmDCA training tunables (consumed by train_adabmdca_in_process when training fires)
    parser.add_argument("--adabmdca-model", default="bmDCA",
                        choices=["bmDCA", "eaDCA", "edDCA", "pseudoDCA"],
                        help="Training algorithm. bmDCA/eaDCA/edDCA are Boltzmann-learning "
                             "variants (high memory); pseudoDCA is pseudolikelihood maximization "
                             "(~2× less peak memory, no MCMC).")
    parser.add_argument("--adabmdca-nepochs", type=int, default=50000)
    parser.add_argument("--adabmdca-target", type=float, default=0.95)
    parser.add_argument("--adabmdca-lr", type=float, default=0.01)
    parser.add_argument("--adabmdca-nchains", type=int, default=10000)
    parser.add_argument("--adabmdca-nsweeps", type=int, default=10)
    parser.add_argument("--adabmdca-device", default="cuda")
    parser.add_argument("--adabmdca-dtype", default="float32",
                        choices=["float32", "float64"])
    parser.add_argument("--adabmdca-seed", type=int, default=0)
    parser.add_argument("--validation-log",
                        help="Validation log from exon_aware_mapping for mutation filtering")
    parser.add_argument("--output", "-o", default=".", help="Output directory")
    parser.add_argument("--quiet", "-q", action="store_true")

    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if fasta_path.is_dir():
        fasta_map = discover_fasta_files(str(fasta_path))
        if not fasta_map:
            print(f"No FASTA files found in {fasta_path}", file=sys.stderr)
            sys.exit(1)

        for gene, fasta_file in sorted(fasta_map.items()):
            print(f"\n== {gene} ==")
            try:
                mutations_file = _find_file_for_gene(gene, args.mutations, ["*.csv"])
                if not mutations_file:
                    print("  No mutations file found -> skipping")
                    continue
                _process_gene(gene, fasta_file, mutations_file, output_dir, args)
            except Exception as exc:
                print(f"  ERROR: {exc}", file=sys.stderr)
        return

    if not fasta_path.is_file():
        print(f"Error: --fasta path not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    gene = args.gene or extract_gene_from_filename(str(fasta_path)) or "GENE"
    print(f"\n== {gene} ==")

    mutations_file = _find_file_for_gene(gene, args.mutations, ["*.csv"])
    if not mutations_file:
        print(f"Error: mutations file not found (gene={gene}, arg={args.mutations})", file=sys.stderr)
        sys.exit(1)

    _process_gene(gene, str(fasta_path), mutations_file, output_dir, args)


if __name__ == "__main__":
    main()
