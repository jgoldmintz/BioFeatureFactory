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
import ...`), NOT the `adabmDCA` console script -- either Boltzmann learning
(bmDCA/eaDCA/edDCA) or pseudolikelihood maximization (pseudoDCA), PyTorch/GPU
native. This module then loads the resulting text-format params and produces
per-mutation TSVs with the same routing logic as evmutation_pipeline.py:

  protein TSV  missense (and synonymous + stop when --skip-codon is set)
  codon TSV    synonymous + stop variants only

Math per mutation (i = position, a = WT token, b = mutant token):
  DeltaH_independent = h_i(b) - h_i(a)
  DeltaH_epistatic   = DeltaH_independent + Sum_{j!=i} [J_ij(b, s_j) - J_ij(a, s_j)]
                   s_j is the WT (focus) token at column j
  pairwise       = DeltaH_epistatic - DeltaH_independent

Zero-sum gauge applied before scoring:
  h_i(a)    -= mean_a h_i(a)
  J_ij(a,b) -= row_mean_b + col_mean_a - global_mean (per (i,j))

Concordance (per-position relative threshold, matches EVmutation convention):
  threshold = 0.5 * std(pairwise) at this position
  |pairwise| <= threshold or std == 0 -> NEUTRAL
  sign(epistatic) == sign(independent) -> CONCORDANT
  otherwise -> DISCORDANT
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
    derive_mutations_root,
    codon_to_aa,
    discover_fasta_files,
    discover_mutation_files,
    find_gene_file,
    extract_gene_from_filename,
    format_aa_token,
    load_validation_failures,
    mint_pkey,
    mutation_class,
    parse_variant,
    protein_consequence,
    read_fasta,
    should_skip_mutation,
    splice_seq,
    split_intronic_tokens,
    trim_muts,
    warn_intronic_unsupported,
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
# successes). We call the training building blocks directly here -- same
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

    # 4. Checkpoint object -- periodic + final save. file_paths["params"] is
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
    # OOM catch: Boltzmann learning's peak memory is ~4-5 x L^2q^2. For codon
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
            raise  # not an OOM -- propagate the original error untouched
        peak_low  = 4 * L * L * q * q * 4 / 1e9
        peak_high = 5 * L * L * q * q * 4 / 1e9
        pl_low    = 2 * L * L * q * q * 4 / 1e9
        pl_high   = 3 * L * L * q * q * 4 / 1e9
        raise RuntimeError(
            f"adabmDCA Boltzmann training OOM'd at L={L}, q={q}.\n"
            f"  Peak memory for Boltzmann is ~4-5 x L^2q^2 x 4 bytes (float32) = "
            f"~{peak_low:.1f}-{peak_high:.1f} GiB.\n"
            f"  Pseudolikelihood mode roughly halves this (~{pl_low:.1f}-{pl_high:.1f} GiB)\n"
            f"  by dropping the MCMC chains + fij_chains buffers. Re-run with:\n"
            f"      --adabmdca-model pseudoDCA\n"
            f"  For genes much larger than ~1200 codons at q=65, even pseudolikelihood\n"
            f"  won't fit on a single GPU -- use the EVmutation/plmc codon backend\n"
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
        # pseudolikelihood save below -- save_params writes every mask-True
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
# buffer, reducing peak GPU memory from ~4-5 x L^2q^2 to ~2-3 x L^2q^2 at the
# same L. The math is per-position conditional likelihood maximization
# (Ekeberg et al. 2013; the approach plmc uses on CPU). GPU-accelerated here
# so it stays fast for the medium-gene regime that doesn't fit Boltzmann.
#
# Note on scaling: the L^2 scaling on the dense J tensor is unchanged. This
# is a constant-factor improvement (~2x), not a scaling-class change. For
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
    tol: float = 1e-3,
    patience: int = 3,
    check_every: int = 10,
) -> str:
    """
    Train an adabmDCA Potts model by pseudolikelihood maximization.

    Parameters mirror Ekeberg 2013 / Hopf 2017 / plmc conventions:
      lambda_h, lambda_J: L2 regularization strengths on h and J.
      nepochs: MAXIMUM number of gradient steps (plmc default: 500). The loop
               stops earlier when the convergence test below fires; nepochs is a
               ceiling, not a target.
      lr: learning rate for plain SGD on h and J.
      chunk_size: positions processed per j-loop chunk in the J gradient
                  accumulator. Larger = faster but more peak intermediate
                  memory (4 x L x q x chunk x q x 4 bytes float32).
      tol, patience, check_every: convergence test. Every `check_every` epochs,
                  ||grad|| is compared against its value at the first check;
                  when the ratio stays below `tol` for `patience` consecutive
                  checks the loop breaks. tol=0 disables the test.

    Convergence: the Boltzmann routines in adabmDCApy self-terminate at
    `target_pearson` (default 0.95), but this pseudolikelihood loop is BFF's own
    reimplementation -- adabmDCApy has no pseudolikelihood routine -- and it ran
    a fixed epoch count with no stopping rule. ||grad|| was already computed at
    the print site and thrown away, so "did 500 epochs converge?" was
    unanswerable by construction. It is now measured.

    The test is RELATIVE to the first checked gradient norm, not absolute: the
    norm scales with L, q and the sequence weights, so any fixed threshold would
    mean something different for every gene. Ratio-to-initial is scale-free and
    comparable across the panel.

    This matters most on rented hardware. `--adabmdca-nepochs` defaults to 50000
    (a Boltzmann-sized number; this function's own signature default is 500), and
    at codon scale the per-epoch cost is the L^2q^2 gradient einsum below, so the
    difference between stopping at convergence and running the ceiling is the
    difference between a run and days of billed compute.

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
    # Memory: B x L x q x 4 bytes. For F9 codon (B=531, L=764, q=65): ~106 MB.
    B = X.shape[0]
    X_oh = torch.nn.functional.one_hot(X, num_classes=q).to(dtype)  # (B, L, q)

    # Heads-up if the user passed a Boltzmann-sized nepochs (default 50000) into
    # pseudolikelihood -- plmc's default is 500. Doesn't fail, just warns.
    if not quiet and nepochs > 5000:
        print(f"    NOTE nepochs={nepochs} is large for pseudolikelihood "
              f"(plmc default: 500). Consider --adabmdca-nepochs 500.")

    # 5. Training loop. Plain SGD on negative-PL + L2.
    #
    # MEMORY-CRITICAL CHANGES vs. the initial version:
    #   (a) torch.no_grad() wraps the loop -- autograd is irrelevant for
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
    # Convergence bookkeeping (see the docstring's Convergence note).
    grad_norm_0 = None      # ||grad|| at the first checked epoch = the baseline
    below = 0               # consecutive checks under tol
    stopped_at = None       # epoch the loop broke at, None if it ran to nepochs
    history = []            # [(epoch, ||grad||)] at each check, for the caller's log
    # Pre-allocate the (B, L, q) E buffer once and reuse via copy_().
    E = torch.empty(B, L, q, device=device, dtype=dtype)
    try:
        with torch.no_grad():
            for epoch in range(nepochs):
                # ----- Forward: compute conditional fields E[s, i, a] -----
                # E[s, i, a] = h[i, a] + Sum_{j != i} J[i, a, j, X[s, j]]
                # Compute the full j-sum, then subtract the diagonal j=i term.

                # Initialize E with the broadcast h (in-place into the reused buffer).
                E.copy_(h.unsqueeze(0).expand(B, -1, -1))

                for j in range(L):
                    # J[:, :, j, X[:, j]] -> (L, q, B) view -> permute to (B, L, q),
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

                # Per-sequence summed log-likelihood: Sum_i log p(X[s, i] | s_{-i})
                log_p_at_data = log_p.gather(-1, X.unsqueeze(-1)).squeeze(-1)  # (B, L)
                PL = (weights.unsqueeze(-1) * log_p_at_data).sum() / weights_sum

                # ----- Gradient of PL (data - model) weighted by sequence weights -----
                diff = X_oh - p                                          # (B, L, q)
                # weights/weights_sum, not weights: PL above is normalized by
                # weights_sum, so the gradient must be too -- otherwise lr and
                # the L2 term below are effectively scaled by M_eff.
                wdiff = diff * (weights / weights_sum).view(-1, 1, 1)    # weighted (B, L, q)
                # Free intermediates we no longer need so PyTorch's allocator can reuse.
                del log_p, p, log_p_at_data, diff, logZ

                # dPL/dh[i, a] = Sum_s w_s * (X_oh[s, i, a] - p[s, i, a]) / weights_sum
                torch.sum(wdiff, dim=0, out=grad_h)                      # in-place into grad_h

                # dPL/dJ[i, a, j, b] -- symmetrized as in Ekeberg 2013 eq. 9.
                # Term1[i, a, j, b] = Sum_s w_s * diff[s, i, a] * X_oh[s, j, b] / weights_sum
                # Symmetrized: grad_J_sym[i,a,j,b] = 0.5*(Term1[i,a,j,b] + Term1[j,b,i,a]).
                #
                # We build the symmetrized gradient DIRECTLY via per-chunk dual
                # accumulation so we never need to materialize the transpose of
                # the full (L, q, L, q) grad_J (which would double peak memory).
                # For each j-chunk:
                #   - 0.5 * Term1[i, a, j_chunk, b]  -> grad_J[:, :, j_chunk, :]
                #   - 0.5 * Term1[j_chunk, b, i, a]  -> grad_J[j_chunk, :, :, :]
                #     (achieved by permuting the chunk axes: (L, q, chunk, q) ->
                #      (chunk, q, L, q) and broadcasting into the j_chunk slice
                #      of dim 0.)
                grad_J.zero_()
                for j_start in range(0, L, chunk_size):
                    j_end = min(j_start + chunk_size, L)
                    X_oh_chunk = X_oh[:, j_start:j_end, :]               # (B, chunk, q)
                    term1_chunk = torch.einsum(
                        "sia,sjb->iajb", wdiff, X_oh_chunk
                    )                                                     # (L, q, chunk, q)
                    # Contribution of 0.5 * Term1 to the (*, *, j_chunk, *) slot:
                    grad_J[:, :, j_start:j_end, :].add_(term1_chunk, alpha=0.5)
                    # Contribution of 0.5 * Term1^T to the (j_chunk, *, *, *) slot:
                    grad_J[j_start:j_end, :, :, :].add_(
                        term1_chunk.permute(2, 3, 0, 1), alpha=0.5
                    )
                    del term1_chunk
                grad_J[L_idx, :, L_idx, :] = 0.0

                # ----- L2 regularization + SGD step (all in-place) -----
                # max (PL - 0.5 lambda ||.||^2)  =>  param <- param*(1 - lr*lambda) + lr*grad
                h.mul_(1.0 - lr * lambda_h).add_(grad_h, alpha=lr)
                J.mul_(1.0 - lr * lambda_J).add_(grad_J, alpha=lr)
                J[L_idx, :, L_idx, :] = 0.0  # keep diagonal at zero after the update

                # ----- Convergence test -----
                # Evaluated on a schedule, not every epoch: ||grad|| is a reduction
                # over the full (L,q,L,q) tensor, which is cheap next to the epoch's
                # einsum but not free at codon scale. check_every amortises it.
                is_check = (tol > 0 and (epoch % check_every == 0 or epoch == nepochs - 1))
                grad_norm = None
                if is_check or (not quiet and epoch % 25 == 0):
                    grad_norm = (grad_h.norm() + grad_J.norm()).item()
                if not quiet and (epoch % 25 == 0 or epoch == nepochs - 1):
                    if grad_norm is None:
                        grad_norm = (grad_h.norm() + grad_J.norm()).item()
                    ratio = "" if grad_norm_0 is None else f", ||grad||/||grad||_0 = {grad_norm / grad_norm_0:.3e}"
                    print(f"    epoch {epoch:4d}: PL = {PL.item():.4f}, "
                          f"||grad|| = {grad_norm:.4e}{ratio}")
                if is_check:
                    if grad_norm_0 is None:
                        grad_norm_0 = grad_norm
                    # PLATEAU detector, not ratio-to-initial. Measured on a toy MSA,
                    # ||grad|| under this plain-SGD update decays roughly as
                    # 1/sqrt(epoch): 1.46 -> 0.36 over 400 epochs, i.e. ratio-to-initial
                    # 0.78 @25, 0.49 @100, 0.34 @200, 0.25 @400. Reaching 1e-3 OF THE
                    # INITIAL NORM would take ~1e5-1e6 epochs, so an absolute-ratio test
                    # can never fire and the early stop would be decorative.
                    # What "converged" actually means here is that the norm has stopped
                    # moving, so the test is the RELATIVE DECREASE PER CHECK.
                    rel_drop = None
                    if history:
                        prev = history[-1][1]
                        if prev > 0:
                            rel_drop = (prev - grad_norm) / prev
                    history.append((epoch, grad_norm))
                    if rel_drop is not None and rel_drop < tol:
                        below += 1
                        if below >= patience:
                            stopped_at = epoch
                            if not quiet:
                                print(f"    converged: ||grad|| improved <{tol:.1e} per "
                                      f"{check_every} epochs for {patience} consecutive "
                                      f"checks; stopping at epoch {epoch} of {nepochs} "
                                      f"(||grad||/||grad||_0 = {grad_norm / grad_norm_0:.3e})")
                            del wdiff
                            break
                    else:
                        # Reset on any excursion -- patience counts CONSECUTIVE checks, so
                        # a single noisy plateau cannot end the run early. rel_drop < 0
                        # (the norm went UP) also resets rather than counting.
                        below = 0
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
            f"  Peak memory ~{pl_low:.1f}-{pl_high:.1f} GiB. GPU couldn't supply.\n"
            f"  Pseudolikelihood is already the more memory-efficient option;\n"
            f"  this gene exceeds the single-GPU dense-Potts ceiling. Use the\n"
            f"  EVmutation/plmc codon backend instead -- it runs the same\n"
            f"  pseudolikelihood math on CPU with system-RAM-scale budgets.\n"
            f"  Original error: {e}"
        ) from e

    if not quiet:
        if tol <= 0:
            print(f"    epochs run: {nepochs} of {nepochs} (early stopping disabled, tol=0)")
        elif stopped_at is None:
            # Ran the ceiling without the test firing. Say so explicitly: it means the
            # model is STILL DESCENDING at the ceiling, which is exactly the thing the
            # fixed-epoch loop could not tell you. Report the observed decay so the
            # next run can be sized from data instead of guessed.
            last_drop = ""
            if len(history) >= 2 and history[-2][1] > 0:
                d = (history[-2][1] - history[-1][1]) / history[-2][1]
                last_drop = f", last decrease {d:.2e} per {check_every} epochs (tol={tol:.1e})"
            ratio = f", ||grad||/||grad||_0 = {history[-1][1] / grad_norm_0:.3e}" \
                    if history and grad_norm_0 else ""
            print(f"    NOT converged: ran all {nepochs} epochs and ||grad|| was still "
                  f"descending{ratio}{last_drop}. Raise --adabmdca-nepochs, or loosen "
                  f"--adabmdca-tol if this decrease is small enough for your purpose.")
        else:
            print(f"    epochs run: {stopped_at + 1} of {nepochs} "
                  f"({(nepochs - stopped_at - 1) / nepochs:.0%} of the ceiling saved)")

    # 6. Save the trained params to disk in adabmDCApy text format so the
    # downstream scoring step can load them the same way as the Boltzmann path.
    final_params = {"bias": h, "coupling_matrix": J}
    save_params(output_params_path, final_params, tokens, mask)
    return output_params_path


# ---------------------------------------------------------------------------
# adabmDCApy text-format params loader (numpy-only, no torch dependency)
# ---------------------------------------------------------------------------

def load_adabmdca_params(path: str, alphabet: str,
                         dtype: np.dtype = np.float32) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse adabmDCApy text-format params.

    Format:
      'J idx0 idx1 aa0 aa1 value'   only upper-triangular i<j entries saved
      'h idx aa value'

    Returns:
      h: (L, q) array
      J: (L, q, L, q) array (symmetrized)

    MEMORY. J is dense (L, L, q, q) and that is the whole cost of this function:

        gene          L     q   one J    old peak   this peak
        SMN2 codon   294    64   2.6 G      5.3 G       1.3 G
        F9   codon   461    64   6.5 G     13.0 G       3.2 G
        PAM  codon   974    64  29.0 G     57.9 G      14.5 G
        BRCA1 codon 1863    64 105.9 G    211.8 G      53.0 G
        BRCA1 aa    1863    21  11.4 G     22.8 G       5.7 G

    Two changes get that 4x:

    1. The symmetrization used to be `J = J + J.transpose(1, 0, 3, 2)`, which
       allocates a SECOND full-size array before the first is released -- peak
       2x, for an operation that writes each value exactly once. The mirrored
       write below fills both halves in the record loop instead, so peak is 1x.
       Exact for the documented format: entries are strictly upper-triangular
       (i < j), so `J[i,j,a,b] = v; J[j,i,b,a] = v` and "assign upper, then add
       the transpose" produce the same array -- the transpose contributes 0 to
       the upper half and v to the lower.
    2. float32 rather than float64. These are DCA couplings consumed by
       _delta_hamiltonian_multi as a sum of differences; float64 was never a
       precision requirement, and params.adabmdca_dtype in bin/main.nf already
       defaults the TRAINING side to float32. Pass dtype=np.float64 to restore
       the old precision.

    NOT fixed here: J is still materialized in full. _delta_hamiltonian_multi
    only ever touches rows for mutated positions, so streaming straight from
    j_records would make this O(mutated sites) instead of O(L^2 q^2) and is what
    a codon-alphabet BRCA1 actually needs. That is a change to the scoring path,
    not to this loader.
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

    h = np.zeros((L, q), dtype=dtype)
    for idx, ai, val in h_records:
        h[idx, ai] = val

    # Mirrored write, not assign-then-add-transpose: see MEMORY in the docstring.
    J = np.zeros((L, L, q, q), dtype=dtype)
    for i, j, ai, bj, val in j_records:
        J[i, j, ai, bj] = val
        J[j, i, bj, ai] = val

    # Reorganize to (L, q, L, q) -- matches EVmutation's J_ij[i, a, j, b] layout.
    # transpose returns a VIEW, so this costs nothing; the copy happens later in
    # apply_zero_sum_gauge, which is the other place peak memory is set.
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
    focus has a gap are skipped -- they don't correspond to a sequence position.

    Returns:
      lookup: dict keyed by (seq_pos, mut_token_char) -> score dict
      col_to_seq: dict mapping alignment column -> seq_pos (1-based)
    """
    L_align = h.shape[0]
    if len(focus_seq) != L_align:
        raise RuntimeError(
            f"Focus sequence length {len(focus_seq)} does not match params L {L_align}"
        )

    token_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    gap_idx = token_to_idx.get("-")

    # Map alignment column -> CDS/AA position (1-based, skipping gap columns in focus)
    col_to_seq: Dict[int, int] = {}
    seq_pos = 0
    for col, ch in enumerate(focus_seq):
        idx = token_to_idx.get(ch)
        if idx is None or idx == gap_idx:
            continue
        seq_pos += 1
        col_to_seq[col] = seq_pos

    # Vectorize WT token indices over all alignment columns (gap -> 0; we'll skip via col_to_seq)
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

        # Sum_{j!=col_i} J[col_i, a, j, s_j]  for each a
        # J_i[:, j, wt_idx_per_col[j]]  -> shape (q, L)
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


# ---------------------------------------------------------------------------
# Multi-site Hamiltonian
# ---------------------------------------------------------------------------

class ModelContext:
    """Everything needed to evaluate the Hamiltonian at arbitrary sites.

    build_lookup already computes all of this and used to discard it: the
    builders returned `[0]`, keeping only the single-mutant dict, so a change
    spanning more than one site had no key and could never acquire one. This
    carries the arrays forward instead.

    Fields:
      h, J           zero-sum-gauged parameters. J is (L, q, L, q) indexed
                     J[i, a, j, b] -- the same layout EVmutation's CouplingsModel
                     uses, which is why _delta_hamiltonian_multi ports across with
                     only the index lookup changed.
      pos_to_col     1-based sequence position -> alignment column. The inverse of
                     build_lookup's col_to_seq. Alignment columns are NOT sequence
                     positions whenever the focus row carries a gap, and handing a
                     position straight to J silently evaluates a different site.
      wt_idx_per_col alphabet index of the focus symbol at each column, i.e. the
                     wild-type background the deltas are taken against.
      column_freqs   (L, q) non-gap column frequencies, for the frequency column.
      alphabet       symbol order; index 0 is the gap in both BFF alphabets.
    """

    __slots__ = ("h", "J", "pos_to_col", "wt_idx_per_col", "column_freqs", "alphabet")

    def __init__(self, h, J, col_to_seq, wt_idx_per_col, column_freqs, alphabet):
        self.h = h
        self.J = J
        self.pos_to_col = {sp: col for col, sp in col_to_seq.items()}
        self.wt_idx_per_col = wt_idx_per_col
        self.column_freqs = column_freqs
        self.alphabet = alphabet


def _delta_hamiltonian_multi(ctx: ModelContext, plan) -> Tuple[float, float, float]:
    """Delta statistical energy for SIMULTANEOUS substitutions at several sites.

    plan is [(seq_pos, wt_symbol, mut_symbol), ...] in 1-based sequence positions.
    Returns (epistatic, independent, pairwise) where
        independent = sum_m [ h_i(b_m) - h_i(a_m) ]
        epistatic   = independent + the coupling terms below
        pairwise    = epistatic - independent

    Same arithmetic as evmutation_pipeline._delta_hamiltonian_multi, which in turn
    reproduces the vendored CouplingsModel: for each moved site, the coupling delta
    is taken against the WILD-TYPE background at every other site, then pairs where
    BOTH sites moved are corrected -- remove the double count and re-add the
    coupling in the new background. Without that correction a two-site change is
    not the same as two one-site changes evaluated independently, which is the
    entire reason a multi-site evaluator is needed rather than summing lookups.

    Raises KeyError for a position outside the model or a symbol outside the
    alphabet; ValueError when the model's own focus symbol disagrees with the ORF.
    """
    token_to_idx = {ch: i for i, ch in enumerate(ctx.alphabet)}
    L = ctx.h.shape[0]
    wt = ctx.wt_idx_per_col

    cols, muts = [], []
    for pos, wt_symbol, mut_symbol in plan:
        col = ctx.pos_to_col[pos]                 # KeyError => not in model
        b = token_to_idx[mut_symbol]              # KeyError => not in alphabet
        a = token_to_idx[wt_symbol]
        if wt[col] != a:
            raise ValueError(
                f"model column {col} (seq pos {pos}) carries "
                f"{ctx.alphabet[wt[col]]!r}, ORF declares {wt_symbol!r}")
        cols.append(col)
        muts.append(b)

    all_j = np.arange(L)
    delta_h = 0.0
    delta_J = 0.0
    for m, (i, b) in enumerate(zip(cols, muts)):
        a = int(wt[i])
        delta_h += float(ctx.h[i, b] - ctx.h[i, a])
        # couplings to every other column, in the wild-type background
        new = ctx.J[i, b, all_j, wt]
        old = ctx.J[i, a, all_j, wt]
        delta_J += float(new.sum() - old.sum() - (new[i] - old[i]))   # drop j == i
        # pairs where BOTH sites moved
        for n in range(m + 1, len(cols)):
            j, b2 = cols[n], muts[n]
            a2 = int(wt[j])
            delta_J -= float(ctx.J[i, b, j, a2])
            delta_J -= float(ctx.J[i, a, j, b2])
            delta_J += float(ctx.J[i, a, j, a2])
            delta_J += float(ctx.J[i, b, j, b2])

    epistatic = delta_J + delta_h
    return epistatic, delta_h, delta_J


def _score_plan_adabm(ctx: Optional[ModelContext], plan):
    """Score a substitution plan. Returns (columns, reason); one is meaningful.

    No DELETION_NO_GAP_SYMBOL_IN_MODEL guard, deliberately: plmc's -g strips the
    gap out of the params file, which is why evmutation_pipeline must refuse
    deletions, but adabmDCA keeps '-' at index 0 of both alphabets and
    build_lookup iterates every symbol including it. A deletion is an ordinary
    substitution to the gap state here.
    """
    if ctx is None:
        return {}, "NO_MODEL"
    try:
        epistatic, independent, pairwise = _delta_hamiltonian_multi(ctx, plan)
    except KeyError:
        return {}, "NOT_IN_MODEL"
    except ValueError:
        return {}, "MODEL_WT_MISMATCH"
    freqs = []
    for pos, _, mut_symbol in plan:
        col = ctx.pos_to_col[pos]
        freqs.append(float(ctx.column_freqs[col, ctx.alphabet.index(mut_symbol)]))
    return {
        "epistatic": epistatic,
        "independent": independent,
        "pairwise": pairwise,
        # Comma-joined per site, matching the mutant/pos columns. The gap column is
        # zeroed by _column_frequencies (gaps are excluded from the denominator),
        # so a deleted site reports empty rather than a fabricated 0.0.
        "frequency": ",".join("" if (f == 0.0 and s == "-") else f"{f}"
                              for f, (_, _, s) in zip(freqs, plan)),
    }, None


# ---------------------------------------------------------------------------
# Non-SNV support
# ---------------------------------------------------------------------------
# Ported from evmutation_pipeline.py so both DCA backends label a token the same
# way. What differs is SCORING: EVmutation evaluates a multi-site plan directly
# off a live CouplingsModel via _delta_hamiltonian_multi, and adabmDCA has no
# equivalent -- build_lookup computes h, J, col_to_seq and column_freqs and then
# _build_protein_lookup_from_params returns only [0], discarding all of it. The
# precomputed lookup holds SINGLE mutants keyed by (seq_pos, token), so nothing
# spanning more than one site can be answered from it.
#
# So these rows are ANNOTATED but not scored: mutation_class, the aa consequence,
# the span, and the model form of the change are all recorded, and every metric
# column stays EMPTY with qc_flags naming why. That is strictly better than the
# status quo, where nt_re rejected the token and it was filed as INVALID_MUTATION
# -- a false label, since 'AAA200GGG' is perfectly well formed.

# protein_consequence's class -> the mutation_class column. The classes the
# original vocabulary cannot express are carried under their own names rather
# than flattened into MISSENSE, which means "one residue swapped". 'mnv' IS a
# substitution (nothing added or removed), so it stays MISSENSE; the multiplicity
# shows up in the comma-joined mutant/pos columns and in qc_flags as AA:mnv.
_AA_CONSEQUENCE_TO_CLASS = {
    "synonymous":     "SYNONYMOUS",
    "stop_gained":    "STOP_GAIN",
    "stop_lost":      "STOP_LOSS",
    "snv":            "MISSENSE",
    "mnv":            "MISSENSE",
    "inframe_del":    "INFRAME_DEL",
    "inframe_ins":    "INFRAME_INS",
    "inframe_delins": "INFRAME_DELINS",
    "frameshift":     "FRAMESHIFT",
}


def _join(values) -> str:
    """Comma-join per-site values. A single-site change produces no comma, so an
    SNV row is unchanged."""
    return ",".join(str(v) for v in values)


def _blank_row(fieldnames: List[str], pkey: str, nt_mut: str, qc_flags: List[str]) -> Dict:
    row = {f: "" for f in fieldnames}
    row.update({"pkey": pkey, "nt_mutant": nt_mut, "qc_flags": ";".join(qc_flags)})
    return row


def _wt_residue(orf_seq: str, aa_pos: int) -> str:
    """The wild-type residue at a 1-based codon position, '' if there is no codon."""
    codon = orf_seq[(aa_pos - 1) * 3:aa_pos * 3]
    if len(codon) != 3:
        return ""
    return aa_symbol(codon_to_aa.get(codon, "X"))


def _anchor_aa_alleles(aa_pos: int, wt_aa: str, mut_aa: str, orf_seq: str):
    """Re-anchor a one-sided aa change on a neighbouring residue.

    protein_consequence returns the MINIMAL span, so a pure insertion has an empty
    wild-type allele. An empty allele does not render as a token and a consumer
    handed one drops the row. Anchor 5' where there is a preceding residue, 3' at
    residue 1 -- the convention infer_aavariant_from_nt uses for the mapping CSVs.
    """
    if wt_aa and mut_aa:
        return aa_pos, wt_aa, mut_aa
    if aa_pos > 1:
        anchor = _wt_residue(orf_seq, aa_pos - 1)
        if anchor:
            return aa_pos - 1, anchor + wt_aa, anchor + mut_aa
    anchor = _wt_residue(orf_seq, aa_pos + len(wt_aa))
    if anchor:
        return aa_pos, wt_aa + anchor, mut_aa + anchor
    return aa_pos, wt_aa, mut_aa


def _substitution_plan(consequence: str, aa_pos: int, wt_aa: str, mut_aa: str):
    """Write a protein change as fixed-L site substitutions, or refuse by name.

    Returns (plan, refusal); exactly one is meaningful. The mutant residues are laid
    down from the 5' end of the changed span and any residues left over on the
    wild-type side become gaps, so a 2-residue deletion written 'MK'->'' and one
    written 'MKL'->'M' land on the same sites.

    A deleted residue is written '-' regardless of what a model turns out to carry.
    Unlike plmc -- whose -g flag strips the gap out of the params file, which is why
    evmutation_pipeline refuses deletions with DELETION_NO_GAP_SYMBOL_IN_MODEL --
    adabmDCA keeps '-' at index 0 of both alphabets (protein default
    '-ACDEFGHIKLMNPQRSTVWY', CODON_ALPHABET 65 chars), and build_lookup:794 iterates
    every b_idx including the gap. So a single-residue deletion is already present in
    the lookup as (aa_pos, '-'), and multi-residue ones become scorable as soon as the
    Hamiltonian evaluator lands. Do NOT port that refusal here.
    """
    if consequence == "frameshift":
        # wt_aa/mut_aa are deliberately empty (protein_consequence declines to pair
        # residues across a frame change) and every downstream site changes identity,
        # so there is nothing to write even in principle.
        return [], "FRAMESHIFT_NOT_REPRESENTABLE_FIXED_L"
    if len(mut_aa) > len(wt_aa):
        return [], "INSERTION_NOT_REPRESENTABLE_FIXED_L"
    plan = [(aa_pos + i, wt_aa[i], mut_aa[i]) for i in range(len(mut_aa))]
    plan += [(aa_pos + i, wt_aa[i], "-") for i in range(len(mut_aa), len(wt_aa))]
    return plan, None


def _non_snv_rows_adabm(pkey, nt_mut, variant, orf_seq, skip_codon,
                        aa_ctx=None, codon_ctx=None):
    """Build the row for one non-SNV token. Returns (protein_row, codon_row).

    Exactly one is not None -- the routing mirrors the SNV path: synonymous and stop
    variants belong to the codon table unless --skip-codon sends them to the protein
    table. Metric columns are always EMPTY (see the module note above); qc_flags
    carries NON_SNV:<kind>, AA:<consequence> and the reason no score exists.
    """
    qc = [f"NON_SNV:{variant.kind}"]

    # Bound the END of the REF span, not its start: a multi-base REF can begin inside
    # the ORF and run off the end, which seq[pos] never catches.
    if variant.pos0 + len(variant.ref) > len(orf_seq):
        return _blank_row(PROTEIN_FIELDNAMES_ADABM, pkey, nt_mut, qc + ["OUT_OF_RANGE"]), None

    first_codon = variant.pos0 // 3
    last_codon = (variant.pos0 + len(variant.ref) - 1) // 3
    # Bound the codon the span ENDS in. A REF that begins in a whole codon can finish
    # inside a ragged final one, and checking only first_codon lets that through as a
    # 4- or 5-base "codon" in wt_codon.
    if last_codon * 3 + 3 > len(orf_seq):
        return _blank_row(PROTEIN_FIELDNAMES_ADABM, pkey, nt_mut, qc + ["PARTIAL_CODON"]), None

    # Guard the WHOLE REF span, not just its first base. Flag and keep going, which is
    # what the SNV path does: the ORF's own bases still describe a real codon, and
    # dropping the row would hide the disagreement.
    if orf_seq[variant.pos0:variant.pos0 + len(variant.ref)].upper() != variant.ref.upper():
        qc.append("REF_MISMATCH")

    cons = protein_consequence(variant, orf_seq)
    consequence = cons["aa_consequence"]
    aa_pos, wt_aa, mut_aa = cons["aa_pos"], cons["wt_aa"], cons["mut_aa"]
    qc.append(f"AA:{consequence}")
    if consequence == "frameshift" and cons["new_stop_aa_pos"] is not None:
        qc.append(f"NEW_STOP_AA:{cons['new_stop_aa_pos']}")

    mut_orf = splice_seq(orf_seq, variant.pos0, variant.ref, variant.alt, validate=False)
    wt_span = orf_seq[first_codon * 3:(last_codon + 1) * 3]
    delta = variant.length_delta
    # A frame-preserving edit shortens or lengthens the mutant span by exactly its
    # length delta. A frameshift does not: the codons after it are re-read rather than
    # removed, so the mutant span is the SAME width and reports what the ribosome now
    # reads over those bases.
    mut_width = len(wt_span) + (delta if delta % 3 == 0 else 0)
    mut_span = mut_orf[first_codon * 3:first_codon * 3 + mut_width]

    shared = {
        "codon_position": first_codon + 1,
        "wt_codon": wt_span,
        "mut_codon": mut_span,
        "mutation_class": _AA_CONSEQUENCE_TO_CLASS.get(consequence, "UNKNOWN"),
    }

    # ---- codon-table classes: synonymous and stop, exactly as for an SNV ----
    if consequence in ("synonymous", "stop_gained", "stop_lost"):
        if skip_codon:
            prow = _blank_row(PROTEIN_FIELDNAMES_ADABM, pkey, nt_mut, qc)
            prow.update(shared)
            if consequence == "synonymous":
                # h_i(M) - h_i(M) = 0 for an unchanged protein: a measured zero, not a
                # coalesced one, which is why it is written rather than left blank.
                prow.update({
                    "prediction_protein_independent_adabm": 0.0,
                    "prediction_protein_epistatic_adabm":   0.0,
                    "protein_pairwise_contribution_adabm":  0.0,
                    "protein_concordance_adabm":            "NEUTRAL",
                })
                qc.append("SYNONYMOUS_PROTEIN_LEVEL")
            else:
                qc.append(_AA_CONSEQUENCE_TO_CLASS[consequence])
            prow["qc_flags"] = ";".join(qc)
            return prow, None

        crow = _blank_row(CODON_FIELDNAMES_ADABM, pkey, nt_mut, qc)
        crow.update({k: v for k, v in shared.items() if k in CODON_FIELDNAMES_ADABM})
        if consequence != "synonymous":
            qc.append(_AA_CONSEQUENCE_TO_CLASS[consequence])
        else:
            codon_plan = [c for c in range(first_codon, last_codon + 1)
                          if orf_seq[c * 3:c * 3 + 3] != mut_orf[c * 3:c * 3 + 3]]
            if not codon_plan:
                # Length-preserving and every codon identical: the token names no change
                # at all. Say so instead of emitting a row of zeros.
                qc.append("NO_CODON_CHANGED")
            else:
                # Both single- and multi-codon changes go through the SAME evaluator,
                # so a token never scores two different ways depending on how many
                # codons it happens to touch. For one codon this reproduces the
                # lookup; for several it is the only route.
                plan = [(c + 1, orf_seq[c * 3:c * 3 + 3], mut_orf[c * 3:c * 3 + 3])
                        for c in codon_plan]
                cols, reason = _score_plan_adabm(codon_ctx, plan)
                if cols:
                    crow.update({
                        "prediction_codon_independent_adabm": cols["independent"],
                        "prediction_codon_epistatic_adabm":   cols["epistatic"],
                        "codon_pairwise_contribution_adabm":  cols["pairwise"],
                        "codon_frequency_adabm":              cols["frequency"],
                        "codon_concordance_adabm": "" if len(plan) > 1 else
                            ("CONCORDANT" if (cols["epistatic"] >= 0) == (cols["independent"] >= 0)
                             else "DISCORDANT"),
                    })
                    qc.append("SYNONYMOUS_SCORED")
                    if len(plan) > 1:
                        qc.append("CONCORDANCE_UNDEFINED_MULTICODON")
                else:
                    qc.append(f"SYNONYMOUS_{reason}")
        crow["qc_flags"] = ";".join(qc)
        return None, crow

    # ---- protein-table classes ----
    prow = _blank_row(PROTEIN_FIELDNAMES_ADABM, pkey, nt_mut, qc)
    prow.update(shared)

    plan, refusal = _substitution_plan(consequence, aa_pos, wt_aa, mut_aa)
    if plan:
        # The model form: one 'M30V' per site moved, comma-joined, a deleted residue
        # written 'K61-'.
        prow.update({"pos":    _join(p for p, _, _ in plan),
                     "wt":     _join(w for _, w, _ in plan),
                     "subs":   _join(s for _, _, s in plan),
                     "mutant": _join(f"{w}{p}{s}" for p, w, s in plan)})
        cols, reason = _score_plan_adabm(aa_ctx, plan)
        if cols:
            prow.update({
                "prediction_protein_independent_adabm": cols["independent"],
                "prediction_protein_epistatic_adabm":   cols["epistatic"],
                "protein_pairwise_contribution_adabm":  cols["pairwise"],
                "frequency_adabm":                      cols["frequency"],
                # Concordance is a comparison against ONE position's pairwise noise
                # floor (build_lookup's per-position std). Several sites moved, so
                # that floor does not exist; left empty and named rather than
                # borrowed from one of them.
                "protein_concordance_adabm": "" if len(plan) > 1 else
                    ("CONCORDANT" if (cols["epistatic"] >= 0) == (cols["independent"] >= 0)
                     else "DISCORDANT"),
            })
            qc.append("PASS")
            if len(plan) > 1:
                qc.append("CONCORDANCE_UNDEFINED_MULTISITE")
        else:
            qc.append(reason)
    else:
        # No model form exists -- that IS the refusal. Fall back to the aa-level token
        # so the row still names the change.
        if consequence == "frameshift":
            wt_res = _wt_residue(orf_seq, aa_pos)
            prow.update({"pos": aa_pos, "wt": wt_res, "subs": "",
                         "mutant": format_aa_token(aa_pos, wt_res, "", "frameshift")})
        else:
            tp, tw, tm = _anchor_aa_alleles(aa_pos, wt_aa, mut_aa, orf_seq)
            prow.update({"pos": tp, "wt": tw, "subs": tm, "mutant": f"{tw}{tp}{tm}"})
        qc.append(refusal)

    prow["qc_flags"] = ";".join(qc)
    return prow, None


def score_nt_mutations_adabm(
    nt_mutations: List[str],
    gene: str,
    orf_seq: str,
    aa_lookup: Optional[Dict[Tuple[int, str], Dict]] = None,
    codon_lookup: Optional[Dict[Tuple[int, str], Dict]] = None,
    failure_map: Optional[Dict] = None,
    skip_codon: bool = False,
    aa_ctx: Optional["ModelContext"] = None,
    codon_ctx: Optional["ModelContext"] = None,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Map and score NT mutations through the adabmDCA backend.

    aa_lookup keys: (aa_pos_1based, mut_aa_symbol) -> {indep, epi, pair, concordance, frequency}
    codon_lookup keys: (aa_pos_1based, mut_codon_3letter) -> {indep, epi, pair, concordance, frequency}

    Routing:
      missense / unknown_codon  -> protein TSV (scored via aa_lookup if present)
      synonymous (skip_codon=False) -> codon TSV
      synonymous (skip_codon=True)  -> protein TSV with score 0 (trivial: M->M)
      stop_gain/loss (skip_codon=False) -> codon TSV (annotation only)
      stop_gain/loss (skip_codon=True)  -> protein TSV with empty scores
    """
    failure_map = failure_map or {}
    nt_re = re.compile(r"^([ACGT])(\d+)([ACGT])$")
    protein_rows: List[Dict] = []
    codon_rows: List[Dict] = []

    for nt_mut in nt_mutations:
        if should_skip_mutation(gene, nt_mut, failure_map):
            continue

        # {GENE}-{sha}. nt_mut comes from trim_muts, which strips only '*' and
        # whitespace, so it is the verbatim token variant_mapping hashed. The old
        # f"{gene}-{nt_mut}" was unbounded in length, which is what broke the SQL key
        # limit on knockout-scale tokens; mint_pkey is len(gene)+13 for every class.
        pkey = mint_pkey(gene, nt_mut)
        qc_flags: List[str] = []

        m = nt_re.match(nt_mut)
        if not m:
            # nt_re is uppercase-ACGT only, so a valid SNV written lowercase or with U
            # (RNA notation) missed it and fell through to INVALID_MUTATION -- a false
            # label, since parse_variant accepts all three spellings. Canonicalise for
            # the GATE only; pkey and nt_mutant keep the user's original spelling, and
            # an already-uppercase token cannot reach here, so SNV output is unchanged.
            m = nt_re.match(nt_mut.upper().replace("U", "T"))
        if not m:
            # nt_re stays the SNV gate and stays FIRST, so every token it accepts takes
            # byte-identical the path it always did. Of the tokens it rejects, only the
            # ones that parse as a genuine multi-base variant are treated differently;
            # an SNV-shaped token it declines still falls through to INVALID_MUTATION.
            #
            # parse_variant is what makes "indel" distinguishable from "garbage";
            # try/except around the legacy parser cannot, because both raise ValueError.
            variant = parse_variant(nt_mut, is_nt=True)
            if variant is not None and not variant.is_snv:
                prow, crow = _non_snv_rows_adabm(pkey, nt_mut, variant, orf_seq, skip_codon,
                                                 aa_ctx=aa_ctx, codon_ctx=codon_ctx)
                if prow is not None:
                    protein_rows.append(prow)
                else:
                    codon_rows.append(crow)
                continue
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
    Build a codon (codon_pos, mut_codon_3letter) -> score lookup.

    codon_msa_encoded_path must be the SINGLE-CHARACTER encoded MSA produced
    by codon_encoding.encode_codon_msa -- same input the params were trained on.
    """
    if not _CODON_ENCODING_AVAILABLE or CODON_ALPHABET is None:
        raise RuntimeError("codon_encoding module not available; cannot score codon side")

    h, J = load_adabmdca_params(codon_params_path, CODON_ALPHABET)
    h, J = apply_zero_sum_gauge(h, J)

    _, focus_seq = _load_focus_sequence(codon_msa_encoded_path, codon_focus_id)
    column_freqs = _column_frequencies(codon_msa_encoded_path, CODON_ALPHABET, encoded=True)

    char_lookup, col_to_seq = build_lookup(h, J, focus_seq, CODON_ALPHABET, column_freqs)

    # Re-key from (seq_pos, mut_char) -> (seq_pos, mut_codon_triplet) for downstream consumption
    out: Dict[Tuple[int, str], Dict] = {}
    for (sp, mut_char), score in char_lookup.items():
        triplet = CHAR_TO_CODON.get(mut_char)
        if triplet is None:
            continue
        out[(sp, triplet)] = score

    # Carry the parameters forward. Previously this returned `out` alone and h/J
    # went out of scope here, which is why nothing spanning more than one codon
    # could ever be scored.
    token_to_idx = {ch: i for i, ch in enumerate(CODON_ALPHABET)}
    wt_idx_per_col = np.array(
        [token_to_idx.get(ch, 0) for ch in focus_seq], dtype=np.int64)
    ctx = ModelContext(h, J, col_to_seq, wt_idx_per_col, column_freqs, CODON_ALPHABET)
    return out, ctx


def _build_protein_lookup_from_params(
    protein_params_path: str,
    protein_msa_path: str,
    protein_focus_id: str,
    alphabet: str = "-ACDEFGHIKLMNPQRSTVWY",
) -> Dict[Tuple[int, str], Dict]:
    """
    Build a protein (aa_pos, mut_aa_symbol) -> score lookup from the protein MSA.
    """
    h, J = load_adabmdca_params(protein_params_path, alphabet)
    h, J = apply_zero_sum_gauge(h, J)

    _, focus_seq = _load_focus_sequence(protein_msa_path, protein_focus_id)
    column_freqs = _column_frequencies(protein_msa_path, alphabet, encoded=False)

    lookup, col_to_seq = build_lookup(h, J, focus_seq, alphabet, column_freqs)
    token_to_idx = {ch: i for i, ch in enumerate(alphabet)}
    wt_idx_per_col = np.array(
        [token_to_idx.get(ch, 0) for ch in focus_seq], dtype=np.int64)
    ctx = ModelContext(h, J, col_to_seq, wt_idx_per_col, column_freqs, alphabet)
    return lookup, ctx


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

    # Intronic gate, before any params are loaded or any row is built. A Potts model
    # has a fixed number of SITES -- residues in the protein model, codons in the codon
    # model -- and an intron occupies neither, so there is no index to evaluate the
    # Hamiltonian at.
    #
    # Unguarded these tokens do not crash: nt_re declines them, parse_variant returns
    # None, and they land in the protein TSV flagged INVALID_MUTATION. That label is
    # false -- 'gd.T5000C' is well formed in a coordinate space this model has no sites
    # for -- and the row's presence implies it was evaluated and found unscorable
    # rather than never being scoreable at all.
    nt_mutations, intronic = split_intronic_tokens(nt_mutations)
    warn_intronic_unsupported(
        'adabmdca', gene, intronic,
        "A Potts model is indexed by residue or codon site; an intron occupies "
        "neither. Score these with RNAfold, miranda, genesplicer or AlphaFold3.")
    # A named row per excluded token, the same contract netNglyc, netphos, netMHC and
    # NetSurfP3 keep. The token is unscorable, not absent: dropping it entirely leaves a
    # hole at adabmDCA for anyone joining the pipelines on pkey, and is
    # indistinguishable from a mutation that was never submitted. Every metric column
    # stays EMPTY -- there is no site index to evaluate.
    intronic_rows = [
        _blank_row(PROTEIN_FIELDNAMES_ADABM, mint_pkey(gene, tok), tok,
                   ["NON_ORF_TOKEN:no_residue_or_codon_site_in_potts_model"])
        for tok in intronic
    ]
    if not nt_mutations:
        print("  (every mutation was intronic)")
        gene_dir = output_dir / gene / "adabmDCA"
        gene_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(intronic_rows, str(gene_dir / f"{gene}.protein.tsv"),
                  PROTEIN_FIELDNAMES_ADABM, extrasaction="ignore")
        if not args.skip_codon:
            write_tsv([], str(gene_dir / f"{gene}.codon.tsv"),
                      CODON_FIELDNAMES_ADABM, extrasaction="ignore")
        return len(intronic_rows), 0

    # Resolve params (training if needed) + the MSA paths the scorer reads.
    protein_params_path, codon_params_path, protein_msa_path, codon_msa_path = \
        _resolve_per_gene_adabm_params(gene, args)

    aa_lookup = None
    aa_ctx = None
    if protein_params_path and protein_msa_path:
        if not args.quiet:
            print(f"  Loading protein adabmDCA params: {protein_params_path}")
        aa_lookup, aa_ctx = _build_protein_lookup_from_params(
            protein_params_path, protein_msa_path, args.focus or gene, args.protein_alphabet
        )

    codon_lookup = None
    codon_ctx = None
    if codon_params_path and codon_msa_path and not args.skip_codon:
        if not args.quiet:
            print(f"  Loading codon adabmDCA params: {codon_params_path}")
        codon_lookup, codon_ctx = _build_codon_lookup_from_params(
            codon_params_path, codon_msa_path, args.codon_focus or "ORF",
        )

    failure_map = load_validation_failures(args.validation_log) if args.validation_log else {}
    protein_rows, codon_rows = score_nt_mutations_adabm(
        nt_mutations, gene, orf_seq,
        aa_lookup=aa_lookup, codon_lookup=codon_lookup,
        failure_map=failure_map, skip_codon=args.skip_codon,
        aa_ctx=aa_ctx, codon_ctx=codon_ctx,
    )
    # Excluded tokens keep their row so a pkey join across pipelines has no hole.
    protein_rows = intronic_rows + protein_rows

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
    """Per-gene file from a direct path or a directory, either layout.

    Delegates to lib/utility.find_gene_file so this and evmutation_pipeline
    resolve the same artifact the same way. The flat glob that stood here was
    non-recursive and returned None at a variant_mapping / MSA-generator root.
    """
    return find_gene_file(path_arg, gene, glob_patterns)


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
    # pseudoDCA uses pseudolikelihood maximization -- different algorithm,
    # different memory footprint (no MCMC chains).
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    # Resolve the backend-aware epoch ceiling. Boltzmann self-terminates at
    # target_pearson so a large ceiling is harmless there; pseudoDCA runs every
    # epoch it is given, so inheriting 50000 was a 100x multiplier on the one path
    # that pays for it. None => per-backend default; an explicit value always wins.
    nepochs = args.adabmdca_nepochs
    if nepochs is None:
        nepochs = 500 if args.adabmdca_model == "pseudoDCA" else 50000
        if not args.quiet:
            print(f"  --adabmdca-nepochs not set; using {nepochs} for "
                  f"{args.adabmdca_model}")
    if args.adabmdca_model == "pseudoDCA":
        train_adabmdca_pseudolikelihood_in_process(
            msa_path=effective_msa,
            alphabet=alphabet,
            output_params_path=output_path,
            nepochs=nepochs,
            lr=args.adabmdca_lr,
            device_str=args.adabmdca_device,
            dtype_str=args.adabmdca_dtype,
            seed=args.adabmdca_seed,
            quiet=args.quiet,
            tol=args.adabmdca_tol,
            patience=args.adabmdca_patience,
            check_every=args.adabmdca_check_every,
        )
    else:
        train_adabmdca_in_process(
            msa_path=effective_msa,
            alphabet=alphabet,
            output_params_path=output_path,
            model=args.adabmdca_model,
            nepochs=nepochs,
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
    parser.add_argument("-f", "--fasta", required=True,
                        help="ORF FASTA file (single gene) or directory")
    parser.add_argument("-m", "--mutations", required=False,
                        help="Mutations CSV file or directory")
    # MSAs are the primary inputs: training is invoked here if params don't exist.
    parser.add_argument("--msa",
                        help="Protein MSA file or directory; triggers `adabmDCA train` if params absent")
    parser.add_argument("-cm", "--codon-msa",
                        help="Codon MSA file or directory (triplets); encoded + trained if params absent")
    # Pre-built params (cache / override) -- skip train if these resolve.
    parser.add_argument("-pp", "--protein-params",
                        help="Pre-built protein params file or directory")
    parser.add_argument("-cp", "--codon-params",
                        help="Pre-built codon params file or directory")
    parser.add_argument("-pa", "--protein-alphabet", default="-ACDEFGHIKLMNPQRSTVWY",
                        help="Protein alphabet (default: standard 20 AA + gap)")
    parser.add_argument("--focus",
                        help="Focus sequence ID in protein MSA (default: gene name)")
    parser.add_argument("-cf", "--codon-focus",
                        help="Focus sequence ID in codon MSA (default: ORF)")
    parser.add_argument("-g", "--gene", help="Gene name override (single-gene mode)")
    parser.add_argument("-sc", "--skip-codon", action="store_true",
                        help="Route synonymous + stop to the protein TSV; do not produce codon TSV")
    parser.add_argument("-st", "--skip-train", action="store_true",
                        help="Skip adabmDCA train invocation (encoding still runs for codon MSAs)")
    # adabmDCA training tunables (consumed by train_adabmdca_in_process when training fires)
    parser.add_argument("-am", "--adabmdca-model", default="bmDCA",
                        choices=["bmDCA", "eaDCA", "edDCA", "pseudoDCA"],
                        help="Training algorithm. bmDCA/eaDCA/edDCA are Boltzmann-learning "
                             "variants (high memory); pseudoDCA is pseudolikelihood maximization "
                             "(~2x less peak memory, no MCMC).")
    # Backend-aware: None resolves in _build_adabmdca_params to 50000 for the
    # Boltzmann routines and 500 for pseudoDCA. A single default cannot serve both
    # -- 50000 is Boltzmann-sized, and the pseudoDCA function's own signature says
    # 500 (plmc's default), so the CLI silently multiplied that path by 100x.
    parser.add_argument("-an", "--adabmdca-nepochs", type=int, default=None,
                        help="Max gradient steps / epochs. Default: 500 for pseudoDCA, "
                             "50000 for bmDCA/eaDCA/edDCA. On pseudoDCA this is a "
                             "CEILING -- the run stops early at convergence.")
    parser.add_argument("-at", "--adabmdca-tol", type=float, default=1e-3,
                        help="pseudoDCA convergence threshold on ||grad||/||grad||_0 "
                             "(default: 1e-3). 0 disables early stopping.")
    parser.add_argument("-ap", "--adabmdca-patience", type=int, default=3,
                        help="Consecutive convergence checks below --adabmdca-tol "
                             "required to stop (default: 3)")
    parser.add_argument("-ace", "--adabmdca-check-every", type=int, default=10,
                        help="Epochs between convergence checks (default: 10)")
    parser.add_argument("-ata", "--adabmdca-target", type=float, default=0.95)
    parser.add_argument("-al", "--adabmdca-lr", type=float, default=0.01)
    parser.add_argument("-anc", "--adabmdca-nchains", type=int, default=10000)
    parser.add_argument("-ans", "--adabmdca-nsweeps", type=int, default=10)
    parser.add_argument("-ad", "--adabmdca-device", default="cuda")
    parser.add_argument("-adt", "--adabmdca-dtype", default="float32",
                        choices=["float32", "float64"])
    parser.add_argument("-as", "--adabmdca-seed", type=int, default=0)
    parser.add_argument("-vl", "--validation-log",
                        help="Validation log from exon_aware_mapping for mutation filtering")
    parser.add_argument("--output", "-o", default=".", help="Output directory")
    parser.add_argument("--quiet", "-q", action="store_true")

    args = parser.parse_args()


    # Directory mode: <root>/<GENE>/mappings/mutations/ sits beside the input,

    # so the root supplies both. Explicit --mutations always wins; FILE MODE and

    # any layout outside the tree are unaffected.

    args.mutations = derive_mutations_root(args.mutations, args.fasta, "adabmdca")

    if not args.mutations:

        parser.error("--mutations is required (no <GENE>/mappings/mutations/ under "

                     f"--fasta {args.fasta})")

    fasta_path = Path(args.fasta)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if fasta_path.is_dir():
        fasta_map = discover_fasta_files(str(fasta_path))
        if not fasta_map:
            print(f"No FASTA files found in {fasta_path}", file=sys.stderr)
            sys.exit(1)

        mutation_index = discover_mutation_files(str(args.mutations)) if args.mutations else {}

        for gene, fasta_file in sorted(fasta_map.items()):
            print(f"\n== {gene} ==")
            try:
                # discover_mutation_files, NOT _find_file_for_gene with "*.csv":
                # a gene tree holds seven CSV types under <GENE>/mappings/ that
                # all extract to the same gene name, so a filename match there is
                # decided by sort order. This selects by directory. The flat
                # fallback keeps FILE MODE and a plain directory working.
                mutations_file = (mutation_index.get(gene)
                                  or _find_file_for_gene(gene, args.mutations, ["*.csv"]))
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
