"""Source-derived RAM admission estimates for the existing EVmutation workflow.

For plmc's gap-reduced L-BFGS path, P = L*q + L*(L-1)*q*q/2 parameters.
The native peak retains x, five optimizer work vectors, twelve history vectors
(history length six), and the marginal frequencies: nineteen P-sized arrays.
Each OpenMP worker additionally retains two L*q*q site blocks and two q vectors.
Native arithmetic is conservatively assumed double precision; this does not
change the executable's compilation or optimizer settings.

EVmutation scoring retains four full float64 L*L*q*q arrays: pair frequencies
and couplings in the loaded model and its independent-model deepcopy. Admission
uses the larger training/scoring phase, plus input/lookup/runtime allowances and
a safety margin. These uncalibrated estimates are not measured RSS guarantees.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct

try:
    from .codon_encoding import CODON_ALPHABET, CODON_TO_CHAR
    from .resource_planner import GIB, blocked_plan, positive, positive_integer
except ImportError:
    from codon_encoding import CODON_ALPHABET, CODON_TO_CHAR
    from resource_planner import GIB, blocked_plan, positive, positive_integer


PROTEIN_ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"
NATIVE_SCALAR_BYTES = 8
LBFGS_HISTORY = 6
NATIVE_INT_MAX = 2 ** 31 - 1
RUNTIME_ALLOWANCE_BYTES = GIB
LOOKUP_BYTES_PER_STATE = 4096


def _fasta_records(filename):
    """Read one alignment record at a time without retaining the whole MSA."""
    header = None
    pieces = []
    with open(filename) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(pieces)
                header = line[1:]
                pieces = []
            else:
                if header is None:
                    raise ValueError(f"Missing FASTA header in {filename}")
                pieces.append(line.strip())
    if header is not None:
        yield header, "".join(pieces)


def alignment_shape(msa, side, focus):
    """Match BFF's explicit custom alphabet, -g, and plmc prefix-focus rules.

    A missing focus leaves all alignment columns modeled, including gaps.
    Explicit -a uses plmc's custom-alphabet branch, not its default protein
    lowercase-insertion branch. Codon triplets use the shared BFF encoding.
    Raw sequence counts are retained conservatively for input storage even
    when plmc excludes invalid rows or codon encoding deduplicates headers.
    """
    if side not in {"protein", "codon"}:
        raise ValueError("EVmutation side must be protein or codon")
    alphabet = CODON_ALPHABET if side == "codon" else PROTEIN_ALPHABET
    sequences = 0
    valid_sequences = 0
    width = None
    focus_header = None
    focus_sites = None
    for header, sequence in _fasta_records(msa):
        if side == "codon":
            sequence = sequence.replace(" ", "").upper()
            if len(sequence) % 3:
                raise ValueError(f"Codon alignment width is not divisible by three: {msa}")
            sequence = "".join(CODON_TO_CHAR.get(sequence[position:position + 3], "-")
                               for position in range(0, len(sequence), 3))
        if not sequence:
            raise ValueError(f"Empty alignment record in {msa}")
        if width is None:
            width = len(sequence)
        elif len(sequence) != width:
            raise ValueError(f"Unequal alignment widths in {msa}")
        sequences += 1
        valid_sequences += all(symbol in alphabet for symbol in sequence)
        if focus_header is None and header.startswith(focus):
            focus_header = header
            focus_sites = sum(symbol != "-" for symbol in sequence)
        elif side == "codon" and header == focus_header:
            focus_sites = sum(symbol != "-" for symbol in sequence)
    if not sequences or not valid_sequences:
        raise ValueError(f"Empty alignment or no valid custom-alphabet sequences: {msa}")
    sites = width if focus_sites is None else focus_sites
    if not sites:
        raise ValueError(f"Focus has no modeled non-gap columns: {msa}")
    return {
        "sequences": sequences, "valid_sequences": valid_sequences,
        "alignment_sites": width, "sites": sites, "states": len(alphabet) - 1,
        "focus": focus, "focus_found": focus_header is not None,
        "focus_header": focus_header,
        "focus_mode": "prefix_non_gap" if focus_header is not None else "missing_focus_full_alignment",
    }


def params_shape(params):
    """Read only the native v2 header and reject missing/truncated model files.

    BFF's EVmutation loader uses the native plmc_v2 format with float32 file
    payloads, regardless of the training executable's arithmetic precision.
    Its dense pair arrays are subsequently allocated as float64.
    """
    filename = Path(params)
    with filename.open("rb") as handle:
        header = handle.read(20)
    if len(header) != 20:
        raise ValueError(f"Invalid or truncated plmc_v2 params header: {params}")
    sites, states, valid, invalid, iterations = struct.unpack("=5i", header)
    if sites < 1 or states < 1 or states > 256 or min(valid, invalid, iterations) < 0:
        raise ValueError(f"Invalid plmc_v2 params dimensions: {params}")
    pairs = sites * (sites - 1) // 2
    required_bytes = 40 + states + 4 * (valid + invalid) + 5 * sites
    required_bytes += 8 * sites * states + 8 * pairs * states * states
    if filename.stat().st_size < required_bytes:
        raise ValueError(f"Truncated plmc_v2 params payload: {params}")
    return {"sites": sites, "states": states, "sequences": valid + invalid}


def estimate_memory(sequences, sites, states, workspace_threads=1, alignment_sites=None,
                    msa_bytes=0, prebuilt=False, margin=1.15):
    """Count the native optimizer and Python scoring allocations without training."""
    sequences = positive_integer(sequences, "sequence count")
    sites = positive_integer(sites, "modeled sites")
    states = positive_integer(states, "model states")
    workspace_threads = min(sites, positive_integer(workspace_threads, "workspace threads"))
    alignment_sites = sites if alignment_sites is None else positive_integer(alignment_sites, "alignment sites")
    margin = positive(margin, "memory margin")
    if margin < 1:
        raise ValueError("Memory margin must be at least one")
    if msa_bytes < 0:
        raise ValueError("MSA byte count must be nonnegative")
    pairs = sites * (sites - 1) // 2
    parameters = sites * states + pairs * states * states
    if not prebuilt and parameters > NATIVE_INT_MAX:
        raise ValueError(
            f"Unsupported plmc dimensions L={sites}, q={states}: {parameters} parameters "
            "exceed the native signed 32-bit parameter count"
        )
    if not prebuilt and sequences * alignment_sites > NATIVE_INT_MAX:
        raise ValueError("Unsupported plmc alignment: native signed 32-bit sequence indexing overflows")
    input_bytes = 12 * sequences * alignment_sites + 4 * msa_bytes
    input_bytes += 256 * sequences + 8 * (sequences + alignment_sites)
    parameter_bytes = parameters * NATIVE_SCALAR_BYTES
    optimizer_bytes = (1 + 5 + 2 * LBFGS_HISTORY) * parameter_bytes
    marginal_bytes = parameter_bytes
    workspace_bytes = workspace_threads * (2 * sites * states * states + 2 * states) * NATIVE_SCALAR_BYTES
    regularization_bytes = (sites + pairs) * NATIVE_SCALAR_BYTES
    training_bytes = optimizer_bytes + marginal_bytes + workspace_bytes + regularization_bytes + input_bytes
    scoring_dense_bytes = 4 * sites * sites * states * states * 8
    scoring_bytes = scoring_dense_bytes + LOOKUP_BYTES_PER_STATE * sites * states + input_bytes
    training_gib = 0.0 if prebuilt else (training_bytes + RUNTIME_ALLOWANCE_BYTES) / GIB * margin
    scoring_gib = (scoring_bytes + RUNTIME_ALLOWANCE_BYTES) / GIB * margin
    return {
        "estimated_memory_gib": max(training_gib, scoring_gib),
        "training_memory_gib": training_gib, "scoring_memory_gib": scoring_gib,
        "parameter_count": parameters, "native_scalar_bytes": NATIVE_SCALAR_BYTES,
        "assumed_native_precision": "float64", "scoring_scalar_bytes": 8,
        "lbfgs_history": LBFGS_HISTORY, "workspace_threads": workspace_threads,
        "optimizer_bytes": 0 if prebuilt else optimizer_bytes,
        "marginal_bytes": 0 if prebuilt else marginal_bytes,
        "thread_workspace_bytes": 0 if prebuilt else workspace_bytes,
        "scoring_dense_bytes": scoring_dense_bytes,
        "input_allowance_bytes": input_bytes,
        "runtime_allowance_bytes": RUNTIME_ALLOWANCE_BYTES,
    }


def plan_evmutation_task(gene, side, msa, config, params=None):
    """Plan one EV side; explicit RAM is a floor, not an estimator bypass.

    Reserve thread-local workspace for the complete usable CPU allocation when
    available, avoiding circular memory/thread-share calculations. The actual
    task still uses config['threads'], which can be a smaller automatic share.
    """
    focus = "ORF" if side == "codon" else gene
    shape = alignment_shape(msa, side, focus)
    scoring_shape = params_shape(params) if params else shape
    hardware = config.get("hardware", {})
    threads = positive_integer(config.get("threads", 1), "threads")
    workspace_threads = positive_integer(hardware.get("cpus", threads), "CPU budget")
    if threads > workspace_threads:
        raise ValueError(f"{gene}/{side}: threads exceed the CPU budget")
    estimates = estimate_memory(
        max(shape["sequences"], scoring_shape["sequences"]), scoring_shape["sites"],
        scoring_shape["states"], workspace_threads=workspace_threads,
        alignment_sites=shape["alignment_sites"], msa_bytes=Path(msa).stat().st_size,
        prebuilt=bool(params), margin=config.get("memory_margin", 1.15),
    )
    memory_floor = config.get("evmutation", {}).get("memory_gib")
    memory_gib = estimates["estimated_memory_gib"]
    if memory_floor is not None:
        memory_gib = max(memory_gib, positive(memory_floor, "EVmutation RAM minimum"))
    if "memory_gib" in hardware and memory_gib > positive(hardware["memory_gib"], "RAM budget"):
        raise ValueError(
            f"{gene}/{side}: EVmutation unschedulable; RAM request {memory_gib:.2f} GiB "
            f"exceeds RAM budget {hardware['memory_gib']:.2f} GiB"
        )
    return {
        "gene": gene, "side": side, "backend": "evmutation", "device": "cpu",
        **shape, **estimates,
        "sites": scoring_shape["sites"], "states": scoring_shape["states"],
        "threads": threads, "memory_gib": memory_gib, "memory_floor_gib": memory_floor,
        "prebuilt": bool(params), "params": str(Path(params).resolve()) if params else None,
        "estimate_basis": "uncalibrated native L-BFGS and dense EV scoring allocation envelope",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan"])
    for name in ("gene", "side", "msa", "output"):
        parser.add_argument(f"--{name}", required=True, choices=("protein", "codon") if name == "side" else None)
    parser.add_argument("--params")
    parser.add_argument("--config")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--memory-gib", type=float)
    parser.add_argument("--defer-errors", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.threads is not None:
            positive_integer(args.threads, "threads")
        config = json.loads(Path(args.config).read_text()) if args.config else {
            "threads": args.threads if args.threads is not None else 1,
            "memory_margin": 1.15, "evmutation": {"memory_gib": args.memory_gib},
        }
        if not isinstance(config, dict):
            raise ValueError("Resource config must be a mapping")
        if not isinstance(config.get("evmutation", {}), dict) or not isinstance(config.get("hardware", {}), dict):
            raise ValueError("EVmutation settings and hardware must be mappings")
        if args.memory_gib is not None:
            cli_floor = positive(args.memory_gib, "EVmutation RAM minimum")
            evmutation = config.setdefault("evmutation", {})
            config_floor = evmutation.get("memory_gib")
            evmutation["memory_gib"] = max(cli_floor, positive(config_floor, "EVmutation RAM minimum")) if config_floor is not None else cli_floor
        positive_integer(config.get("threads", 1), "threads")
        if positive(config.get("memory_margin", 1.15), "memory margin") < 1:
            raise ValueError("Memory margin must be at least one")
        memory_floor = config.get("evmutation", {}).get("memory_gib")
        if memory_floor is not None:
            positive(memory_floor, "EVmutation RAM minimum")
        hardware = config.get("hardware", {})
        if "cpus" in hardware:
            positive_integer(hardware["cpus"], "CPU budget")
        if "memory_gib" in hardware:
            positive(hardware["memory_gib"], "RAM budget")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"EVmutation resource planning failed: {error}\n")
    try:
        plan = plan_evmutation_task(args.gene, args.side, args.msa, config, args.params)
    except (OSError, ValueError, KeyError, TypeError) as error:
        if not args.defer_errors:
            parser.exit(1, f"EVmutation resource planning failed: {error}\n")
        plan = blocked_plan(args.gene, args.side, "evmutation", error)
    try:
        Path(args.output).write_text(json.dumps(plan, indent=2) + "\n")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"EVmutation resource planning failed: {error}\n")
    if "resource_error" not in plan:
        print(f"[resources] {args.gene}/{args.side} EVmutation: CPU; RAM {plan['memory_gib']:.2f} GiB, L={plan['sites']}, q={plan['states']}")


if __name__ == "__main__":
    main()
