"""Verify EVmutation artifacts against routing settings and input contents."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

from biofeaturefactory.mutation_effects.bin.resource_planner import file_digest


SCHEMA_VERSION = 1


def routing_fingerprint(gene, side, fasta, mutations, validation_log=None,
                        skip_codon=False, score_missense_codon=False):
    """Hash scoring inputs and effective routing independently of file locations."""
    backend_directory = Path(__file__).resolve().parent.parent
    payload = {
        "schema": SCHEMA_VERSION,
        "gene": gene,
        "side": side,
        "inputs": {
            "fasta": file_digest(fasta),
            "mutations": file_digest(mutations),
            "validation_log": file_digest(validation_log) if validation_log else None,
        },
        "settings": {
            "skip_codon": bool(skip_codon),
            "score_missense_codon": bool(score_missense_codon),
        },
        "implementation": {
            filename: file_digest(backend_directory / filename)
            for filename in ("evmutation_pipeline.py", "bin/codon_encoding.py")
        },
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _file_record(filename):
    path = Path(filename)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty EVmutation artifact: {path}")
    return {"size_bytes": path.stat().st_size, "sha256": file_digest(path)}


def verify_completion(marker, fingerprint, tsv, params, msa=None):
    """Return False for old, incomplete, changed, or unreadable artifacts."""
    try:
        if not isinstance(fingerprint, str) or not fingerprint.strip():
            return False
        manifest = json.loads(Path(marker).read_text())
        if manifest["schema"] != SCHEMA_VERSION or manifest["success"] is not True:
            return False
        if manifest["fingerprint"] != fingerprint:
            return False
        artifacts = {"tsv": tsv, "params": params}
        if msa is not None:
            artifacts["msa"] = msa
        if set(manifest["files"]) != set(artifacts):
            return False
        return all(manifest["files"][role] == _file_record(path) for role, path in artifacts.items())
    except (OSError, ValueError, KeyError, TypeError):
        return False


def write_completion(marker, fingerprint, tsv, params, msa=None):
    """Write the completion marker atomically after hashing complete artifacts."""
    if not isinstance(fingerprint, str) or not fingerprint.strip():
        raise ValueError("A nonempty routing fingerprint is required")
    artifacts = {"tsv": tsv, "params": params}
    if msa is not None:
        artifacts["msa"] = msa
    manifest = {
        "schema": SCHEMA_VERSION,
        "success": True,
        "fingerprint": fingerprint,
        "files": {role: _file_record(path) for role, path in artifacts.items()},
    }
    marker = Path(marker)
    marker.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{marker.name}.", suffix=".tmp", dir=marker.parent)
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(manifest, handle, sort_keys=True, indent=2)
            handle.write("\n")
        os.replace(temporary, marker)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["write"])
    for name in ("marker", "fingerprint", "tsv", "params"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--msa")
    arguments = parser.parse_args()
    try:
        write_completion(arguments.marker, arguments.fingerprint, arguments.tsv,
                         arguments.params, arguments.msa)
    except (OSError, ValueError, TypeError) as error:
        parser.exit(1, f"EVmutation completion failed: {error}\n")


if __name__ == "__main__":
    main()
