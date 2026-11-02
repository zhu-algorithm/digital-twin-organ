"""Integrity manifests for mechanistic sensitivity analyses."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file


def write_sensitivity_audit(
    output_dir: str | Path,
    filename: str,
    stage: str,
    inputs: dict,
    settings: dict,
    generated_filenames: tuple[str, ...],
) -> Path:
    """Bind sensitivity inputs, settings and exact result tables."""
    output = Path(output_dir)
    generated = {}
    for name in generated_filenames:
        path = output / name
        if not path.is_file():
            raise ValueError(f"sensitivity output missing before audit: {name}")
        generated[name] = {"file": name, "sha256": sha256_file(path)}
    payload = {
        "schema_version": 1,
        "stage": stage,
        "research_only": True,
        "inputs": inputs,
        "settings": settings,
        "outputs": generated,
        "limitation": (
            "SHA-256 integrity record only; it does not validate biological assumptions, "
            "cohort transferability, or clinical utility."
        ),
    }
    audit_path = output / filename
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit_path


def verify_sensitivity_audit(audit_path: str | Path) -> dict:
    """Fail closed when a sensitivity result no longer matches its audit."""
    path = Path(audit_path)
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid sensitivity audit: {path}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") not in {
                "systemic_sensitivity", "spatial_sensitivity",
            }):
        raise ValueError("unsupported sensitivity audit schema or stage")
    outputs = payload.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("sensitivity audit has no bound outputs")
    verified = []
    for key, record in outputs.items():
        if (not isinstance(record, dict) or record.get("file") != key
                or not isinstance(record.get("sha256"), str)):
            raise ValueError("invalid sensitivity output record")
        target = path.parent / key
        if not target.is_file() or sha256_file(target) != record["sha256"]:
            raise ValueError(f"sensitivity output SHA-256 mismatch: {key}")
        verified.append(key)
    return {
        "valid": True,
        "stage": payload["stage"],
        "verified_outputs": sorted(verified),
        "audit_sha256": sha256_file(path),
        "payload": payload,
    }
