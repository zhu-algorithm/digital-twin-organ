"""Integrity chain for label-free mechanistic feature augmentation."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file


def write_augmentation_audit(
    output_path: str | Path,
    inputs: dict,
    settings: dict,
) -> Path:
    output = Path(output_path)
    if not output.is_file():
        raise ValueError("augmented feature output must exist before audit")
    payload = {
        "schema_version": 1,
        "stage": "label_free_twin_augmentation",
        "research_only": True,
        "labels_used": False,
        "inputs": inputs,
        "settings": settings,
        "output": {"file": output.name, "sha256": sha256_file(output)},
        "limitation": (
            "Integrity record only; it does not validate assay units, mechanistic "
            "assumptions, treatment effects, or clinical utility."
        ),
    }
    audit = Path(f"{output}.twin-audit.json")
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit


def verify_augmentation_audit(
    audit_path: str | Path,
    output_path: str | Path,
    input_paths: dict[str, str | Path],
) -> dict:
    audit = Path(audit_path)
    try:
        payload = json.loads(audit.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid augmentation audit: {audit}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "label_free_twin_augmentation"
            or payload.get("labels_used") is not False):
        raise ValueError("unsupported or non-label-free augmentation audit")
    output = Path(output_path)
    recorded_output = payload.get("output", {})
    if (recorded_output.get("file") != output.name
            or not output.is_file()
            or recorded_output.get("sha256") != sha256_file(output)):
        raise ValueError("augmented feature output SHA-256 mismatch")
    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("augmentation audit has no input records")
    verified_inputs = []
    for name, record in inputs.items():
        if record is None:
            continue
        if (not isinstance(record, dict) or not isinstance(record.get("file"), str)
                or not isinstance(record.get("sha256"), str)):
            raise ValueError(f"invalid augmentation input record: {name}")
        if name not in input_paths:
            raise ValueError(f"augmentation verification requires input path: {name}")
        target = Path(input_paths[name])
        if (not target.is_file() or target.name != record["file"]
                or sha256_file(target) != record["sha256"]):
            raise ValueError(f"augmentation input SHA-256 mismatch: {name}")
        derivation = record.get("derivation_audit")
        if derivation is not None:
            suffixes = {
                "features": ".feature-audit.json",
                "spatial_initial_zones": ".spatial-audit.json",
            }
            suffix = suffixes.get(name)
            derived_audit = Path(f"{target}{suffix}") if suffix else None
            if (derived_audit is None or not isinstance(derivation, dict)
                    or derivation.get("file") != derived_audit.name
                    or not derived_audit.is_file()
                    or derivation.get("sha256") != sha256_file(derived_audit)):
                raise ValueError(f"{name} derivation audit SHA-256 mismatch")
        verified_inputs.append(name)
    return {
        "valid": True,
        "verified_inputs": sorted(verified_inputs),
        "output_sha256": recorded_output["sha256"],
        "audit_sha256": sha256_file(audit),
        "payload": payload,
    }
