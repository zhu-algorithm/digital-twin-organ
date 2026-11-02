"""Integrity manifest for calibrated mechanistic uncertainty propagation."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file


def write_uncertainty_audit(output_dir, inputs: dict, settings: dict) -> Path:
    output = Path(output_dir)
    outputs = {}
    for name in ("uncertainty_samples.csv", "uncertainty_summary.json"):
        path = output / name
        if not path.is_file():
            raise ValueError(f"uncertainty output missing before audit: {name}")
        outputs[name] = {"file": name, "sha256": sha256_file(path)}
    payload = {
        "schema_version": 1,
        "stage": "systemic_calibration_uncertainty_propagation",
        "research_only": True,
        "inputs": inputs,
        "settings": settings,
        "outputs": outputs,
        "limitation": (
            "Local asymptotic parameter uncertainty only; excludes structural, "
            "measurement-model, cohort-shift and clinical uncertainty."
        ),
    }
    audit = output / "uncertainty_provenance.json"
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit


def verify_uncertainty_audit(audit_path, input_paths: dict) -> dict:
    audit = Path(audit_path)
    try:
        payload = json.loads(audit.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid uncertainty audit: {audit}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "systemic_calibration_uncertainty_propagation"):
        raise ValueError("unsupported uncertainty audit")
    verified_inputs = []
    for name, record in payload.get("inputs", {}).items():
        if name not in input_paths:
            raise ValueError(f"uncertainty verification requires input path: {name}")
        target = Path(input_paths[name])
        if (not target.is_file() or target.name != record.get("file")
                or sha256_file(target) != record.get("sha256")):
            raise ValueError(f"uncertainty input SHA-256 mismatch: {name}")
        verified_inputs.append(name)
    verified_outputs = []
    for name, record in payload.get("outputs", {}).items():
        target = audit.parent / name
        if (record.get("file") != name or not target.is_file()
                or sha256_file(target) != record.get("sha256")):
            raise ValueError(f"uncertainty output SHA-256 mismatch: {name}")
        verified_outputs.append(name)
    if not verified_inputs or not verified_outputs:
        raise ValueError("uncertainty audit must bind inputs and outputs")
    return {
        "valid": True,
        "verified_inputs": sorted(verified_inputs),
        "verified_outputs": sorted(verified_outputs),
        "audit_sha256": sha256_file(audit),
        "payload": payload,
    }
