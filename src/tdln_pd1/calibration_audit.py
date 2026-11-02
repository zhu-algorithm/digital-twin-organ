"""Cryptographic input/output manifests for mechanistic calibration runs."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file

CALIBRATION_STAGES = {"systemic_calibration", "spatial_calibration"}


def write_calibration_audit(
    output_dir: str | Path,
    filename: str,
    stage: str,
    inputs: dict,
    settings: dict,
    generated_filenames: tuple[str, ...],
) -> Path:
    if stage not in CALIBRATION_STAGES:
        raise ValueError(f"unsupported calibration audit stage: {stage}")
    output = Path(output_dir)
    generated = {}
    for name in generated_filenames:
        target = output / name
        if not target.is_file():
            raise ValueError(f"calibration output missing before audit: {name}")
        generated[name] = {"file": name, "sha256": sha256_file(target)}
    payload = {
        "schema_version": 1, "stage": stage, "research_only": True,
        "inputs": inputs, "settings": settings, "outputs": generated,
        "limitation": (
            "Integrity manifest only; it does not establish identifiability, biological "
            "validity, cohort transportability, or clinical utility."
        ),
    }
    path = output / filename
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def verify_calibration_audit(audit_path: str | Path) -> dict:
    path = Path(audit_path)
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid calibration audit: {path}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") not in CALIBRATION_STAGES):
        raise ValueError("unsupported calibration audit schema or stage")
    outputs = payload.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("calibration audit has no bound outputs")
    verified = []
    for name, record in outputs.items():
        if (not isinstance(record, dict) or record.get("file") != name
                or not isinstance(record.get("sha256"), str)):
            raise ValueError("invalid calibration output record")
        target = path.parent / name
        if not target.is_file() or sha256_file(target) != record["sha256"]:
            raise ValueError(f"calibration output SHA-256 mismatch: {name}")
        verified.append(name)
    return {"valid": True, "stage": payload["stage"],
            "verified_outputs": sorted(verified), "audit_sha256": sha256_file(path),
            "payload": payload}
