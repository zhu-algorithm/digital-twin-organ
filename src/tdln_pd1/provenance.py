"""Cohort provenance manifest and cryptographic file verification."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

REQUIRED_MANIFEST_FIELDS = {
    "cohort_id", "source", "license", "modalities", "outcome_definition", "files"
}


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def attach_derivation_audit(
    record: dict, source_path: str | Path, suffix: str, expected_stage: str,
) -> dict:
    """Attach a validated transform sidecar fingerprint when one is present."""
    source = Path(source_path)
    sidecar = Path(f"{source}{suffix}")
    if not sidecar.exists():
        return record
    try:
        payload = json.loads(sidecar.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid derivation audit: {sidecar}") from error
    output = payload.get("output", {}) if isinstance(payload, dict) else {}
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != expected_stage
            or payload.get("labels_used") is not False
            or output.get("file") != source.name
            or output.get("sha256") != sha256_file(source)):
        raise ValueError("invalid or mismatched derivation audit")
    bound = dict(record)
    bound["derivation_audit"] = {
        "file": sidecar.name, "sha256": sha256_file(sidecar),
        "stage": expected_stage,
    }
    return bound


def verify_bound_derivation_audit(record: dict, source_path: str | Path) -> bool:
    """Verify an attached derivation sidecar and its binding to the source file."""
    derivation = record.get("derivation_audit")
    if derivation is None:
        return False
    source = Path(source_path)
    stage = derivation.get("stage") if isinstance(derivation, dict) else None
    suffixes = {
        "label_free_single_cell_featurization": ".feature-audit.json",
        "label_free_spatial_state_featurization": ".spatial-audit.json",
    }
    suffix = suffixes.get(stage)
    sidecar = Path(f"{source}{suffix}") if suffix else None
    if (sidecar is None or derivation.get("file") != sidecar.name
            or not sidecar.is_file() or derivation.get("sha256") != sha256_file(sidecar)):
        raise ValueError("derivation audit SHA-256 mismatch")
    payload = json.loads(sidecar.read_text())
    output = payload.get("output", {}) if isinstance(payload, dict) else {}
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != stage
            or payload.get("labels_used") is not False
            or output.get("file") != source.name
            or output.get("sha256") != sha256_file(source)):
        raise ValueError("derivation audit no longer binds source output")
    return True


def verify_manifest(path: str | Path) -> dict:
    """Validate metadata and verify every declared local file hash."""
    path = Path(path)
    manifest = json.loads(path.read_text())
    missing = REQUIRED_MANIFEST_FIELDS - set(manifest)
    errors: list[str] = []
    if missing:
        errors.append(f"missing manifest fields: {sorted(missing)}")
    for field in ("cohort_id", "source", "license", "outcome_definition"):
        if field in manifest and not str(manifest[field]).strip():
            errors.append(f"manifest field is empty: {field}")
    if not isinstance(manifest.get("modalities", []), list) or not manifest.get("modalities"):
        errors.append("modalities must be a non-empty list")
    files = manifest.get("files", [])
    if not isinstance(files, list) or not files:
        errors.append("files must be a non-empty list")
        files = []
    checked = []
    for item in files:
        if not isinstance(item, dict) or not {"path", "sha256", "role"}.issubset(item):
            errors.append("each file needs path, sha256 and role")
            continue
        target = (path.parent / item["path"]).resolve()
        if not target.is_file():
            errors.append(f"missing file: {item['path']}")
            continue
        observed = sha256_file(target)
        match = observed.lower() == str(item["sha256"]).lower()
        checked.append({"path": item["path"], "role": item["role"],
                        "expected_sha256": item["sha256"], "observed_sha256": observed,
                        "match": match})
        if not match:
            errors.append(f"hash mismatch: {item['path']}")
    return {"valid": not errors, "cohort_id": manifest.get("cohort_id"),
            "checked_files": checked, "errors": errors,
            "manifest": str(path.resolve())}
