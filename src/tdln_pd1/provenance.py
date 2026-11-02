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
