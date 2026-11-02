"""Integrity chain for label-free single-cell patient featurization."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file


def write_feature_audit(
    cells_path: str | Path,
    output_path: str | Path,
    settings: dict,
) -> Path:
    cells, output = Path(cells_path), Path(output_path)
    if not cells.is_file() or not output.is_file():
        raise ValueError("cell input and feature output must exist before audit")
    payload = {
        "schema_version": 1,
        "stage": "label_free_single_cell_featurization",
        "research_only": True,
        "labels_used": False,
        "input": {"file": cells.name, "sha256": sha256_file(cells)},
        "settings": settings,
        "output": {"file": output.name, "sha256": sha256_file(output)},
        "limitation": (
            "Integrity record only; it does not establish annotation quality, "
            "assay comparability, biological validity, or clinical utility."
        ),
    }
    audit = Path(f"{output}.feature-audit.json")
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit


def verify_feature_audit(
    audit_path: str | Path,
    cells_path: str | Path,
    output_path: str | Path,
) -> dict:
    audit, cells, output = Path(audit_path), Path(cells_path), Path(output_path)
    try:
        payload = json.loads(audit.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid feature audit: {audit}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "label_free_single_cell_featurization"
            or payload.get("labels_used") is not False):
        raise ValueError("unsupported or non-label-free feature audit")
    for name, path, record in (
        ("cell input", cells, payload.get("input", {})),
        ("feature output", output, payload.get("output", {})),
    ):
        if (not path.is_file() or record.get("file") != path.name
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"{name} SHA-256 mismatch")
    return {
        "valid": True,
        "input_sha256": payload["input"]["sha256"],
        "output_sha256": payload["output"]["sha256"],
        "audit_sha256": sha256_file(audit),
        "payload": payload,
    }
