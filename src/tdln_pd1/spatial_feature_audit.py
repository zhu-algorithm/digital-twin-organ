"""Integrity chain for spot-to-zone TDLN state aggregation."""

from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256_file


def write_spatial_feature_audit(
    spots_path: str | Path,
    output_path: str | Path,
    weight_column: str | None,
) -> Path:
    spots, output = Path(spots_path), Path(output_path)
    if not spots.is_file() or not output.is_file():
        raise ValueError("spatial spot input and state output must exist before audit")
    payload = {
        "schema_version": 1,
        "stage": "label_free_spatial_state_featurization",
        "research_only": True,
        "labels_used": False,
        "input": {"file": spots.name, "sha256": sha256_file(spots)},
        "settings": {
            "weight_column": weight_column,
            "anatomy_policy": "curated_zone_labels_no_coordinate_inference",
            "missing_zone_policy": "explicit_zero_with_zone_observed_false",
            "labels_used": False,
        },
        "output": {"file": output.name, "sha256": sha256_file(output)},
        "limitation": (
            "Integrity record only; it does not validate segmentation, anatomy, "
            "deconvolution, abundance units, or clinical utility."
        ),
    }
    audit = Path(f"{output}.spatial-audit.json")
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit


def verify_spatial_feature_audit(
    audit_path: str | Path,
    spots_path: str | Path,
    output_path: str | Path,
) -> dict:
    audit, spots, output = Path(audit_path), Path(spots_path), Path(output_path)
    try:
        payload = json.loads(audit.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid spatial feature audit: {audit}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "label_free_spatial_state_featurization"
            or payload.get("labels_used") is not False):
        raise ValueError("unsupported or non-label-free spatial feature audit")
    for name, path, record in (
        ("spatial spot input", spots, payload.get("input", {})),
        ("spatial state output", output, payload.get("output", {})),
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
