"""Integrity chain for spot-to-zone TDLN state aggregation."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
import pandas as pd

from .provenance import sha256_file
from .spatial_data import aggregate_spatial_spots


def write_spatial_feature_audit(
    spots_path: str | Path,
    output_path: str | Path,
    weight_column: str | None,
    time_column: str | None = None,
    baseline_values: tuple[str, ...] | None = None,
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
            "time_column": time_column,
            "baseline_values": list(baseline_values or ()),
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
    settings = payload.get("settings")
    if (not isinstance(settings, dict)
            or settings.get("anatomy_policy") != "curated_zone_labels_no_coordinate_inference"
            or settings.get("missing_zone_policy") != "explicit_zero_with_zone_observed_false"
            or settings.get("labels_used") is not False
            or not isinstance(settings.get("baseline_values"), list)):
        raise ValueError("invalid spatial feature derivation settings")
    for name, path, record in (
        ("spatial spot input", spots, payload.get("input", {})),
        ("spatial state output", output, payload.get("output", {})),
    ):
        if (not path.is_file() or record.get("file") != path.name
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"{name} SHA-256 mismatch")
    replayed = aggregate_spatial_spots(
        pd.read_csv(spots), settings.get("weight_column"), settings.get("time_column"),
        tuple(str(value) for value in settings["baseline_values"]) or None,
    )
    replayed_sha = hashlib.sha256(
        replayed.to_csv(index=False).encode("utf-8")
    ).hexdigest()
    if replayed_sha != payload["output"]["sha256"]:
        raise ValueError("spatial feature derivation replay mismatch")
    return {
        "valid": True,
        "input_sha256": payload["input"]["sha256"],
        "output_sha256": payload["output"]["sha256"],
        "audit_sha256": sha256_file(audit),
        "derivation_replayed": True,
        "payload": payload,
    }
