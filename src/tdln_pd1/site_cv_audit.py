"""Replayable chain of custody for leave-one-site-out validation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .augmentation_audit import verify_augmentation_audit
from .model import leave_one_site_out
from .provenance import sha256_file


AUDIT_FILENAME = "site_cv_provenance.json"
METRICS_FILENAME = "site_cv_metrics.json"
PREDICTIONS_FILENAME = "site_cv_predictions.csv"


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_site_cv_audit(
    features_path: str | Path,
    labels_path: str | Path,
    output: str | Path,
    *,
    site_column: str,
    bootstrap: int,
    augmentation_input_paths: dict[str, str | Path] | None = None,
) -> Path:
    """Bind exact cohort inputs, settings and both site-CV result files."""
    features_path, labels_path, output = Path(features_path), Path(labels_path), Path(output)
    metrics_path = output / METRICS_FILENAME
    predictions_path = output / PREDICTIONS_FILENAME
    twin_audit_path = Path(f"{features_path}.twin-audit.json")
    twin_audit_record = None
    if twin_audit_path.is_file():
        replay = verify_augmentation_audit(
            twin_audit_path, features_path, augmentation_input_paths or {},
        )
        twin_audit_record = {
            "file": twin_audit_path.name,
            "sha256": sha256_file(twin_audit_path),
            "derivation_replayed": True,
            "verified_inputs": replay["verified_inputs"],
            "nested_derivations_replayed": replay["nested_derivations_replayed"],
        }
    payload = {
        "schema_version": 1,
        "stage": "leave_one_site_out_validation",
        "research_only": True,
        "settings": {"site_column": site_column, "bootstrap_draws": int(bootstrap)},
        "inputs": {
            "features": {"file": features_path.name, "sha256": sha256_file(features_path)},
            "labels": {"file": labels_path.name, "sha256": sha256_file(labels_path)},
            "augmentation_audit": twin_audit_record,
        },
        "outputs": {
            "metrics": {"file": metrics_path.name, "sha256": sha256_file(metrics_path)},
            "predictions": {
                "file": predictions_path.name,
                "sha256": sha256_file(predictions_path),
            },
        },
        "limitation": (
            "Deterministic integrity and internal centre-transportability evidence only; "
            "does not establish independent external or clinical validation."
        ),
    }
    audit_path = output / AUDIT_FILENAME
    _write_json(audit_path, payload)
    return audit_path


def verify_site_cv_audit(
    audit_path: str | Path,
    features_path: str | Path,
    labels_path: str | Path,
    metrics_path: str | Path,
    predictions_path: str | Path,
    augmentation_input_paths: dict[str, str | Path] | None = None,
) -> dict:
    """Re-run leave-one-site-out fitting and require exact bound results."""
    audit_path = Path(audit_path)
    features_path, labels_path = Path(features_path), Path(labels_path)
    metrics_path, predictions_path = Path(metrics_path), Path(predictions_path)
    try:
        payload = json.loads(audit_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("invalid site-CV audit") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "leave_one_site_out_validation"):
        raise ValueError("unsupported site-CV audit schema")
    settings = payload.get("settings", {})
    site_column, bootstrap = settings.get("site_column"), settings.get("bootstrap_draws")
    if (not isinstance(site_column, str) or not site_column.strip()
            or not isinstance(bootstrap, int) or isinstance(bootstrap, bool) or bootstrap < 0):
        raise ValueError("invalid site-CV audit settings")
    for name, path in (("features", features_path), ("labels", labels_path)):
        record = payload.get("inputs", {}).get(name, {})
        if (record.get("file") != path.name or not path.is_file()
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"site-CV {name} SHA-256 mismatch")
    twin_record = payload.get("inputs", {}).get("augmentation_audit")
    twin_audit_path = Path(f"{features_path}.twin-audit.json")
    feature_columns = pd.read_csv(features_path, nrows=0).columns
    mechanistic = any(
        column.startswith(("twin_", "liver_cf_")) for column in feature_columns
    )
    if twin_record is not None:
        if (not isinstance(twin_record, dict)
                or twin_record.get("file") != twin_audit_path.name
                or not twin_audit_path.is_file()
                or twin_record.get("sha256") != sha256_file(twin_audit_path)):
            raise ValueError("site-CV augmentation audit SHA-256 mismatch")
        replay = verify_augmentation_audit(
            twin_audit_path, features_path, augmentation_input_paths or {},
        )
        if (twin_record.get("derivation_replayed") is not True
                or twin_record.get("verified_inputs") != replay["verified_inputs"]
                or twin_record.get("nested_derivations_replayed")
                != replay["nested_derivations_replayed"]):
            raise ValueError("site-CV augmentation replay contract mismatch")
    elif mechanistic:
        raise ValueError("mechanistic site-CV replay requires a bound augmentation audit")
    for name, path in (("metrics", metrics_path), ("predictions", predictions_path)):
        record = payload.get("outputs", {}).get(name, {})
        if (record.get("file") != path.name or not path.is_file()
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"site-CV {name} SHA-256 mismatch")
    replayed_metrics, replayed_predictions = leave_one_site_out(
        pd.read_csv(features_path), pd.read_csv(labels_path), site_column, bootstrap,
    )
    if replayed_metrics != json.loads(metrics_path.read_text()):
        raise ValueError("site-CV metric derivation replay mismatch")
    expected_predictions = replayed_predictions.to_csv(index=False)
    if expected_predictions != predictions_path.read_text():
        raise ValueError("site-CV prediction derivation replay mismatch")
    return {
        "valid": True,
        "stage": payload["stage"],
        "site_column": site_column,
        "bootstrap_draws": bootstrap,
        "patients_replayed": int(len(replayed_predictions)),
        "augmentation_audit_verified": twin_record is not None,
        "augmentation_derivation_replayed": twin_record is not None,
        "derivation_replayed": True,
        "audit_sha256": sha256_file(audit_path),
    }
