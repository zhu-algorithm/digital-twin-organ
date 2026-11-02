"""Cryptographic chain of custody for locked external validation files."""

from __future__ import annotations

import json
from pathlib import Path

import joblib

from .provenance import sha256_file


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _load_json(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid validation audit: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"validation audit must contain one JSON object: {path}")
    return payload


def write_prediction_audit(
    model_path: str | Path,
    features_path: str | Path,
    predictions_path: str | Path,
    patient_count: int,
) -> Path:
    """Bind one locked artifact and feature table to its prediction file."""
    model_path, features_path = Path(model_path), Path(features_path)
    predictions_path = Path(predictions_path)
    artifact = joblib.load(model_path)
    required = {"artifact_version", "training_cohort_fingerprint",
                "training_feature_cohort_fingerprint", "point_prediction_method"}
    if not isinstance(artifact, dict) or not required.issubset(artifact):
        raise ValueError("locked model lacks versioned training fingerprints")
    if (artifact["artifact_version"] != 8
            or artifact["point_prediction_method"] != "full_cohort_locked_model"):
        raise ValueError("locked model lacks the version 8 point-prediction contract")
    payload = {
        "schema_version": 2,
        "stage": "locked_prediction",
        "research_only": True,
        "model_refit": False,
        "patient_count": int(patient_count),
        "artifact_version": int(artifact["artifact_version"]),
        "point_prediction_method": artifact["point_prediction_method"],
        "decision_threshold": float(artifact.get("threshold", 0.5)),
        "training_cohort_fingerprint": artifact["training_cohort_fingerprint"],
        "training_feature_cohort_fingerprint": artifact[
            "training_feature_cohort_fingerprint"
        ],
        "inputs": {
            "model": {"file": model_path.name, "sha256": sha256_file(model_path)},
            "features": {"file": features_path.name, "sha256": sha256_file(features_path)},
        },
        "output": {
            "predictions": {
                "file": predictions_path.name,
                "sha256": sha256_file(predictions_path),
            }
        },
        "limitation": (
            "Integrity chain only; does not prove cohort independence, prospective "
            "registration, or correct outcome ascertainment."
        ),
    }
    audit_path = Path(f"{predictions_path}.audit.json")
    _write_json(audit_path, payload)
    return audit_path


def verify_prediction_audit(
    predictions_path: str | Path,
    audit_path: str | Path | None = None,
) -> dict:
    """Fail closed if a prediction file no longer matches its locked audit."""
    predictions_path = Path(predictions_path)
    audit_path = Path(audit_path) if audit_path else Path(f"{predictions_path}.audit.json")
    payload = _load_json(audit_path)
    if payload.get("schema_version") != 2 or payload.get("stage") != "locked_prediction":
        raise ValueError("unsupported locked-prediction audit schema")
    if payload.get("model_refit") is not False:
        raise ValueError("prediction audit does not attest locked inference")
    if payload.get("point_prediction_method") != "full_cohort_locked_model":
        raise ValueError("prediction audit lacks the locked point-prediction contract")
    recorded = payload.get("output", {}).get("predictions", {})
    if recorded.get("sha256") != sha256_file(predictions_path):
        raise ValueError("prediction file SHA-256 does not match locked-prediction audit")
    if recorded.get("file") != predictions_path.name:
        raise ValueError("prediction filename does not match locked-prediction audit")
    return {"payload": payload, "audit_sha256": sha256_file(audit_path),
            "audit_file": audit_path.name}


def write_evaluation_audit(
    predictions_path: str | Path,
    prediction_audit: dict,
    labels_path: str | Path,
    metrics_path: str | Path,
    bootstrap: int,
) -> Path:
    """Bind verified predictions and outcomes to the emitted metric report."""
    predictions_path, labels_path = Path(predictions_path), Path(labels_path)
    metrics_path = Path(metrics_path)
    source = prediction_audit["payload"]
    payload = {
        "schema_version": 1,
        "stage": "locked_external_evaluation",
        "research_only": True,
        "settings": {"bootstrap_draws": int(bootstrap), "primary_denominator": "all_patients"},
        "inputs": {
            "model_sha256": source["inputs"]["model"]["sha256"],
            "features_sha256": source["inputs"]["features"]["sha256"],
            "predictions": {"file": predictions_path.name,
                            "sha256": sha256_file(predictions_path)},
            "prediction_audit": {"file": prediction_audit["audit_file"],
                                 "sha256": prediction_audit["audit_sha256"]},
            "labels": {"file": labels_path.name, "sha256": sha256_file(labels_path)},
        },
        "output": {"metrics": {"file": metrics_path.name,
                                "sha256": sha256_file(metrics_path)}},
        "limitation": (
            "Integrity chain only; preregistration and independent-cohort status require "
            "external governance evidence."
        ),
    }
    audit_path = Path(f"{metrics_path}.audit.json")
    _write_json(audit_path, payload)
    return audit_path
