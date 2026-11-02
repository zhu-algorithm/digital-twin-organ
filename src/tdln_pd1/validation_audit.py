"""Cryptographic chain of custody for locked external validation files."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import joblib
import pandas as pd

from .data import load_labels
from .model import evaluate_external, predict_locked
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
        "schema_version": 3,
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
    model_path: str | Path | None = None,
    features_path: str | Path | None = None,
) -> dict:
    """Re-run locked inference and fail if any bound input or output differs."""
    predictions_path = Path(predictions_path)
    audit_path = Path(audit_path) if audit_path else Path(f"{predictions_path}.audit.json")
    payload = _load_json(audit_path)
    if payload.get("schema_version") != 3 or payload.get("stage") != "locked_prediction":
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
    if model_path is None or features_path is None:
        raise ValueError("prediction replay requires the locked model and feature table")
    model_path, features_path = Path(model_path), Path(features_path)
    for name, path in (("model", model_path), ("features", features_path)):
        record = payload.get("inputs", {}).get(name, {})
        if (record.get("file") != path.name or not path.is_file()
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"prediction replay {name} SHA-256 mismatch")
    replayed = predict_locked(model_path, pd.read_csv(features_path))
    if len(replayed) != payload.get("patient_count"):
        raise ValueError("prediction replay patient count mismatch")
    replayed_sha = hashlib.sha256(
        replayed.to_csv(index=False).encode("utf-8")
    ).hexdigest()
    if replayed_sha != recorded["sha256"]:
        raise ValueError("locked prediction derivation replay mismatch")
    return {"payload": payload, "audit_sha256": sha256_file(audit_path),
            "audit_file": audit_path.name, "derivation_replayed": True}


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
        "schema_version": 2,
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


def verify_evaluation_audit(
    audit_path: str | Path,
    predictions_path: str | Path,
    model_path: str | Path,
    features_path: str | Path,
    labels_path: str | Path,
    metrics_path: str | Path,
    prediction_audit_path: str | Path | None = None,
) -> dict:
    """Replay locked inference and evaluation, then verify the exact metric report."""
    audit_path = Path(audit_path)
    predictions_path, labels_path = Path(predictions_path), Path(labels_path)
    metrics_path = Path(metrics_path)
    payload = _load_json(audit_path)
    if payload.get("schema_version") != 2 or payload.get("stage") != "locked_external_evaluation":
        raise ValueError("unsupported locked-evaluation audit schema")
    settings = payload.get("settings", {})
    bootstrap = settings.get("bootstrap_draws")
    if (not isinstance(bootstrap, int) or isinstance(bootstrap, bool) or bootstrap < 0
            or settings.get("primary_denominator") != "all_patients"):
        raise ValueError("invalid locked-evaluation settings")
    verified_prediction = verify_prediction_audit(
        predictions_path, prediction_audit_path, model_path, features_path,
    )
    source_inputs = payload.get("inputs", {})
    expected_inputs = {
        "model_sha256": verified_prediction["payload"]["inputs"]["model"]["sha256"],
        "features_sha256": verified_prediction["payload"]["inputs"]["features"]["sha256"],
    }
    for name, expected in expected_inputs.items():
        if source_inputs.get(name) != expected:
            raise ValueError(f"evaluation audit {name} mismatch")
    bound_files = {
        "predictions": predictions_path,
        "prediction_audit": Path(prediction_audit_path) if prediction_audit_path else Path(
            f"{predictions_path}.audit.json"
        ),
        "labels": labels_path,
    }
    for name, path in bound_files.items():
        record = source_inputs.get(name, {})
        if (record.get("file") != path.name or not path.is_file()
                or record.get("sha256") != sha256_file(path)):
            raise ValueError(f"evaluation replay {name} SHA-256 mismatch")
    output = payload.get("output", {}).get("metrics", {})
    if (output.get("file") != metrics_path.name or not metrics_path.is_file()
            or output.get("sha256") != sha256_file(metrics_path)):
        raise ValueError("evaluation metric report SHA-256 mismatch")
    replayed = evaluate_external(
        pd.read_csv(predictions_path), load_labels(labels_path), bootstrap,
    )
    replayed["validation_chain"] = {
        "prediction_audit_sha256": verified_prediction["audit_sha256"],
        "model_sha256": verified_prediction["payload"]["inputs"]["model"]["sha256"],
        "features_sha256": verified_prediction["payload"]["inputs"]["features"]["sha256"],
        "integrity_verified": True,
        "prediction_derivation_replayed": True,
    }
    replayed_bytes = (json.dumps(replayed, indent=2) + "\n").encode("utf-8")
    if hashlib.sha256(replayed_bytes).hexdigest() != output["sha256"]:
        raise ValueError("locked evaluation derivation replay mismatch")
    return {
        "valid": True,
        "schema_version": 2,
        "prediction_derivation_replayed": True,
        "evaluation_derivation_replayed": True,
        "bootstrap_draws": bootstrap,
        "metrics_sha256": output["sha256"],
    }
