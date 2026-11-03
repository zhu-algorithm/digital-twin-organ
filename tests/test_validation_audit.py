import json
import hashlib
from pathlib import Path

import joblib
import pandas as pd
import pytest

from tdln_pd1.validation_audit import (verify_evaluation_audit, verify_prediction_audit,
                                       write_evaluation_audit,
                                       write_prediction_audit)


def _files(tmp_path: Path):
    model = tmp_path / "model.joblib"
    joblib.dump({
        "artifact_version": 9,
        "training_cohort_fingerprint": "a" * 64,
        "training_feature_cohort_fingerprint": "b" * 64,
        "point_prediction_method": "full_cohort_locked_model",
    }, model)
    features = tmp_path / "features.csv"
    features.write_text("patient_id,x\np1,1\n")
    predictions = tmp_path / "predictions.csv"
    predictions.write_text("patient_id,probability\np1,0.7\n")
    return model, features, predictions


def test_prediction_and_evaluation_audits_bind_every_file(tmp_path: Path, monkeypatch):
    model, features, predictions = _files(tmp_path)
    monkeypatch.setattr(
        "tdln_pd1.validation_audit.predict_locked",
        lambda *_: pd.DataFrame({"patient_id": ["p1"], "probability": [0.7]}),
    )
    audit_path = write_prediction_audit(model, features, predictions, 1)
    verified = verify_prediction_audit(predictions, model_path=model, features_path=features)
    assert verified["payload"]["model_refit"] is False
    assert verified["payload"]["patient_count"] == 1
    assert verified["payload"]["decision_threshold"] == 0.5
    assert verified["payload"]["schema_version"] == 3
    assert verified["payload"]["point_prediction_method"] == "full_cohort_locked_model"
    assert verified["derivation_replayed"] is True

    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,response\np1,1\n")
    metrics = tmp_path / "metrics.json"
    metrics.write_text('{"roc_auc": 0.5}\n')
    evaluation_audit = write_evaluation_audit(
        predictions, verified, labels, metrics, bootstrap=500,
    )
    payload = json.loads(evaluation_audit.read_text())
    assert payload["stage"] == "locked_external_evaluation"
    assert payload["schema_version"] == 2
    assert payload["settings"]["bootstrap_draws"] == 500
    assert payload["inputs"]["prediction_audit"]["file"] == audit_path.name


def test_prediction_audit_rejects_posthoc_file_replacement(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    write_prediction_audit(model, features, predictions, 1)
    predictions.write_text("patient_id,probability\np1,0.1\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify_prediction_audit(predictions, model_path=model, features_path=features)


def test_prediction_audit_requires_versioned_training_fingerprints(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    joblib.dump({"artifact_version": 1}, model)
    with pytest.raises(ValueError, match="fingerprints"):
        write_prediction_audit(model, features, predictions, 1)


def test_prediction_audit_requires_replay_inputs(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    write_prediction_audit(model, features, predictions, 1)
    with pytest.raises(ValueError, match="requires the locked model"):
        verify_prediction_audit(predictions)


def test_evaluation_audit_replays_exact_metric_report(tmp_path: Path, monkeypatch):
    model, features, predictions = _files(tmp_path)
    replayed_predictions = pd.DataFrame({"patient_id": ["p1"], "probability": [0.7]})
    monkeypatch.setattr("tdln_pd1.validation_audit.predict_locked", lambda *_: replayed_predictions)
    prediction_audit = write_prediction_audit(model, features, predictions, 1)
    verified = verify_prediction_audit(predictions, model_path=model, features_path=features)
    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,response\np1,1\n")
    expected = {"roc_auc": 0.5, "validation_chain": {
        "prediction_audit_sha256": verified["audit_sha256"],
        "model_sha256": verified["payload"]["inputs"]["model"]["sha256"],
        "features_sha256": verified["payload"]["inputs"]["features"]["sha256"],
        "integrity_verified": True,
        "prediction_derivation_replayed": True,
    }}
    monkeypatch.setattr("tdln_pd1.validation_audit.evaluate_external", lambda *_: {"roc_auc": 0.5})
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(expected, indent=2) + "\n")
    audit = write_evaluation_audit(predictions, verified, labels, metrics, bootstrap=500)
    report = verify_evaluation_audit(
        audit, predictions, model, features, labels, metrics, prediction_audit,
    )
    assert report["valid"] is True
    assert report["evaluation_derivation_replayed"] is True

    metrics.write_text(json.dumps({"roc_auc": 0.9}, indent=2) + "\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify_evaluation_audit(
            audit, predictions, model, features, labels, metrics, prediction_audit,
        )

    forged = json.loads(audit.read_text())
    forged["output"]["metrics"]["sha256"] = hashlib.sha256(metrics.read_bytes()).hexdigest()
    audit.write_text(json.dumps(forged, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="derivation replay"):
        verify_evaluation_audit(
            audit, predictions, model, features, labels, metrics, prediction_audit,
        )
