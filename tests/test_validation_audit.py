import json
from pathlib import Path

import joblib
import pytest

from tdln_pd1.validation_audit import (verify_prediction_audit,
                                       write_evaluation_audit,
                                       write_prediction_audit)


def _files(tmp_path: Path):
    model = tmp_path / "model.joblib"
    joblib.dump({
        "artifact_version": 7,
        "training_cohort_fingerprint": "a" * 64,
        "training_feature_cohort_fingerprint": "b" * 64,
    }, model)
    features = tmp_path / "features.csv"
    features.write_text("patient_id,x\np1,1\n")
    predictions = tmp_path / "predictions.csv"
    predictions.write_text("patient_id,probability\np1,0.7\n")
    return model, features, predictions


def test_prediction_and_evaluation_audits_bind_every_file(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    audit_path = write_prediction_audit(model, features, predictions, 1)
    verified = verify_prediction_audit(predictions)
    assert verified["payload"]["model_refit"] is False
    assert verified["payload"]["patient_count"] == 1

    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,response\np1,1\n")
    metrics = tmp_path / "metrics.json"
    metrics.write_text('{"roc_auc": 0.5}\n')
    evaluation_audit = write_evaluation_audit(
        predictions, verified, labels, metrics, bootstrap=500,
    )
    payload = json.loads(evaluation_audit.read_text())
    assert payload["stage"] == "locked_external_evaluation"
    assert payload["settings"]["bootstrap_draws"] == 500
    assert payload["inputs"]["prediction_audit"]["file"] == audit_path.name


def test_prediction_audit_rejects_posthoc_file_replacement(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    write_prediction_audit(model, features, predictions, 1)
    predictions.write_text("patient_id,probability\np1,0.1\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify_prediction_audit(predictions)


def test_prediction_audit_requires_versioned_training_fingerprints(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    joblib.dump({"artifact_version": 1}, model)
    with pytest.raises(ValueError, match="fingerprints"):
        write_prediction_audit(model, features, predictions, 1)
