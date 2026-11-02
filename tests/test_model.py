from pathlib import Path
import pandas as pd
import pytest
from tdln_pd1.model import (evaluate_external, evaluate_probabilities,
                            leave_one_site_out, predict_locked, train_evaluate)


def test_train_evaluate_writes_reproducible_artifacts(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4],
        "frac_cd8_t": [0.2, 0.3, 0.2, 0.3, 0.7, 0.8, 0.7, 0.8],
    })
    labels = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)], "response": [0]*4 + [1]*4})
    metrics = train_evaluate(features, labels, tmp_path, folds=2)
    assert metrics["n_patients"] == 8
    assert (tmp_path / "model.joblib").exists()
    assert (tmp_path / "cv_predictions.csv").exists()


def test_probability_evaluation_reports_calibration_and_intervals():
    metrics = evaluate_probabilities([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=20)
    assert metrics["roc_auc"] == 1.0
    assert metrics["brier_score"] < 0.1
    assert "calibration_slope" in metrics
    assert metrics["bootstrap"]["roc_auc"]["valid_draws"] > 0


def test_locked_prediction_and_external_evaluation(tmp_path: Path):
    features = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)],
                             "x": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2)
    predictions = predict_locked(tmp_path / "model.joblib", features)
    assert predictions.probability.between(0, 1).all()
    metrics = evaluate_external(predictions, labels, bootstrap=10)
    assert metrics["evaluation_type"] == "external_locked_model"
    with pytest.raises(ValueError):
        predict_locked(tmp_path / "model.joblib", pd.DataFrame({"patient_id": ["z"]}))


def test_leave_one_site_out_predicts_each_evaluable_patient_once():
    patient = [f"p{i}" for i in range(12)]
    features = pd.DataFrame({"patient_id": patient, "x": [0.1, 1.1, 0.2, 1.2] * 3})
    labels = pd.DataFrame({"patient_id": patient, "response": [0, 1, 0, 1] * 3,
                           "center": ["a"] * 4 + ["b"] * 4 + ["c"] * 4})
    metrics, predictions = leave_one_site_out(features, labels, bootstrap=10)
    assert metrics["n_sites_evaluated"] == 3
    assert metrics["evaluation_type"] == "leave_one_site_out"
    assert predictions.patient_id.nunique() == 12
    assert not predictions.patient_id.duplicated().any()
