from pathlib import Path
import joblib
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
    metrics = train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=8)
    assert metrics["n_patients"] == 8
    assert (tmp_path / "model.joblib").exists()
    assert (tmp_path / "cv_predictions.csv").exists()
    predictions = pd.read_csv(tmp_path / "cv_predictions.csv")
    assert predictions.patient_id.nunique() == len(predictions)
    assert predictions.outer_fold.ge(1).all()
    assert predictions.fold_selected_C.isin([0.01, 0.1, 1.0, 10.0, 100.0]).all()
    artifact = joblib.load(tmp_path / "model.joblib")
    assert artifact["artifact_version"] == 4
    assert artifact["selected_C"] in artifact["C_grid"]
    assert len(artifact["training_patient_hashes"]) == 8
    assert len(artifact["training_cohort_fingerprint"]) == 64
    assert len(artifact["bootstrap_models"]) == 8
    assert metrics["prediction_uncertainty"]["members"] == 8
    assert metrics["selection_metric"] == "negative_log_loss"


def test_probability_evaluation_reports_calibration_and_intervals():
    metrics = evaluate_probabilities([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=20)
    assert metrics["roc_auc"] == 1.0
    assert metrics["brier_score"] < 0.1
    assert "calibration_slope" in metrics
    assert metrics["bootstrap"]["roc_auc"]["valid_draws"] > 0
    assert len(metrics["decision_curve"]) == 9
    assert metrics["decision_curve"][4]["threshold"] == 0.5


def test_locked_prediction_and_external_evaluation(tmp_path: Path):
    features = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)],
                             "x": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=8)
    training_predictions = predict_locked(tmp_path / "model.joblib", features)
    assert training_predictions.training_overlap.all()
    with pytest.raises(ValueError, match="overlaps"):
        evaluate_external(training_predictions, labels, bootstrap=10)
    external_features = features.assign(patient_id=[f"q{i}" for i in range(8)])
    external_labels = labels.assign(patient_id=external_features.patient_id)
    predictions = predict_locked(tmp_path / "model.joblib", external_features)
    assert predictions.probability.between(0, 1).all()
    assert predictions.uncertainty_available.all()
    assert predictions.probability_lower_95.le(predictions.probability).all()
    assert predictions.probability_upper_95.ge(predictions.probability).all()
    assert predictions.probability_interval_width.ge(0).all()
    assert not predictions.training_overlap.any()
    metrics = evaluate_external(predictions, external_labels, bootstrap=10)
    assert metrics["evaluation_type"] == "external_locked_model"
    assert metrics["training_overlap_checked"] is True
    assert metrics["training_overlap_count"] == 0
    assert predictions.decision_threshold.eq(0.5).all()
    with pytest.raises(ValueError):
        predict_locked(tmp_path / "model.joblib", pd.DataFrame({"patient_id": ["z"]}))


def test_patient_alignment_fails_closed_instead_of_silent_inner_join(tmp_path: Path):
    features = pd.DataFrame({"patient_id": ["p0", "p1", "p2", "p3"], "x": [0, 1, 2, 3]})
    labels = pd.DataFrame({"patient_id": ["p0", "p1", "p2", "other"], "response": [0, 0, 1, 1]})
    with pytest.raises(ValueError, match="patient sets differ"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=4)


def test_negative_ensemble_size_is_rejected(tmp_path: Path):
    features = pd.DataFrame({"patient_id": ["a", "b", "c", "d"], "x": [0, 1, 2, 3]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0, 0, 1, 1]})
    with pytest.raises(ValueError, match="ensemble_members"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=-1)


def test_external_audit_parses_csv_boole_but_rejects_ambiguous_values():
    predictions = pd.DataFrame({
        "patient_id": ["a", "b"], "probability": [0.2, 0.8],
        "training_overlap": ["false", "False"],
        "overlap_check_available": ["true", "True"],
    })
    labels = pd.DataFrame({"patient_id": ["a", "b"], "response": [0, 1]})
    assert evaluate_external(predictions, labels, bootstrap=0)["training_overlap_count"] == 0
    predictions.loc[0, "training_overlap"] = "unknown"
    with pytest.raises(ValueError, match="true/false"):
        evaluate_external(predictions, labels, bootstrap=0)


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
