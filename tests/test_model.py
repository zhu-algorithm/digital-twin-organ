from pathlib import Path
import joblib
import pandas as pd
import pytest
from tdln_pd1.model import (evaluate_external, evaluate_probabilities,
                            leave_one_site_out, predict_locked, train_evaluate,
                            _safe_group_splits)


def test_train_evaluate_writes_reproducible_artifacts(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4],
        "frac_cd8_t": [0.2, 0.3, 0.2, 0.3, 0.7, 0.8, 0.7, 0.8],
        "n_cells": [100, 200, 300, 400, 500, 600, 700, 800],
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
    assert artifact["artifact_version"] == 8
    assert artifact["selected_C"] in artifact["C_grid"]
    assert len(artifact["training_patient_hashes"]) == 8
    assert len(artifact["training_cohort_fingerprint"]) == 64
    assert len(artifact["training_feature_row_hashes"]) == 8
    assert len(artifact["training_feature_cohort_fingerprint"]) == 64
    assert len(artifact["bootstrap_models"]) == 8
    assert metrics["prediction_uncertainty"]["members"] == 8
    assert metrics["selection_metric"] == "negative_log_loss"
    assert "n_cells" not in artifact["features"]
    assert artifact["excluded_features"]["n_cells"] == "technical cell-count/QC feature"
    assert metrics["feature_schema"]["included"] == ["sig_exhaustion", "frac_cd8_t"]
    assert set(artifact["training_domain"]["features"]) == {"sig_exhaustion", "frac_cd8_t"}


def test_probability_evaluation_reports_calibration_and_intervals():
    metrics = evaluate_probabilities([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=20)
    assert metrics["roc_auc"] == 1.0
    assert metrics["brier_score"] < 0.1
    assert "calibration_slope" in metrics
    assert metrics["bootstrap"]["roc_auc"]["valid_draws"] > 0
    assert metrics["bootstrap"]["log_loss"]["valid_draws"] > 0
    assert metrics["bootstrap"]["sensitivity"]["valid_draws"] > 0
    assert metrics["bootstrap"]["specificity"]["valid_draws"] > 0
    assert metrics["bootstrap_design"]["unit"] == "patient"
    assert len(metrics["decision_curve"]) == 9
    assert metrics["decision_curve"][4]["threshold"] == 0.5


def test_probability_evaluation_uses_whole_centre_cluster_bootstrap():
    metrics = evaluate_probabilities(
        [0, 1, 0, 1, 0, 1],
        [0.1, 0.8, 0.2, 0.7, 0.3, 0.9],
        bootstrap=20,
        bootstrap_groups=["a", "a", "b", "b", "c", "c"],
    )
    assert metrics["bootstrap_design"]["unit"] == "centre"
    assert metrics["bootstrap_design"]["groups"] == 3
    assert metrics["bootstrap"]["roc_auc"]["valid_draws"] == 20
    assert metrics["bootstrap"]["balanced_accuracy"]["valid_draws"] == 20
    assert metrics["bootstrap"]["sensitivity"]["lower_95"] <= 1.0
    assert metrics["bootstrap"]["specificity"]["upper_95"] >= 1.0


def test_cluster_bootstrap_rejects_unusable_groups():
    with pytest.raises(ValueError, match="at least two"):
        evaluate_probabilities(
            [0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=10,
            bootstrap_groups=["only"] * 4,
        )


def test_locked_prediction_and_external_evaluation(tmp_path: Path):
    features = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)],
                             "x": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=8)
    training_predictions = predict_locked(tmp_path / "model.joblib", features)
    assert training_predictions.training_overlap.all()
    with pytest.raises(ValueError, match="overlaps"):
        evaluate_external(training_predictions, labels, bootstrap=10)
    external_features = pd.DataFrame({
        "patient_id": [f"q{i}" for i in range(8)],
        "x": [0.15, 0.25, 0.35, 0.39, 1.11, 1.19, 1.29, 1.39],
    })
    external_labels = labels.assign(patient_id=external_features.patient_id)
    predictions = predict_locked(tmp_path / "model.joblib", external_features)
    artifact = joblib.load(tmp_path / "model.joblib")
    expected_point = artifact["model"].predict_proba(external_features[["x"]])[:, 1]
    assert predictions.probability.tolist() == pytest.approx(expected_point)
    assert predictions.point_prediction_method.eq("full_cohort_locked_model").all()
    assert predictions.probability.between(0, 1).all()
    assert predictions.uncertainty_available.all()
    assert predictions.probability_lower_95.le(predictions.bootstrap_probability_median).all()
    assert predictions.probability_upper_95.ge(predictions.bootstrap_probability_median).all()
    assert predictions.probability_interval_width.ge(0).all()
    assert not predictions.training_overlap.any()
    assert not predictions.training_feature_reuse.any()
    assert predictions.domain_check_available.all()
    assert not predictions.domain_shift_flag.any()
    metrics = evaluate_external(predictions, external_labels, bootstrap=10)
    assert metrics["evaluation_type"] == "external_locked_model"
    assert metrics["training_overlap_checked"] is True
    assert metrics["training_overlap_count"] == 0
    assert metrics["training_feature_reuse_count"] == 0
    assert metrics["applicability_checked"] is True
    assert metrics["domain_shift_count"] == 0
    assert metrics["in_domain_sensitivity_analysis"]["status"] == "evaluated"
    assert metrics["in_domain_sensitivity_analysis"]["n_patients"] == 8
    assert predictions.decision_threshold.eq(0.5).all()
    with pytest.raises(ValueError):
        predict_locked(tmp_path / "model.joblib", pd.DataFrame({"patient_id": ["z"]}))


def test_patient_alignment_fails_closed_instead_of_silent_inner_join(tmp_path: Path):
    features = pd.DataFrame({"patient_id": ["p0", "p1", "p2", "p3"], "x": [0, 1, 2, 3]})
    labels = pd.DataFrame({"patient_id": ["p0", "p1", "p2", "other"], "response": [0, 0, 1, 1]})
    with pytest.raises(ValueError, match="patient sets differ"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=4)


def test_all_model_entry_points_reject_missing_response_labels(tmp_path: Path):
    features = pd.DataFrame({"patient_id": ["p0", "p1", "p2", "p3"],
                             "x": [0.1, 0.2, 0.8, 0.9]})
    labels = pd.DataFrame({"patient_id": features.patient_id,
                           "response": [0, 0, 1, None]})
    with pytest.raises(ValueError, match="complete"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)

    predictions = pd.DataFrame({
        "patient_id": features.patient_id, "probability": [0.1, 0.2, 0.8, 0.9],
        "training_overlap": [False] * 4, "overlap_check_available": [True] * 4,
        "training_feature_reuse": [False] * 4,
        "feature_reuse_check_available": [True] * 4,
        "domain_check_available": [True] * 4,
        "domain_shift_flag": [False] * 4, "decision_threshold": [0.5] * 4,
    })
    with pytest.raises(ValueError, match="complete"):
        evaluate_external(predictions, labels, bootstrap=0)


def test_negative_ensemble_size_is_rejected(tmp_path: Path):
    features = pd.DataFrame({"patient_id": ["a", "b", "c", "d"], "x": [0, 1, 2, 3]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0, 0, 1, 1]})
    with pytest.raises(ValueError, match="ensemble_members"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=-1)


def test_outcome_named_and_exact_label_proxy_features_are_rejected(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    named = pd.DataFrame({"patient_id": patients, "response_score": range(8)})
    with pytest.raises(ValueError, match="outcome-like"):
        train_evaluate(named, labels, tmp_path, folds=2, ensemble_members=2)
    copied = pd.DataFrame({"patient_id": patients, "perfect_proxy": labels.response,
                           "safe_feature": range(8)})
    with pytest.raises(ValueError, match="exactly copy"):
        train_evaluate(copied, labels, tmp_path, folds=2, ensemble_members=2)


def test_non_numeric_metadata_must_not_enter_model_matrix(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({"patient_id": patients, "x": range(8), "batch": ["a", "b"] * 4})
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    with pytest.raises(ValueError, match="must be numeric"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)


def test_model_matrix_rejects_infinite_measurements_but_allows_missing_values(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    features = pd.DataFrame({
        "patient_id": patients,
        "x": [0.1, 0.2, float("nan"), 0.4, 1.1, 1.2, 1.3, 1.4],
    })
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    missing = pd.DataFrame({"patient_id": ["external"], "x": [float("nan")]})
    assert predict_locked(tmp_path / "model.joblib", missing).probability.notna().all()

    infinite = features.copy()
    infinite.loc[0, "x"] = float("inf")
    with pytest.raises(ValueError, match="not infinite"):
        train_evaluate(infinite, labels, tmp_path / "inf", folds=2, ensemble_members=2)
    with pytest.raises(ValueError, match="not infinite"):
        predict_locked(
            tmp_path / "model.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": [float("-inf")]}),
        )
    with pytest.raises(ValueError, match="must be numeric"):
        predict_locked(
            tmp_path / "model.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": ["unknown"]}),
        )


def test_external_audit_parses_csv_boole_but_rejects_ambiguous_values():
    predictions = pd.DataFrame({
        "patient_id": ["a", "b"], "probability": [0.2, 0.8],
        "training_overlap": ["false", "False"],
        "overlap_check_available": ["true", "True"],
        "domain_check_available": ["true", "True"],
        "domain_shift_flag": ["false", "False"],
        "feature_reuse_check_available": ["true", "True"],
        "training_feature_reuse": ["false", "False"],
        "decision_threshold": [0.5, 0.5],
    })
    labels = pd.DataFrame({"patient_id": ["a", "b"], "response": [0, 1]})
    assert evaluate_external(predictions, labels, bootstrap=0)["training_overlap_count"] == 0
    predictions.loc[0, "training_overlap"] = "unknown"
    with pytest.raises(ValueError, match="true/false"):
        evaluate_external(predictions, labels, bootstrap=0)


def test_locked_prediction_flags_out_of_domain_patients(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({"patient_id": patients,
                             "x": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4]})
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    external = pd.DataFrame({"patient_id": ["near", "far"], "x": [0.5, 100.0]})
    predictions = predict_locked(tmp_path / "model.joblib", external)
    assert predictions.domain_check_available.all()
    assert not predictions.loc[predictions.patient_id.eq("near"), "domain_shift_flag"].item()
    assert predictions.loc[predictions.patient_id.eq("far"), "domain_shift_flag"].item()
    assert predictions.loc[predictions.patient_id.eq("far"), "domain_max_abs_robust_z"].item() > 4


def test_external_evaluation_rejects_missing_domain_audit():
    predictions = pd.DataFrame({
        "patient_id": ["a", "b"], "probability": [0.2, 0.8],
        "training_overlap": [False, False], "overlap_check_available": [True, True],
        "feature_reuse_check_available": [True, True],
        "training_feature_reuse": [False, False],
    })
    labels = pd.DataFrame({"patient_id": ["a", "b"], "response": [0, 1]})
    with pytest.raises(ValueError, match="applicability"):
        evaluate_external(predictions, labels, bootstrap=0)


def test_renaming_training_patients_cannot_fake_external_validation(tmp_path: Path):
    features = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)],
                             "x": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4]})
    labels = pd.DataFrame({"patient_id": features.patient_id, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    renamed = features.assign(patient_id=[f"renamed_{i}" for i in range(8)])
    renamed_labels = labels.assign(patient_id=renamed.patient_id)
    predictions = predict_locked(tmp_path / "model.joblib", renamed)
    assert not predictions.training_overlap.any()
    assert predictions.training_feature_reuse.all()
    with pytest.raises(ValueError, match="reuses exact training feature rows"):
        evaluate_external(predictions, renamed_labels, bootstrap=0)


def test_external_primary_analysis_keeps_domain_shifted_patients():
    predictions = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(6)],
        "probability": [0.1, 0.2, 0.3, 0.7, 0.8, 0.9],
        "training_overlap": [False] * 6,
        "overlap_check_available": [True] * 6,
        "training_feature_reuse": [False] * 6,
        "feature_reuse_check_available": [True] * 6,
        "domain_check_available": [True] * 6,
        "domain_shift_flag": [True, False, False, False, False, False],
        "decision_threshold": [0.5] * 6,
    })
    labels = pd.DataFrame({"patient_id": predictions.patient_id,
                           "response": [0, 0, 0, 1, 1, 1]})
    metrics = evaluate_external(predictions, labels, bootstrap=0)
    assert metrics["n_patients"] == 6
    assert metrics["domain_shift_count"] == 1
    sensitivity = metrics["in_domain_sensitivity_analysis"]
    assert sensitivity["role"] == "prespecified_sensitivity_analysis_not_primary"
    assert sensitivity["status"] == "evaluated"
    assert sensitivity["n_patients"] == 5
    assert sensitivity["coverage"] == pytest.approx(5 / 6)


def test_external_evaluation_uses_the_locked_decision_threshold():
    predictions = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(6)],
        "probability": [0.1, 0.2, 0.3, 0.55, 0.65, 0.9],
        "training_overlap": [False] * 6,
        "overlap_check_available": [True] * 6,
        "training_feature_reuse": [False] * 6,
        "feature_reuse_check_available": [True] * 6,
        "domain_check_available": [True] * 6,
        "domain_shift_flag": [False] * 6,
        "decision_threshold": [0.7] * 6,
    })
    labels = pd.DataFrame({"patient_id": predictions.patient_id,
                           "response": [0, 0, 0, 1, 1, 1]})
    metrics = evaluate_external(predictions, labels, bootstrap=0)
    assert metrics["decision_threshold"] == 0.7
    assert metrics["sensitivity"] == pytest.approx(1 / 3)
    assert metrics["specificity"] == 1.0


def test_external_evaluation_exposes_center_heterogeneity_without_dropping_sites():
    predictions = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(10)],
        "probability": [0.1, 0.8, 0.2, 0.9, 0.15, 0.75, 0.25, 0.85, 0.2, 0.3],
        "training_overlap": [False] * 10,
        "overlap_check_available": [True] * 10,
        "training_feature_reuse": [False] * 10,
        "feature_reuse_check_available": [True] * 10,
        "domain_check_available": [True] * 10,
        "domain_shift_flag": [False] * 10,
        "decision_threshold": [0.5] * 10,
    })
    labels = pd.DataFrame({
        "patient_id": predictions.patient_id,
        "response": [0, 1, 0, 1, 0, 1, 0, 1, 0, 0],
        "center": ["a"] * 4 + ["b"] * 4 + ["single_class"] * 2,
    })
    metrics = evaluate_external(predictions, labels, bootstrap=20)
    stratified = metrics["center_stratified_evaluation"]
    assert metrics["n_patients"] == 10
    assert stratified["role"] == "descriptive_center_heterogeneity_not_primary"
    assert stratified["n_centers"] == 3
    assert stratified["n_centers_evaluable"] == 2
    by_center = {row["center"]: row for row in stratified["centers"]}
    assert by_center["a"]["status"] == "evaluated"
    assert by_center["single_class"]["status"] == "not_evaluable"
    assert metrics["bootstrap_design"]["unit"] == "centre"


def test_evaluation_rejects_invalid_probabilities_bootstrap_and_threshold():
    with pytest.raises(ValueError, match="finite"):
        evaluate_probabilities([0, 0, 1, 1], [0.1, float("nan"), 0.8, 0.9])
    with pytest.raises(ValueError, match="non-negative integer"):
        evaluate_probabilities([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=-1)
    with pytest.raises(ValueError, match="strictly"):
        evaluate_probabilities(
            [0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], decision_threshold=1.0,
        )


def test_leave_one_site_out_predicts_each_evaluable_patient_once():
    patient = [f"p{i}" for i in range(12)]
    features = pd.DataFrame({"patient_id": patient, "x": [0.1, 1.1, 0.2, 1.2] * 3})
    labels = pd.DataFrame({"patient_id": patient, "response": [0, 1, 0, 1] * 3,
                           "center": ["a"] * 4 + ["b"] * 4 + ["c"] * 4})
    metrics, predictions = leave_one_site_out(features, labels, bootstrap=10)
    assert metrics["n_sites_evaluated"] == 3
    assert metrics["evaluation_type"] == "leave_one_site_out"
    assert metrics["bootstrap_design"]["unit"] == "centre"
    assert metrics["bootstrap_design"]["groups"] == 3
    assert len(metrics["site_hyperparameter_selection"]) == 3
    assert all(row["mode"] == "nested_grid_search"
               for row in metrics["site_hyperparameter_selection"])
    assert predictions.selected_C.isin([0.01, 0.1, 1.0, 10.0, 100.0]).all()
    assert predictions.patient_id.nunique() == 12
    assert not predictions.patient_id.duplicated().any()


def test_site_cv_records_fixed_c_when_inner_tuning_is_not_evaluable():
    features = pd.DataFrame({"patient_id": ["a0", "a1", "b0", "b1"],
                             "x": [0.1, 0.9, 0.2, 0.8]})
    labels = pd.DataFrame({"patient_id": features.patient_id,
                           "response": [0, 1, 0, 1],
                           "center": ["a", "a", "b", "b"]})
    metrics, predictions = leave_one_site_out(features, labels, bootstrap=0)
    assert len(predictions) == 4
    assert predictions.selected_C.eq(1.0).all()
    assert all(row["mode"] == "fixed_C_inner_cv_not_evaluable"
               for row in metrics["site_hyperparameter_selection"])


def test_site_cv_keeps_single_class_held_out_centers_in_primary_analysis():
    features = pd.DataFrame({
        "patient_id": ["a0", "a1", "b0", "b1", "c0", "c1"],
        "x": [0.1, 0.9, 0.2, 0.8, 0.15, 0.25],
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id,
        "response": [0, 1, 0, 1, 0, 0],
        "center": ["a", "a", "b", "b", "single_class", "single_class"],
    })
    metrics, predictions = leave_one_site_out(features, labels, bootstrap=10)
    assert predictions.patient_id.nunique() == 6
    assert metrics["n_patients"] == 6
    assert metrics["n_sites_evaluated"] == 3
    assert metrics["skipped_sites"] == []
    by_center = {
        row["center"]: row
        for row in metrics["center_stratified_evaluation"]["centers"]
    }
    assert by_center["single_class"]["status"] == "not_evaluable"
    assert by_center["a"]["status"] == "evaluated"


def test_safe_group_splits_never_emit_single_class_partitions():
    y = pd.Series([0, 1, 0, 1, 0, 1, 0, 1])
    groups = pd.Series([f"p{i}" for i in range(8)])
    splits, n_splits = _safe_group_splits(y, groups, maximum_splits=3, seed=101)
    assert splits is not None and n_splits >= 2
    for train, test in splits:
        assert set(y.iloc[train]) == {0, 1}
        assert set(y.iloc[test]) == {0, 1}
