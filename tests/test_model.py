from pathlib import Path
import joblib
import pandas as pd
import pytest
from tdln_pd1.hybrid import augment_with_twin
from tdln_pd1.model import (evaluate_external, evaluate_probabilities,
                            leave_one_site_out, predict_locked, train_evaluate,
                            build_model, _safe_group_splits)
from tdln_pd1.spatial_dynamics import SPECIES, ZONE_NAMES
from tdln_pd1.simulator import TwinParameters


def _zero_member_uncertainty(predictions: pd.DataFrame) -> pd.DataFrame:
    """Attach a valid explicit v10 no-ensemble contract to hand-built predictions."""
    return predictions.assign(
        uncertainty_method="patient_bootstrap_refit",
        uncertainty_bootstrap_unit="patient",
        uncertainty_bootstrap_groups=1,
        uncertainty_bootstrap_members=0,
        uncertainty_available=False,
        bootstrap_probability_median=float("nan"),
        probability_lower_95=float("nan"),
        probability_upper_95=float("nan"),
        probability_interval_width=float("nan"),
    )


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
    assert artifact["artifact_version"] == 10
    assert artifact["selected_C"] in artifact["C_grid"]
    assert len(artifact["training_patient_hashes"]) == 8
    assert len(artifact["training_cohort_fingerprint"]) == 64
    assert len(artifact["training_feature_row_hashes"]) == 8
    assert len(artifact["training_feature_cohort_fingerprint"]) == 64
    assert len(artifact["bootstrap_models"]) == 8
    assert metrics["prediction_uncertainty"]["members"] == 8
    assert metrics["prediction_uncertainty"]["unit"] == "patient"
    assert artifact["uncertainty_method"] == "patient_bootstrap_refit"
    assert artifact["uncertainty_bootstrap_design"] == {
        "unit": "patient", "groups": 8, "members": 8, "seed": 42,
    }
    assert metrics["selection_metric"] == "negative_log_loss"
    assert metrics["probability_estimand"].startswith("unweighted_logistic_likelihood")
    assert "n_cells" not in artifact["features"]
    assert artifact["excluded_features"]["n_cells"] == "technical cell-count/QC feature"
    assert metrics["feature_schema"]["included"] == ["sig_exhaustion", "frac_cd8_t"]
    assert set(artifact["training_domain"]["features"]) == {"sig_exhaustion", "frac_cd8_t"}
    assert metrics["mechanism_data_fusion_ablation"]["status"] == "not_evaluable"


def test_probability_model_preserves_observed_prevalence_without_class_reweighting():
    features = pd.DataFrame({"constant": [1.0] * 20})
    outcomes = pd.Series([0] * 16 + [1] * 4)
    model = build_model(["constant"], c=100.0).fit(features, outcomes)
    classifier = model.named_steps["classifier"]
    assert classifier.class_weight is None
    probabilities = model.predict_proba(features)[:, 1]
    assert probabilities.mean() == pytest.approx(outcomes.mean(), abs=0.02)


def test_training_reports_paired_mechanism_data_fusion_ablation(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(12)],
        "sig_exhaustion": [0.10, 0.21, 0.31, 0.39, 0.48, 0.59,
                           0.62, 0.73, 0.79, 0.91, 1.02, 1.08],
        "twin_effector_gain": [0.3, 0.8, 0.4, 0.9, 0.5, 1.0,
                               1.4, 1.1, 1.7, 1.3, 1.9, 1.5],
        "liver_cf_tumour_relative_reduction": [
            0.01, 0.04, 0.02, 0.05, 0.03, 0.06,
            0.10, 0.08, 0.13, 0.09, 0.15, 0.12,
        ],
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id,
        "response": [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    })
    metrics = train_evaluate(features, labels, tmp_path, folds=3, ensemble_members=4)
    ablation = metrics["mechanism_data_fusion_ablation"]
    assert ablation["status"] == "evaluated"
    assert ablation["design"] == "paired_nested_patient_level_cross_validation_same_outer_folds"
    assert ablation["data_only"]["features"] == ["sig_exhaustion"]
    assert ablation["mechanism_only"]["features"] == [
        "twin_effector_gain", "liver_cf_tumour_relative_reduction",
    ]
    assert set(ablation["hybrid_minus_data_only"]) == {
        "roc_auc", "average_precision", "brier_score_improvement", "log_loss_improvement",
    }
    paired = ablation["paired_bootstrap"]
    assert paired["draws"] == 500
    assert paired["seed"] == 42
    assert paired["unit"] == "patient"
    assert paired["groups"] == 12
    assert paired["estimand"] == "paired_difference_of_fixed_out_of_fold_predictions"
    assert set(paired["intervals"]) == set(ablation["hybrid_minus_data_only"])
    assert all(value["valid_draws"] > 0 for value in paired["intervals"].values())
    predictions = pd.read_csv(tmp_path / "cv_predictions.csv")
    assert predictions.data_only_probability.notna().all()
    assert predictions.mechanism_only_probability.notna().all()


def test_training_and_fusion_ablation_bootstrap_whole_centres(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(12)],
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                           0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
        "twin_effector_gain": [0.2, 0.4, 0.3, 0.5, 0.4, 0.6,
                               1.0, 0.8, 1.2, 0.9, 1.4, 1.1],
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id,
        "response": [0, 1, 0, 1] * 3,
        "center": ["a"] * 2 + ["b"] * 4 + ["c"] * 6,
    })
    metrics = train_evaluate(features, labels, tmp_path, folds=3, ensemble_members=2)
    assert metrics["bootstrap_design"]["unit"] == "centre"
    assert metrics["bootstrap_design"]["groups"] == 3
    paired = metrics["mechanism_data_fusion_ablation"]["paired_bootstrap"]
    assert paired["unit"] == "centre"
    assert paired["groups"] == 3
    assert metrics["prediction_uncertainty"]["method"] == "centre_cluster_bootstrap_refit"
    assert metrics["prediction_uncertainty"]["unit"] == "centre"
    artifact = joblib.load(tmp_path / "model.joblib")
    assert artifact["uncertainty_method"] == "centre_cluster_bootstrap_refit"
    assert artifact["uncertainty_bootstrap_design"] == {
        "unit": "centre", "groups": 3, "members": 2, "seed": 42,
    }
    second_scaler = artifact["bootstrap_models"][1].named_steps["prep"].named_transformers_[
        "numeric"
    ].named_steps["scale"]
    assert second_scaler.n_samples_seen_ == 14


def test_training_rejects_empty_center_identifiers(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1.0],
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id, "response": [0] * 4 + [1] * 4,
        "center": ["a", "a", "b", "b", "c", "c", "d", ""],
    })
    with pytest.raises(ValueError, match="center must be non-empty"):
        train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)


def test_locked_prediction_enforces_embedded_preprocessing_contract(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    vocabulary_hash = "a" * 64
    features = pd.DataFrame({
        "patient_id": patients,
        "frac_cd8_t": [0.2, 0.3, 0.25, 0.35, 0.6, 0.7, 0.65, 0.75],
        "qc_cell_type_vocabulary_sha256": [vocabulary_hash] * 8,
    })
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    artifact = joblib.load(tmp_path / "model.joblib")
    assert artifact["preprocessing_contracts"] == {
        "qc_cell_type_vocabulary_sha256": vocabulary_hash,
    }
    assert artifact["preprocessing_contract_sources"] == {
        "qc_cell_type_vocabulary_sha256": "embedded_feature_qc",
    }
    external = features.iloc[[0]].copy()
    external["patient_id"] = "external"
    assert predict_locked(tmp_path / "model.joblib", external).probability.notna().all()
    external["qc_cell_type_vocabulary_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="preprocessing contract mismatch"):
        predict_locked(tmp_path / "model.joblib", external)
    with pytest.raises(ValueError, match="missing preprocessing contract"):
        predict_locked(
            tmp_path / "model.joblib",
            external.drop(columns="qc_cell_type_vocabulary_sha256"),
        )


def test_locked_prediction_rejects_cell_quality_policy_drift(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({
        "patient_id": patients,
        "sig_cytotoxic": [0.2, 0.3, 0.25, 0.35, 0.6, 0.7, 0.65, 0.75],
        "qc_cell_quality_policy_sha256": ["a" * 64] * 8,
    })
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=0)
    artifact = joblib.load(tmp_path / "model.joblib")
    assert artifact["preprocessing_contracts"][
        "qc_cell_quality_policy_sha256"] == "a" * 64
    external = features.iloc[[0]].copy()
    external["patient_id"] = "external"
    external["qc_cell_quality_policy_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="preprocessing contract mismatch"):
        predict_locked(tmp_path / "model.joblib", external)


def test_locked_prediction_rejects_spatial_fusion_semantic_drift(tmp_path: Path):
    base = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "sig_tcf7_memory": [0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1.0],
    })
    training = augment_with_twin(base, days=1, dose_days=(0,))
    labels = pd.DataFrame({
        "patient_id": base.patient_id, "response": [0] * 4 + [1] * 4,
    })
    train_evaluate(training, labels, tmp_path, folds=2, ensemble_members=0)
    artifact = joblib.load(tmp_path / "model.joblib")
    contract = training.qc_spatial_fusion_contract_sha256.unique().tolist()
    assert artifact["preprocessing_contracts"][
        "qc_spatial_fusion_contract_sha256"
    ] == contract[0]

    external_base = pd.DataFrame({
        "patient_id": ["external"], "sig_tcf7_memory": [0.5],
    })
    measured_states = pd.DataFrame([
        {"patient_id": "external", "species": species, "zone": zone,
         "amount": 0.0, "zone_observed": zone == "paracortex"}
        for species in SPECIES for zone in ZONE_NAMES
    ])
    external = augment_with_twin(
        external_base, days=1, dose_days=(0,), spatial_states=measured_states,
    )
    with pytest.raises(ValueError, match="qc_spatial_fusion_contract_sha256"):
        predict_locked(tmp_path / "model.joblib", external)
    different_horizon = augment_with_twin(
        external_base, days=2, dose_days=(0,),
    )
    with pytest.raises(ValueError, match="qc_mechanistic_simulation_contract_sha256"):
        predict_locked(tmp_path / "model.joblib", different_horizon)
    different_parameters = augment_with_twin(
        external_base, days=1, dose_days=(0,),
        systemic_parameters=TwinParameters(tumour_growth=0.2),
    )
    with pytest.raises(ValueError, match="qc_mechanistic_parameter_contract_sha256"):
        predict_locked(tmp_path / "model.joblib", different_parameters)


def test_spatial_evidence_support_marks_lower_coverage_as_domain_shift(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "x": [0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1.0],
        "qc_spatial_zone_coverage": [0.8, 0.9, 1.0, 0.8, 0.9, 1.0, 0.9, 0.8],
        "qc_spatial_feature_prior_imputed_fraction": [
            0.2, 0.1, 0.0, 0.2, 0.1, 0.0, 0.1, 0.2,
        ],
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id, "response": [0] * 4 + [1] * 4,
    })
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=0)
    artifact = joblib.load(tmp_path / "model.joblib")
    assert artifact["spatial_evidence_support"]["features"][
        "qc_spatial_zone_coverage"
    ] == {"minimum": 0.8, "maximum": 1.0}

    external = pd.DataFrame({
        "patient_id": ["external"], "x": [0.5],
        "qc_spatial_zone_coverage": [1 / 6],
        "qc_spatial_feature_prior_imputed_fraction": [5 / 6],
    })
    prediction = predict_locked(tmp_path / "model.joblib", external)
    assert prediction.spatial_evidence_check_available.item()
    assert prediction.spatial_evidence_outside_training_range.item()
    assert prediction.domain_shift_flag.item()


def test_spatial_evidence_support_requires_external_qc(tmp_path: Path):
    features = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(8)],
        "x": [0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1.0],
        "qc_spatial_zone_coverage": [1.0] * 8,
    })
    labels = pd.DataFrame({
        "patient_id": features.patient_id, "response": [0] * 4 + [1] * 4,
    })
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=0)
    with pytest.raises(ValueError, match="missing spatial evidence QC"):
        predict_locked(
            tmp_path / "model.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": [0.5]}),
        )


def test_locked_prediction_rejects_unreviewed_model_features_but_allows_audit_columns(
    tmp_path: Path,
):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({
        "patient_id": patients,
        "sig_exhaustion": [0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4],
        "n_cells": [100] * 8,
        "qc_feature_audit_sha256": ["a" * 64] * 8,
    })
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    external = pd.DataFrame({
        "patient_id": ["external"],
        "sig_exhaustion": [0.5],
        "n_cells": [90],
        "n_cells_tdln": [80],
        "qc_feature_audit_sha256": ["a" * 64],
        "qc_external_batch": ["recorded_not_modelled"],
    })
    assert predict_locked(tmp_path / "model.joblib", external).probability.notna().all()
    for column in ("sig_unreviewed", "twin_new_endpoint"):
        drifted = external.assign(**{column: 1.0})
        with pytest.raises(ValueError, match="unexpected unlocked model features"):
            predict_locked(tmp_path / "model.joblib", drifted)
    with pytest.raises(ValueError, match="preprocessing contract mismatch"):
        predict_locked(
            tmp_path / "model.joblib", external.assign(frac_new_annotation=1.0)
        )


def test_probability_evaluation_reports_calibration_and_intervals():
    metrics = evaluate_probabilities([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], bootstrap=20)
    assert metrics["roc_auc"] == 1.0
    assert metrics["brier_score"] < 0.1
    assert "calibration_slope" in metrics
    assert len(metrics["calibration_curve"]) == 2
    assert sum(row["n_patients"] for row in metrics["calibration_curve"]) == 4
    assert metrics["expected_calibration_error"] == pytest.approx(0.15)
    assert metrics["calibration_in_the_large"] == pytest.approx(0.0)
    assert metrics["bootstrap"]["expected_calibration_error"]["valid_draws"] > 0
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
    assert predictions.uncertainty_method.eq("patient_bootstrap_refit").all()
    assert predictions.uncertainty_bootstrap_unit.eq("patient").all()
    assert predictions.uncertainty_bootstrap_groups.eq(8).all()
    assert predictions.uncertainty_bootstrap_members.eq(8).all()
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
        "point_prediction_method": ["full_cohort_locked_model"] * 4,
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
        "point_prediction_method": ["full_cohort_locked_model"] * 2,
    })
    predictions = _zero_member_uncertainty(predictions)
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


def test_locked_prediction_rejects_legacy_or_ambiguous_point_estimate_artifacts(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({"patient_id": patients, "x": range(8)})
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    artifact = joblib.load(tmp_path / "model.joblib")
    artifact["artifact_version"] = 8
    joblib.dump(artifact, tmp_path / "legacy.joblib")
    with pytest.raises(ValueError, match="version 10"):
        predict_locked(
            tmp_path / "legacy.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": [0.5]}),
        )


def test_locked_prediction_rejects_tampered_uncertainty_bootstrap_contract(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({"patient_id": patients, "x": range(8)})
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    artifact = joblib.load(tmp_path / "model.joblib")
    artifact["uncertainty_bootstrap_design"]["unit"] = "centre"
    joblib.dump(artifact, tmp_path / "tampered.joblib")
    with pytest.raises(ValueError, match="invalid v10 uncertainty-bootstrap"):
        predict_locked(
            tmp_path / "tampered.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": [0.5]}),
        )
    artifact = joblib.load(tmp_path / "model.joblib")
    artifact["bootstrap_models"] = artifact["bootstrap_models"][:-1]
    joblib.dump(artifact, tmp_path / "missing_member.joblib")
    with pytest.raises(ValueError, match="invalid v10 uncertainty-bootstrap"):
        predict_locked(
            tmp_path / "missing_member.joblib",
            pd.DataFrame({"patient_id": ["external"], "x": [0.5]}),
        )


def test_external_evaluation_rejects_tampered_v10_uncertainty_rows(tmp_path: Path):
    patients = [f"p{i}" for i in range(8)]
    features = pd.DataFrame({"patient_id": patients, "x": range(8)})
    labels = pd.DataFrame({"patient_id": patients, "response": [0] * 4 + [1] * 4})
    train_evaluate(features, labels, tmp_path, folds=2, ensemble_members=2)
    external = pd.DataFrame({
        "patient_id": [f"q{i}" for i in range(8)],
        "x": [value + 0.1 for value in range(8)],
    })
    external_labels = labels.assign(patient_id=external.patient_id)
    predictions = predict_locked(tmp_path / "model.joblib", external)
    mixed = predictions.copy()
    mixed.loc[0, "uncertainty_method"] = "centre_cluster_bootstrap_refit"
    with pytest.raises(ValueError, match="one valid v10 uncertainty method"):
        evaluate_external(mixed, external_labels, bootstrap=0)
    unavailable = predictions.assign(uncertainty_available=False)
    with pytest.raises(ValueError, match="availability conflicts"):
        evaluate_external(unavailable, external_labels, bootstrap=0)
    inverted = predictions.copy()
    inverted.loc[0, "probability_lower_95"] = 1.0
    with pytest.raises(ValueError, match="invalid locked-model bootstrap"):
        evaluate_external(inverted, external_labels, bootstrap=0)


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
        "spatial_evidence_check_available": [True] * 6,
        "spatial_evidence_outside_training_range": [True, False, False, False, False, False],
        "decision_threshold": [0.5] * 6,
        "point_prediction_method": ["full_cohort_locked_model"] * 6,
    })
    predictions = _zero_member_uncertainty(predictions)
    labels = pd.DataFrame({"patient_id": predictions.patient_id,
                           "response": [0, 0, 0, 1, 1, 1]})
    metrics = evaluate_external(predictions, labels, bootstrap=0)
    assert metrics["n_patients"] == 6
    assert metrics["domain_shift_count"] == 1
    assert metrics["spatial_evidence_applicability"] == {
        "checked": True, "outside_count": 1, "outside_fraction": pytest.approx(1 / 6),
        "role": "evidence_support_not_fitted_biological_predictor",
    }
    sensitivity = metrics["in_domain_sensitivity_analysis"]
    assert sensitivity["role"] == "prespecified_sensitivity_analysis_not_primary"
    assert sensitivity["status"] == "evaluated"
    assert sensitivity["n_patients"] == 5
    assert sensitivity["coverage"] == pytest.approx(5 / 6)


def test_external_evaluation_rejects_spatial_domain_flag_contradiction():
    predictions = pd.DataFrame({
        "patient_id": ["a", "b"], "probability": [0.2, 0.8],
        "training_overlap": [False, False], "overlap_check_available": [True, True],
        "training_feature_reuse": [False, False],
        "feature_reuse_check_available": [True, True],
        "domain_check_available": [True, True], "domain_shift_flag": [False, False],
        "spatial_evidence_check_available": [True, True],
        "spatial_evidence_outside_training_range": [True, False],
        "decision_threshold": [0.5, 0.5],
        "point_prediction_method": ["full_cohort_locked_model"] * 2,
    })
    predictions = _zero_member_uncertainty(predictions)
    labels = pd.DataFrame({"patient_id": ["a", "b"], "response": [0, 1]})
    with pytest.raises(ValueError, match="conflict with domain-shift"):
        evaluate_external(predictions, labels, bootstrap=0)


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
        "point_prediction_method": ["full_cohort_locked_model"] * 6,
    })
    predictions = _zero_member_uncertainty(predictions)
    labels = pd.DataFrame({"patient_id": predictions.patient_id,
                           "response": [0, 0, 0, 1, 1, 1]})
    metrics = evaluate_external(predictions, labels, bootstrap=0)
    assert metrics["decision_threshold"] == 0.7
    assert metrics["sensitivity"] == pytest.approx(1 / 3)
    assert metrics["specificity"] == 1.0


def test_external_evaluation_requires_one_locked_point_prediction_method():
    predictions = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"],
        "probability": [0.1, 0.2, 0.8, 0.9],
        "training_overlap": [False] * 4,
        "overlap_check_available": [True] * 4,
        "training_feature_reuse": [False] * 4,
        "feature_reuse_check_available": [True] * 4,
        "domain_check_available": [True] * 4,
        "domain_shift_flag": [False] * 4,
        "decision_threshold": [0.5] * 4,
        "point_prediction_method": ["bootstrap_median"] * 4,
    })
    labels = pd.DataFrame({
        "patient_id": predictions.patient_id, "response": [0, 0, 1, 1],
    })
    with pytest.raises(ValueError, match="full_cohort_locked_model"):
        evaluate_external(predictions, labels, bootstrap=0)
    with pytest.raises(ValueError, match="point-prediction method"):
        evaluate_external(predictions.drop(columns="point_prediction_method"), labels, bootstrap=0)


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
        "point_prediction_method": ["full_cohort_locked_model"] * 10,
    })
    predictions = _zero_member_uncertainty(predictions)
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
    for invalid in (float("inf"), -0.01, 1.01):
        with pytest.raises(ValueError, match=r"within \[0, 1\]"):
            evaluate_probabilities(
                [0, 0, 1, 1], [0.1, invalid, 0.8, 0.9], bootstrap=0,
            )
    for invalid_outcomes in ([0, 0.5, 1, 1], [0, "no", 1, 1]):
        with pytest.raises(ValueError, match="only 0 and 1"):
            evaluate_probabilities(
                invalid_outcomes, [0.1, 0.2, 0.8, 0.9], bootstrap=0,
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
