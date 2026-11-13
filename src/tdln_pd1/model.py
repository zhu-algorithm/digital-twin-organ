"""Training and evaluation without cell-level leakage."""

from __future__ import annotations

import json
import hashlib
import re
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, log_loss, roc_auc_score,
                             confusion_matrix)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


DEFAULT_C_GRID = (0.01, 0.1, 1.0, 10.0, 100.0)
PROBABILITY_ESTIMAND = "unweighted_logistic_likelihood_observed_development_prevalence"
FUSION_ABLATION_BOOTSTRAP_DRAWS = 500
FORBIDDEN_FEATURE_TOKENS = (
    "response", "responder", "outcome", "label", "recist",
    "progression_free", "overall_survival", "pfs", "os_months",
)
DOMAIN_MAX_ABS_Z = 4.0
DOMAIN_MAX_OUTSIDE_FRACTION = 0.2
DOMAIN_MAX_MISSINGNESS_INCREASE = 0.2
SPATIAL_EVIDENCE_COLUMNS = (
    "qc_spatial_zone_coverage",
    "qc_spatial_sampling_coverage_mean",
    "qc_spatial_sampling_coverage_min",
    "qc_spatial_feature_prior_imputed_fraction",
)
ARTIFACT_VERSION = 10


def _normalise_patient_ids(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    if "patient_id" not in frame:
        raise ValueError(f"{name} need patient_id")
    if frame.patient_id.isna().any():
        raise ValueError(f"{name} contain missing patient_id")
    out = frame.copy()
    out["patient_id"] = out.patient_id.astype(str).str.strip()
    if out.patient_id.eq("").any() or out.patient_id.duplicated().any():
        raise ValueError(f"{name} need unique, non-empty patient_id")
    return out


def _align_patients(features: pd.DataFrame, labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    features = _normalise_patient_ids(features, "features")
    labels = _normalise_patient_ids(labels, "labels")
    if "response" not in labels:
        raise ValueError("labels need response")
    response = pd.to_numeric(labels.response, errors="coerce")
    if response.isna().any() or not response.isin([0, 1]).all():
        raise ValueError("response must be complete and contain only 0 and 1")
    labels["response"] = response.astype(int)
    feature_ids, label_ids = set(features.patient_id), set(labels.patient_id)
    if feature_ids != label_ids:
        missing_labels = sorted(feature_ids - label_ids)
        missing_features = sorted(label_ids - feature_ids)
        raise ValueError(
            "patient sets differ: "
            f"missing_labels={missing_labels[:10]}, missing_features={missing_features[:10]}"
        )
    return features, labels


def _patient_hash(patient_id: str) -> str:
    return hashlib.sha256(patient_id.encode("utf-8")).hexdigest()


def _feature_row_hash(row: pd.Series, columns: list[str]) -> str:
    """Hash one locked-schema feature vector independently of patient identity."""
    encoded = []
    for column in columns:
        value = float(row[column])
        if np.isnan(value):
            canonical = "nan"
        elif np.isposinf(value):
            canonical = "+inf"
        elif np.isneginf(value):
            canonical = "-inf"
        elif value == 0.0:
            canonical = "0x0.0p+0"
        else:
            canonical = value.hex()
        encoded.append(f"{column}={canonical}")
    return hashlib.sha256("\x1f".join(encoded).encode("utf-8")).hexdigest()


def _boolean_audit_column(series: pd.Series, name: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    if not normalized.isin({"true", "false"}).all():
        raise ValueError(f"{name} must contain only true/false values")
    return normalized.eq("true")


def _validate_model_matrix(
    frame: pd.DataFrame,
    columns: list[str],
    context: str,
) -> None:
    """Require numeric, finite-or-missing measurements before sklearn sees them."""
    non_numeric = [
        column for column in columns
        if not pd.api.types.is_numeric_dtype(frame[column])
    ]
    if non_numeric:
        raise ValueError(f"{context} model features must be numeric: {sorted(non_numeric)}")
    non_finite = []
    for column in columns:
        values = frame[column].to_numpy(dtype=float)
        if np.isinf(values).any():
            non_finite.append(column)
    if non_finite:
        raise ValueError(
            f"{context} model features may be finite or missing, but not infinite: "
            f"{sorted(non_finite)}"
        )


def _technical_audit_feature(column: str) -> bool:
    """Return whether a column is explicitly excluded technical/QC metadata."""
    return column == "n_cells" or column.startswith("n_cells_") or column.startswith("qc_")


def _select_model_features(
    features: pd.DataFrame,
    labels: pd.DataFrame,
) -> tuple[list[str], dict[str, str]]:
    """Audit feature roles and reject common forms of target leakage."""
    candidates = [column for column in features.columns if column != "patient_id"]
    forbidden = [
        column for column in candidates
        if any(token in re.sub(r"[^a-z0-9]+", "_", column.lower())
               for token in FORBIDDEN_FEATURE_TOKENS)
    ]
    if forbidden:
        raise ValueError(f"forbidden outcome-like feature columns: {sorted(forbidden)}")
    excluded = {
        column: "technical cell-count/QC feature"
        for column in candidates if _technical_audit_feature(column)
    }
    included = [column for column in candidates if column not in excluded]
    non_numeric = [column for column in included
                   if not pd.api.types.is_numeric_dtype(features[column])]
    if non_numeric:
        raise ValueError(
            "model features must be numeric; keep metadata in the label/cohort table: "
            f"{sorted(non_numeric)}"
        )
    if not included:
        raise ValueError("no eligible biological or mechanistic model features remain")
    _validate_model_matrix(features, included, "training")
    response = labels.set_index("patient_id").loc[features.patient_id, "response"].to_numpy(dtype=float)
    copied = []
    for column in included:
        values = pd.to_numeric(features[column], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(values).all() and (
            np.allclose(values, response, rtol=0.0, atol=1e-12)
            or np.allclose(values, 1.0 - response, rtol=0.0, atol=1e-12)
        ):
            copied.append(column)
    if copied:
        raise ValueError(f"features exactly copy or invert the response labels: {sorted(copied)}")
    return included, excluded


def _preprocessing_contracts(features: pd.DataFrame) -> dict[str, str]:
    """Freeze constant cryptographic QC contracts without entering them into fitting."""
    contracts = {}
    for column in sorted(
        name for name in features if name.startswith("qc_") and name.endswith("_sha256")
    ):
        values = features[column].dropna().astype(str).str.strip().unique().tolist()
        if (len(values) != 1 or len(values[0]) != 64
                or any(character not in "0123456789abcdef" for character in values[0].lower())):
            raise ValueError(f"invalid or non-constant preprocessing contract: {column}")
        contracts[column] = values[0].lower()
    if "qc_cell_type_vocabulary_sha256" not in contracts:
        fraction_schema = tuple(sorted(
            column for column in features if column.startswith("frac_")
        ))
        contracts["qc_cell_type_vocabulary_sha256"] = hashlib.sha256(
            json.dumps(fraction_schema, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
    return contracts


def build_model(feature_columns: list[str], c: float = 1.0) -> Pipeline:
    prep = ColumnTransformer([
        ("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
        ]), feature_columns)
    ])
    return Pipeline([("prep", prep), ("classifier", LogisticRegression(
        C=c, max_iter=2000, class_weight=None, solver="liblinear", random_state=42,
    ))])


def validate_uncertainty_bootstrap_contract(artifact: dict) -> dict:
    """Require one internally consistent locked refit-uncertainty design."""
    method = artifact.get("uncertainty_method")
    design = artifact.get("uncertainty_bootstrap_design")
    expected_unit = {
        "patient_bootstrap_refit": "patient",
        "centre_cluster_bootstrap_refit": "centre",
    }.get(method)
    if expected_unit is None or not isinstance(design, dict):
        raise ValueError("locked model lacks the v10 uncertainty-bootstrap contract")
    groups, seed, members = design.get("groups"), design.get("seed"), design.get("members")
    models = artifact.get("bootstrap_models")
    training_patients = artifact.get("training_patients")
    if (design.get("unit") != expected_unit
            or not isinstance(groups, (int, np.integer)) or isinstance(groups, bool)
            or groups < (2 if expected_unit == "centre" else 1)
            or not isinstance(members, (int, np.integer)) or isinstance(members, bool)
            or members < 0 or not isinstance(models, list) or len(models) != members
            or seed != 42
            or not isinstance(training_patients, (int, np.integer))
            or isinstance(training_patients, bool) or training_patients < 1
            or (expected_unit == "patient" and groups != training_patients)
            or (expected_unit == "centre" and groups > training_patients)):
        raise ValueError("invalid v10 uncertainty-bootstrap contract")
    return {"method": method, "unit": expected_unit, "groups": int(groups),
            "members": int(members), "seed": 42}


def _nested_oof_for_columns(
    data: pd.DataFrame,
    y: pd.Series,
    columns: list[str],
    outer_cv: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, list[dict]]:
    """Generate group-safe nested-CV probabilities for one prespecified feature block."""
    probabilities = np.full(len(data), np.nan)
    tuning = []
    X = data[columns]
    for fold, (train_index, test_index) in enumerate(outer_cv, start=1):
        y_train = y.iloc[train_index]
        inner_cv, inner_splits = _safe_group_splits(
            y_train, data.patient_id.iloc[train_index], 3, fold,
        )
        if inner_cv is not None:
            search = GridSearchCV(
                build_model(columns), {"classifier__C": list(DEFAULT_C_GRID)},
                scoring="neg_log_loss", cv=inner_cv, refit=True, error_score="raise",
            )
            search.fit(X.iloc[train_index], y_train)
            fitted = search.best_estimator_
            best_c = float(search.best_params_["classifier__C"])
            mode = "nested_grid_search"
        else:
            fitted = build_model(columns).fit(X.iloc[train_index], y_train)
            best_c = 1.0
            mode = "fixed_C_inner_cv_not_evaluable"
        probabilities[test_index] = fitted.predict_proba(X.iloc[test_index])[:, 1]
        tuning.append({
            "fold": fold, "inner_splits": inner_splits, "selected_C": best_c,
            "mode": mode, "n_train": int(len(train_index)),
            "n_test": int(len(test_index)),
        })
    if np.isnan(probabilities).any():
        raise RuntimeError("feature-block ablation did not predict every patient exactly once")
    return probabilities, tuning


def _fusion_ablation(
    data: pd.DataFrame,
    y: pd.Series,
    columns: list[str],
    outer_cv: list[tuple[np.ndarray, np.ndarray]],
    hybrid_probability: np.ndarray,
    bootstrap_groups=None,
) -> tuple[dict | None, dict[str, np.ndarray]]:
    """Compare data-only, mechanism-only and fused models on identical outer folds."""
    mechanism_prefixes = ("twin_", "liver_cf_")
    mechanistic = [column for column in columns if column.startswith(mechanism_prefixes)]
    measured = [column for column in columns if not column.startswith(mechanism_prefixes)]
    if not mechanistic or not measured:
        return None, {}
    block_probabilities = {"hybrid": hybrid_probability}
    block_reports = {}
    for name, block in (("data_only", measured), ("mechanism_only", mechanistic)):
        probability, tuning = _nested_oof_for_columns(data, y, block, outer_cv)
        block_probabilities[name] = probability
        block_reports[name] = {
            "features": block,
            "metrics": evaluate_probabilities(y, probability, bootstrap=0),
            "hyperparameter_selection": tuning,
        }
    hybrid_metrics = evaluate_probabilities(y, hybrid_probability, bootstrap=0)
    data_metrics = block_reports["data_only"]["metrics"]
    rng = np.random.default_rng(42)
    paired_draws = {name: [] for name in (
        "roc_auc", "average_precision", "brier_score_improvement", "log_loss_improvement",
    )}
    data_probability = block_probabilities["data_only"]
    outcome = y.to_numpy(dtype=int)
    groups = None if bootstrap_groups is None else np.asarray(bootstrap_groups)
    unique_groups = np.unique(groups) if groups is not None else None
    for _ in range(FUSION_ABLATION_BOOTSTRAP_DRAWS):
        if unique_groups is None:
            index = rng.integers(0, len(outcome), len(outcome))
        else:
            sampled_groups = rng.choice(unique_groups, len(unique_groups), replace=True)
            index = np.concatenate([
                np.flatnonzero(groups == group) for group in sampled_groups
            ])
        if len(np.unique(outcome[index])) < 2:
            continue
        sampled_y = outcome[index]
        sampled_hybrid = hybrid_probability[index]
        sampled_data = data_probability[index]
        paired_draws["roc_auc"].append(
            roc_auc_score(sampled_y, sampled_hybrid)
            - roc_auc_score(sampled_y, sampled_data)
        )
        paired_draws["average_precision"].append(
            average_precision_score(sampled_y, sampled_hybrid)
            - average_precision_score(sampled_y, sampled_data)
        )
        paired_draws["brier_score_improvement"].append(
            brier_score_loss(sampled_y, sampled_data)
            - brier_score_loss(sampled_y, sampled_hybrid)
        )
        paired_draws["log_loss_improvement"].append(
            log_loss(sampled_y, sampled_data, labels=[0, 1])
            - log_loss(sampled_y, sampled_hybrid, labels=[0, 1])
        )
    paired_intervals = {
        name: {
            "lower_95": float(np.percentile(values, 2.5)),
            "upper_95": float(np.percentile(values, 97.5)),
            "valid_draws": len(values),
        }
        for name, values in paired_draws.items()
    }
    report = {
        "status": "evaluated",
        "design": "paired_nested_patient_level_cross_validation_same_outer_folds",
        "hybrid": {"features": columns, "metrics": hybrid_metrics},
        **block_reports,
        "hybrid_minus_data_only": {
            "roc_auc": float(hybrid_metrics["roc_auc"] - data_metrics["roc_auc"]),
            "average_precision": float(
                hybrid_metrics["average_precision"] - data_metrics["average_precision"]
            ),
            "brier_score_improvement": float(
                data_metrics["brier_score"] - hybrid_metrics["brier_score"]
            ),
            "log_loss_improvement": float(
                data_metrics["log_loss"] - hybrid_metrics["log_loss"]
            ),
        },
        "paired_bootstrap": {
            "draws": FUSION_ABLATION_BOOTSTRAP_DRAWS,
            "seed": 42,
            "unit": "centre" if groups is not None else "patient",
            "groups": int(len(unique_groups)) if unique_groups is not None else int(len(outcome)),
            "intervals": paired_intervals,
            "estimand": "paired_difference_of_fixed_out_of_fold_predictions",
            "limitation": (
                "Resamples whole centres with their paired predictions; does not repeat model "
                "training or capture cross-validation split uncertainty."
                if groups is not None else
                "Resamples patients and their paired predictions; does not repeat model "
                "training or capture cross-validation split uncertainty."
            ),
        },
        "interpretation": (
            "Exploratory internal ablation only; positive values do not establish clinical "
            "utility, causality, or external transportability."
        ),
    }
    spatial_features = [
        column for column in mechanistic
        if column.startswith("twin_spatial_")
        or (column.startswith("liver_cf_")
            and ("spatial_" in column or "efferent_" in column))
    ]
    without_spatial = [column for column in columns if column not in spatial_features]
    if spatial_features and without_spatial:
        reference_probability, reference_tuning = _nested_oof_for_columns(
            data, y, without_spatial, outer_cv)
        block_probabilities["without_spatial_mechanism"] = reference_probability
        reference_metrics = evaluate_probabilities(
            y, reference_probability, bootstrap=0)
        spatial_differences = {
            "roc_auc": float(hybrid_metrics["roc_auc"] - reference_metrics["roc_auc"]),
            "average_precision": float(
                hybrid_metrics["average_precision"]
                - reference_metrics["average_precision"]),
            "brier_score_improvement": float(
                reference_metrics["brier_score"] - hybrid_metrics["brier_score"]),
            "log_loss_improvement": float(
                reference_metrics["log_loss"] - hybrid_metrics["log_loss"]),
        }
        spatial_draws = {name: [] for name in spatial_differences}
        spatial_rng = np.random.default_rng(42)
        for _ in range(FUSION_ABLATION_BOOTSTRAP_DRAWS):
            if unique_groups is None:
                index = spatial_rng.integers(0, len(outcome), len(outcome))
            else:
                sampled_groups = spatial_rng.choice(
                    unique_groups, len(unique_groups), replace=True)
                index = np.concatenate([
                    np.flatnonzero(groups == group) for group in sampled_groups
                ])
            if len(np.unique(outcome[index])) < 2:
                continue
            sampled_y = outcome[index]
            sampled_hybrid = hybrid_probability[index]
            sampled_reference = reference_probability[index]
            spatial_draws["roc_auc"].append(
                roc_auc_score(sampled_y, sampled_hybrid)
                - roc_auc_score(sampled_y, sampled_reference))
            spatial_draws["average_precision"].append(
                average_precision_score(sampled_y, sampled_hybrid)
                - average_precision_score(sampled_y, sampled_reference))
            spatial_draws["brier_score_improvement"].append(
                brier_score_loss(sampled_y, sampled_reference)
                - brier_score_loss(sampled_y, sampled_hybrid))
            spatial_draws["log_loss_improvement"].append(
                log_loss(sampled_y, sampled_reference, labels=[0, 1])
                - log_loss(sampled_y, sampled_hybrid, labels=[0, 1]))
        report["spatial_mechanism_increment"] = {
            "status": "evaluated",
            "design": "paired_nested_patient_level_cross_validation_same_outer_folds",
            "spatial_features": spatial_features,
            "without_spatial_mechanism": {
                "features": without_spatial,
                "metrics": reference_metrics,
                "hyperparameter_selection": reference_tuning,
            },
            "hybrid_minus_without_spatial": spatial_differences,
            "paired_bootstrap": {
                "draws": FUSION_ABLATION_BOOTSTRAP_DRAWS,
                "seed": 42,
                "unit": "centre" if groups is not None else "patient",
                "groups": (int(len(unique_groups)) if unique_groups is not None
                           else int(len(outcome))),
                "intervals": {
                    name: {
                        "lower_95": float(np.percentile(values, 2.5)),
                        "upper_95": float(np.percentile(values, 97.5)),
                        "valid_draws": len(values),
                    }
                    for name, values in spatial_draws.items()
                },
                "estimand": "paired_difference_of_fixed_out_of_fold_predictions",
            },
            "interpretation": (
                "Exploratory spatial-mechanism increment on the development cohort; "
                "not causal or external validation."
            ),
        }
    else:
        report["spatial_mechanism_increment"] = {
            "status": "not_evaluable",
            "reason": "requires at least one explicit spatial TDLN feature and one retained reference feature",
            "spatial_features": spatial_features,
        }
    return report, block_probabilities


def _safe_group_splits(
    y: pd.Series,
    groups: pd.Series,
    maximum_splits: int,
    seed: int,
) -> tuple[list[tuple[np.ndarray, np.ndarray]] | None, int]:
    """Return grouped stratified splits only when every partition has both classes."""
    labels = np.asarray(y, dtype=int)
    group_values = np.asarray(groups)
    for n_splits in range(min(maximum_splits, int(pd.Series(labels).value_counts().min())), 1, -1):
        if len(np.unique(group_values)) == len(group_values):
            # One patient-level row per group: ordinary stratification is exactly
            # group-safe and avoids version-dependent singleton-group heuristics.
            splitter = StratifiedKFold(
                n_splits=n_splits, shuffle=True, random_state=seed,
            )
            splits = list(splitter.split(np.zeros((len(labels), 1)), labels))
        else:
            splitter = StratifiedGroupKFold(
                n_splits=n_splits, shuffle=True, random_state=seed,
            )
            splits = list(splitter.split(np.zeros((len(labels), 1)), labels, group_values))
        if all(len(np.unique(labels[train])) == 2 and len(np.unique(labels[test])) == 2
               for train, test in splits):
            return splits, n_splits
    return None, 0


def _training_domain_reference(features: pd.DataFrame, columns: list[str]) -> dict:
    """Persist label-free training support for locked-model applicability checks."""
    reference = {}
    for column in columns:
        values = pd.to_numeric(features[column], errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if not len(finite):
            raise ValueError(f"model feature has no finite training values: {column}")
        q1, median, q3 = np.percentile(finite, [25, 50, 75])
        robust_scale = float((q3 - q1) / 1.349)
        if robust_scale <= 1e-12:
            robust_scale = float(np.std(finite))
        reference[column] = {
            "median": float(median), "robust_scale": max(robust_scale, 1e-12),
            "minimum": float(np.min(finite)), "maximum": float(np.max(finite)),
            "missing_fraction": float(1.0 - len(finite) / len(values)),
        }
    return {
        "method": "robust_z_and_training_range", "features": reference,
        "thresholds": {
            "max_abs_robust_z": DOMAIN_MAX_ABS_Z,
            "max_outside_training_range_fraction": DOMAIN_MAX_OUTSIDE_FRACTION,
            "max_missingness_increase": DOMAIN_MAX_MISSINGNESS_INCREASE,
        },
        "note": "Applicability diagnostic only; not proof of cohort transportability.",
    }


def _domain_diagnostics(features: pd.DataFrame, required: list[str], reference: dict | None) -> dict:
    if not reference or set(reference.get("features", {})) != set(required):
        n = len(features)
        return {"available": np.zeros(n, dtype=bool), "max_abs_z": np.full(n, np.nan),
                "outside_fraction": np.full(n, np.nan), "missing_fraction": np.full(n, np.nan),
                "flag": np.zeros(n, dtype=bool)}
    matrix = features[required].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    missing = ~np.isfinite(matrix)
    medians = np.array([reference["features"][c]["median"] for c in required])
    scales = np.array([reference["features"][c]["robust_scale"] for c in required])
    minima = np.array([reference["features"][c]["minimum"] for c in required])
    maxima = np.array([reference["features"][c]["maximum"] for c in required])
    training_missing = np.array([reference["features"][c]["missing_fraction"] for c in required])
    filled = np.where(missing, medians, matrix)
    max_abs_z = np.max(np.abs((filled - medians) / scales), axis=1)
    outside_fraction = ((((matrix < minima) | (matrix > maxima)) & ~missing).mean(axis=1))
    missing_fraction = missing.mean(axis=1)
    thresholds = reference["thresholds"]
    flag = ((max_abs_z > thresholds["max_abs_robust_z"])
            | (outside_fraction > thresholds["max_outside_training_range_fraction"])
            | (missing_fraction > float(training_missing.mean())
               + thresholds["max_missingness_increase"]))
    return {"available": np.ones(len(features), dtype=bool), "max_abs_z": max_abs_z,
            "outside_fraction": outside_fraction, "missing_fraction": missing_fraction,
            "flag": flag}


def _spatial_evidence_reference(features: pd.DataFrame) -> dict | None:
    """Freeze development support for model-excluded spatial evidence QC."""
    support = {}
    for column in SPATIAL_EVIDENCE_COLUMNS:
        if column not in features:
            continue
        values = pd.to_numeric(features[column], errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if not len(finite):
            continue
        if ((finite < 0) | (finite > 1)).any() or len(finite) != len(values):
            raise ValueError(
                f"spatial evidence QC must be complete and within [0,1]: {column}"
            )
        support[column] = {
            "minimum": float(np.min(finite)), "maximum": float(np.max(finite)),
        }
    if not support:
        return None
    return {
        "method": "closed_development_range_for_spatial_evidence_qc",
        "features": support,
        "note": (
            "Applicability evidence only; coverage support is not a biological confidence score."
        ),
    }


def _spatial_evidence_diagnostics(features: pd.DataFrame, reference: dict | None) -> dict:
    n = len(features)
    if not isinstance(reference, dict) or not isinstance(reference.get("features"), dict):
        return {"available": np.zeros(n, dtype=bool), "outside": np.zeros(n, dtype=bool)}
    columns = list(reference["features"])
    missing = [column for column in columns if column not in features]
    if missing:
        raise ValueError(f"missing spatial evidence QC: {missing}")
    outside = np.zeros(n, dtype=bool)
    for column in columns:
        values = pd.to_numeric(features[column], errors="coerce").to_numpy(dtype=float)
        bounds = reference["features"][column]
        if (not np.isfinite(values).all() or ((values < 0) | (values > 1)).any()
                or not isinstance(bounds, dict)):
            raise ValueError(f"invalid spatial evidence QC: {column}")
        lower, upper = bounds.get("minimum"), bounds.get("maximum")
        if (not isinstance(lower, (int, float)) or not isinstance(upper, (int, float))
                or not np.isfinite([lower, upper]).all() or lower > upper):
            raise ValueError("locked model has invalid spatial evidence support")
        outside |= (values < lower - 1e-12) | (values > upper + 1e-12)
    return {"available": np.ones(n, dtype=bool), "outside": outside}


def _decision_curve(y: np.ndarray, prob: np.ndarray) -> list[dict[str, float]]:
    """Decision-curve net benefit across prespecified probability thresholds."""
    rows = []
    prevalence = float(y.mean())
    for threshold in np.arange(0.1, 1.0, 0.1):
        prediction = prob >= threshold
        tp = float(np.sum(prediction & (y == 1)))
        fp = float(np.sum(prediction & (y == 0)))
        odds = threshold / (1.0 - threshold)
        rows.append({
            "threshold": round(float(threshold), 1),
            "model_net_benefit": float(tp / len(y) - fp / len(y) * odds),
            "treat_all_net_benefit": float(prevalence - (1.0 - prevalence) * odds),
            "treat_none_net_benefit": 0.0,
        })
    return rows


def _calibration_curve(y: np.ndarray, prob: np.ndarray) -> tuple[list[dict], float]:
    """Return deterministic equal-frequency bins and their weighted absolute gap."""
    n_bins = min(5, max(2, int(np.sqrt(len(y)))))
    ordered = np.argsort(prob, kind="stable")
    rows = []
    weighted_gap = 0.0
    for number, index in enumerate(np.array_split(ordered, n_bins), start=1):
        mean_probability = float(prob[index].mean())
        observed = float(y[index].mean())
        gap = abs(mean_probability - observed)
        weighted_gap += len(index) / len(y) * gap
        rows.append({
            "bin": number,
            "n_patients": int(len(index)),
            "mean_probability": mean_probability,
            "observed_response_rate": observed,
            "absolute_gap": float(gap),
        })
    return rows, float(weighted_gap)


def evaluate_probabilities(
    y_true,
    probability,
    bootstrap: int = 500,
    seed: int = 42,
    bootstrap_groups=None,
    decision_threshold: float = 0.5,
) -> dict:
    """Discrimination, calibration and patient- or centre-bootstrap intervals."""
    raw_outcome = np.asarray(y_true)
    try:
        numeric_outcome = raw_outcome.astype(float)
    except (TypeError, ValueError) as error:
        raise ValueError("evaluation outcomes must contain only 0 and 1") from error
    if (not np.isfinite(numeric_outcome).all()
            or not np.isin(numeric_outcome, (0.0, 1.0)).all()):
        raise ValueError("evaluation outcomes must contain only 0 and 1")
    y = numeric_outcome.astype(int)
    raw_probability = np.asarray(probability, dtype=float)
    if len(y) != len(raw_probability) or len(y) == 0 or set(np.unique(y)) != {0, 1}:
        raise ValueError("evaluation requires aligned probabilities and both response classes")
    if (not np.isfinite(raw_probability).all()
            or (raw_probability < 0).any() or (raw_probability > 1).any()):
        raise ValueError("evaluation probabilities must be finite and within [0, 1]")
    prob = np.clip(raw_probability, 1e-7, 1 - 1e-7)
    if not isinstance(bootstrap, (int, np.integer)) or bootstrap < 0:
        raise ValueError("bootstrap draws must be a non-negative integer")
    if not np.isfinite(decision_threshold) or not 0 < decision_threshold < 1:
        raise ValueError("decision threshold must fall strictly between 0 and 1")
    groups = None
    if bootstrap_groups is not None:
        groups = np.asarray(bootstrap_groups)
        if len(groups) != len(y):
            raise ValueError("bootstrap groups must align one-to-one with outcomes")
        if pd.isna(groups).any() or len(np.unique(groups)) < 2:
            raise ValueError("cluster bootstrap requires at least two non-missing groups")
    pred = prob >= decision_threshold
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    logits = np.log(prob / (1 - prob)).reshape(-1, 1)
    calibrator = LogisticRegression(C=1e6, solver="lbfgs").fit(logits, y)
    calibration_curve, calibration_error = _calibration_curve(y, prob)
    point = {
        "roc_auc": float(roc_auc_score(y, prob)),
        "average_precision": float(average_precision_score(y, prob)),
        "brier_score": float(brier_score_loss(y, prob)),
        "log_loss": float(log_loss(y, prob)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "calibration_intercept": float(calibrator.intercept_[0]),
        "calibration_slope": float(calibrator.coef_[0, 0]),
        "calibration_in_the_large": float(prob.mean() - y.mean()),
        "expected_calibration_error": calibration_error,
        "calibration_curve": calibration_curve,
        "calibration_curve_design": {
            "method": "equal_frequency",
            "bins": len(calibration_curve),
            "maximum_bins": 5,
            "role": "descriptive_small_cohort_diagnostic",
        },
        "decision_threshold": float(decision_threshold),
        "decision_curve": _decision_curve(y, prob),
    }
    if bootstrap:
        rng = np.random.default_rng(seed)
        draws = {name: [] for name in (
            "roc_auc", "average_precision", "brier_score", "log_loss",
            "balanced_accuracy", "sensitivity", "specificity",
            "expected_calibration_error", "calibration_in_the_large",
        )}
        unique_groups = np.unique(groups) if groups is not None else None
        for _ in range(bootstrap):
            if unique_groups is None:
                index = rng.integers(0, len(y), len(y))
            else:
                sampled_groups = rng.choice(unique_groups, len(unique_groups), replace=True)
                index = np.concatenate([
                    np.flatnonzero(groups == group) for group in sampled_groups
                ])
            if len(np.unique(y[index])) < 2:
                continue
            draws["roc_auc"].append(roc_auc_score(y[index], prob[index]))
            draws["average_precision"].append(average_precision_score(y[index], prob[index]))
            draws["brier_score"].append(brier_score_loss(y[index], prob[index]))
            draws["log_loss"].append(log_loss(y[index], prob[index], labels=[0, 1]))
            _, sampled_calibration_error = _calibration_curve(y[index], prob[index])
            draws["expected_calibration_error"].append(sampled_calibration_error)
            draws["calibration_in_the_large"].append(
                float(prob[index].mean() - y[index].mean())
            )
            sampled_prediction = prob[index] >= decision_threshold
            sampled_tn, sampled_fp, sampled_fn, sampled_tp = confusion_matrix(
                y[index], sampled_prediction, labels=[0, 1]
            ).ravel()
            draws["balanced_accuracy"].append(
                balanced_accuracy_score(y[index], sampled_prediction)
            )
            draws["sensitivity"].append(
                sampled_tp / (sampled_tp + sampled_fn)
            )
            draws["specificity"].append(
                sampled_tn / (sampled_tn + sampled_fp)
            )
        point["bootstrap"] = {
            name: {"lower_95": float(np.percentile(values, 2.5)),
                   "upper_95": float(np.percentile(values, 97.5)),
                   "valid_draws": len(values)}
            for name, values in draws.items() if values
        }
        point["bootstrap_design"] = {
            "unit": "centre" if groups is not None else "patient",
            "requested_draws": int(bootstrap),
            "groups": int(len(unique_groups)) if unique_groups is not None else int(len(y)),
            "seed": int(seed),
            "note": (
                "Entire centres were sampled with replacement; all patients in each sampled "
                "centre were retained, including repeated centres."
                if groups is not None else
                "Patients were sampled with replacement."
            ),
        }
    return point


def _available_cluster_groups(frame: pd.DataFrame, column: str = "center"):
    """Use a declared cluster only when at least two complete groups are present."""
    if column not in frame:
        return None
    if frame[column].isna().any():
        raise ValueError(f"{column} must be complete when supplied for clustered evaluation")
    normalized = frame[column].astype(str).str.strip()
    if normalized.eq("").any():
        raise ValueError(f"{column} must be non-empty when supplied for clustered evaluation")
    return normalized if normalized.nunique() >= 2 else None


def _center_stratified_evaluation(
    data: pd.DataFrame,
    decision_threshold: float,
    center_column: str = "center",
) -> dict | None:
    """Expose centre heterogeneity without dropping centres from pooled metrics."""
    if center_column not in data:
        return None
    rows = []
    normalized_centers = data[center_column].astype(str)
    for center, group in data.groupby(normalized_centers, sort=True):
        row = {
            "center": center,
            "n_patients": int(len(group)),
            "responders": int(group.response.sum()),
            "response_prevalence": float(group.response.mean()),
        }
        if group.response.nunique() == 2:
            row.update({"status": "evaluated", "metrics": evaluate_probabilities(
                group.response, group.probability, bootstrap=0,
                decision_threshold=decision_threshold,
            )})
        else:
            row.update({
                "status": "not_evaluable",
                "reason": "centre contains only one response class",
            })
        rows.append(row)
    evaluated = [row for row in rows if row["status"] == "evaluated"]
    summary = {
        "role": "descriptive_center_heterogeneity_not_primary",
        "n_centers": len(rows),
        "n_centers_evaluable": len(evaluated),
        "centers": rows,
    }
    if evaluated:
        auc = [row["metrics"]["roc_auc"] for row in evaluated]
        brier = [row["metrics"]["brier_score"] for row in evaluated]
        summary["evaluable_center_range"] = {
            "roc_auc_min": float(min(auc)), "roc_auc_max": float(max(auc)),
            "brier_score_min": float(min(brier)), "brier_score_max": float(max(brier)),
        }
    return summary


def predict_locked(model_path: str | Path, features: pd.DataFrame) -> pd.DataFrame:
    """Apply an immutable fitted artifact to a new cohort without refitting."""
    artifact = joblib.load(model_path)
    if not isinstance(artifact, dict) or artifact.get("artifact_version") != ARTIFACT_VERSION:
        raise ValueError(f"locked prediction requires model artifact version {ARTIFACT_VERSION}")
    point_method = artifact.get("point_prediction_method")
    if point_method != "full_cohort_locked_model":
        raise ValueError("locked model lacks the version 10 point-prediction contract")
    if artifact.get("probability_estimand") != PROBABILITY_ESTIMAND:
        raise ValueError("locked model lacks the v10 probability-estimand contract")
    uncertainty_contract = validate_uncertainty_bootstrap_contract(artifact)
    features = _normalise_patient_ids(features, "features")
    contracts = artifact.get("preprocessing_contracts")
    contract_sources = artifact.get("preprocessing_contract_sources")
    if (not isinstance(contracts, dict)
            or "qc_cell_type_vocabulary_sha256" not in contracts
            or not isinstance(contract_sources, dict)):
        raise ValueError("locked model lacks preprocessing contracts")
    for column, expected in contracts.items():
        if column not in features:
            if (column == "qc_cell_type_vocabulary_sha256"
                    and contract_sources.get(column) == "fraction_column_schema_fallback"):
                fraction_schema = tuple(sorted(
                    name for name in features if name.startswith("frac_")
                ))
                observed = hashlib.sha256(
                    json.dumps(fraction_schema, separators=(",", ":")).encode("utf-8")
                ).hexdigest()
                values = [observed]
            else:
                raise ValueError(f"missing preprocessing contract: {column}")
        else:
            values = features[column].dropna().astype(str).str.strip().str.lower().unique().tolist()
        if values != [expected]:
            raise ValueError(f"preprocessing contract mismatch: {column}")
    required = artifact["features"]
    missing = set(required) - set(features.columns)
    if missing:
        raise ValueError(f"missing locked-model features: {sorted(missing)}")
    unexpected = sorted(
        column for column in features.columns
        if column != "patient_id"
        and column not in required
        and not _technical_audit_feature(column)
    )
    if unexpected:
        raise ValueError(
            "unexpected unlocked model features; external biological and mechanistic "
            f"schema must exactly match the artifact: {unexpected}"
        )
    _validate_model_matrix(features, required, "prediction")
    # The primary point estimate comes from the model fitted to the complete
    # development cohort with the locked selected C. Bootstrap refits describe
    # resampling instability; their median is a separate estimand.
    probability = artifact["model"].predict_proba(features[required])[:, 1]
    ensemble = artifact.get("bootstrap_models", [])
    if ensemble:
        ensemble_probability = np.column_stack([
            member.predict_proba(features[required])[:, 1] for member in ensemble
        ])
        bootstrap_median = np.median(ensemble_probability, axis=1)
        lower = np.percentile(ensemble_probability, 2.5, axis=1)
        upper = np.percentile(ensemble_probability, 97.5, axis=1)
        uncertainty_available = True
    else:
        bootstrap_median = np.full(len(features), np.nan)
        lower = np.full(len(features), np.nan)
        upper = np.full(len(features), np.nan)
        uncertainty_available = False
    threshold = float(artifact.get("threshold", 0.5))
    if not 0.0 < threshold < 1.0:
        raise ValueError("locked-model threshold must fall strictly between 0 and 1")
    training_hashes = set(artifact.get("training_patient_hashes", []))
    overlap_check_available = bool(training_hashes)
    overlap = features.patient_id.map(lambda value: _patient_hash(value) in training_hashes)
    training_feature_hashes = set(artifact.get("training_feature_row_hashes", []))
    feature_reuse_check_available = bool(training_feature_hashes)
    feature_reuse = features[required].apply(
        lambda row: _feature_row_hash(row, required) in training_feature_hashes, axis=1
    )
    domain = _domain_diagnostics(features, required, artifact.get("training_domain"))
    spatial_evidence = _spatial_evidence_diagnostics(
        features, artifact.get("spatial_evidence_support")
    )
    combined_domain_flag = domain["flag"] | spatial_evidence["outside"]
    return pd.DataFrame({"patient_id": features.patient_id, "probability": probability,
                         "bootstrap_probability_median": bootstrap_median,
                         "probability_lower_95": lower, "probability_upper_95": upper,
                         "probability_interval_width": upper - lower,
                         "uncertainty_available": uncertainty_available,
                         "point_prediction_method": point_method,
                         "uncertainty_method": uncertainty_contract["method"],
                         "uncertainty_bootstrap_unit": uncertainty_contract["unit"],
                         "uncertainty_bootstrap_groups": uncertainty_contract["groups"],
                         "uncertainty_bootstrap_members": uncertainty_contract["members"],
                         "predicted_response": (probability >= threshold).astype(int),
                         "decision_threshold": threshold,
                         "training_overlap": overlap,
                         "overlap_check_available": overlap_check_available,
                         "training_feature_reuse": feature_reuse,
                         "feature_reuse_check_available": feature_reuse_check_available,
                         "domain_check_available": domain["available"],
                         "domain_max_abs_robust_z": domain["max_abs_z"],
                         "domain_outside_range_fraction": domain["outside_fraction"],
                         "domain_missing_feature_fraction": domain["missing_fraction"],
                         "spatial_evidence_check_available": spatial_evidence["available"],
                         "spatial_evidence_outside_training_range": spatial_evidence["outside"],
                         "domain_shift_flag": combined_domain_flag})


def _prediction_uncertainty_contract(predictions: pd.DataFrame) -> dict:
    """Validate one v10 uncertainty contract and its row-level interval semantics."""
    required = {
        "uncertainty_method", "uncertainty_bootstrap_unit",
        "uncertainty_bootstrap_groups", "uncertainty_bootstrap_members",
        "uncertainty_available", "bootstrap_probability_median",
        "probability_lower_95", "probability_upper_95", "probability_interval_width",
    }
    missing = required.difference(predictions.columns)
    if missing:
        raise ValueError(f"external evaluation missing v10 uncertainty columns: {sorted(missing)}")
    methods = predictions.uncertainty_method.astype(str).str.strip().unique().tolist()
    units = predictions.uncertainty_bootstrap_unit.astype(str).str.strip().unique().tolist()
    expected = {"patient_bootstrap_refit": "patient",
                "centre_cluster_bootstrap_refit": "centre"}
    if len(methods) != 1 or len(units) != 1 or expected.get(methods[0]) != units[0]:
        raise ValueError("external evaluation requires one valid v10 uncertainty method and unit")
    numeric = {}
    for column in ("uncertainty_bootstrap_groups", "uncertainty_bootstrap_members"):
        values = pd.to_numeric(predictions[column], errors="coerce")
        if (not np.isfinite(values).all() or values.nunique() != 1
                or not float(values.iloc[0]).is_integer()):
            raise ValueError("external evaluation requires one integer v10 uncertainty design")
        numeric[column] = int(values.iloc[0])
    groups, members = (numeric["uncertainty_bootstrap_groups"],
                       numeric["uncertainty_bootstrap_members"])
    if groups < (2 if units[0] == "centre" else 1) or members < 0:
        raise ValueError("external evaluation requires a valid v10 uncertainty design")
    available = _boolean_audit_column(
        predictions.uncertainty_available, "uncertainty_available"
    )
    if available.nunique() != 1 or bool(available.iloc[0]) != (members > 0):
        raise ValueError("uncertainty availability conflicts with bootstrap member count")
    interval_columns = ["bootstrap_probability_median", "probability_lower_95",
                        "probability_upper_95", "probability_interval_width"]
    interval = predictions[interval_columns].apply(pd.to_numeric, errors="coerce")
    if members:
        if (not np.isfinite(interval.to_numpy()).all()
                or not interval.bootstrap_probability_median.between(0, 1).all()
                or not interval.probability_lower_95.between(0, 1).all()
                or not interval.probability_upper_95.between(0, 1).all()
                or (interval.probability_lower_95
                    > interval.bootstrap_probability_median).any()
                or (interval.bootstrap_probability_median
                    > interval.probability_upper_95).any()
                or not np.allclose(
                    interval.probability_interval_width,
                    interval.probability_upper_95 - interval.probability_lower_95,
                    rtol=1e-10, atol=1e-12,
                )):
            raise ValueError("invalid locked-model bootstrap probability intervals")
    elif not interval.isna().all().all():
        raise ValueError("zero-member uncertainty intervals must be unavailable")
    return {"method": methods[0], "unit": units[0], "groups": groups,
            "members": members}


def evaluate_external(predictions: pd.DataFrame, labels: pd.DataFrame, bootstrap: int = 500) -> dict:
    predictions, labels = _align_patients(predictions, labels)
    if "overlap_check_available" not in predictions:
        raise ValueError("external evaluation requires a model artifact with training-overlap fingerprints")
    check_available = _boolean_audit_column(
        predictions.overlap_check_available, "overlap_check_available"
    )
    if not check_available.all():
        raise ValueError("external evaluation requires a model artifact with training-overlap fingerprints")
    if "training_overlap" not in predictions:
        raise ValueError("predictions need training_overlap audit column")
    overlap_flags = _boolean_audit_column(predictions.training_overlap, "training_overlap")
    overlap = predictions.loc[overlap_flags, "patient_id"].tolist()
    if overlap:
        raise ValueError(f"external cohort overlaps model training patients: {overlap[:10]}")
    if "feature_reuse_check_available" not in predictions:
        raise ValueError("external evaluation requires training-feature reuse fingerprints")
    reuse_available = _boolean_audit_column(
        predictions.feature_reuse_check_available, "feature_reuse_check_available"
    )
    if not reuse_available.all() or "training_feature_reuse" not in predictions:
        raise ValueError("external evaluation requires training-feature reuse fingerprints")
    reuse_flags = _boolean_audit_column(
        predictions.training_feature_reuse, "training_feature_reuse"
    )
    reused = predictions.loc[reuse_flags, "patient_id"].tolist()
    if reused:
        raise ValueError(f"external cohort reuses exact training feature rows: {reused[:10]}")
    if "domain_check_available" not in predictions:
        raise ValueError("external evaluation requires locked-model applicability diagnostics")
    domain_available = _boolean_audit_column(
        predictions.domain_check_available, "domain_check_available"
    )
    if not domain_available.all() or "domain_shift_flag" not in predictions:
        raise ValueError("external evaluation requires locked-model applicability diagnostics")
    domain_flags = _boolean_audit_column(predictions.domain_shift_flag, "domain_shift_flag")
    spatial_columns = {
        "spatial_evidence_check_available",
        "spatial_evidence_outside_training_range",
    }
    supplied_spatial_columns = spatial_columns.intersection(predictions.columns)
    if supplied_spatial_columns and supplied_spatial_columns != spatial_columns:
        raise ValueError("external evaluation has incomplete spatial evidence diagnostics")
    if supplied_spatial_columns:
        spatial_available = _boolean_audit_column(
            predictions.spatial_evidence_check_available,
            "spatial_evidence_check_available",
        )
        spatial_outside = _boolean_audit_column(
            predictions.spatial_evidence_outside_training_range,
            "spatial_evidence_outside_training_range",
        )
        if spatial_available.nunique() != 1:
            raise ValueError("external evaluation requires one spatial evidence audit mode")
        spatial_checked = bool(spatial_available.iloc[0])
        if ((~spatial_available & spatial_outside).any()
                or (spatial_outside & ~domain_flags).any()):
            raise ValueError("spatial evidence diagnostics conflict with domain-shift flags")
    else:
        spatial_checked = False
        spatial_outside = pd.Series(False, index=predictions.index)
    if "decision_threshold" not in predictions:
        raise ValueError("locked predictions need decision_threshold")
    thresholds = pd.to_numeric(predictions.decision_threshold, errors="coerce")
    if (not np.isfinite(thresholds).all() or thresholds.nunique() != 1
            or not 0 < float(thresholds.iloc[0]) < 1):
        raise ValueError("locked predictions need one valid decision_threshold")
    decision_threshold = float(thresholds.iloc[0])
    if "point_prediction_method" not in predictions:
        raise ValueError("external evaluation requires a locked point-prediction method")
    point_methods = predictions.point_prediction_method.astype(str).str.strip()
    if (point_methods.nunique() != 1
            or point_methods.iloc[0] != "full_cohort_locked_model"):
        raise ValueError(
            "external evaluation requires point_prediction_method=full_cohort_locked_model"
        )
    uncertainty_contract = _prediction_uncertainty_contract(predictions)
    audited_predictions = predictions.copy()
    audited_predictions["_domain_shift_flag"] = domain_flags.to_numpy(dtype=bool)
    data = audited_predictions.merge(labels, on="patient_id", validate="one_to_one")
    centre_groups = _available_cluster_groups(data)
    metrics = evaluate_probabilities(
        data.response, data.probability, bootstrap=bootstrap,
        bootstrap_groups=centre_groups,
        decision_threshold=decision_threshold,
    )
    metrics.update({"n_patients": int(len(data)), "evaluation_type": "external_locked_model",
                    "training_overlap_checked": True, "training_overlap_count": 0,
                    "uncertainty_bootstrap_contract": uncertainty_contract,
                    "training_feature_reuse_checked": True,
                    "training_feature_reuse_count": 0,
                    "applicability_checked": True,
                    "domain_shift_count": int(domain_flags.sum()),
                    "domain_shift_fraction": float(domain_flags.mean()),
                    "in_domain_fraction": float((~domain_flags).mean()),
                    "spatial_evidence_applicability": {
                        "checked": spatial_checked,
                        "outside_count": (int(spatial_outside.sum())
                                          if spatial_checked else None),
                        "outside_fraction": (float(spatial_outside.mean())
                                             if spatial_checked else None),
                        "role": "evidence_support_not_fitted_biological_predictor",
                    },
                    "note": "External only if this cohort was never used for fitting, tuning, or feature design."})
    in_domain = data[~data._domain_shift_flag]
    class_counts = in_domain.response.value_counts()
    sensitivity = {
        "role": "prespecified_sensitivity_analysis_not_primary",
        "selection_rule": "locked artifact domain_shift_flag == false",
        "n_patients": int(len(in_domain)),
        "coverage": float(len(in_domain) / len(data)),
        "minimum_requirement": "at least 4 patients and at least 2 per response class",
    }
    if len(in_domain) >= 4 and len(class_counts) == 2 and int(class_counts.min()) >= 2:
        sensitivity["status"] = "evaluated"
        sensitivity["metrics"] = evaluate_probabilities(
            in_domain.response, in_domain.probability, bootstrap=bootstrap,
            bootstrap_groups=_available_cluster_groups(in_domain),
            decision_threshold=decision_threshold,
        )
    else:
        sensitivity["status"] = "not_evaluable"
        sensitivity["reason"] = "insufficient in-domain class support"
    metrics["in_domain_sensitivity_analysis"] = sensitivity
    center_stratified = _center_stratified_evaluation(data, decision_threshold)
    if center_stratified is not None:
        metrics["center_stratified_evaluation"] = center_stratified
    return metrics


def leave_one_site_out(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    site_column: str = "center",
    bootstrap: int = 500,
) -> tuple[dict, pd.DataFrame]:
    """Train on all other sites and predict each held-out site exactly once."""
    features, labels = _align_patients(features, labels)
    if site_column not in labels:
        raise ValueError(f"labels need site column: {site_column}")
    data = features.merge(labels[["patient_id", "response", site_column]], on="patient_id", validate="one_to_one")
    if data[site_column].isna().any() or data[site_column].nunique() < 2:
        raise ValueError("site evaluation needs at least two non-missing sites")
    cols, excluded = _select_model_features(features, labels)
    predictions = []
    skipped = []
    tuning = []
    for site_index, site in enumerate(sorted(data[site_column].astype(str).unique()), start=1):
        test = data[data[site_column].astype(str) == site]
        train = data[data[site_column].astype(str) != site]
        # A single-class held-out centre cannot support a centre-specific AUC,
        # but its patients remain valid transportability tests for proper
        # scoring rules and the pooled analysis. Dropping that centre would
        # condition evaluation on its outcomes and can make performance look
        # artificially favourable. Only an unusable training partition makes
        # prediction impossible.
        if train.response.nunique() < 2:
            skipped.append(site)
            continue
        y_train = train.response.astype(int)
        inner_cv, inner_splits = _safe_group_splits(
            y_train, train.patient_id, 3, 100 + site_index,
        )
        if inner_cv is not None:
            search = GridSearchCV(
                build_model(cols), {"classifier__C": list(DEFAULT_C_GRID)},
                scoring="neg_log_loss", cv=inner_cv, refit=True, error_score="raise",
            )
            search.fit(train[cols], y_train)
            model = search.best_estimator_
            selected_c = float(search.best_params_["classifier__C"])
            mode = "nested_grid_search"
        else:
            model = build_model(cols).fit(train[cols], y_train)
            selected_c = 1.0
            mode = "fixed_C_inner_cv_not_evaluable"
        probability = model.predict_proba(test[cols])[:, 1]
        predictions.append(pd.DataFrame({"patient_id": test.patient_id, site_column: site,
                                         "response": test.response.astype(int),
                                         "probability": probability,
                                         "selected_C": selected_c,
                                         "predicted_response": (probability >= 0.5).astype(int)}))
        tuning.append({"held_out_site": site, "inner_splits": inner_splits,
                       "selected_C": selected_c, "mode": mode,
                       "n_train": int(len(train)), "n_test": int(len(test))})
    if not predictions:
        raise ValueError("no site could be evaluated with both classes in train and test")
    out = pd.concat(predictions, ignore_index=True)
    metrics = evaluate_probabilities(
        out.response, out.probability, bootstrap=bootstrap,
        bootstrap_groups=out[site_column],
    )
    metrics.update({"n_patients": int(len(out)), "n_sites_evaluated": int(out[site_column].nunique()),
                    "skipped_sites": skipped, "site_column": site_column,
                    "evaluation_type": "leave_one_site_out",
                    "site_hyperparameter_selection": tuning,
                    "selection_metric": "negative_log_loss",
                    "C_grid": list(DEFAULT_C_GRID),
                    "feature_schema": {"included": cols, "excluded": excluded},
                    "note": "Stronger internal transportability check; not independent external validation."})
    metrics["center_stratified_evaluation"] = _center_stratified_evaluation(
        decision_threshold=0.5,
        data=out,
        center_column=site_column,
    )
    return metrics, out


def train_evaluate(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    output: str | Path,
    folds: int = 5,
    ensemble_members: int = 50,
    input_provenance: dict | None = None,
) -> dict:
    features, labels = _align_patients(features, labels)
    preprocessing_contracts = _preprocessing_contracts(features)
    preprocessing_contract_sources = {
        name: ("embedded_feature_qc" if name in features else "fraction_column_schema_fallback")
        for name in preprocessing_contracts
    }
    data = features.merge(labels, on="patient_id", validate="one_to_one")
    if data["response"].nunique() != 2:
        raise ValueError("training needs both response classes")
    if ensemble_members < 0:
        raise ValueError("ensemble_members must be non-negative")
    counts = data["response"].value_counts()
    n_splits = min(folds, int(counts.min()))
    if n_splits < 2:
        raise ValueError("each class needs at least two patients")
    cols, excluded = _select_model_features(features, labels)
    X, y = data[cols], data["response"].astype(int)
    outer_cv, n_splits = _safe_group_splits(y, data.patient_id, n_splits, 42)
    if outer_cv is None:
        raise ValueError("could not construct grouped outer folds with both classes")
    prob = np.full(len(data), np.nan)
    fold_id = np.full(len(data), -1, dtype=int)
    selected_c = np.full(len(data), np.nan)
    tuning = []
    for fold, (train_index, test_index) in enumerate(
        outer_cv, start=1
    ):
        y_train = y.iloc[train_index]
        inner_cv, inner_splits = _safe_group_splits(
            y_train, data.patient_id.iloc[train_index], 3, fold,
        )
        if inner_cv is not None:
            search = GridSearchCV(
                build_model(cols), {"classifier__C": list(DEFAULT_C_GRID)},
                scoring="neg_log_loss", cv=inner_cv, refit=True, error_score="raise",
            )
            search.fit(X.iloc[train_index], y_train)
            fitted = search.best_estimator_
            best_c = float(search.best_params_["classifier__C"])
            mode = "nested_grid_search"
        else:
            fitted = build_model(cols).fit(X.iloc[train_index], y_train)
            best_c = 1.0
            mode = "fixed_C_inner_cv_not_evaluable"
        prob[test_index] = fitted.predict_proba(X.iloc[test_index])[:, 1]
        fold_id[test_index] = fold
        selected_c[test_index] = best_c
        tuning.append({"fold": fold, "inner_splits": inner_splits, "selected_C": best_c,
                       "mode": mode, "n_train": int(len(train_index)), "n_test": int(len(test_index))})
    if np.isnan(prob).any() or (fold_id < 1).any():
        raise RuntimeError("outer cross-validation did not predict every patient exactly once")
    evaluation_groups = _available_cluster_groups(data)
    metrics = evaluate_probabilities(y, prob, bootstrap_groups=evaluation_groups)
    metrics.update({"n_patients": int(len(data)), "n_splits": n_splits,
                    "hyperparameter_selection": tuning,
                    "selection_metric": "negative_log_loss",
                    "probability_estimand": PROBABILITY_ESTIMAND,
                    "development_response_prevalence": float(y.mean()),
                    "C_grid": list(DEFAULT_C_GRID),
                    "feature_schema": {"included": cols, "excluded": excluded},
                    "evaluation_type": "internal_cross_validation",
                    "note": "Nested patient-level cross-validation; not external or clinical validation."})
    ablation, ablation_probabilities = _fusion_ablation(
        data, y, cols, outer_cv, prob, bootstrap_groups=evaluation_groups,
    )
    metrics["mechanism_data_fusion_ablation"] = (
        ablation if ablation is not None else {
            "status": "not_evaluable",
            "reason": (
                "requires at least one twin_ or liver_cf_ mechanistic feature "
                "and one measured feature"
            ),
        }
    )
    metrics["training_input_provenance"] = input_provenance
    final_cv, final_splits = _safe_group_splits(y, data.patient_id, min(5, int(counts.min())), 42)
    if final_cv is None:
        raise ValueError("could not construct grouped final-selection folds with both classes")
    final_search = GridSearchCV(
        build_model(cols), {"classifier__C": list(DEFAULT_C_GRID)},
        scoring="neg_log_loss", cv=final_cv, refit=True, error_score="raise",
    )
    final_search.fit(X, y)
    model = final_search.best_estimator_
    final_c = float(final_search.best_params_["classifier__C"])
    rng = np.random.default_rng(42)
    bootstrap_models = []
    ensemble_groups = (
        np.asarray(evaluation_groups) if evaluation_groups is not None else None
    )
    unique_ensemble_groups = (
        np.unique(ensemble_groups) if ensemble_groups is not None else None
    )
    uncertainty_method = (
        "centre_cluster_bootstrap_refit"
        if ensemble_groups is not None else "patient_bootstrap_refit"
    )
    attempts = 0
    while len(bootstrap_models) < ensemble_members and attempts < max(100, ensemble_members * 20):
        attempts += 1
        if unique_ensemble_groups is None:
            index = rng.integers(0, len(data), len(data))
        else:
            sampled_groups = rng.choice(
                unique_ensemble_groups, len(unique_ensemble_groups), replace=True,
            )
            index = np.concatenate([
                np.flatnonzero(ensemble_groups == group) for group in sampled_groups
            ])
        y_bootstrap = y.iloc[index]
        if y_bootstrap.nunique() < 2:
            continue
        member = build_model(cols, c=final_c).fit(X.iloc[index], y_bootstrap)
        bootstrap_models.append(member)
    if len(bootstrap_models) != ensemble_members:
        raise RuntimeError("could not construct requested two-class bootstrap ensemble")
    metrics["prediction_uncertainty"] = {
        "method": uncertainty_method,
        "unit": "centre" if ensemble_groups is not None else "patient",
        "groups": (
            int(len(unique_ensemble_groups))
            if unique_ensemble_groups is not None else int(len(data))
        ),
        "members": ensemble_members,
        "interval": "2.5th to 97.5th percentile",
        "fixed_selected_C": final_c,
        "point_estimate": "full_cohort_locked_model",
        "bootstrap_summary": "median_and_2.5th_to_97.5th_percentiles_reported_separately",
        "seed": 42,
        "note": "Model-resampling uncertainty only; excludes cohort shift and structural uncertainty.",
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    training_hashes = sorted(_patient_hash(value) for value in data.patient_id)
    cohort_fingerprint = hashlib.sha256("\n".join(training_hashes).encode("utf-8")).hexdigest()
    training_feature_hashes = sorted(
        _feature_row_hash(row, cols) for _, row in X.iterrows()
    )
    training_feature_cohort_fingerprint = hashlib.sha256(
        "\n".join(training_feature_hashes).encode("utf-8")
    ).hexdigest()
    training_domain = _training_domain_reference(X, cols)
    spatial_evidence_support = _spatial_evidence_reference(features)
    joblib.dump({"model": model, "features": cols, "threshold": 0.5,
                 "excluded_features": excluded,
                 "training_patients": int(len(data)), "selected_C": final_c,
                 "selection_metric": "negative_log_loss", "C_grid": list(DEFAULT_C_GRID),
                 "probability_estimand": PROBABILITY_ESTIMAND,
                 "development_response_prevalence": float(y.mean()),
                 "training_patient_hashes": training_hashes,
                 "training_cohort_fingerprint": cohort_fingerprint,
                 "training_feature_row_hashes": training_feature_hashes,
                 "training_feature_cohort_fingerprint": training_feature_cohort_fingerprint,
                 "training_domain": training_domain,
                 "spatial_evidence_support": spatial_evidence_support,
                 "preprocessing_contracts": preprocessing_contracts,
                 "preprocessing_contract_sources": preprocessing_contract_sources,
                 "training_input_provenance": input_provenance,
                 "bootstrap_models": bootstrap_models,
                 "uncertainty_method": uncertainty_method,
                 "uncertainty_bootstrap_design": {
                     "unit": "centre" if ensemble_groups is not None else "patient",
                     "groups": (
                         int(len(unique_ensemble_groups))
                         if unique_ensemble_groups is not None else int(len(data))
                     ),
                     "members": int(ensemble_members),
                     "seed": 42,
                 },
                 "point_prediction_method": "full_cohort_locked_model",
                 "artifact_version": ARTIFACT_VERSION}, output / "model.joblib")
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    prediction_columns = {
        "patient_id": data["patient_id"], "response": y, "probability": prob,
        "outer_fold": fold_id, "fold_selected_C": selected_c,
    }
    if ablation_probabilities:
        prediction_columns["data_only_probability"] = ablation_probabilities["data_only"]
        prediction_columns["mechanism_only_probability"] = ablation_probabilities[
            "mechanism_only"
        ]
        if "without_spatial_mechanism" in ablation_probabilities:
            prediction_columns["without_spatial_mechanism_probability"] = (
                ablation_probabilities["without_spatial_mechanism"])
    pd.DataFrame(prediction_columns).to_csv(
        output / "cv_predictions.csv", index=False
    )
    return metrics
