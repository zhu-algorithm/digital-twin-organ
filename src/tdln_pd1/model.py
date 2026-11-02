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
FORBIDDEN_FEATURE_TOKENS = (
    "response", "responder", "outcome", "label", "recist",
    "progression_free", "overall_survival", "pfs", "os_months",
)
DOMAIN_MAX_ABS_Z = 4.0
DOMAIN_MAX_OUTSIDE_FRACTION = 0.2
DOMAIN_MAX_MISSINGNESS_INCREASE = 0.2


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
    non_numeric = [column for column in candidates
                   if not pd.api.types.is_numeric_dtype(features[column])]
    if non_numeric:
        raise ValueError(
            "model features must be numeric; keep metadata in the label/cohort table: "
            f"{sorted(non_numeric)}"
        )
    excluded = {
        column: "technical cell-count/QC feature"
        for column in candidates
        if column == "n_cells" or column.startswith("n_cells_") or column.startswith("qc_")
    }
    included = [column for column in candidates if column not in excluded]
    if not included:
        raise ValueError("no eligible biological or mechanistic model features remain")
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


def build_model(feature_columns: list[str], c: float = 1.0) -> Pipeline:
    prep = ColumnTransformer([
        ("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
        ]), feature_columns)
    ])
    return Pipeline([("prep", prep), ("classifier", LogisticRegression(
        C=c, max_iter=2000, class_weight="balanced", solver="liblinear", random_state=42,
    ))])


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


def evaluate_probabilities(
    y_true,
    probability,
    bootstrap: int = 500,
    seed: int = 42,
    bootstrap_groups=None,
    decision_threshold: float = 0.5,
) -> dict:
    """Discrimination, calibration and patient- or centre-bootstrap intervals."""
    y = np.asarray(y_true, dtype=int)
    prob = np.clip(np.asarray(probability, dtype=float), 1e-7, 1 - 1e-7)
    if len(y) != len(prob) or len(y) == 0 or set(np.unique(y)) != {0, 1}:
        raise ValueError("evaluation requires aligned probabilities and both response classes")
    if not np.isfinite(prob).all():
        raise ValueError("evaluation probabilities must be finite")
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
        "decision_threshold": float(decision_threshold),
        "decision_curve": _decision_curve(y, prob),
    }
    if bootstrap:
        rng = np.random.default_rng(seed)
        draws = {name: [] for name in ("roc_auc", "average_precision", "brier_score")}
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
    return frame[column] if frame[column].nunique() >= 2 else None


def predict_locked(model_path: str | Path, features: pd.DataFrame) -> pd.DataFrame:
    """Apply an immutable fitted artifact to a new cohort without refitting."""
    artifact = joblib.load(model_path)
    features = _normalise_patient_ids(features, "features")
    required = artifact["features"]
    missing = set(required) - set(features.columns)
    if missing:
        raise ValueError(f"missing locked-model features: {sorted(missing)}")
    ensemble = artifact.get("bootstrap_models", [])
    if ensemble:
        ensemble_probability = np.column_stack([
            member.predict_proba(features[required])[:, 1] for member in ensemble
        ])
        probability = np.median(ensemble_probability, axis=1)
        lower = np.percentile(ensemble_probability, 2.5, axis=1)
        upper = np.percentile(ensemble_probability, 97.5, axis=1)
        uncertainty_available = True
    else:
        probability = artifact["model"].predict_proba(features[required])[:, 1]
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
    return pd.DataFrame({"patient_id": features.patient_id, "probability": probability,
                         "probability_lower_95": lower, "probability_upper_95": upper,
                         "probability_interval_width": upper - lower,
                         "uncertainty_available": uncertainty_available,
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
                         "domain_shift_flag": domain["flag"]})


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
    if "decision_threshold" not in predictions:
        raise ValueError("locked predictions need decision_threshold")
    thresholds = pd.to_numeric(predictions.decision_threshold, errors="coerce")
    if (not np.isfinite(thresholds).all() or thresholds.nunique() != 1
            or not 0 < float(thresholds.iloc[0]) < 1):
        raise ValueError("locked predictions need one valid decision_threshold")
    decision_threshold = float(thresholds.iloc[0])
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
                    "training_feature_reuse_checked": True,
                    "training_feature_reuse_count": 0,
                    "applicability_checked": True,
                    "domain_shift_count": int(domain_flags.sum()),
                    "domain_shift_fraction": float(domain_flags.mean()),
                    "in_domain_fraction": float((~domain_flags).mean()),
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
        if train.response.nunique() < 2 or test.response.nunique() < 2:
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
    metrics = evaluate_probabilities(y, prob)
    metrics.update({"n_patients": int(len(data)), "n_splits": n_splits,
                    "hyperparameter_selection": tuning,
                    "selection_metric": "negative_log_loss",
                    "C_grid": list(DEFAULT_C_GRID),
                    "feature_schema": {"included": cols, "excluded": excluded},
                    "evaluation_type": "internal_cross_validation",
                    "note": "Nested patient-level cross-validation; not external or clinical validation."})
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
    attempts = 0
    while len(bootstrap_models) < ensemble_members and attempts < max(100, ensemble_members * 20):
        attempts += 1
        index = rng.integers(0, len(data), len(data))
        y_bootstrap = y.iloc[index]
        if y_bootstrap.nunique() < 2:
            continue
        member = build_model(cols, c=final_c).fit(X.iloc[index], y_bootstrap)
        bootstrap_models.append(member)
    if len(bootstrap_models) != ensemble_members:
        raise RuntimeError("could not construct requested two-class patient bootstrap ensemble")
    metrics["prediction_uncertainty"] = {
        "method": "patient_bootstrap_refit",
        "members": ensemble_members,
        "interval": "2.5th to 97.5th percentile",
        "fixed_selected_C": final_c,
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
    joblib.dump({"model": model, "features": cols, "threshold": 0.5,
                 "excluded_features": excluded,
                 "training_patients": int(len(data)), "selected_C": final_c,
                 "selection_metric": "negative_log_loss", "C_grid": list(DEFAULT_C_GRID),
                 "training_patient_hashes": training_hashes,
                 "training_cohort_fingerprint": cohort_fingerprint,
                 "training_feature_row_hashes": training_feature_hashes,
                 "training_feature_cohort_fingerprint": training_feature_cohort_fingerprint,
                 "training_domain": training_domain,
                 "training_input_provenance": input_provenance,
                 "bootstrap_models": bootstrap_models,
                 "uncertainty_method": "patient_bootstrap_refit",
                 "artifact_version": 7}, output / "model.joblib")
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame({"patient_id": data["patient_id"], "response": y, "probability": prob,
                  "outer_fold": fold_id, "fold_selected_C": selected_c}).to_csv(
        output / "cv_predictions.csv", index=False
    )
    return metrics
