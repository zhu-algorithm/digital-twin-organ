"""Training and evaluation without cell-level leakage."""

from __future__ import annotations

import json
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
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def build_model(feature_columns: list[str]) -> Pipeline:
    prep = ColumnTransformer([
        ("numeric", Pipeline([
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
        ]), feature_columns)
    ])
    return Pipeline([("prep", prep), ("classifier", LogisticRegression(max_iter=2000, class_weight="balanced"))])


def evaluate_probabilities(y_true, probability, bootstrap: int = 500, seed: int = 42) -> dict:
    """Discrimination, calibration and patient-level bootstrap intervals."""
    y = np.asarray(y_true, dtype=int)
    prob = np.clip(np.asarray(probability, dtype=float), 1e-7, 1 - 1e-7)
    if len(y) != len(prob) or len(y) == 0 or set(np.unique(y)) != {0, 1}:
        raise ValueError("evaluation requires aligned probabilities and both response classes")
    pred = prob >= 0.5
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
    }
    if bootstrap:
        rng = np.random.default_rng(seed)
        draws = {name: [] for name in ("roc_auc", "average_precision", "brier_score")}
        for _ in range(bootstrap):
            index = rng.integers(0, len(y), len(y))
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
    return point


def predict_locked(model_path: str | Path, features: pd.DataFrame) -> pd.DataFrame:
    """Apply an immutable fitted artifact to a new cohort without refitting."""
    artifact = joblib.load(model_path)
    required = artifact["features"]
    missing = set(required) - set(features.columns)
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    if missing:
        raise ValueError(f"missing locked-model features: {sorted(missing)}")
    probability = artifact["model"].predict_proba(features[required])[:, 1]
    return pd.DataFrame({"patient_id": features.patient_id, "probability": probability,
                         "predicted_response": (probability >= 0.5).astype(int)})


def evaluate_external(predictions: pd.DataFrame, labels: pd.DataFrame, bootstrap: int = 500) -> dict:
    data = predictions.merge(labels, on="patient_id", validate="one_to_one")
    metrics = evaluate_probabilities(data.response, data.probability, bootstrap=bootstrap)
    metrics.update({"n_patients": int(len(data)), "evaluation_type": "external_locked_model",
                    "note": "External only if this cohort was never used for fitting, tuning, or feature design."})
    return metrics


def leave_one_site_out(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    site_column: str = "center",
    bootstrap: int = 500,
) -> tuple[dict, pd.DataFrame]:
    """Train on all other sites and predict each held-out site exactly once."""
    if site_column not in labels:
        raise ValueError(f"labels need site column: {site_column}")
    data = features.merge(labels[["patient_id", "response", site_column]], on="patient_id", validate="one_to_one")
    if data[site_column].isna().any() or data[site_column].nunique() < 2:
        raise ValueError("site evaluation needs at least two non-missing sites")
    cols = [column for column in features.columns if column != "patient_id"]
    predictions = []
    skipped = []
    for site in sorted(data[site_column].astype(str).unique()):
        test = data[data[site_column].astype(str) == site]
        train = data[data[site_column].astype(str) != site]
        if train.response.nunique() < 2 or test.response.nunique() < 2:
            skipped.append(site)
            continue
        model = build_model(cols).fit(train[cols], train.response.astype(int))
        probability = model.predict_proba(test[cols])[:, 1]
        predictions.append(pd.DataFrame({"patient_id": test.patient_id, site_column: site,
                                         "response": test.response.astype(int),
                                         "probability": probability,
                                         "predicted_response": (probability >= 0.5).astype(int)}))
    if not predictions:
        raise ValueError("no site could be evaluated with both classes in train and test")
    out = pd.concat(predictions, ignore_index=True)
    metrics = evaluate_probabilities(out.response, out.probability, bootstrap=bootstrap)
    metrics.update({"n_patients": int(len(out)), "n_sites_evaluated": int(out[site_column].nunique()),
                    "skipped_sites": skipped, "site_column": site_column,
                    "evaluation_type": "leave_one_site_out",
                    "note": "Stronger internal transportability check; not independent external validation."})
    return metrics, out


def train_evaluate(features: pd.DataFrame, labels: pd.DataFrame, output: str | Path, folds: int = 5) -> dict:
    data = features.merge(labels, on="patient_id", validate="one_to_one")
    if data["response"].nunique() != 2:
        raise ValueError("training needs both response classes")
    counts = data["response"].value_counts()
    n_splits = min(folds, int(counts.min()))
    if n_splits < 2:
        raise ValueError("each class needs at least two patients")
    cols = [c for c in data.columns if c not in {"patient_id", "response"}]
    X, y = data[cols], data["response"].astype(int)
    model = build_model(cols)
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    prob = cross_val_predict(model, X, y, groups=data["patient_id"], cv=cv, method="predict_proba")[:, 1]
    metrics = evaluate_probabilities(y, prob)
    metrics.update({"n_patients": int(len(data)), "n_splits": n_splits,
                    "evaluation_type": "internal_cross_validation",
                    "note": "Cross-validated baseline metrics; not external or clinical validation."})
    model.fit(X, y)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "features": cols, "threshold": 0.5,
                 "training_patients": int(len(data)), "artifact_version": 1}, output / "model.joblib")
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame({"patient_id": data["patient_id"], "response": y, "probability": prob}).to_csv(output / "cv_predictions.csv", index=False)
    return metrics
