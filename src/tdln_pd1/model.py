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
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
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
    metrics = {
        "n_patients": int(len(data)),
        "n_splits": n_splits,
        "roc_auc": float(roc_auc_score(y, prob)),
        "average_precision": float(average_precision_score(y, prob)),
        "balanced_accuracy": float(balanced_accuracy_score(y, prob >= 0.5)),
        "note": "Cross-validated baseline metrics; not a clinical validation.",
    }
    model.fit(X, y)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "features": cols}, output / "model.joblib")
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    pd.DataFrame({"patient_id": data["patient_id"], "response": y, "probability": prob}).to_csv(output / "cv_predictions.csv", index=False)
    return metrics
