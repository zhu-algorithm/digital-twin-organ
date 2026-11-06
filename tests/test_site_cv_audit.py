import json

import pandas as pd
import pytest

from tdln_pd1.model import leave_one_site_out
from tdln_pd1.site_cv_audit import verify_site_cv_audit, write_site_cv_audit


def _run(tmp_path):
    features = tmp_path / "features.csv"
    labels = tmp_path / "labels.csv"
    output = tmp_path / "site_cv"
    output.mkdir()
    pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(12)],
        "x": [0.1, 0.9, 0.2, 0.8] * 3,
    }).to_csv(features, index=False)
    pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(12)],
        "response": [0, 1, 0, 1] * 3,
        "center": ["a"] * 4 + ["b"] * 4 + ["c"] * 4,
    }).to_csv(labels, index=False)
    metrics, predictions = leave_one_site_out(
        pd.read_csv(features), pd.read_csv(labels), bootstrap=10,
    )
    predictions_path = output / "site_cv_predictions.csv"
    metrics_path = output / "site_cv_metrics.json"
    predictions.to_csv(predictions_path, index=False)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    audit = write_site_cv_audit(
        features, labels, output, site_column="center", bootstrap=10,
    )
    return features, labels, metrics_path, predictions_path, audit


def test_site_cv_audit_replays_complete_validation(tmp_path):
    features, labels, metrics, predictions, audit = _run(tmp_path)
    report = verify_site_cv_audit(audit, features, labels, metrics, predictions)
    assert report["valid"] is True
    assert report["derivation_replayed"] is True
    assert report["patients_replayed"] == 12
    assert report["site_column"] == "center"
    assert report["bootstrap_draws"] == 10


@pytest.mark.parametrize("target", ["features", "labels", "metrics", "predictions"])
def test_site_cv_audit_rejects_any_bound_file_replacement(tmp_path, target):
    features, labels, metrics, predictions, audit = _run(tmp_path)
    paths = {"features": features, "labels": labels,
             "metrics": metrics, "predictions": predictions}
    paths[target].write_text(paths[target].read_text() + "\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_site_cv_audit(audit, features, labels, metrics, predictions)


def test_site_cv_audit_rejects_forged_output_even_with_updated_hash(tmp_path):
    features, labels, metrics, predictions, audit = _run(tmp_path)
    payload = json.loads(audit.read_text())
    report = json.loads(metrics.read_text())
    report["roc_auc"] = 0.123
    metrics.write_text(json.dumps(report, indent=2) + "\n")
    from tdln_pd1.provenance import sha256_file
    payload["outputs"]["metrics"]["sha256"] = sha256_file(metrics)
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="metric derivation replay"):
        verify_site_cv_audit(audit, features, labels, metrics, predictions)
