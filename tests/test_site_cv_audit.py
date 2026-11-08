import hashlib
import json

import pandas as pd
import pytest

from tdln_pd1.augmentation_audit import write_augmentation_audit
from tdln_pd1.hybrid import augment_with_twin
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


def test_mechanistic_site_cv_audit_binds_twin_augmentation_audit(tmp_path):
    features, labels, metrics, predictions, _ = _run(tmp_path)
    raw_features = tmp_path / "raw_features.csv"
    pd.read_csv(features).to_csv(raw_features, index=False)
    augment_with_twin(
        pd.read_csv(raw_features), days=1, dose_days=(0,),
    ).to_csv(features, index=False)
    twin_audit = write_augmentation_audit(
        features,
        inputs={"features": {
            "file": raw_features.name,
            "sha256": hashlib.sha256(raw_features.read_bytes()).hexdigest(),
        }},
        settings={
            "days": 1, "dose": 1, "dose_days": [0], "uncertainty_samples": 0,
            "uncertainty_seed": 42, "labels_used": False,
            "cross_layer_sampling": (
                "independent_systemic_spatial_mapping_no_joint_covariance"
            ),
            "spatial_missing_zone_policy": "feature_prior",
        },
    )
    augmentation_inputs = {"features": raw_features}
    output = metrics.parent
    replayed_metrics, replayed_predictions = leave_one_site_out(
        pd.read_csv(features), pd.read_csv(labels), bootstrap=10,
    )
    predictions.write_text(replayed_predictions.to_csv(index=False))
    metrics.write_text(json.dumps(replayed_metrics, indent=2) + "\n")
    audit = write_site_cv_audit(
        features, labels, output, site_column="center", bootstrap=10,
        augmentation_input_paths=augmentation_inputs,
    )
    report = verify_site_cv_audit(
        audit, features, labels, metrics, predictions, augmentation_inputs,
    )
    assert report["augmentation_audit_verified"] is True
    assert report["augmentation_derivation_replayed"] is True
    original_raw = raw_features.read_text()
    raw_features.write_text(original_raw + "\n")
    with pytest.raises(ValueError, match="augmentation input SHA-256 mismatch"):
        verify_site_cv_audit(
            audit, features, labels, metrics, predictions, augmentation_inputs,
        )
    raw_features.write_text(original_raw)
    twin_audit.unlink()
    with pytest.raises(ValueError, match="augmentation audit SHA-256 mismatch"):
        verify_site_cv_audit(
            audit, features, labels, metrics, predictions, augmentation_inputs,
        )
    twin_audit = write_augmentation_audit(
        features,
        inputs={"features": {
            "file": raw_features.name,
            "sha256": hashlib.sha256(raw_features.read_bytes()).hexdigest(),
        }},
        settings={
            "days": 1, "dose": 1, "dose_days": [0], "uncertainty_samples": 0,
            "uncertainty_seed": 42, "labels_used": False,
            "cross_layer_sampling": (
                "independent_systemic_spatial_mapping_no_joint_covariance"
            ),
            "spatial_missing_zone_policy": "feature_prior",
        },
    )
    payload = json.loads(audit.read_text())
    payload["inputs"]["augmentation_audit"] = None
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="requires a bound augmentation audit"):
        verify_site_cv_audit(
            audit, features, labels, metrics, predictions, augmentation_inputs,
        )


def test_mechanistic_site_cv_requires_full_augmentation_replay_inputs(tmp_path):
    features, labels, metrics, predictions, _ = _run(tmp_path)
    raw_features = tmp_path / "raw_features.csv"
    pd.read_csv(features).to_csv(raw_features, index=False)
    augment_with_twin(pd.read_csv(raw_features), days=1, dose_days=(0,)).to_csv(
        features, index=False,
    )
    write_augmentation_audit(
        features,
        inputs={"features": {
            "file": raw_features.name,
            "sha256": hashlib.sha256(raw_features.read_bytes()).hexdigest(),
        }},
        settings={
            "days": 1, "dose": 1, "dose_days": [0], "uncertainty_samples": 0,
            "uncertainty_seed": 42, "labels_used": False,
            "cross_layer_sampling": (
                "independent_systemic_spatial_mapping_no_joint_covariance"
            ),
            "spatial_missing_zone_policy": "feature_prior",
        },
    )
    replayed_metrics, replayed_predictions = leave_one_site_out(
        pd.read_csv(features), pd.read_csv(labels), bootstrap=10,
    )
    predictions.write_text(replayed_predictions.to_csv(index=False))
    metrics.write_text(json.dumps(replayed_metrics, indent=2) + "\n")
    with pytest.raises(ValueError, match="requires input path: features"):
        write_site_cv_audit(
            features, labels, metrics.parent, site_column="center", bootstrap=10,
        )
