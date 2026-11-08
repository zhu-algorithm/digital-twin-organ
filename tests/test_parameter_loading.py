from dataclasses import asdict
import hashlib
import json

import pandas as pd
import pytest

from tdln_pd1.augmentation_audit import write_augmentation_audit
from tdln_pd1.cli import (_load_liver_replay_manifest, _load_parameter_file,
                          _training_input_provenance)
from tdln_pd1.hybrid import FeatureStateMapping, augment_with_twin
from tdln_pd1.simulator import TwinParameters


def _audited_twin_features(tmp_path):
    raw = tmp_path / "raw.csv"
    features = tmp_path / "hybrid.csv"
    pd.DataFrame({"patient_id": ["p0"], "x": [1.0]}).to_csv(raw, index=False)
    augment_with_twin(pd.read_csv(raw), days=1, dose_days=(0,)).to_csv(
        features, index=False,
    )
    write_augmentation_audit(
        features,
        inputs={"features": {
            "file": raw.name,
            "sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
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
    return features, raw, {"features": raw}


def test_parameter_file_loading_is_exact_and_fingerprinted(tmp_path):
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(asdict(TwinParameters())))
    parameters, provenance = _load_parameter_file(str(path), TwinParameters)
    assert parameters == TwinParameters()
    assert len(provenance["sha256"]) == 64
    assert provenance["parameter_class"] == "TwinParameters"


def test_liver_replay_manifest_resolves_relative_paths_and_rejects_unknown_keys(tmp_path):
    manifest = tmp_path / "replay.json"
    manifest.write_text(json.dumps({
        "features": "raw/features.csv",
        "treatment_manifest": "/locked/treatment.csv",
    }))
    loaded = _load_liver_replay_manifest(manifest)
    assert loaded["features"] == str(tmp_path / "raw/features.csv")
    assert loaded["treatment_manifest"] == "/locked/treatment.csv"
    manifest.write_text(json.dumps({"features": "raw.csv", "labels": "labels.csv"}))
    with pytest.raises(ValueError, match="manifest schema"):
        _load_liver_replay_manifest(manifest)


def test_parameter_file_rejects_schema_drift_and_nonpositive_values(tmp_path):
    payload = asdict(TwinParameters())
    payload.pop("tumour_growth")
    path = tmp_path / "missing.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="schema mismatch"):
        _load_parameter_file(str(path), TwinParameters)
    payload = asdict(TwinParameters())
    payload["tumour_growth"] = 0
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="strictly positive"):
        _load_parameter_file(str(path), TwinParameters)


def test_feature_state_mapping_file_is_exact_and_fingerprinted(tmp_path):
    path = tmp_path / "mapping.json"
    path.write_text(json.dumps(asdict(FeatureStateMapping())))
    mapping, provenance = _load_parameter_file(str(path), FeatureStateMapping)
    assert mapping == FeatureStateMapping()
    assert provenance["parameter_class"] == "FeatureStateMapping"
    assert len(provenance["sha256"]) == 64


def test_training_input_binds_feature_and_parameter_sidecar_hashes(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,x\np0,1\n")
    sidecar = tmp_path / "hybrid.csv.parameters.json"
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "systemic": {"file": "fit.json", "sha256": "a" * 64,
                     "parameter_class": "TwinParameters"},
    }, "research_only": True}))
    provenance = _training_input_provenance(str(features))
    assert len(provenance["feature_sha256"]) == 64
    assert len(provenance["parameter_sidecar_sha256"]) == 64
    assert provenance["parameter_provenance"]["sources"]["systemic"]["sha256"] == "a" * 64


def test_training_input_accepts_and_validates_mapping_provenance(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,x\np0,1\n")
    sidecar = tmp_path / "hybrid.csv.parameters.json"
    source = {"file": "mapping.json", "sha256": "b" * 64,
              "parameter_class": "FeatureStateMapping"}
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "feature_state_mapping": source,
    }, "research_only": True}))
    provenance = _training_input_provenance(str(features))
    assert provenance["parameter_provenance"]["sources"]["feature_state_mapping"] == source
    source["parameter_class"] = "TwinParameters"
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "feature_state_mapping": source,
    }, "research_only": True}))
    with pytest.raises(ValueError, match="invalid parameter source fingerprint"):
        _training_input_provenance(str(features))


def test_training_input_binds_and_checks_label_free_augmentation_audit(tmp_path):
    features, _, inputs = _audited_twin_features(tmp_path)
    provenance = _training_input_provenance(str(features), inputs)
    assert provenance["augmentation_audit"]["labels_used"] is False
    assert len(provenance["augmentation_audit_sha256"]) == 64
    assert provenance["augmentation_derivation_replayed"] is True
    assert provenance["augmentation_verified_inputs"] == ["features"]

    features.write_text(features.read_text() + "\n")
    with pytest.raises(ValueError, match="mismatched feature augmentation audit"):
        _training_input_provenance(str(features), inputs)


def test_mechanistic_training_features_require_augmentation_lineage(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,sig_x,twin_effector_gain\np0,1,2\n")
    with pytest.raises(ValueError, match="mechanistic training features require"):
        _training_input_provenance(str(features))

    features, _, inputs = _audited_twin_features(tmp_path)
    provenance = _training_input_provenance(str(features), inputs)
    assert provenance["mechanistic_feature_lineage_required"] is True


def test_training_rejects_claimed_augmentation_without_replay_inputs(tmp_path):
    features, _, _ = _audited_twin_features(tmp_path)
    with pytest.raises(ValueError, match="requires input path: features"):
        _training_input_provenance(str(features))
