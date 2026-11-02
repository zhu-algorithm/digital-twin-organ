import hashlib
import json

import pandas as pd
import pytest

from tdln_pd1.augmentation_audit import (
    verify_augmentation_audit, write_augmentation_audit,
)
from tdln_pd1.hybrid import augment_with_twin
from tdln_pd1.spatial_dynamics import SPECIES, ZONE_NAMES


SETTINGS = {
    "days": 1.0, "dose_days": [0.0], "dose": 1.0, "labels_used": False,
    "uncertainty_samples": 0, "uncertainty_seed": 42,
    "cross_layer_sampling": "independent_systemic_spatial_mapping_no_joint_covariance",
}


def _spatial_states() -> pd.DataFrame:
    return pd.DataFrame([
        {"patient_id": "p1", "species": species, "zone": zone, "amount": 0.0}
        for species in SPECIES for zone in ZONE_NAMES
    ])


def _case(tmp_path, spatial=False, feature_derivation=False, spatial_derivation=False):
    features = tmp_path / "features.csv"
    feature_frame = pd.DataFrame({"patient_id": ["p1"], "sig_tcf7_memory": [0.5]})
    feature_frame.to_csv(features, index=False)
    inputs = {"features": {
        "file": features.name,
        "sha256": hashlib.sha256(features.read_bytes()).hexdigest(),
    }}
    if feature_derivation:
        sidecar = tmp_path / "features.csv.feature-audit.json"
        sidecar.write_text("{\"stage\": \"label_free_single_cell_featurization\"}\n")
        inputs["features"]["derivation_audit"] = {
            "file": sidecar.name,
            "sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
        }
    states = None
    paths = {"features": features}
    if spatial:
        states = tmp_path / "states.csv"
        state_frame = _spatial_states()
        state_frame.to_csv(states, index=False)
        inputs["spatial_initial_zones"] = {
            "file": states.name,
            "sha256": hashlib.sha256(states.read_bytes()).hexdigest(),
        }
        if spatial_derivation:
            sidecar = tmp_path / "states.csv.spatial-audit.json"
            sidecar.write_text("{\"stage\": \"label_free_spatial_state_featurization\"}\n")
            inputs["spatial_initial_zones"]["derivation_audit"] = {
                "file": sidecar.name,
                "sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
            }
        paths["spatial_initial_zones"] = states
    output = tmp_path / "hybrid.csv"
    augment_with_twin(
        feature_frame, days=1, dose_days=(0,),
        spatial_states=pd.read_csv(states) if states else None,
    ).to_csv(output, index=False)
    audit = write_augmentation_audit(output, inputs=inputs, settings=SETTINGS)
    return features, states, output, audit, paths


def test_augmentation_audit_binds_all_inputs_replays_and_output(tmp_path):
    _, _, output, audit, paths = _case(tmp_path, spatial=True)
    verified = verify_augmentation_audit(audit, output, paths)
    assert verified["valid"] is True
    assert verified["derivation_replayed"] is True
    assert verified["verified_inputs"] == ["features", "spatial_initial_zones"]
    assert len(verified["output_sha256"]) == 64


def test_augmentation_audit_detects_input_and_output_replacement(tmp_path):
    features, _, output, audit, paths = _case(tmp_path)
    original_features = features.read_text()
    features.write_text("patient_id,sig_tcf7_memory\np1,999\n")
    with pytest.raises(ValueError, match="input SHA-256 mismatch"):
        verify_augmentation_audit(audit, output, paths)
    features.write_text(original_features)
    output.write_text("patient_id,x\np1,999\n")
    with pytest.raises(ValueError, match="output SHA-256 mismatch"):
        verify_augmentation_audit(audit, output, paths)


def test_augmentation_audit_binds_upstream_feature_derivation_audit(tmp_path):
    _, _, output, audit, paths = _case(tmp_path, feature_derivation=True)
    assert verify_augmentation_audit(audit, output, paths)["valid"] is True
    derivation = tmp_path / "features.csv.feature-audit.json"
    derivation.write_text("replaced\n")
    with pytest.raises(ValueError, match="derivation audit SHA-256 mismatch"):
        verify_augmentation_audit(audit, output, paths)


def test_augmentation_audit_binds_upstream_spatial_derivation_audit(tmp_path):
    _, _, output, audit, paths = _case(
        tmp_path, spatial=True, spatial_derivation=True
    )
    assert verify_augmentation_audit(audit, output, paths)["valid"] is True
    derivation = tmp_path / "states.csv.spatial-audit.json"
    derivation.write_text("replaced\n")
    with pytest.raises(ValueError, match="spatial_initial_zones derivation"):
        verify_augmentation_audit(audit, output, paths)


def test_augmentation_replay_detects_setting_tampering(tmp_path):
    _, _, output, audit, paths = _case(tmp_path)
    payload = json.loads(audit.read_text())
    payload["settings"]["dose"] = 0.0
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="augmentation derivation replay mismatch"):
        verify_augmentation_audit(audit, output, paths)
