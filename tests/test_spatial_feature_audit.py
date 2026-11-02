import json

import pandas as pd
import pytest

from tdln_pd1.spatial_feature_audit import (
    verify_spatial_feature_audit, write_spatial_feature_audit,
)
from tdln_pd1.spatial_data import aggregate_spatial_spots


def _write_valid_spots_and_states(spots, output, zone="paracortex"):
    frame = pd.DataFrame({
        "patient_id": ["p1"], "spot_id": ["s1"], "zone": [zone],
        "sampling_day": [0], "area": [2.0], "antigen": [0.0],
        "mature_dc": [1.0], "tpex": [2.0], "effector": [3.0],
    })
    frame.to_csv(spots, index=False)
    aggregate_spatial_spots(frame, "area", "sampling_day", ("0",)).to_csv(
        output, index=False
    )


def test_spatial_feature_audit_binds_spots_weight_and_state_output(tmp_path):
    spots = tmp_path / "spots.csv"
    output = tmp_path / "states.csv"
    _write_valid_spots_and_states(spots, output)
    audit = write_spatial_feature_audit(
        spots, output, "area", "sampling_day", ("0",)
    )
    verified = verify_spatial_feature_audit(audit, spots, output)
    assert verified["valid"] is True
    assert verified["payload"]["settings"]["weight_column"] == "area"
    assert verified["payload"]["settings"]["time_column"] == "sampling_day"
    assert verified["payload"]["settings"]["baseline_values"] == ["0"]
    assert verified["payload"]["labels_used"] is False
    assert verified["derivation_replayed"] is True
    assert len(verified["audit_sha256"]) == 64


def test_spatial_feature_audit_detects_replacement_and_label_drift(tmp_path):
    spots = tmp_path / "spots.csv"
    output = tmp_path / "states.csv"
    _write_valid_spots_and_states(spots, output)
    audit = write_spatial_feature_audit(spots, output, "area", "sampling_day", ("0",))
    _write_valid_spots_and_states(spots, tmp_path / "replacement.csv", zone="medulla")
    with pytest.raises(ValueError, match="spot input SHA-256 mismatch"):
        verify_spatial_feature_audit(audit, spots, output)
    _write_valid_spots_and_states(spots, output)
    output.write_text("replaced\n")
    with pytest.raises(ValueError, match="state output SHA-256 mismatch"):
        verify_spatial_feature_audit(audit, spots, output)
    _write_valid_spots_and_states(spots, output)
    payload = json.loads(audit.read_text())
    payload["labels_used"] = True
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="non-label-free"):
        verify_spatial_feature_audit(audit, spots, output)


def test_spatial_feature_audit_replay_detects_setting_tampering(tmp_path):
    spots, output = tmp_path / "spots.csv", tmp_path / "states.csv"
    _write_valid_spots_and_states(spots, output)
    audit = write_spatial_feature_audit(
        spots, output, "area", "sampling_day", ("0",)
    )
    payload = json.loads(audit.read_text())
    payload["settings"]["weight_column"] = None
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="derivation replay mismatch"):
        verify_spatial_feature_audit(audit, spots, output)
