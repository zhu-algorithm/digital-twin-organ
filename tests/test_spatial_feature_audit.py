import json

import pytest

from tdln_pd1.spatial_feature_audit import (
    verify_spatial_feature_audit, write_spatial_feature_audit,
)


def test_spatial_feature_audit_binds_spots_weight_and_state_output(tmp_path):
    spots = tmp_path / "spots.csv"
    spots.write_text("patient_id,zone,area\np1,paracortex,2\n")
    output = tmp_path / "states.csv"
    output.write_text("patient_id,zone,species,amount\np1,paracortex,tpex,2\n")
    audit = write_spatial_feature_audit(spots, output, "area")
    verified = verify_spatial_feature_audit(audit, spots, output)
    assert verified["valid"] is True
    assert verified["payload"]["settings"]["weight_column"] == "area"
    assert verified["payload"]["labels_used"] is False
    assert len(verified["audit_sha256"]) == 64


def test_spatial_feature_audit_detects_replacement_and_label_drift(tmp_path):
    spots = tmp_path / "spots.csv"
    spots.write_text("patient_id,zone\np1,paracortex\n")
    output = tmp_path / "states.csv"
    output.write_text("patient_id,zone,species,amount\np1,paracortex,tpex,1\n")
    audit = write_spatial_feature_audit(spots, output, None)
    spots.write_text("patient_id,zone\np1,medulla\n")
    with pytest.raises(ValueError, match="spot input SHA-256 mismatch"):
        verify_spatial_feature_audit(audit, spots, output)
    spots.write_text("patient_id,zone\np1,paracortex\n")
    output.write_text("replaced\n")
    with pytest.raises(ValueError, match="state output SHA-256 mismatch"):
        verify_spatial_feature_audit(audit, spots, output)
    output.write_text("patient_id,zone,species,amount\np1,paracortex,tpex,1\n")
    payload = json.loads(audit.read_text())
    payload["labels_used"] = True
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="non-label-free"):
        verify_spatial_feature_audit(audit, spots, output)
