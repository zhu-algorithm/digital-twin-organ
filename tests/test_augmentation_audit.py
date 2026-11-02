import hashlib

import pytest

from tdln_pd1.augmentation_audit import (
    verify_augmentation_audit, write_augmentation_audit,
)


def test_augmentation_audit_binds_all_inputs_and_output(tmp_path):
    features = tmp_path / "features.csv"
    features.write_text("patient_id,x\np1,1\n")
    states = tmp_path / "states.csv"
    states.write_text("patient_id,zone\np1,paracortex\n")
    output = tmp_path / "hybrid.csv"
    output.write_text("patient_id,x,twin\np1,1,2\n")
    audit = write_augmentation_audit(
        output,
        inputs={
            "features": {"file": features.name,
                         "sha256": hashlib.sha256(features.read_bytes()).hexdigest()},
            "spatial_initial_zones": {"file": states.name,
                                      "sha256": hashlib.sha256(states.read_bytes()).hexdigest()},
        },
        settings={"days": 7, "labels_used": False},
    )
    verified = verify_augmentation_audit(
        audit, output, {"features": features, "spatial_initial_zones": states},
    )
    assert verified["valid"] is True
    assert verified["verified_inputs"] == ["features", "spatial_initial_zones"]
    assert len(verified["output_sha256"]) == 64


def test_augmentation_audit_detects_input_and_output_replacement(tmp_path):
    features = tmp_path / "features.csv"
    features.write_text("patient_id,x\np1,1\n")
    output = tmp_path / "hybrid.csv"
    output.write_text("patient_id,x\np1,1\n")
    audit = write_augmentation_audit(
        output,
        inputs={"features": {"file": features.name,
                             "sha256": hashlib.sha256(features.read_bytes()).hexdigest()}},
        settings={},
    )
    features.write_text("patient_id,x\np1,999\n")
    with pytest.raises(ValueError, match="input SHA-256 mismatch"):
        verify_augmentation_audit(audit, output, {"features": features})
    features.write_text("patient_id,x\np1,1\n")
    output.write_text("patient_id,x\np1,999\n")
    with pytest.raises(ValueError, match="output SHA-256 mismatch"):
        verify_augmentation_audit(audit, output, {"features": features})
