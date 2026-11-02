import hashlib

import pytest

from tdln_pd1.calibration_audit import (
    verify_calibration_audit, write_calibration_audit,
)


def test_calibration_audit_binds_exact_outputs_and_detects_replacement(tmp_path):
    observations = tmp_path / "observations.csv"
    observations.write_text("day,value\n0,1\n")
    parameters = tmp_path / "parameters.json"
    diagnostics = tmp_path / "diagnostics.json"
    residuals = tmp_path / "residuals.csv"
    parameters.write_text("{}\n")
    diagnostics.write_text("{}\n")
    residuals.write_text("day,residual\n0,0\n")
    audit = write_calibration_audit(
        tmp_path, "audit.json", "systemic_calibration",
        inputs={"observations": {
            "file": observations.name,
            "sha256": hashlib.sha256(observations.read_bytes()).hexdigest(),
        }},
        settings={"fit_parameters": ["tumour_growth"]},
        generated_filenames=(parameters.name, diagnostics.name, residuals.name),
    )
    verified = verify_calibration_audit(audit, {"observations": observations})
    assert verified["valid"] is True
    assert verified["stage"] == "systemic_calibration"
    assert verified["verified_inputs"] == ["observations"]
    assert verified["verified_outputs"] == [
        "diagnostics.json", "parameters.json", "residuals.csv",
    ]
    parameters.write_text('{"changed": true}\n')
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_calibration_audit(audit, {"observations": observations})


def test_calibration_audit_requires_and_verifies_declared_input(tmp_path):
    observations = tmp_path / "observations.csv"
    observations.write_text("day,value\n0,1\n")
    output = tmp_path / "parameters.json"
    output.write_text("{}\n")
    audit = write_calibration_audit(
        tmp_path, "audit.json", "systemic_calibration",
        inputs={"observations": {
            "file": observations.name,
            "sha256": hashlib.sha256(observations.read_bytes()).hexdigest(),
        }}, settings={}, generated_filenames=(output.name,),
    )
    with pytest.raises(ValueError, match="requires input path"):
        verify_calibration_audit(audit)
    observations.write_text("day,value\n0,999\n")
    with pytest.raises(ValueError, match="input SHA-256 mismatch"):
        verify_calibration_audit(audit, {"observations": observations})


def test_calibration_audit_rejects_unknown_stage_and_missing_output(tmp_path):
    with pytest.raises(ValueError, match="unsupported"):
        write_calibration_audit(tmp_path, "audit.json", "unknown", {}, {}, ())
    with pytest.raises(ValueError, match="missing before audit"):
        write_calibration_audit(
            tmp_path, "audit.json", "spatial_calibration", {}, {}, ("missing.csv",),
        )
