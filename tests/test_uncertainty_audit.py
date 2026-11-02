import hashlib

import pytest

from tdln_pd1.uncertainty_audit import (
    verify_uncertainty_audit, write_uncertainty_audit,
)


def test_uncertainty_audit_binds_inputs_settings_and_outputs(tmp_path):
    parameters = tmp_path / "parameters.json"
    parameters.write_text("{}\n")
    diagnostics = tmp_path / "diagnostics.json"
    diagnostics.write_text("{}\n")
    (tmp_path / "uncertainty_samples.csv").write_text("draw,x\n0,1\n")
    (tmp_path / "uncertainty_summary.json").write_text("{}\n")
    inputs = {
        "parameters": {"file": parameters.name,
                       "sha256": hashlib.sha256(parameters.read_bytes()).hexdigest()},
        "diagnostics": {"file": diagnostics.name,
                        "sha256": hashlib.sha256(diagnostics.read_bytes()).hexdigest()},
    }
    audit = write_uncertainty_audit(tmp_path, inputs, {"samples": 1})
    verified = verify_uncertainty_audit(audit, {
        "parameters": parameters, "diagnostics": diagnostics,
    })
    assert verified["valid"] is True
    assert verified["verified_outputs"] == [
        "uncertainty_samples.csv", "uncertainty_summary.json",
    ]


def test_uncertainty_audit_detects_input_and_output_replacement(tmp_path):
    parameters = tmp_path / "parameters.json"
    parameters.write_text("{}\n")
    diagnostics = tmp_path / "diagnostics.json"
    diagnostics.write_text("{}\n")
    samples = tmp_path / "uncertainty_samples.csv"
    samples.write_text("draw,x\n0,1\n")
    summary = tmp_path / "uncertainty_summary.json"
    summary.write_text("{}\n")
    audit = write_uncertainty_audit(tmp_path, {
        "parameters": {"file": parameters.name,
                       "sha256": hashlib.sha256(parameters.read_bytes()).hexdigest()},
        "diagnostics": {"file": diagnostics.name,
                        "sha256": hashlib.sha256(diagnostics.read_bytes()).hexdigest()},
    }, {})
    parameters.write_text("replaced\n")
    with pytest.raises(ValueError, match="input SHA-256 mismatch"):
        verify_uncertainty_audit(audit, {
            "parameters": parameters, "diagnostics": diagnostics,
        })
    parameters.write_text("{}\n")
    summary.write_text("replaced\n")
    with pytest.raises(ValueError, match="output SHA-256 mismatch"):
        verify_uncertainty_audit(audit, {
            "parameters": parameters, "diagnostics": diagnostics,
        })
