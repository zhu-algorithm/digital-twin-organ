import hashlib

import pytest

from tdln_pd1.uncertainty_audit import (
    verify_uncertainty_audit, write_uncertainty_audit,
)
from tdln_pd1.spatial_feature_audit import write_spatial_feature_audit


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


def test_spatial_uncertainty_stage_supports_spatial_output_contract(tmp_path):
    parameters = tmp_path / "spatial_parameters.json"
    parameters.write_text("{}\n")
    diagnostics = tmp_path / "spatial_diagnostics.json"
    diagnostics.write_text("{}\n")
    (tmp_path / "spatial_uncertainty_samples.csv").write_text("draw,x\n0,1\n")
    (tmp_path / "spatial_uncertainty_summary.json").write_text("{}\n")
    audit = write_uncertainty_audit(
        tmp_path,
        inputs={
            "parameters": {"file": parameters.name,
                           "sha256": hashlib.sha256(parameters.read_bytes()).hexdigest()},
            "diagnostics": {"file": diagnostics.name,
                            "sha256": hashlib.sha256(diagnostics.read_bytes()).hexdigest()},
            "optional": None,
        },
        settings={"paired_control": True},
        stage="spatial_calibration_uncertainty_propagation",
        output_filenames=(
            "spatial_uncertainty_samples.csv", "spatial_uncertainty_summary.json",
        ),
    )
    verified = verify_uncertainty_audit(audit, {
        "parameters": parameters, "diagnostics": diagnostics,
    })
    assert verified["stage"] == "spatial_calibration_uncertainty_propagation"
    assert "optional" not in verified["payload"]["inputs"]


def test_uncertainty_verifier_rechecks_bound_spatial_derivation(tmp_path):
    spots = tmp_path / "spots.csv"
    spots.write_text("patient_id,zone\np1,paracortex\n")
    states = tmp_path / "states.csv"
    states.write_text("patient_id,zone,species,amount\np1,paracortex,tpex,1\n")
    sidecar = write_spatial_feature_audit(spots, states, None)
    parameters = tmp_path / "parameters.json"
    parameters.write_text("{}\n")
    (tmp_path / "uncertainty_samples.csv").write_text("draw,x\n0,1\n")
    (tmp_path / "uncertainty_summary.json").write_text("{}\n")
    inputs = {
        "parameters": {"file": parameters.name,
                       "sha256": hashlib.sha256(parameters.read_bytes()).hexdigest()},
        "initial_zones": {
            "file": states.name,
            "sha256": hashlib.sha256(states.read_bytes()).hexdigest(),
            "derivation_audit": {
                "file": sidecar.name,
                "sha256": hashlib.sha256(sidecar.read_bytes()).hexdigest(),
                "stage": "label_free_spatial_state_featurization",
            },
        },
    }
    audit = write_uncertainty_audit(tmp_path, inputs, {})
    paths = {"parameters": parameters, "initial_zones": states}
    assert verify_uncertainty_audit(audit, paths)["valid"] is True
    sidecar.write_text("replaced\n")
    with pytest.raises(ValueError, match="derivation audit SHA-256 mismatch"):
        verify_uncertainty_audit(audit, paths)
