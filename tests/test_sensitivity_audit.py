import pytest

from tdln_pd1.sensitivity_audit import (
    verify_sensitivity_audit, write_sensitivity_audit,
)


def test_sensitivity_audit_binds_and_verifies_exact_outputs(tmp_path):
    samples = tmp_path / "samples.csv"
    ranking = tmp_path / "ranking.csv"
    samples.write_text("draw,value\n0,1\n")
    ranking.write_text("parameter,rho\nx,0.5\n")
    audit = write_sensitivity_audit(
        tmp_path, "audit.json", "spatial_sensitivity",
        inputs={"parameter_source": None}, settings={"seed": 42},
        generated_filenames=(samples.name, ranking.name),
    )
    verified = verify_sensitivity_audit(audit)
    assert verified["valid"] is True
    assert verified["stage"] == "spatial_sensitivity"
    assert verified["verified_outputs"] == ["ranking.csv", "samples.csv"]
    assert len(verified["audit_sha256"]) == 64


def test_sensitivity_audit_detects_posthoc_result_replacement(tmp_path):
    samples = tmp_path / "samples.csv"
    samples.write_text("draw,value\n0,1\n")
    audit = write_sensitivity_audit(
        tmp_path, "audit.json", "systemic_sensitivity",
        inputs={}, settings={}, generated_filenames=(samples.name,),
    )
    samples.write_text("draw,value\n0,999\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_sensitivity_audit(audit)


def test_sensitivity_audit_rejects_missing_output_before_writing(tmp_path):
    with pytest.raises(ValueError, match="missing before audit"):
        write_sensitivity_audit(
            tmp_path, "audit.json", "systemic_sensitivity",
            inputs={}, settings={}, generated_filenames=("missing.csv",),
        )
