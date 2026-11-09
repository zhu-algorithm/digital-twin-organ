import json

import pandas as pd
import pytest

from tdln_pd1.regimen_sweep import (sweep_patient_regimens,
                                    verify_regimen_sweep_audit,
                                    write_regimen_sweep_audit)


def test_regimen_sweep_is_label_free_finite_and_marks_pareto_front():
    row = pd.Series({"patient_id": "p0", "sig_cytotoxic": 0.5,
                     "sig_exhaustion": 0.2, "frac_cd8_t": 0.3})
    result = sweep_patient_regimens(
        row, dose_amounts=(0.5, 1.0), intervals=(2, 3), days=4, step=1)
    assert len(result) == 4
    assert result.pareto_front.any()
    assert result.interpretation.nunique() == 1
    assert (result.dose_count > 0).all()
    assert result.drop(columns=["dose_days", "interpretation", "pareto_front"]).notna().all().all()
    with pytest.raises(ValueError, match="outcome-free"):
        sweep_patient_regimens(pd.Series({"response": 1}), dose_amounts=(1,),
                               intervals=(2,), days=3, step=1)


def test_regimen_sweep_validates_grid():
    row = pd.Series({"sig_cytotoxic": 0.5})
    with pytest.raises(ValueError, match="dose_amounts"):
        sweep_patient_regimens(row, dose_amounts=(0,), intervals=(2,), days=3)
    with pytest.raises(ValueError, match="intervals"):
        sweep_patient_regimens(row, dose_amounts=(1,), intervals=(2, 2), days=3)


def test_regimen_sweep_audit_replays_and_detects_tampering(tmp_path):
    features = tmp_path / "features.csv"
    output = tmp_path / "sweep.csv"
    audit = tmp_path / "sweep.audit.json"
    pd.DataFrame({"patient_id": ["p0"], "sig_cytotoxic": [0.5]}).to_csv(
        features, index=False)
    settings = {"dose_amounts": [0.5, 1.0], "intervals": [2.0],
                "days": 3.0, "step": 1.0,
                "portal_flow_multiplier": None,
                "arterial_flow_multiplier": None}
    sweep_patient_regimens(pd.read_csv(features).iloc[0], **{
        "dose_amounts": settings["dose_amounts"], "intervals": settings["intervals"],
        "days": settings["days"], "step": settings["step"],
    }).to_csv(output, index=False)
    write_regimen_sweep_audit(audit, features_path=features, output_path=output,
                              patient_id="p0", settings=settings)
    verified = verify_regimen_sweep_audit(audit, features, output)
    assert verified["labels_used"] is False
    frame = pd.read_csv(output)
    frame.loc[0, "tumour_relative_reduction"] += 0.1
    frame.to_csv(output, index=False)
    with pytest.raises(ValueError, match="output hash mismatch"):
        verify_regimen_sweep_audit(audit, features, output)


def test_regimen_sweep_audit_rejects_forged_hash_with_wrong_derivation(tmp_path):
    from tdln_pd1.provenance import sha256_file

    features = tmp_path / "features.csv"
    output = tmp_path / "sweep.csv"
    audit = tmp_path / "sweep.audit.json"
    pd.DataFrame({"patient_id": ["p0"], "sig_cytotoxic": [0.5]}).to_csv(features, index=False)
    settings = {"dose_amounts": [1.0], "intervals": [2.0], "days": 3.0,
                "step": 1.0, "portal_flow_multiplier": None,
                "arterial_flow_multiplier": None}
    sweep_patient_regimens(pd.read_csv(features).iloc[0], dose_amounts=(1,),
                           intervals=(2,), days=3, step=1).to_csv(output, index=False)
    write_regimen_sweep_audit(audit, features_path=features, output_path=output,
                              patient_id="p0", settings=settings)
    frame = pd.read_csv(output)
    frame.loc[0, "pareto_front"] = False
    frame.to_csv(output, index=False)
    payload = json.loads(audit.read_text())
    payload["output"]["sha256"] = sha256_file(output)
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="derivation replay mismatch"):
        verify_regimen_sweep_audit(audit, features, output)


def test_regimen_sweep_binds_parameter_source_and_loaded_content(tmp_path):
    from tdln_pd1.hybrid import FeatureStateMapping

    features = tmp_path / "features.csv"
    output = tmp_path / "sweep.csv"
    audit = tmp_path / "sweep.audit.json"
    mapping = tmp_path / "mapping.json"
    pd.DataFrame({"patient_id": ["p0"], "sig_cytotoxic": [0.5]}).to_csv(
        features, index=False)
    mapping.write_text(json.dumps(FeatureStateMapping().__dict__))
    settings = {"dose_amounts": [1.0], "intervals": [2.0], "days": 3.0,
                "step": 1.0, "portal_flow_multiplier": None,
                "arterial_flow_multiplier": None}
    sweep_patient_regimens(pd.read_csv(features).iloc[0], dose_amounts=(1,),
                           intervals=(2,), days=3, step=1).to_csv(output, index=False)
    write_regimen_sweep_audit(
        audit, features_path=features, output_path=output, patient_id="p0",
        settings=settings, source_paths={"feature_state_mapping": mapping})
    assert verify_regimen_sweep_audit(
        audit, features, output,
        source_paths={"feature_state_mapping": mapping})["research_only"] is True
    with pytest.raises(ValueError, match="source paths"):
        verify_regimen_sweep_audit(audit, features, output)
    mapping.write_text(json.dumps({**FeatureStateMapping().__dict__,
                                   "tcf7_tpex_coefficient": 1.1}))
    with pytest.raises(ValueError, match="source mismatch"):
        verify_regimen_sweep_audit(
            audit, features, output,
            source_paths={"feature_state_mapping": mapping})
