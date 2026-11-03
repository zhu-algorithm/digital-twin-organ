import json

import pandas as pd
import pytest

from tdln_pd1.liver_bridge_features import (patient_liver_bridge_from_features,
                                             simulate_patient_liver_tdln,
                                             simulate_patient_liver_tdln_counterfactual,
                                             verify_patient_liver_tdln_audit)
from tdln_pd1.simulator import Dose
from tdln_pd1.provenance import sha256_file
from tdln_pd1.liver_tdln_coupling import LiverTDLNCouplingParameters
from tdln_pd1.liver_tdln_spatial import LiverTDLNSpatialBridge
from tdln_pd1.spatial_dynamics import SpatialParameters


def test_single_cell_states_change_only_linked_bridge_mechanisms():
    base_coupling = LiverTDLNCouplingParameters()
    base_bridge = LiverTDLNSpatialBridge()
    base_spatial = SpatialParameters()
    low = patient_liver_bridge_from_features(pd.Series({
        "sig_dendritic_activation_dc": 0.1, "sig_cytotoxic_cd8_t": 0.1,
        "sig_exhaustion_cd8_t": 1.0, "frac_cd8_t": 0.1,
    }))
    high = patient_liver_bridge_from_features(pd.Series({
        "sig_dendritic_activation_dc": 1.0, "sig_cytotoxic_cd8_t": 1.0,
        "sig_exhaustion_cd8_t": 0.1, "frac_cd8_t": 0.5,
    }))
    assert high.spatial_bridge.antigen_drainage_gain > low.spatial_bridge.antigen_drainage_gain
    assert high.spatial_parameters.dc_yield > low.spatial_parameters.dc_yield
    assert high.coupling.immune_killing_gain > low.coupling.immune_killing_gain
    assert high.coupling.necrotic_antigen_gain == base_coupling.necrotic_antigen_gain
    assert high.spatial_parameters.antigen_decay == base_spatial.antigen_decay
    assert base_bridge.antigen_drainage_gain > 0


def test_missing_features_preserve_population_bases_and_are_audited():
    patient = patient_liver_bridge_from_features(pd.Series(dtype=float))
    assert patient.coupling == LiverTDLNCouplingParameters()
    assert patient.spatial_bridge == LiverTDLNSpatialBridge()
    assert patient.spatial_parameters == SpatialParameters()
    assert patient.audit["labels_used"] is False
    assert len(patient.audit["fallback_features"]) == 5
    assert set(patient.audit["multipliers"].values()) == {1.0}


def test_patient_simulation_uses_feature_specific_bridge_and_states():
    low = pd.Series({"sig_dendritic_activation": 0.0, "sig_cytotoxic": 0.0,
                     "sig_exhaustion": 1.0, "frac_cd8_t": 0.1})
    high = pd.Series({"sig_dendritic_activation": 1.0, "sig_cytotoxic": 1.0,
                      "sig_exhaustion": 0.0, "frac_cd8_t": 0.5})
    _, low_spatial, low_audit = simulate_patient_liver_tdln(low, days=5, step=0.5)
    _, high_spatial, high_audit = simulate_patient_liver_tdln(high, days=5, step=0.5)
    selector = lambda frame: frame[frame.species.eq("mature_dc")
                                   & frame.zone.eq("paracortex")].amount.iloc[-1]
    assert selector(high_spatial) > selector(low_spatial)
    assert high_audit["multipliers"]["immune_killing"] > low_audit["multipliers"]["immune_killing"]


def test_feature_bridge_rejects_outcomes_and_invalid_measurements():
    with pytest.raises(ValueError, match="outcome-free"):
        patient_liver_bridge_from_features(pd.Series({"response": 1, "sig_cytotoxic": 1.0}))
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        patient_liver_bridge_from_features(pd.Series({"frac_cd8_t": 1.2}))
    with pytest.raises(ValueError, match="finite"):
        patient_liver_bridge_from_features(pd.Series({"sig_cytotoxic": float("inf")}))
    with pytest.raises(ValueError, match="non-negative"):
        patient_liver_bridge_from_features(pd.Series({"tumour_antigen_proxy": -0.1}))


def test_patient_audit_verifier_detects_tampering(tmp_path):
    features = tmp_path / "features.csv"
    output = tmp_path / "output"
    output.mkdir()
    features.write_text("patient_id,sig_cytotoxic\np0,0.5\n")
    trajectory = output / "patient_liver_tdln_systemic.csv"
    trajectory.write_text("day,tumour\n0,0.1\n")
    audit = {
        "schema_version": 1, "labels_used": False,
        "input": {"sha256": sha256_file(features)},
        "outputs": {trajectory.name: sha256_file(trajectory)},
    }
    (output / "patient_liver_tdln_audit.json").write_text(json.dumps(audit))
    assert verify_patient_liver_tdln_audit(output, features) == audit
    trajectory.write_text("day,tumour\n0,0.9\n")
    with pytest.raises(ValueError, match="output hash mismatch"):
        verify_patient_liver_tdln_audit(output, features)


def test_patient_counterfactual_is_matched_and_reports_mechanistic_effect():
    row = pd.Series({"sig_dendritic_activation": 0.5, "sig_cytotoxic": 0.5,
                     "sig_exhaustion": 0.2, "frac_cd8_t": 0.3})
    systemic, spatial, effect, audit = simulate_patient_liver_tdln_counterfactual(
        row, doses=(Dose(0, 1.0), Dose(7, 1.0)), days=14, step=0.5)
    control = systemic[systemic.arm.eq("control")].reset_index(drop=True)
    treated = systemic[systemic.arm.eq("anti_pd1")].reset_index(drop=True)
    assert control.iloc[0].liver_tumour_fraction == pytest.approx(
        treated.iloc[0].liver_tumour_fraction)
    assert control.pd1_occupancy.max() == pytest.approx(0.0)
    assert treated.pd1_occupancy.max() > 0
    assert effect["anti_pd1_final_liver_tumour_fraction"] < effect[
        "control_final_liver_tumour_fraction"]
    assert effect["effector_tumour_auc_fold_change"] > 1
    assert effect["efferent_effector_endpoint_fold_change"] > 1
    assert effect["interpretation"].endswith("not_clinical_response")
    assert set(spatial.arm) == {"control", "anti_pd1"}
    assert audit["labels_used"] is False


def test_patient_counterfactual_requires_active_treatment():
    with pytest.raises(ValueError, match="positive dose schedule"):
        simulate_patient_liver_tdln_counterfactual(
            pd.Series(dtype=float), doses=(), days=2)
    with pytest.raises(ValueError, match="positive dose schedule"):
        simulate_patient_liver_tdln_counterfactual(
            pd.Series(dtype=float), doses=(Dose(0, 0.0),), days=2)
