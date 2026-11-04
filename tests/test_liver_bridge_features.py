import json

import pandas as pd
import pytest

from tdln_pd1.liver_bridge_features import (FeatureBridgeMapping,
                                             patient_liver_bridge_from_features,
                                             augment_with_liver_tdln_counterfactuals,
                                             liver_bridge_parameter_draws,
                                             load_liver_bridge_parameters,
                                             validate_perfusion_manifest,
                                             validate_treatment_manifest,
                                             simulate_patient_liver_tdln,
                                             simulate_patient_liver_tdln_counterfactual,
                                             verify_liver_counterfactual_augmentation,
                                             verify_patient_liver_tdln_audit)
from tdln_pd1.simulator import Dose
from tdln_pd1.provenance import sha256_file
from tdln_pd1.liver_tdln_coupling import LiverTDLNCouplingParameters
from tdln_pd1.liver_tdln_spatial import LiverTDLNSpatialBridge
from tdln_pd1.liver_tumor import LiverTumourParameters
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


def test_patient_audit_replays_dual_inflow_settings(tmp_path):
    features = tmp_path / "features.csv"
    output = tmp_path / "output"
    output.mkdir()
    features.write_text("patient_id,sig_cytotoxic\np0,0.5\n")
    systemic, spatial, bridge_audit = simulate_patient_liver_tdln(
        pd.Series({"patient_id": "p0", "sig_cytotoxic": 0.5}), days=2, step=0.5,
        portal_flow_multiplier=0.7, arterial_flow_multiplier=0.4)
    systemic_path = output / "patient_liver_tdln_systemic.csv"
    spatial_path = output / "patient_liver_tdln_spatial.csv"
    systemic.to_csv(systemic_path, index=False)
    spatial.to_csv(spatial_path, index=False)
    audit = {**bridge_audit, "patient_id": "p0",
             "input": {"sha256": sha256_file(features)},
             "settings": {"days": 2, "step": 0.5, "dose_days": [], "dose": 1.0,
                          "portal_flow_multiplier": 0.7,
                          "arterial_flow_multiplier": 0.4,
                          "paired_counterfactual": False},
             "outputs": {systemic_path.name: sha256_file(systemic_path),
                         spatial_path.name: sha256_file(spatial_path)}}
    audit_path = output / "patient_liver_tdln_audit.json"
    audit_path.write_text(json.dumps(audit))
    assert verify_patient_liver_tdln_audit(output, features) == audit
    audit["settings"]["portal_flow_multiplier"] = 0.8
    audit_path.write_text(json.dumps(audit))
    with pytest.raises(ValueError, match="derivation replay mismatch"):
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


def test_patient_counterfactual_preserves_same_dual_inflow_in_both_arms():
    systemic, _, _, _ = simulate_patient_liver_tdln_counterfactual(
        pd.Series({"sig_cytotoxic": 0.4}), doses=(Dose(0, 1.0),), days=4,
        portal_flow_multiplier=0.6, arterial_flow_multiplier=0.3)
    control = systemic[systemic.arm.eq("control")].reset_index(drop=True)
    treated = systemic[systemic.arm.eq("anti_pd1")].reset_index(drop=True)
    assert control.liver_portal_flow.iloc[0] == pytest.approx(
        treated.liver_portal_flow.iloc[0])
    assert control.liver_arterial_flow.iloc[0] == pytest.approx(
        treated.liver_arterial_flow.iloc[0])


def test_patient_counterfactual_requires_active_treatment():
    with pytest.raises(ValueError, match="positive dose schedule"):
        simulate_patient_liver_tdln_counterfactual(
            pd.Series(dtype=float), doses=(), days=2)
    with pytest.raises(ValueError, match="positive dose schedule"):
        simulate_patient_liver_tdln_counterfactual(
            pd.Series(dtype=float), doses=(Dose(0, 0.0),), days=2)


def test_cohort_counterfactual_augmentation_is_aligned_and_outcome_free():
    features = pd.DataFrame({
        "patient_id": ["p1", "p2"],
        "sig_cytotoxic": [0.1, 1.0], "sig_exhaustion": [0.8, 0.1],
        "frac_cd8_t": [0.1, 0.5],
    })
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1.0),), days=5, step=0.5)
    assert out.patient_id.tolist() == ["p1", "p2"]
    assert "response" not in out
    assert out.liver_cf_peak_pd1_occupancy.gt(0).all()
    assert out.liver_cf_effector_tumour_auc_fold_change.gt(1).all()
    assert out.liver_cf_tumour_relative_reduction.ge(0).all()
    with pytest.raises(ValueError, match="outcome-free"):
        augment_with_liver_tdln_counterfactuals(
            features.assign(response=[0, 1]), doses=(Dose(0, 1),), days=2)


def test_cohort_counterfactual_audit_replays_and_rejects_tampering(tmp_path):
    features_path = tmp_path / "features.csv"
    output_path = tmp_path / "augmented.csv"
    audit_path = tmp_path / "augmented.csv.liver-counterfactual-audit.json"
    features = pd.DataFrame({"patient_id": ["p1"], "sig_cytotoxic": [0.5]})
    features.to_csv(features_path, index=False)
    mapping = FeatureBridgeMapping()
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1.0),), days=2, step=0.5, mapping=mapping)
    out.to_csv(output_path, index=False)
    audit = {
        "schema_version": 1,
        "stage": "label_free_liver_tdln_counterfactual_augmentation",
        "labels_used": False,
        "input": {"sha256": sha256_file(features_path)},
        "output": {"sha256": sha256_file(output_path)},
        "settings": {"days": 2, "step": 0.5, "dose_days": [0],
                     "dose_amounts": [1.0]},
        "mapping": FeatureBridgeMapping().__dict__,
        "base_parameters": {
            "coupling": LiverTDLNCouplingParameters().__dict__,
            "spatial_bridge": LiverTDLNSpatialBridge().__dict__,
        },
    }
    audit_path.write_text(json.dumps(audit))
    assert verify_liver_counterfactual_augmentation(
        audit_path, features_path, output_path) == audit
    output_path.write_text(output_path.read_text().replace("p1", "changed"))
    with pytest.raises(ValueError, match="output hash mismatch"):
        verify_liver_counterfactual_augmentation(
            audit_path, features_path, output_path)


def test_saved_bridge_parameters_change_cohort_features_and_fail_closed(tmp_path):
    parameter_path = tmp_path / "bridge.json"
    parameter_path.write_text(json.dumps({
        "coupling": LiverTDLNCouplingParameters(immune_killing_gain=0.02).__dict__,
        "spatial_bridge": LiverTDLNSpatialBridge(antigen_drainage_gain=0.8).__dict__,
        "research_only": True,
    }))
    liver, coupling, bridge = load_liver_bridge_parameters(parameter_path)
    assert liver.arterial_buffer_fraction > 0
    features = pd.DataFrame({"patient_id": ["p1"], "sig_cytotoxic": [0.5]})
    default = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1),), days=5, step=0.5)
    calibrated = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1),), days=5, step=0.5,
        base_coupling=coupling, base_bridge=bridge)
    assert calibrated.liver_cf_tumour_relative_reduction.item() > default[
        "liver_cf_tumour_relative_reduction"].item()
    parameter_path.write_text(json.dumps({
        "coupling": {}, "spatial_bridge": {}, "research_only": False,
    }))
    with pytest.raises(ValueError, match="research_only"):
        load_liver_bridge_parameters(parameter_path)


def test_correlated_bridge_uncertainty_produces_ordered_reproducible_intervals():
    liver = LiverTumourParameters()
    coupling, bridge = LiverTDLNCouplingParameters(), LiverTDLNSpatialBridge()
    diagnostics = {
        "locally_identifiable": True,
        "fit_parameters": {"bridge.antigen_drainage_gain": bridge.antigen_drainage_gain},
        "local_log_parameter_covariance": {
            "parameter_order": ["bridge.antigen_drainage_gain"],
            "log_scale_covariance": [[0.04]],
        },
    }
    draws = liver_bridge_parameter_draws(liver, coupling, bridge, diagnostics, 8, 19)
    repeated = liver_bridge_parameter_draws(liver, coupling, bridge, diagnostics, 8, 19)
    assert draws == repeated
    assert all(draw_bridge.antigen_drainage_gain > 0
               for _, _, draw_bridge in draws)
    features = pd.DataFrame({"patient_id": ["p1"], "sig_cytotoxic": [0.5]})
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1),), days=3, step=0.5,
        base_coupling=coupling, base_bridge=bridge, parameter_draws=draws)
    for base in ("liver_cf_tumour_relative_reduction",
                 "liver_cf_efferent_effector_fold_change"):
        assert out[f"{base}_q025"].item() <= out[f"{base}_q50"].item()
        assert out[f"{base}_q50"].item() <= out[f"{base}_q975"].item()
    assert out.qc_liver_cf_uncertainty_samples.item() == 8
    bad = {**diagnostics, "fit_parameters": {"bridge.antigen_drainage_gain": 0.8}}
    with pytest.raises(ValueError, match="different centers"):
        liver_bridge_parameter_draws(liver, coupling, bridge, bad, 4, 1)


def test_patient_treatment_manifest_controls_schedules_and_fails_mismatch():
    manifest = pd.DataFrame({
        "patient_id": ["p1", "p2"], "regimen_id": ["d0", "d0_d2"],
        "agent_class": ["anti_pd1", "anti_pd1"], "dose_days": ["0", "0,2"],
        "model_dose_amount": [1.0, 0.5],
        "dose_scale": ["arbitrary_model_exposure"] * 2,
        "response_definition": ["RECIST_1_1_best_overall_response"] * 2,
    })
    schedules, contract = validate_treatment_manifest(manifest, ["p1", "p2"], days=5)
    assert [dose.day for dose in schedules["p2"]] == [0, 2]
    assert [dose.amount for dose in schedules["p2"]] == [0.5, 0.5]
    assert contract["response_definition"] == "RECIST_1_1_best_overall_response"
    features = pd.DataFrame({"patient_id": ["p1", "p2"], "sig_cytotoxic": [0.2, 0.2]})
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(), patient_doses=schedules, days=5, step=0.5)
    assert out.loc[1, "liver_cf_peak_pd1_occupancy"] != out.loc[
        0, "liver_cf_peak_pd1_occupancy"]
    with pytest.raises(ValueError, match="exactly match"):
        validate_treatment_manifest(manifest.iloc[:1], ["p1", "p2"], days=5)
    with pytest.raises(ValueError, match="only agent_class"):
        validate_treatment_manifest(
            manifest.assign(agent_class=["anti_pd1", "anti_ctla4"]),
            ["p1", "p2"], days=5)
    with pytest.raises(ValueError, match="shared response_definition"):
        validate_treatment_manifest(
            manifest.assign(response_definition=["RECIST", "irRECIST"]),
            ["p1", "p2"], days=5)


def test_patient_perfusion_manifest_controls_cohort_dual_inflow_and_fails_closed():
    manifest = pd.DataFrame({
        "patient_id": ["p1", "p2"],
        "portal_flow_multiplier": [1.0, 0.55],
        "arterial_flow_multiplier": [1.0, 0.35],
        "measurement_definition": ["synthetic_duplex_ratio"] * 2,
        "normalization_reference": ["synthetic_population_baseline"] * 2,
    })
    flows, contract = validate_perfusion_manifest(manifest, ["p1", "p2"])
    assert flows["p2"] == (0.55, 0.35)
    assert contract["scale"] == "dimensionless_multiplier"
    features = pd.DataFrame({"patient_id": ["p1", "p2"], "sig_cytotoxic": [0.2, 0.2]})
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1.0),), patient_flows=flows, days=5, step=0.5)
    assert out.loc[1, "liver_cf_control_final_tumour_fraction"] != pytest.approx(
        out.loc[0, "liver_cf_control_final_tumour_fraction"])
    with pytest.raises(ValueError, match="exactly match"):
        validate_perfusion_manifest(manifest.iloc[:1], ["p1", "p2"])
    with pytest.raises(ValueError, match="strictly positive"):
        validate_perfusion_manifest(
            manifest.assign(portal_flow_multiplier=[1.0, 0.0]), ["p1", "p2"])
    with pytest.raises(ValueError, match="shared measurement_definition"):
        validate_perfusion_manifest(
            manifest.assign(measurement_definition=["duplex", "dce_mri"]),
            ["p1", "p2"])


def test_cohort_audit_binds_and_replays_perfusion_manifest(tmp_path):
    features_path = tmp_path / "features.csv"
    output_path = tmp_path / "augmented.csv"
    perfusion_path = tmp_path / "perfusion.csv"
    audit_path = tmp_path / "augmented.csv.liver-counterfactual-audit.json"
    features = pd.DataFrame({"patient_id": ["p1"], "sig_cytotoxic": [0.5]})
    perfusion = pd.DataFrame({
        "patient_id": ["p1"], "portal_flow_multiplier": [0.7],
        "arterial_flow_multiplier": [0.4],
        "measurement_definition": ["synthetic_duplex_ratio"],
        "normalization_reference": ["synthetic_baseline"],
    })
    features.to_csv(features_path, index=False)
    perfusion.to_csv(perfusion_path, index=False)
    flows, contract = validate_perfusion_manifest(perfusion, ["p1"])
    mapping = FeatureBridgeMapping()
    out = augment_with_liver_tdln_counterfactuals(
        features, doses=(Dose(0, 1.0),), days=2, step=0.5,
        mapping=mapping, patient_flows=flows)
    out.to_csv(output_path, index=False)
    audit = {
        "schema_version": 1,
        "stage": "label_free_liver_tdln_counterfactual_augmentation",
        "labels_used": False,
        "input": {"sha256": sha256_file(features_path)},
        "output": {"sha256": sha256_file(output_path)},
        "settings": {"days": 2, "step": 0.5, "dose_days": [0],
                     "dose_amounts": [1.0]},
        "mapping": mapping.__dict__,
        "base_parameters": {
            "coupling": LiverTDLNCouplingParameters().__dict__,
            "spatial_bridge": LiverTDLNSpatialBridge().__dict__,
        },
        "perfusion_manifest": {"sha256": sha256_file(perfusion_path)},
        "perfusion_contract": contract,
    }
    audit_path.write_text(json.dumps(audit))
    assert verify_liver_counterfactual_augmentation(
        audit_path, features_path, output_path,
        perfusion_manifest_path=perfusion_path) == audit
    perfusion.loc[0, "portal_flow_multiplier"] = 0.8
    perfusion.to_csv(perfusion_path, index=False)
    with pytest.raises(ValueError, match="perfusion manifest hash mismatch"):
        verify_liver_counterfactual_augmentation(
            audit_path, features_path, output_path,
            perfusion_manifest_path=perfusion_path)
