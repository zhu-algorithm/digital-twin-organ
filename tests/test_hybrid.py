import pandas as pd
import pytest

from tdln_pd1.hybrid import (
    FeatureStateMapping, augment_with_twin, initial_state_from_features,
    spatial_initial_state_from_features,
)
from tdln_pd1.simulator import TwinParameters
from tdln_pd1.spatial_dynamics import SPECIES, ZONE_NAMES, SpatialParameters


def _complete_spatial_state(patient_id: str, tpex: float, effector: float) -> pd.DataFrame:
    rows = []
    for species in SPECIES:
        for zone in ZONE_NAMES:
            amount = 0.0
            if species == "tpex" and zone == "paracortex":
                amount = tpex
            elif species == "effector" and zone == "paracortex":
                amount = effector
            rows.append({"patient_id": patient_id, "species": species, "zone": zone,
                         "amount": amount, "zone_observed": zone == "paracortex"})
    return pd.DataFrame(rows)


def test_feature_mapping_changes_biologically_linked_states():
    low = initial_state_from_features(pd.Series({"sig_tcf7_memory": 0.1, "sig_cytotoxic": 0.1}))
    high = initial_state_from_features(pd.Series({"sig_tcf7_memory": 1.5, "sig_cytotoxic": 1.5}))
    assert high.tpex_tdln > low.tpex_tdln
    assert high.effector_tumour > low.effector_tumour


def test_feature_mapping_prefers_lineage_state_and_falls_back_to_global():
    lineage_low = initial_state_from_features(pd.Series({
        "sig_cytotoxic_cd8_t": 0.1, "sig_cytotoxic": 2.0,
        "sig_tcf7_memory_cd8_t": 0.1, "sig_tcf7_memory": 2.0,
    }))
    lineage_high = initial_state_from_features(pd.Series({
        "sig_cytotoxic_cd8_t": 1.0, "sig_cytotoxic": 0.0,
        "sig_tcf7_memory_cd8_t": 1.0, "sig_tcf7_memory": 0.0,
    }))
    fallback = initial_state_from_features(pd.Series({
        "sig_cytotoxic_cd8_t": None, "sig_cytotoxic": 1.0,
        "sig_tcf7_memory_cd8_t": None, "sig_tcf7_memory": 1.0,
    }))
    assert lineage_high.effector_tumour > lineage_low.effector_tumour
    assert lineage_high.tpex_tdln > lineage_low.tpex_tdln
    assert fallback.effector_tumour == lineage_high.effector_tumour
    assert fallback.tpex_tdln == lineage_high.tpex_tdln


def test_feature_mapping_initializes_patient_specific_spatial_states():
    low = spatial_initial_state_from_features(pd.Series({
        "sig_tcf7_memory": 0.1, "sig_cytotoxic": 0.1,
        "sig_dendritic_activation": 0.1, "tumour_antigen_proxy": 0.1,
    }))
    high = spatial_initial_state_from_features(pd.Series({
        "sig_tcf7_memory": 1.5, "sig_cytotoxic": 1.5,
        "sig_dendritic_activation": 1.5, "tumour_antigen_proxy": 1.5,
    }))
    assert high.afferent_antigen > low.afferent_antigen
    assert high.afferent_dc > low.afferent_dc
    assert high.paracortex_tpex > low.paracortex_tpex
    assert high.paracortex_effector > low.paracortex_effector


def test_canonical_cd8_fraction_alias_and_missing_antigen_proxy_are_supported():
    low = spatial_initial_state_from_features(pd.Series({"frac_cd8_t_cell": 0.05}))
    high = spatial_initial_state_from_features(pd.Series({"frac_cd8_t_cell": 0.80}))
    assert high.paracortex_effector > low.paracortex_effector
    assert high.afferent_antigen == low.afferent_antigen


def test_patient_feature_mapping_rejects_supplied_invalid_measurements():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        initial_state_from_features(pd.Series({"frac_cd8_t": 1.2}))
    with pytest.raises(ValueError, match="finite"):
        spatial_initial_state_from_features(pd.Series({"tumour_antigen_proxy": float("inf")}))
    with pytest.raises(ValueError, match="numeric"):
        spatial_initial_state_from_features(pd.Series({"sig_tcf7_memory": "high"}))


def test_signed_finite_signatures_map_monotonically_to_bounded_positive_states():
    negative = initial_state_from_features(pd.Series({
        "sig_cytotoxic": -2.0, "sig_tcf7_memory": -2.0,
    }))
    centered = initial_state_from_features(pd.Series({
        "sig_cytotoxic": 0.0, "sig_tcf7_memory": 0.0,
    }))
    positive = initial_state_from_features(pd.Series({
        "sig_cytotoxic": 2.0, "sig_tcf7_memory": 2.0,
    }))
    assert 0 < negative.effector_tumour < centered.effector_tumour < positive.effector_tumour
    assert 0 < negative.tpex_tdln < centered.tpex_tdln < positive.tpex_tdln

    spatial_negative = spatial_initial_state_from_features(pd.Series({
        "sig_dendritic_activation": -2.0,
    }))
    assert spatial_negative.afferent_dc > 0
    with pytest.raises(ValueError, match="non-negative"):
        spatial_initial_state_from_features(pd.Series({"tumour_antigen_proxy": -0.1}))


def test_truly_missing_optional_patient_features_keep_documented_defaults():
    absent = spatial_initial_state_from_features(pd.Series(dtype=float))
    missing = spatial_initial_state_from_features(pd.Series({"sig_cytotoxic": None}))
    assert absent == missing


def test_versioned_feature_state_mapping_controls_both_model_layers():
    row = pd.Series({"sig_tcf7_memory": 0.5, "sig_cytotoxic": 0.5})
    default_systemic = initial_state_from_features(row)
    default_spatial = spatial_initial_state_from_features(row)
    stronger = FeatureStateMapping(tcf7_coefficient=2.0, cytotoxic_coefficient=2.0)
    mapped_systemic = initial_state_from_features(row, mapping=stronger)
    mapped_spatial = spatial_initial_state_from_features(row, mapping=stronger)
    assert mapped_systemic.tpex_tdln > default_systemic.tpex_tdln
    assert mapped_systemic.effector_tumour > default_systemic.effector_tumour
    assert mapped_spatial.paracortex_tpex > default_spatial.paracortex_tpex
    assert mapped_spatial.paracortex_effector > default_spatial.paracortex_effector


def test_feature_state_mapping_validation_fails_closed():
    row = pd.Series(dtype=float)
    with pytest.raises(ValueError, match="strictly positive"):
        initial_state_from_features(row, mapping=FeatureStateMapping(tcf7_coefficient=0))
    with pytest.raises(ValueError, match="below upper bounds"):
        spatial_initial_state_from_features(
            row, mapping=FeatureStateMapping(signature_lower_multiplier=6.0)
        )


def test_twin_augmentation_is_label_free_and_patient_aligned():
    features = pd.DataFrame({
        "patient_id": ["a", "b"], "sig_cytotoxic": [0.2, 1.0],
        "sig_exhaustion": [0.3, 1.2], "sig_tcf7_memory": [1.0, 0.2],
        "sig_dendritic_activation": [0.3, 0.8], "frac_treg": [0.1, 0.2],
        "frac_cd8_t": [0.3, 0.5],
    })
    out = augment_with_twin(features, days=3, dose_days=(0,), dose=1)
    assert out.patient_id.tolist() == ["a", "b"]
    assert "response" not in out
    assert out.twin_peak_pd1_occupancy.gt(0).all()
    assert out.twin_treatment_tumour_delta.lt(0).all()
    assert out.twin_spatial_efferent_effector_gain.gt(0).all()
    assert out.twin_spatial_treated_paracortex_tpex.gt(0).all()


def test_augmentation_records_explicit_mapping_without_leaking_qc_into_model_role():
    features = pd.DataFrame({"patient_id": ["a"], "sig_tcf7_memory": [0.5]})
    out = augment_with_twin(
        features, days=1, dose_days=(0,),
        feature_state_mapping=FeatureStateMapping(tcf7_coefficient=2.0)
    )
    assert out.qc_feature_state_mapping_supplied.item() == 1


def test_twin_augmentation_rejects_outcomes_in_patient_features_before_simulation():
    base = pd.DataFrame({"patient_id": ["a", "b"], "sig_cytotoxic": [0.1, 0.2]})
    for column in ("response", "RECIST_category", "pfs_months", "overall_survival_days"):
        contaminated = base.assign(**{column: [0, 1]})
        with pytest.raises(ValueError, match="outcome-free"):
            augment_with_twin(contaminated, days=1)


def test_duplicate_patient_is_rejected():
    with pytest.raises(ValueError):
        augment_with_twin(pd.DataFrame({"patient_id": ["a", "a"]}), days=1)


def test_measured_spatial_states_drive_hybrid_endpoints_without_labels():
    features = pd.DataFrame({"patient_id": ["a", "b"], "sig_tcf7_memory": [0.2, 0.2]})
    states = pd.concat([
        _complete_spatial_state("a", 1.0e4, 5.0e3),
        _complete_spatial_state("b", 2.0e5, 8.0e4),
    ], ignore_index=True)
    out = augment_with_twin(
        features, days=3, dose_days=(0,), spatial_states=states,
        systemic_parameters=TwinParameters(tumour_growth=0.15),
        spatial_parameters=SpatialParameters(effector_to_medulla=0.2),
    )
    assert out.qc_spatial_initialization_measured.eq(1).all()
    assert out.qc_spatial_zone_coverage.eq(1 / 6).all()
    assert out.qc_spatial_sampling_coverage_mean.isna().all()
    assert out.qc_systemic_parameters_supplied.eq(1).all()
    assert out.qc_spatial_parameters_supplied.eq(1).all()
    assert "response" not in out
    assert out.loc[out.patient_id.eq("b"), "twin_spatial_treated_paracortex_tpex"].item() > (
        out.loc[out.patient_id.eq("a"), "twin_spatial_treated_paracortex_tpex"].item()
    )


def test_hybrid_reports_sampling_support_without_using_it_as_state():
    features = pd.DataFrame({"patient_id": ["a"], "sig_tcf7_memory": [0.2]})
    states = _complete_spatial_state("a", 1.0e4, 5.0e3)
    coverage = {zone: (0.5 if zone == "paracortex" else 0.0)
                for zone in states.zone.unique()}
    states["sampling_unit_coverage"] = states.zone.map(coverage)
    out = augment_with_twin(features, days=1, dose_days=(0,), spatial_states=states)
    assert out.qc_spatial_sampling_coverage_mean.item() == pytest.approx(0.5 / 6)
    assert out.qc_spatial_sampling_coverage_min.item() == 0.0


def test_spatial_hybrid_fails_closed_on_patient_mismatch_or_outcome_column():
    features = pd.DataFrame({"patient_id": ["a", "b"]})
    states = pd.DataFrame({"patient_id": ["a"], "species": ["tpex"],
                           "zone": ["paracortex"], "amount": [1.0]})
    with pytest.raises(ValueError, match="patient sets differ"):
        augment_with_twin(features, days=1, spatial_states=states)
    states = pd.concat([states, states.assign(patient_id="b")], ignore_index=True)
    states["response"] = [0, 1]
    with pytest.raises(ValueError, match="outcome-free"):
        augment_with_twin(features, days=1, spatial_states=states)


def test_hybrid_exposes_reproducible_label_free_mechanistic_uncertainty():
    features = pd.DataFrame({"patient_id": ["a"], "sig_tcf7_memory": [0.5]})
    systemic = (
        TwinParameters(tumour_growth=0.12),
        TwinParameters(tumour_growth=0.18),
        TwinParameters(tumour_growth=0.24),
    )
    spatial = (
        SpatialParameters(effector_to_medulla=0.12),
        SpatialParameters(effector_to_medulla=0.18),
        SpatialParameters(effector_to_medulla=0.24),
    )
    first = augment_with_twin(
        features, days=2, dose_days=(0,),
        systemic_parameter_draws=systemic, spatial_parameter_draws=spatial,
    )
    second = augment_with_twin(
        features, days=2, dose_days=(0,),
        systemic_parameter_draws=systemic, spatial_parameter_draws=spatial,
    )
    uncertainty = [name for name in first
                   if name.startswith("twin_") and "_uncertainty_" in name]
    assert len(uncertainty) == 16
    assert first[uncertainty].equals(second[uncertainty])
    assert first.filter(like="_uncertainty_width_95").ge(0).all().all()
    assert first.qc_mechanistic_uncertainty_samples.item() == 3
    assert not any("response" in name or "outcome" in name for name in uncertainty)


def test_hybrid_uncertainty_draw_contract_fails_closed():
    features = pd.DataFrame({"patient_id": ["a"]})
    with pytest.raises(ValueError, match="non-empty"):
        augment_with_twin(features, days=1, systemic_parameter_draws=())
    with pytest.raises(ValueError, match="draw counts must match"):
        augment_with_twin(
            features, days=1,
            systemic_parameter_draws=(TwinParameters(), TwinParameters()),
            spatial_parameter_draws=(SpatialParameters(),),
        )


def test_mapping_calibration_draws_propagate_into_patient_mechanistic_intervals():
    features = pd.DataFrame({"patient_id": ["a"], "sig_tcf7_memory_cd8_t": [0.8]})
    mappings = tuple(
        FeatureStateMapping(tcf7_coefficient=value) for value in (0.6, 1.0, 1.8)
    )
    out = augment_with_twin(
        features, days=2, dose_days=(0,),
        feature_state_mapping=FeatureStateMapping(),
        feature_state_mapping_draws=mappings,
    )
    assert out.qc_mapping_uncertainty_supplied.item() == 1
    assert out.qc_mechanistic_uncertainty_samples.item() == 3
    assert out.twin_treated_tpex_tdln_uncertainty_width_95.item() > 0
    assert out.twin_spatial_efferent_effector_gain_uncertainty_width_95.item() >= 0
