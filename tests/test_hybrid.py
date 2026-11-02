import pandas as pd
import pytest

from tdln_pd1.hybrid import (
    augment_with_twin, initial_state_from_features, spatial_initial_state_from_features,
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
    assert out.qc_systemic_parameters_supplied.eq(1).all()
    assert out.qc_spatial_parameters_supplied.eq(1).all()
    assert "response" not in out
    assert out.loc[out.patient_id.eq("b"), "twin_spatial_treated_paracortex_tpex"].item() > (
        out.loc[out.patient_id.eq("a"), "twin_spatial_treated_paracortex_tpex"].item()
    )


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
