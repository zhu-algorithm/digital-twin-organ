import pandas as pd
import pytest

from tdln_pd1.hybrid import augment_with_twin, initial_state_from_features


def test_feature_mapping_changes_biologically_linked_states():
    low = initial_state_from_features(pd.Series({"sig_tcf7_memory": 0.1, "sig_cytotoxic": 0.1}))
    high = initial_state_from_features(pd.Series({"sig_tcf7_memory": 1.5, "sig_cytotoxic": 1.5}))
    assert high.tpex_tdln > low.tpex_tdln
    assert high.effector_tumour > low.effector_tumour


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


def test_duplicate_patient_is_rejected():
    with pytest.raises(ValueError):
        augment_with_twin(pd.DataFrame({"patient_id": ["a", "a"]}), days=1)
