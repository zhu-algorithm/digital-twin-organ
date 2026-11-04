import json

import numpy as np
import pandas as pd
import pytest

from tdln_pd1.liver_bridge_features import (FeatureBridgeMapping,
                                             patient_liver_bridge_from_features)
from tdln_pd1.liver_mapping_assimilation import (
    calibrate_feature_bridge_mapping, feature_bridge_mapping_draws,
    load_feature_bridge_mapping,
    save_feature_bridge_mapping_calibration, validate_feature_bridge_observations)


def _observations(coefficient=1.4):
    truth = FeatureBridgeMapping(killing_cytotoxic_coefficient=coefficient)
    rows = []
    for index, score in enumerate((0.1, 0.2, 0.35, 0.5, 0.7, 0.9)):
        feature = pd.Series({"patient_id": f"p{index}", "sig_cytotoxic_cd8_t": score})
        value = patient_liver_bridge_from_features(
            feature, mapping=truth).audit["multipliers"]["immune_killing"]
        rows.append({"patient_id": f"p{index}",
                     "target": "immune_killing_multiplier", "value": value,
                     "sd": value * 0.01, "sig_cytotoxic_cd8_t": score,
                     "source": "synthetic_unit_test"})
    return pd.DataFrame(rows)


def test_feature_bridge_calibration_recovers_synthetic_coefficient(tmp_path):
    mapping, diagnostics, residuals = calibrate_feature_bridge_mapping(
        _observations(), fit_parameters=("killing_cytotoxic_coefficient",),
        base_mapping=FeatureBridgeMapping(killing_cytotoxic_coefficient=0.4))
    assert mapping.killing_cytotoxic_coefficient == pytest.approx(1.4, rel=2e-3)
    assert diagnostics["success"] is True
    assert diagnostics["labels_used"] is False
    assert diagnostics["weighted_rmse"] < 1e-3
    save_feature_bridge_mapping_calibration(tmp_path, mapping, diagnostics, residuals)
    loaded = load_feature_bridge_mapping(
        tmp_path / "fitted_liver_feature_bridge_mapping.json")
    assert loaded == mapping
    payload = json.loads((tmp_path / "fitted_liver_feature_bridge_mapping.json").read_text())
    assert payload["research_only"] is True and payload["labels_used"] is False


def test_feature_bridge_observation_contract_fails_closed(tmp_path):
    observations = _observations()
    with pytest.raises(ValueError, match="unknown"):
        validate_feature_bridge_observations(observations.assign(target="response"))
    with pytest.raises(ValueError, match="unique"):
        validate_feature_bridge_observations(
            pd.concat([observations, observations.iloc[[0]]]))
    with pytest.raises(ValueError, match="strictly positive"):
        validate_feature_bridge_observations(observations.assign(sd=0))
    with pytest.raises(ValueError, match="outcome-free"):
        validate_feature_bridge_observations(observations.assign(response_label=1))
    with pytest.raises(ValueError, match="fit parameters"):
        calibrate_feature_bridge_mapping(
            observations, fit_parameters=("reference_cd8_fraction",))
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"mapping": {}, "research_only": True,
                               "labels_used": True}))
    with pytest.raises(ValueError, match="label-free"):
        load_feature_bridge_mapping(bad)


def test_feature_bridge_mapping_uncertainty_is_correlated_and_reproducible():
    mapping, diagnostics, _ = calibrate_feature_bridge_mapping(
        _observations(), fit_parameters=("killing_cytotoxic_coefficient",),
        base_mapping=FeatureBridgeMapping(killing_cytotoxic_coefficient=0.4))
    draws = feature_bridge_mapping_draws(mapping, diagnostics, 12, 7)
    assert draws == feature_bridge_mapping_draws(mapping, diagnostics, 12, 7)
    assert len(draws) == 12
    assert all(draw.killing_cytotoxic_coefficient > 0 for draw in draws)
    altered = {**diagnostics, "fit_parameters": {
        "killing_cytotoxic_coefficient": 0.5}}
    with pytest.raises(ValueError, match="different fit centers"):
        feature_bridge_mapping_draws(mapping, altered, 4, 7)
