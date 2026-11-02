import pandas as pd
import pytest

from tdln_pd1.hybrid import FeatureStateMapping, initial_state_from_features
from tdln_pd1.mapping_assimilation import (
    calibrate_feature_state_mapping, validate_mapping_observations,
)


def _synthetic_mapping_observations() -> pd.DataFrame:
    truth = FeatureStateMapping(tcf7_coefficient=1.7)
    rows = []
    for index, value in enumerate((0.1, 0.2, 0.4, 0.6, 0.8, 1.0)):
        features = pd.Series({"sig_tcf7_memory_cd8_t": value})
        target = initial_state_from_features(features, mapping=truth).tpex_tdln
        rows.append({
            "patient_id": f"p{index}", "target": "tpex_tdln",
            "value": target, "sd": 100.0, "sig_tcf7_memory_cd8_t": value,
        })
    return pd.DataFrame(rows)


def test_mapping_calibration_recovers_synthetic_truth_and_covariance():
    fitted, diagnostics, compared = calibrate_feature_state_mapping(
        _synthetic_mapping_observations(),
        fit_parameters=("tcf7_coefficient",),
    )
    assert fitted.tcf7_coefficient == pytest.approx(1.7, rel=1e-5)
    assert diagnostics["success"] is True
    assert diagnostics["locally_identifiable"] is True
    assert diagnostics["jacobian_rank"] == 1
    assert diagnostics["local_log_parameter_covariance"]["parameter_order"] == [
        "tcf7_coefficient"
    ]
    assert compared.weighted_residual.abs().max() < 1e-4


def test_mapping_calibration_supports_spatial_targets():
    observations = _synthetic_mapping_observations().copy()
    observations["target"] = "paracortex_tpex"
    observations["value"] = observations["value"] / 80000.0 * 80000.0
    fitted, diagnostics, _ = calibrate_feature_state_mapping(
        observations, fit_parameters=("tcf7_coefficient",),
    )
    assert fitted.tcf7_coefficient == pytest.approx(1.7, rel=1e-5)
    assert diagnostics["locally_identifiable"] is True


def test_mapping_observation_contract_fails_closed():
    base = _synthetic_mapping_observations()
    with pytest.raises(ValueError, match="unknown mapping targets"):
        validate_mapping_observations(base.assign(target="response"))
    with pytest.raises(ValueError, match="unique patient-target"):
        validate_mapping_observations(pd.concat([base, base.iloc[[0]]], ignore_index=True))
    with pytest.raises(ValueError, match="sd finite/positive"):
        validate_mapping_observations(base.assign(sd=0))
    with pytest.raises(ValueError, match="invalid mapping fit parameters"):
        calibrate_feature_state_mapping(base, fit_parameters=("not_a_parameter",))
    with pytest.raises(ValueError, match="unknown mapping observation columns"):
        validate_mapping_observations(base.assign(response=1))
    repeated = pd.concat([
        base.iloc[[0]], base.iloc[[0]].assign(target="paracortex_tpex",
                                             sig_tcf7_memory_cd8_t=0.9),
    ], ignore_index=True)
    with pytest.raises(ValueError, match="consistent across targets"):
        validate_mapping_observations(repeated)
    with pytest.raises(ValueError, match="unsupported mapping calibration parameters"):
        calibrate_feature_state_mapping(
            base, fit_parameters=("signature_lower_multiplier",)
        )
