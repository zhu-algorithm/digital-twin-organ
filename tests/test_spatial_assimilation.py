import numpy as np
import pandas as pd
import pytest

from tdln_pd1.spatial_assimilation import (
    calibrate_spatial, spatial_calibration_ensemble, validate_spatial_observations,
)
from tdln_pd1.spatial_dynamics import SpatialParameters, simulate_spatial


def test_spatial_calibration_recovers_synthetic_antigen_decay():
    truth = SpatialParameters(antigen_decay=0.2)
    trajectory = simulate_spatial(days=5, step=0.1, parameters=truth)
    series = trajectory[
        trajectory.species.eq("antigen") & trajectory.zone.eq("afferent_lymph")
    ]
    days = np.arange(0, 6, dtype=float)
    values = np.interp(days, series.day, series.amount)
    observations = pd.DataFrame({
        "day": days, "species": "antigen", "zone": "afferent_lymph",
        "amount": values, "sd": np.maximum(values * 0.01, 1e-6),
    })
    fitted, diagnostics, compared = calibrate_spatial(
        observations, fit_parameters=("antigen_decay",), max_nfev=80,
        profile_points=5, profile_fold_range=2.0,
    )
    assert diagnostics["success"]
    assert fitted.antigen_decay == pytest.approx(0.2, rel=0.03)
    assert diagnostics["weighted_rmse"] < 0.1
    assert diagnostics["jacobian_rank"] == 1
    assert diagnostics["locally_identifiable"] is True
    assert diagnostics["local_parameter_uncertainty"]["antigen_decay"]
    profile = diagnostics["profile_likelihood"]["profiles"]["antigen_decay"]
    assert len(profile["points"]) == 5
    assert profile["points"][2]["delta_weighted_sse"] == pytest.approx(0.0, abs=1e-9)
    assert len(compared) == len(observations)


def test_spatial_observation_contract_rejects_unknown_or_duplicate_rows():
    bad = pd.DataFrame({"day": [0], "species": ["magic"], "zone": ["paracortex"],
                        "amount": [1.0]})
    with pytest.raises(ValueError, match="unknown spatial"):
        validate_spatial_observations(bad)
    duplicate = pd.DataFrame({"day": [0, 0], "species": ["tpex", "tpex"],
                              "zone": ["paracortex", "paracortex"],
                              "amount": [1.0, 1.0]})
    with pytest.raises(ValueError, match="unique"):
        validate_spatial_observations(duplicate)
    with pytest.raises(ValueError, match="must not be empty"):
        validate_spatial_observations(
            pd.DataFrame(columns=["day", "species", "zone", "amount"])
        )
    missing_coordinate = pd.DataFrame({
        "day": [0], "species": [None], "zone": ["paracortex"], "amount": [1.0],
    })
    with pytest.raises(ValueError, match="non-missing"):
        validate_spatial_observations(missing_coordinate)


def test_spatial_log_space_fit_requires_unique_positive_base_parameters():
    observations = pd.DataFrame({
        "day": [0, 1], "species": ["antigen"] * 2,
        "zone": ["afferent_lymph"] * 2, "amount": [1.0, 0.5], "sd": [0.1, 0.1],
    })
    with pytest.raises(ValueError, match="must be unique"):
        calibrate_spatial(
            observations, fit_parameters=("antigen_decay", "antigen_decay")
        )
    with pytest.raises(ValueError, match="strictly positive"):
        calibrate_spatial(
            observations, base=SpatialParameters(antigen_decay=0.0),
            fit_parameters=("antigen_decay",),
        )


def test_spatial_calibration_can_fit_transport_rate():
    truth = SpatialParameters(dc_to_sinus=0.35)
    trajectory = simulate_spatial(days=4, step=0.1, parameters=truth)
    series = trajectory[
        trajectory.species.eq("mature_dc") & trajectory.zone.eq("afferent_lymph")
    ]
    days = np.arange(0, 5, dtype=float)
    values = np.interp(days, series.day, series.amount)
    observations = pd.DataFrame({
        "day": days, "species": "mature_dc", "zone": "afferent_lymph",
        "amount": values, "sd": np.maximum(values * 0.01, 1e-6),
    })
    fitted, diagnostics, _ = calibrate_spatial(
        observations, fit_parameters=("dc_to_sinus",), max_nfev=80,
    )
    assert diagnostics["success"]
    assert fitted.dc_to_sinus == pytest.approx(0.35, rel=0.03)


def test_spatial_calibration_ensemble_propagates_correlated_paired_effects():
    parameters = SpatialParameters(tpex_growth=0.1, checkpoint_growth_gain=4.0)
    diagnostics = {
        "locally_identifiable": True,
        "fit_parameters": {"tpex_growth": 0.1, "checkpoint_growth_gain": 4.0},
        "local_log_parameter_covariance": {
            "parameter_order": ["tpex_growth", "checkpoint_growth_gain"],
            "log_scale_covariance": [[0.01, 0.008], [0.008, 0.01]],
        },
    }
    a = spatial_calibration_ensemble(parameters, diagnostics, days=3, n=12, seed=4)
    b = spatial_calibration_ensemble(parameters, diagnostics, days=3, n=12, seed=4)
    pd.testing.assert_frame_equal(a, b)
    assert a.parameter_tpex_growth.corr(a.parameter_checkpoint_growth_gain) > 0.5
    assert a.pd1_efferent_effector_delta.gt(0).all()
    assert a.pd1_paracortex_tpex_delta.gt(0).all()
