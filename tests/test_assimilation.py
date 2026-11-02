import numpy as np
import pandas as pd
import pytest

from tdln_pd1.assimilation import (
    calibrate, calibration_ensemble, ensemble, validate_observations,
)
from tdln_pd1.simulator import TwinParameters, simulate


def test_calibration_recovers_synthetic_growth_rate():
    truth = TwinParameters(tumour_growth=0.11)
    trajectory = simulate(days=14, parameters=truth)
    days = np.array([0, 2, 4, 7, 10, 14], dtype=float)
    values = np.interp(days, trajectory.day, trajectory.tumour_cells)
    observations = pd.DataFrame({"day": days, "variable": "tumour_cells", "value": values, "sd": values * 0.01})
    fitted, diagnostics, compared = calibrate(observations, fit_parameters=("tumour_growth",), max_nfev=80)
    assert diagnostics["success"]
    assert fitted.tumour_growth == pytest.approx(0.11, rel=0.03)
    assert diagnostics["weighted_rmse"] < 0.1
    assert diagnostics["jacobian_rank"] == 1
    assert diagnostics["locally_identifiable"] is True
    interval = diagnostics["local_parameter_uncertainty"]["tumour_growth"]
    assert interval["approximate_lower_95"] <= truth.tumour_growth * (1 + 1e-12)
    assert interval["approximate_upper_95"] >= truth.tumour_growth * (1 - 1e-12)
    assert len(compared) == len(observations)


def test_systemic_calibration_flags_a_structurally_uninformative_parameter():
    trajectory = simulate(days=5)
    days = np.arange(0, 6, dtype=float)
    values = np.interp(days, trajectory.day, trajectory.tumour_cells)
    observations = pd.DataFrame({
        "day": days, "variable": "tumour_cells", "value": values,
        "sd": np.maximum(values * 0.01, 1.0),
    })
    _, diagnostics, _ = calibrate(
        observations, fit_parameters=("occupancy_on",), max_nfev=20,
    )
    assert diagnostics["jacobian_rank"] == 0
    assert diagnostics["locally_identifiable"] is False
    assert diagnostics["identifiability_warning"] is True
    assert diagnostics["local_parameter_uncertainty"] is None


def test_systemic_profile_likelihood_scans_fitted_parameter():
    truth = TwinParameters(tumour_growth=0.11)
    trajectory = simulate(days=8, parameters=truth)
    days = np.array([0, 2, 4, 6, 8], dtype=float)
    values = np.interp(days, trajectory.day, trajectory.tumour_cells)
    observations = pd.DataFrame({
        "day": days, "variable": "tumour_cells", "value": values,
        "sd": values * 0.02,
    })
    _, diagnostics, _ = calibrate(
        observations, fit_parameters=("tumour_growth",), max_nfev=50,
        profile_points=5, profile_fold_range=2.0,
    )
    profile = diagnostics["profile_likelihood"]
    rows = profile["profiles"]["tumour_growth"]["points"]
    assert len(rows) == 5
    assert rows[2]["delta_weighted_sse"] == pytest.approx(0.0, abs=1e-9)
    fitted_value = profile["gaussian_profile_optimum"]["tumour_growth"]
    assert rows[0]["parameter_value"] < fitted_value < rows[-1]["parameter_value"]


def test_profile_likelihood_rejects_even_grid_size():
    observations = pd.DataFrame({
        "day": [0, 1], "variable": ["tumour_cells"] * 2,
        "value": [2e7, 2.2e7], "sd": [1e6, 1e6],
    })
    with pytest.raises(ValueError, match="odd integer"):
        calibrate(observations, fit_parameters=("tumour_growth",), profile_points=4)


def test_observation_contract_rejects_unknown_state():
    bad = pd.DataFrame({"day": [0], "variable": ["magic"], "value": [1]})
    with pytest.raises(ValueError):
        validate_observations(bad)


def test_systemic_calibration_rejects_ambiguous_or_empty_observations():
    duplicate = pd.DataFrame({
        "day": [0, 0], "variable": ["tumour_cells", "tumour_cells"],
        "value": [2e7, 2e7], "sd": [1e6, 1e6],
    })
    with pytest.raises(ValueError, match="unique day-variable"):
        validate_observations(duplicate)
    with pytest.raises(ValueError, match="must not be empty"):
        validate_observations(pd.DataFrame(columns=["day", "variable", "value"]))
    missing_name = pd.DataFrame({"day": [0], "variable": [None], "value": [1]})
    with pytest.raises(ValueError, match="non-missing"):
        validate_observations(missing_name)


def test_systemic_log_space_fit_requires_unique_positive_base_parameters():
    observations = pd.DataFrame({
        "day": [0, 1], "variable": ["tumour_cells"] * 2,
        "value": [2e7, 2.2e7], "sd": [1e6, 1e6],
    })
    with pytest.raises(ValueError, match="must be unique"):
        calibrate(observations, fit_parameters=("tumour_growth", "tumour_growth"))
    with pytest.raises(ValueError, match="strictly positive"):
        calibrate(
            observations, base=TwinParameters(tumour_growth=0.0),
            fit_parameters=("tumour_growth",),
        )


def test_uncertainty_ensemble_is_reproducible():
    a = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    b = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    pd.testing.assert_frame_equal(a, b)


def test_calibration_ensemble_uses_saved_correlated_log_covariance():
    parameters = TwinParameters(tumour_growth=0.2, killing_rate=2e-8)
    diagnostics = {"locally_identifiable": True, "local_log_parameter_covariance": {
        "parameter_order": ["tumour_growth", "killing_rate"],
        "log_scale_covariance": [[0.01, 0.008], [0.008, 0.01]],
    }}
    a = calibration_ensemble(parameters, diagnostics, (), days=2, n=20, seed=7)
    b = calibration_ensemble(parameters, diagnostics, (), days=2, n=20, seed=7)
    pd.testing.assert_frame_equal(a, b)
    assert {"parameter_tumour_growth", "parameter_killing_rate"}.issubset(a)
    assert a.parameter_tumour_growth.corr(a.parameter_killing_rate) > 0.5
    assert a.tumour_cells.nunique() > 1


def test_calibration_ensemble_rejects_missing_or_invalid_covariance():
    with pytest.raises(ValueError, match="estimable covariance"):
        calibration_ensemble(TwinParameters(), {}, (), days=2, n=4)
    invalid = {"locally_identifiable": True, "local_log_parameter_covariance": {
        "parameter_order": ["tumour_growth"],
        "log_scale_covariance": [[-1.0]],
    }}
    with pytest.raises(ValueError, match="positive semidefinite"):
        calibration_ensemble(TwinParameters(), invalid, (), days=2, n=4)
    nonidentifiable = {
        "locally_identifiable": False,
        "local_log_parameter_covariance": {
            "parameter_order": ["tumour_growth"],
            "log_scale_covariance": [[0.01]],
        },
    }
    with pytest.raises(ValueError, match="locally identifiable"):
        calibration_ensemble(TwinParameters(), nonidentifiable, (), days=2, n=4)
