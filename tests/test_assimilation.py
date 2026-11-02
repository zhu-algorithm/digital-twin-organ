import numpy as np
import pandas as pd
import pytest

from tdln_pd1.assimilation import calibrate, ensemble, validate_observations
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


def test_uncertainty_ensemble_is_reproducible():
    a = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    b = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    pd.testing.assert_frame_equal(a, b)
