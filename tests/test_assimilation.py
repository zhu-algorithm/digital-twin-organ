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
    assert len(compared) == len(observations)


def test_observation_contract_rejects_unknown_state():
    bad = pd.DataFrame({"day": [0], "variable": ["magic"], "value": [1]})
    with pytest.raises(ValueError):
        validate_observations(bad)


def test_uncertainty_ensemble_is_reproducible():
    a = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    b = ensemble(TwinParameters(), (), days=2, n=4, seed=7)
    pd.testing.assert_frame_equal(a, b)
