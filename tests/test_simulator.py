import numpy as np
from tdln_pd1.simulator import Dose, STATE_NAMES, simulate


def test_simulator_is_nonnegative_and_has_expected_schema():
    out = simulate(days=5, step=0.5)
    assert set(STATE_NAMES).issubset(out.columns)
    assert np.isfinite(out[list(STATE_NAMES)].to_numpy()).all()
    assert (out[list(STATE_NAMES)] >= 0).all().all()


def test_antibody_reaches_tdln_and_binds_pd1():
    treated = simulate(days=10, doses=(Dose(0, 1.0),), step=0.5)
    control = simulate(days=10, step=0.5)
    assert treated.anti_pd1_tdln.max() > 0
    assert treated.pd1_occupancy.max() > control.pd1_occupancy.max()
    assert treated.effector_tumour.iloc[-1] > control.effector_tumour.iloc[-1]


def test_duplicate_doses_on_same_day_are_combined():
    a = simulate(days=2, doses=(Dose(0, 1), Dose(0, 2)))
    b = simulate(days=2, doses=(Dose(0, 3),))
    assert np.allclose(a.anti_pd1_plasma, b.anti_pd1_plasma)
