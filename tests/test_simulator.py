import numpy as np
import pytest
from tdln_pd1.simulator import (Dose, InitialState, STATE_NAMES, TwinParameters,
                                simulate)


def test_simulator_is_nonnegative_and_has_expected_schema():
    out = simulate(days=5, step=0.5)
    assert set(STATE_NAMES).issubset(out.columns)
    assert np.isfinite(out[list(STATE_NAMES)].to_numpy()).all()
    assert (out[list(STATE_NAMES)] >= 0).all().all()


def test_simulator_rejects_nonphysical_parameters_and_initial_states():
    with pytest.raises(ValueError, match="tumour_growth"):
        simulate(days=1, parameters=TwinParameters(tumour_growth=-0.1))
    with pytest.raises(ValueError, match="tumour_capacity"):
        simulate(days=1, parameters=TwinParameters(tumour_capacity=0))
    with pytest.raises(ValueError, match="tumour_cells"):
        simulate(days=1, initial=InitialState(tumour_cells=-1))
    with pytest.raises(ValueError, match="pd1_occupancy"):
        simulate(days=1, initial=InitialState(pd1_occupancy=1.1))


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


def test_bolus_boundary_is_post_dose_and_auditable():
    trajectory = simulate(days=2, doses=(Dose(1, 1.5),), step=0.25)
    boundary = trajectory.loc[trajectory.day.eq(1.0)].iloc[0]
    assert boundary.administered_dose == 1.5
    assert boundary.anti_pd1_plasma == pytest.approx(1.5)
    assert trajectory.day.is_unique


def test_dose_at_simulation_endpoint_is_rejected_instead_of_ignored():
    with pytest.raises(ValueError, match="before the simulation endpoint"):
        simulate(days=2, doses=(Dose(2, 1.0),))
