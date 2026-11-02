import numpy as np

import pandas as pd
import pytest

from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.spatial_dynamics import (
    SPECIES, simulate_spatial, simulate_spatial_coupled, spatial_summary,
)


def test_transport_only_conserves_each_species():
    trajectory = simulate_spatial(days=5, reactions=False)
    totals = trajectory.groupby(["day", "species"]).amount.sum().unstack()
    for species in SPECIES:
        assert np.allclose(totals[species], totals[species].iloc[0], rtol=1e-7, atol=1e-7)


def test_antigen_reaches_sinus_and_dc_reaches_paracortex():
    trajectory = simulate_spatial(days=4)
    summary = spatial_summary(trajectory)
    assert summary[("antigen", "subcapsular_sinus")].max() > 0
    assert summary[("mature_dc", "paracortex")].max() > 0


def test_effectors_exit_through_efferent_lymph_and_states_are_nonnegative():
    trajectory = simulate_spatial(days=6)
    efferent = trajectory[(trajectory.species == "effector") & (trajectory.zone == "efferent_lymph")]
    assert efferent.amount.iloc[-1] > 0
    assert trajectory.amount.ge(0).all()


def test_zero_occupancy_coupling_is_exactly_the_uncoupled_baseline():
    systemic = simulate(days=6)
    baseline = simulate_spatial(days=6)
    coupled = simulate_spatial_coupled(systemic, days=6)
    assert np.allclose(coupled.amount, baseline.amount, rtol=1e-10, atol=1e-10)


def test_anti_pd1_occupancy_increases_spatial_effector_output():
    control = simulate_spatial_coupled(simulate(days=14), days=14)
    treated = simulate_spatial_coupled(simulate(days=14, doses=(Dose(0, 1.0),)), days=14)
    selector = lambda frame: frame[(frame.species == "effector") & (frame.zone == "efferent_lymph")]
    assert selector(treated).amount.iloc[-1] > selector(control).amount.iloc[-1]
    tpex = lambda frame: frame[(frame.species == "tpex") & (frame.zone == "paracortex")]
    assert tpex(treated).amount.iloc[-1] > tpex(control).amount.iloc[-1]


def test_coupling_rejects_invalid_or_incomplete_occupancy():
    with pytest.raises(ValueError, match="cover"):
        simulate_spatial_coupled(pd.DataFrame({"day": [0, 1], "pd1_occupancy": [0, 0.5]}), days=2)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        simulate_spatial(days=2, checkpoint_occupancy=1.1)
