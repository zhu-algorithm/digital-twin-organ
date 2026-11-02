import numpy as np

import pandas as pd
import pytest

from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.spatial_dynamics import (
    SPECIES, simulate_spatial, simulate_spatial_coupled,
    simulate_spatial_counterfactual, spatial_summary,
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


def test_matched_spatial_counterfactual_has_auditable_occupancy_and_effect():
    paired = simulate_spatial_counterfactual((Dose(0, 1.0),), days=7, step=1)
    assert set(paired.arm) == {"control", "anti_pd1"}
    assert paired[paired.arm.eq("control")].pd1_occupancy.eq(0).all()
    assert paired[paired.arm.eq("anti_pd1")].pd1_occupancy.max() > 0
    endpoint = paired[(paired.species == "effector") & (paired.zone == "efferent_lymph")]
    values = endpoint.groupby("arm").amount.last()
    assert values["anti_pd1"] > values["control"]


def test_counterfactual_requires_explicit_schedule_and_zero_dose_is_equivalent():
    with pytest.raises(ValueError, match="explicit dose"):
        simulate_spatial_counterfactual((), days=2)
    paired = simulate_spatial_counterfactual((Dose(0, 0.0),), days=3, step=1)
    control = paired[paired.arm.eq("control")].reset_index(drop=True)
    treated = paired[paired.arm.eq("anti_pd1")].reset_index(drop=True)
    assert np.allclose(control.amount, treated.amount, rtol=1e-6, atol=1e-6)
    assert treated.pd1_occupancy.eq(0).all()
