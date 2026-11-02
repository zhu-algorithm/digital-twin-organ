import numpy as np

import pandas as pd
import pytest
from types import SimpleNamespace
import tdln_pd1.spatial_dynamics as spatial_module

from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.spatial_dynamics import (
    SPECIES, SpatialInitialState, SpatialParameters, simulate_spatial, simulate_spatial_coupled,
    simulate_spatial_counterfactual, spatial_summary,
)


def test_transport_only_conserves_each_species():
    trajectory = simulate_spatial(days=5, reactions=False)
    totals = trajectory.groupby(["day", "species"]).amount.sum().unstack()
    for species in SPECIES:
        assert np.allclose(totals[species], totals[species].iloc[0], rtol=1e-7, atol=1e-7)


def test_spatial_model_rejects_nonphysical_parameters_and_initial_states():
    with pytest.raises(ValueError, match="antigen_decay"):
        simulate_spatial(days=1, parameters=SpatialParameters(antigen_decay=-0.1))
    with pytest.raises(ValueError, match="afferent_dc"):
        simulate_spatial(days=1, initial=SpatialInitialState(afferent_dc=-1))
    with pytest.raises(ValueError, match="tpex_growth"):
        simulate_spatial(days=1, parameters=SpatialParameters(tpex_growth=float("nan")))
    with pytest.raises(ValueError, match="antigen_half_saturation"):
        simulate_spatial(days=1, parameters=SpatialParameters(antigen_half_saturation=0))
    with pytest.raises(ValueError, match="tpex_carrying_capacity"):
        simulate_spatial(days=1, parameters=SpatialParameters(tpex_carrying_capacity=0))


def test_antigen_half_saturation_is_an_explicit_mechanistic_control():
    baseline = simulate_spatial(days=4)
    weak_antigen_sensitivity = simulate_spatial(
        days=4, parameters=SpatialParameters(antigen_half_saturation=2.0)
    )
    selector = lambda frame: frame[(frame.species == "tpex") & (frame.zone == "paracortex")]
    assert selector(weak_antigen_sensitivity).amount.iloc[-1] < selector(baseline).amount.iloc[-1]


def test_paracortical_tpex_expansion_respects_finite_niche_capacity():
    parameters = SpatialParameters(
        antigen_decay=0.0, antigen_uptake=0.0, dc_death=0.0,
        priming_rate=0.0, tpex_growth=2.0, tpex_carrying_capacity=1.0e5,
        tpex_differentiation=0.0, tpex_death=0.0,
        tpex_to_medulla=0.0, tpex_to_follicle=0.0,
    )
    rows = [
        {"species": species, "zone": zone, "amount": (
            1.0 if species == "antigen" and zone == "paracortex" else
            8.0e4 if species == "tpex" and zone == "paracortex" else 0.0
        )}
        for species in spatial_module.SPECIES for zone in spatial_module.ZONE_NAMES
    ]
    trajectory = simulate_spatial(
        days=20, parameters=parameters, initial=pd.DataFrame(rows), step=0.2,
    )
    tpex = trajectory[
        trajectory.species.eq("tpex") & trajectory.zone.eq("paracortex")
    ].amount
    assert tpex.iloc[-1] == pytest.approx(parameters.tpex_carrying_capacity, rel=2e-4)
    assert tpex.max() <= parameters.tpex_carrying_capacity * 1.0002


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
    with pytest.raises(ValueError, match="finite"):
        simulate_spatial(days=2, checkpoint_occupancy=lambda _: float("nan"))
    with pytest.raises(ValueError, match=r"within \[0, 1\]"):
        simulate_spatial(days=2, checkpoint_occupancy=lambda _: 1.01)
    with pytest.raises(ValueError, match="numeric"):
        simulate_spatial(days=2, checkpoint_occupancy=lambda _: "invalid")
    with pytest.raises(ValueError, match="numeric"):
        simulate_spatial(days=2, checkpoint_occupancy="invalid")


def test_coupling_never_silently_drops_invalid_systemic_timeline_rows():
    missing = pd.DataFrame({
        "day": [0.0, 1.0, 2.0],
        "pd1_occupancy": [0.0, float("nan"), 0.2],
    })
    with pytest.raises(ValueError, match="complete and finite"):
        simulate_spatial_coupled(missing, days=2)
    non_numeric = pd.DataFrame({
        "day": [0.0, "invalid", 2.0],
        "pd1_occupancy": [0.0, 0.1, 0.2],
    })
    with pytest.raises(ValueError, match="complete and finite"):
        simulate_spatial_coupled(non_numeric, days=2)
    duplicate = pd.DataFrame({
        "day": [0.0, 1.0, 1.0, 2.0],
        "pd1_occupancy": [0.0, 0.1, 0.2, 0.3],
    })
    with pytest.raises(ValueError, match="unique"):
        simulate_spatial_coupled(duplicate, days=2)


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


def test_spatial_solver_cannot_hide_material_negative_states(monkeypatch):
    def invalid_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(SPECIES) * len(spatial_module.ZONE_NAMES), len(t_eval)))
        values[0, -1] = -1.0
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(spatial_module, "solve_ivp", invalid_solver)
    with pytest.raises(RuntimeError, match="non-negative state bounds"):
        simulate_spatial(days=1)


def test_spatial_solver_clips_only_tolerance_scale_negative_drift(monkeypatch):
    def tiny_drift_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(SPECIES) * len(spatial_module.ZONE_NAMES), len(t_eval)))
        values[0, -1] = -1e-10
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(spatial_module, "solve_ivp", tiny_drift_solver)
    assert simulate_spatial(days=1).amount.min() == 0.0
