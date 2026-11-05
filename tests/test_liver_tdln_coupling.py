import numpy as np
import pytest
from scipy.integrate import trapezoid
from types import SimpleNamespace

import tdln_pd1.liver_tdln_coupling as coupling_module
from tdln_pd1.liver_tdln_coupling import (LiverTDLNCouplingParameters,
                                          COMBINED_STATE_NAMES,
                                          simulate_liver_tdln)
from tdln_pd1.simulator import Dose


def test_closed_loop_is_finite_nonnegative_and_exposes_both_organs():
    frame = simulate_liver_tdln(days=14, step=0.5)
    assert np.isfinite(frame.select_dtypes("number").to_numpy()).all()
    assert (frame.drop(columns="day").select_dtypes("number") >= 0).all().all()
    for column in ("liver_normal_cells", "liver_oxygen_core", "tumour_antigen",
                   "mature_dc_tdln", "tpex_tdln", "effector_tumour", "tumour_cells"):
        assert column in frame


def test_liver_burden_drives_tdln_antigen_and_dendritic_response():
    frame = simulate_liver_tdln(days=28)
    assert frame.tumour_antigen.max() > frame.tumour_antigen.iloc[0]
    assert frame.mature_dc_tdln.max() > 0
    assert frame.tumour_cells.iloc[-1] != pytest.approx(frame.tumour_cells.iloc[0])


def test_anti_pd1_closes_the_loop_back_to_lower_viable_liver_tumour():
    control = simulate_liver_tdln(days=56)
    treated = simulate_liver_tdln(days=56, doses=(Dose(0, 1.0), Dose(21, 1.0), Dose(42, 1.0)))
    assert treated.pd1_occupancy.max() > control.pd1_occupancy.max()
    assert trapezoid(treated.effector_tumour, treated.day) > trapezoid(
        control.effector_tumour, control.day)
    assert treated.tumour_cells.iloc[-1] < control.tumour_cells.iloc[-1]


def test_zero_killing_gain_removes_treatment_effect_on_liver_burden():
    coupling = LiverTDLNCouplingParameters(immune_killing_gain=0.0)
    control = simulate_liver_tdln(days=28, coupling_parameters=coupling)
    treated = simulate_liver_tdln(days=28, doses=(Dose(0, 1.0),), coupling_parameters=coupling)
    assert np.allclose(control.tumour_cells, treated.tumour_cells, rtol=1e-7, atol=1e-6)


def test_independent_dual_inflow_propagates_through_closed_loop():
    portal_loss = simulate_liver_tdln(days=28, portal_flow_multiplier=2 / 3)
    arterial_loss = simulate_liver_tdln(days=28, arterial_flow_multiplier=0.001)
    assert arterial_loss.liver_oxygenated_inflow.iloc[0] < (
        portal_loss.liver_oxygenated_inflow.iloc[0])
    assert arterial_loss.liver_oxygen_core.iloc[-1] < portal_loss.liver_oxygen_core.iloc[-1]
    assert not np.allclose(arterial_loss.tumour_antigen, portal_loss.tumour_antigen)


def test_invalid_coupling_contract_fails_closed():
    with pytest.raises(ValueError, match="tumour_cells_per_fraction"):
        simulate_liver_tdln(days=2, coupling_parameters=LiverTDLNCouplingParameters(
            tumour_cells_per_fraction=0))
    with pytest.raises(ValueError, match="portal_flow_multiplier"):
        simulate_liver_tdln(days=2, portal_flow_multiplier=0)
    with pytest.raises(ValueError, match="solver tolerances"):
        simulate_liver_tdln(days=2, rtol=0)


def test_coupled_solver_cannot_hide_material_state_bound_violations(monkeypatch):
    def invalid_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(COMBINED_STATE_NAMES), len(t_eval)))
        values[COMBINED_STATE_NAMES.index("liver_normal_cells"), -1] = -1.0
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(coupling_module, "solve_ivp", invalid_solver)
    with pytest.raises(RuntimeError, match="physical state bounds"):
        simulate_liver_tdln(days=1)


def test_coupled_solver_clips_only_tolerance_scale_drift(monkeypatch):
    def tiny_drift_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(COMBINED_STATE_NAMES), len(t_eval)))
        values[COMBINED_STATE_NAMES.index("liver_normal_cells"), -1] = -1e-10
        values[COMBINED_STATE_NAMES.index("pd1_occupancy"), -1] = 1.0 + 1e-10
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(coupling_module, "solve_ivp", tiny_drift_solver)
    out = simulate_liver_tdln(days=1)
    assert out.liver_normal_cells.min() == 0.0
    assert out.pd1_occupancy.max() == 1.0
