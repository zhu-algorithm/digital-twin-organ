import numpy as np
import pytest
from types import SimpleNamespace

import tdln_pd1.liver_tumor as liver_module
from tdln_pd1.liver_tumor import (LiverTumourInitialState, LiverTumourParameters,
                                  STATE_NAMES, simulate_liver_tumour)


def test_liver_tumour_simulation_is_finite_nonnegative_and_grows():
    frame = simulate_liver_tumour(days=60)
    assert np.isfinite(frame.to_numpy()).all()
    assert (frame.drop(columns="day") >= 0).all().all()
    assert frame.total_tumour.iloc[-1] > frame.total_tumour.iloc[0]
    assert frame.normal_cells.iloc[-1] < frame.normal_cells.iloc[0]


def test_reduced_blood_flow_reduces_oxygen_and_increases_necrosis():
    normal = simulate_liver_tumour(days=90)
    low = simulate_liver_tumour(days=90, flow_multiplier=0.55)
    assert low.oxygen_core.iloc[-1] < normal.oxygen_core.iloc[-1]
    assert low.necrotic_tumour.iloc[-1] > normal.necrotic_tumour.iloc[-1]


def test_portal_and_arterial_inflow_are_independently_executable():
    # Portal loss triggers partial arterial buffering, whereas an explicitly
    # measured arterial loss is never inferred away.
    portal_loss = simulate_liver_tumour(days=60, portal_flow_multiplier=2 / 3)
    arterial_loss = simulate_liver_tumour(days=60, arterial_flow_multiplier=0.001)
    assert portal_loss.portal_flow.iloc[0] < simulate_liver_tumour(days=1).portal_flow.iloc[0]
    assert arterial_loss.arterial_flow.iloc[0] < portal_loss.arterial_flow.iloc[0]
    assert arterial_loss.oxygenated_inflow.iloc[0] < portal_loss.oxygenated_inflow.iloc[0]
    assert arterial_loss.oxygen_core.iloc[-1] < portal_loss.oxygen_core.iloc[-1]
    assert portal_loss.arterial_buffer_flow.iloc[0] > 0
    assert np.allclose(portal_loss.portal_flow + portal_loss.arterial_flow,
                       (0.75 * 2 / 3 + 0.25 + 0.25 * 0.75 / 3)
                       * portal_loss.functional_vessels)


def test_hepatic_arterial_buffer_is_partial_bounded_and_not_double_counted():
    buffered = simulate_liver_tumour(days=2, portal_flow_multiplier=0.4)
    disabled = simulate_liver_tumour(
        days=2, portal_flow_multiplier=0.4,
        parameters=LiverTumourParameters(arterial_buffer_fraction=0.0))
    both_measured = simulate_liver_tumour(
        days=2, portal_flow_multiplier=0.4, arterial_flow_multiplier=1.0)
    assert buffered.arterial_flow.iloc[0] > disabled.arterial_flow.iloc[0]
    assert buffered.perfusion.iloc[0] < simulate_liver_tumour(days=2).perfusion.iloc[0]
    assert both_measured.arterial_buffer_flow.max() == pytest.approx(0.0)
    capped = simulate_liver_tumour(
        days=2, portal_flow_multiplier=0.01,
        parameters=LiverTumourParameters(
            arterial_buffer_fraction=1.0, arterial_buffer_max_multiplier=1.2))
    assert capped.arterial_flow.iloc[0] <= (
        0.25 * 1.2 * capped.functional_vessels.iloc[0] + 1e-12)


def test_legacy_total_flow_multiplier_equals_two_explicit_inflow_multipliers():
    legacy = simulate_liver_tumour(days=10, flow_multiplier=0.65)
    explicit = simulate_liver_tumour(
        days=10, portal_flow_multiplier=0.65, arterial_flow_multiplier=0.65)
    assert np.allclose(legacy.to_numpy(), explicit.to_numpy(), rtol=1e-10, atol=1e-12)


def test_zero_tumour_does_not_spontaneously_create_tumour():
    state = LiverTumourInitialState(tumour_rim=0, tumour_core=0, necrotic_tumour=0)
    frame = simulate_liver_tumour(days=30, initial=state)
    assert frame.total_tumour.max() == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("kwargs", ({"days": 0}, {"step": -1}, {"flow_multiplier": 0},
                                    {"portal_flow_multiplier": 0},
                                    {"arterial_flow_multiplier": float("nan")},
                                    {"angiogenesis_inhibition": 1.2},
                                    {"parameters": LiverTumourParameters(portal_flow=0)},
                                    {"parameters": LiverTumourParameters(
                                        arterial_buffer_fraction=1.2)},
                                    {"parameters": LiverTumourParameters(
                                        arterial_buffer_max_multiplier=0.9)}))
def test_invalid_simulation_contracts_fail(kwargs):
    with pytest.raises(ValueError):
        simulate_liver_tumour(**kwargs)


def test_liver_solver_cannot_hide_material_state_bound_violations(monkeypatch):
    def invalid_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(STATE_NAMES), len(t_eval)))
        values[STATE_NAMES.index("normal_cells"), -1] = -1.0
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(liver_module, "solve_ivp", invalid_solver)
    with pytest.raises(RuntimeError, match="physical state bounds"):
        simulate_liver_tumour(days=1)


def test_liver_solver_clips_only_tolerance_scale_drift(monkeypatch):
    def tiny_drift_solver(*args, t_eval, **kwargs):
        values = np.zeros((len(STATE_NAMES), len(t_eval)))
        values[STATE_NAMES.index("normal_cells"), -1] = -1e-11
        return SimpleNamespace(success=True, message="", t=np.asarray(t_eval), y=values)

    monkeypatch.setattr(liver_module, "solve_ivp", tiny_drift_solver)
    out = simulate_liver_tumour(days=1)
    assert out.normal_cells.min() == 0.0
