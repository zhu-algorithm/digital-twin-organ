"""Mechanistic liver-tumour microenvironment simulator.

The model is deliberately compact: it is a research/demo digital-organ layer,
not a clinical growth or treatment predictor.  States are normalised tissue
fractions (cells and vessels) or oxygen saturation fractions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class LiverTumourParameters:
    portal_flow: float = 0.75
    arterial_flow: float = 0.25
    portal_oxygen_content: float = 0.75
    arterial_oxygen_content: float = 1.0
    arterial_buffer_fraction: float = 0.25
    arterial_buffer_max_multiplier: float = 2.5
    normal_turnover: float = 0.045
    normal_capacity: float = 1.0
    tumour_growth: float = 0.085
    tumour_capacity: float = 0.72
    core_conversion: float = 0.030
    hypoxic_necrosis: float = 0.050
    necrotic_clearance: float = 0.018
    angiogenesis: float = 0.060
    vessel_remodelling: float = 0.035
    vessel_compression: float = 1.10
    oxygen_delivery: float = 0.90
    oxygen_turnover: float = 0.65
    normal_oxygen_use: float = 0.24
    tumour_oxygen_use: float = 0.48


@dataclass(frozen=True)
class LiverTumourInitialState:
    normal_cells: float = 0.88
    tumour_rim: float = 0.045
    tumour_core: float = 0.012
    necrotic_tumour: float = 0.003
    functional_vessels: float = 0.92
    oxygen_normal: float = 0.92
    oxygen_rim: float = 0.74
    oxygen_core: float = 0.42


STATE_NAMES = tuple(field.name for field in fields(LiverTumourInitialState))


def _validate(values: object, label: str) -> None:
    array = np.asarray(list(asdict(values).values()), dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"{label} must contain finite, non-negative values")


def _inflow(p: LiverTumourParameters, flow_multiplier: float,
            portal_flow_multiplier: float | None,
            arterial_flow_multiplier: float | None) -> tuple[float, float, float, float]:
    """Return portal, arterial, buffer and normalised oxygenated inflow."""
    portal_scale = flow_multiplier if portal_flow_multiplier is None else portal_flow_multiplier
    arterial_scale = flow_multiplier if arterial_flow_multiplier is None else arterial_flow_multiplier
    portal = p.portal_flow * portal_scale
    unbuffered_arterial = p.arterial_flow * arterial_scale
    buffer_flow = 0.0
    # Infer HABR only for an isolated portal perturbation. If both flows were
    # measured/supplied, the arterial measurement already contains any response.
    if portal_flow_multiplier is not None and arterial_flow_multiplier is None:
        portal_deficit = max(0.0, p.portal_flow * flow_multiplier - portal)
        arterial_cap = p.arterial_flow * flow_multiplier * p.arterial_buffer_max_multiplier
        buffered = min(arterial_cap,
                       unbuffered_arterial + p.arterial_buffer_fraction * portal_deficit)
        buffer_flow = max(0.0, buffered - unbuffered_arterial)
    arterial = unbuffered_arterial + buffer_flow
    baseline_oxygenated = (p.portal_flow * p.portal_oxygen_content
                           + p.arterial_flow * p.arterial_oxygen_content)
    oxygenated = (portal * p.portal_oxygen_content
                  + arterial * p.arterial_oxygen_content) / baseline_oxygenated
    return portal, arterial, buffer_flow, oxygenated


def _rhs(_time: float, y: np.ndarray, p: LiverTumourParameters,
         flow_multiplier: float, angiogenesis_inhibition: float,
         portal_flow_multiplier: float | None = None,
         arterial_flow_multiplier: float | None = None) -> np.ndarray:
    h, rim, core, nec, vessel, o_h, o_r, o_c = np.maximum(y, 0.0)
    tumour = rim + core
    portal, arterial, _buffer_flow, oxygenated_inflow = _inflow(
        p, flow_multiplier, portal_flow_multiplier, arterial_flow_multiplier)
    total_flow = portal + arterial
    perfusion = total_flow * vessel / (1.0 + p.vessel_compression * tumour)
    oxygen_perfusion = oxygenated_inflow * vessel / (
        1.0 + p.vessel_compression * tumour)
    rim_perf = oxygen_perfusion * (0.68 + 0.32 * vessel)
    core_perf = oxygen_perfusion * np.exp(-2.4 * core)

    space_h = max(0.0, 1.0 - h / p.normal_capacity - 0.55 * tumour)
    space_t = max(0.0, 1.0 - tumour / p.tumour_capacity - 0.32 * h)
    hypoxia_r = max(0.0, 0.55 - o_r)
    hypoxia_c = max(0.0, 0.42 - o_c)

    dh = p.normal_turnover * h * space_h * o_h - 0.022 * tumour * h
    drim = p.tumour_growth * rim * space_t * (0.35 + 0.65 * o_r) - p.core_conversion * rim * (0.4 + hypoxia_r)
    dcore = p.core_conversion * rim * (0.4 + hypoxia_r) - p.hypoxic_necrosis * core * (0.15 + 2.0 * hypoxia_c)
    dnec = p.hypoxic_necrosis * core * (0.15 + 2.0 * hypoxia_c) - p.necrotic_clearance * nec

    angiogenic_signal = tumour * (0.25 + hypoxia_r + hypoxia_c)
    dvessel = (p.angiogenesis * (1.0 - angiogenesis_inhibition) * angiogenic_signal * (1.0 - vessel)
               - p.vessel_remodelling * tumour * vessel)
    do_h = p.oxygen_turnover * (p.oxygen_delivery * oxygen_perfusion - o_h) - p.normal_oxygen_use * h * o_h
    do_r = p.oxygen_turnover * (p.oxygen_delivery * rim_perf - o_r) - p.tumour_oxygen_use * rim * o_r
    do_c = p.oxygen_turnover * (p.oxygen_delivery * core_perf - o_c) - 1.35 * p.tumour_oxygen_use * core * o_c
    return np.asarray((dh, drim, dcore, dnec, dvessel, do_h, do_r, do_c))


def simulate_liver_tumour(*, days: float = 120.0, step: float = 0.5,
                          parameters: LiverTumourParameters | None = None,
                          initial: LiverTumourInitialState | None = None,
                          flow_multiplier: float = 1.0,
                          portal_flow_multiplier: float | None = None,
                          arterial_flow_multiplier: float | None = None,
                          angiogenesis_inhibition: float = 0.0) -> pd.DataFrame:
    """Simulate liver parenchyma, tumour, perfusion and oxygen trajectories."""
    p = parameters or LiverTumourParameters()
    state = initial or LiverTumourInitialState()
    _validate(p, "parameters")
    _validate(state, "initial state")
    if (p.portal_flow <= 0 or p.arterial_flow <= 0
            or p.portal_oxygen_content <= 0 or p.arterial_oxygen_content <= 0):
        raise ValueError("dual-inflow rates and oxygen contents must be strictly positive")
    if not 0 <= p.arterial_buffer_fraction <= 1:
        raise ValueError("arterial_buffer_fraction must lie in [0, 1]")
    if p.arterial_buffer_max_multiplier < 1:
        raise ValueError("arterial_buffer_max_multiplier must be at least 1")
    if not np.isfinite(days) or days <= 0 or not np.isfinite(step) or step <= 0:
        raise ValueError("days and step must be finite and strictly positive")
    if not np.isfinite(flow_multiplier) or flow_multiplier <= 0:
        raise ValueError("flow_multiplier must be finite and strictly positive")
    for name, value in (("portal_flow_multiplier", portal_flow_multiplier),
                        ("arterial_flow_multiplier", arterial_flow_multiplier)):
        if value is not None and (not np.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be finite and strictly positive")
    if not np.isfinite(angiogenesis_inhibition) or not 0 <= angiogenesis_inhibition <= 1:
        raise ValueError("angiogenesis_inhibition must lie in [0, 1]")
    t_eval = np.append(np.arange(0.0, days, step), days)
    solved = solve_ivp(_rhs, (0.0, days), list(asdict(state).values()), t_eval=t_eval,
                       args=(p, flow_multiplier, angiogenesis_inhibition,
                             portal_flow_multiplier, arterial_flow_multiplier),
                       method="LSODA", rtol=1e-8, atol=1e-10)
    if not solved.success or not np.all(np.isfinite(solved.y)):
        raise RuntimeError(f"liver tumour integration failed: {solved.message}")
    values = np.maximum(solved.y.T, 0.0)
    frame = pd.DataFrame(values, columns=STATE_NAMES)
    frame.insert(0, "day", solved.t)
    frame["viable_tumour"] = frame.tumour_rim + frame.tumour_core
    frame["total_tumour"] = frame.viable_tumour + frame.necrotic_tumour
    portal, arterial, buffer_flow, oxygenated_inflow = _inflow(
        p, flow_multiplier, portal_flow_multiplier, arterial_flow_multiplier)
    total_flow = portal + arterial
    frame["perfusion"] = total_flow * frame.functional_vessels / (
        1.0 + p.vessel_compression * frame.viable_tumour)
    frame["portal_flow"] = portal * frame.functional_vessels
    frame["arterial_flow"] = arterial * frame.functional_vessels
    frame["arterial_buffer_flow"] = buffer_flow * frame.functional_vessels
    frame["oxygenated_inflow"] = oxygenated_inflow * frame.functional_vessels
    frame["arterial_flow_fraction"] = arterial / max(total_flow, 1e-12)
    frame["tumour_fraction"] = frame.total_tumour / (
        frame.normal_cells + frame.total_tumour).clip(lower=1e-12)
    return frame


def parameter_manifest() -> dict[str, float]:
    return asdict(LiverTumourParameters())

