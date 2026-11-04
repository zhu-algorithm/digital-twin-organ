"""Closed-loop liver-tumour ↔ TDLN/anti-PD-1 simulation.

The liver module supplies viable tumour burden and antigen to the systemic
immune model.  TDLN-derived tumour effectors feed back as killing at the viable
tumour rim.  This is an explicit research hypothesis layer: the dimensional
bridges are declared parameters and require cohort calibration.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .liver_tumor import (LiverTumourInitialState, LiverTumourParameters,
                          STATE_NAMES as LIVER_STATES, _rhs as _liver_rhs,
                          _inflow as _liver_inflow,
                          _validate as _validate_liver,
                          _validate_parameters as _validate_liver_parameters)
from .simulator import (Dose, InitialState, STATE_NAMES as IMMUNE_STATES,
                        TwinParameters, _rhs as _immune_rhs,
                        _validate_model_inputs)

IMMUNE_DYNAMIC_STATES = IMMUNE_STATES[1:]


@dataclass(frozen=True)
class LiverTDLNCouplingParameters:
    tumour_cells_per_fraction: float = 3.0e8
    immune_killing_gain: float = 0.004
    necrotic_antigen_gain: float = 0.35
    hypoxic_antigen_gain: float = 0.20


def _validate_coupling(p: LiverTDLNCouplingParameters) -> None:
    values = np.asarray(list(asdict(p).values()), dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("coupling parameters must be finite and non-negative")
    if p.tumour_cells_per_fraction <= 0:
        raise ValueError("tumour_cells_per_fraction must be strictly positive")


def _combined_rhs(time: float, y: np.ndarray, liver_p: LiverTumourParameters,
                  immune_p: TwinParameters, coupling: LiverTDLNCouplingParameters,
                  flow_multiplier: float, angiogenesis_inhibition: float,
                  portal_flow_multiplier: float | None,
                  arterial_flow_multiplier: float | None) -> np.ndarray:
    liver = np.maximum(y[:len(LIVER_STATES)], 0.0)
    immune_dynamic = np.maximum(y[len(LIVER_STATES):], 0.0)
    rim = liver[LIVER_STATES.index("tumour_rim")]
    core = liver[LIVER_STATES.index("tumour_core")]
    necrotic = liver[LIVER_STATES.index("necrotic_tumour")]
    oxygen_core = liver[LIVER_STATES.index("oxygen_core")]
    viable_cells = (rim + core) * coupling.tumour_cells_per_fraction

    full_immune = np.concatenate(([viable_cells], immune_dynamic))
    immune_dydt = _immune_rhs(time, full_immune, immune_p)
    liver_dydt = _liver_rhs(time, liver, liver_p, flow_multiplier,
                            angiogenesis_inhibition, portal_flow_multiplier,
                            arterial_flow_multiplier)

    occupancy = np.clip(full_immune[IMMUNE_STATES.index("pd1_occupancy")], 0.0, 1.0)
    effector = full_immune[IMMUNE_STATES.index("effector_tumour")]
    checkpoint_relief = 0.25 + 0.75 * occupancy
    killed_cells = immune_p.killing_rate * effector * viable_cells * checkpoint_relief
    killed_fraction = coupling.immune_killing_gain * killed_cells / coupling.tumour_cells_per_fraction
    liver_dydt[LIVER_STATES.index("tumour_rim")] -= min(killed_fraction, rim)
    liver_dydt[LIVER_STATES.index("necrotic_tumour")] += min(killed_fraction, rim)

    # Damaged and hypoxic liver tumour contributes additional drainable antigen.
    antigen_index = IMMUNE_STATES.index("tumour_antigen") - 1
    immune_dydt[1:][antigen_index] += (
        coupling.necrotic_antigen_gain * necrotic
        + coupling.hypoxic_antigen_gain * (rim + core) * max(0.0, 0.45 - oxygen_core)
    )
    out = np.concatenate((liver_dydt, immune_dydt[1:]))
    return np.where((y <= 0.0) & (out < 0.0), 0.0, out)


def simulate_liver_tdln(*, days: float = 84.0, step: float = 0.25,
                        doses: tuple[Dose, ...] = (),
                        liver_parameters: LiverTumourParameters | None = None,
                        immune_parameters: TwinParameters | None = None,
                        coupling_parameters: LiverTDLNCouplingParameters | None = None,
                        liver_initial: LiverTumourInitialState | None = None,
                        immune_initial: InitialState | None = None,
                        flow_multiplier: float = 1.0,
                        portal_flow_multiplier: float | None = None,
                        arterial_flow_multiplier: float | None = None,
                        angiogenesis_inhibition: float = 0.0) -> pd.DataFrame:
    """Simulate the simultaneous liver tumour–TDLN–anti-PD-1 feedback loop."""
    lp, ip = liver_parameters or LiverTumourParameters(), immune_parameters or TwinParameters()
    cp = coupling_parameters or LiverTDLNCouplingParameters()
    li, ii = liver_initial or LiverTumourInitialState(), immune_initial or InitialState()
    _validate_liver_parameters(lp)
    _validate_liver(li, "liver initial state")
    _validate_model_inputs(ip, ii)
    _validate_coupling(cp)
    if not np.isfinite(days) or not np.isfinite(step) or days <= 0 or step <= 0:
        raise ValueError("days and step must be finite and strictly positive")
    if not np.isfinite(flow_multiplier) or flow_multiplier <= 0:
        raise ValueError("flow_multiplier must be finite and strictly positive")
    for name, value in (("portal_flow_multiplier", portal_flow_multiplier),
                        ("arterial_flow_multiplier", arterial_flow_multiplier)):
        if value is not None and (not np.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be finite and strictly positive")
    if not np.isfinite(angiogenesis_inhibition) or not 0 <= angiogenesis_inhibition <= 1:
        raise ValueError("angiogenesis_inhibition must lie in [0, 1]")
    if any(not np.isfinite(d.day) or not np.isfinite(d.amount) or d.day < 0
           or d.day >= days or d.amount < 0 for d in doses):
        raise ValueError("doses must have non-negative amount and occur before the simulation endpoint")

    events: dict[float, float] = {}
    for dose in doses:
        events[dose.day] = events.get(dose.day, 0.0) + dose.amount
    boundaries = sorted({0.0, days, *events})
    immune_values = np.asarray([getattr(ii, name) for name in IMMUNE_DYNAMIC_STATES])
    y = np.concatenate((np.asarray(list(asdict(li).values())), immune_values))
    plasma_index = len(LIVER_STATES) + IMMUNE_DYNAMIC_STATES.index("anti_pd1_plasma")
    frames = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if start in events:
            y[plasma_index] += events[start]
        times = np.unique(np.append(np.arange(start, end, step), end))
        solved = solve_ivp(_combined_rhs, (start, end), y, t_eval=times,
                           # Both inflows remain independently executable inside
                           # the simultaneous tumour–immune feedback solve.
                           args=(lp, ip, cp, flow_multiplier, angiogenesis_inhibition,
                                 portal_flow_multiplier, arterial_flow_multiplier),
                           rtol=1e-7, atol=1e-9)
        if not solved.success or not np.isfinite(solved.y).all():
            raise RuntimeError(f"coupled liver–TDLN integration failed: {solved.message}")
        values = np.maximum(solved.y.T, 0.0)
        frame = pd.DataFrame(values, columns=[f"liver_{x}" for x in LIVER_STATES]
                             + list(IMMUNE_DYNAMIC_STATES))
        frame.insert(0, "day", solved.t)
        if frames:
            frames[-1] = frames[-1].iloc[:-1] if start in events else frames[-1]
            if start not in events:
                frame = frame.iloc[1:]
        frames.append(frame)
        y = values[-1]
    out = pd.concat(frames, ignore_index=True)
    out.insert(1, "administered_dose", out.day.map(events).fillna(0.0))
    out["tumour_cells"] = ((out.liver_tumour_rim + out.liver_tumour_core)
                           * cp.tumour_cells_per_fraction)
    out["liver_total_tumour"] = (out.liver_tumour_rim + out.liver_tumour_core
                                  + out.liver_necrotic_tumour)
    out["liver_tumour_fraction"] = out.liver_total_tumour / (
        out.liver_normal_cells + out.liver_total_tumour).clip(lower=1e-12)
    portal, arterial, buffer_flow, oxygenated = _liver_inflow(
        lp, flow_multiplier, portal_flow_multiplier, arterial_flow_multiplier)
    out["liver_perfusion"] = ((portal + arterial)
                               * out.liver_functional_vessels
                               / (1.0 + lp.vessel_compression
                                  * (out.liver_tumour_rim + out.liver_tumour_core)))
    out["liver_portal_flow"] = portal * out.liver_functional_vessels
    out["liver_arterial_flow"] = arterial * out.liver_functional_vessels
    out["liver_arterial_buffer_flow"] = buffer_flow * out.liver_functional_vessels
    out["liver_oxygenated_inflow"] = oxygenated * out.liver_functional_vessels
    return out

