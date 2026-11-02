"""Mechanistic tumour–blood–TDLN digital twin.

The model is deliberately compact enough to audit.  It is not a clinical
physiology claim: parameters are configurable hypotheses and outputs require
calibration against a matched cohort before patient interpretation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

STATE_NAMES = (
    "tumour_cells", "tumour_antigen", "mature_dc_tdln", "naive_cd8_tdln",
    "tpex_tdln", "effector_tdln", "effector_tumour", "exhausted_tumour",
    "treg_tdln", "anti_pd1_plasma", "anti_pd1_tdln", "pd1_occupancy",
)


@dataclass(frozen=True)
class TwinParameters:
    # Time is days; cell states are scaled cell counts; drug is arbitrary exposure.
    tumour_growth: float = 0.18
    tumour_capacity: float = 1.0e9
    antigen_release: float = 2.0e-8
    antigen_clearance: float = 0.8
    dc_migration: float = 0.35
    dc_turnover: float = 0.25
    priming_rate: float = 1.8e-8
    tpex_self_renewal: float = 0.12
    tpex_to_effector: float = 0.18
    effector_egress: float = 0.35
    effector_death: float = 0.16
    tumour_infiltration: float = 0.28
    exhaustion_rate: float = 0.22
    reinvigoration_rate: float = 0.16
    killing_rate: float = 1.2e-8
    treg_suppression: float = 2.0e-7
    antibody_clearance: float = 0.09
    antibody_tdln_exchange: float = 0.35
    antibody_tdln_loss: float = 0.12
    occupancy_on: float = 1.4
    occupancy_off: float = 0.18


@dataclass(frozen=True)
class InitialState:
    tumour_cells: float = 2.0e7
    tumour_antigen: float = 0.2
    mature_dc_tdln: float = 2.0e4
    naive_cd8_tdln: float = 2.0e6
    tpex_tdln: float = 8.0e4
    effector_tdln: float = 3.0e4
    effector_tumour: float = 2.0e4
    exhausted_tumour: float = 1.5e5
    treg_tdln: float = 5.0e4
    anti_pd1_plasma: float = 0.0
    anti_pd1_tdln: float = 0.0
    pd1_occupancy: float = 0.0


@dataclass(frozen=True)
class Dose:
    day: float
    amount: float


def _validate_model_inputs(parameters: TwinParameters, initial: InitialState) -> None:
    """Reject non-physical parameter and initial-state domains before integration."""
    parameter_values = asdict(parameters)
    invalid_parameters = [
        name for name, value in parameter_values.items()
        if not np.isfinite(value) or value < 0
    ]
    if invalid_parameters:
        raise ValueError(
            "systemic parameters must be finite and non-negative: "
            f"{sorted(invalid_parameters)}"
        )
    if parameters.tumour_capacity <= 0:
        raise ValueError("tumour_capacity must be strictly positive")

    initial_values = asdict(initial)
    invalid_initial = [
        name for name, value in initial_values.items()
        if not np.isfinite(value) or value < 0
    ]
    if invalid_initial:
        raise ValueError(
            "systemic initial states must be finite and non-negative: "
            f"{sorted(invalid_initial)}"
        )
    if initial.pd1_occupancy > 1:
        raise ValueError("initial pd1_occupancy must fall within [0, 1]")


def _rhs(_: float, y: np.ndarray, p: TwinParameters) -> np.ndarray:
    (cancer, antigen, dc, naive, tpex, eff_ln, eff_tum, tex, treg,
     ab_plasma, ab_ln, occupancy) = np.maximum(y, 0.0)
    suppression = 1.0 / (1.0 + p.treg_suppression * treg)
    checkpoint_relief = 0.25 + 0.75 * np.clip(occupancy, 0.0, 1.0)
    priming = p.priming_rate * dc * naive * suppression
    tpex_diff = p.tpex_to_effector * tpex * checkpoint_relief
    infiltration = p.tumour_infiltration * eff_tum
    killing = p.killing_rate * eff_tum * cancer * checkpoint_relief
    dydt = np.array([
        p.tumour_growth * cancer * (1.0 - cancer / p.tumour_capacity) - killing,
        p.antigen_release * cancer - p.antigen_clearance * antigen - p.dc_migration * antigen,
        p.dc_migration * antigen * 1.0e5 - p.dc_turnover * dc,
        -priming + 0.015 * (2.0e6 - naive),
        priming + p.tpex_self_renewal * tpex * suppression - tpex_diff - 0.05 * tpex,
        tpex_diff - p.effector_egress * eff_ln - p.effector_death * eff_ln,
        p.effector_egress * eff_ln + p.reinvigoration_rate * occupancy * tex
        - infiltration - p.exhaustion_rate * (1.0 - occupancy) * eff_tum - p.effector_death * eff_tum,
        p.exhaustion_rate * (1.0 - occupancy) * eff_tum
        - p.reinvigoration_rate * occupancy * tex - 0.04 * tex,
        0.03 * dc - 0.05 * treg,
        -p.antibody_clearance * ab_plasma - p.antibody_tdln_exchange * (ab_plasma - ab_ln),
        p.antibody_tdln_exchange * (ab_plasma - ab_ln) - p.antibody_tdln_loss * ab_ln,
        p.occupancy_on * ab_ln * (1.0 - occupancy) - p.occupancy_off * occupancy,
    ])
    # Numerical guard: a state at zero cannot have a negative derivative.
    return np.where((y <= 0.0) & (dydt < 0.0), 0.0, dydt)


def _validated_solution_states(values: np.ndarray, rtol: float, atol: float) -> np.ndarray:
    """Clip only solver-scale drift; reject material violations of state bounds."""
    states = np.asarray(values, dtype=float)
    if not np.isfinite(states).all():
        raise RuntimeError("systemic solver returned non-finite state values")
    scales = np.maximum(np.max(np.abs(states), axis=1), 1.0)
    tolerance = 10.0 * (atol + rtol * scales)
    minima = np.min(states, axis=1)
    invalid_negative = [
        STATE_NAMES[index] for index in np.flatnonzero(minima < -tolerance)
    ]
    occupancy_index = STATE_NAMES.index("pd1_occupancy")
    occupancy_max = float(np.max(states[occupancy_index]))
    invalid_occupancy = occupancy_max > 1.0 + tolerance[occupancy_index]
    if invalid_negative or invalid_occupancy:
        raise RuntimeError(
            "systemic solver violated physical state bounds: "
            f"negative={invalid_negative}, pd1_occupancy_max={occupancy_max}"
        )
    states = np.maximum(states, 0.0)
    states[occupancy_index] = np.minimum(states[occupancy_index], 1.0)
    return states


def simulate(
    days: float = 42.0,
    doses: tuple[Dose, ...] = (),
    parameters: TwinParameters | None = None,
    initial: InitialState | None = None,
    step: float = 0.25,
    rtol: float = 1e-7,
    atol: float = 1e-9,
) -> pd.DataFrame:
    """Simulate piecewise across bolus doses and return a tidy trajectory."""
    p, init = parameters or TwinParameters(), initial or InitialState()
    if not np.isfinite(days) or not np.isfinite(step) or days <= 0 or step <= 0:
        raise ValueError("days and step must be finite and positive")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol <= 0 or atol <= 0:
        raise ValueError("solver tolerances must be finite and positive")
    _validate_model_inputs(p, init)
    if any(not np.isfinite(d.day) or not np.isfinite(d.amount)
           or d.day < 0 or d.day >= days or d.amount < 0 for d in doses):
        raise ValueError(
            "doses must have non-negative amount and occur before the simulation endpoint"
        )
    events: dict[float, float] = {}
    for dose in doses:
        events[dose.day] = events.get(dose.day, 0.0) + dose.amount
    boundaries = sorted({0.0, days, *events})
    y = np.array([getattr(init, name) for name in STATE_NAMES], dtype=float)
    frames: list[pd.DataFrame] = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if start in events:
            y[STATE_NAMES.index("anti_pd1_plasma")] += events[start]
        times = np.arange(start, end, step)
        times = np.unique(np.append(times, end))
        sol = solve_ivp(_rhs, (start, end), y, args=(p,), t_eval=times, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        physical_states = _validated_solution_states(sol.y, rtol, atol)
        frame = pd.DataFrame(physical_states.T, columns=STATE_NAMES)
        frame.insert(0, "day", sol.t)
        if frames:
            if start in events:
                # Report a right-continuous trajectory: at a bolus boundary the
                # retained row is the post-dose state, not the pre-dose endpoint.
                frames[-1] = frames[-1].iloc[:-1]
            else:
                frame = frame.iloc[1:]
        frames.append(frame)
        y = physical_states[:, -1]
    out = pd.concat(frames, ignore_index=True)
    out.insert(1, "administered_dose", out.day.map(events).fillna(0.0))
    out["tumour_change_pct"] = 100.0 * (out.tumour_cells / out.tumour_cells.iloc[0] - 1.0)
    out["tpex_effector_ratio"] = out.tpex_tdln / out.effector_tdln.clip(lower=1.0)
    return out


def parameter_manifest(parameters: TwinParameters | None = None) -> dict[str, float]:
    return asdict(parameters or TwinParameters())
