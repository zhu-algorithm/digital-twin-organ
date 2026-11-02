"""Zone-resolved TDLN transport and immune reaction model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .spatial import ZONES
from .simulator import Dose, TwinParameters, simulate

ZONE_NAMES = tuple(zone.name for zone in ZONES)
SPECIES = ("antigen", "mature_dc", "tpex", "effector")


@dataclass(frozen=True)
class SpatialParameters:
    antigen_to_sinus: float = 0.90
    antigen_sinus_to_paracortex: float = 0.12
    antigen_sinus_to_efferent: float = 0.18
    dc_to_sinus: float = 0.55
    dc_sinus_to_paracortex: float = 0.65
    tpex_to_medulla: float = 0.05
    tpex_medulla_to_efferent: float = 0.35
    tpex_to_follicle: float = 0.02
    tpex_follicle_to_paracortex: float = 0.04
    effector_to_medulla: float = 0.18
    effector_medulla_to_efferent: float = 0.65
    antigen_decay: float = 0.35
    antigen_uptake: float = 0.12
    dc_yield: float = 2.0e4
    dc_death: float = 0.18
    priming_rate: float = 1.0e-7
    tpex_growth: float = 0.10
    tpex_differentiation: float = 0.12
    checkpoint_growth_gain: float = 4.0
    checkpoint_differentiation_gain: float = 0.2
    tpex_death: float = 0.04
    effector_death: float = 0.12


@dataclass(frozen=True)
class SpatialInitialState:
    afferent_antigen: float = 1.0
    afferent_dc: float = 2.0e4
    paracortex_tpex: float = 8.0e4
    paracortex_effector: float = 2.0e4


# Directed first-order transport rates per day.  Movement is conservative:
# material removed from the source is added to the destination.
TRANSPORT = {
    "antigen": (
        ("afferent_lymph", "subcapsular_sinus", "antigen_to_sinus"),
        ("subcapsular_sinus", "paracortex", "antigen_sinus_to_paracortex"),
        ("subcapsular_sinus", "efferent_lymph", "antigen_sinus_to_efferent"),
    ),
    "mature_dc": (
        ("afferent_lymph", "subcapsular_sinus", "dc_to_sinus"),
        ("subcapsular_sinus", "paracortex", "dc_sinus_to_paracortex"),
    ),
    "tpex": (
        ("paracortex", "medulla", "tpex_to_medulla"),
        ("medulla", "efferent_lymph", "tpex_medulla_to_efferent"),
        ("paracortex", "follicle", "tpex_to_follicle"),
        ("follicle", "paracortex", "tpex_follicle_to_paracortex"),
    ),
    "effector": (
        ("paracortex", "medulla", "effector_to_medulla"),
        ("medulla", "efferent_lymph", "effector_medulla_to_efferent"),
    ),
}


def _index(species: str, zone: str) -> int:
    return SPECIES.index(species) * len(ZONE_NAMES) + ZONE_NAMES.index(zone)


def initial_vector(initial: SpatialInitialState | pd.DataFrame | None = None) -> np.ndarray:
    if isinstance(initial, pd.DataFrame):
        required = {"species", "zone", "amount"}
        missing = required.difference(initial.columns)
        if missing:
            raise ValueError(f"zone initial state missing columns: {sorted(missing)}")
        x = initial.loc[:, ["species", "zone", "amount"]].copy()
        if x[["species", "zone"]].isna().any().any() or x.duplicated(["species", "zone"]).any():
            raise ValueError("zone initial state needs unique, non-missing species-zone pairs")
        invalid_species = sorted(set(x.species) - set(SPECIES))
        invalid_zones = sorted(set(x.zone) - set(ZONE_NAMES))
        if invalid_species or invalid_zones:
            raise ValueError(f"invalid initial species={invalid_species}, zones={invalid_zones}")
        x["amount"] = pd.to_numeric(x.amount, errors="coerce")
        if not np.isfinite(x.amount).all() or x.amount.lt(0).any():
            raise ValueError("zone initial amounts must be finite and non-negative")
        y = np.zeros(len(SPECIES) * len(ZONE_NAMES), dtype=float)
        for row in x.itertuples(index=False):
            y[_index(row.species, row.zone)] = row.amount
        return y
    initial = initial or SpatialInitialState()
    y = np.zeros(len(SPECIES) * len(ZONE_NAMES), dtype=float)
    y[_index("antigen", "afferent_lymph")] = initial.afferent_antigen
    y[_index("mature_dc", "afferent_lymph")] = initial.afferent_dc
    y[_index("tpex", "paracortex")] = initial.paracortex_tpex
    y[_index("effector", "paracortex")] = initial.paracortex_effector
    return y


def _transport(y: np.ndarray, parameters: SpatialParameters) -> np.ndarray:
    dydt = np.zeros_like(y)
    for species, edges in TRANSPORT.items():
        for source, target, parameter_name in edges:
            rate = getattr(parameters, parameter_name)
            flux = rate * max(y[_index(species, source)], 0.0)
            dydt[_index(species, source)] -= flux
            dydt[_index(species, target)] += flux
    return dydt


def _rhs(
    time: float,
    y: np.ndarray,
    p: SpatialParameters,
    reactions: bool,
    checkpoint_occupancy: float | Callable[[float], float],
) -> np.ndarray:
    dydt = _transport(y, p)
    if reactions:
        # Antigen capture occurs in the subcapsular sinus; priming occurs in paracortex.
        scs_ag = max(y[_index("antigen", "subcapsular_sinus")], 0.0)
        uptake = p.antigen_uptake * scs_ag
        dydt[_index("antigen", "subcapsular_sinus")] -= uptake
        dydt[_index("mature_dc", "subcapsular_sinus")] += p.dc_yield * uptake
        for zone in ZONE_NAMES:
            dydt[_index("antigen", zone)] -= p.antigen_decay * max(y[_index("antigen", zone)], 0.0)
            dydt[_index("mature_dc", zone)] -= p.dc_death * max(y[_index("mature_dc", zone)], 0.0)
            dydt[_index("tpex", zone)] -= p.tpex_death * max(y[_index("tpex", zone)], 0.0)
            dydt[_index("effector", zone)] -= p.effector_death * max(y[_index("effector", zone)], 0.0)
        dc = max(y[_index("mature_dc", "paracortex")], 0.0)
        antigen = max(y[_index("antigen", "paracortex")], 0.0)
        tpex = max(y[_index("tpex", "paracortex")], 0.0)
        priming = p.priming_rate * dc * antigen
        occupancy = checkpoint_occupancy(time) if callable(checkpoint_occupancy) else checkpoint_occupancy
        occupancy = float(np.clip(occupancy, 0.0, 1.0))
        growth = (p.tpex_growth * tpex * antigen / (0.2 + antigen)
                  * (1.0 + p.checkpoint_growth_gain * occupancy))
        differentiation = (p.tpex_differentiation * tpex
                           * (1.0 + p.checkpoint_differentiation_gain * occupancy))
        dydt[_index("tpex", "paracortex")] += priming + growth - differentiation
        dydt[_index("effector", "paracortex")] += differentiation
    return np.where((y <= 0) & (dydt < 0), 0.0, dydt)


def simulate_spatial(
    days: float = 14.0,
    step: float = 0.1,
    parameters: SpatialParameters | None = None,
    initial: SpatialInitialState | pd.DataFrame | None = None,
    reactions: bool = True,
    checkpoint_occupancy: float | Callable[[float], float] = 0.0,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> pd.DataFrame:
    """Return tidy species-by-zone trajectories.

    ``checkpoint_occupancy`` is either a bounded constant or a time-varying
    callable, normally interpolated from the systemic tumour–TDLN model.  It
    modulates only reaction terms, so transport-only conservation remains an
    independently auditable invariant.
    """
    if not np.isfinite(days) or not np.isfinite(step) or days <= 0 or step <= 0:
        raise ValueError("days and step must be finite and positive")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol <= 0 or atol <= 0:
        raise ValueError("solver tolerances must be finite and positive")
    p = parameters or SpatialParameters()
    times = np.unique(np.append(np.arange(0.0, days, step), days))
    if not callable(checkpoint_occupancy) and not 0.0 <= checkpoint_occupancy <= 1.0:
        raise ValueError("checkpoint_occupancy must fall within [0, 1]")
    solution = solve_ivp(_rhs, (0.0, days), initial_vector(initial),
                         args=(p, reactions, checkpoint_occupancy),
                         t_eval=times, rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(solution.message)
    rows = []
    for time_index, day in enumerate(solution.t):
        for species in SPECIES:
            for zone in ZONE_NAMES:
                rows.append({"day": day, "species": species, "zone": zone,
                             "amount": max(solution.y[_index(species, zone), time_index], 0.0)})
    return pd.DataFrame(rows)


def simulate_spatial_coupled(
    system_trajectory: pd.DataFrame,
    days: float = 14.0,
    step: float = 0.1,
    parameters: SpatialParameters | None = None,
    initial: SpatialInitialState | pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Couple systemic PD-1 occupancy into the six-zone reaction network."""
    required = {"day", "pd1_occupancy"}
    missing = required.difference(system_trajectory.columns)
    if missing:
        raise ValueError(f"system trajectory missing columns: {sorted(missing)}")
    timeline = system_trajectory.loc[:, ["day", "pd1_occupancy"]].dropna().sort_values("day")
    if timeline.empty or timeline.day.duplicated().any():
        raise ValueError("system trajectory must contain unique, finite time points")
    if timeline.day.iloc[0] > 0 or timeline.day.iloc[-1] < days:
        raise ValueError("system trajectory must cover the requested spatial interval")
    if not timeline.pd1_occupancy.between(0.0, 1.0).all():
        raise ValueError("system PD-1 occupancy must fall within [0, 1]")

    occupancy = lambda time: float(np.interp(time, timeline.day, timeline.pd1_occupancy))
    return simulate_spatial(days, step, parameters, initial, True, occupancy)


def simulate_spatial_counterfactual(
    doses: tuple[Dose, ...],
    days: float = 14.0,
    step: float = 0.1,
    parameters: SpatialParameters | None = None,
    initial: SpatialInitialState | pd.DataFrame | None = None,
    systemic_parameters: TwinParameters | None = None,
) -> pd.DataFrame:
    """Run matched control and anti-PD-1 spatial arms from one initial state."""
    if not doses:
        raise ValueError("counterfactual simulation requires an explicit dose schedule")
    systemic_control = simulate(days=days, parameters=systemic_parameters)
    systemic_treated = simulate(days=days, doses=doses, parameters=systemic_parameters)
    frames = []
    for arm, systemic in (("control", systemic_control), ("anti_pd1", systemic_treated)):
        spatial = simulate_spatial_coupled(
            systemic, days=days, step=step, parameters=parameters, initial=initial,
        )
        spatial["pd1_occupancy"] = np.interp(
            spatial.day, systemic.day, systemic.pd1_occupancy
        )
        spatial["arm"] = arm
        frames.append(spatial)
    return pd.concat(frames, ignore_index=True)


def spatial_summary(trajectory: pd.DataFrame) -> pd.DataFrame:
    """Wide summary for visualization and invariant audits."""
    return trajectory.pivot_table(index="day", columns=["species", "zone"], values="amount").sort_index(axis=1)
