"""Zone-resolved TDLN transport and immune reaction model."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .spatial import ZONES

ZONE_NAMES = tuple(zone.name for zone in ZONES)
SPECIES = ("antigen", "mature_dc", "tpex", "effector")


@dataclass(frozen=True)
class SpatialParameters:
    antigen_decay: float = 0.35
    antigen_uptake: float = 0.12
    dc_yield: float = 2.0e4
    dc_death: float = 0.18
    priming_rate: float = 1.0e-7
    tpex_growth: float = 0.10
    tpex_differentiation: float = 0.12
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
        ("afferent_lymph", "subcapsular_sinus", 0.90),
        ("subcapsular_sinus", "paracortex", 0.12),
        ("subcapsular_sinus", "efferent_lymph", 0.18),
    ),
    "mature_dc": (
        ("afferent_lymph", "subcapsular_sinus", 0.55),
        ("subcapsular_sinus", "paracortex", 0.65),
    ),
    "tpex": (
        ("paracortex", "medulla", 0.05),
        ("medulla", "efferent_lymph", 0.35),
        ("paracortex", "follicle", 0.02),
        ("follicle", "paracortex", 0.04),
    ),
    "effector": (
        ("paracortex", "medulla", 0.18),
        ("medulla", "efferent_lymph", 0.65),
    ),
}


def _index(species: str, zone: str) -> int:
    return SPECIES.index(species) * len(ZONE_NAMES) + ZONE_NAMES.index(zone)


def initial_vector(initial: SpatialInitialState | None = None) -> np.ndarray:
    initial = initial or SpatialInitialState()
    y = np.zeros(len(SPECIES) * len(ZONE_NAMES), dtype=float)
    y[_index("antigen", "afferent_lymph")] = initial.afferent_antigen
    y[_index("mature_dc", "afferent_lymph")] = initial.afferent_dc
    y[_index("tpex", "paracortex")] = initial.paracortex_tpex
    y[_index("effector", "paracortex")] = initial.paracortex_effector
    return y


def _transport(y: np.ndarray) -> np.ndarray:
    dydt = np.zeros_like(y)
    for species, edges in TRANSPORT.items():
        for source, target, rate in edges:
            flux = rate * max(y[_index(species, source)], 0.0)
            dydt[_index(species, source)] -= flux
            dydt[_index(species, target)] += flux
    return dydt


def _rhs(_: float, y: np.ndarray, p: SpatialParameters, reactions: bool) -> np.ndarray:
    dydt = _transport(y)
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
        growth = p.tpex_growth * tpex * antigen / (0.2 + antigen)
        differentiation = p.tpex_differentiation * tpex
        dydt[_index("tpex", "paracortex")] += priming + growth - differentiation
        dydt[_index("effector", "paracortex")] += differentiation
    return np.where((y <= 0) & (dydt < 0), 0.0, dydt)


def simulate_spatial(
    days: float = 14.0,
    step: float = 0.1,
    parameters: SpatialParameters | None = None,
    initial: SpatialInitialState | None = None,
    reactions: bool = True,
) -> pd.DataFrame:
    """Return tidy species-by-zone trajectories."""
    if days <= 0 or step <= 0:
        raise ValueError("days and step must be positive")
    p = parameters or SpatialParameters()
    times = np.unique(np.append(np.arange(0.0, days, step), days))
    solution = solve_ivp(_rhs, (0.0, days), initial_vector(initial), args=(p, reactions),
                         t_eval=times, rtol=1e-8, atol=1e-10)
    if not solution.success:
        raise RuntimeError(solution.message)
    rows = []
    for time_index, day in enumerate(solution.t):
        for species in SPECIES:
            for zone in ZONE_NAMES:
                rows.append({"day": day, "species": species, "zone": zone,
                             "amount": max(solution.y[_index(species, zone), time_index], 0.0)})
    return pd.DataFrame(rows)


def spatial_summary(trajectory: pd.DataFrame) -> pd.DataFrame:
    """Wide summary for visualization and invariant audits."""
    return trajectory.pivot_table(index="day", columns=["species", "zone"], values="amount").sort_index(axis=1)
