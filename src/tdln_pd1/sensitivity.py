"""Global parameter sensitivity for the mechanistic TDLN model."""

from __future__ import annotations

from dataclasses import asdict, replace
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import qmc, spearmanr

from .simulator import Dose, TwinParameters, simulate
from .spatial_dynamics import (
    SpatialInitialState, SpatialParameters, simulate_spatial,
)

DEFAULT_PARAMETERS = (
    "tumour_growth", "dc_migration", "priming_rate", "tpex_self_renewal",
    "tpex_to_effector", "exhaustion_rate", "reinvigoration_rate",
    "killing_rate", "treg_suppression", "antibody_tdln_exchange",
)
OUTCOMES = ("tumour_cells", "effector_tumour", "tpex_tdln", "peak_pd1_occupancy")
DEFAULT_SPATIAL_PARAMETERS = (
    "antigen_to_sinus", "dc_sinus_to_paracortex", "antigen_decay",
    "antigen_uptake", "priming_rate", "antigen_half_saturation",
    "tpex_growth", "tpex_differentiation", "checkpoint_growth_gain",
    "checkpoint_differentiation_gain", "effector_medulla_to_efferent",
)
SPATIAL_OUTCOMES = (
    "control_efferent_effector", "occupied_efferent_effector",
    "pd1_efferent_effector_delta", "control_paracortex_tpex",
    "occupied_paracortex_tpex", "pd1_paracortex_tpex_delta",
)


def global_sensitivity(
    parameters: TwinParameters | None = None,
    names: Iterable[str] = DEFAULT_PARAMETERS,
    days: float = 42.0,
    doses: tuple[Dose, ...] = (Dose(0, 1), Dose(14, 1), Dose(28, 1)),
    n: int = 128,
    fold_range: float = 4.0,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Latin-hypercube sample log-scaled parameters and report rank correlations.

    Correlations are screening statistics, not variance-decomposition Sobol indices.
    They expose monotonic influence and help prioritize calibration measurements.
    """
    base = parameters or TwinParameters()
    names = tuple(names)
    invalid = set(names) - set(asdict(base))
    if invalid or not names:
        raise ValueError(f"invalid sensitivity parameters: {sorted(invalid)}")
    if n < 8 or fold_range <= 1 or days <= 0:
        raise ValueError("n must be >=8, fold_range >1, and days positive")
    unit = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    log_fold = np.log(fold_range)
    multipliers = np.exp(-log_fold + unit * 2 * log_fold)
    rows = []
    active_doses = tuple(dose for dose in doses if 0 <= dose.day <= days)
    for draw, vector in enumerate(multipliers):
        values = {name: getattr(base, name) * multiple for name, multiple in zip(names, vector)}
        trajectory = simulate(days=days, doses=active_doses, parameters=replace(base, **values))
        final = trajectory.iloc[-1]
        row = {"draw": draw, **values, "tumour_cells": final.tumour_cells,
               "effector_tumour": final.effector_tumour, "tpex_tdln": final.tpex_tdln,
               "peak_pd1_occupancy": trajectory.pd1_occupancy.max()}
        rows.append(row)
    samples = pd.DataFrame(rows)
    effects = []
    for name in names:
        for outcome in OUTCOMES:
            rho, p_value = spearmanr(samples[name], samples[outcome])
            effects.append({"parameter": name, "outcome": outcome,
                            "spearman_rho": float(rho), "p_value": float(p_value),
                            "abs_rho": float(abs(rho)), "n": n})
    ranking = pd.DataFrame(effects).sort_values(["outcome", "abs_rho"], ascending=[True, False])
    return samples, ranking.reset_index(drop=True)


def global_spatial_sensitivity(
    parameters: SpatialParameters | None = None,
    initial: SpatialInitialState | pd.DataFrame | None = None,
    names: Iterable[str] = DEFAULT_SPATIAL_PARAMETERS,
    days: float = 14.0,
    checkpoint_occupancy: float = 0.8,
    n: int = 128,
    fold_range: float = 4.0,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Screen six-zone parameters against baseline and PD-1 occupancy effects.

    Every Latin-hypercube draw is evaluated twice from the same initial state:
    once without checkpoint occupancy and once at the prespecified occupancy.
    This isolates the model-implied spatial PD-1 contrast from parameter effects
    on the untreated trajectory. Correlations remain screening statistics.
    """
    base = parameters or SpatialParameters()
    names = tuple(names)
    schema = asdict(base)
    invalid = set(names) - set(schema)
    nonpositive = [name for name in names if name in schema and schema[name] <= 0]
    if invalid or not names:
        raise ValueError(f"invalid spatial sensitivity parameters: {sorted(invalid)}")
    if nonpositive:
        raise ValueError(
            "log-scaled spatial sensitivity parameters must be strictly positive: "
            f"{sorted(nonpositive)}"
        )
    if (n < 8 or not np.isfinite(fold_range) or fold_range <= 1
            or not np.isfinite(days) or days <= 0):
        raise ValueError("n must be >=8, fold_range finite and >1, and days finite and positive")
    if (not np.isfinite(checkpoint_occupancy)
            or not 0 < checkpoint_occupancy <= 1):
        raise ValueError("checkpoint_occupancy must be finite and fall within (0, 1]")

    unit = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    log_fold = np.log(fold_range)
    multipliers = np.exp(-log_fold + unit * 2 * log_fold)
    rows = []
    for draw, vector in enumerate(multipliers):
        values = {name: getattr(base, name) * multiple for name, multiple in zip(names, vector)}
        sampled = replace(base, **values)
        control = simulate_spatial(
            days=days, parameters=sampled, initial=initial, checkpoint_occupancy=0.0,
        )
        occupied = simulate_spatial(
            days=days, parameters=sampled, initial=initial,
            checkpoint_occupancy=checkpoint_occupancy,
        )

        def endpoint(frame: pd.DataFrame, species: str, zone: str) -> float:
            selected = frame[(frame.species == species) & (frame.zone == zone)]
            return float(selected.amount.iloc[-1])

        control_effector = endpoint(control, "effector", "efferent_lymph")
        occupied_effector = endpoint(occupied, "effector", "efferent_lymph")
        control_tpex = endpoint(control, "tpex", "paracortex")
        occupied_tpex = endpoint(occupied, "tpex", "paracortex")
        rows.append({
            "draw": draw, **values,
            "control_efferent_effector": control_effector,
            "occupied_efferent_effector": occupied_effector,
            "pd1_efferent_effector_delta": occupied_effector - control_effector,
            "control_paracortex_tpex": control_tpex,
            "occupied_paracortex_tpex": occupied_tpex,
            "pd1_paracortex_tpex_delta": occupied_tpex - control_tpex,
        })
    samples = pd.DataFrame(rows)
    effects = []
    for name in names:
        for outcome in SPATIAL_OUTCOMES:
            rho, p_value = spearmanr(samples[name], samples[outcome])
            effects.append({
                "parameter": name, "outcome": outcome,
                "spearman_rho": float(rho), "p_value": float(p_value),
                "abs_rho": float(abs(rho)), "n": n,
                "checkpoint_occupancy": float(checkpoint_occupancy),
            })
    ranking = pd.DataFrame(effects).sort_values(
        ["outcome", "abs_rho"], ascending=[True, False]
    )
    return samples, ranking.reset_index(drop=True)
