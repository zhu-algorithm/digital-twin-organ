"""Global parameter sensitivity for the mechanistic TDLN model."""

from __future__ import annotations

from dataclasses import asdict, replace
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import qmc, spearmanr

from .simulator import Dose, TwinParameters, simulate

DEFAULT_PARAMETERS = (
    "tumour_growth", "dc_migration", "priming_rate", "tpex_self_renewal",
    "tpex_to_effector", "exhaustion_rate", "reinvigoration_rate",
    "killing_rate", "treg_suppression", "antibody_tdln_exchange",
)
OUTCOMES = ("tumour_cells", "effector_tumour", "tpex_tdln", "peak_pd1_occupancy")


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
