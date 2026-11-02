"""Patient/cohort data assimilation for the mechanistic TDLN twin."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .simulator import Dose, InitialState, STATE_NAMES, TwinParameters, simulate

DEFAULT_FIT_PARAMETERS = (
    "tumour_growth", "dc_migration", "priming_rate", "tpex_to_effector",
    "exhaustion_rate", "killing_rate", "antibody_tdln_exchange",
)


def validate_observations(observations: pd.DataFrame) -> pd.DataFrame:
    """Validate long-format observations: day, variable, value, optional sd."""
    required = {"day", "variable", "value"}
    missing = required - set(observations.columns)
    if missing:
        raise ValueError(f"missing observation columns: {sorted(missing)}")
    out = observations.copy()
    out["day"] = pd.to_numeric(out["day"], errors="raise")
    out["value"] = pd.to_numeric(out["value"], errors="raise")
    if "sd" not in out:
        out["sd"] = np.maximum(np.abs(out["value"]) * 0.1, 1.0)
    out["sd"] = pd.to_numeric(out["sd"], errors="raise")
    if (out.day < 0).any() or (out.sd <= 0).any() or not np.isfinite(out[["day", "value", "sd"]]).all().all():
        raise ValueError("days must be non-negative; values finite; sd positive")
    allowed = set(STATE_NAMES) | {"tumour_change_pct", "tpex_effector_ratio"}
    unknown = set(out.variable) - allowed
    if unknown:
        raise ValueError(f"unknown observed variables: {sorted(unknown)}")
    return out.sort_values(["day", "variable"]).reset_index(drop=True)


def _predictions(trajectory: pd.DataFrame, observations: pd.DataFrame) -> np.ndarray:
    return np.array([
        np.interp(row.day, trajectory.day, trajectory[row.variable])
        for row in observations.itertuples()
    ])


def calibrate(
    observations: pd.DataFrame,
    doses: tuple[Dose, ...] = (),
    base: TwinParameters | None = None,
    initial: InitialState | None = None,
    fit_parameters: Iterable[str] = DEFAULT_FIT_PARAMETERS,
    max_nfev: int = 200,
) -> tuple[TwinParameters, dict, pd.DataFrame]:
    """Fit positive parameters in log space using weighted robust residuals."""
    obs = validate_observations(observations)
    base = base or TwinParameters()
    names = tuple(fit_parameters)
    invalid = set(names) - set(asdict(base))
    if invalid or not names:
        raise ValueError(f"invalid fit parameters: {sorted(invalid)}")
    x0 = np.log([getattr(base, name) for name in names])
    days = max(float(obs.day.max()), 1e-6)

    def unpack(log_values: np.ndarray) -> TwinParameters:
        return replace(base, **dict(zip(names, np.exp(log_values))))

    def residual(log_values: np.ndarray) -> np.ndarray:
        trajectory = simulate(days=days, doses=doses, parameters=unpack(log_values), initial=initial)
        return (_predictions(trajectory, obs) - obs.value.to_numpy()) / obs.sd.to_numpy()

    result = least_squares(
        residual, x0, bounds=(x0 - np.log(100.0), x0 + np.log(100.0)),
        loss="soft_l1", max_nfev=max_nfev,
    )
    fitted = unpack(result.x)
    trajectory = simulate(days=days, doses=doses, parameters=fitted, initial=initial)
    prediction = _predictions(trajectory, obs)
    weighted = (prediction - obs.value.to_numpy()) / obs.sd.to_numpy()
    diagnostics = {
        "success": bool(result.success), "message": result.message,
        "n_observations": int(len(obs)), "n_parameters": len(names), "nfev": int(result.nfev),
        "weighted_rmse": float(np.sqrt(np.mean(weighted ** 2))),
        "aic_like": float(len(obs) * np.log(max(np.mean(weighted ** 2), 1e-15)) + 2 * len(names)),
        "fit_parameters": {name: float(getattr(fitted, name)) for name in names},
        "identifiability_warning": len(obs) < 5 * len(names),
        "research_only": True,
    }
    compared = obs.assign(prediction=prediction, weighted_residual=weighted)
    return fitted, diagnostics, compared


def save_calibration(output: str | Path, fitted: TwinParameters, diagnostics: dict, compared: pd.DataFrame) -> None:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "fitted_parameters.json").write_text(json.dumps(asdict(fitted), indent=2) + "\n")
    (output / "calibration_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    compared.to_csv(output / "calibration_residuals.csv", index=False)


def ensemble(
    parameters: TwinParameters,
    doses: tuple[Dose, ...],
    days: float,
    n: int = 100,
    relative_sd: float = 0.15,
    seed: int = 42,
) -> pd.DataFrame:
    """Propagate independent log-normal parameter uncertainty into outcome intervals."""
    if n < 2 or relative_sd <= 0:
        raise ValueError("n must be >=2 and relative_sd positive")
    rng = np.random.default_rng(seed)
    base = asdict(parameters)
    rows = []
    sigma = np.sqrt(np.log1p(relative_sd ** 2))
    for draw in range(n):
        sampled = {key: value * np.exp(rng.normal(-0.5 * sigma ** 2, sigma)) for key, value in base.items()}
        trajectory = simulate(days=days, doses=doses, parameters=TwinParameters(**sampled))
        final = trajectory.iloc[-1]
        rows.append({"draw": draw, "tumour_cells": final.tumour_cells,
                     "tumour_change_pct": final.tumour_change_pct,
                     "effector_tumour": final.effector_tumour,
                     "peak_pd1_occupancy": trajectory.pd1_occupancy.max()})
    return pd.DataFrame(rows)
