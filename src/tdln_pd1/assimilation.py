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


def _local_calibration_diagnostics(result, names, lower, upper, n_observations) -> dict:
    """Return local Jacobian, boundary and approximate log-scale uncertainty diagnostics."""
    jacobian = np.asarray(result.jac, dtype=float)
    singular = np.linalg.svd(jacobian, compute_uv=False) if jacobian.size else np.array([])
    rank = int(np.linalg.matrix_rank(jacobian)) if jacobian.size else 0
    condition = (float(singular[0] / singular[-1])
                 if len(singular) and singular[-1] > 0 else float("inf"))
    span = np.asarray(upper) - np.asarray(lower)
    tolerance = 1e-6 * np.maximum(span, 1.0)
    boundary = [
        name for name, value, low, high, tol in
        zip(names, result.x, lower, upper, tolerance)
        if value - low <= tol or high - value <= tol
    ]
    uncertainty = None
    dof = int(n_observations - len(names))
    if rank == len(names) and dof > 0:
        residual_variance = float(np.sum(np.asarray(result.fun) ** 2) / dof)
        covariance = np.linalg.pinv(jacobian.T @ jacobian) * residual_variance
        standard_error = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        uncertainty = {
            name: {
                "log_scale_standard_error": float(error),
                "approximate_lower_95": float(np.exp(value - 1.96 * error)),
                "approximate_upper_95": float(np.exp(value + 1.96 * error)),
            }
            for name, value, error in zip(names, result.x, standard_error)
        }
    return {
        "jacobian_rank": rank,
        "jacobian_condition_number": condition if np.isfinite(condition) else None,
        "parameters_at_search_boundary": boundary,
        "local_parameter_uncertainty": uncertainty,
        "uncertainty_note": (
            "Approximate local intervals from the weighted residual Jacobian in log-parameter "
            "space; unavailable when rank deficient or residual degrees of freedom are absent."
        ),
        "locally_identifiable": bool(
            rank == len(names) and np.isfinite(condition) and condition <= 1e8 and not boundary
        ),
    }


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

    lower, upper = x0 - np.log(100.0), x0 + np.log(100.0)
    result = least_squares(
        residual, x0, bounds=(lower, upper),
        loss="soft_l1", max_nfev=max_nfev,
    )
    fitted = unpack(result.x)
    trajectory = simulate(days=days, doses=doses, parameters=fitted, initial=initial)
    prediction = _predictions(trajectory, obs)
    weighted = (prediction - obs.value.to_numpy()) / obs.sd.to_numpy()
    local = _local_calibration_diagnostics(result, names, lower, upper, len(obs))
    diagnostics = {
        "success": bool(result.success), "message": result.message,
        "n_observations": int(len(obs)), "n_parameters": len(names), "nfev": int(result.nfev),
        "weighted_rmse": float(np.sqrt(np.mean(weighted ** 2))),
        "aic_like": float(len(obs) * np.log(max(np.mean(weighted ** 2), 1e-15)) + 2 * len(names)),
        "fit_parameters": {name: float(getattr(fitted, name)) for name in names},
        **local,
        "identifiability_warning": bool(
            len(obs) < 5 * len(names) or not local["locally_identifiable"]
        ),
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
