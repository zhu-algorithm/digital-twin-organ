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
    covariance_payload = None
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
        covariance_payload = {
            "parameter_order": list(names),
            "log_scale_covariance": covariance.tolist(),
        }
    return {
        "jacobian_rank": rank,
        "jacobian_condition_number": condition if np.isfinite(condition) else None,
        "parameters_at_search_boundary": boundary,
        "local_parameter_uncertainty": uncertainty,
        "local_log_parameter_covariance": covariance_payload,
        "uncertainty_note": (
            "Approximate local intervals from the weighted residual Jacobian in log-parameter "
            "space; unavailable when rank deficient or residual degrees of freedom are absent."
        ),
        "locally_identifiable": bool(
            rank == len(names) and np.isfinite(condition) and condition <= 1e8 and not boundary
        ),
    }


def _profile_likelihood(
    residual,
    result,
    names,
    lower,
    upper,
    points: int,
    fold_range: float,
    max_nfev: int,
) -> dict | None:
    """Profile each fitted parameter while re-optimizing all remaining parameters."""
    if points == 0:
        return None
    if points < 3 or points % 2 == 0:
        raise ValueError("profile_points must be zero or an odd integer >= 3")
    if not np.isfinite(fold_range) or fold_range <= 1:
        raise ValueError("profile_fold_range must be finite and > 1")
    names = tuple(names)
    lower, upper = np.asarray(lower), np.asarray(upper)
    gaussian_fit = least_squares(
        residual, np.asarray(result.x, dtype=float), bounds=(lower, upper),
        loss="linear", max_nfev=max_nfev,
    )
    optimum = np.asarray(gaussian_fit.x, dtype=float)
    baseline_sse = float(np.sum(np.asarray(residual(optimum)) ** 2))
    half = points // 2
    radius = np.log(fold_range)
    profiles = {}
    for fixed_index, name in enumerate(names):
        left = np.linspace(max(lower[fixed_index], optimum[fixed_index] - radius),
                           optimum[fixed_index], half + 1)
        right = np.linspace(optimum[fixed_index],
                            min(upper[fixed_index], optimum[fixed_index] + radius),
                            half + 1)[1:]
        grid = np.concatenate([left, right])
        rows = []
        free = np.array([index for index in range(len(names)) if index != fixed_index])
        for fixed_log_value in grid:
            if len(free):
                def profiled_residual(free_values):
                    candidate = optimum.copy()
                    candidate[fixed_index] = fixed_log_value
                    candidate[free] = free_values
                    return residual(candidate)

                refit = least_squares(
                    profiled_residual, optimum[free],
                    bounds=(lower[free], upper[free]), loss="linear",
                    max_nfev=max_nfev,
                )
                candidate = optimum.copy()
                candidate[fixed_index] = fixed_log_value
                candidate[free] = refit.x
                refit_success = bool(refit.success)
            else:
                candidate = np.array([fixed_log_value])
                refit_success = True
            sse = float(np.sum(np.asarray(residual(candidate)) ** 2))
            rows.append({
                "parameter_value": float(np.exp(fixed_log_value)),
                "weighted_sse": sse,
                "delta_weighted_sse": float(max(0.0, sse - baseline_sse)),
                "nuisance_refit_success": refit_success,
            })
        supported = [row for row in rows if row["delta_weighted_sse"] <= 3.841]
        profiles[name] = {
            "points": rows,
            "scanned_support_lower": min(row["parameter_value"] for row in supported),
            "scanned_support_upper": max(row["parameter_value"] for row in supported),
            "lower_limit_resolved": rows[0]["delta_weighted_sse"] > 3.841,
            "upper_limit_resolved": rows[-1]["delta_weighted_sse"] > 3.841,
        }
    return {
        "method": "one_parameter_profile_with_nuisance_refit",
        "threshold_delta_weighted_sse": 3.841,
        "threshold_reference": "asymptotic chi-square, 1 degree of freedom",
        "points_per_parameter": points,
        "fold_range": float(fold_range),
        "baseline_weighted_sse": baseline_sse,
        "gaussian_baseline_refit_success": bool(gaussian_fit.success),
        "gaussian_profile_optimum": {
            name: float(np.exp(value)) for name, value in zip(names, optimum)
        },
        "profiles": profiles,
        "note": (
            "Scanned support is not a resolved confidence limit when the corresponding "
            "limit_resolved flag is false; asymptotic coverage is not guaranteed."
        ),
    }


def validate_observations(observations: pd.DataFrame) -> pd.DataFrame:
    """Validate long-format observations: day, variable, value, optional sd."""
    required = {"day", "variable", "value"}
    missing = required - set(observations.columns)
    if missing:
        raise ValueError(f"missing observation columns: {sorted(missing)}")
    out = observations.copy()
    if out.empty:
        raise ValueError("calibration observations must not be empty")
    if out["variable"].isna().any():
        raise ValueError("observed variable names must be non-missing and non-empty")
    out["variable"] = out["variable"].astype(str).str.strip()
    if out["variable"].eq("").any():
        raise ValueError("observed variable names must be non-missing and non-empty")
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
    if out.duplicated(["day", "variable"]).any():
        raise ValueError("calibration observations need unique day-variable rows")
    return out.sort_values(["day", "variable"]).reset_index(drop=True)


def _validate_fit_parameters(base, fit_parameters: Iterable[str], label: str) -> tuple[str, ...]:
    """Validate a log-space parameter selection against its dataclass instance."""
    names = tuple(fit_parameters)
    available = asdict(base)
    invalid = set(names) - set(available)
    if invalid or not names:
        raise ValueError(f"invalid {label}fit parameters: {sorted(invalid)}")
    if len(set(names)) != len(names):
        raise ValueError(f"{label}fit parameters must be unique")
    nonpositive = [
        name for name in names
        if not np.isfinite(available[name]) or available[name] <= 0
    ]
    if nonpositive:
        raise ValueError(
            f"{label}fit parameters need finite, strictly positive base values for log-space fitting: "
            f"{sorted(nonpositive)}"
        )
    return names


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
    profile_points: int = 0,
    profile_fold_range: float = 4.0,
) -> tuple[TwinParameters, dict, pd.DataFrame]:
    """Fit positive parameters in log space using weighted robust residuals."""
    obs = validate_observations(observations)
    base = base or TwinParameters()
    names = _validate_fit_parameters(base, fit_parameters, "")
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
    profile = _profile_likelihood(
        residual, result, names, lower, upper, profile_points,
        profile_fold_range, max_nfev,
    )
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
        "profile_likelihood": profile,
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


def calibration_ensemble(
    parameters: TwinParameters,
    diagnostics: dict,
    doses: tuple[Dose, ...],
    days: float,
    n: int = 100,
    seed: int = 42,
) -> pd.DataFrame:
    """Propagate the fitted local log-parameter covariance into endpoints.

    Only fitted coordinates in the saved covariance are sampled; all other
    mechanistic parameters remain fixed at their calibrated-file values.
    This is local asymptotic parameter uncertainty, not structural uncertainty.
    """
    if not np.isfinite(days) or days <= 0:
        raise ValueError("calibration ensemble days must be finite and positive")
    names, draws = calibration_parameter_draws(parameters, diagnostics, n, seed)
    rows = []
    for draw, sampled_values in enumerate(draws):
        sampled_parameters = replace(parameters, **sampled_values)
        trajectory = simulate(days=days, doses=doses, parameters=sampled_parameters)
        final = trajectory.iloc[-1]
        rows.append({
            "draw": draw,
            **{f"parameter_{name}": sampled_values[name] for name in names},
            "tumour_cells": float(final.tumour_cells),
            "tumour_change_pct": float(final.tumour_change_pct),
            "effector_tumour": float(final.effector_tumour),
            "tpex_tdln": float(final.tpex_tdln),
            "peak_pd1_occupancy": float(trajectory.pd1_occupancy.max()),
        })
    return pd.DataFrame(rows)


def calibration_parameter_draws(parameters, diagnostics: dict, n: int, seed: int):
    """Draw correlated positive dataclass parameters from saved log covariance."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("calibration ensemble requires n >= 2")
    payload = diagnostics.get("local_log_parameter_covariance")
    if not isinstance(payload, dict):
        raise ValueError("calibration diagnostics do not contain estimable covariance")
    if diagnostics.get("locally_identifiable") is not True:
        raise ValueError("calibration uncertainty requires locally identifiable parameters")
    names = payload.get("parameter_order")
    covariance = np.asarray(payload.get("log_scale_covariance"), dtype=float)
    available = asdict(parameters)
    if (not isinstance(names, list) or not names or len(set(names)) != len(names)
            or any(name not in available for name in names)
            or covariance.shape != (len(names), len(names))
            or not np.isfinite(covariance).all()
            or not np.allclose(covariance, covariance.T, rtol=1e-8, atol=1e-12)):
        raise ValueError("invalid calibration log-parameter covariance")
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues.min() < -1e-10 * max(1.0, float(np.abs(eigenvalues).max())):
        raise ValueError("calibration covariance must be positive semidefinite")
    rng = np.random.default_rng(seed)
    log_center = np.log([available[name] for name in names])
    vectors = rng.multivariate_normal(log_center, covariance, size=n, check_valid="raise")
    return names, [dict(zip(names, np.exp(vector))) for vector in vectors]
