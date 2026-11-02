"""Longitudinal calibration for the zone-resolved TDLN model."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .spatial_dynamics import (SPECIES, ZONE_NAMES, SpatialInitialState,
                               SpatialParameters, simulate_spatial,
                               simulate_spatial_coupled)


DEFAULT_SPATIAL_FIT_PARAMETERS = (
    "antigen_to_sinus", "dc_sinus_to_paracortex", "effector_to_medulla",
    "antigen_decay", "antigen_uptake", "dc_death", "priming_rate",
    "tpex_growth", "tpex_differentiation", "effector_death",
)


def validate_spatial_observations(observations: pd.DataFrame) -> pd.DataFrame:
    required = {"day", "species", "zone", "amount"}
    missing = required.difference(observations.columns)
    if missing:
        raise ValueError(f"missing spatial observation columns: {sorted(missing)}")
    out = observations.copy()
    out["day"] = pd.to_numeric(out.day, errors="raise")
    out["amount"] = pd.to_numeric(out.amount, errors="raise")
    if "sd" not in out:
        out["sd"] = np.maximum(np.abs(out.amount) * 0.1, 1e-6)
    out["sd"] = pd.to_numeric(out.sd, errors="raise")
    if (out.day < 0).any() or (out.amount < 0).any() or (out.sd <= 0).any():
        raise ValueError("spatial days and amounts must be non-negative; sd must be positive")
    if not np.isfinite(out[["day", "amount", "sd"]]).all().all():
        raise ValueError("spatial observations must be finite")
    invalid_species = sorted(set(out.species) - set(SPECIES))
    invalid_zones = sorted(set(out.zone) - set(ZONE_NAMES))
    if invalid_species or invalid_zones:
        raise ValueError(
            f"unknown spatial species={invalid_species}, zones={invalid_zones}"
        )
    if out.duplicated(["day", "species", "zone"]).any():
        raise ValueError("spatial observations need unique day-species-zone rows")
    return out.sort_values(["day", "species", "zone"]).reset_index(drop=True)


def _spatial_predictions(trajectory: pd.DataFrame, observations: pd.DataFrame) -> np.ndarray:
    values = []
    for row in observations.itertuples(index=False):
        series = trajectory[
            trajectory.species.eq(row.species) & trajectory.zone.eq(row.zone)
        ]
        values.append(np.interp(row.day, series.day, series.amount))
    return np.asarray(values)


def calibrate_spatial(
    observations: pd.DataFrame,
    base: SpatialParameters | None = None,
    initial: SpatialInitialState | pd.DataFrame | None = None,
    systemic_trajectory: pd.DataFrame | None = None,
    fit_parameters: Iterable[str] = DEFAULT_SPATIAL_FIT_PARAMETERS,
    max_nfev: int = 200,
) -> tuple[SpatialParameters, dict, pd.DataFrame]:
    """Fit selected positive spatial parameters using robust weighted residuals."""
    obs = validate_spatial_observations(observations)
    base = base or SpatialParameters()
    names = tuple(fit_parameters)
    invalid = set(names) - set(asdict(base))
    if invalid or not names:
        raise ValueError(f"invalid spatial fit parameters: {sorted(invalid)}")
    x0 = np.log([getattr(base, name) for name in names])
    days = max(float(obs.day.max()), 1e-6)

    def unpack(log_values: np.ndarray) -> SpatialParameters:
        return replace(base, **dict(zip(names, np.exp(log_values))))

    def run(parameters: SpatialParameters) -> pd.DataFrame:
        if systemic_trajectory is None:
            return simulate_spatial(days=days, parameters=parameters, initial=initial)
        return simulate_spatial_coupled(
            systemic_trajectory, days=days, parameters=parameters, initial=initial,
        )

    def residual(log_values: np.ndarray) -> np.ndarray:
        prediction = _spatial_predictions(run(unpack(log_values)), obs)
        return (prediction - obs.amount.to_numpy()) / obs.sd.to_numpy()

    result = least_squares(
        residual, x0, bounds=(x0 - np.log(100.0), x0 + np.log(100.0)),
        loss="soft_l1", max_nfev=max_nfev,
    )
    fitted = unpack(result.x)
    prediction = _spatial_predictions(run(fitted), obs)
    weighted = (prediction - obs.amount.to_numpy()) / obs.sd.to_numpy()
    singular = np.linalg.svd(result.jac, compute_uv=False) if result.jac.size else np.array([])
    rank = int(np.linalg.matrix_rank(result.jac)) if result.jac.size else 0
    condition = (float(singular[0] / singular[-1])
                 if len(singular) and singular[-1] > 0 else float("inf"))
    warning = len(obs) < 5 * len(names) or rank < len(names) or condition > 1e8
    diagnostics = {
        "success": bool(result.success), "message": result.message,
        "n_observations": int(len(obs)), "n_parameters": len(names),
        "nfev": int(result.nfev),
        "weighted_rmse": float(np.sqrt(np.mean(weighted ** 2))),
        "aic_like": float(len(obs) * np.log(max(np.mean(weighted ** 2), 1e-15))
                          + 2 * len(names)),
        "fit_parameters": {name: float(getattr(fitted, name)) for name in names},
        "jacobian_rank": rank,
        "jacobian_condition_number": condition if np.isfinite(condition) else None,
        "identifiability_warning": bool(warning), "research_only": True,
    }
    compared = obs.assign(prediction=prediction, weighted_residual=weighted)
    return fitted, diagnostics, compared


def save_spatial_calibration(
    output: str | Path,
    fitted: SpatialParameters,
    diagnostics: dict,
    compared: pd.DataFrame,
) -> None:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "fitted_spatial_parameters.json").write_text(
        json.dumps(asdict(fitted), indent=2) + "\n"
    )
    (output / "spatial_calibration_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )
    compared.to_csv(output / "spatial_calibration_residuals.csv", index=False)
