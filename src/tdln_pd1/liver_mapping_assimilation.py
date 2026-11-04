"""Calibrate label-free single-cell features to liver–TDLN bridge multipliers."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .assimilation import (_local_calibration_diagnostics,
                           calibration_parameter_draws)
from .hybrid import _outcome_like_columns
from .liver_bridge_features import (FeatureBridgeMapping, _validate_mapping,
                                    patient_liver_bridge_from_features)


TARGETS = {
    "antigen_drainage_multiplier",
    "dc_yield_multiplier",
    "immune_killing_multiplier",
}
FIT_PARAMETERS = {
    "drainage_dendritic_coefficient",
    "drainage_antigen_proxy_coefficient",
    "dc_yield_dendritic_coefficient",
    "killing_cytotoxic_coefficient",
    "killing_exhaustion_coefficient",
    "killing_cd8_fraction_coefficient",
}
FEATURE_COLUMNS = {
    "sig_dendritic_activation_dc", "sig_dendritic_activation",
    "sig_cytotoxic_cd8_t", "sig_cytotoxic",
    "sig_exhaustion_cd8_t", "sig_exhaustion",
    "frac_cd8_t", "frac_cd8_t_cell", "tumour_antigen_proxy",
}
AUDIT_MULTIPLIERS = {
    "antigen_drainage_multiplier": "antigen_drainage",
    "dc_yield_multiplier": "dc_yield",
    "immune_killing_multiplier": "immune_killing",
}


def validate_feature_bridge_observations(observations: pd.DataFrame) -> pd.DataFrame:
    """Validate long-form mechanistic targets without accepting outcome labels."""
    required = {"patient_id", "target", "value", "sd"}
    missing = required - set(observations)
    if missing:
        raise ValueError(f"missing feature-bridge observation columns: {sorted(missing)}")
    out = observations.copy()
    if out.empty:
        raise ValueError("feature-bridge observations must not be empty")
    out["patient_id"] = out.patient_id.astype(str).str.strip()
    out["target"] = out.target.astype(str).str.strip()
    if out.patient_id.eq("").any() or out.target.eq("").any():
        raise ValueError("feature-bridge patient IDs and targets must be non-empty")
    for column in ("value", "sd"):
        out[column] = pd.to_numeric(out[column], errors="raise")
    if not np.isfinite(out[["value", "sd"]]).all().all():
        raise ValueError("feature-bridge values and uncertainty must be finite")
    if out.value.le(0).any() or out.sd.le(0).any():
        raise ValueError("feature-bridge values and sd must be strictly positive")
    unknown = sorted(set(out.target) - TARGETS)
    if unknown:
        raise ValueError(f"unknown feature-bridge targets: {unknown}")
    if out.duplicated(["patient_id", "target"]).any():
        raise ValueError("feature-bridge observations need unique patient-target rows")
    forbidden = _outcome_like_columns(out[list(set(out) - {"target", "value", "sd"})])
    if forbidden:
        raise ValueError(f"feature-bridge calibration must be outcome-free: {forbidden}")
    feature_columns = sorted(FEATURE_COLUMNS & set(out))
    if not feature_columns:
        raise ValueError("feature-bridge observations contain no supported features")
    for column in feature_columns:
        out[column] = pd.to_numeric(out[column], errors="raise")
        if not np.isfinite(out[column]).all():
            raise ValueError(f"feature-bridge feature must be finite: {column}")
    for _, group in out.groupby("patient_id", sort=False):
        if any(group[column].nunique(dropna=False) != 1 for column in feature_columns):
            raise ValueError("feature values must be consistent within each patient")
    return out.sort_values(["patient_id", "target"]).reset_index(drop=True)


def _predictions(observations: pd.DataFrame,
                 mapping: FeatureBridgeMapping) -> np.ndarray:
    values = []
    for row in observations.itertuples(index=False):
        series = pd.Series(row._asdict())
        patient = patient_liver_bridge_from_features(series, mapping=mapping)
        values.append(patient.audit["multipliers"][AUDIT_MULTIPLIERS[row.target]])
    return np.asarray(values, dtype=float)


def calibrate_feature_bridge_mapping(
    observations: pd.DataFrame, *,
    fit_parameters: Iterable[str],
    base_mapping: FeatureBridgeMapping | None = None,
    max_nfev: int = 200,
) -> tuple[FeatureBridgeMapping, dict, pd.DataFrame]:
    """Fit selected positive mapping coefficients using weighted residuals."""
    obs = validate_feature_bridge_observations(observations)
    names = tuple(fit_parameters)
    invalid = sorted(set(names) - FIT_PARAMETERS)
    if not names or invalid:
        raise ValueError(f"invalid feature-bridge fit parameters: {invalid}")
    if len(set(names)) != len(names):
        raise ValueError("feature-bridge fit parameters must be unique")
    if not isinstance(max_nfev, (int, np.integer)) or max_nfev < 1:
        raise ValueError("max_nfev must be a positive integer")
    mapping = base_mapping or FeatureBridgeMapping()
    _validate_mapping(mapping)
    initial = np.asarray([getattr(mapping, name) for name in names], dtype=float)
    if not np.isfinite(initial).all() or np.any(initial <= 0):
        raise ValueError("fitted feature-bridge coefficients must start strictly positive")
    x0 = np.log(initial)
    lower, upper = x0 - np.log(100.0), x0 + np.log(100.0)

    def unpack(log_values: np.ndarray) -> FeatureBridgeMapping:
        return replace(mapping, **dict(zip(names, np.exp(log_values))))

    def residual(log_values: np.ndarray) -> np.ndarray:
        return ((_predictions(obs, unpack(log_values)) - obs.value.to_numpy())
                / obs.sd.to_numpy())

    result = least_squares(residual, x0, bounds=(lower, upper), loss="soft_l1",
                           diff_step=1e-3, max_nfev=max_nfev)
    fitted = unpack(result.x)
    prediction = _predictions(obs, fitted)
    weighted = (prediction - obs.value.to_numpy()) / obs.sd.to_numpy()
    local = _local_calibration_diagnostics(result, names, lower, upper, len(obs))
    diagnostics = {
        "success": bool(result.success), "message": result.message,
        "n_observations": int(len(obs)), "n_patients": int(obs.patient_id.nunique()),
        "n_parameters": len(names), "nfev": int(result.nfev),
        "weighted_rmse": float(np.sqrt(np.mean(weighted ** 2))),
        "fit_parameters": {name: float(getattr(fitted, name)) for name in names},
        **local,
        "identifiability_warning": bool(
            len(obs) < 5 * len(names) or not local["locally_identifiable"]),
        "labels_used": False, "research_only": True,
        "limitation": (
            "Local calibration to mechanistic bridge measurements; it is not "
            "calibration to PD-1 response labels or clinical validation."
        ),
    }
    compared = obs.assign(prediction=prediction, weighted_residual=weighted)
    return fitted, diagnostics, compared


def save_feature_bridge_mapping_calibration(
    output: str | Path, mapping: FeatureBridgeMapping,
    diagnostics: dict, compared: pd.DataFrame,
) -> None:
    target = Path(output)
    target.mkdir(parents=True, exist_ok=True)
    (target / "fitted_liver_feature_bridge_mapping.json").write_text(json.dumps({
        "mapping": asdict(mapping), "labels_used": False, "research_only": True,
    }, indent=2) + "\n")
    (target / "liver_feature_bridge_mapping_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n")
    compared.to_csv(target / "liver_feature_bridge_mapping_residuals.csv", index=False)


def load_feature_bridge_mapping(path: str | Path) -> FeatureBridgeMapping:
    payload = json.loads(Path(path).read_text())
    if payload.get("research_only") is not True or payload.get("labels_used") is not False:
        raise ValueError("feature-bridge mapping must be research-only and label-free")
    try:
        mapping = FeatureBridgeMapping(**payload["mapping"])
    except (KeyError, TypeError) as error:
        raise ValueError("invalid feature-bridge mapping schema") from error
    _validate_mapping(mapping)
    return mapping


def feature_bridge_mapping_draws(
    mapping: FeatureBridgeMapping, diagnostics: dict, n: int, seed: int,
) -> tuple[FeatureBridgeMapping, ...]:
    """Draw a reproducible correlated ensemble around a fitted mapping."""
    names, samples = calibration_parameter_draws(mapping, diagnostics, n, seed)
    if any(name not in FIT_PARAMETERS for name in names):
        raise ValueError("feature-bridge covariance contains a non-fittable parameter")
    draws = tuple(replace(mapping, **sample) for sample in samples)
    for draw in draws:
        _validate_mapping(draw)
    return draws
