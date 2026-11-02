"""Calibrate the measured-feature to mechanistic-state bridge."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .assimilation import _local_calibration_diagnostics, _validate_fit_parameters
from .hybrid import (
    FeatureStateMapping, initial_state_from_features,
    spatial_initial_state_from_features,
)


MAPPING_TARGETS = (
    "mature_dc_tdln", "tpex_tdln", "effector_tumour", "exhausted_tumour",
    "treg_tdln", "afferent_antigen", "afferent_dc", "paracortex_tpex",
    "paracortex_effector",
)
MAPPING_FEATURES = (
    "sig_cytotoxic_cd8_t", "sig_cytotoxic", "sig_exhaustion_cd8_t",
    "sig_exhaustion", "sig_tcf7_memory_cd8_t", "sig_tcf7_memory",
    "sig_dendritic_activation_dc", "sig_dendritic_activation", "frac_treg",
    "frac_regulatory_t_cell", "frac_cd8_t", "frac_cd8_t_cell",
    "tumour_antigen_proxy",
)
CALIBRATABLE_MAPPING_PARAMETERS = (
    "signature_intercept", "cytotoxic_coefficient", "cd8_fraction_coefficient",
    "exhaustion_coefficient", "tcf7_coefficient", "dendritic_coefficient",
    "treg_fraction_coefficient", "antigen_proxy_coefficient",
)


def validate_mapping_observations(observations: pd.DataFrame) -> pd.DataFrame:
    """Validate long-form paired assay/state measurements without using outcomes."""
    required = {"patient_id", "target", "value", "sd"}
    missing = required - set(observations.columns)
    if missing:
        raise ValueError(f"mapping observations missing columns: {sorted(missing)}")
    out = observations.copy()
    unknown_columns = sorted(set(out.columns) - required - set(MAPPING_FEATURES))
    if unknown_columns:
        raise ValueError(f"unknown mapping observation columns: {unknown_columns}")
    if out.empty or out[list(required)].isna().any().any():
        raise ValueError("mapping observations require complete patient, target, value and sd")
    out["patient_id"] = out.patient_id.astype(str).str.strip()
    out["target"] = out.target.astype(str).str.strip()
    if out.patient_id.eq("").any() or out.target.eq("").any():
        raise ValueError("mapping observations require non-empty patient and target")
    unknown = sorted(set(out.target) - set(MAPPING_TARGETS))
    if unknown:
        raise ValueError(f"unknown mapping targets: {unknown}")
    if out.duplicated(["patient_id", "target"]).any():
        raise ValueError("mapping observations need unique patient-target rows")
    feature_columns = [name for name in MAPPING_FEATURES if name in out]
    for name in feature_columns:
        out[name] = pd.to_numeric(out[name], errors="raise")
        if out[name].notna().any() and not np.isfinite(out.loc[out[name].notna(), name]).all():
            raise ValueError(f"mapping feature must be finite when supplied: {name}")
    if feature_columns:
        inconsistent = [
            name for name in feature_columns
            if out.groupby("patient_id")[name].nunique(dropna=False).gt(1).any()
        ]
        if inconsistent:
            raise ValueError(
                "mapping features must be consistent across targets for each patient: "
                f"{sorted(inconsistent)}"
            )
    for name in ("value", "sd"):
        out[name] = pd.to_numeric(out[name], errors="raise")
    if (not np.isfinite(out[["value", "sd"]]).all().all()
            or out.value.lt(0).any() or out.sd.le(0).any()):
        raise ValueError("mapping values must be finite/non-negative and sd finite/positive")
    return out.sort_values(["patient_id", "target"]).reset_index(drop=True)


def _predict_row(row: pd.Series, mapping: FeatureStateMapping) -> float:
    systemic = initial_state_from_features(row, mapping=mapping)
    spatial = spatial_initial_state_from_features(row, mapping=mapping)
    target = str(row.target)
    if hasattr(systemic, target):
        return float(getattr(systemic, target))
    return float(getattr(spatial, target))


def calibrate_feature_state_mapping(
    observations: pd.DataFrame,
    base: FeatureStateMapping | None = None,
    fit_parameters: Iterable[str] = ("tcf7_coefficient",),
    max_nfev: int = 200,
) -> tuple[FeatureStateMapping, dict, pd.DataFrame]:
    """Fit positive mapping coefficients against quantitative immune states."""
    obs = validate_mapping_observations(observations)
    base = base or FeatureStateMapping()
    names = _validate_fit_parameters(base, fit_parameters, "mapping ")
    unsupported = sorted(set(names) - set(CALIBRATABLE_MAPPING_PARAMETERS))
    if unsupported:
        raise ValueError(f"unsupported mapping calibration parameters: {unsupported}")
    x0 = np.log([getattr(base, name) for name in names])

    def unpack(log_values: np.ndarray) -> FeatureStateMapping:
        return replace(base, **dict(zip(names, np.exp(log_values))))

    def predictions(log_values: np.ndarray) -> np.ndarray:
        mapping = unpack(log_values)
        return np.asarray([_predict_row(row, mapping) for _, row in obs.iterrows()])

    def residual(log_values: np.ndarray) -> np.ndarray:
        return (predictions(log_values) - obs.value.to_numpy()) / obs.sd.to_numpy()

    lower, upper = x0 - np.log(100.0), x0 + np.log(100.0)
    result = least_squares(
        residual, x0, bounds=(lower, upper), loss="soft_l1", max_nfev=max_nfev,
    )
    fitted = unpack(result.x)
    prediction = predictions(result.x)
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
            len(obs) < 5 * len(names) or not local["locally_identifiable"]
        ),
        "research_only": True,
        "limitation": (
            "Local calibration of an assay-to-state bridge; validity depends on quantitative "
            "target measurements and does not establish clinical utility."
        ),
    }
    compared = obs.assign(prediction=prediction, weighted_residual=weighted)
    return fitted, diagnostics, compared


def save_mapping_calibration(
    output: str | Path,
    fitted: FeatureStateMapping,
    diagnostics: dict,
    compared: pd.DataFrame,
) -> None:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "fitted_feature_state_mapping.json").write_text(
        json.dumps(asdict(fitted), indent=2) + "\n"
    )
    (output / "mapping_calibration_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )
    compared.to_csv(output / "mapping_calibration_residuals.csv", index=False)
