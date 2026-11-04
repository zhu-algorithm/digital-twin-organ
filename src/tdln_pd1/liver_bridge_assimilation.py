"""Weighted longitudinal calibration of liver-tumour ↔ spatial-TDLN bridges."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .assimilation import _local_calibration_diagnostics
from .liver_tdln_coupling import LiverTDLNCouplingParameters
from .liver_tumor import LiverTumourParameters
from .liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                 simulate_liver_tdln_spatial)
from .simulator import Dose
from .spatial_dynamics import SPECIES, ZONE_NAMES

SYSTEM_TARGETS = {
    "liver_tumour_fraction", "tumour_antigen", "effector_tumour", "pd1_occupancy",
    "liver_portal_flow", "liver_arterial_flow", "liver_arterial_buffer_flow",
    "liver_oxygenated_inflow", "liver_perfusion",
}
SPATIAL_TARGETS = {f"{species}@{zone}" for species in SPECIES for zone in ZONE_NAMES}
TARGETS = SYSTEM_TARGETS | SPATIAL_TARGETS
FIT_PARAMETERS = {
    "coupling.immune_killing_gain",
    "coupling.necrotic_antigen_gain",
    "coupling.hypoxic_antigen_gain",
    "bridge.antigen_drainage_gain",
    "liver.arterial_buffer_fraction",
    "liver.arterial_buffer_max_multiplier",
    "liver.portal_oxygen_content",
    "liver.arterial_oxygen_content",
}


def validate_bridge_observations(observations: pd.DataFrame) -> pd.DataFrame:
    required = {"day", "target", "value", "sd"}
    missing = required.difference(observations.columns)
    if missing:
        raise ValueError(f"missing liver bridge observation columns: {sorted(missing)}")
    out = observations.copy()
    if out.empty:
        raise ValueError("liver bridge observations must not be empty")
    if out.target.isna().any():
        raise ValueError("liver bridge targets must be non-missing and non-empty")
    out["target"] = out.target.astype(str).str.strip()
    if out.target.eq("").any():
        raise ValueError("liver bridge targets must be non-missing and non-empty")
    for column in ("day", "value", "sd"):
        out[column] = pd.to_numeric(out[column], errors="raise")
    if not np.isfinite(out[["day", "value", "sd"]]).all().all():
        raise ValueError("liver bridge observations must be finite")
    if out.day.lt(0).any() or out.value.lt(0).any() or out.sd.le(0).any():
        raise ValueError("observation days/values must be non-negative and sd strictly positive")
    unknown = sorted(set(out.target) - TARGETS)
    if unknown:
        raise ValueError(f"unknown liver bridge targets: {unknown}")
    if out.duplicated(["day", "target"]).any():
        raise ValueError("liver bridge observations need unique day-target rows")
    return out.sort_values(["day", "target"]).reset_index(drop=True)


def _predictions(systemic: pd.DataFrame, spatial: pd.DataFrame,
                 observations: pd.DataFrame) -> np.ndarray:
    values = []
    for row in observations.itertuples(index=False):
        if row.target in SYSTEM_TARGETS:
            values.append(np.interp(row.day, systemic.day, systemic[row.target]))
        else:
            species, zone = row.target.split("@", 1)
            series = spatial[spatial.species.eq(species) & spatial.zone.eq(zone)]
            values.append(np.interp(row.day, series.day, series.amount))
    return np.asarray(values, dtype=float)


def calibrate_liver_bridge(
    observations: pd.DataFrame,
    *,
    fit_parameters: Iterable[str] = ("bridge.antigen_drainage_gain",),
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_liver: LiverTumourParameters | None = None,
    doses: tuple[Dose, ...] = (),
    portal_flow_multiplier: float | None = None,
    arterial_flow_multiplier: float | None = None,
    step: float = 0.5,
    max_nfev: int = 100,
) -> tuple[LiverTumourParameters, LiverTDLNCouplingParameters,
           LiverTDLNSpatialBridge, dict, pd.DataFrame]:
    """Fit positive bridge parameters in log space to weighted observations."""
    obs = validate_bridge_observations(observations)
    names = tuple(fit_parameters)
    invalid = sorted(set(names) - FIT_PARAMETERS)
    if not names or invalid:
        raise ValueError(f"invalid liver bridge fit parameters: {invalid}")
    if len(set(names)) != len(names):
        raise ValueError("liver bridge fit parameters must be unique")
    buffer_names = {"liver.arterial_buffer_fraction",
                    "liver.arterial_buffer_max_multiplier"}
    if set(names) & buffer_names and (
            portal_flow_multiplier is None or arterial_flow_multiplier is not None):
        raise ValueError(
            "arterial buffer calibration requires an isolated portal perturbation: "
            "set portal_flow_multiplier and leave arterial_flow_multiplier unspecified")
    if not np.isfinite(step) or step <= 0 or max_nfev < 1:
        raise ValueError("calibration step must be positive and max_nfev >= 1")
    coupling = base_coupling or LiverTDLNCouplingParameters()
    bridge = base_bridge or LiverTDLNSpatialBridge()
    liver = base_liver or LiverTumourParameters()

    def get(name):
        group, field = name.split(".", 1)
        source = {"coupling": coupling, "bridge": bridge, "liver": liver}[group]
        return getattr(source, field)

    initial_values = np.asarray([get(name) for name in names], dtype=float)
    if not np.isfinite(initial_values).all() or np.any(initial_values <= 0):
        raise ValueError("fitted bridge parameters need finite, strictly positive base values")
    x0 = np.log(initial_values)
    lower, upper = x0 - np.log(100.0), x0 + np.log(100.0)
    for index, name in enumerate(names):
        if name == "liver.arterial_buffer_fraction":
            lower[index], upper[index] = np.log(1e-6), np.log(1.0)
        elif name == "liver.arterial_buffer_max_multiplier":
            lower[index], upper[index] = np.log(1.0), np.log(10.0)
    days = max(float(obs.day.max()), step)

    def unpack(log_values):
        coupling_updates, bridge_updates, liver_updates = {}, {}, {}
        for name, value in zip(names, np.exp(log_values)):
            group, field = name.split(".", 1)
            {"coupling": coupling_updates, "bridge": bridge_updates,
             "liver": liver_updates}[group][field] = value
        return (replace(liver, **liver_updates), replace(coupling, **coupling_updates),
                replace(bridge, **bridge_updates))

    def run(log_values):
        fitted_liver, fitted_coupling, fitted_bridge = unpack(log_values)
        return simulate_liver_tdln_spatial(
            days=days, step=step, doses=doses, bridge=fitted_bridge,
            coupling_parameters=fitted_coupling, liver_parameters=fitted_liver,
            portal_flow_multiplier=portal_flow_multiplier,
            arterial_flow_multiplier=arterial_flow_multiplier)

    def residual(log_values):
        systemic, spatial = run(log_values)
        return (_predictions(systemic, spatial, obs) - obs.value.to_numpy()) / obs.sd.to_numpy()

    result = least_squares(residual, x0, bounds=(lower, upper), loss="soft_l1",
                           diff_step=1e-3, max_nfev=max_nfev)
    fitted_liver, fitted_coupling, fitted_bridge = unpack(result.x)
    systemic, spatial = run(result.x)
    prediction = _predictions(systemic, spatial, obs)
    weighted = (prediction - obs.value.to_numpy()) / obs.sd.to_numpy()
    local = _local_calibration_diagnostics(result, names, lower, upper, len(obs))
    diagnostics = {
        "success": bool(result.success), "message": result.message,
        "n_observations": int(len(obs)), "n_parameters": len(names),
        "nfev": int(result.nfev),
        "weighted_rmse": float(np.sqrt(np.mean(weighted ** 2))),
        "fit_parameters": {name: float(np.exp(value)) for name, value in zip(names, result.x)},
        **local,
        "identifiability_warning": bool(
            len(obs) < 5 * len(names) or not local["locally_identifiable"]),
        "research_only": True,
        "limitation": (
            "Local fit of dimensional bridge hypotheses; external biological and clinical "
            "validation remains required."
        ),
    }
    compared = obs.assign(prediction=prediction, weighted_residual=weighted)
    return fitted_liver, fitted_coupling, fitted_bridge, diagnostics, compared


def save_liver_bridge_calibration(output: str | Path,
                                  liver: LiverTumourParameters,
                                  coupling: LiverTDLNCouplingParameters,
                                  bridge: LiverTDLNSpatialBridge,
                                  diagnostics: dict,
                                  compared: pd.DataFrame) -> None:
    target = Path(output)
    target.mkdir(parents=True, exist_ok=True)
    (target / "fitted_liver_tdln_bridge.json").write_text(json.dumps({
        "liver_tumour": asdict(liver), "coupling": asdict(coupling),
        "spatial_bridge": asdict(bridge),
        "research_only": True,
    }, indent=2) + "\n")
    (target / "liver_tdln_bridge_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n")
    compared.to_csv(target / "liver_tdln_bridge_residuals.csv", index=False)

