"""Outcome-free anti-PD-1 regimen sweeps for the liver-tumour--TDLN twin."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .hybrid import FeatureStateMapping
from .liver_bridge_features import (FeatureBridgeMapping,
                                    simulate_patient_liver_tdln_counterfactual)
from .liver_tdln_coupling import LiverTDLNCouplingParameters
from .liver_tdln_spatial import LiverTDLNSpatialBridge
from .liver_tumor import LiverTumourParameters
from .provenance import sha256_file
from .simulator import Dose, TwinParameters
from .spatial_dynamics import SpatialParameters


ENDPOINT_COLUMNS = (
    "tumour_relative_reduction",
    "normal_cell_preservation",
    "normal_oxygen_change",
    "perfusion_relative_change",
)


def _positive_grid(values, name: str) -> tuple[float, ...]:
    try:
        grid = tuple(float(value) for value in values)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must contain finite positive values") from error
    if not grid or not np.isfinite(grid).all() or any(value <= 0 for value in grid):
        raise ValueError(f"{name} must contain finite positive values")
    if len(set(grid)) != len(grid):
        raise ValueError(f"{name} must not contain duplicates")
    return tuple(sorted(grid))


def _pareto_front(frame: pd.DataFrame) -> np.ndarray:
    """Maximise four organ endpoints while minimising total model exposure."""
    objectives = frame[list(ENDPOINT_COLUMNS)].to_numpy(dtype=float)
    exposure = frame.total_model_exposure.to_numpy(dtype=float)
    front = np.ones(len(frame), dtype=bool)
    for candidate in range(len(frame)):
        no_worse = ((objectives >= objectives[candidate]).all(axis=1)
                    & (exposure <= exposure[candidate]))
        strictly_better = ((objectives > objectives[candidate]).any(axis=1)
                           | (exposure < exposure[candidate]))
        no_worse[candidate] = False
        if np.any(no_worse & strictly_better):
            front[candidate] = False
    return front


def sweep_patient_regimens(
    row: pd.Series, *, dose_amounts=(0.5, 1.0, 2.0), intervals=(7.0, 14.0),
    days: float = 28.0, step: float = 0.25,
    mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    portal_flow_multiplier: float | None = None,
    arterial_flow_multiplier: float | None = None,
) -> pd.DataFrame:
    """Evaluate a prespecified regimen grid without using response labels.

    Dose amounts are dimensionless model exposure units.  Pareto membership is
    descriptive and must not be interpreted as a clinical dose recommendation.
    """
    amounts = _positive_grid(dose_amounts, "dose_amounts")
    gaps = _positive_grid(intervals, "intervals")
    if not np.isfinite(days) or not np.isfinite(step) or days <= 0 or step <= 0:
        raise ValueError("days and step must be finite and positive")
    records = []
    for interval in gaps:
        dose_days = tuple(float(day) for day in np.arange(0.0, days, interval))
        for amount in amounts:
            doses = tuple(Dose(day, amount) for day in dose_days)
            _, _, effect, _ = simulate_patient_liver_tdln_counterfactual(
                row, doses=doses, days=days, step=step, mapping=mapping,
                base_coupling=base_coupling, base_bridge=base_bridge,
                base_spatial=base_spatial, base_liver=base_liver,
                base_immune=base_immune,
                feature_state_mapping=feature_state_mapping,
                portal_flow_multiplier=portal_flow_multiplier,
                arterial_flow_multiplier=arterial_flow_multiplier,
            )
            records.append({
                "dose_amount_model_units": amount,
                "dose_interval_days": interval,
                "dose_count": len(doses),
                "total_model_exposure": amount * len(doses),
                "dose_days": ",".join(format(day, ".12g") for day in dose_days),
                "tumour_relative_reduction": effect[
                    "final_liver_tumour_fraction_relative_reduction"],
                "normal_cell_preservation": effect[
                    "final_liver_normal_cell_preservation"],
                "normal_oxygen_change": effect["min_liver_normal_oxygen_change"],
                "perfusion_relative_change": effect[
                    "final_liver_perfusion_relative_change"],
                "peak_pd1_occupancy": effect["anti_pd1_peak_pd1_occupancy"],
                "effector_tumour_auc_fold_change": effect[
                    "effector_tumour_auc_fold_change"],
                "efferent_effector_endpoint_fold_change": effect[
                    "efferent_effector_endpoint_fold_change"],
            })
    result = pd.DataFrame(records).sort_values(
        ["total_model_exposure", "dose_interval_days", "dose_amount_model_units"],
        kind="stable",
    ).reset_index(drop=True)
    if not np.isfinite(result.drop(columns="dose_days").to_numpy(dtype=float)).all():
        raise RuntimeError("regimen sweep produced non-finite endpoints")
    result["pareto_front"] = _pareto_front(result)
    result["interpretation"] = (
        "mechanistic_research_grid_not_clinical_dose_recommendation"
    )
    return result


def write_regimen_sweep_audit(
    path: str | Path, *, features_path: str | Path, output_path: str | Path,
    patient_id: str, settings: dict, mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
) -> Path:
    audit = {
        "schema_version": 1,
        "stage": "label_free_patient_pd1_regimen_sweep",
        "labels_used": False,
        "research_only": True,
        "patient_id": str(patient_id),
        "input": {"file": Path(features_path).name,
                  "sha256": sha256_file(features_path)},
        "output": {"file": Path(output_path).name,
                   "sha256": sha256_file(output_path)},
        "settings": settings,
        "base_parameters": {
            "mapping": asdict(mapping or FeatureBridgeMapping()),
            "coupling": asdict(base_coupling or LiverTDLNCouplingParameters()),
            "spatial_bridge": asdict(base_bridge or LiverTDLNSpatialBridge()),
            "spatial": asdict(base_spatial or SpatialParameters()),
            "liver": asdict(base_liver or LiverTumourParameters()),
            "immune": asdict(base_immune or TwinParameters()),
            "feature_state_mapping": asdict(
                feature_state_mapping or FeatureStateMapping()),
        },
        "pareto_contract": {
            "maximise": list(ENDPOINT_COLUMNS),
            "minimise": ["total_model_exposure"],
            "role": "descriptive_multiobjective_screen",
        },
        "limitation": (
            "Dose amounts are model units; the grid is not pharmacologically calibrated "
            "and Pareto membership is not a clinical recommendation."
        ),
    }
    target = Path(path)
    target.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    return target


def verify_regimen_sweep_audit(
    audit_path: str | Path, features_path: str | Path, output_path: str | Path,
) -> dict:
    audit = json.loads(Path(audit_path).read_text())
    if (audit.get("schema_version") != 1
            or audit.get("stage") != "label_free_patient_pd1_regimen_sweep"
            or audit.get("labels_used") is not False
            or audit.get("research_only") is not True):
        raise ValueError("invalid regimen sweep audit contract")
    if audit.get("input", {}).get("sha256") != sha256_file(features_path):
        raise ValueError("regimen sweep input hash mismatch")
    if audit.get("output", {}).get("sha256") != sha256_file(output_path):
        raise ValueError("regimen sweep output hash mismatch")
    features = pd.read_csv(features_path)
    selected = features[features.patient_id.astype(str).eq(str(audit["patient_id"]))]
    if len(selected) != 1:
        raise ValueError("regimen sweep patient must select exactly one feature row")
    try:
        p, s = audit["base_parameters"], audit["settings"]
        replay = sweep_patient_regimens(
            selected.iloc[0], dose_amounts=s["dose_amounts"], intervals=s["intervals"],
            days=float(s["days"]), step=float(s["step"]),
            mapping=FeatureBridgeMapping(**p["mapping"]),
            base_coupling=LiverTDLNCouplingParameters(**p["coupling"]),
            base_bridge=LiverTDLNSpatialBridge(**p["spatial_bridge"]),
            base_spatial=SpatialParameters(**p["spatial"]),
            base_liver=LiverTumourParameters(**p["liver"]),
            base_immune=TwinParameters(**p["immune"]),
            feature_state_mapping=FeatureStateMapping(**p["feature_state_mapping"]),
            portal_flow_multiplier=s.get("portal_flow_multiplier"),
            arterial_flow_multiplier=s.get("arterial_flow_multiplier"),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid regimen sweep replay settings") from error
    observed = pd.read_csv(output_path)
    try:
        pd.testing.assert_frame_equal(observed, replay, check_exact=False,
                                      rtol=1e-12, atol=1e-12)
    except AssertionError as error:
        raise ValueError("regimen sweep derivation replay mismatch") from error
    return audit
