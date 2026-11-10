"""Outcome-free anti-PD-1 regimen sweeps for the liver-tumour--TDLN twin."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .hybrid import FeatureStateMapping
from .liver_bridge_features import (FeatureBridgeMapping, feature_state_mapping_draws,
                                    liver_bridge_parameter_draws,
                                    load_feature_state_mapping,
                                    load_liver_bridge_parameters,
                                    load_spatial_parameters,
                                    load_systemic_parameters,
                                    spatial_parameter_draws,
                                    systemic_parameter_draws,
                                    simulate_patient_liver_tdln_counterfactual)
from .liver_mapping_assimilation import (feature_bridge_mapping_draws,
                                         load_feature_bridge_mapping)
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
    uncertainty_scenarios: tuple[dict, ...] | None = None,
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
    if uncertainty_scenarios:
        sampled = []
        for scenario in uncertainty_scenarios:
            if not isinstance(scenario, dict):
                raise ValueError("uncertainty scenarios must be parameter dictionaries")
            sampled.append(sweep_patient_regimens(
                row, dose_amounts=amounts, intervals=gaps, days=days, step=step,
                mapping=scenario.get("mapping", mapping),
                base_coupling=scenario.get("base_coupling", base_coupling),
                base_bridge=scenario.get("base_bridge", base_bridge),
                base_spatial=scenario.get("base_spatial", base_spatial),
                base_liver=scenario.get("base_liver", base_liver),
                base_immune=scenario.get("base_immune", base_immune),
                feature_state_mapping=scenario.get(
                    "feature_state_mapping", feature_state_mapping),
                portal_flow_multiplier=portal_flow_multiplier,
                arterial_flow_multiplier=arterial_flow_multiplier,
            ))
        for column in ENDPOINT_COLUMNS:
            values = np.stack([frame[column].to_numpy(dtype=float) for frame in sampled])
            result[f"{column}_lower_95"] = np.quantile(values, 0.025, axis=0)
            result[f"{column}_median"] = np.quantile(values, 0.5, axis=0)
            result[f"{column}_upper_95"] = np.quantile(values, 0.975, axis=0)
        result["pareto_frequency"] = np.mean(np.stack([
            frame.pareto_front.to_numpy(dtype=float) for frame in sampled
        ]), axis=0)
        result["uncertainty_samples"] = len(sampled)
        result["uncertainty_interpretation"] = (
            "calibration_covariance_propagation_not_predictive_coverage"
        )
    result["interpretation"] = (
        "mechanistic_research_grid_not_clinical_dose_recommendation"
    )
    return result


def regimen_uncertainty_scenarios(
    *, n: int, seed: int,
    mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    bridge_diagnostics: dict | None = None,
    mapping_diagnostics: dict | None = None,
    systemic_diagnostics: dict | None = None,
    spatial_diagnostics: dict | None = None,
    state_mapping_diagnostics: dict | None = None,
) -> tuple[dict, ...]:
    """Build aligned deterministic draws from independent calibration layers."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("uncertainty sample count must be a positive integer")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("uncertainty seed must be a non-negative integer")
    diagnostics = (bridge_diagnostics, mapping_diagnostics, systemic_diagnostics,
                   spatial_diagnostics, state_mapping_diagnostics)
    if not any(item is not None for item in diagnostics):
        raise ValueError("uncertainty propagation requires calibration diagnostics")
    liver, coupling, bridge = (base_liver or LiverTumourParameters(),
                               base_coupling or LiverTDLNCouplingParameters(),
                               base_bridge or LiverTDLNSpatialBridge())
    mapping = mapping or FeatureBridgeMapping()
    spatial = base_spatial or SpatialParameters()
    immune = base_immune or TwinParameters()
    state_mapping = feature_state_mapping or FeatureStateMapping()
    bridge_draws = (liver_bridge_parameter_draws(
        liver, coupling, bridge, bridge_diagnostics, n, seed)
        if bridge_diagnostics is not None else ((liver, coupling, bridge),) * n)
    mapping_draws = (feature_bridge_mapping_draws(
        mapping, mapping_diagnostics, n, seed + 1_000_003)
        if mapping_diagnostics is not None else (mapping,) * n)
    immune_draws = (systemic_parameter_draws(
        immune, systemic_diagnostics, n, seed + 2_000_003)
        if systemic_diagnostics is not None else (immune,) * n)
    spatial_draws = (spatial_parameter_draws(
        spatial, spatial_diagnostics, n, seed + 3_000_017)
        if spatial_diagnostics is not None else (spatial,) * n)
    state_draws = (feature_state_mapping_draws(
        state_mapping, state_mapping_diagnostics, n, seed + 4_000_037)
        if state_mapping_diagnostics is not None else (state_mapping,) * n)
    return tuple({
        "base_liver": bridge_draws[index][0],
        "base_coupling": bridge_draws[index][1],
        "base_bridge": bridge_draws[index][2],
        "mapping": mapping_draws[index],
        "base_immune": immune_draws[index],
        "base_spatial": spatial_draws[index],
        "feature_state_mapping": state_draws[index],
    } for index in range(n))


def write_regimen_sweep_audit(
    path: str | Path, *, features_path: str | Path, output_path: str | Path,
    patient_id: str, settings: dict, mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    source_paths: dict[str, str | Path | None] | None = None,
) -> Path:
    sources = {
        name: {"file": Path(source).name, "sha256": sha256_file(source)}
        for name, source in (source_paths or {}).items() if source is not None
    }
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
        "parameter_sources": sources,
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
    source_paths: dict[str, str | Path | None] | None = None,
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
    declared_sources = audit.get("parameter_sources")
    if not isinstance(declared_sources, dict):
        raise ValueError("regimen sweep audit lacks parameter source contract")
    supplied_sources = {
        name: Path(source) for name, source in (source_paths or {}).items()
        if source is not None
    }
    if set(supplied_sources) != set(declared_sources):
        raise ValueError("regimen sweep parameter source paths do not match audit")
    for name, record in declared_sources.items():
        source = supplied_sources[name]
        if (not isinstance(record, dict) or record.get("file") != source.name
                or record.get("sha256") != sha256_file(source)):
            raise ValueError(f"regimen sweep parameter source mismatch: {name}")
    features = pd.read_csv(features_path)
    selected = features[features.patient_id.astype(str).eq(str(audit["patient_id"]))]
    if len(selected) != 1:
        raise ValueError("regimen sweep patient must select exactly one feature row")
    try:
        p, s = audit["base_parameters"], audit["settings"]
        loaded = {}
        if "systemic_parameters" in supplied_sources:
            loaded["immune"] = asdict(load_systemic_parameters(
                supplied_sources["systemic_parameters"]))
        if "spatial_parameters" in supplied_sources:
            loaded["spatial"] = asdict(load_spatial_parameters(
                supplied_sources["spatial_parameters"]))
        if "liver_bridge_parameters" in supplied_sources:
            liver, coupling, bridge = load_liver_bridge_parameters(
                supplied_sources["liver_bridge_parameters"])
            loaded.update({"liver": asdict(liver), "coupling": asdict(coupling),
                           "spatial_bridge": asdict(bridge)})
        if "feature_bridge_mapping" in supplied_sources:
            loaded["mapping"] = asdict(load_feature_bridge_mapping(
                supplied_sources["feature_bridge_mapping"]))
        if "feature_state_mapping" in supplied_sources:
            loaded["feature_state_mapping"] = asdict(load_feature_state_mapping(
                supplied_sources["feature_state_mapping"]))
        if any(p.get(name) != values for name, values in loaded.items()):
            raise ValueError("parameter source content differs from replay parameters")
        diagnostic_keys = {
            "bridge_diagnostics": "bridge_diagnostics",
            "feature_bridge_diagnostics": "mapping_diagnostics",
            "systemic_diagnostics": "systemic_diagnostics",
            "spatial_diagnostics": "spatial_diagnostics",
            "feature_state_diagnostics": "state_mapping_diagnostics",
        }
        diagnostic_payloads = {
            argument: json.loads(supplied_sources[source].read_text())
            for source, argument in diagnostic_keys.items()
            if source in supplied_sources
        }
        samples = int(s.get("uncertainty_samples", 0))
        if samples:
            scenarios = regimen_uncertainty_scenarios(
                n=samples, seed=int(s.get("uncertainty_seed", 42)),
                mapping=FeatureBridgeMapping(**p["mapping"]),
                base_coupling=LiverTDLNCouplingParameters(**p["coupling"]),
                base_bridge=LiverTDLNSpatialBridge(**p["spatial_bridge"]),
                base_spatial=SpatialParameters(**p["spatial"]),
                base_liver=LiverTumourParameters(**p["liver"]),
                base_immune=TwinParameters(**p["immune"]),
                feature_state_mapping=FeatureStateMapping(**p["feature_state_mapping"]),
                **diagnostic_payloads,
            )
        else:
            if diagnostic_payloads:
                raise ValueError("diagnostics supplied without uncertainty samples")
            scenarios = None
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
            uncertainty_scenarios=scenarios,
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
