"""Label-free single-cell features → liver–TDLN bridge hypotheses."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .hybrid import (_outcome_like_columns, initial_state_from_features,
                     spatial_initial_state_from_features)
from .liver_tdln_coupling import LiverTDLNCouplingParameters
from .liver_tdln_coupling import _validate_coupling
from .liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                 _validate_bridge,
                                 simulate_liver_tdln_spatial)
from .provenance import sha256_file
from .simulator import Dose
from .spatial_dynamics import SpatialParameters


@dataclass(frozen=True)
class FeatureBridgeMapping:
    lower_multiplier: float = 0.25
    upper_multiplier: float = 4.0
    drainage_dendritic_coefficient: float = 0.50
    drainage_antigen_proxy_coefficient: float = 0.50
    dc_yield_dendritic_coefficient: float = 0.70
    killing_cytotoxic_coefficient: float = 0.80
    killing_exhaustion_coefficient: float = 0.60
    killing_cd8_fraction_coefficient: float = 1.00
    reference_cd8_fraction: float = 0.15


@dataclass(frozen=True)
class PatientLiverBridge:
    coupling: LiverTDLNCouplingParameters
    spatial_bridge: LiverTDLNSpatialBridge
    spatial_parameters: SpatialParameters
    audit: dict


def verify_patient_liver_tdln_audit(
    output_directory: str | Path, features_path: str | Path,
) -> dict:
    """Recompute every recorded file hash and reject incomplete audits."""
    output_directory = Path(output_directory)
    features_path = Path(features_path)
    audit_path = output_directory / "patient_liver_tdln_audit.json"
    audit = json.loads(audit_path.read_text())
    if audit.get("schema_version") != 1 or audit.get("labels_used") is not False:
        raise ValueError("invalid patient liver-TDLN audit contract")
    if audit.get("input", {}).get("sha256") != sha256_file(features_path):
        raise ValueError("patient liver-TDLN input hash mismatch")
    outputs = audit.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("patient liver-TDLN audit has no outputs")
    for name, expected in outputs.items():
        path = output_directory / Path(name).name
        if sha256_file(path) != expected:
            raise ValueError(f"patient liver-TDLN output hash mismatch: {path.name}")
    return audit


def _validate_mapping(mapping: FeatureBridgeMapping) -> None:
    values = np.asarray(list(asdict(mapping).values()), dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("feature-bridge mapping values must be finite and non-negative")
    if mapping.lower_multiplier <= 0 or mapping.upper_multiplier <= mapping.lower_multiplier:
        raise ValueError("feature-bridge multiplier bounds must be positive and ordered")
    if not 0 <= mapping.reference_cd8_fraction <= 1:
        raise ValueError("reference_cd8_fraction must lie within [0, 1]")


def _read(row: pd.Series, aliases: tuple[str, ...], default: float,
          used: dict, missing: list[str]) -> float:
    for name in aliases:
        if name in row and pd.notna(row[name]):
            try:
                value = float(row[name])
            except (TypeError, ValueError) as error:
                raise ValueError(f"patient bridge feature must be numeric: {name}") from error
            if not np.isfinite(value):
                raise ValueError(f"patient bridge feature must be finite: {name}")
            used[name] = value
            return value
    missing.append(aliases[0])
    return default


def _multiplier(score: float, mapping: FeatureBridgeMapping) -> float:
    return float(np.clip(np.exp(np.clip(score, -20.0, 20.0)),
                         mapping.lower_multiplier, mapping.upper_multiplier))


def patient_liver_bridge_from_features(
    row: pd.Series,
    *,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    mapping: FeatureBridgeMapping | None = None,
) -> PatientLiverBridge:
    """Create bounded bridge parameters without using treatment outcomes."""
    mapping = mapping or FeatureBridgeMapping()
    _validate_mapping(mapping)
    forbidden = _outcome_like_columns(pd.DataFrame([row]))
    if forbidden:
        raise ValueError(f"patient bridge features must be outcome-free: {forbidden}")
    coupling = base_coupling or LiverTDLNCouplingParameters()
    bridge = base_bridge or LiverTDLNSpatialBridge()
    spatial = base_spatial or SpatialParameters()
    used, missing = {}, []
    dendritic = _read(row, ("sig_dendritic_activation_dc", "sig_dendritic_activation"),
                      0.0, used, missing)
    cytotoxic = _read(row, ("sig_cytotoxic_cd8_t", "sig_cytotoxic"),
                      0.0, used, missing)
    exhaustion = _read(row, ("sig_exhaustion_cd8_t", "sig_exhaustion"),
                       0.0, used, missing)
    cd8_fraction = _read(row, ("frac_cd8_t", "frac_cd8_t_cell"),
                         mapping.reference_cd8_fraction, used, missing)
    if not 0 <= cd8_fraction <= 1:
        raise ValueError("patient bridge CD8 fraction must lie within [0, 1]")
    antigen_proxy = _read(row, ("tumour_antigen_proxy",), 0.0, used, missing)
    if antigen_proxy < 0:
        raise ValueError("tumour_antigen_proxy must be non-negative")

    drainage_multiplier = _multiplier(
        mapping.drainage_dendritic_coefficient * dendritic
        + mapping.drainage_antigen_proxy_coefficient * antigen_proxy, mapping)
    dc_multiplier = _multiplier(
        mapping.dc_yield_dendritic_coefficient * dendritic, mapping)
    killing_multiplier = _multiplier(
        mapping.killing_cytotoxic_coefficient * cytotoxic
        - mapping.killing_exhaustion_coefficient * exhaustion
        + mapping.killing_cd8_fraction_coefficient
        * (cd8_fraction - mapping.reference_cd8_fraction), mapping)
    mapped_coupling = replace(
        coupling, immune_killing_gain=coupling.immune_killing_gain * killing_multiplier)
    mapped_bridge = replace(
        bridge, antigen_drainage_gain=bridge.antigen_drainage_gain * drainage_multiplier)
    mapped_spatial = replace(spatial, dc_yield=spatial.dc_yield * dc_multiplier)
    audit = {
        "schema_version": 1,
        "labels_used": False,
        "used_features": used,
        "fallback_features": sorted(set(missing)),
        "multipliers": {"antigen_drainage": drainage_multiplier,
                        "dc_yield": dc_multiplier,
                        "immune_killing": killing_multiplier},
        "mapping": asdict(mapping),
        "research_only": True,
        "limitation": "Bounded feature-to-parameter hypothesis; coefficients require calibration.",
    }
    return PatientLiverBridge(mapped_coupling, mapped_bridge, mapped_spatial, audit)


def simulate_patient_liver_tdln(
    row: pd.Series, *, days: float = 28.0, step: float = 0.25,
    doses: tuple[Dose, ...] = (), mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run the full organ/system/six-zone model from one outcome-free feature row."""
    patient = patient_liver_bridge_from_features(
        row, mapping=mapping, base_coupling=base_coupling,
        base_bridge=base_bridge, base_spatial=base_spatial)
    immune_initial = initial_state_from_features(row)
    spatial_initial = spatial_initial_state_from_features(row)
    systemic, spatial = simulate_liver_tdln_spatial(
        days=days, step=step, doses=doses, bridge=patient.spatial_bridge,
        spatial_parameters=patient.spatial_parameters,
        spatial_initial=spatial_initial, coupling_parameters=patient.coupling,
        immune_initial=immune_initial,
    )
    return systemic, spatial, patient.audit


def _time_integral(frame: pd.DataFrame, column: str) -> float:
    time = frame.day.to_numpy(dtype=float)
    values = frame[column].to_numpy(dtype=float)
    return float(np.sum(np.diff(time) * (values[:-1] + values[1:]) / 2.0))


def _spatial_endpoint(frame: pd.DataFrame, species: str, zone: str) -> float:
    selected = frame[frame.species.eq(species) & frame.zone.eq(zone)]
    if selected.empty:
        raise ValueError(f"missing spatial endpoint: {species}@{zone}")
    return float(selected.sort_values("day").amount.iloc[-1])


def simulate_patient_liver_tdln_counterfactual(
    row: pd.Series, *, doses: tuple[Dose, ...], days: float = 28.0,
    step: float = 0.25, mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    """Run matched natural-history and anti-PD-1 arms for one patient.

    The returned effect summary is a mechanistic model contrast, not a
    calibrated probability or clinical response classification.
    """
    if not doses or not any(dose.amount > 0 for dose in doses):
        raise ValueError("patient counterfactual requires an explicit positive dose schedule")
    frames_systemic, frames_spatial = [], []
    arm_runs: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    audit: dict | None = None
    for arm, schedule in (("control", ()), ("anti_pd1", doses)):
        systemic, spatial, arm_audit = simulate_patient_liver_tdln(
            row, days=days, step=step, doses=schedule, mapping=mapping,
            base_coupling=base_coupling, base_bridge=base_bridge,
            base_spatial=base_spatial)
        systemic, spatial = systemic.copy(), spatial.copy()
        systemic["arm"], spatial["arm"] = arm, arm
        frames_systemic.append(systemic)
        frames_spatial.append(spatial)
        arm_runs[arm] = systemic, spatial
        audit = arm_audit

    control, control_spatial = arm_runs["control"]
    treated, treated_spatial = arm_runs["anti_pd1"]
    control_endpoint = float(control.liver_tumour_fraction.iloc[-1])
    treated_endpoint = float(treated.liver_tumour_fraction.iloc[-1])
    relative_reduction = ((control_endpoint - treated_endpoint)
                          / max(control_endpoint, 1e-12))
    effect = {
        "schema_version": 1,
        "interpretation": "mechanistic_matched_counterfactual_not_clinical_response",
        "control_final_liver_tumour_fraction": control_endpoint,
        "anti_pd1_final_liver_tumour_fraction": treated_endpoint,
        "final_liver_tumour_fraction_absolute_change": treated_endpoint - control_endpoint,
        "final_liver_tumour_fraction_relative_reduction": relative_reduction,
        "control_effector_tumour_auc": _time_integral(control, "effector_tumour"),
        "anti_pd1_effector_tumour_auc": _time_integral(treated, "effector_tumour"),
        "anti_pd1_peak_pd1_occupancy": float(treated.pd1_occupancy.max()),
        "control_efferent_effector_endpoint": _spatial_endpoint(
            control_spatial, "effector", "efferent_lymph"),
        "anti_pd1_efferent_effector_endpoint": _spatial_endpoint(
            treated_spatial, "effector", "efferent_lymph"),
    }
    effect["effector_tumour_auc_fold_change"] = (
        effect["anti_pd1_effector_tumour_auc"]
        / max(effect["control_effector_tumour_auc"], 1e-12))
    effect["efferent_effector_endpoint_fold_change"] = (
        effect["anti_pd1_efferent_effector_endpoint"]
        / max(effect["control_efferent_effector_endpoint"], 1e-12))
    return (pd.concat(frames_systemic, ignore_index=True),
            pd.concat(frames_spatial, ignore_index=True), effect, audit or {})


COUNTERFACTUAL_FEATURES = {
    "control_final_liver_tumour_fraction": "liver_cf_control_final_tumour_fraction",
    "anti_pd1_final_liver_tumour_fraction": "liver_cf_anti_pd1_final_tumour_fraction",
    "final_liver_tumour_fraction_absolute_change": "liver_cf_tumour_absolute_change",
    "final_liver_tumour_fraction_relative_reduction": "liver_cf_tumour_relative_reduction",
    "control_effector_tumour_auc": "liver_cf_control_effector_tumour_auc",
    "anti_pd1_effector_tumour_auc": "liver_cf_anti_pd1_effector_tumour_auc",
    "effector_tumour_auc_fold_change": "liver_cf_effector_tumour_auc_fold_change",
    "anti_pd1_peak_pd1_occupancy": "liver_cf_peak_pd1_occupancy",
    "control_efferent_effector_endpoint": "liver_cf_control_efferent_effector_endpoint",
    "anti_pd1_efferent_effector_endpoint": "liver_cf_anti_pd1_efferent_effector_endpoint",
    "efferent_effector_endpoint_fold_change": "liver_cf_efferent_effector_fold_change",
}


def augment_with_liver_tdln_counterfactuals(
    features: pd.DataFrame, *, doses: tuple[Dose, ...], days: float = 28.0,
    step: float = 0.25, mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
) -> pd.DataFrame:
    """Append label-free matched liver–TDLN contrasts for every patient."""
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    forbidden = _outcome_like_columns(features)
    if forbidden:
        raise ValueError(f"counterfactual augmentation must be outcome-free: {forbidden}")
    out = features.copy()
    if out.patient_id.isna().any():
        raise ValueError("features need non-missing patient_id")
    out["patient_id"] = out.patient_id.astype(str).str.strip()
    if out.patient_id.eq("").any() or out.patient_id.duplicated().any():
        raise ValueError("features must have one unique, non-empty row per patient")
    rows = []
    for _, row in out.iterrows():
        _, _, effect, _ = simulate_patient_liver_tdln_counterfactual(
            row, doses=doses, days=days, step=step, mapping=mapping,
            base_coupling=base_coupling, base_bridge=base_bridge,
            base_spatial=base_spatial)
        values = {target: float(effect[source])
                  for source, target in COUNTERFACTUAL_FEATURES.items()}
        if not np.isfinite(list(values.values())).all():
            raise RuntimeError(f"non-finite counterfactual feature for patient {row.patient_id}")
        rows.append(values)
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def verify_liver_counterfactual_augmentation(
    audit_path: str | Path, features_path: str | Path, output_path: str | Path,
    bridge_parameters_path: str | Path | None = None,
) -> dict:
    """Verify hashes and replay a cohort counterfactual augmentation exactly."""
    audit_path, features_path, output_path = map(
        Path, (audit_path, features_path, output_path))
    audit = json.loads(audit_path.read_text())
    if (audit.get("schema_version") != 1 or audit.get("labels_used") is not False
            or audit.get("stage") != "label_free_liver_tdln_counterfactual_augmentation"):
        raise ValueError("invalid liver counterfactual augmentation audit contract")
    if audit.get("input", {}).get("sha256") != sha256_file(features_path):
        raise ValueError("liver counterfactual augmentation input hash mismatch")
    if audit.get("output", {}).get("sha256") != sha256_file(output_path):
        raise ValueError("liver counterfactual augmentation output hash mismatch")
    parameter_record = audit.get("bridge_parameters")
    if parameter_record is not None:
        if bridge_parameters_path is None:
            raise ValueError("bridge parameter file is required to replay this augmentation")
        if parameter_record.get("sha256") != sha256_file(bridge_parameters_path):
            raise ValueError("liver counterfactual bridge parameter hash mismatch")
    settings = audit.get("settings", {})
    try:
        doses = tuple(Dose(float(day), float(amount)) for day, amount in zip(
            settings["dose_days"], settings["dose_amounts"], strict=True))
        days, step = float(settings["days"]), float(settings["step"])
        mapping = FeatureBridgeMapping(**audit["mapping"])
        base_coupling = LiverTDLNCouplingParameters(**audit["base_parameters"]["coupling"])
        base_bridge = LiverTDLNSpatialBridge(**audit["base_parameters"]["spatial_bridge"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid liver counterfactual replay settings") from error
    replay = augment_with_liver_tdln_counterfactuals(
        pd.read_csv(features_path), doses=doses, days=days, step=step, mapping=mapping,
        base_coupling=base_coupling, base_bridge=base_bridge)
    replay_sha = hashlib.sha256(
        replay.to_csv(index=False).encode("utf-8")).hexdigest()
    if replay_sha != audit["output"]["sha256"]:
        raise ValueError("liver counterfactual augmentation derivation replay mismatch")
    return audit


def load_liver_bridge_parameters(
    path: str | Path,
) -> tuple[LiverTDLNCouplingParameters, LiverTDLNSpatialBridge]:
    """Load and fail-closed validate a saved liver–TDLN calibration."""
    payload = json.loads(Path(path).read_text())
    if payload.get("research_only") is not True:
        raise ValueError("liver bridge parameter file must declare research_only=true")
    try:
        coupling = LiverTDLNCouplingParameters(**payload["coupling"])
        bridge = LiverTDLNSpatialBridge(**payload["spatial_bridge"])
    except (KeyError, TypeError) as error:
        raise ValueError("invalid liver bridge parameter schema") from error
    _validate_coupling(coupling)
    _validate_bridge(bridge)
    return coupling, bridge

