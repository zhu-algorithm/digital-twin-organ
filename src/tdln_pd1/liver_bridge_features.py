"""Label-free single-cell features → liver–TDLN bridge hypotheses."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .hybrid import (_outcome_like_columns, initial_state_from_features,
                     spatial_initial_state_from_features)
from .liver_tdln_coupling import LiverTDLNCouplingParameters
from .liver_tdln_spatial import (LiverTDLNSpatialBridge,
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
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run the full organ/system/six-zone model from one outcome-free feature row."""
    patient = patient_liver_bridge_from_features(row, mapping=mapping)
    immune_initial = initial_state_from_features(row)
    spatial_initial = spatial_initial_state_from_features(row)
    systemic, spatial = simulate_liver_tdln_spatial(
        days=days, step=step, doses=doses, bridge=patient.spatial_bridge,
        spatial_parameters=patient.spatial_parameters,
        spatial_initial=spatial_initial, coupling_parameters=patient.coupling,
        immune_initial=immune_initial,
    )
    return systemic, spatial, patient.audit

