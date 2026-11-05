"""Label-free single-cell features → liver–TDLN bridge hypotheses."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .assimilation import calibration_parameter_draws
from .hybrid import (FeatureStateMapping, _outcome_like_columns,
                     _validate_mapping as _validate_state_mapping,
                     initial_state_from_features,
                     spatial_initial_state_from_features)
from .liver_tdln_coupling import LiverTDLNCouplingParameters
from .liver_tdln_coupling import _validate_coupling
from .liver_tumor import (LiverTumourParameters,
                          _validate_parameters as _validate_liver_parameters)
from .liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                 _validate_bridge,
                                 simulate_liver_tdln_spatial)
from .provenance import sha256_file
from .simulator import (Dose, InitialState, TwinParameters,
                        _validate_model_inputs)
from .spatial_dynamics import (SpatialParameters,
                               _validate_parameters as _validate_spatial_parameters)


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
    systemic_parameters_path: str | Path | None = None,
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
    systemic_record = audit.get("systemic_parameters")
    if systemic_record is not None:
        if systemic_parameters_path is None:
            raise ValueError("patient liver-TDLN replay requires systemic parameters")
        if systemic_record.get("sha256") != sha256_file(systemic_parameters_path):
            raise ValueError("patient liver-TDLN systemic parameter hash mismatch")
    outputs = audit.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError("patient liver-TDLN audit has no outputs")
    for name, expected in outputs.items():
        path = output_directory / Path(name).name
        if sha256_file(path) != expected:
            raise ValueError(f"patient liver-TDLN output hash mismatch: {path.name}")
    settings = audit.get("settings")
    if settings is not None:
        try:
            table = pd.read_csv(features_path)
            selected = table[table.patient_id.astype(str).eq(str(audit["patient_id"]))]
            if len(selected) != 1:
                raise ValueError("patient audit must select exactly one feature row")
            doses = tuple(Dose(float(day), float(settings["dose"]))
                          for day in settings["dose_days"])
            kwargs = {
                "days": float(settings["days"]), "step": float(settings["step"]),
                "doses": doses,
                "portal_flow_multiplier": settings.get("portal_flow_multiplier"),
                "arterial_flow_multiplier": settings.get("arterial_flow_multiplier"),
                "base_immune": TwinParameters(**audit.get("base_parameters", {}).get(
                    "immune", asdict(TwinParameters()))),
            }
            if settings.get("paired_counterfactual") is True:
                systemic, spatial, effect, _ = simulate_patient_liver_tdln_counterfactual(
                    selected.iloc[0], **kwargs)
            else:
                systemic, spatial, _ = simulate_patient_liver_tdln(selected.iloc[0], **kwargs)
                effect = None
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("invalid patient liver-TDLN replay settings") from error
        replay = {
            "patient_liver_tdln_systemic.csv": hashlib.sha256(
                systemic.to_csv(index=False).encode("utf-8")).hexdigest(),
            "patient_liver_tdln_spatial.csv": hashlib.sha256(
                spatial.to_csv(index=False).encode("utf-8")).hexdigest(),
        }
        if effect is not None:
            replay["patient_liver_tdln_counterfactual.json"] = hashlib.sha256(
                (json.dumps(effect, indent=2, sort_keys=True) + "\n").encode("utf-8")
            ).hexdigest()
        if replay != outputs:
            raise ValueError("patient liver-TDLN derivation replay mismatch")
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
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    portal_flow_multiplier: float | None = None,
    arterial_flow_multiplier: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run the full organ/system/six-zone model from one outcome-free feature row."""
    patient = patient_liver_bridge_from_features(
        row, mapping=mapping, base_coupling=base_coupling,
        base_bridge=base_bridge, base_spatial=base_spatial)
    immune_initial = initial_state_from_features(row, mapping=feature_state_mapping)
    spatial_initial = spatial_initial_state_from_features(
        row, mapping=feature_state_mapping)
    systemic, spatial = simulate_liver_tdln_spatial(
        days=days, step=step, doses=doses, bridge=patient.spatial_bridge,
        spatial_parameters=patient.spatial_parameters,
        spatial_initial=spatial_initial, coupling_parameters=patient.coupling,
        immune_initial=immune_initial, immune_parameters=base_immune,
        liver_parameters=base_liver,
        portal_flow_multiplier=portal_flow_multiplier,
        arterial_flow_multiplier=arterial_flow_multiplier,
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
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    portal_flow_multiplier: float | None = None,
    arterial_flow_multiplier: float | None = None,
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
            base_spatial=base_spatial,
            base_liver=base_liver,
            base_immune=base_immune,
            feature_state_mapping=feature_state_mapping,
            portal_flow_multiplier=portal_flow_multiplier,
            arterial_flow_multiplier=arterial_flow_multiplier)
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
UNCERTAINTY_QUANTILES = ((0.025, "q025"), (0.5, "q50"), (0.975, "q975"))


def liver_bridge_parameter_draws(
    liver: LiverTumourParameters, coupling: LiverTDLNCouplingParameters,
    bridge: LiverTDLNSpatialBridge,
    diagnostics: dict, n: int, seed: int,
) -> tuple[tuple[LiverTumourParameters, LiverTDLNCouplingParameters,
                 LiverTDLNSpatialBridge], ...]:
    """Draw a physiologically truncated correlated ensemble in log space."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("liver bridge uncertainty requires at least two samples")
    covariance_payload = diagnostics.get("local_log_parameter_covariance")
    if not isinstance(covariance_payload, dict) or diagnostics.get(
            "locally_identifiable") is not True:
        raise ValueError("liver bridge diagnostics lack identifiable covariance")
    _validate_liver_parameters(liver)
    _validate_coupling(coupling)
    _validate_bridge(bridge)
    names = covariance_payload.get("parameter_order")
    covariance = np.asarray(covariance_payload.get("log_scale_covariance"), dtype=float)
    available = {
        **{f"liver.{key}": value for key, value in asdict(liver).items()},
        **{f"coupling.{key}": value for key, value in asdict(coupling).items()},
        **{f"bridge.{key}": value for key, value in asdict(bridge).items()},
    }
    fitted = diagnostics.get("fit_parameters")
    if (not isinstance(names, list) or not names or len(set(names)) != len(names)
            or any(name not in available for name in names)
            or covariance.shape != (len(names), len(names))
            or not np.isfinite(covariance).all()
            or not np.allclose(covariance, covariance.T, rtol=1e-8, atol=1e-12)
            or not isinstance(fitted, dict) or set(fitted) != set(names)):
        raise ValueError("invalid namespaced liver bridge covariance")
    centers = np.asarray([float(fitted[name]) for name in names], dtype=float)
    loaded = np.asarray([available[name] for name in names], dtype=float)
    if (not np.isfinite(centers).all() or np.any(centers <= 0)
            or not np.allclose(centers, loaded, rtol=1e-10, atol=0.0)):
        raise ValueError("liver bridge diagnostics and parameter file have different centers")
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues.min() < -1e-10 * max(1.0, float(np.abs(eigenvalues).max())):
        raise ValueError("liver bridge covariance must be positive semidefinite")
    rng = np.random.default_rng(seed)
    draws = []
    attempts, max_attempts = 0, max(1000, 1000 * n)
    while len(draws) < n and attempts < max_attempts:
        batch_size = min(max(4, n - len(draws)), max_attempts - attempts)
        vectors = rng.multivariate_normal(
            np.log(loaded), covariance, size=batch_size, check_valid="raise")
        attempts += batch_size
        with np.errstate(over="ignore", invalid="ignore"):
            candidates = np.exp(vectors)
        for vector in candidates:
            liver_updates, coupling_updates, bridge_updates = {}, {}, {}
            for name, value in zip(names, vector):
                group, field = name.split(".", 1)
                {"liver": liver_updates, "coupling": coupling_updates,
                 "bridge": bridge_updates}[group][field] = float(value)
            candidate = (replace(liver, **liver_updates),
                         replace(coupling, **coupling_updates),
                         replace(bridge, **bridge_updates))
            try:
                _validate_liver_parameters(candidate[0])
                _validate_coupling(candidate[1])
                _validate_bridge(candidate[2])
            except ValueError:
                continue
            draws.append(candidate)
            if len(draws) == n:
                break
    if len(draws) != n:
        raise ValueError(
            "unable to draw a physiological liver bridge ensemble from the "
            "calibration covariance; inspect boundary contact and identifiability")
    return tuple(draws)


def systemic_parameter_draws(
    parameters: TwinParameters, diagnostics: dict, n: int, seed: int,
) -> tuple[TwinParameters, ...]:
    """Draw and validate systemic immune/PK parameters for liver coupling."""
    _, samples = calibration_parameter_draws(parameters, diagnostics, n, seed)
    draws = tuple(replace(parameters, **sample) for sample in samples)
    for draw in draws:
        _validate_model_inputs(draw, InitialState())
    return draws


def spatial_parameter_draws(
    parameters: SpatialParameters, diagnostics: dict, n: int, seed: int,
) -> tuple[SpatialParameters, ...]:
    """Draw and validate six-zone TDLN transport/reaction parameters."""
    _, samples = calibration_parameter_draws(parameters, diagnostics, n, seed)
    draws = tuple(replace(parameters, **sample) for sample in samples)
    for draw in draws:
        _validate_spatial_parameters(draw)
    return draws


def feature_state_mapping_draws(
    mapping: FeatureStateMapping, diagnostics: dict, n: int, seed: int,
) -> tuple[FeatureStateMapping, ...]:
    """Draw and validate assay-to-initial-state mapping coefficients."""
    _, samples = calibration_parameter_draws(mapping, diagnostics, n, seed)
    draws = tuple(replace(mapping, **sample) for sample in samples)
    for draw in draws:
        _validate_state_mapping(draw)
    return draws


def augment_with_liver_tdln_counterfactuals(
    features: pd.DataFrame, *, doses: tuple[Dose, ...], days: float = 28.0,
    step: float = 0.25, mapping: FeatureBridgeMapping | None = None,
    base_coupling: LiverTDLNCouplingParameters | None = None,
    base_bridge: LiverTDLNSpatialBridge | None = None,
    base_spatial: SpatialParameters | None = None,
    base_liver: LiverTumourParameters | None = None,
    base_immune: TwinParameters | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    parameter_draws: tuple[
        tuple[LiverTumourParameters, LiverTDLNCouplingParameters,
              LiverTDLNSpatialBridge], ...
    ] | None = None,
    mapping_draws: tuple[FeatureBridgeMapping, ...] | None = None,
    immune_draws: tuple[TwinParameters, ...] | None = None,
    spatial_draws: tuple[SpatialParameters, ...] | None = None,
    state_mapping_draws: tuple[FeatureStateMapping, ...] | None = None,
    patient_doses: dict[str, tuple[Dose, ...]] | None = None,
    patient_flows: dict[str, tuple[float, float]] | None = None,
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
    if patient_doses is not None and set(patient_doses) != set(out.patient_id):
        raise ValueError("patient treatment schedules must exactly match feature patients")
    if patient_flows is not None and set(patient_flows) != set(out.patient_id):
        raise ValueError("patient perfusion settings must exactly match feature patients")
    rows = []
    if parameter_draws is not None and len(parameter_draws) < 2:
        raise ValueError("parameter_draws must contain at least two samples")
    if mapping_draws is not None and len(mapping_draws) < 2:
        raise ValueError("mapping_draws must contain at least two samples")
    if immune_draws is not None and len(immune_draws) < 2:
        raise ValueError("immune_draws must contain at least two samples")
    if spatial_draws is not None and len(spatial_draws) < 2:
        raise ValueError("spatial_draws must contain at least two samples")
    if state_mapping_draws is not None and len(state_mapping_draws) < 2:
        raise ValueError("state_mapping_draws must contain at least two samples")
    draw_lengths = {len(draws) for draws in (
        parameter_draws, mapping_draws, immune_draws, spatial_draws,
        state_mapping_draws)
                    if draws is not None}
    if len(draw_lengths) > 1:
        raise ValueError("all uncertainty draws must have equal sample counts")
    for _, row in out.iterrows():
        schedule = patient_doses[str(row.patient_id)] if patient_doses is not None else doses
        portal_flow, arterial_flow = (
            patient_flows[str(row.patient_id)] if patient_flows is not None
            else (None, None))
        _, _, effect, _ = simulate_patient_liver_tdln_counterfactual(
            row, doses=schedule, days=days, step=step, mapping=mapping,
            base_coupling=base_coupling, base_bridge=base_bridge,
            base_spatial=base_spatial, base_liver=base_liver,
            base_immune=base_immune,
            feature_state_mapping=feature_state_mapping,
            portal_flow_multiplier=portal_flow,
            arterial_flow_multiplier=arterial_flow)
        values = {target: float(effect[source])
                  for source, target in COUNTERFACTUAL_FEATURES.items()}
        if not np.isfinite(list(values.values())).all():
            raise RuntimeError(f"non-finite counterfactual feature for patient {row.patient_id}")
        if any(draws is not None for draws in (
                parameter_draws, mapping_draws, immune_draws, spatial_draws,
                state_mapping_draws)):
            sampled = {source: [] for source in COUNTERFACTUAL_FEATURES}
            sample_count = next(iter(draw_lengths))
            for sample_index in range(sample_count):
                sampled_liver, sampled_coupling, sampled_bridge = (
                    parameter_draws[sample_index] if parameter_draws is not None
                    else (base_liver, base_coupling, base_bridge))
                sampled_mapping = (mapping_draws[sample_index]
                                   if mapping_draws is not None else mapping)
                sampled_immune = (immune_draws[sample_index]
                                  if immune_draws is not None else base_immune)
                sampled_spatial = (spatial_draws[sample_index]
                                   if spatial_draws is not None else base_spatial)
                sampled_state_mapping = (
                    state_mapping_draws[sample_index]
                    if state_mapping_draws is not None else feature_state_mapping)
                _, _, draw_effect, _ = simulate_patient_liver_tdln_counterfactual(
                    row, doses=schedule, days=days, step=step, mapping=sampled_mapping,
                    base_coupling=sampled_coupling, base_bridge=sampled_bridge,
                    base_spatial=sampled_spatial, base_liver=sampled_liver,
                    base_immune=sampled_immune,
                    feature_state_mapping=sampled_state_mapping,
                    portal_flow_multiplier=portal_flow,
                    arterial_flow_multiplier=arterial_flow)
                for source in sampled:
                    sampled[source].append(float(draw_effect[source]))
            for source, target in COUNTERFACTUAL_FEATURES.items():
                for quantile, suffix in UNCERTAINTY_QUANTILES:
                    values[f"{target}_{suffix}"] = float(np.quantile(
                        sampled[source], quantile))
            values["qc_liver_cf_uncertainty_samples"] = sample_count
        rows.append(values)
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def verify_liver_counterfactual_augmentation(
    audit_path: str | Path, features_path: str | Path, output_path: str | Path,
    bridge_parameters_path: str | Path | None = None,
    bridge_diagnostics_path: str | Path | None = None,
    treatment_manifest_path: str | Path | None = None,
    perfusion_manifest_path: str | Path | None = None,
    feature_bridge_mapping_path: str | Path | None = None,
    feature_bridge_diagnostics_path: str | Path | None = None,
    systemic_parameters_path: str | Path | None = None,
    systemic_diagnostics_path: str | Path | None = None,
    spatial_parameters_path: str | Path | None = None,
    spatial_diagnostics_path: str | Path | None = None,
    feature_state_mapping_path: str | Path | None = None,
    feature_state_diagnostics_path: str | Path | None = None,
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
    mapping_record = audit.get("feature_bridge_mapping")
    if mapping_record is not None:
        if feature_bridge_mapping_path is None:
            raise ValueError("feature-bridge mapping file is required to replay this augmentation")
        if mapping_record.get("sha256") != sha256_file(feature_bridge_mapping_path):
            raise ValueError("liver counterfactual feature-bridge mapping hash mismatch")
        mapping_payload = json.loads(Path(feature_bridge_mapping_path).read_text())
        if (mapping_payload.get("research_only") is not True
                or mapping_payload.get("labels_used") is not False
                or mapping_payload.get("mapping") != audit.get("mapping")):
            raise ValueError("liver counterfactual feature-bridge mapping contract mismatch")
    systemic_record = audit.get("systemic_parameters")
    if systemic_record is not None:
        if systemic_parameters_path is None:
            raise ValueError("systemic parameter file is required to replay this augmentation")
        if systemic_record.get("sha256") != sha256_file(systemic_parameters_path):
            raise ValueError("liver counterfactual systemic parameter hash mismatch")
    systemic_diagnostic_record = audit.get("systemic_diagnostics")
    if systemic_diagnostic_record is not None:
        if systemic_diagnostics_path is None:
            raise ValueError("systemic diagnostics are required to replay uncertainty")
        if systemic_diagnostic_record.get("sha256") != sha256_file(
                systemic_diagnostics_path):
            raise ValueError("liver counterfactual systemic diagnostics hash mismatch")
    spatial_record = audit.get("spatial_parameters")
    if spatial_record is not None:
        if spatial_parameters_path is None:
            raise ValueError("spatial parameter file is required to replay this augmentation")
        if spatial_record.get("sha256") != sha256_file(spatial_parameters_path):
            raise ValueError("liver counterfactual spatial parameter hash mismatch")
    spatial_diagnostic_record = audit.get("spatial_diagnostics")
    if spatial_diagnostic_record is not None:
        if spatial_diagnostics_path is None:
            raise ValueError("spatial diagnostics are required to replay uncertainty")
        if spatial_diagnostic_record.get("sha256") != sha256_file(spatial_diagnostics_path):
            raise ValueError("liver counterfactual spatial diagnostics hash mismatch")
    state_mapping_record = audit.get("feature_state_mapping")
    if state_mapping_record is not None:
        if feature_state_mapping_path is None:
            raise ValueError("feature-state mapping is required to replay this augmentation")
        if state_mapping_record.get("sha256") != sha256_file(feature_state_mapping_path):
            raise ValueError("liver counterfactual feature-state mapping hash mismatch")
    state_diagnostic_record = audit.get("feature_state_diagnostics")
    if state_diagnostic_record is not None:
        if feature_state_diagnostics_path is None:
            raise ValueError("feature-state diagnostics are required to replay uncertainty")
        if state_diagnostic_record.get("sha256") != sha256_file(
                feature_state_diagnostics_path):
            raise ValueError("liver counterfactual feature-state diagnostics hash mismatch")
    settings = audit.get("settings", {})
    try:
        doses = tuple(Dose(float(day), float(amount)) for day, amount in zip(
            settings["dose_days"], settings["dose_amounts"], strict=True))
        days, step = float(settings["days"]), float(settings["step"])
        mapping = FeatureBridgeMapping(**audit["mapping"])
        base_coupling = LiverTDLNCouplingParameters(**audit["base_parameters"]["coupling"])
        base_bridge = LiverTDLNSpatialBridge(**audit["base_parameters"]["spatial_bridge"])
        base_liver = LiverTumourParameters(**audit["base_parameters"].get(
            "liver_tumour", asdict(LiverTumourParameters())))
        base_immune = TwinParameters(**audit["base_parameters"].get(
            "immune", asdict(TwinParameters())))
        base_spatial = SpatialParameters(**audit["base_parameters"].get(
            "spatial", asdict(SpatialParameters())))
        feature_state_mapping = FeatureStateMapping(**audit["base_parameters"].get(
            "feature_state_mapping", asdict(FeatureStateMapping())))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid liver counterfactual replay settings") from error
    samples, seed = int(settings.get("uncertainty_samples", 0)), int(
        settings.get("uncertainty_seed", 42))
    parameter_draws = None
    diagnostic_record = audit.get("bridge_diagnostics")
    if samples and diagnostic_record is not None:
        if bridge_diagnostics_path is None:
            raise ValueError("bridge diagnostics are required to replay bridge uncertainty")
        if diagnostic_record.get("sha256") != sha256_file(bridge_diagnostics_path):
            raise ValueError("liver counterfactual bridge diagnostics hash mismatch")
        parameter_draws = liver_bridge_parameter_draws(
            base_liver, base_coupling, base_bridge,
            json.loads(Path(bridge_diagnostics_path).read_text()), samples, seed)
    elif diagnostic_record is not None:
        raise ValueError("bridge diagnostics supplied without uncertainty samples")
    mapping_draws = None
    mapping_diagnostic_record = audit.get("feature_bridge_diagnostics")
    if samples and mapping_diagnostic_record is not None:
        if feature_bridge_diagnostics_path is None or feature_bridge_mapping_path is None:
            raise ValueError("feature-bridge mapping and diagnostics are required for replay")
        if mapping_diagnostic_record.get("sha256") != sha256_file(
                feature_bridge_diagnostics_path):
            raise ValueError("liver counterfactual feature-bridge diagnostics hash mismatch")
        from .liver_mapping_assimilation import feature_bridge_mapping_draws
        mapping_draws = feature_bridge_mapping_draws(
            mapping, json.loads(Path(feature_bridge_diagnostics_path).read_text()),
            samples, seed + 1_000_003)
    elif mapping_diagnostic_record is not None:
        raise ValueError("feature-bridge diagnostics supplied without uncertainty samples")
    immune_draws = None
    if samples and systemic_diagnostic_record is not None:
        immune_draws = systemic_parameter_draws(
            base_immune, json.loads(Path(systemic_diagnostics_path).read_text()),
            samples, seed + 2_000_003)
    elif systemic_diagnostic_record is not None:
        raise ValueError("systemic diagnostics supplied without uncertainty samples")
    spatial_draws = None
    if samples and spatial_diagnostic_record is not None:
        spatial_draws = spatial_parameter_draws(
            base_spatial, json.loads(Path(spatial_diagnostics_path).read_text()),
            samples, seed + 3_000_017)
    elif spatial_diagnostic_record is not None:
        raise ValueError("spatial diagnostics supplied without uncertainty samples")
    state_mapping_draws = None
    if samples and state_diagnostic_record is not None:
        state_mapping_draws = feature_state_mapping_draws(
            feature_state_mapping,
            json.loads(Path(feature_state_diagnostics_path).read_text()),
            samples, seed + 4_000_037)
    elif state_diagnostic_record is not None:
        raise ValueError("feature-state diagnostics supplied without uncertainty samples")
    if (samples and parameter_draws is None and mapping_draws is None
            and immune_draws is None and spatial_draws is None
            and state_mapping_draws is None):
        raise ValueError("uncertainty samples require at least one diagnostics file")
    patient_doses = None
    treatment_record = audit.get("treatment_manifest")
    if treatment_record is not None:
        if treatment_manifest_path is None:
            raise ValueError("treatment manifest is required to replay this augmentation")
        if treatment_record.get("sha256") != sha256_file(treatment_manifest_path):
            raise ValueError("liver counterfactual treatment manifest hash mismatch")
        patient_doses, manifest_summary = validate_treatment_manifest(
            pd.read_csv(treatment_manifest_path),
            pd.read_csv(features_path).patient_id.astype(str).tolist(), days)
        if manifest_summary != audit.get("treatment_contract"):
            raise ValueError("treatment manifest contract mismatch")
    patient_flows = None
    perfusion_record = audit.get("perfusion_manifest")
    if perfusion_record is not None:
        if perfusion_manifest_path is None:
            raise ValueError("perfusion manifest is required to replay this augmentation")
        if perfusion_record.get("sha256") != sha256_file(perfusion_manifest_path):
            raise ValueError("liver counterfactual perfusion manifest hash mismatch")
        patient_flows, perfusion_summary = validate_perfusion_manifest(
            pd.read_csv(perfusion_manifest_path),
            pd.read_csv(features_path).patient_id.astype(str).tolist())
        if perfusion_summary != audit.get("perfusion_contract"):
            raise ValueError("perfusion manifest contract mismatch")
    replay = augment_with_liver_tdln_counterfactuals(
        pd.read_csv(features_path), doses=doses, days=days, step=step, mapping=mapping,
        base_coupling=base_coupling, base_bridge=base_bridge,
        base_spatial=base_spatial,
        feature_state_mapping=feature_state_mapping,
        base_liver=base_liver,
        base_immune=base_immune,
        parameter_draws=parameter_draws, mapping_draws=mapping_draws,
        immune_draws=immune_draws,
        spatial_draws=spatial_draws,
        state_mapping_draws=state_mapping_draws,
        patient_doses=patient_doses,
        patient_flows=patient_flows)
    replay_sha = hashlib.sha256(
        replay.to_csv(index=False).encode("utf-8")).hexdigest()
    if replay_sha != audit["output"]["sha256"]:
        raise ValueError("liver counterfactual augmentation derivation replay mismatch")
    return audit


def load_liver_bridge_parameters(
    path: str | Path,
) -> tuple[LiverTumourParameters, LiverTDLNCouplingParameters,
           LiverTDLNSpatialBridge]:
    """Load and fail-closed validate a saved liver–TDLN calibration."""
    payload = json.loads(Path(path).read_text())
    if payload.get("research_only") is not True:
        raise ValueError("liver bridge parameter file must declare research_only=true")
    try:
        liver = LiverTumourParameters(**payload.get("liver_tumour", {}))
        coupling = LiverTDLNCouplingParameters(**payload["coupling"])
        bridge = LiverTDLNSpatialBridge(**payload["spatial_bridge"])
    except (KeyError, TypeError) as error:
        raise ValueError("invalid liver bridge parameter schema") from error
    _validate_coupling(coupling)
    _validate_bridge(bridge)
    _validate_liver_parameters(liver)
    return liver, coupling, bridge


def load_systemic_parameters(path: str | Path) -> TwinParameters:
    """Load an exact calibrated systemic parameter artifact for liver coupling."""
    payload = json.loads(Path(path).read_text())
    expected = {field.name for field in fields(TwinParameters)}
    if not isinstance(payload, dict) or set(payload) != expected:
        raise ValueError("invalid systemic parameter schema for liver counterfactuals")
    try:
        parameters = TwinParameters(**{name: float(payload[name]) for name in expected})
    except (TypeError, ValueError) as error:
        raise ValueError("systemic parameters must be numeric") from error
    _validate_model_inputs(parameters, InitialState())
    return parameters


def load_spatial_parameters(path: str | Path) -> SpatialParameters:
    """Load an exact calibrated six-zone parameter artifact for liver coupling."""
    payload = json.loads(Path(path).read_text())
    expected = {field.name for field in fields(SpatialParameters)}
    if not isinstance(payload, dict) or set(payload) != expected:
        raise ValueError("invalid spatial parameter schema for liver counterfactuals")
    try:
        parameters = SpatialParameters(**{name: float(payload[name]) for name in expected})
    except (TypeError, ValueError) as error:
        raise ValueError("spatial parameters must be numeric") from error
    _validate_spatial_parameters(parameters)
    return parameters


def load_feature_state_mapping(path: str | Path) -> FeatureStateMapping:
    """Load an exact calibrated assay-to-state mapping for liver coupling."""
    payload = json.loads(Path(path).read_text())
    expected = {field.name for field in fields(FeatureStateMapping)}
    if not isinstance(payload, dict) or set(payload) != expected:
        raise ValueError("invalid feature-state mapping schema for liver counterfactuals")
    try:
        mapping = FeatureStateMapping(**{name: float(payload[name]) for name in expected})
    except (TypeError, ValueError) as error:
        raise ValueError("feature-state mapping coefficients must be numeric") from error
    _validate_state_mapping(mapping)
    return mapping


def validate_treatment_manifest(
    manifest: pd.DataFrame, patient_ids: list[str], days: float,
) -> tuple[dict[str, tuple[Dose, ...]], dict]:
    """Validate patient-specific anti-PD-1 schedules and outcome semantics."""
    required = {"patient_id", "regimen_id", "agent_class", "dose_days",
                "model_dose_amount", "dose_scale", "response_definition"}
    missing = required - set(manifest)
    if missing:
        raise ValueError(f"treatment manifest missing columns: {sorted(missing)}")
    table = manifest.copy()
    if table[list(required)].isna().any().any():
        raise ValueError("treatment manifest required fields must be complete")
    text_fields = ("patient_id", "regimen_id", "agent_class", "dose_days",
                   "dose_scale", "response_definition")
    if any(table[field].astype(str).str.strip().eq("").any() for field in text_fields):
        raise ValueError("treatment manifest required text fields must be non-empty")
    table["patient_id"] = table.patient_id.astype(str).str.strip()
    expected = {str(value).strip() for value in patient_ids}
    if (table.patient_id.eq("").any() or table.patient_id.duplicated().any()
            or set(table.patient_id) != expected):
        raise ValueError("treatment manifest must exactly match unique feature patients")
    agent_classes = set(table.agent_class.astype(str).str.strip().str.lower())
    if agent_classes != {"anti_pd1"}:
        raise ValueError("treatment manifest supports only agent_class=anti_pd1")
    scales = set(table.dose_scale.astype(str).str.strip())
    if scales != {"arbitrary_model_exposure"}:
        raise ValueError("treatment manifest dose_scale must be arbitrary_model_exposure")
    response_definitions = set(table.response_definition.astype(str).str.strip())
    if len(response_definitions) != 1 or "" in response_definitions:
        raise ValueError("treatment manifest requires one shared response_definition")
    schedules = {}
    for row in table.itertuples(index=False):
        try:
            amount = float(row.model_dose_amount)
            dose_days = tuple(float(value.strip()) for value in str(row.dose_days).split(",")
                              if value.strip())
        except (TypeError, ValueError) as error:
            raise ValueError("treatment manifest doses must be numeric") from error
        if (not np.isfinite(amount) or amount <= 0 or not dose_days
                or len(set(dose_days)) != len(dose_days)
                or any(not np.isfinite(day) or day < 0 or day >= days for day in dose_days)):
            raise ValueError("treatment manifest has invalid positive in-window schedule")
        schedules[row.patient_id] = tuple(Dose(day, amount) for day in sorted(dose_days))
    summary = {
        "agent_class": "anti_pd1",
        "dose_scale": "arbitrary_model_exposure",
        "response_definition": next(iter(response_definitions)),
        "regimen_ids": sorted(set(table.regimen_id.astype(str).str.strip())),
        "patient_specific_schedules": True,
        "schedule_support": {
            "model_dose_amount": [
                float(min(schedule[0].amount for schedule in schedules.values())),
                float(max(schedule[0].amount for schedule in schedules.values())),
            ],
            "dose_count": [
                int(min(len(schedule) for schedule in schedules.values())),
                int(max(len(schedule) for schedule in schedules.values())),
            ],
            "total_model_exposure": [
                float(min(sum(dose.amount for dose in schedule)
                          for schedule in schedules.values())),
                float(max(sum(dose.amount for dose in schedule)
                          for schedule in schedules.values())),
            ],
            "first_dose_day": [
                float(min(schedule[0].day for schedule in schedules.values())),
                float(max(schedule[0].day for schedule in schedules.values())),
            ],
            "last_dose_day": [
                float(min(schedule[-1].day for schedule in schedules.values())),
                float(max(schedule[-1].day for schedule in schedules.values())),
            ],
        },
    }
    return schedules, summary


def validate_perfusion_manifest(
    manifest: pd.DataFrame, patient_ids: list[str],
) -> tuple[dict[str, tuple[float, float]], dict]:
    """Validate patient-aligned dimensionless dual-inflow measurements."""
    required = {"patient_id", "portal_flow_multiplier", "arterial_flow_multiplier",
                "measurement_definition", "normalization_reference"}
    missing = required - set(manifest)
    if missing:
        raise ValueError(f"perfusion manifest missing columns: {sorted(missing)}")
    table = manifest.copy()
    if table[list(required)].isna().any().any():
        raise ValueError("perfusion manifest requires complete values")
    table["patient_id"] = table.patient_id.astype(str).str.strip()
    if table.patient_id.eq("").any() or table.patient_id.duplicated().any():
        raise ValueError("perfusion manifest requires unique non-empty patient_id")
    if set(table.patient_id) != {str(value).strip() for value in patient_ids}:
        raise ValueError("perfusion manifest patients must exactly match feature patients")
    definitions = set(table.measurement_definition.astype(str).str.strip())
    references = set(table.normalization_reference.astype(str).str.strip())
    if len(definitions) != 1 or "" in definitions:
        raise ValueError("perfusion manifest requires one shared measurement_definition")
    if len(references) != 1 or "" in references:
        raise ValueError("perfusion manifest requires one shared normalization_reference")
    flows = {}
    for row in table.itertuples(index=False):
        try:
            portal = float(row.portal_flow_multiplier)
            arterial = float(row.arterial_flow_multiplier)
        except (TypeError, ValueError) as error:
            raise ValueError("perfusion multipliers must be numeric") from error
        if not np.isfinite((portal, arterial)).all() or portal <= 0 or arterial <= 0:
            raise ValueError("perfusion multipliers must be finite and strictly positive")
        flows[str(row.patient_id)] = (portal, arterial)
    return flows, {
        "scale": "dimensionless_multiplier",
        "measurement_definition": next(iter(definitions)),
        "normalization_reference": next(iter(references)),
        "patient_count": len(flows),
    }

