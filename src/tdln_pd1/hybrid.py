"""Fuse patient-level single-cell features with mechanistic counterfactuals."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import re
import numpy as np
import pandas as pd

from .simulator import Dose, InitialState, TwinParameters, simulate
from .spatial_dynamics import (SpatialInitialState, SpatialParameters,
                               SPECIES, ZONE_NAMES, simulate_spatial_coupled)


OUTCOME_TOKENS = (
    "response", "responder", "outcome", "label", "recist",
    "progression_free", "overall_survival", "pfs", "os_months",
)


@dataclass(frozen=True)
class FeatureStateMapping:
    """Versionable coefficients connecting assay summaries to model states."""

    signature_intercept: float = 0.5
    signature_lower_multiplier: float = 0.25
    immune_upper_multiplier: float = 4.0
    effector_upper_multiplier: float = 5.0
    cytotoxic_coefficient: float = 1.0
    cd8_fraction_coefficient: float = 1.0
    exhaustion_coefficient: float = 1.0
    tcf7_coefficient: float = 1.0
    dendritic_coefficient: float = 1.0
    treg_fraction_coefficient: float = 4.0
    antigen_proxy_coefficient: float = 1.0


def _validate_mapping(mapping: FeatureStateMapping) -> None:
    values = asdict(mapping)
    invalid = [name for name, value in values.items()
               if not np.isfinite(value) or value <= 0]
    if invalid:
        raise ValueError(
            "feature-state mapping coefficients must be finite and strictly positive: "
            f"{sorted(invalid)}"
        )
    if mapping.signature_lower_multiplier >= min(
        mapping.immune_upper_multiplier, mapping.effector_upper_multiplier
    ):
        raise ValueError("feature-state mapping lower multiplier must be below upper bounds")


def _outcome_like_columns(frame: pd.DataFrame) -> list[str]:
    """Return columns that could carry treatment outcome into the twin transform."""
    return sorted([
        column for column in frame.columns
        if any(
            token in re.sub(r"[^a-z0-9]+", "_", str(column).lower())
            for token in OUTCOME_TOKENS
        )
    ])


def _finite(row: pd.Series, name: str, default: float) -> float:
    if name not in row or pd.isna(row[name]):
        return default
    value = row[name]
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"patient feature must be numeric when supplied: {name}") from error
    if not np.isfinite(value):
        raise ValueError(f"patient feature must be finite when supplied: {name}")
    return value


def _finite_alias(row: pd.Series, names: tuple[str, ...], default: float) -> float:
    """Prefer a lineage-restricted state estimate, then a legacy global one."""
    for name in names:
        if name in row and pd.notna(row[name]):
            return _finite(row, name, default)
    return default


def _nonnegative(row: pd.Series, name: str, default: float) -> float:
    value = _finite(row, name, default)
    if value < 0:
        raise ValueError(f"patient feature must be non-negative: {name}")
    return value


def _fraction(row: pd.Series, names: tuple[str, ...], default: float) -> float:
    """Read the first available canonicalized cell-fraction alias."""
    for name in names:
        if name in row and pd.notna(row[name]):
            value = _finite(row, name, default)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"cell-fraction feature must fall within [0, 1]: {name}")
            return value
    return default


def initial_state_from_features(
    row: pd.Series,
    base: InitialState | None = None,
    mapping: FeatureStateMapping | None = None,
) -> InitialState:
    """Map measured immune features to initial states with explicit bounded rules.

    Signature inputs are finite signed summaries so centered/scaled assays remain
    admissible; fractions are in [0,1]. Bounded multipliers keep signed deviations
    from producing negative mechanistic states. The mapping is intentionally simple
    and inspectable so a calibrated mapping can replace it without changing the simulator.
    """
    base = base or InitialState()
    mapping = mapping or FeatureStateMapping()
    _validate_mapping(mapping)
    cytotoxic = _finite_alias(row, ("sig_cytotoxic_cd8_t", "sig_cytotoxic"), 0.0)
    exhaustion = _finite_alias(row, ("sig_exhaustion_cd8_t", "sig_exhaustion"), 0.0)
    tcf7 = _finite_alias(row, ("sig_tcf7_memory_cd8_t", "sig_tcf7_memory"), 0.0)
    dendritic = _finite_alias(
        row, ("sig_dendritic_activation_dc", "sig_dendritic_activation"), 0.0
    )
    treg_fraction = _fraction(row, ("frac_treg", "frac_regulatory_t_cell"), 0.05)
    cd8_fraction = _fraction(row, ("frac_cd8_t", "frac_cd8_t_cell"), 0.15)
    return replace(
        base,
        mature_dc_tdln=base.mature_dc_tdln * np.clip(
            mapping.signature_intercept + mapping.dendritic_coefficient * dendritic,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        ),
        tpex_tdln=base.tpex_tdln * np.clip(
            mapping.signature_intercept + mapping.tcf7_coefficient * tcf7,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        ),
        effector_tumour=base.effector_tumour * np.clip(
            mapping.signature_intercept + mapping.cytotoxic_coefficient * cytotoxic
            + mapping.cd8_fraction_coefficient * cd8_fraction,
            mapping.signature_lower_multiplier, mapping.effector_upper_multiplier,
        ),
        exhausted_tumour=base.exhausted_tumour * np.clip(
            mapping.signature_intercept + mapping.exhaustion_coefficient * exhaustion,
            mapping.signature_lower_multiplier, mapping.effector_upper_multiplier,
        ),
        treg_tdln=base.treg_tdln * np.clip(
            mapping.signature_intercept + mapping.treg_fraction_coefficient * treg_fraction,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        ),
    )


def spatial_initial_state_from_features(
    row: pd.Series,
    base: SpatialInitialState | None = None,
    mapping: FeatureStateMapping | None = None,
) -> SpatialInitialState:
    """Map the same measured immune features to bounded six-zone inputs.

    No outcome label enters this mapping. ``tumour_antigen_proxy`` is optional
    because transcriptomes do not directly measure lymph-borne antigen mass; if
    absent, the auditable population default is retained.
    """
    base = base or SpatialInitialState()
    mapping = mapping or FeatureStateMapping()
    _validate_mapping(mapping)
    dendritic = _finite_alias(
        row, ("sig_dendritic_activation_dc", "sig_dendritic_activation"), 0.0
    )
    tcf7 = _finite_alias(row, ("sig_tcf7_memory_cd8_t", "sig_tcf7_memory"), 0.0)
    cytotoxic = _finite_alias(row, ("sig_cytotoxic_cd8_t", "sig_cytotoxic"), 0.0)
    cd8_fraction = _fraction(row, ("frac_cd8_t", "frac_cd8_t_cell"), 0.15)
    antigen_proxy = _nonnegative(row, "tumour_antigen_proxy", 0.5)
    return replace(
        base,
        afferent_antigen=float(base.afferent_antigen * np.clip(
            mapping.signature_intercept + mapping.antigen_proxy_coefficient * antigen_proxy,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        )),
        afferent_dc=float(base.afferent_dc * np.clip(
            mapping.signature_intercept + mapping.dendritic_coefficient * dendritic,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        )),
        paracortex_tpex=float(base.paracortex_tpex * np.clip(
            mapping.signature_intercept + mapping.tcf7_coefficient * tcf7,
            mapping.signature_lower_multiplier, mapping.immune_upper_multiplier,
        )),
        paracortex_effector=float(
            base.paracortex_effector * np.clip(
                mapping.signature_intercept + mapping.cytotoxic_coefficient * cytotoxic
                + mapping.cd8_fraction_coefficient * cd8_fraction,
                mapping.signature_lower_multiplier, mapping.effector_upper_multiplier,
            )
        ),
    )


def augment_with_twin(
    features: pd.DataFrame,
    days: float = 42.0,
    dose_days: tuple[float, ...] = (0.0, 14.0, 28.0),
    dose: float = 1.0,
    spatial_states: pd.DataFrame | None = None,
    systemic_parameters: TwinParameters | None = None,
    spatial_parameters: SpatialParameters | None = None,
    systemic_parameter_draws: tuple[TwinParameters, ...] | None = None,
    spatial_parameter_draws: tuple[SpatialParameters, ...] | None = None,
    feature_state_mapping: FeatureStateMapping | None = None,
    feature_state_mapping_draws: tuple[FeatureStateMapping, ...] | None = None,
    spatial_missing_zone_policy: str = "feature_prior",
) -> pd.DataFrame:
    """Add label-free point and optional calibrated-uncertainty endpoints."""
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    features = features.copy()
    mapping = feature_state_mapping or FeatureStateMapping()
    _validate_mapping(mapping)
    if spatial_missing_zone_policy not in {"feature_prior", "explicit_zero"}:
        raise ValueError(
            "spatial_missing_zone_policy must be feature_prior or explicit_zero"
        )
    forbidden_features = _outcome_like_columns(features)
    if forbidden_features:
        raise ValueError(
            f"patient features must be outcome-free before twin augmentation: {forbidden_features}"
        )
    if features.patient_id.isna().any():
        raise ValueError("features need non-missing patient_id")
    features["patient_id"] = features.patient_id.astype(str).str.strip()
    if features.patient_id.eq("").any() or features.patient_id.duplicated().any():
        raise ValueError("features must have one row per patient")
    for name, draws in (
        ("systemic_parameter_draws", systemic_parameter_draws),
        ("spatial_parameter_draws", spatial_parameter_draws),
        ("feature_state_mapping_draws", feature_state_mapping_draws),
    ):
        if draws is not None and len(draws) == 0:
            raise ValueError(f"{name} must be non-empty when supplied")
    supplied_draw_counts = {
        len(draws) for draws in (
            systemic_parameter_draws, spatial_parameter_draws, feature_state_mapping_draws,
        ) if draws is not None
    }
    if len(supplied_draw_counts) > 1:
        raise ValueError("systemic, spatial and mapping parameter draw counts must match")
    uncertainty_samples = next(iter(supplied_draw_counts), 0)
    measured_states = None
    if spatial_states is not None:
        if "patient_id" not in spatial_states:
            raise ValueError("spatial initial states need patient_id")
        forbidden = _outcome_like_columns(spatial_states)
        if forbidden:
            raise ValueError(f"spatial initial states must be outcome-free: {sorted(forbidden)}")
        measured_states = spatial_states.copy()
        if measured_states.patient_id.isna().any():
            raise ValueError("spatial initial states need non-missing patient_id")
        measured_states["patient_id"] = measured_states.patient_id.astype(str).str.strip()
        if measured_states.patient_id.eq("").any():
            raise ValueError("spatial initial states need non-empty patient_id")
        feature_ids, spatial_ids = set(features.patient_id), set(measured_states.patient_id)
        if feature_ids != spatial_ids:
            raise ValueError(
                "spatial and feature patient sets differ: "
                f"missing_spatial={sorted(feature_ids - spatial_ids)[:10]}, "
                f"extra_spatial={sorted(spatial_ids - feature_ids)[:10]}"
            )
    spatial_fusion_contract = {
        "schema_version": 1,
        "spatial_initialization_source": (
            "measured_spatial_with_patient_feature_prior"
            if measured_states is not None else "patient_feature_prior_only"
        ),
        "missing_zone_policy": spatial_missing_zone_policy,
        "observed_zone_rule": "retain_measured_amounts",
        "unobserved_zone_rule": (
            "patient_feature_state_mapping"
            if spatial_missing_zone_policy == "feature_prior" else "explicit_zero"
        ),
    }
    spatial_fusion_contract_sha256 = hashlib.sha256(
        json.dumps(
            spatial_fusion_contract, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    doses = tuple(Dose(day, dose) for day in dose_days)
    simulation_contract = {
        "schema_version": 1,
        "systemic_endpoint_day": float(days),
        "spatial_endpoint_day": float(days),
        "spatial_horizon_rule": "same_as_systemic_endpoint",
        "dose_days": [float(day) for day in dose_days],
        "dose_amount": float(dose),
        "counterfactual_design": "paired_same_patient_initial_state",
    }
    mechanistic_simulation_contract_sha256 = hashlib.sha256(
        json.dumps(
            simulation_contract, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    rows = []
    for row in features.to_dict(orient="records"):
        measured_spatial_initial = (
            measured_states[measured_states.patient_id.eq(row["patient_id"])]
            if measured_states is not None
            else None
        )
        if measured_states is not None and "zone_observed" in measured_spatial_initial:
            zone_audit = measured_spatial_initial[["zone", "zone_observed"]].drop_duplicates()
            if zone_audit.zone.duplicated().any():
                raise ValueError("zone_observed must be consistent across species within each zone")
            if pd.api.types.is_bool_dtype(zone_audit.zone_observed):
                observed_flag = zone_audit.zone_observed
            else:
                normalized = zone_audit.zone_observed.astype(str).str.strip().str.lower()
                if not normalized.isin({"true", "false"}).all():
                    raise ValueError("zone_observed must contain only true/false values")
                observed_flag = normalized.eq("true")
            spatial_zone_coverage = float(observed_flag.mean())
            observed_zones = set(zone_audit.loc[observed_flag, "zone"])
        else:
            spatial_zone_coverage = np.nan
            observed_zones = set(ZONE_NAMES) if measured_states is not None else set()
        if measured_states is not None and "sampling_unit_coverage" in measured_spatial_initial:
            sampling_audit = measured_spatial_initial[
                ["zone", "sampling_unit_coverage"]
            ].drop_duplicates()
            if sampling_audit.zone.duplicated().any():
                raise ValueError(
                    "sampling_unit_coverage must be consistent across species within each zone"
                )
            sampling_audit["sampling_unit_coverage"] = pd.to_numeric(
                sampling_audit.sampling_unit_coverage, errors="coerce"
            )
            if (not np.isfinite(sampling_audit.sampling_unit_coverage).all()
                    or not sampling_audit.sampling_unit_coverage.between(0, 1).all()):
                raise ValueError("sampling_unit_coverage must be finite and between zero and one")
            spatial_sampling_coverage_mean = float(
                sampling_audit.sampling_unit_coverage.mean()
            )
            spatial_sampling_coverage_min = float(
                sampling_audit.sampling_unit_coverage.min()
            )
        else:
            spatial_sampling_coverage_mean = np.nan
            spatial_sampling_coverage_min = np.nan
        def endpoint(frame: pd.DataFrame, species: str, zone: str) -> float:
            selected = frame[(frame.species == species) & (frame.zone == zone)]
            return float(selected.amount.iloc[-1])

        def run_pair(
            systemic: TwinParameters | None,
            spatial: SpatialParameters | None,
            state_mapping: FeatureStateMapping,
        ) -> dict[str, float]:
            initial = initial_state_from_features(pd.Series(row), mapping=state_mapping)
            feature_spatial_initial = spatial_initial_state_from_features(
                pd.Series(row), mapping=state_mapping
            )
            if measured_spatial_initial is None:
                spatial_initial = feature_spatial_initial
            elif (spatial_missing_zone_policy == "feature_prior"
                  and observed_zones != set(ZONE_NAMES)):
                prior_amounts = {
                    (species, zone): 0.0
                    for species in SPECIES for zone in ZONE_NAMES
                }
                prior_amounts[("antigen", "afferent_lymph")] = (
                    feature_spatial_initial.afferent_antigen
                )
                prior_amounts[("mature_dc", "afferent_lymph")] = (
                    feature_spatial_initial.afferent_dc
                )
                prior_amounts[("tpex", "paracortex")] = (
                    feature_spatial_initial.paracortex_tpex
                )
                prior_amounts[("effector", "paracortex")] = (
                    feature_spatial_initial.paracortex_effector
                )
                spatial_initial = measured_spatial_initial.copy()
                unobserved = ~spatial_initial.zone.isin(observed_zones)
                spatial_initial.loc[unobserved, "amount"] = [
                    prior_amounts[(species, zone)]
                    for species, zone in zip(
                        spatial_initial.loc[unobserved, "species"],
                        spatial_initial.loc[unobserved, "zone"],
                    )
                ]
            else:
                spatial_initial = measured_spatial_initial
            control = simulate(days=days, initial=initial, parameters=systemic)
            treated = simulate(days=days, doses=doses, initial=initial, parameters=systemic)
            spatial_control = simulate_spatial_coupled(
                control, days=days, initial=spatial_initial, parameters=spatial,
            )
            spatial_treated = simulate_spatial_coupled(
                treated, days=days, initial=spatial_initial, parameters=spatial,
            )
            c, t = control.iloc[-1], treated.iloc[-1]
            control_exit = endpoint(spatial_control, "effector", "efferent_lymph")
            treated_exit = endpoint(spatial_treated, "effector", "efferent_lymph")
            return {
                "twin_control_tumour_cells": c.tumour_cells,
                "twin_treated_tumour_cells": t.tumour_cells,
                "twin_treatment_tumour_delta": t.tumour_cells - c.tumour_cells,
                "twin_treated_effector_tumour": t.effector_tumour,
                "twin_treated_tpex_tdln": t.tpex_tdln,
                "twin_treated_tpex_effector_ratio": t.tpex_effector_ratio,
                "twin_peak_pd1_occupancy": treated.pd1_occupancy.max(),
                "twin_effector_gain": t.effector_tumour - c.effector_tumour,
                "twin_spatial_control_efferent_effector": control_exit,
                "twin_spatial_treated_efferent_effector": treated_exit,
                "twin_spatial_efferent_effector_gain": treated_exit - control_exit,
                "twin_spatial_treated_paracortex_tpex": endpoint(
                    spatial_treated, "tpex", "paracortex"
                ),
            }

        point = run_pair(systemic_parameters, spatial_parameters, mapping)
        result = {
            "patient_id": row["patient_id"],
            **point,
            "qc_spatial_initialization_measured": int(measured_states is not None),
            "qc_spatial_fusion_contract_sha256": spatial_fusion_contract_sha256,
            "qc_mechanistic_simulation_contract_sha256": (
                mechanistic_simulation_contract_sha256
            ),
            "qc_systemic_endpoint_day": float(days),
            "qc_spatial_endpoint_day": float(days),
            "qc_spatial_zone_coverage": spatial_zone_coverage,
            "qc_spatial_feature_prior_imputed_fraction": (
                1.0 - spatial_zone_coverage
                if measured_states is not None
                and spatial_missing_zone_policy == "feature_prior"
                and np.isfinite(spatial_zone_coverage)
                else 0.0
            ),
            "qc_spatial_sampling_coverage_mean": spatial_sampling_coverage_mean,
            "qc_spatial_sampling_coverage_min": spatial_sampling_coverage_min,
            "qc_systemic_parameters_supplied": int(systemic_parameters is not None),
            "qc_spatial_parameters_supplied": int(spatial_parameters is not None),
            "qc_mechanistic_uncertainty_samples": uncertainty_samples,
            "qc_systemic_uncertainty_supplied": int(systemic_parameter_draws is not None),
            "qc_spatial_uncertainty_supplied": int(spatial_parameter_draws is not None),
            "qc_feature_state_mapping_supplied": int(feature_state_mapping is not None),
            "qc_mapping_uncertainty_supplied": int(feature_state_mapping_draws is not None),
        }
        if uncertainty_samples:
            systemic_draws = systemic_parameter_draws or (systemic_parameters,) * uncertainty_samples
            spatial_draws = spatial_parameter_draws or (spatial_parameters,) * uncertainty_samples
            mapping_draws = feature_state_mapping_draws or (mapping,) * uncertainty_samples
            ensemble = [run_pair(systemic, spatial, state_mapping)
                        for systemic, spatial, state_mapping
                        in zip(systemic_draws, spatial_draws, mapping_draws)]
            for name in (
                "twin_treatment_tumour_delta", "twin_effector_gain",
                "twin_spatial_efferent_effector_gain", "twin_treated_tpex_tdln",
            ):
                values = np.asarray([draw[name] for draw in ensemble], dtype=float)
                lower, median, upper = np.quantile(values, (0.025, 0.5, 0.975))
                result.update({
                    f"{name}_uncertainty_median": median,
                    f"{name}_uncertainty_lower_95": lower,
                    f"{name}_uncertainty_upper_95": upper,
                    f"{name}_uncertainty_width_95": upper - lower,
                })
        rows.append(result)
    return features.merge(pd.DataFrame(rows), on="patient_id", validate="one_to_one")
