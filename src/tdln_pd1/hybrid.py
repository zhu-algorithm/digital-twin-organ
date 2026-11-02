"""Fuse patient-level single-cell features with mechanistic counterfactuals."""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pandas as pd

from .simulator import Dose, InitialState, TwinParameters, simulate
from .spatial_dynamics import (SpatialInitialState, SpatialParameters,
                               simulate_spatial_coupled)


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


def initial_state_from_features(row: pd.Series, base: InitialState | None = None) -> InitialState:
    """Map measured immune features to initial states with explicit bounded rules.

    Signature inputs are finite signed summaries so centered/scaled assays remain
    admissible; fractions are in [0,1]. Bounded multipliers keep signed deviations
    from producing negative mechanistic states. The mapping is intentionally simple
    and inspectable so a calibrated mapping can replace it without changing the simulator.
    """
    base = base or InitialState()
    cytotoxic = _finite(row, "sig_cytotoxic", 0.0)
    exhaustion = _finite(row, "sig_exhaustion", 0.0)
    tcf7 = _finite(row, "sig_tcf7_memory", 0.0)
    dendritic = _finite(row, "sig_dendritic_activation", 0.0)
    treg_fraction = _fraction(row, ("frac_treg", "frac_regulatory_t_cell"), 0.05)
    cd8_fraction = _fraction(row, ("frac_cd8_t", "frac_cd8_t_cell"), 0.15)
    return replace(
        base,
        mature_dc_tdln=base.mature_dc_tdln * np.clip(0.5 + dendritic, 0.25, 4.0),
        tpex_tdln=base.tpex_tdln * np.clip(0.5 + tcf7, 0.25, 4.0),
        effector_tumour=base.effector_tumour * np.clip(0.5 + cytotoxic + cd8_fraction, 0.25, 5.0),
        exhausted_tumour=base.exhausted_tumour * np.clip(0.5 + exhaustion, 0.25, 5.0),
        treg_tdln=base.treg_tdln * np.clip(0.5 + 4.0 * treg_fraction, 0.25, 4.0),
    )


def spatial_initial_state_from_features(
    row: pd.Series,
    base: SpatialInitialState | None = None,
) -> SpatialInitialState:
    """Map the same measured immune features to bounded six-zone inputs.

    No outcome label enters this mapping. ``tumour_antigen_proxy`` is optional
    because transcriptomes do not directly measure lymph-borne antigen mass; if
    absent, the auditable population default is retained.
    """
    base = base or SpatialInitialState()
    dendritic = _finite(row, "sig_dendritic_activation", 0.0)
    tcf7 = _finite(row, "sig_tcf7_memory", 0.0)
    cytotoxic = _finite(row, "sig_cytotoxic", 0.0)
    cd8_fraction = _fraction(row, ("frac_cd8_t", "frac_cd8_t_cell"), 0.15)
    antigen_proxy = _nonnegative(row, "tumour_antigen_proxy", 0.5)
    return replace(
        base,
        afferent_antigen=float(base.afferent_antigen * np.clip(0.5 + antigen_proxy, 0.25, 4.0)),
        afferent_dc=float(base.afferent_dc * np.clip(0.5 + dendritic, 0.25, 4.0)),
        paracortex_tpex=float(base.paracortex_tpex * np.clip(0.5 + tcf7, 0.25, 4.0)),
        paracortex_effector=float(
            base.paracortex_effector * np.clip(0.5 + cytotoxic + cd8_fraction, 0.25, 5.0)
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
) -> pd.DataFrame:
    """Add per-patient treated/control mechanistic endpoints without using labels."""
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    features = features.copy()
    if features.patient_id.isna().any():
        raise ValueError("features need non-missing patient_id")
    features["patient_id"] = features.patient_id.astype(str).str.strip()
    if features.patient_id.eq("").any() or features.patient_id.duplicated().any():
        raise ValueError("features must have one row per patient")
    measured_states = None
    if spatial_states is not None:
        if "patient_id" not in spatial_states:
            raise ValueError("spatial initial states need patient_id")
        forbidden = [column for column in spatial_states.columns
                     if any(token in column.lower()
                            for token in ("response", "outcome", "label", "recist"))]
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
    doses = tuple(Dose(day, dose) for day in dose_days)
    rows = []
    for row in features.to_dict(orient="records"):
        initial = initial_state_from_features(pd.Series(row))
        control = simulate(days=days, initial=initial, parameters=systemic_parameters)
        treated = simulate(days=days, doses=doses, initial=initial,
                           parameters=systemic_parameters)
        spatial_initial = (
            measured_states[measured_states.patient_id.eq(row["patient_id"])]
            if measured_states is not None
            else spatial_initial_state_from_features(pd.Series(row))
        )
        if measured_states is not None and "zone_observed" in spatial_initial:
            zone_audit = spatial_initial[["zone", "zone_observed"]].drop_duplicates()
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
        else:
            spatial_zone_coverage = np.nan
        spatial_days = min(days, 14.0)
        spatial_control = simulate_spatial_coupled(
            control, days=spatial_days, initial=spatial_initial,
            parameters=spatial_parameters,
        )
        spatial_treated = simulate_spatial_coupled(
            treated, days=spatial_days, initial=spatial_initial,
            parameters=spatial_parameters,
        )
        c, t = control.iloc[-1], treated.iloc[-1]
        def endpoint(frame: pd.DataFrame, species: str, zone: str) -> float:
            selected = frame[(frame.species == species) & (frame.zone == zone)]
            return float(selected.amount.iloc[-1])

        control_exit = endpoint(spatial_control, "effector", "efferent_lymph")
        treated_exit = endpoint(spatial_treated, "effector", "efferent_lymph")
        rows.append({
            "patient_id": row["patient_id"],
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
            "qc_spatial_initialization_measured": int(measured_states is not None),
            "qc_spatial_zone_coverage": spatial_zone_coverage,
            "qc_systemic_parameters_supplied": int(systemic_parameters is not None),
            "qc_spatial_parameters_supplied": int(spatial_parameters is not None),
        })
    return features.merge(pd.DataFrame(rows), on="patient_id", validate="one_to_one")
