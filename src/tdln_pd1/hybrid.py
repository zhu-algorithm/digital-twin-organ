"""Fuse patient-level single-cell features with mechanistic counterfactuals."""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pandas as pd

from .simulator import Dose, InitialState, simulate
from .spatial_dynamics import SpatialInitialState, simulate_spatial_coupled


def _finite(row: pd.Series, name: str, default: float) -> float:
    value = row.get(name, default)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if np.isfinite(value) else default


def _fraction(row: pd.Series, names: tuple[str, ...], default: float) -> float:
    """Read the first available canonicalized cell-fraction alias."""
    for name in names:
        if name in row and pd.notna(row[name]):
            return float(np.clip(_finite(row, name, default), 0.0, 1.0))
    return default


def initial_state_from_features(row: pd.Series, base: InitialState | None = None) -> InitialState:
    """Map measured immune features to initial states with explicit bounded rules.

    Signature inputs are assumed non-negative expression summaries; fractions are
    in [0,1].  The mapping is intentionally simple and inspectable so a calibrated
    mapping can replace it without changing the simulator.
    """
    base = base or InitialState()
    cytotoxic = max(_finite(row, "sig_cytotoxic", 0.0), 0.0)
    exhaustion = max(_finite(row, "sig_exhaustion", 0.0), 0.0)
    tcf7 = max(_finite(row, "sig_tcf7_memory", 0.0), 0.0)
    dendritic = max(_finite(row, "sig_dendritic_activation", 0.0), 0.0)
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
    dendritic = max(_finite(row, "sig_dendritic_activation", 0.0), 0.0)
    tcf7 = max(_finite(row, "sig_tcf7_memory", 0.0), 0.0)
    cytotoxic = max(_finite(row, "sig_cytotoxic", 0.0), 0.0)
    cd8_fraction = _fraction(row, ("frac_cd8_t", "frac_cd8_t_cell"), 0.15)
    antigen_proxy = max(_finite(row, "tumour_antigen_proxy", 0.5), 0.0)
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
) -> pd.DataFrame:
    """Add per-patient treated/control mechanistic endpoints without using labels."""
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    if features.patient_id.duplicated().any():
        raise ValueError("features must have one row per patient")
    doses = tuple(Dose(day, dose) for day in dose_days)
    rows = []
    for row in features.to_dict(orient="records"):
        initial = initial_state_from_features(pd.Series(row))
        control = simulate(days=days, initial=initial)
        treated = simulate(days=days, doses=doses, initial=initial)
        spatial_initial = spatial_initial_state_from_features(pd.Series(row))
        spatial_days = min(days, 14.0)
        spatial_control = simulate_spatial_coupled(
            control, days=spatial_days, initial=spatial_initial,
        )
        spatial_treated = simulate_spatial_coupled(
            treated, days=spatial_days, initial=spatial_initial,
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
        })
    return features.merge(pd.DataFrame(rows), on="patient_id", validate="one_to_one")
