"""Fuse patient-level single-cell features with mechanistic counterfactuals."""

from __future__ import annotations

from dataclasses import replace
import numpy as np
import pandas as pd

from .simulator import Dose, InitialState, simulate


def _finite(row: pd.Series, name: str, default: float) -> float:
    value = row.get(name, default)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return default
    return value if np.isfinite(value) else default


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
    treg_fraction = np.clip(_finite(row, "frac_treg", 0.05), 0.0, 1.0)
    cd8_fraction = np.clip(_finite(row, "frac_cd8_t", 0.15), 0.0, 1.0)
    return replace(
        base,
        mature_dc_tdln=base.mature_dc_tdln * np.clip(0.5 + dendritic, 0.25, 4.0),
        tpex_tdln=base.tpex_tdln * np.clip(0.5 + tcf7, 0.25, 4.0),
        effector_tumour=base.effector_tumour * np.clip(0.5 + cytotoxic + cd8_fraction, 0.25, 5.0),
        exhausted_tumour=base.exhausted_tumour * np.clip(0.5 + exhaustion, 0.25, 5.0),
        treg_tdln=base.treg_tdln * np.clip(0.5 + 4.0 * treg_fraction, 0.25, 4.0),
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
        c, t = control.iloc[-1], treated.iloc[-1]
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
        })
    return features.merge(pd.DataFrame(rows), on="patient_id", validate="one_to_one")
