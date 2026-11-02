"""Validated spatial-omics bridge into the six-zone TDLN state model."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .spatial import ZONES
from .spatial_dynamics import SPECIES


def aggregate_spatial_spots(spots: pd.DataFrame, weight_column: str | None = None) -> pd.DataFrame:
    """Aggregate segmented/deconvolved spots into patient-zone-species amounts.

    Inputs must already use the model's anatomical zone vocabulary and contain
    non-negative abundance estimates in the four mechanistic state columns.
    The function deliberately does not infer anatomy from coordinates alone.
    """
    required = {"patient_id", "zone", *SPECIES}
    missing = required.difference(spots.columns)
    if missing:
        raise ValueError(f"spatial spots missing columns: {sorted(missing)}")
    x = spots.copy()
    if x.empty:
        raise ValueError("spatial spot table must not be empty")
    if x.patient_id.isna().any() or x.zone.isna().any():
        raise ValueError("spatial patient_id and zone must be non-missing")
    x["patient_id"] = x.patient_id.astype(str).str.strip()
    x["zone"] = x.zone.astype(str).str.strip().str.lower()
    if x.patient_id.eq("").any():
        raise ValueError("spatial patient_id must be non-empty")
    valid_zones = {zone.name for zone in ZONES}
    invalid_zones = sorted(set(x.zone) - valid_zones)
    if invalid_zones:
        raise ValueError(f"unknown TDLN zones: {invalid_zones}")
    if "spot_id" in x and (
        x.spot_id.isna().any() or x.assign(spot_id=x.spot_id.astype(str)).duplicated(
            ["patient_id", "spot_id"]
        ).any()
    ):
        raise ValueError("spot_id must be unique within patient and non-missing")
    if weight_column in required:
        raise ValueError("weight column must be distinct from identifiers, zone and species")
    columns = list(SPECIES) + ([weight_column] if weight_column else [])
    if weight_column and weight_column not in x:
        raise ValueError(f"spatial spots missing weight column: {weight_column}")
    for column in columns:
        x[column] = pd.to_numeric(x[column], errors="coerce")
        if not np.isfinite(x[column]).all() or x[column].lt(0).any():
            raise ValueError(f"spatial column must be finite and non-negative: {column}")
    weight = x[weight_column] if weight_column else pd.Series(1.0, index=x.index)
    weighted = x.loc[:, list(SPECIES)].multiply(weight, axis=0)
    weighted["patient_id"] = x.patient_id
    weighted["zone"] = x.zone
    observed = weighted.groupby(["patient_id", "zone"], sort=True, as_index=False).sum()
    patients = sorted(observed.patient_id.unique())
    zones = [zone.name for zone in ZONES]
    complete_index = pd.MultiIndex.from_product(
        [patients, zones], names=["patient_id", "zone"]
    )
    observed_index = pd.MultiIndex.from_frame(observed[["patient_id", "zone"]])
    out = observed.set_index(["patient_id", "zone"]).reindex(
        complete_index, fill_value=0.0,
    ).reset_index()
    out["zone_observed"] = pd.MultiIndex.from_frame(
        out[["patient_id", "zone"]]
    ).isin(observed_index)
    return out.melt(id_vars=["patient_id", "zone", "zone_observed"],
                    value_vars=list(SPECIES), var_name="species", value_name="amount").sort_values(
        ["patient_id", "species", "zone"], ignore_index=True
    )
