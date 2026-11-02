"""Patient-level immune features from cell-level expression tables."""

from __future__ import annotations

import numpy as np
import pandas as pd

SIGNATURES = {
    "cytotoxic": ("NKG7", "GNLY", "GZMB", "PRF1"),
    "exhaustion": ("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"),
    "interferon": ("IFNG", "STAT1", "CXCL9", "CXCL10", "ISG15"),
    "treg": ("FOXP3", "IL2RA", "CTLA4"),
    "tcf7_memory": ("TCF7", "CCR7", "IL7R", "LTB"),
    "dendritic_activation": ("CD80", "CD86", "CCR7", "HLA-DRA"),
}
TEMPORAL_COLUMNS = ("sampling_day", "timepoint", "visit", "collection_timepoint")
UNASSIGNED_CLONOTYPES = {"", "nan", "none", "unassigned", "na", "n/a"}


def select_baseline_cells(
    cells: pd.DataFrame,
    time_column: str | None = None,
    baseline_values: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Select an explicitly declared pre-treatment slice, failing on ambiguity."""
    detected = [column for column in TEMPORAL_COLUMNS if column in cells.columns]
    if time_column is None:
        if detected:
            raise ValueError(
                f"temporal columns detected {detected}; specify time_column and baseline_values"
            )
        return cells.copy()
    if time_column not in cells:
        raise ValueError(f"time column not found: {time_column}")
    if not baseline_values:
        raise ValueError("baseline_values are required when time_column is specified")
    values = cells[time_column]
    text = values.astype(str).str.strip().str.lower()
    mask = pd.Series(False, index=cells.index)
    numeric = pd.to_numeric(values, errors="coerce")
    for baseline in baseline_values:
        baseline_text = str(baseline).strip().lower()
        mask |= text.eq(baseline_text)
        try:
            numeric_value = float(baseline_text)
        except ValueError:
            continue
        mask |= np.isclose(numeric, numeric_value, rtol=0.0, atol=1e-12)
    selected = cells.loc[mask].copy()
    if selected.empty:
        raise ValueError("baseline selection retained zero cells")
    all_patients = set(cells.patient_id.astype(str)) if "patient_id" in cells else set()
    selected_patients = set(selected.patient_id.astype(str)) if "patient_id" in selected else set()
    missing = sorted(all_patients - selected_patients)
    if missing:
        raise ValueError(f"baseline selection removed all cells for patients: {missing[:10]}")
    return selected


def _mean_present(frame: pd.DataFrame, genes: tuple[str, ...]) -> pd.Series:
    present = [g for g in genes if g in frame.columns]
    if not present:
        return pd.Series(np.nan, index=frame.index)
    return frame[present].apply(pd.to_numeric, errors="coerce").mean(axis=1)


def _assigned_clonotypes(series: pd.Series) -> pd.Series:
    clones = series.dropna().astype(str).str.strip()
    return clones[~clones.str.lower().isin(UNASSIGNED_CLONOTYPES)]


def aggregate_cells(
    cells: pd.DataFrame,
    time_column: str | None = None,
    baseline_values: tuple[str, ...] | None = None,
    feature_tissue_pattern: str | None = r"tdln|lymph",
) -> pd.DataFrame:
    """Aggregate long cell rows into one leakage-safe row per patient.

    Required columns are patient_id and cell_type. Gene columns are optional;
    absent signatures are retained as NaN for downstream imputation.
    """
    required = {"patient_id", "cell_type"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    x = select_baseline_cells(cells, time_column, baseline_values)
    for name, genes in SIGNATURES.items():
        x[f"sig_{name}"] = _mean_present(x, genes)
    rows = []
    for patient, group in x.groupby("patient_id", sort=True):
        feature_group = group
        if "tissue" in group and feature_tissue_pattern and feature_tissue_pattern.lower() != "all":
            mask = group.tissue.astype(str).str.contains(feature_tissue_pattern, case=False, regex=True, na=False)
            feature_group = group.loc[mask]
            if feature_group.empty:
                raise ValueError(
                    f"patient {patient} has no cells matching feature_tissue_pattern={feature_tissue_pattern!r}"
                )
        row: dict[str, float | str] = {
            "patient_id": patient, "n_cells": len(feature_group), "n_cells_total_selected": len(group),
        }
        fractions = feature_group["cell_type"].astype(str).value_counts(normalize=True)
        for cell_type, value in fractions.items():
            safe = cell_type.lower().replace(" ", "_").replace("/", "_")
            row[f"frac_{safe}"] = float(value)
        for name in SIGNATURES:
            row[f"sig_{name}"] = float(feature_group[f"sig_{name}"].mean())
        if "clonotype_id" in group:
            clones = _assigned_clonotypes(group["clonotype_id"])
            counts = clones.value_counts()
            if len(counts):
                probabilities = counts.to_numpy(dtype=float) / counts.sum()
                entropy = float(-(probabilities * np.log(probabilities)).sum())
                row["tcr_clone_richness"] = float(len(counts))
                row["tcr_simpson_concentration"] = float(np.square(probabilities).sum())
                row["tcr_normalized_entropy"] = float(entropy / np.log(len(counts))) if len(counts) > 1 else 0.0
                row["tcr_top_clone_fraction"] = float(probabilities.max())
            else:
                row.update({"tcr_clone_richness": 0.0, "tcr_simpson_concentration": np.nan,
                            "tcr_normalized_entropy": np.nan, "tcr_top_clone_fraction": np.nan})
            if "tissue" in group:
                tissue = group["tissue"].astype(str).str.lower()
                tdln_clones = _assigned_clonotypes(
                    group.loc[tissue.str.contains("tdln|lymph"), "clonotype_id"]
                )
                tumour_clones = _assigned_clonotypes(
                    group.loc[tissue.str.contains("tumou|tumor"), "clonotype_id"]
                )
                tdln, tumour = set(tdln_clones), set(tumour_clones)
                shared = tdln & tumour
                union = tdln | tumour
                row["tcr_tdln_tumour_shared_fraction"] = (
                    float(len(shared) / len(tdln)) if tdln else np.nan
                )
                row["tcr_tdln_tumour_shared_cell_fraction"] = (
                    float(tdln_clones.isin(shared).mean()) if len(tdln_clones) else np.nan
                )
                row["tcr_tdln_tumour_jaccard"] = (
                    float(len(shared) / len(union)) if union else np.nan
                )
        rows.append(row)
    return pd.DataFrame(rows).fillna({"n_cells": 0})
