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


def _mean_present(frame: pd.DataFrame, genes: tuple[str, ...]) -> pd.Series:
    present = [g for g in genes if g in frame.columns]
    if not present:
        return pd.Series(np.nan, index=frame.index)
    return frame[present].apply(pd.to_numeric, errors="coerce").mean(axis=1)


def aggregate_cells(cells: pd.DataFrame) -> pd.DataFrame:
    """Aggregate long cell rows into one leakage-safe row per patient.

    Required columns are patient_id and cell_type. Gene columns are optional;
    absent signatures are retained as NaN for downstream imputation.
    """
    required = {"patient_id", "cell_type"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    x = cells.copy()
    for name, genes in SIGNATURES.items():
        x[f"sig_{name}"] = _mean_present(x, genes)
    rows = []
    for patient, group in x.groupby("patient_id", sort=True):
        row: dict[str, float | str] = {"patient_id": patient, "n_cells": len(group)}
        fractions = group["cell_type"].astype(str).value_counts(normalize=True)
        for cell_type, value in fractions.items():
            safe = cell_type.lower().replace(" ", "_").replace("/", "_")
            row[f"frac_{safe}"] = float(value)
        for name in SIGNATURES:
            row[f"sig_{name}"] = float(group[f"sig_{name}"].mean())
        if "clonotype_id" in group:
            clones = group["clonotype_id"].dropna().astype(str)
            clones = clones[~clones.str.lower().isin({"", "nan", "none", "unassigned"})]
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
                tdln = set(group.loc[tissue.str.contains("tdln|lymph"), "clonotype_id"].dropna().astype(str))
                tumour = set(group.loc[tissue.str.contains("tumou|tumor"), "clonotype_id"].dropna().astype(str))
                tdln.discard(""); tumour.discard("")
                row["tcr_tdln_tumour_shared_fraction"] = float(len(tdln & tumour) / len(tdln)) if tdln else np.nan
        rows.append(row)
    return pd.DataFrame(rows).fillna({"n_cells": 0})
