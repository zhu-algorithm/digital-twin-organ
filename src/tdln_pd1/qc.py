"""Pre-training quality gates for TDLN single-cell cohorts."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import SIGNATURES


def cohort_qc(cells: pd.DataFrame, labels: pd.DataFrame) -> dict:
    required_cells = {"patient_id", "cell_type"}
    required_labels = {"patient_id", "response"}
    errors, warnings = [], []
    if missing := required_cells - set(cells):
        errors.append(f"cell table missing: {sorted(missing)}")
    if missing := required_labels - set(labels):
        errors.append(f"label table missing: {sorted(missing)}")
    if errors:
        return {"passed": False, "errors": errors, "warnings": warnings}
    if labels.patient_id.duplicated().any():
        errors.append("label table has duplicate patient_id")
    if not set(labels.response.dropna().unique()).issubset({0, 1}):
        errors.append("response must be binary 0/1")
    cell_patients = set(cells.patient_id.astype(str))
    label_patients = set(labels.patient_id.astype(str))
    matched = cell_patients & label_patients
    missing_labels = sorted(cell_patients - label_patients)
    missing_cells = sorted(label_patients - cell_patients)
    if missing_labels:
        errors.append(f"patients without labels: {missing_labels}")
    if missing_cells:
        errors.append(f"labels without cells: {missing_cells}")
    counts = cells.groupby(cells.patient_id.astype(str)).size()
    response_counts = labels[labels.patient_id.astype(str).isin(matched)].response.value_counts().to_dict()
    if len(response_counts) < 2:
        errors.append("matched cohort does not contain both response classes")
    elif min(response_counts.values()) < 2:
        warnings.append("fewer than two patients in one response class")
    signature_genes = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    present_genes = sorted(set(signature_genes) & set(cells.columns))
    missing_genes = sorted(set(signature_genes) - set(cells.columns))
    if len(present_genes) < 3:
        warnings.append("fewer than three configured signature genes are present")
    tissue_report = {"available": "tissue" in cells, "paired_tdln_tumour_patients": 0}
    if "tissue" in cells:
        tissue = cells.assign(_tissue=cells.tissue.astype(str).str.lower())
        flags = tissue.groupby(tissue.patient_id.astype(str))._tissue.agg(
            lambda x: {"tdln": x.str.contains("tdln|lymph").any(),
                       "tumour": x.str.contains("tumou|tumor").any()})
        tissue_report["paired_tdln_tumour_patients"] = int(sum(v["tdln"] and v["tumour"] for v in flags))
    clonotype_report = {"available": "clonotype_id" in cells, "assigned_fraction": None}
    if "clonotype_id" in cells:
        clone = cells.clonotype_id.dropna().astype(str)
        assigned = ~clone.str.lower().isin({"", "nan", "none", "unassigned"})
        clonotype_report["assigned_fraction"] = float(assigned.mean()) if len(clone) else 0.0
    return {
        "passed": not errors, "errors": errors, "warnings": warnings,
        "n_cells": int(len(cells)), "n_cell_patients": len(cell_patients),
        "n_label_patients": len(label_patients), "n_matched_patients": len(matched),
        "cells_per_patient": {"min": int(counts.min()) if len(counts) else 0,
                              "median": float(np.median(counts)) if len(counts) else 0.0,
                              "max": int(counts.max()) if len(counts) else 0},
        "response_counts": {str(key): int(value) for key, value in response_counts.items()},
        "signature_genes_present": present_genes, "signature_genes_missing": missing_genes,
        "tissue_pairing": tissue_report, "clonotypes": clonotype_report,
    }
