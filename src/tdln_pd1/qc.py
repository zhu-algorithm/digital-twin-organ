"""Pre-training quality gates for TDLN single-cell cohorts."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .features import (MIN_SIGNATURE_GENE_COVERAGE, SIGNATURES, _assigned_clonotypes,
                       _validate_signature_coverage_threshold, aggregate_cells,
                       select_baseline_cells, signature_gene_coverage)


def cohort_qc(
    cells: pd.DataFrame,
    labels: pd.DataFrame,
    time_column: str | None = None,
    baseline_values: tuple[str, ...] | None = None,
    feature_tissue_pattern: str | None = r"tdln|lymph",
    min_signature_gene_coverage: float = MIN_SIGNATURE_GENE_COVERAGE,
    sample_column: str | None = None,
) -> dict:
    try:
        min_signature_gene_coverage = _validate_signature_coverage_threshold(
            min_signature_gene_coverage
        )
    except ValueError as error:
        return {"passed": False, "errors": [str(error)], "warnings": []}
    required_cells = {"patient_id", "cell_type"}
    required_labels = {"patient_id", "response"}
    errors, warnings = [], []
    if missing := required_cells - set(cells):
        errors.append(f"cell table missing: {sorted(missing)}")
    if missing := required_labels - set(labels):
        errors.append(f"label table missing: {sorted(missing)}")
    if errors:
        return {"passed": False, "errors": errors, "warnings": warnings}
    if cells[["patient_id", "cell_type"]].isna().any().any():
        errors.append("cell patient_id and cell_type must be non-missing")
    if labels.patient_id.isna().any():
        errors.append("label patient_id must be non-missing")
    raw_cell_count = len(cells)
    derived_features = None
    derivation_error = None
    try:
        # Execute the exact downstream transformation as part of the gate. This
        # prevents QC from passing data that later fail expression, identifier,
        # cell-type naming, temporal or tissue-compartment contracts.
        derived_features = aggregate_cells(
            cells, time_column, baseline_values, feature_tissue_pattern,
            min_signature_gene_coverage, sample_column=sample_column,
        )
    except ValueError as error:
        derivation_error = f"feature derivation failed: {error}"
    try:
        cells = select_baseline_cells(cells, time_column, baseline_values)
    except ValueError as error:
        errors.append(str(error))
        return {"passed": False, "errors": errors, "warnings": warnings,
                "n_cells_raw": raw_cell_count, "n_cells_selected": 0}
    normalized_label_ids = labels.patient_id.astype(str).str.strip()
    if normalized_label_ids.eq("").any() or normalized_label_ids.duplicated().any():
        errors.append("label table needs unique, non-empty patient_id")
    response = pd.to_numeric(labels.response, errors="coerce")
    if response.isna().any() or not response.isin([0, 1]).all():
        errors.append("response must be complete binary 0/1")
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
    feature_cell_count = len(cells)
    if "tissue" in cells and feature_tissue_pattern and feature_tissue_pattern.lower() != "all":
        feature_mask = cells.tissue.astype(str).str.contains(
            feature_tissue_pattern, case=False, regex=True, na=False,
        )
        feature_cells = cells.loc[feature_mask]
        feature_cell_count = len(feature_cells)
        feature_patients = set(feature_cells.patient_id.astype(str))
        missing_feature_tissue = sorted(cell_patients - feature_patients)
        if missing_feature_tissue:
            errors.append(
                "patients without feature-compartment cells "
                f"for pattern {feature_tissue_pattern!r}: {missing_feature_tissue[:10]}"
            )
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
    signature_coverage = signature_gene_coverage(
        cells.columns, min_signature_gene_coverage,
    )
    unusable_signatures = sorted(
        name for name, record in signature_coverage.items() if not record["usable"]
    )
    if unusable_signatures:
        warnings.append(
            "signatures below minimum gene coverage are emitted as missing: "
            f"{unusable_signatures}"
        )
    tissue_report = {"available": "tissue" in cells, "paired_tdln_tumour_patients": 0}
    if "tissue" in cells:
        tissue = cells.assign(_tissue=cells.tissue.astype(str).str.lower())
        flags = tissue.groupby(tissue.patient_id.astype(str))._tissue.agg(
            lambda x: {"tdln": x.str.contains("tdln|lymph").any(),
                       "tumour": x.str.contains("tumou|tumor").any()})
        tissue_report["paired_tdln_tumour_patients"] = int(sum(v["tdln"] and v["tumour"] for v in flags))
    clonotype_report = {"available": "clonotype_id" in cells, "assigned_fraction": None}
    if "clonotype_id" in cells:
        assigned = _assigned_clonotypes(cells.clonotype_id)
        clonotype_report["assigned_fraction"] = float(len(assigned) / len(cells)) if len(cells) else 0.0
    if derivation_error is not None:
        errors.append(derivation_error)
    return {
        "passed": not errors, "errors": errors, "warnings": warnings,
        "n_cells": int(len(cells)), "n_cells_raw": int(raw_cell_count),
        "n_cells_selected": int(len(cells)),
        "temporal_selection": {"time_column": time_column,
                               "baseline_values": list(baseline_values or ())},
        "sample_aggregation": {"sample_column": sample_column,
                               "weighting": ("equal_samples" if sample_column else "cells")},
        "feature_compartment": {"pattern": feature_tissue_pattern,
                                "n_cells": int(feature_cell_count)},
        "feature_derivation": {
            "validated": derived_features is not None,
            "n_patient_rows": int(len(derived_features)) if derived_features is not None else 0,
            "shared_with_featurize": True,
        },
        "n_cell_patients": len(cell_patients),
        "n_label_patients": len(label_patients), "n_matched_patients": len(matched),
        "cells_per_patient": {"min": int(counts.min()) if len(counts) else 0,
                              "median": float(np.median(counts)) if len(counts) else 0.0,
                              "max": int(counts.max()) if len(counts) else 0},
        "response_counts": {str(key): int(value) for key, value in response_counts.items()},
        "signature_genes_present": present_genes, "signature_genes_missing": missing_genes,
        "signature_gene_coverage": signature_coverage,
        "minimum_signature_gene_coverage": float(min_signature_gene_coverage),
        "tissue_pairing": tissue_report, "clonotypes": clonotype_report,
    }
