"""Patient-level immune features from cell-level expression tables."""

from __future__ import annotations

import hashlib
import json
import re
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
# Cell-state summaries are most interpretable within the lineage in which the
# underlying programme is expected.  Global summaries remain available for
# backwards compatibility and exploratory use.
LINEAGE_SIGNATURES = {
    "sig_cytotoxic_cd8_t": ("cytotoxic", r"\bcd8\b|cytotoxic|tpex|exhausted[ _-]*t"),
    "sig_exhaustion_cd8_t": ("exhaustion", r"\bcd8\b|cytotoxic|tpex|exhausted[ _-]*t"),
    "sig_tcf7_memory_cd8_t": ("tcf7_memory", r"\bcd8\b|cytotoxic|tpex|exhausted[ _-]*t"),
    "sig_treg_treg": ("treg", r"\btreg\b|regulatory[ _-]*t"),
    "sig_dendritic_activation_dc": ("dendritic_activation", r"\bdc\b|dendritic"),
}
TEMPORAL_COLUMNS = ("sampling_day", "timepoint", "visit", "collection_timepoint")
UNASSIGNED_CLONOTYPES = {"", "nan", "none", "unassigned", "na", "n/a"}
MIN_SIGNATURE_GENE_COVERAGE = 0.5


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


def _validate_signature_coverage_threshold(value: float) -> float:
    try:
        threshold = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("minimum signature gene coverage must be numeric") from error
    if not np.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError("minimum signature gene coverage must fall within (0, 1]")
    return threshold


def _mean_present(
    frame: pd.DataFrame, genes: tuple[str, ...], minimum_coverage: float,
) -> pd.Series:
    present = [g for g in genes if g in frame.columns]
    if len(present) / len(genes) < minimum_coverage:
        return pd.Series(np.nan, index=frame.index)
    return frame[present].mean(axis=1)


def signature_gene_coverage(
    columns, minimum_coverage: float = MIN_SIGNATURE_GENE_COVERAGE,
) -> dict[str, dict]:
    minimum_coverage = _validate_signature_coverage_threshold(minimum_coverage)
    available = set(columns)
    return {
        name: {
            "present": [gene for gene in genes if gene in available],
            "missing": [gene for gene in genes if gene not in available],
            "fraction": float(sum(gene in available for gene in genes) / len(genes)),
            "usable": bool(
                sum(gene in available for gene in genes) / len(genes)
                >= minimum_coverage
            ),
        }
        for name, genes in SIGNATURES.items()
    }


def _validate_expression_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Preserve true missingness but reject corrupt supplied expression values."""
    out = frame.copy()
    signature_genes = {gene for genes in SIGNATURES.values() for gene in genes}
    for gene in sorted(signature_genes.intersection(out.columns)):
        raw = out[gene]
        numeric = pd.to_numeric(raw, errors="coerce")
        invalid = raw.notna() & numeric.isna()
        if invalid.any():
            raise ValueError(f"expression values must be numeric or missing: {gene}")
        if np.isinf(numeric.to_numpy(dtype=float)).any():
            raise ValueError(f"expression values must be finite or missing: {gene}")
        out[gene] = numeric
    return out


def _cell_type_feature_names(cell_types: pd.Series) -> dict[str, str]:
    """Build stable fraction names and reject lossy label canonicalization."""
    labels = sorted(set(cell_types.astype(str)))
    mapping = {
        label: re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
        for label in labels
    }
    empty = [label for label, canonical in mapping.items() if not canonical]
    if empty:
        raise ValueError(f"cell_type cannot form a model feature name: {empty[:10]}")
    inverse: dict[str, list[str]] = {}
    for label, canonical in mapping.items():
        inverse.setdefault(canonical, []).append(label)
    collisions = {name: values for name, values in inverse.items() if len(values) > 1}
    if collisions:
        raise ValueError(f"cell_type labels collide after feature normalization: {collisions}")
    return mapping


def _assigned_clonotypes(series: pd.Series) -> pd.Series:
    clones = series.dropna().astype(str).str.strip()
    return clones[~clones.str.lower().isin(UNASSIGNED_CLONOTYPES)]


def aggregate_cells(
    cells: pd.DataFrame,
    time_column: str | None = None,
    baseline_values: tuple[str, ...] | None = None,
    feature_tissue_pattern: str | None = r"tdln|lymph",
    min_signature_gene_coverage: float = MIN_SIGNATURE_GENE_COVERAGE,
    cell_type_vocabulary: tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Aggregate long cell rows into one leakage-safe row per patient.

    Required columns are patient_id and cell_type. Gene columns are optional;
    absent signatures are retained as NaN for downstream imputation.
    """
    min_signature_gene_coverage = _validate_signature_coverage_threshold(
        min_signature_gene_coverage
    )
    required = {"patient_id", "cell_type"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    if cells[["patient_id", "cell_type"]].isna().any().any():
        raise ValueError("patient_id and cell_type must be non-missing")
    x = cells.copy()
    for column in ("patient_id", "cell_type"):
        x[column] = x[column].astype(str).str.strip()
        if x[column].eq("").any():
            raise ValueError(f"{column} must be non-empty")
    x = select_baseline_cells(x, time_column, baseline_values)
    x = _validate_expression_columns(x)
    vocabulary_cells = x
    if ("tissue" in x and feature_tissue_pattern
            and feature_tissue_pattern.lower() != "all"):
        vocabulary_cells = x.loc[
            x.tissue.astype(str).str.contains(
                feature_tissue_pattern, case=False, regex=True, na=False,
            )
        ]
    if cell_type_vocabulary is not None:
        vocabulary = tuple(str(value).strip() for value in cell_type_vocabulary)
        if (not vocabulary or any(not value for value in vocabulary)
                or len(set(vocabulary)) != len(vocabulary)):
            raise ValueError("cell_type_vocabulary must contain unique non-empty labels")
        unknown = sorted(set(vocabulary_cells["cell_type"]) - set(vocabulary))
        if unknown:
            raise ValueError(f"cell types are outside the locked vocabulary: {unknown[:10]}")
        cell_type_names = _cell_type_feature_names(pd.Series(vocabulary, dtype=str))
    else:
        vocabulary = tuple(sorted(set(vocabulary_cells["cell_type"])))
        cell_type_names = _cell_type_feature_names(vocabulary_cells["cell_type"])
    vocabulary_sha256 = hashlib.sha256(json.dumps(
        vocabulary, ensure_ascii=False, separators=(",", ":"),
    ).encode("utf-8")).hexdigest()
    coverage = signature_gene_coverage(x.columns, min_signature_gene_coverage)
    for name, genes in SIGNATURES.items():
        x[f"sig_{name}"] = _mean_present(x, genes, min_signature_gene_coverage)
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
            "qc_min_signature_gene_coverage": min_signature_gene_coverage,
            "qc_cell_type_vocabulary_sha256": vocabulary_sha256,
            "qc_cell_type_vocabulary_locked": int(cell_type_vocabulary is not None),
        }
        for cell_type in vocabulary:
            row[f"frac_{cell_type_names[cell_type]}"] = 0.0
        fractions = feature_group["cell_type"].astype(str).value_counts(normalize=True)
        for cell_type, value in fractions.items():
            safe = cell_type_names[cell_type]
            row[f"frac_{safe}"] = float(value)
        for name in SIGNATURES:
            row[f"sig_{name}"] = float(feature_group[f"sig_{name}"].mean())
            row[f"qc_sig_{name}_gene_coverage"] = coverage[name]["fraction"]
        normalized_cell_type = feature_group["cell_type"].astype(str).str.lower()
        for feature_name, (signature_name, lineage_pattern) in LINEAGE_SIGNATURES.items():
            lineage = feature_group.loc[
                normalized_cell_type.str.contains(lineage_pattern, regex=True, na=False)
            ]
            row[feature_name] = (
                float(lineage[f"sig_{signature_name}"].mean())
                if not lineage.empty else np.nan
            )
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
    result = pd.DataFrame(rows)
    # Cell-type fractions form a compositional vector over the selected feature
    # compartment.  A type absent from one patient's selected cells is an
    # observed zero, not an unmeasured value for median imputation.  Signature
    # and clonotype NaNs remain untouched because they can represent genuinely
    # unavailable assays or insufficient gene coverage.
    fraction_columns = [column for column in result if column.startswith("frac_")]
    if fraction_columns:
        result[fraction_columns] = result[fraction_columns].fillna(0.0)
    return result.fillna({"n_cells": 0})
