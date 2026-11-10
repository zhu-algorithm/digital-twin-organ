"""Input validation and optional AnnData conversion."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

from .features import SAMPLE_COLUMNS, TEMPORAL_COLUMNS


def _finite_nonnegative_column(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        raise ValueError(f"cell-quality column not found: {column}")
    values = pd.to_numeric(frame[column], errors="coerce")
    if (values.isna().any() or not np.isfinite(values).all() or (values < 0).any()):
        raise ValueError(f"cell-quality column must be complete, finite and non-negative: {column}")
    return values


def filter_cells_by_quality(cells: pd.DataFrame, policy: dict | None = None) -> pd.DataFrame:
    """Apply an explicit, assay-unit-aware cell QC policy without guessing columns."""
    if policy is None:
        return cells
    if not isinstance(policy, dict):
        raise ValueError("cell-quality policy must be a mapping")
    allowed = {
        "n_genes_column", "min_genes", "max_genes",
        "total_counts_column", "min_total_counts", "max_total_counts",
        "mitochondrial_fraction_column", "max_mitochondrial_fraction",
        "doublet_column",
    }
    unknown = set(policy) - allowed
    if unknown:
        raise ValueError(f"unknown cell-quality policy fields: {sorted(unknown)}")
    normalized = {key: value for key, value in policy.items() if value is not None}
    if not normalized:
        return cells
    mask = pd.Series(True, index=cells.index)
    for prefix, lower_name, upper_name in (
        ("n_genes", "min_genes", "max_genes"),
        ("total_counts", "min_total_counts", "max_total_counts"),
    ):
        column = normalized.get(f"{prefix}_column")
        lower, upper = normalized.get(lower_name), normalized.get(upper_name)
        if (lower is not None or upper is not None) and column is None:
            raise ValueError(f"{prefix}_column is required when its bounds are supplied")
        if column is not None and lower is None and upper is None:
            raise ValueError(f"{prefix} column requires at least one bound")
        if column is not None:
            values = _finite_nonnegative_column(cells, str(column))
            for name, bound in ((lower_name, lower), (upper_name, upper)):
                if bound is not None and (not np.isfinite(float(bound)) or float(bound) < 0):
                    raise ValueError(f"{name} must be finite and non-negative")
            if lower is not None and upper is not None and float(lower) > float(upper):
                raise ValueError(f"{lower_name} cannot exceed {upper_name}")
            if lower is not None:
                mask &= values >= float(lower)
            if upper is not None:
                mask &= values <= float(upper)
    mito_column = normalized.get("mitochondrial_fraction_column")
    mito_max = normalized.get("max_mitochondrial_fraction")
    if (mito_column is None) != (mito_max is None):
        raise ValueError("mitochondrial fraction column and maximum must be supplied together")
    if mito_column is not None:
        values = _finite_nonnegative_column(cells, str(mito_column))
        maximum = float(mito_max)
        if not np.isfinite(maximum) or not 0 <= maximum <= 1:
            raise ValueError("max_mitochondrial_fraction must lie in [0, 1]")
        if (values > 1).any():
            raise ValueError("mitochondrial fraction values must use the [0, 1] scale")
        mask &= values <= maximum
    doublet_column = normalized.get("doublet_column")
    if doublet_column is not None:
        if doublet_column not in cells:
            raise ValueError(f"cell-quality column not found: {doublet_column}")
        raw = cells[doublet_column]
        if raw.isna().any():
            raise ValueError("doublet indicator must be complete")
        text = raw.astype(str).str.strip().str.lower()
        if not text.isin({"true", "false", "1", "0"}).all():
            raise ValueError("doublet indicator must contain only true/false or 1/0")
        mask &= text.isin({"false", "0"})
    retained = cells.loc[mask].copy()
    original_counts = cells.groupby(cells.patient_id.astype(str)).size()
    retained_counts = retained.groupby(retained.patient_id.astype(str)).size()
    removed_patients = sorted(set(original_counts.index) - set(retained_counts.index))
    if removed_patients:
        raise ValueError(f"cell-quality filtering removed all cells for patients: {removed_patients[:10]}")
    retained.attrs["cell_quality_summary"] = {
        patient: {"input": int(count), "retained": int(retained_counts.get(patient, 0))}
        for patient, count in original_counts.items()
    }
    retained.attrs["cell_quality_policy"] = normalized
    return retained


def _validate_cell_identifiers(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"patient_id", "cell_type"}
    if not required.issubset(frame.columns):
        raise ValueError("cell data must contain patient_id and cell_type")
    if frame[["patient_id", "cell_type"]].isna().any().any():
        raise ValueError("patient_id and cell_type must be non-missing")
    out = frame.copy()
    for column in ("patient_id", "cell_type"):
        out[column] = out[column].astype(str).str.strip()
        if out[column].eq("").any():
            raise ValueError(f"{column} must be non-empty")
    return out


def load_cells(
    path: str | Path,
    metadata_columns: tuple[str, ...] = (),
    gene_symbol_column: str | None = None,
    expression_layer: str | None = None,
) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        if gene_symbol_column or expression_layer:
            raise ValueError(
                "gene_symbol_column and expression_layer apply only to .h5ad input"
            )
        return _validate_cell_identifiers(pd.read_csv(path))
    if path.suffix.lower() == ".h5ad":
        try:
            import anndata as ad
        except ImportError as exc:
            raise ImportError("install tdln-pd1[single-cell] for .h5ad input") from exc
        adata = ad.read_h5ad(path)
        required = {"patient_id", "cell_type"}
        if not required.issubset(adata.obs.columns):
            raise ValueError("AnnData obs must contain patient_id and cell_type")
        if adata.obs[["patient_id", "cell_type"]].isna().any().any():
            raise ValueError("AnnData patient_id and cell_type must be non-missing")
        if gene_symbol_column:
            if gene_symbol_column not in adata.var.columns:
                raise ValueError(
                    f"AnnData var gene-symbol column not found: {gene_symbol_column}"
                )
            raw_symbols = adata.var[gene_symbol_column]
            if raw_symbols.isna().any():
                raise ValueError("AnnData gene symbols must be complete")
            symbols = raw_symbols.astype(str).str.strip().tolist()
            if any(not symbol for symbol in symbols):
                raise ValueError("AnnData gene symbols must be non-empty")
        else:
            symbols = [str(value).strip() for value in adata.var_names]
        selected = [index for index, symbol in enumerate(symbols)
                    if symbol in _signature_genes()]
        genes = [symbols[index] for index in selected]
        duplicates = sorted({gene for gene in genes if genes.count(gene) > 1})
        if duplicates:
            raise ValueError(
                f"AnnData signature gene symbols must be unique: {duplicates}"
            )
        if expression_layer:
            if expression_layer not in adata.layers:
                raise ValueError(f"AnnData expression layer not found: {expression_layer}")
            values = adata.layers[expression_layer][:, selected]
            if sparse.issparse(values):
                values = values.toarray()
            matrix = pd.DataFrame(np.asarray(values), columns=genes)
        elif gene_symbol_column:
            values = adata.X[:, selected]
            if sparse.issparse(values):
                values = values.toarray()
            matrix = pd.DataFrame(np.asarray(values), columns=genes)
        else:
            # Preserve AnnData's own backed/sparse slicing behavior while only
            # materializing the small prespecified signature-gene panel.
            matrix = adata[:, genes].to_df()
        matrix.insert(0, "cell_type", adata.obs["cell_type"].values)
        matrix.insert(0, "patient_id", adata.obs["patient_id"].values)
        metadata = dict.fromkeys(("clonotype_id", "tissue", *TEMPORAL_COLUMNS,
                                  *SAMPLE_COLUMNS,
                                  *metadata_columns))
        for column in metadata:
            if column in adata.obs:
                if column in matrix:
                    raise ValueError(f"AnnData metadata collides with gene column: {column}")
                matrix.insert(len(matrix.columns), column, adata.obs[column].astype(str).values)
        return _validate_cell_identifiers(matrix)
    raise ValueError("input must be .csv or .h5ad")


def _signature_genes() -> set[str]:
    from .features import SIGNATURES
    return {gene for genes in SIGNATURES.values() for gene in genes}


def load_labels(path: str | Path) -> pd.DataFrame:
    labels = pd.read_csv(path)
    required = {"patient_id", "response"}
    if not required.issubset(labels.columns):
        raise ValueError("labels need patient_id,response; response must be 0/1")
    if labels.patient_id.isna().any():
        raise ValueError("label patient_id must be non-missing")
    labels = labels.copy()
    labels["patient_id"] = labels.patient_id.astype(str).str.strip()
    if labels.patient_id.eq("").any() or labels.patient_id.duplicated().any():
        raise ValueError("label patient_id must be unique and non-empty")
    response = pd.to_numeric(labels.response, errors="coerce")
    if response.isna().any() or not response.isin([0, 1]).all():
        raise ValueError("response must be complete and contain only 0 and 1")
    labels["response"] = response.astype(int)
    return labels
