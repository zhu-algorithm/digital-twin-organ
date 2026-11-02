"""Input validation and optional AnnData conversion."""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from .features import TEMPORAL_COLUMNS


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
) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
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
        genes = [g for g in adata.var_names if g in _signature_genes()]
        matrix = adata[:, genes].to_df()
        matrix.insert(0, "cell_type", adata.obs["cell_type"].values)
        matrix.insert(0, "patient_id", adata.obs["patient_id"].values)
        metadata = dict.fromkeys(("clonotype_id", "tissue", *TEMPORAL_COLUMNS,
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
