"""Input validation and optional AnnData conversion."""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def load_cells(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".h5ad":
        try:
            import anndata as ad
        except ImportError as exc:
            raise ImportError("install tdln-pd1[single-cell] for .h5ad input") from exc
        adata = ad.read_h5ad(path)
        required = {"patient_id", "cell_type"}
        if not required.issubset(adata.obs.columns):
            raise ValueError("AnnData obs must contain patient_id and cell_type")
        genes = [g for g in adata.var_names if g in _signature_genes()]
        matrix = adata[:, genes].to_df()
        matrix.insert(0, "cell_type", adata.obs["cell_type"].astype(str).values)
        matrix.insert(0, "patient_id", adata.obs["patient_id"].astype(str).values)
        for column in ("clonotype_id", "tissue"):
            if column in adata.obs:
                matrix.insert(len(matrix.columns), column, adata.obs[column].astype(str).values)
        return matrix
    raise ValueError("input must be .csv or .h5ad")


def _signature_genes() -> set[str]:
    from .features import SIGNATURES
    return {gene for genes in SIGNATURES.values() for gene in genes}


def load_labels(path: str | Path) -> pd.DataFrame:
    labels = pd.read_csv(path)
    required = {"patient_id", "response"}
    if not required.issubset(labels.columns):
        raise ValueError("labels need patient_id,response; response must be 0/1")
    if not set(labels["response"].dropna().unique()).issubset({0, 1}):
        raise ValueError("response must contain only 0 and 1")
    return labels
