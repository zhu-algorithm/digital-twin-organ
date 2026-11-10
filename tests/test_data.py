import sys
from types import SimpleNamespace

import pandas as pd
import pytest
import numpy as np
from scipy import sparse

from tdln_pd1.data import filter_cells_by_quality, load_cells, load_labels
from tdln_pd1.features import aggregate_cells


class _FakeAnnData:
    def __init__(self):
        self.obs = pd.DataFrame({
            "patient_id": ["p1", "p1", "p2", "p2"],
            "cell_type": ["CD8 T"] * 4,
            "tissue": ["TDLN"] * 4,
            "clonotype_id": ["a", "b", "c", "d"],
            "timepoint": ["baseline", "post", "baseline", "post"],
            "custom_visit": ["screen", "cycle2", "screen", "cycle2"],
        })
        self.var_names = ["NKG7", "GZMB", "irrelevant"]
        self._expression = pd.DataFrame({
            "NKG7": [1.0, 100.0, 2.0, 200.0],
            "GZMB": [1.0, 100.0, 2.0, 200.0],
            "irrelevant": [0.0] * 4,
        })

    def __getitem__(self, key):
        _, genes = key
        parent = self

        class _View:
            def to_df(self):
                return parent._expression.loc[:, genes].copy()

        return _View()


def test_h5ad_preserves_temporal_tissue_and_clonotype_metadata(monkeypatch, tmp_path):
    fake = _FakeAnnData()
    monkeypatch.setitem(sys.modules, "anndata", SimpleNamespace(read_h5ad=lambda _: fake))
    path = tmp_path / "cohort.h5ad"
    path.touch()
    cells = load_cells(path, metadata_columns=("custom_visit",))
    assert {"timepoint", "custom_visit", "tissue", "clonotype_id"}.issubset(cells)
    with pytest.raises(ValueError, match="temporal columns detected"):
        aggregate_cells(cells)
    features = aggregate_cells(cells, "timepoint", ("baseline",))
    assert features.set_index("patient_id").loc["p1", "sig_cytotoxic"] == 1.0
    assert features.set_index("patient_id").loc["p2", "sig_cytotoxic"] == 2.0


def test_h5ad_explicit_gene_symbols_and_sparse_expression_layer(monkeypatch, tmp_path):
    fake = _FakeAnnData()
    fake.var_names = ["ENSG1", "ENSG2", "ENSG3"]
    fake.var = pd.DataFrame({
        "gene_symbol": ["NKG7", "GZMB", "irrelevant"],
    }, index=fake.var_names)
    fake.X = np.zeros((4, 3))
    fake.layers = {
        "lognorm": sparse.csr_matrix([
            [1.0, 1.0, 0.0], [100.0, 100.0, 0.0],
            [2.0, 2.0, 0.0], [200.0, 200.0, 0.0],
        ])
    }
    monkeypatch.setitem(sys.modules, "anndata", SimpleNamespace(read_h5ad=lambda _: fake))
    path = tmp_path / "cohort.h5ad"
    path.touch()
    cells = load_cells(
        path, gene_symbol_column="gene_symbol", expression_layer="lognorm")
    assert list(cells[["NKG7", "GZMB"]].iloc[0]) == [1.0, 1.0]
    features = aggregate_cells(cells, "timepoint", ("baseline",))
    assert features.set_index("patient_id").loc["p2", "sig_cytotoxic"] == 2.0


def test_h5ad_gene_mapping_fails_closed_on_missing_or_duplicate_symbols(
    monkeypatch, tmp_path,
):
    fake = _FakeAnnData()
    fake.var = pd.DataFrame({"symbol": ["NKG7", "NKG7", "other"]})
    fake.X = np.zeros((4, 3))
    fake.layers = {}
    monkeypatch.setitem(sys.modules, "anndata", SimpleNamespace(read_h5ad=lambda _: fake))
    path = tmp_path / "cohort.h5ad"
    path.touch()
    with pytest.raises(ValueError, match="must be unique"):
        load_cells(path, gene_symbol_column="symbol")
    with pytest.raises(ValueError, match="not found"):
        load_cells(path, gene_symbol_column="missing")
    fake.var["symbol"] = ["NKG7", "GZMB", "other"]
    with pytest.raises(ValueError, match="layer not found"):
        load_cells(path, gene_symbol_column="symbol", expression_layer="missing")


def test_csv_rejects_anndata_only_expression_options(tmp_path):
    path = tmp_path / "cells.csv"
    pd.DataFrame({"patient_id": ["p1"], "cell_type": ["T"]}).to_csv(
        path, index=False)
    with pytest.raises(ValueError, match="only to .h5ad"):
        load_cells(path, gene_symbol_column="symbol")


def test_cell_loader_rejects_missing_identifiers_for_csv_and_h5ad(monkeypatch, tmp_path):
    csv_path = tmp_path / "bad.csv"
    pd.DataFrame({"patient_id": ["p1", None], "cell_type": ["T", "T"]}).to_csv(
        csv_path, index=False,
    )
    with pytest.raises(ValueError, match="non-missing"):
        load_cells(csv_path)

    fake = _FakeAnnData()
    fake.obs.loc[0, "patient_id"] = None
    monkeypatch.setitem(sys.modules, "anndata", SimpleNamespace(read_h5ad=lambda _: fake))
    h5ad_path = tmp_path / "bad.h5ad"
    h5ad_path.touch()
    with pytest.raises(ValueError, match="non-missing"):
        load_cells(h5ad_path)


def test_label_loader_requires_complete_binary_outcomes_and_unique_patients(tmp_path):
    missing = tmp_path / "missing.csv"
    pd.DataFrame({"patient_id": ["p1", "p2"], "response": [0, None]}).to_csv(
        missing, index=False,
    )
    with pytest.raises(ValueError, match="complete"):
        load_labels(missing)

    duplicate = tmp_path / "duplicate.csv"
    pd.DataFrame({"patient_id": ["p1", "p1"], "response": [0, 1]}).to_csv(
        duplicate, index=False,
    )
    with pytest.raises(ValueError, match="unique"):
        load_labels(duplicate)


def test_explicit_cell_quality_policy_filters_and_audits_retention():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 4 + ["p2"] * 2,
        "cell_type": ["CD8 T"] * 6,
        "n_genes": [100, 500, 500, 500, 500, 500],
        "counts": [1000, 1000, 1000, 1000, 1000, 1000],
        "mito_fraction": [0.05, 0.05, 0.4, 0.05, 0.05, 0.05],
        "is_doublet": [False, False, False, True, False, False],
        "NKG7": [1, 2, 3, 4, 5, 7], "GZMB": [1, 2, 3, 4, 5, 7],
    })
    policy = {"n_genes_column": "n_genes", "min_genes": 200,
              "total_counts_column": "counts", "min_total_counts": 500,
              "mitochondrial_fraction_column": "mito_fraction",
              "max_mitochondrial_fraction": 0.2,
              "doublet_column": "is_doublet"}
    filtered = filter_cells_by_quality(cells, policy)
    assert filtered.groupby("patient_id").size().to_dict() == {"p1": 1, "p2": 2}
    features = aggregate_cells(filtered, feature_tissue_pattern="all").set_index("patient_id")
    assert features.loc["p1", "qc_n_cells_input_quality"] == 4
    assert features.loc["p1", "qc_n_cells_retained_quality"] == 1
    assert features.loc["p1", "qc_cell_quality_retained_fraction"] == pytest.approx(0.25)
    assert features.qc_cell_quality_policy_applied.eq(1).all()
    assert features.qc_cell_quality_policy_sha256.str.len().eq(64).all()
    changed = aggregate_cells(filter_cells_by_quality(
        cells, {**policy, "min_genes": 50}), feature_tissue_pattern="all")
    assert changed.qc_cell_quality_policy_sha256.iloc[0] != (
        features.qc_cell_quality_policy_sha256.iloc[0])


def test_cell_quality_policy_fails_closed_on_units_and_patient_loss():
    cells = pd.DataFrame({"patient_id": ["p1", "p2"], "cell_type": ["T", "T"],
                          "n_genes": [100, 500], "mito": [5.0, 0.1]})
    with pytest.raises(ValueError, match=r"\[0, 1\] scale"):
        filter_cells_by_quality(cells, {
            "mitochondrial_fraction_column": "mito",
            "max_mitochondrial_fraction": 0.2})
    with pytest.raises(ValueError, match="removed all cells"):
        filter_cells_by_quality(cells, {
            "n_genes_column": "n_genes", "min_genes": 200})
    with pytest.raises(ValueError, match="requires at least one bound"):
        filter_cells_by_quality(cells, {"n_genes_column": "n_genes"})
