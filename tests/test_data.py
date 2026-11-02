import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from tdln_pd1.data import load_cells
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
