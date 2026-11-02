import pandas as pd
import pytest
from tdln_pd1.features import aggregate_cells


def test_aggregate_cells_builds_fractions_and_signatures():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2"],
        "cell_type": ["CD8 T", "DC", "CD8 T"],
        "PDCD1": [2.0, 0.0, 1.0],
        "GZMB": [4.0, 0.0, 2.0],
    })
    out = aggregate_cells(cells).set_index("patient_id")
    assert out.loc["p1", "frac_cd8_t"] == 0.5
    assert out.loc["p1", "sig_exhaustion"] == 1.0
    assert out.loc["p2", "sig_cytotoxic"] == 2.0


def test_required_columns_are_checked():
    with pytest.raises(ValueError):
        aggregate_cells(pd.DataFrame({"patient_id": ["p1"]}))


def test_tcr_clonality_and_cross_tissue_sharing():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 6,
        "cell_type": ["Tpex"] * 6,
        "clonotype_id": ["a", "a", "a", "b", "c", "d"],
        "tissue": ["TDLN", "TDLN", "tumour", "TDLN", "tumour", "TDLN"],
    })
    row = aggregate_cells(cells).iloc[0]
    assert row.tcr_clone_richness == 4
    assert row.tcr_top_clone_fraction == 0.5
    assert row.tcr_tdln_tumour_shared_fraction == pytest.approx(1 / 3)
    assert 0 < row.tcr_simpson_concentration < 1
