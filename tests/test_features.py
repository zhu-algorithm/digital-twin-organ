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
    assert row.tcr_tdln_tumour_shared_cell_fraction == pytest.approx(1 / 2)
    assert row.tcr_tdln_tumour_jaccard == pytest.approx(1 / 4)
    assert 0 < row.tcr_simpson_concentration < 1
    assert row.n_cells == 4
    assert row.n_cells_total_selected == 6


def test_unassigned_clonotype_placeholders_never_create_false_sharing():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 8,
        "cell_type": ["T"] * 8,
        "tissue": ["TDLN"] * 4 + ["tumour"] * 4,
        "clonotype_id": ["a", "a", "unassigned", "NA",
                         "b", "None", "nan", "N/A"],
    })
    row = aggregate_cells(cells).iloc[0]
    assert row.tcr_clone_richness == 2
    assert row.tcr_tdln_tumour_shared_fraction == 0.0
    assert row.tcr_tdln_tumour_shared_cell_fraction == 0.0
    assert row.tcr_tdln_tumour_jaccard == 0.0


def test_temporal_input_requires_explicit_baseline_and_excludes_post_treatment():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2", "p2"],
        "cell_type": ["Tpex"] * 4,
        "timepoint": ["baseline", "post", "baseline", "post"],
        "TCF7": [1.0, 100.0, 2.0, 200.0],
    })
    with pytest.raises(ValueError, match="temporal columns detected"):
        aggregate_cells(cells)
    out = aggregate_cells(cells, "timepoint", ("baseline",)).set_index("patient_id")
    assert out.loc["p1", "sig_tcf7_memory"] == 1.0
    assert out.loc["p2", "sig_tcf7_memory"] == 2.0
    assert out.n_cells.eq(1).all()


def test_signatures_use_tdln_while_tcr_sharing_uses_all_tissues():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 4,
        "cell_type": ["Tpex"] * 4,
        "tissue": ["TDLN", "TDLN", "tumour", "tumour"],
        "TCF7": [2.0, 4.0, 100.0, 100.0],
        "clonotype_id": ["a", "b", "a", "c"],
    })
    row = aggregate_cells(cells).iloc[0]
    assert row.sig_tcf7_memory == 3.0
    assert row.tcr_tdln_tumour_shared_fraction == pytest.approx(0.5)


def test_each_patient_must_have_selected_baseline_and_tdln_cells():
    temporal = pd.DataFrame({"patient_id": ["p1", "p2"], "cell_type": ["T", "T"],
                             "timepoint": ["baseline", "post"]})
    with pytest.raises(ValueError, match="removed all cells"):
        aggregate_cells(temporal, "timepoint", ("baseline",))
    tissue = pd.DataFrame({"patient_id": ["p1"], "cell_type": ["T"], "tissue": ["tumour"]})
    with pytest.raises(ValueError, match="no cells matching"):
        aggregate_cells(tissue)
