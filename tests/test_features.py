import pandas as pd
import pytest
from tdln_pd1.features import aggregate_cells


def test_aggregate_cells_builds_fractions_and_signatures():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2"],
        "cell_type": ["CD8 T", "DC", "CD8 T"],
        "PDCD1": [2.0, 0.0, 1.0],
        "CTLA4": [2.0, 0.0, 1.0],
        "LAG3": [2.0, 0.0, 1.0],
        "NKG7": [4.0, 0.0, 2.0],
        "GZMB": [4.0, 0.0, 2.0],
    })
    out = aggregate_cells(cells).set_index("patient_id")
    assert out.loc["p1", "frac_cd8_t"] == 0.5
    assert out.loc["p1", "sig_exhaustion"] == 1.0
    assert out.loc["p2", "sig_cytotoxic"] == 2.0
    assert out.loc["p1", "sig_exhaustion_cd8_t"] == 2.0
    assert out.loc["p1", "sig_cytotoxic_cd8_t"] == 4.0
    assert out.loc["p2", "frac_dc"] == 0.0


def test_cell_fraction_matrix_uses_observed_zero_not_missing_imputation():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2", "p2"],
        "cell_type": ["CD8 T", "DC", "CD8 T", "Treg"],
        "tissue": ["TDLN"] * 4,
    })
    out = aggregate_cells(cells).set_index("patient_id")
    assert out.loc["p1", "frac_treg"] == 0.0
    assert out.loc["p2", "frac_dc"] == 0.0
    assert out.loc["p1", ["frac_cd8_t", "frac_dc", "frac_treg"]].sum() == 1.0
    assert out.loc["p2", ["frac_cd8_t", "frac_dc", "frac_treg"]].sum() == 1.0
    assert not out.filter(like="frac_").isna().any().any()


def test_locked_cell_type_vocabulary_stabilizes_external_feature_schema():
    vocabulary = ("CD8 T", "DC", "Treg")
    development = pd.DataFrame({
        "patient_id": ["p1", "p1", "p1"],
        "cell_type": list(vocabulary),
    })
    external = pd.DataFrame({
        "patient_id": ["e1", "e1"], "cell_type": ["CD8 T", "CD8 T"],
    })
    dev = aggregate_cells(
        development, feature_tissue_pattern="all", cell_type_vocabulary=vocabulary,
    )
    ext = aggregate_cells(
        external, feature_tissue_pattern="all", cell_type_vocabulary=vocabulary,
    )
    fraction_columns = [column for column in dev if column.startswith("frac_")]
    assert fraction_columns == ["frac_cd8_t", "frac_dc", "frac_treg"]
    assert [column for column in ext if column.startswith("frac_")] == fraction_columns
    assert ext.loc[0, "frac_dc"] == 0.0
    assert ext.loc[0, "frac_treg"] == 0.0
    unknown = external.assign(cell_type=["CD8 T", "B cell"])
    with pytest.raises(ValueError, match="outside the locked vocabulary"):
        aggregate_cells(
            unknown, feature_tissue_pattern="all", cell_type_vocabulary=vocabulary,
        )


def test_lineage_restricted_signatures_separate_state_from_composition():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 4,
        "cell_type": ["CD8 T", "CD8 T", "DC", "Treg"],
        "PDCD1": [1.0, 3.0, 100.0, 50.0],
        "CTLA4": [1.0, 3.0, 100.0, 50.0],
        "LAG3": [1.0, 3.0, 100.0, 50.0],
        "NKG7": [2.0, 4.0, 80.0, 0.0],
        "GZMB": [2.0, 4.0, 80.0, 0.0],
        "TCF7": [5.0, 7.0, 90.0, 1.0],
        "CCR7": [5.0, 7.0, 90.0, 1.0],
        "FOXP3": [0.0, 0.0, 0.0, 8.0],
        "IL2RA": [0.0, 0.0, 0.0, 8.0],
        "CD80": [0.0, 0.0, 9.0, 0.0],
        "CD86": [0.0, 0.0, 9.0, 0.0],
    })
    row = aggregate_cells(cells, feature_tissue_pattern="all").iloc[0]
    assert row.sig_exhaustion_cd8_t == 2.0
    assert row.sig_cytotoxic_cd8_t == 3.0
    assert row.sig_tcf7_memory_cd8_t == 6.0
    assert row.sig_treg_treg == 22.0
    assert row.sig_dendritic_activation_dc == 36.0
    assert row.sig_exhaustion != row.sig_exhaustion_cd8_t


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
    assert row.tcr_tumour_tdln_shared_cell_fraction == pytest.approx(1 / 2)
    assert row.tcr_tdln_tumour_jaccard == pytest.approx(1 / 4)
    assert row.tcr_tdln_clone_richness == 3
    assert row.tcr_tumour_clone_richness == 2
    assert row.tcr_tdln_top_clone_fraction == 0.5
    assert row.tcr_tumour_top_clone_fraction == 0.5
    assert row.tcr_tdln_tumour_morisita_horn == pytest.approx(4 / 7)
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
    assert row.tcr_tumour_tdln_shared_cell_fraction == 0.0
    assert row.tcr_tdln_tumour_jaccard == 0.0
    assert row.tcr_tdln_tumour_morisita_horn == 0.0


def test_temporal_input_requires_explicit_baseline_and_excludes_post_treatment():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2", "p2"],
        "cell_type": ["Tpex"] * 4,
        "timepoint": ["baseline", "post", "baseline", "post"],
        "TCF7": [1.0, 100.0, 2.0, 200.0],
        "CCR7": [1.0, 100.0, 2.0, 200.0],
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
        "CCR7": [2.0, 4.0, 100.0, 100.0],
        "clonotype_id": ["a", "b", "a", "c"],
    })
    row = aggregate_cells(cells).iloc[0]
    assert row.sig_tcf7_memory == 3.0
    assert row.tcr_tdln_tumour_shared_fraction == pytest.approx(0.5)


def test_inferred_composition_vocabulary_excludes_tumour_only_annotations():
    cells = pd.DataFrame({
        "patient_id": ["p1"] * 4,
        "cell_type": ["Tpex", "DC", "Tumour-infiltrating T", "Malignant"],
        "tissue": ["TDLN", "TDLN", "tumour", "tumour"],
        "clonotype_id": ["a", None, "a", None],
    })
    row = aggregate_cells(cells).iloc[0]
    assert row.frac_tpex == 0.5
    assert row.frac_dc == 0.5
    assert "frac_tumour_infiltrating_t" not in row.index
    assert "frac_malignant" not in row.index
    assert row.tcr_tdln_tumour_shared_fraction == 1.0


def test_each_patient_must_have_selected_baseline_and_tdln_cells():
    temporal = pd.DataFrame({"patient_id": ["p1", "p2"], "cell_type": ["T", "T"],
                             "timepoint": ["baseline", "post"]})
    with pytest.raises(ValueError, match="removed all cells"):
        aggregate_cells(temporal, "timepoint", ("baseline",))
    tissue = pd.DataFrame({"patient_id": ["p1"], "cell_type": ["T"], "tissue": ["tumour"]})
    with pytest.raises(ValueError, match="no cells matching"):
        aggregate_cells(tissue)


def test_expression_contract_rejects_corrupt_values_but_preserves_missingness():
    valid = pd.DataFrame({
        "patient_id": ["p1", "p1"], "cell_type": ["T", "T"],
        "NKG7": [1.0, float("nan")], "GZMB": [1.0, float("nan")],
    })
    row = aggregate_cells(valid, feature_tissue_pattern="all").iloc[0]
    assert row.sig_cytotoxic == 1.0

    non_numeric = valid.assign(GZMB=pd.Series([1.0, "bad"], dtype=object))
    with pytest.raises(ValueError, match="numeric or missing: GZMB"):
        aggregate_cells(non_numeric, feature_tissue_pattern="all")
    infinite = valid.copy()
    infinite.loc[1, "GZMB"] = float("inf")
    with pytest.raises(ValueError, match="finite or missing: GZMB"):
        aggregate_cells(infinite, feature_tissue_pattern="all")


def test_incomplete_signature_gene_panels_are_explicitly_missing_and_audited():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1"], "cell_type": ["CD8 T", "CD8 T"],
        "GZMB": [2.0, 4.0],
    })
    row = aggregate_cells(cells, feature_tissue_pattern="all").iloc[0]
    assert pd.isna(row.sig_cytotoxic)
    assert pd.isna(row.sig_cytotoxic_cd8_t)
    assert row.qc_sig_cytotoxic_gene_coverage == 0.25
    permissive = aggregate_cells(
        cells, feature_tissue_pattern="all", min_signature_gene_coverage=0.25,
    ).iloc[0]
    assert permissive.sig_cytotoxic == 3.0
    assert permissive.qc_min_signature_gene_coverage == 0.25
    for invalid in (0, 1.1, float("nan"), "invalid"):
        with pytest.raises(ValueError, match="minimum signature gene coverage"):
            aggregate_cells(
                cells, feature_tissue_pattern="all",
                min_signature_gene_coverage=invalid,
            )


def test_cell_type_fraction_names_cannot_silently_collide():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1"],
        "cell_type": ["CD8 T", "CD8/T"],
    })
    with pytest.raises(ValueError, match="collide"):
        aggregate_cells(cells, feature_tissue_pattern="all")
    with pytest.raises(ValueError, match="non-missing"):
        aggregate_cells(
            pd.DataFrame({"patient_id": ["p1"], "cell_type": [None]}),
            feature_tissue_pattern="all",
        )
