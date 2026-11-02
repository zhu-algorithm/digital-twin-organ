import hashlib
import json
from pathlib import Path

import pandas as pd

from tdln_pd1.provenance import verify_manifest
from tdln_pd1.qc import cohort_qc


def test_qc_reports_pairing_clonotypes_and_response_balance():
    cells = pd.DataFrame({
        "patient_id": ["a", "a", "b", "b", "c", "c", "d", "d"],
        "cell_type": ["T cell"] * 8, "tissue": ["TDLN", "tumour"] * 4,
        "clonotype_id": ["x", "x", "y", None, "z", "z", "w", "w"],
        "PDCD1": [1.0] * 8, "TCF7": [1.0] * 8, "GZMB": [1.0] * 8,
    })
    labels = pd.DataFrame({"patient_id": ["a", "b", "c", "d"], "response": [0, 0, 1, 1]})
    report = cohort_qc(cells, labels)
    assert report["passed"]
    assert report["tissue_pairing"]["paired_tdln_tumour_patients"] == 4
    assert report["clonotypes"]["assigned_fraction"] > 0
    assert report["feature_derivation"] == {
        "validated": True, "n_patient_rows": 4, "shared_with_featurize": True,
    }
    assert report["minimum_signature_gene_coverage"] == 0.5
    assert report["signature_gene_coverage"]["cytotoxic"]["fraction"] == 0.25
    assert report["signature_gene_coverage"]["cytotoxic"]["usable"] is False
    assert any("emitted as missing" in warning for warning in report["warnings"])


def test_qc_clonotype_assignment_uses_same_placeholder_filter_as_features():
    cells = pd.DataFrame({
        "patient_id": ["p1", "p1", "p2", "p2"],
        "cell_type": ["T"] * 4,
        "clonotype_id": ["a", "NA", "b", "unassigned"],
    })
    labels = pd.DataFrame({"patient_id": ["p1", "p2"], "response": [0, 1]})
    report = cohort_qc(cells, labels, feature_tissue_pattern="all")
    assert report["clonotypes"]["assigned_fraction"] == 0.5


def test_qc_and_featurization_share_configurable_signature_coverage():
    cells = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"], "cell_type": ["CD8 T"] * 4,
        "GZMB": [1.0, 2.0, 3.0, 4.0],
    })
    labels = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"], "response": [0, 0, 1, 1],
    })
    report = cohort_qc(
        cells, labels, feature_tissue_pattern="all",
        min_signature_gene_coverage=0.25,
    )
    assert report["passed"]
    assert report["minimum_signature_gene_coverage"] == 0.25
    assert report["signature_gene_coverage"]["cytotoxic"]["usable"] is True


def test_qc_fails_on_missing_outcome_or_patient_identity():
    cells = pd.DataFrame({"patient_id": ["p1", "p2"], "cell_type": ["T", "T"]})
    missing_outcome = pd.DataFrame({"patient_id": ["p1", "p2"],
                                    "response": [0, None]})
    report = cohort_qc(cells, missing_outcome)
    assert not report["passed"]
    assert any("complete binary" in error for error in report["errors"])

    missing_patient = cells.copy()
    missing_patient.loc[0, "patient_id"] = None
    labels = pd.DataFrame({"patient_id": ["p1", "p2"], "response": [0, 1]})
    report = cohort_qc(missing_patient, labels)
    assert not report["passed"]
    assert any("non-missing" in error for error in report["errors"])


def test_qc_fails_on_ambiguous_timepoints_and_reports_explicit_baseline_selection():
    cells = pd.DataFrame({
        "patient_id": ["a", "a", "b", "b", "c", "c", "d", "d"],
        "cell_type": ["T cell"] * 8,
        "timepoint": ["baseline", "post"] * 4,
        "PDCD1": [1.0] * 8, "TCF7": [1.0] * 8, "GZMB": [1.0] * 8,
    })
    labels = pd.DataFrame({"patient_id": ["a", "b", "c", "d"], "response": [0, 0, 1, 1]})
    rejected = cohort_qc(cells, labels)
    assert not rejected["passed"]
    assert "temporal columns detected" in rejected["errors"][0]
    selected = cohort_qc(cells, labels, "timepoint", ("baseline",))
    assert selected["passed"]
    assert selected["n_cells_raw"] == 8
    assert selected["n_cells_selected"] == 4
    assert selected["temporal_selection"] == {
        "time_column": "timepoint", "baseline_values": ["baseline"],
    }


def test_qc_and_featurization_share_the_tdln_compartment_gate():
    cells = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"], "cell_type": ["T cell"] * 4,
        "tissue": ["TDLN", "TDLN", "tumour", "TDLN"],
        "PDCD1": [1.0] * 4, "TCF7": [1.0] * 4, "GZMB": [1.0] * 4,
    })
    labels = pd.DataFrame({"patient_id": ["a", "b", "c", "d"], "response": [0, 0, 1, 1]})
    report = cohort_qc(cells, labels)
    assert not report["passed"]
    assert "patients without feature-compartment cells" in report["errors"][0]
    all_tissues = cohort_qc(cells, labels, feature_tissue_pattern="all")
    assert all_tissues["passed"]
    assert all_tissues["feature_compartment"]["n_cells"] == 4


def test_qc_executes_expression_and_cell_type_feature_contracts():
    labels = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"], "response": [0, 0, 1, 1],
    })
    invalid_expression = pd.DataFrame({
        "patient_id": ["a", "b", "c", "d"], "cell_type": ["T"] * 4,
        "GZMB": pd.Series([1.0, 2.0, "bad", 4.0], dtype=object),
    })
    report = cohort_qc(invalid_expression, labels, feature_tissue_pattern="all")
    assert not report["passed"]
    assert any("feature derivation failed" in error and "numeric" in error
               for error in report["errors"])
    assert report["feature_derivation"]["validated"] is False

    colliding_types = pd.DataFrame({
        "patient_id": ["a", "a", "b", "c", "d"],
        "cell_type": ["CD8 T", "CD8/T", "T", "T", "T"],
    })
    report = cohort_qc(colliding_types, labels, feature_tissue_pattern="all")
    assert not report["passed"]
    assert any("collide" in error for error in report["errors"])


def test_manifest_hash_verification_detects_tampering(tmp_path: Path):
    data = tmp_path / "cells.csv"
    data.write_text("patient_id,cell_type\na,T cell\n")
    digest = hashlib.sha256(data.read_bytes()).hexdigest()
    manifest = tmp_path / "manifest.json"
    payload = {"cohort_id": "demo", "source": "synthetic", "license": "CC0",
               "modalities": ["scRNA-seq"], "outcome_definition": "synthetic binary",
               "files": [{"path": "cells.csv", "sha256": digest, "role": "cells"}]}
    manifest.write_text(json.dumps(payload))
    assert verify_manifest(manifest)["valid"]
    data.write_text("tampered")
    result = verify_manifest(manifest)
    assert not result["valid"]
    assert "hash mismatch" in result["errors"][0]
