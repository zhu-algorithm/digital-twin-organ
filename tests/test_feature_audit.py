import json

import pandas as pd
import pytest

from tdln_pd1.feature_audit import verify_feature_audit, write_feature_audit
from tdln_pd1.features import aggregate_cells


DEFAULT_SETTINGS = {
    "time_column": None, "baseline_values": [],
    "feature_tissue_pattern": "tdln|lymph",
    "minimum_signature_gene_coverage": 0.5,
    "gene_symbol_column": None, "expression_layer": None,
    "labels_used": False,
}


def _write_cells_and_features(cells, output, cell_type="CD8 T"):
    frame = pd.DataFrame({"patient_id": ["p1"], "cell_type": [cell_type]})
    frame.to_csv(cells, index=False)
    aggregate_cells(frame).to_csv(output, index=False)


def test_feature_audit_binds_cells_settings_and_output(tmp_path):
    cells = tmp_path / "cells.csv"
    frame = pd.DataFrame({
        "patient_id": ["p1"], "cell_type": ["CD8 T"], "visit": ["baseline"],
    })
    frame.to_csv(cells, index=False)
    output = tmp_path / "features.csv"
    settings = {
        "time_column": "visit", "baseline_values": ["baseline"],
        "feature_tissue_pattern": "tdln|lymph",
        "minimum_signature_gene_coverage": 0.75, "labels_used": False,
    }
    aggregate_cells(frame, "visit", ("baseline",), "tdln|lymph", 0.75).to_csv(
        output, index=False
    )
    audit = write_feature_audit(cells, output, settings)
    verified = verify_feature_audit(audit, cells, output)
    assert verified["valid"] is True
    assert verified["payload"]["settings"] == settings
    assert verified["payload"]["labels_used"] is False
    assert verified["derivation_replayed"] is True
    assert len(verified["audit_sha256"]) == 64


def test_feature_audit_detects_input_output_and_manifest_replacement(tmp_path):
    cells = tmp_path / "cells.csv"
    output = tmp_path / "features.csv"
    _write_cells_and_features(cells, output, "T")
    audit = write_feature_audit(cells, output, DEFAULT_SETTINGS)
    replacement = tmp_path / "replacement.csv"
    _write_cells_and_features(cells, replacement, "DC")
    with pytest.raises(ValueError, match="cell input SHA-256 mismatch"):
        verify_feature_audit(audit, cells, output)
    _write_cells_and_features(cells, output, "T")
    output.write_text("patient_id,x\np1,2\n")
    with pytest.raises(ValueError, match="feature output SHA-256 mismatch"):
        verify_feature_audit(audit, cells, output)
    _write_cells_and_features(cells, output, "T")
    payload = json.loads(audit.read_text())
    payload["labels_used"] = True
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="non-label-free"):
        verify_feature_audit(audit, cells, output)


def test_feature_audit_replay_detects_setting_tampering(tmp_path):
    cells, output = tmp_path / "cells.csv", tmp_path / "features.csv"
    _write_cells_and_features(cells, output)
    audit = write_feature_audit(cells, output, DEFAULT_SETTINGS)
    payload = json.loads(audit.read_text())
    payload["settings"]["minimum_signature_gene_coverage"] = 1.0
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="derivation replay mismatch"):
        verify_feature_audit(audit, cells, output)


def test_feature_audit_rejects_invalid_anndata_matrix_contract(tmp_path):
    cells, output = tmp_path / "cells.csv", tmp_path / "features.csv"
    _write_cells_and_features(cells, output)
    audit = write_feature_audit(cells, output, DEFAULT_SETTINGS)
    payload = json.loads(audit.read_text())
    payload["settings"]["expression_layer"] = ""
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="expression_layer"):
        verify_feature_audit(audit, cells, output)
