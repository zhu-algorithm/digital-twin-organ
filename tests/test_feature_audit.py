import json

import pytest

from tdln_pd1.feature_audit import verify_feature_audit, write_feature_audit


def test_feature_audit_binds_cells_settings_and_output(tmp_path):
    cells = tmp_path / "cells.csv"
    cells.write_text("patient_id,cell_type\np1,CD8 T\n")
    output = tmp_path / "features.csv"
    output.write_text("patient_id,sig\np1,1\n")
    settings = {
        "time_column": "visit", "baseline_values": ["baseline"],
        "feature_tissue_pattern": "tdln|lymph",
        "minimum_signature_gene_coverage": 0.75, "labels_used": False,
    }
    audit = write_feature_audit(cells, output, settings)
    verified = verify_feature_audit(audit, cells, output)
    assert verified["valid"] is True
    assert verified["payload"]["settings"] == settings
    assert verified["payload"]["labels_used"] is False
    assert len(verified["audit_sha256"]) == 64


def test_feature_audit_detects_input_output_and_manifest_replacement(tmp_path):
    cells = tmp_path / "cells.csv"
    cells.write_text("patient_id,cell_type\np1,T\n")
    output = tmp_path / "features.csv"
    output.write_text("patient_id,x\np1,1\n")
    audit = write_feature_audit(cells, output, {})
    cells.write_text("patient_id,cell_type\np1,DC\n")
    with pytest.raises(ValueError, match="cell input SHA-256 mismatch"):
        verify_feature_audit(audit, cells, output)
    cells.write_text("patient_id,cell_type\np1,T\n")
    output.write_text("patient_id,x\np1,2\n")
    with pytest.raises(ValueError, match="feature output SHA-256 mismatch"):
        verify_feature_audit(audit, cells, output)
    output.write_text("patient_id,x\np1,1\n")
    payload = json.loads(audit.read_text())
    payload["labels_used"] = True
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="non-label-free"):
        verify_feature_audit(audit, cells, output)
