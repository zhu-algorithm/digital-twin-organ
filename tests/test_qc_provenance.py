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
