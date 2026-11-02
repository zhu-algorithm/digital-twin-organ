from dataclasses import asdict
import json

import pytest

from tdln_pd1.cli import _load_parameter_file, _training_input_provenance
from tdln_pd1.simulator import TwinParameters


def test_parameter_file_loading_is_exact_and_fingerprinted(tmp_path):
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(asdict(TwinParameters())))
    parameters, provenance = _load_parameter_file(str(path), TwinParameters)
    assert parameters == TwinParameters()
    assert len(provenance["sha256"]) == 64
    assert provenance["parameter_class"] == "TwinParameters"


def test_parameter_file_rejects_schema_drift_and_nonpositive_values(tmp_path):
    payload = asdict(TwinParameters())
    payload.pop("tumour_growth")
    path = tmp_path / "missing.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="schema mismatch"):
        _load_parameter_file(str(path), TwinParameters)
    payload = asdict(TwinParameters())
    payload["tumour_growth"] = 0
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="strictly positive"):
        _load_parameter_file(str(path), TwinParameters)


def test_training_input_binds_feature_and_parameter_sidecar_hashes(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,x\np0,1\n")
    sidecar = tmp_path / "hybrid.csv.parameters.json"
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "systemic": {"file": "fit.json", "sha256": "a" * 64,
                     "parameter_class": "TwinParameters"},
    }, "research_only": True}))
    provenance = _training_input_provenance(str(features))
    assert len(provenance["feature_sha256"]) == 64
    assert len(provenance["parameter_sidecar_sha256"]) == 64
    assert provenance["parameter_provenance"]["sources"]["systemic"]["sha256"] == "a" * 64
