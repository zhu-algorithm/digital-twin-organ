from dataclasses import asdict
import hashlib
import json

import pytest

from tdln_pd1.cli import _load_parameter_file, _training_input_provenance
from tdln_pd1.hybrid import FeatureStateMapping
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


def test_feature_state_mapping_file_is_exact_and_fingerprinted(tmp_path):
    path = tmp_path / "mapping.json"
    path.write_text(json.dumps(asdict(FeatureStateMapping())))
    mapping, provenance = _load_parameter_file(str(path), FeatureStateMapping)
    assert mapping == FeatureStateMapping()
    assert provenance["parameter_class"] == "FeatureStateMapping"
    assert len(provenance["sha256"]) == 64


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


def test_training_input_accepts_and_validates_mapping_provenance(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,x\np0,1\n")
    sidecar = tmp_path / "hybrid.csv.parameters.json"
    source = {"file": "mapping.json", "sha256": "b" * 64,
              "parameter_class": "FeatureStateMapping"}
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "feature_state_mapping": source,
    }, "research_only": True}))
    provenance = _training_input_provenance(str(features))
    assert provenance["parameter_provenance"]["sources"]["feature_state_mapping"] == source
    source["parameter_class"] = "TwinParameters"
    sidecar.write_text(json.dumps({"schema_version": 1, "sources": {
        "feature_state_mapping": source,
    }, "research_only": True}))
    with pytest.raises(ValueError, match="invalid parameter source fingerprint"):
        _training_input_provenance(str(features))


def test_training_input_binds_and_checks_label_free_augmentation_audit(tmp_path):
    features = tmp_path / "hybrid.csv"
    features.write_text("patient_id,x\np0,1\n")
    feature_sha = hashlib.sha256(features.read_bytes()).hexdigest()
    audit = tmp_path / "hybrid.csv.twin-audit.json"
    audit.write_text(json.dumps({
        "schema_version": 1,
        "stage": "label_free_twin_augmentation",
        "labels_used": False,
        "output": {"file": features.name, "sha256": feature_sha},
        "inputs": {},
    }))
    provenance = _training_input_provenance(str(features))
    assert provenance["augmentation_audit"]["labels_used"] is False
    assert len(provenance["augmentation_audit_sha256"]) == 64

    features.write_text("patient_id,x\np0,2\n")
    with pytest.raises(ValueError, match="mismatched feature augmentation audit"):
        _training_input_provenance(str(features))
