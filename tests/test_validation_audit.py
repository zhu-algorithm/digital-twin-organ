import json
import hashlib
from pathlib import Path

import joblib
import pandas as pd
import pytest

from tdln_pd1.validation_audit import (verify_evaluation_audit, verify_prediction_audit,
                                       write_evaluation_audit,
                                       write_prediction_audit,
                                       validate_prediction_feature_lineage)


def _files(tmp_path: Path):
    model = tmp_path / "model.joblib"
    joblib.dump({
        "artifact_version": 9,
        "training_cohort_fingerprint": "a" * 64,
        "training_feature_cohort_fingerprint": "b" * 64,
        "point_prediction_method": "full_cohort_locked_model",
    }, model)
    features = tmp_path / "features.csv"
    features.write_text("patient_id,x\np1,1\n")
    predictions = tmp_path / "predictions.csv"
    predictions.write_text("patient_id,probability\np1,0.7\n")
    return model, features, predictions


def _add_treatment_lineage(model: Path, features: Path, response="RECIST 1.1"):
    contract = {"agent_class": "anti_pd1", "dose_scale": "arbitrary_model_exposure",
                "response_definition": response, "regimen_ids": ["q3w"]}
    perfusion = {"scale": "dimensionless_multiplier",
                 "measurement_definition": "duplex_ratio",
                 "normalization_reference": "population_baseline", "patient_count": 1}
    artifact = joblib.load(model)
    artifact["training_input_provenance"] = {"augmentation_audit": {"inputs": {
        "features": {"derivation_payload": {"treatment_contract": contract,
                                               "perfusion_contract": perfusion}}
    }}}
    joblib.dump(artifact, model)
    payload = {
        "schema_version": 1, "stage": "label_free_twin_augmentation",
        "labels_used": False,
        "inputs": {"features": {"derivation_payload": {
            "stage": "label_free_liver_tdln_counterfactual_augmentation",
            "labels_used": False, "treatment_contract": contract,
            "perfusion_contract": perfusion,
        }}},
        "output": {"file": features.name,
                   "sha256": hashlib.sha256(features.read_bytes()).hexdigest()},
    }
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    return payload


def _add_mechanistic_contract(model: Path, payload: dict):
    contract = {
        "systemic_parameters": {"sha256": "0" * 64},
        "systemic_diagnostics": {"sha256": "9" * 64},
        "spatial_parameters": {"sha256": "8" * 64},
        "spatial_diagnostics": {"sha256": "7" * 64},
        "bridge_parameters": {"sha256": "1" * 64},
        "bridge_diagnostics": {"sha256": "2" * 64},
        "feature_bridge_mapping": {"sha256": "3" * 64},
        "feature_bridge_diagnostics": {"sha256": "4" * 64},
        "mapping": {"killing_cytotoxic_coefficient": 1.4},
        "base_parameters": {"coupling": {"immune_killing_gain": 0.004}},
        "mechanistic_feature_columns": ["liver_cf_tumour_relative_reduction"],
        "settings": {"days": 28.0, "step": 0.25, "paired_control": True,
                     "uncertainty_samples": 100, "uncertainty_seed": 42},
    }
    artifact = joblib.load(model)
    training = artifact["training_input_provenance"]["augmentation_audit"][
        "inputs"]["features"]["derivation_payload"]
    training.update(json.loads(json.dumps(contract)))
    joblib.dump(artifact, model)
    external = payload["inputs"]["features"]["derivation_payload"]
    external.update(json.loads(json.dumps(contract)))
    return external


def test_treatment_aware_prediction_requires_compatible_feature_lineage(tmp_path: Path):
    model, features, _ = _files(tmp_path)
    payload = _add_treatment_lineage(model, features)
    result = validate_prediction_feature_lineage(model, features)
    assert result["required"] is result["verified"] is True
    assert result["compatible_semantics"]["response_definition"] == "RECIST 1.1"
    assert result["compatible_perfusion_semantics"]["measurement_definition"] == (
        "duplex_ratio")

    Path(f"{features}.twin-audit.json").unlink()
    with pytest.raises(ValueError, match="requires a feature twin audit"):
        validate_prediction_feature_lineage(model, features)

    audit = Path(f"{features}.twin-audit.json")
    payload["inputs"]["features"]["derivation_payload"]["treatment_contract"][
        "response_definition"] = "iRECIST"
    audit.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="response_definition"):
        validate_prediction_feature_lineage(model, features)


def test_prediction_requires_compatible_perfusion_semantics(tmp_path: Path):
    model, features, _ = _files(tmp_path)
    payload = _add_treatment_lineage(model, features)
    derivation = payload["inputs"]["features"]["derivation_payload"]
    derivation["perfusion_contract"]["normalization_reference"] = "liver_volume"
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="normalization_reference"):
        validate_prediction_feature_lineage(model, features)
    del derivation["perfusion_contract"]
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="lack a perfusion contract"):
        validate_prediction_feature_lineage(model, features)


def test_prediction_rejects_treatment_schedule_extrapolation(tmp_path: Path):
    model, features, _ = _files(tmp_path)
    payload = _add_treatment_lineage(model, features)
    support = {
        "model_dose_amount": [0.5, 1.0], "dose_count": [1, 3],
        "total_model_exposure": [1.0, 2.0], "first_dose_day": [0.0, 1.0],
        "last_dose_day": [0.0, 14.0],
    }
    artifact = joblib.load(model)
    training = artifact["training_input_provenance"]["augmentation_audit"][
        "inputs"]["features"]["derivation_payload"]["treatment_contract"]
    training["schedule_support"] = json.loads(json.dumps(support))
    joblib.dump(artifact, model)
    external = payload["inputs"]["features"]["derivation_payload"][
        "treatment_contract"]
    external["schedule_support"] = json.loads(json.dumps(support))
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    result = validate_prediction_feature_lineage(model, features)
    assert result["compatible_treatment_schedule_support"] == support
    external["schedule_support"]["total_model_exposure"] = [1.0, 3.0]
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="total_model_exposure"):
        validate_prediction_feature_lineage(model, features)


@pytest.mark.parametrize("mutation,field", [
    (lambda contract: contract["systemic_parameters"].update(sha256="e" * 64),
     "systemic_parameters_sha256"),
    (lambda contract: contract["systemic_diagnostics"].update(sha256="d" * 64),
     "systemic_diagnostics_sha256"),
    (lambda contract: contract["spatial_parameters"].update(sha256="c" * 64),
     "spatial_parameters_sha256"),
    (lambda contract: contract["spatial_diagnostics"].update(sha256="b" * 64),
     "spatial_diagnostics_sha256"),
    (lambda contract: contract["feature_bridge_mapping"].update(sha256="f" * 64),
     "feature_bridge_mapping_sha256"),
    (lambda contract: contract["settings"].update(days=21.0), "settings"),
    (lambda contract: contract["mapping"].update(killing_cytotoxic_coefficient=0.8),
     "mapping"),
    (lambda contract: contract["base_parameters"]["coupling"].update(
        immune_killing_gain=0.02), "base_parameters"),
])
def test_prediction_locks_mechanistic_feature_semantics(tmp_path: Path, mutation, field):
    model, features, _ = _files(tmp_path)
    payload = _add_treatment_lineage(model, features)
    external = _add_mechanistic_contract(model, payload)
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    verified = validate_prediction_feature_lineage(model, features)
    assert verified["compatible_mechanistic_feature_contract"][
        "feature_bridge_mapping_sha256"] == "3" * 64
    mutation(external)
    Path(f"{features}.twin-audit.json").write_text(json.dumps(payload))
    with pytest.raises(ValueError, match=field):
        validate_prediction_feature_lineage(model, features)


def test_prediction_lineage_rejects_feature_or_audit_tampering(tmp_path: Path):
    model, features, _ = _files(tmp_path)
    _add_treatment_lineage(model, features)
    features.write_text("patient_id,x\np1,2\n")
    with pytest.raises(ValueError, match="mismatched prediction feature twin audit"):
        validate_prediction_feature_lineage(model, features)


def test_prediction_and_evaluation_audits_bind_every_file(tmp_path: Path, monkeypatch):
    model, features, predictions = _files(tmp_path)
    monkeypatch.setattr(
        "tdln_pd1.validation_audit.predict_locked",
        lambda *_: pd.DataFrame({"patient_id": ["p1"], "probability": [0.7]}),
    )
    audit_path = write_prediction_audit(model, features, predictions, 1)
    verified = verify_prediction_audit(predictions, model_path=model, features_path=features)
    assert verified["payload"]["model_refit"] is False
    assert verified["payload"]["patient_count"] == 1
    assert verified["payload"]["decision_threshold"] == 0.5
    assert verified["payload"]["schema_version"] == 3
    assert verified["payload"]["point_prediction_method"] == "full_cohort_locked_model"
    assert verified["derivation_replayed"] is True

    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,response\np1,1\n")
    metrics = tmp_path / "metrics.json"
    metrics.write_text('{"roc_auc": 0.5}\n')
    evaluation_audit = write_evaluation_audit(
        predictions, verified, labels, metrics, bootstrap=500,
    )
    payload = json.loads(evaluation_audit.read_text())
    assert payload["stage"] == "locked_external_evaluation"
    assert payload["schema_version"] == 2
    assert payload["settings"]["bootstrap_draws"] == 500
    assert payload["inputs"]["prediction_audit"]["file"] == audit_path.name


def test_prediction_audit_rejects_posthoc_file_replacement(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    write_prediction_audit(model, features, predictions, 1)
    predictions.write_text("patient_id,probability\np1,0.1\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify_prediction_audit(predictions, model_path=model, features_path=features)


def test_prediction_audit_requires_versioned_training_fingerprints(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    joblib.dump({"artifact_version": 1}, model)
    with pytest.raises(ValueError, match="fingerprints"):
        write_prediction_audit(model, features, predictions, 1)


def test_prediction_audit_requires_replay_inputs(tmp_path: Path):
    model, features, predictions = _files(tmp_path)
    write_prediction_audit(model, features, predictions, 1)
    with pytest.raises(ValueError, match="requires the locked model"):
        verify_prediction_audit(predictions)


def test_evaluation_audit_replays_exact_metric_report(tmp_path: Path, monkeypatch):
    model, features, predictions = _files(tmp_path)
    replayed_predictions = pd.DataFrame({"patient_id": ["p1"], "probability": [0.7]})
    monkeypatch.setattr("tdln_pd1.validation_audit.predict_locked", lambda *_: replayed_predictions)
    prediction_audit = write_prediction_audit(model, features, predictions, 1)
    verified = verify_prediction_audit(predictions, model_path=model, features_path=features)
    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,response\np1,1\n")
    expected = {"roc_auc": 0.5, "validation_chain": {
        "prediction_audit_sha256": verified["audit_sha256"],
        "model_sha256": verified["payload"]["inputs"]["model"]["sha256"],
        "features_sha256": verified["payload"]["inputs"]["features"]["sha256"],
        "integrity_verified": True,
        "prediction_derivation_replayed": True,
    }}
    monkeypatch.setattr("tdln_pd1.validation_audit.evaluate_external", lambda *_: {"roc_auc": 0.5})
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(expected, indent=2) + "\n")
    audit = write_evaluation_audit(predictions, verified, labels, metrics, bootstrap=500)
    report = verify_evaluation_audit(
        audit, predictions, model, features, labels, metrics, prediction_audit,
    )
    assert report["valid"] is True
    assert report["evaluation_derivation_replayed"] is True

    metrics.write_text(json.dumps({"roc_auc": 0.9}, indent=2) + "\n")
    with pytest.raises(ValueError, match="SHA-256"):
        verify_evaluation_audit(
            audit, predictions, model, features, labels, metrics, prediction_audit,
        )

    forged = json.loads(audit.read_text())
    forged["output"]["metrics"]["sha256"] = hashlib.sha256(metrics.read_bytes()).hexdigest()
    audit.write_text(json.dumps(forged, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="derivation replay"):
        verify_evaluation_audit(
            audit, predictions, model, features, labels, metrics, prediction_audit,
        )
