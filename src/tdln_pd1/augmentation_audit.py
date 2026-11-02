"""Integrity chain for label-free mechanistic feature augmentation."""

from __future__ import annotations

from dataclasses import fields, replace
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from .assimilation import calibration_parameter_draws
from .hybrid import FeatureStateMapping, augment_with_twin
from .provenance import sha256_file
from .simulator import TwinParameters
from .spatial_dynamics import SpatialParameters


def _load_exact_parameters(path: str | Path | None, parameter_class):
    if path is None:
        return None
    payload = json.loads(Path(path).read_text())
    expected = {field.name for field in fields(parameter_class)}
    if not isinstance(payload, dict) or set(payload) != expected:
        raise ValueError(f"augmentation replay parameter schema mismatch: {parameter_class.__name__}")
    try:
        values = {name: float(payload[name]) for name in expected}
    except (TypeError, ValueError) as error:
        raise ValueError("augmentation replay parameters must be numeric") from error
    array = np.asarray(list(values.values()), dtype=float)
    if not np.isfinite(array).all() or (array <= 0).any():
        raise ValueError("augmentation replay parameters must be finite and positive")
    return parameter_class(**values)


def write_augmentation_audit(
    output_path: str | Path,
    inputs: dict,
    settings: dict,
) -> Path:
    output = Path(output_path)
    if not output.is_file():
        raise ValueError("augmented feature output must exist before audit")
    payload = {
        "schema_version": 1,
        "stage": "label_free_twin_augmentation",
        "research_only": True,
        "labels_used": False,
        "inputs": inputs,
        "settings": settings,
        "output": {"file": output.name, "sha256": sha256_file(output)},
        "limitation": (
            "Integrity record only; it does not validate assay units, mechanistic "
            "assumptions, treatment effects, or clinical utility."
        ),
    }
    audit = Path(f"{output}.twin-audit.json")
    audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit


def verify_augmentation_audit(
    audit_path: str | Path,
    output_path: str | Path,
    input_paths: dict[str, str | Path],
) -> dict:
    audit = Path(audit_path)
    try:
        payload = json.loads(audit.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid augmentation audit: {audit}") from error
    if (not isinstance(payload, dict) or payload.get("schema_version") != 1
            or payload.get("stage") != "label_free_twin_augmentation"
            or payload.get("labels_used") is not False):
        raise ValueError("unsupported or non-label-free augmentation audit")
    output = Path(output_path)
    recorded_output = payload.get("output", {})
    if (recorded_output.get("file") != output.name
            or not output.is_file()
            or recorded_output.get("sha256") != sha256_file(output)):
        raise ValueError("augmented feature output SHA-256 mismatch")
    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("augmentation audit has no input records")
    verified_inputs = []
    for name, record in inputs.items():
        if record is None:
            continue
        if (not isinstance(record, dict) or not isinstance(record.get("file"), str)
                or not isinstance(record.get("sha256"), str)):
            raise ValueError(f"invalid augmentation input record: {name}")
        if name not in input_paths:
            raise ValueError(f"augmentation verification requires input path: {name}")
        target = Path(input_paths[name])
        if (not target.is_file() or target.name != record["file"]
                or sha256_file(target) != record["sha256"]):
            raise ValueError(f"augmentation input SHA-256 mismatch: {name}")
        derivation = record.get("derivation_audit")
        if derivation is not None:
            suffixes = {
                "features": ".feature-audit.json",
                "spatial_initial_zones": ".spatial-audit.json",
            }
            suffix = suffixes.get(name)
            derived_audit = Path(f"{target}{suffix}") if suffix else None
            if (derived_audit is None or not isinstance(derivation, dict)
                    or derivation.get("file") != derived_audit.name
                    or not derived_audit.is_file()
                    or derivation.get("sha256") != sha256_file(derived_audit)):
                raise ValueError(f"{name} derivation audit SHA-256 mismatch")
        verified_inputs.append(name)
    settings = payload.get("settings")
    if (not isinstance(settings, dict) or settings.get("labels_used") is not False
            or settings.get("cross_layer_sampling")
            != "independent_systemic_spatial_mapping_no_joint_covariance"):
        raise ValueError("invalid augmentation replay settings")
    try:
        days, dose = float(settings["days"]), float(settings["dose"])
        dose_days = tuple(float(value) for value in settings["dose_days"])
        samples, seed = int(settings["uncertainty_samples"]), int(settings["uncertainty_seed"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid numeric augmentation replay settings") from error
    if (not np.isfinite(days) or days <= 0 or not np.isfinite(dose) or dose < 0
            or samples < 0 or samples == 1):
        raise ValueError("augmentation replay settings are outside supported domains")

    def supplied(name: str):
        return input_paths.get(name) if inputs.get(name) is not None else None

    systemic = _load_exact_parameters(supplied("systemic_parameters"), TwinParameters)
    spatial = _load_exact_parameters(supplied("spatial_parameters"), SpatialParameters)
    mapping = _load_exact_parameters(supplied("feature_state_mapping"), FeatureStateMapping)
    systemic_draws = spatial_draws = mapping_draws = None
    for diagnostics_name, parameters, offset, parameter_class, target_name in (
        ("systemic_diagnostics", systemic, 0, TwinParameters, "systemic"),
        ("spatial_diagnostics", spatial, 1, SpatialParameters, "spatial"),
        ("mapping_diagnostics", mapping, 2, FeatureStateMapping, "mapping"),
    ):
        diagnostics_path = supplied(diagnostics_name)
        if diagnostics_path is None:
            continue
        if parameters is None or samples < 2:
            raise ValueError(f"augmentation replay cannot sample {target_name} diagnostics")
        diagnostics = json.loads(Path(diagnostics_path).read_text())
        _, draws = calibration_parameter_draws(parameters, diagnostics, samples, seed + offset)
        converted = tuple(replace(parameters, **draw) for draw in draws)
        if parameter_class is TwinParameters:
            systemic_draws = converted
        elif parameter_class is SpatialParameters:
            spatial_draws = converted
        else:
            mapping_draws = converted
    replayed = augment_with_twin(
        pd.read_csv(input_paths["features"]), days, dose_days, dose,
        pd.read_csv(supplied("spatial_initial_zones"))
        if supplied("spatial_initial_zones") is not None else None,
        systemic_parameters=systemic, spatial_parameters=spatial,
        systemic_parameter_draws=systemic_draws,
        spatial_parameter_draws=spatial_draws,
        feature_state_mapping=mapping,
        feature_state_mapping_draws=mapping_draws,
    )
    replayed_sha = hashlib.sha256(
        replayed.to_csv(index=False).encode("utf-8")
    ).hexdigest()
    if replayed_sha != recorded_output["sha256"]:
        raise ValueError("augmentation derivation replay mismatch")
    return {
        "valid": True,
        "verified_inputs": sorted(verified_inputs),
        "output_sha256": recorded_output["sha256"],
        "audit_sha256": sha256_file(audit),
        "derivation_replayed": True,
        "payload": payload,
    }
