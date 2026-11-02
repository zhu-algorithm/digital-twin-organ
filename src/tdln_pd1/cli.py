from __future__ import annotations

import argparse
from dataclasses import asdict, fields, replace
import hashlib
import json
from pathlib import Path
import numpy as np
from .data import load_cells, load_labels
from .features import aggregate_cells
from .model import evaluate_external, leave_one_site_out, predict_locked, train_evaluate
from .simulator import Dose, TwinParameters, simulate
from .assimilation import (calibrate, calibration_ensemble, calibration_parameter_draws,
                           save_calibration)
from .hybrid import FeatureStateMapping, augment_with_twin
from .mapping_assimilation import (
    calibrate_feature_state_mapping, save_mapping_calibration,
)
from .spatial_dynamics import (SpatialParameters, simulate_spatial, simulate_spatial_coupled,
                               simulate_spatial_counterfactual)
from .spatial_data import aggregate_spatial_spots
from .spatial_assimilation import (
    calibrate_spatial, save_spatial_calibration, spatial_calibration_ensemble,
)
from .sensitivity import global_sensitivity, global_spatial_sensitivity
from .sensitivity_audit import verify_sensitivity_audit, write_sensitivity_audit
from .provenance import attach_derivation_audit, verify_manifest
from .qc import cohort_qc
from .verification import verify_twin
from .validation_audit import (verify_prediction_audit, write_evaluation_audit,
                               write_prediction_audit)
from .calibration_audit import verify_calibration_audit, write_calibration_audit
from .augmentation_audit import verify_augmentation_audit, write_augmentation_audit
from .feature_audit import verify_feature_audit, write_feature_audit
from .spatial_feature_audit import (
    verify_spatial_feature_audit, write_spatial_feature_audit,
)
from .uncertainty_audit import verify_uncertainty_audit, write_uncertainty_audit


def _load_parameter_file(path: str | None, parameter_class):
    if not path:
        return None, None
    source = Path(path)
    payload = json.loads(source.read_text())
    if not isinstance(payload, dict):
        raise ValueError("parameter file must contain one JSON object")
    expected = {field.name for field in fields(parameter_class)}
    if set(payload) != expected:
        raise ValueError(
            f"parameter schema mismatch for {parameter_class.__name__}: "
            f"missing={sorted(expected - set(payload))}, extra={sorted(set(payload) - expected)}"
        )
    converted = {name: float(payload[name]) for name in expected}
    values = np.asarray([converted[name] for name in sorted(expected)], dtype=float)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("parameter files require finite, strictly positive values")
    fingerprint = hashlib.sha256(source.read_bytes()).hexdigest()
    return parameter_class(**converted), {
        "file": source.name, "sha256": fingerprint,
        "parameter_class": parameter_class.__name__,
    }


def _write_parameter_provenance(output: str, sources: dict) -> None:
    used = {name: value for name, value in sources.items() if value is not None}
    if used:
        Path(f"{output}.parameters.json").write_text(json.dumps({
            "schema_version": 1, "sources": used, "research_only": True,
        }, indent=2) + "\n")


def _load_patient_spatial_initial(path: str | None, patient_id: str | None):
    """Load one patient's complete zone state and bind it to its source file."""
    if not path:
        if patient_id is not None:
            raise ValueError("--patient-id requires --initial-zones")
        return None, None
    import pandas as pd
    source = Path(path)
    states = pd.read_csv(source)
    if "patient_id" not in states:
        raise ValueError("initial zone states need patient_id")
    patients = states.patient_id.astype(str).unique().tolist()
    selected_id = str(patient_id) if patient_id is not None else (
        patients[0] if len(patients) == 1 else None
    )
    if selected_id is None:
        raise ValueError("--patient-id is required when initial states contain multiple patients")
    selected = states[states.patient_id.astype(str).eq(selected_id)].copy()
    if selected.empty:
        raise ValueError(f"patient not found in initial zone states: {selected_id}")
    provenance = {
        "file": source.name,
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "patient_id": selected_id,
        "rows": int(len(selected)),
    }
    provenance = attach_derivation_audit(
        provenance, source, ".spatial-audit.json",
        "label_free_spatial_state_featurization",
    )
    return selected, provenance


def _source_record(path: str) -> dict:
    source = Path(path)
    return {"file": source.name,
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}


def _training_input_provenance(feature_path: str) -> dict:
    source = Path(feature_path)
    provenance = {
        "feature_file": source.name,
        "feature_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    sidecar = Path(f"{feature_path}.parameters.json")
    if sidecar.exists():
        payload = json.loads(sidecar.read_text())
        if payload.get("schema_version") != 1 or not isinstance(payload.get("sources"), dict):
            raise ValueError("invalid feature parameter-provenance sidecar")
        expected_classes = {
            "systemic": "TwinParameters",
            "spatial": "SpatialParameters",
            "feature_state_mapping": "FeatureStateMapping",
        }
        for name, entry in payload["sources"].items():
            digest = entry.get("sha256") if isinstance(entry, dict) else None
            parameter_class = entry.get("parameter_class") if isinstance(entry, dict) else None
            if (name not in expected_classes or parameter_class != expected_classes.get(name)
                    or not isinstance(digest, str)
                    or len(digest) != 64
                    or any(character not in "0123456789abcdef" for character in digest.lower())):
                raise ValueError("invalid parameter source fingerprint in feature sidecar")
        provenance["parameter_provenance"] = payload
        provenance["parameter_sidecar_sha256"] = hashlib.sha256(sidecar.read_bytes()).hexdigest()
    else:
        provenance["parameter_provenance"] = None
    augmentation_audit = Path(f"{feature_path}.twin-audit.json")
    if augmentation_audit.exists():
        payload = json.loads(augmentation_audit.read_text())
        recorded_output = payload.get("output", {}) if isinstance(payload, dict) else {}
        if (payload.get("schema_version") != 1
                or payload.get("stage") != "label_free_twin_augmentation"
                or payload.get("labels_used") is not False
                or recorded_output.get("file") != source.name
                or recorded_output.get("sha256") != provenance["feature_sha256"]):
            raise ValueError("invalid or mismatched feature augmentation audit")
        provenance["augmentation_audit"] = payload
        provenance["augmentation_audit_sha256"] = hashlib.sha256(
            augmentation_audit.read_bytes()
        ).hexdigest()
    else:
        provenance["augmentation_audit"] = None
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description="TDLN PD-1 response baseline")
    sub = parser.add_subparsers(dest="command", required=True)
    featurize = sub.add_parser("featurize")
    featurize.add_argument("--cells", required=True)
    featurize.add_argument("--output", required=True)
    featurize.add_argument("--time-column")
    featurize.add_argument("--baseline-values", default="")
    featurize.add_argument("--feature-tissue-pattern", default="tdln|lymph")
    featurize.add_argument("--min-signature-gene-coverage", type=float, default=0.5)
    train = sub.add_parser("train")
    train.add_argument("--features", required=True)
    train.add_argument("--labels", required=True)
    train.add_argument("--output", required=True)
    train.add_argument("--folds", type=int, default=5)
    train.add_argument("--ensemble-members", type=int, default=50)
    predict = sub.add_parser("predict")
    predict.add_argument("--model", required=True)
    predict.add_argument("--features", required=True)
    predict.add_argument("--output", required=True)
    evaluate = sub.add_parser("evaluate")
    evaluate.add_argument("--predictions", required=True)
    evaluate.add_argument("--labels", required=True)
    evaluate.add_argument("--output", required=True)
    evaluate.add_argument("--bootstrap", type=int, default=500)
    evaluate.add_argument("--prediction-audit",
                          help="defaults to <predictions>.audit.json")
    site_cv = sub.add_parser("site-cv")
    site_cv.add_argument("--features", required=True)
    site_cv.add_argument("--labels", required=True)
    site_cv.add_argument("--site-column", default="center")
    site_cv.add_argument("--bootstrap", type=int, default=500)
    site_cv.add_argument("--output", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--days", type=float, default=42.0)
    sim.add_argument("--dose-days", default="0,14,28")
    sim.add_argument("--dose", type=float, default=1.0)
    sim.add_argument("--parameters", help="fitted TwinParameters JSON")
    sim.add_argument("--output", required=True)
    cal = sub.add_parser("calibrate")
    cal.add_argument("--observations", required=True)
    cal.add_argument("--fit", required=True, help="comma-separated TwinParameters fields")
    cal.add_argument("--dose-days", default="")
    cal.add_argument("--dose", type=float, default=1.0)
    cal.add_argument("--output", required=True)
    cal.add_argument("--profile-points", type=int, default=0,
                     help="zero disables; otherwise use an odd integer >=3")
    cal.add_argument("--profile-fold-range", type=float, default=4.0)
    mapping_cal = sub.add_parser("calibrate-feature-state-mapping")
    mapping_cal.add_argument("--observations", required=True)
    mapping_cal.add_argument("--fit", required=True,
                             help="comma-separated FeatureStateMapping fields")
    mapping_cal.add_argument("--base-mapping", help="starting FeatureStateMapping JSON")
    mapping_cal.add_argument("--max-nfev", type=int, default=200)
    mapping_cal.add_argument("--output", required=True)
    spatial_cal = sub.add_parser("calibrate-spatial")
    spatial_cal.add_argument("--observations", required=True)
    spatial_cal.add_argument("--fit", required=True,
                             help="comma-separated SpatialParameters fields")
    spatial_cal.add_argument("--initial-zones")
    spatial_cal.add_argument("--patient-id")
    spatial_cal.add_argument("--dose-days", default="")
    spatial_cal.add_argument("--dose", type=float, default=1.0)
    spatial_cal.add_argument("--max-nfev", type=int, default=200)
    spatial_cal.add_argument("--profile-points", type=int, default=0,
                             help="zero disables; otherwise use an odd integer >=3")
    spatial_cal.add_argument("--profile-fold-range", type=float, default=4.0)
    spatial_cal.add_argument("--output", required=True)
    uncertainty = sub.add_parser("propagate-uncertainty")
    uncertainty.add_argument("--parameters", required=True)
    uncertainty.add_argument("--diagnostics", required=True)
    uncertainty.add_argument("--days", type=float, default=42.0)
    uncertainty.add_argument("--dose-days", default="0,14,28")
    uncertainty.add_argument("--dose", type=float, default=1.0)
    uncertainty.add_argument("--samples", type=int, default=200)
    uncertainty.add_argument("--seed", type=int, default=42)
    uncertainty.add_argument("--output", required=True)
    spatial_uncertainty = sub.add_parser("propagate-spatial-uncertainty")
    spatial_uncertainty.add_argument("--parameters", required=True)
    spatial_uncertainty.add_argument("--diagnostics", required=True)
    spatial_uncertainty.add_argument("--systemic-parameters")
    spatial_uncertainty.add_argument("--initial-zones")
    spatial_uncertainty.add_argument("--patient-id")
    spatial_uncertainty.add_argument("--days", type=float, default=14.0)
    spatial_uncertainty.add_argument("--dose-days", default="0")
    spatial_uncertainty.add_argument("--dose", type=float, default=1.0)
    spatial_uncertainty.add_argument("--samples", type=int, default=200)
    spatial_uncertainty.add_argument("--seed", type=int, default=42)
    spatial_uncertainty.add_argument("--output", required=True)
    aug = sub.add_parser("augment-twin")
    aug.add_argument("--features", required=True)
    aug.add_argument("--days", type=float, default=42.0)
    aug.add_argument("--dose-days", default="0,14,28")
    aug.add_argument("--dose", type=float, default=1.0)
    aug.add_argument("--spatial-initial-zones",
                     help="patient-zone-species states from spatial-featurize")
    aug.add_argument("--systemic-parameters", help="fitted TwinParameters JSON")
    aug.add_argument("--spatial-parameters", help="fitted SpatialParameters JSON")
    aug.add_argument("--feature-state-mapping",
                     help="versioned FeatureStateMapping JSON")
    aug.add_argument("--systemic-diagnostics", help="systemic calibration diagnostics JSON")
    aug.add_argument("--spatial-diagnostics", help="spatial calibration diagnostics JSON")
    aug.add_argument("--uncertainty-samples", type=int, default=0,
                     help="0 disables; otherwise requires at least 2 correlated draws")
    aug.add_argument("--uncertainty-seed", type=int, default=42)
    aug.add_argument("--output", required=True)
    spatial = sub.add_parser("simulate-spatial")
    spatial.add_argument("--days", type=float, default=14.0)
    spatial.add_argument("--step", type=float, default=0.1)
    spatial.add_argument("--transport-only", action="store_true")
    spatial.add_argument("--dose-days", default="", help="comma-separated anti-PD-1 dose days")
    spatial.add_argument("--dose", type=float, default=1.0)
    spatial.add_argument("--paired-counterfactual", action="store_true")
    spatial.add_argument("--systemic-parameters", help="fitted TwinParameters JSON")
    spatial.add_argument("--spatial-parameters", help="fitted SpatialParameters JSON")
    spatial.add_argument("--initial-zones", help="long patient-zone-species state CSV")
    spatial.add_argument("--patient-id", help="patient to select from --initial-zones")
    spatial.add_argument("--output", required=True)
    spatial_features = sub.add_parser("spatial-featurize")
    spatial_features.add_argument("--spots", required=True)
    spatial_features.add_argument("--weight-column")
    spatial_features.add_argument("--output", required=True)
    sens = sub.add_parser("sensitivity")
    sens.add_argument("--parameters", default="", help="comma-separated fields; empty uses defaults")
    sens.add_argument("--base-parameters", help="fitted TwinParameters JSON")
    sens.add_argument("--days", type=float, default=42.0)
    sens.add_argument("--dose-days", default=None,
                      help="comma-separated schedule; omitted uses adaptive days 0/14/28")
    sens.add_argument("--dose", type=float, default=1.0)
    sens.add_argument("--samples", type=int, default=128)
    sens.add_argument("--fold-range", type=float, default=4.0)
    sens.add_argument("--seed", type=int, default=42)
    sens.add_argument("--output", required=True)
    spatial_sens = sub.add_parser("spatial-sensitivity")
    spatial_sens.add_argument("--parameters", default="", help="comma-separated fields; empty uses defaults")
    spatial_sens.add_argument("--spatial-parameters", help="fitted SpatialParameters JSON")
    spatial_sens.add_argument("--systemic-parameters", help="fitted TwinParameters JSON")
    spatial_sens.add_argument("--initial-zones", help="long patient-zone-species state CSV")
    spatial_sens.add_argument("--patient-id", help="patient to select from --initial-zones")
    spatial_sens.add_argument("--days", type=float, default=14.0)
    spatial_sens.add_argument("--checkpoint-occupancy", type=float, default=0.8)
    spatial_sens.add_argument("--dose-days", default="",
                              help="when supplied, use systemic time-varying occupancy")
    spatial_sens.add_argument("--dose", type=float, default=1.0)
    spatial_sens.add_argument("--samples", type=int, default=128)
    spatial_sens.add_argument("--fold-range", type=float, default=4.0)
    spatial_sens.add_argument("--seed", type=int, default=42)
    spatial_sens.add_argument("--output", required=True)
    qc = sub.add_parser("qc")
    qc.add_argument("--cells", required=True)
    qc.add_argument("--labels", required=True)
    qc.add_argument("--manifest")
    qc.add_argument("--output", required=True)
    qc.add_argument("--time-column")
    qc.add_argument("--baseline-values", default="")
    qc.add_argument("--feature-tissue-pattern", default="tdln|lymph")
    qc.add_argument("--min-signature-gene-coverage", type=float, default=0.5)
    verify = sub.add_parser("verify")
    verify.add_argument("--output", required=True)
    verify_sensitivity = sub.add_parser("verify-sensitivity")
    verify_sensitivity.add_argument("--audit", required=True)
    verify_sensitivity.add_argument("--base-parameters")
    verify_sensitivity.add_argument("--spatial-parameters")
    verify_sensitivity.add_argument("--systemic-parameters")
    verify_sensitivity.add_argument("--initial-zones")
    verify_calibration = sub.add_parser("verify-calibration")
    verify_calibration.add_argument("--audit", required=True)
    verify_calibration.add_argument("--observations", required=True)
    verify_calibration.add_argument("--initial-zones")
    verify_calibration.add_argument("--base-mapping")
    verify_augmentation = sub.add_parser("verify-augmentation")
    verify_augmentation.add_argument("--audit", required=True)
    verify_augmentation.add_argument("--output", required=True)
    verify_augmentation.add_argument("--features", required=True)
    verify_augmentation.add_argument("--spatial-initial-zones")
    verify_augmentation.add_argument("--systemic-parameters")
    verify_augmentation.add_argument("--spatial-parameters")
    verify_augmentation.add_argument("--feature-state-mapping")
    verify_augmentation.add_argument("--systemic-diagnostics")
    verify_augmentation.add_argument("--spatial-diagnostics")
    verify_features = sub.add_parser("verify-features")
    verify_features.add_argument("--audit", required=True)
    verify_features.add_argument("--cells", required=True)
    verify_features.add_argument("--output", required=True)
    verify_spatial_features = sub.add_parser("verify-spatial-features")
    verify_spatial_features.add_argument("--audit", required=True)
    verify_spatial_features.add_argument("--spots", required=True)
    verify_spatial_features.add_argument("--output", required=True)
    verify_uncertainty = sub.add_parser("verify-uncertainty")
    verify_uncertainty.add_argument("--audit", required=True)
    verify_uncertainty.add_argument("--parameters", required=True)
    verify_uncertainty.add_argument("--diagnostics", required=True)
    verify_spatial_uncertainty = sub.add_parser("verify-spatial-uncertainty")
    verify_spatial_uncertainty.add_argument("--audit", required=True)
    verify_spatial_uncertainty.add_argument("--parameters", required=True)
    verify_spatial_uncertainty.add_argument("--diagnostics", required=True)
    verify_spatial_uncertainty.add_argument("--systemic-parameters")
    verify_spatial_uncertainty.add_argument("--initial-zones")
    args = parser.parse_args()
    if args.command == "featurize":
        baseline_values = tuple(value.strip() for value in args.baseline_values.split(",") if value.strip())
        aggregate_cells(load_cells(args.cells, (args.time_column,) if args.time_column else ()),
                        args.time_column, baseline_values or None,
                        args.feature_tissue_pattern,
                        args.min_signature_gene_coverage).to_csv(args.output, index=False)
        write_feature_audit(args.cells, args.output, {
            "time_column": args.time_column,
            "baseline_values": list(baseline_values),
            "feature_tissue_pattern": args.feature_tissue_pattern,
            "minimum_signature_gene_coverage": args.min_signature_gene_coverage,
            "labels_used": False,
        })
    elif args.command == "train":
        import pandas as pd
        metrics = train_evaluate(
            pd.read_csv(args.features), load_labels(args.labels), Path(args.output),
            args.folds, args.ensemble_members, _training_input_provenance(args.features),
        )
        print(json.dumps(metrics, indent=2))
    elif args.command == "predict":
        import pandas as pd
        predictions = predict_locked(args.model, pd.read_csv(args.features))
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        predictions.to_csv(args.output, index=False)
        audit_path = write_prediction_audit(
            args.model, args.features, args.output, len(predictions),
        )
        print(json.dumps({"n_patients": len(predictions), "output": args.output,
                          "model_refit": False, "audit": str(audit_path)}, indent=2))
    elif args.command == "evaluate":
        import pandas as pd
        prediction_audit = verify_prediction_audit(
            args.predictions, args.prediction_audit,
        )
        metrics = evaluate_external(pd.read_csv(args.predictions), load_labels(args.labels), args.bootstrap)
        metrics["validation_chain"] = {
            "prediction_audit_sha256": prediction_audit["audit_sha256"],
            "model_sha256": prediction_audit["payload"]["inputs"]["model"]["sha256"],
            "features_sha256": prediction_audit["payload"]["inputs"]["features"]["sha256"],
            "integrity_verified": True,
        }
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(metrics, indent=2) + "\n")
        write_evaluation_audit(
            args.predictions, prediction_audit, args.labels, args.output, args.bootstrap,
        )
        print(json.dumps(metrics, indent=2))
    elif args.command == "site-cv":
        import pandas as pd
        metrics, predictions = leave_one_site_out(pd.read_csv(args.features), pd.read_csv(args.labels),
                                                  args.site_column, args.bootstrap)
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        predictions.to_csv(output / "site_cv_predictions.csv", index=False)
        (output / "site_cv_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
        print(json.dumps(metrics, indent=2))
    elif args.command == "simulate":
        doses = tuple(Dose(float(day), args.dose)
                      for day in args.dose_days.split(",") if day.strip())
        parameters, provenance = _load_parameter_file(args.parameters, TwinParameters)
        trajectory = simulate(days=args.days, doses=doses, parameters=parameters)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        trajectory.to_csv(args.output, index=False)
        _write_parameter_provenance(args.output, {"systemic": provenance})
        print(json.dumps({
            "final_tumour_change_pct": float(trajectory.tumour_change_pct.iloc[-1]),
            "peak_pd1_occupancy": float(trajectory.pd1_occupancy.max()),
            "output": args.output,
            "research_only": True,
        }, indent=2))
    elif args.command == "calibrate":
        import pandas as pd
        doses = tuple(Dose(float(day), args.dose) for day in args.dose_days.split(",") if day.strip())
        fitted, diagnostics, compared = calibrate(
            pd.read_csv(args.observations), doses=doses,
            fit_parameters=tuple(name.strip() for name in args.fit.split(",") if name.strip()),
            profile_points=args.profile_points,
            profile_fold_range=args.profile_fold_range,
        )
        save_calibration(args.output, fitted, diagnostics, compared)
        fit_names = tuple(name.strip() for name in args.fit.split(",") if name.strip())
        write_calibration_audit(
            args.output, "calibration_provenance.json", "systemic_calibration",
            inputs={"observations": _source_record(args.observations)},
            settings={
                "fit_parameters": list(fit_names),
                "dose_days": [dose.day for dose in doses],
                "dose_amounts": [dose.amount for dose in doses],
                "profile_points": args.profile_points,
                "profile_fold_range": args.profile_fold_range,
            },
            generated_filenames=(
                "fitted_parameters.json", "calibration_diagnostics.json",
                "calibration_residuals.csv",
            ),
        )
        print(json.dumps(diagnostics, indent=2))
    elif args.command == "calibrate-feature-state-mapping":
        import pandas as pd
        base, base_provenance = _load_parameter_file(
            args.base_mapping, FeatureStateMapping
        )
        fit_names = tuple(name.strip() for name in args.fit.split(",") if name.strip())
        fitted, diagnostics, compared = calibrate_feature_state_mapping(
            pd.read_csv(args.observations), base=base,
            fit_parameters=fit_names, max_nfev=args.max_nfev,
        )
        save_mapping_calibration(args.output, fitted, diagnostics, compared)
        write_calibration_audit(
            args.output, "mapping_calibration_provenance.json",
            "feature_state_mapping_calibration",
            inputs={
                "observations": _source_record(args.observations),
                "base_mapping": base_provenance,
            },
            settings={"fit_parameters": list(fit_names), "max_nfev": args.max_nfev},
            generated_filenames=(
                "fitted_feature_state_mapping.json",
                "mapping_calibration_diagnostics.json",
                "mapping_calibration_residuals.csv",
            ),
        )
        print(json.dumps(diagnostics, indent=2))
    elif args.command == "calibrate-spatial":
        import pandas as pd
        observations = pd.read_csv(args.observations)
        initial, initial_provenance = _load_patient_spatial_initial(
            args.initial_zones, args.patient_id
        )
        doses = tuple(Dose(float(day), args.dose)
                      for day in args.dose_days.split(",") if day.strip())
        systemic = None
        if doses:
            calibration_days = max(float(pd.to_numeric(observations.day).max()), 1e-6)
            systemic = simulate(days=calibration_days, doses=doses)
        fitted, diagnostics, compared = calibrate_spatial(
            observations, initial=initial, systemic_trajectory=systemic,
            fit_parameters=tuple(name.strip() for name in args.fit.split(",") if name.strip()),
            max_nfev=args.max_nfev,
            profile_points=args.profile_points,
            profile_fold_range=args.profile_fold_range,
        )
        save_spatial_calibration(args.output, fitted, diagnostics, compared)
        fit_names = tuple(name.strip() for name in args.fit.split(",") if name.strip())
        write_calibration_audit(
            args.output, "spatial_calibration_provenance.json", "spatial_calibration",
            inputs={"observations": _source_record(args.observations),
                    "initial_state_source": initial_provenance},
            settings={
                "fit_parameters": list(fit_names),
                "dose_days": [dose.day for dose in doses],
                "dose_amounts": [dose.amount for dose in doses],
                "systemic_parameters": asdict(TwinParameters()) if doses else None,
                "max_nfev": args.max_nfev,
                "profile_points": args.profile_points,
                "profile_fold_range": args.profile_fold_range,
            },
            generated_filenames=(
                "fitted_spatial_parameters.json",
                "spatial_calibration_diagnostics.json",
                "spatial_calibration_residuals.csv",
            ),
        )
        print(json.dumps(diagnostics, indent=2))
    elif args.command == "propagate-uncertainty":
        parameters, parameter_source = _load_parameter_file(
            args.parameters, TwinParameters,
        )
        diagnostics_path = Path(args.diagnostics)
        diagnostics = json.loads(diagnostics_path.read_text())
        doses = tuple(Dose(float(day), args.dose)
                      for day in args.dose_days.split(",") if day.strip())
        samples = calibration_ensemble(
            parameters, diagnostics, doses, args.days, args.samples, args.seed,
        )
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        samples.to_csv(output / "uncertainty_samples.csv", index=False)
        endpoints = (
            "tumour_cells", "tumour_change_pct", "effector_tumour",
            "tpex_tdln", "peak_pd1_occupancy",
        )
        summary = {
            "method": "local_multivariate_log_parameter_normal",
            "samples": int(len(samples)), "seed": args.seed,
            "endpoints": {
                name: {
                    "median": float(samples[name].median()),
                    "lower_95": float(samples[name].quantile(0.025)),
                    "upper_95": float(samples[name].quantile(0.975)),
                }
                for name in endpoints
            },
            "research_only": True,
            "limitation": "Local parameter uncertainty only; not a clinical interval.",
        }
        (output / "uncertainty_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        diagnostics_source = _source_record(args.diagnostics)
        write_uncertainty_audit(output, {
            "parameters": parameter_source,
            "diagnostics": diagnostics_source,
        }, {
            "days": args.days, "dose_days": [dose.day for dose in doses],
            "dose": args.dose, "samples": args.samples, "seed": args.seed,
            "sampling": "multivariate_normal_in_fitted_log_parameter_space",
        })
        print(json.dumps(summary, indent=2))
    elif args.command == "propagate-spatial-uncertainty":
        parameters, parameter_source = _load_parameter_file(
            args.parameters, SpatialParameters,
        )
        systemic_parameters, systemic_source = _load_parameter_file(
            args.systemic_parameters, TwinParameters,
        )
        initial, initial_source = _load_patient_spatial_initial(
            args.initial_zones, args.patient_id,
        )
        diagnostics_path = Path(args.diagnostics)
        diagnostics = json.loads(diagnostics_path.read_text())
        doses = tuple(Dose(float(day), args.dose)
                      for day in args.dose_days.split(",") if day.strip())
        if not doses or not np.isfinite(args.dose) or args.dose <= 0:
            raise ValueError(
                "spatial uncertainty treated arm requires an explicit positive dose schedule"
            )
        systemic_trajectory = simulate(
            days=args.days, doses=doses, parameters=systemic_parameters,
        )
        samples = spatial_calibration_ensemble(
            parameters, diagnostics, args.days, initial, systemic_trajectory,
            args.samples, args.seed,
        )
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        samples.to_csv(output / "spatial_uncertainty_samples.csv", index=False)
        endpoints = (
            "control_efferent_effector", "treated_efferent_effector",
            "pd1_efferent_effector_delta", "control_paracortex_tpex",
            "treated_paracortex_tpex", "pd1_paracortex_tpex_delta",
        )
        summary = {
            "method": "paired_spatial_local_multivariate_log_parameter_normal",
            "samples": int(len(samples)), "seed": args.seed,
            "checkpoint_mode": "dose_driven_systemic_trajectory",
            "peak_pd1_occupancy": float(systemic_trajectory.pd1_occupancy.max()),
            "endpoints": {
                name: {
                    "median": float(samples[name].median()),
                    "lower_95": float(samples[name].quantile(0.025)),
                    "upper_95": float(samples[name].quantile(0.975)),
                }
                for name in endpoints
            },
            "research_only": True,
            "limitation": "Local spatial parameter uncertainty only; not a clinical interval.",
        }
        (output / "spatial_uncertainty_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        write_uncertainty_audit(
            output,
            inputs={
                "parameters": parameter_source,
                "diagnostics": _source_record(args.diagnostics),
                "systemic_parameters": systemic_source,
                "initial_zones": initial_source,
            },
            settings={
                "days": args.days, "dose_days": [dose.day for dose in doses],
                "dose": args.dose, "samples": args.samples, "seed": args.seed,
                "sampling": "multivariate_normal_in_fitted_log_parameter_space",
                "paired_control": True,
            },
            stage="spatial_calibration_uncertainty_propagation",
            output_filenames=(
                "spatial_uncertainty_samples.csv",
                "spatial_uncertainty_summary.json",
            ),
        )
        print(json.dumps(summary, indent=2))
    elif args.command == "augment-twin":
        import pandas as pd
        dose_days = tuple(float(day) for day in args.dose_days.split(",") if day.strip())
        spatial_states = (pd.read_csv(args.spatial_initial_zones)
                          if args.spatial_initial_zones else None)
        systemic_parameters, systemic_provenance = _load_parameter_file(
            args.systemic_parameters, TwinParameters
        )
        spatial_parameters, spatial_provenance = _load_parameter_file(
            args.spatial_parameters, SpatialParameters
        )
        feature_state_mapping, mapping_provenance = _load_parameter_file(
            args.feature_state_mapping, FeatureStateMapping
        )
        if args.uncertainty_samples < 0 or args.uncertainty_samples == 1:
            raise ValueError("--uncertainty-samples must be 0 or at least 2")
        if (args.systemic_diagnostics or args.spatial_diagnostics) and not args.uncertainty_samples:
            raise ValueError("diagnostics require --uncertainty-samples of at least 2")
        if args.systemic_diagnostics and systemic_parameters is None:
            raise ValueError("--systemic-diagnostics requires --systemic-parameters")
        if args.spatial_diagnostics and spatial_parameters is None:
            raise ValueError("--spatial-diagnostics requires --spatial-parameters")
        if args.uncertainty_samples and not (
            args.systemic_diagnostics or args.spatial_diagnostics
        ):
            raise ValueError("uncertainty sampling requires at least one diagnostics file")
        systemic_draws = spatial_draws = None
        if args.systemic_diagnostics:
            diagnostics = json.loads(Path(args.systemic_diagnostics).read_text())
            _, draws = calibration_parameter_draws(
                systemic_parameters, diagnostics, args.uncertainty_samples,
                args.uncertainty_seed,
            )
            systemic_draws = tuple(replace(systemic_parameters, **draw) for draw in draws)
        if args.spatial_diagnostics:
            diagnostics = json.loads(Path(args.spatial_diagnostics).read_text())
            _, draws = calibration_parameter_draws(
                spatial_parameters, diagnostics, args.uncertainty_samples,
                args.uncertainty_seed + 1,
            )
            spatial_draws = tuple(replace(spatial_parameters, **draw) for draw in draws)
        augmented = augment_with_twin(pd.read_csv(args.features), args.days, dose_days,
                                      args.dose, spatial_states,
                                      systemic_parameters=systemic_parameters,
                                      spatial_parameters=spatial_parameters,
                                      systemic_parameter_draws=systemic_draws,
                                      spatial_parameter_draws=spatial_draws,
                                      feature_state_mapping=feature_state_mapping)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        augmented.to_csv(args.output, index=False)
        _write_parameter_provenance(args.output, {
            "systemic": systemic_provenance, "spatial": spatial_provenance,
            "feature_state_mapping": mapping_provenance,
        })
        feature_record = _source_record(args.features)
        feature_audit = Path(f"{args.features}.feature-audit.json")
        if feature_audit.exists():
            feature_payload = json.loads(feature_audit.read_text())
            if (feature_payload.get("stage") != "label_free_single_cell_featurization"
                    or feature_payload.get("labels_used") is not False
                    or feature_payload.get("output", {}).get("file") != Path(args.features).name
                    or feature_payload.get("output", {}).get("sha256") != feature_record["sha256"]):
                raise ValueError("invalid or mismatched feature derivation audit")
            feature_record["derivation_audit"] = _source_record(str(feature_audit))
        spatial_record = None
        if args.spatial_initial_zones:
            spatial_record = _source_record(args.spatial_initial_zones)
            spatial_audit = Path(f"{args.spatial_initial_zones}.spatial-audit.json")
            if spatial_audit.exists():
                spatial_payload = json.loads(spatial_audit.read_text())
                if (spatial_payload.get("stage") != "label_free_spatial_state_featurization"
                        or spatial_payload.get("labels_used") is not False
                        or spatial_payload.get("output", {}).get("file")
                        != Path(args.spatial_initial_zones).name
                        or spatial_payload.get("output", {}).get("sha256")
                        != spatial_record["sha256"]):
                    raise ValueError("invalid or mismatched spatial derivation audit")
                spatial_record["derivation_audit"] = _source_record(str(spatial_audit))
        write_augmentation_audit(
            args.output,
            inputs={
                "features": feature_record,
                "spatial_initial_zones": spatial_record,
                "systemic_parameters": systemic_provenance,
                "spatial_parameters": spatial_provenance,
                "feature_state_mapping": mapping_provenance,
                "systemic_diagnostics": (_source_record(args.systemic_diagnostics)
                                           if args.systemic_diagnostics else None),
                "spatial_diagnostics": (_source_record(args.spatial_diagnostics)
                                          if args.spatial_diagnostics else None),
            },
            settings={
                "days": args.days, "dose_days": list(dose_days),
                "dose": args.dose, "labels_used": False,
                "uncertainty_samples": args.uncertainty_samples,
                "uncertainty_seed": args.uncertainty_seed,
                "cross_layer_sampling": "independent_no_joint_covariance",
            },
        )
        print(json.dumps({"n_patients": len(augmented), "n_features": len(augmented.columns) - 1,
                          "output": args.output, "labels_used": False}, indent=2))
    elif args.command == "simulate-spatial":
        systemic_parameters, systemic_provenance = _load_parameter_file(
            args.systemic_parameters, TwinParameters
        )
        spatial_parameters, spatial_provenance = _load_parameter_file(
            args.spatial_parameters, SpatialParameters
        )
        initial = None
        if args.initial_zones:
            import pandas as pd
            states = pd.read_csv(args.initial_zones)
            if "patient_id" not in states:
                raise ValueError("initial zone states need patient_id")
            patients = states.patient_id.astype(str).unique().tolist()
            patient_id = args.patient_id or (patients[0] if len(patients) == 1 else None)
            if patient_id is None:
                raise ValueError("--patient-id is required when initial states contain multiple patients")
            initial = states[states.patient_id.astype(str) == str(patient_id)]
            if initial.empty:
                raise ValueError(f"patient not found in initial zone states: {patient_id}")
        doses = tuple(Dose(float(day), args.dose)
                      for day in args.dose_days.split(",") if day.strip())
        if args.transport_only and doses:
            raise ValueError("anti-PD-1 dosing is incompatible with --transport-only")
        if args.paired_counterfactual and not doses:
            raise ValueError("--paired-counterfactual requires --dose-days")
        if args.paired_counterfactual:
            trajectory = simulate_spatial_counterfactual(
                doses, args.days, args.step, parameters=spatial_parameters,
                initial=initial, systemic_parameters=systemic_parameters,
            )
        elif doses:
            systemic = simulate(days=args.days, doses=doses,
                                parameters=systemic_parameters)
            trajectory = simulate_spatial_coupled(
                systemic, args.days, args.step, parameters=spatial_parameters,
                initial=initial,
            )
            trajectory["pd1_occupancy"] = np.interp(
                trajectory.day, systemic.day, systemic.pd1_occupancy
            )
            trajectory["arm"] = "anti_pd1"
        else:
            trajectory = simulate_spatial(args.days, args.step, parameters=spatial_parameters,
                                          initial=initial,
                                          reactions=not args.transport_only)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        trajectory.to_csv(args.output, index=False)
        _write_parameter_provenance(args.output, {
            "systemic": systemic_provenance, "spatial": spatial_provenance,
        })
        efferent = trajectory[(trajectory.species == "effector") & (trajectory.zone == "efferent_lymph")]
        report = {"rows": len(trajectory), "zones": int(trajectory.zone.nunique()),
                  "species": int(trajectory.species.nunique()), "output": args.output}
        if "arm" in efferent:
            report["final_efferent_effector_by_arm"] = (
                efferent.groupby("arm").amount.last().to_dict()
            )
        else:
            report["final_efferent_effector"] = float(efferent.amount.iloc[-1])
        print(json.dumps(report, indent=2))
    elif args.command == "spatial-featurize":
        import pandas as pd
        states = aggregate_spatial_spots(pd.read_csv(args.spots), args.weight_column)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        states.to_csv(args.output, index=False)
        write_spatial_feature_audit(args.spots, args.output, args.weight_column)
        observed_zones = states.loc[states.zone_observed].groupby("patient_id").zone.nunique()
        print(json.dumps({"patients": int(states.patient_id.nunique()),
                          "zones_in_state_grid": int(states.zone.nunique()),
                          "minimum_observed_zones": int(observed_zones.min()),
                          "maximum_observed_zones": int(observed_zones.max()),
                          "species": int(states.species.nunique()),
                          "output": args.output, "labels_used": False}, indent=2))
    elif args.command == "sensitivity":
        from .sensitivity import DEFAULT_PARAMETERS
        names = tuple(name.strip() for name in args.parameters.split(",") if name.strip()) or DEFAULT_PARAMETERS
        doses = (None if args.dose_days is None else tuple(
            Dose(float(day), args.dose)
            for day in args.dose_days.split(",") if day.strip()
        ))
        effective_doses = (doses if doses is not None else tuple(
            Dose(day, 1.0) for day in (0.0, 14.0, 28.0) if day < args.days
        ))
        base_parameters, parameter_provenance = _load_parameter_file(
            args.base_parameters, TwinParameters
        )
        samples, ranking = global_sensitivity(parameters=base_parameters, names=names,
                                              days=args.days, doses=doses,
                                              n=args.samples, fold_range=args.fold_range, seed=args.seed)
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        samples.to_csv(output / "sensitivity_samples.csv", index=False)
        ranking.to_csv(output / "sensitivity_ranking.csv", index=False)
        write_sensitivity_audit(
            output, "sensitivity_provenance.json", "systemic_sensitivity",
            inputs={"parameter_source": parameter_provenance},
            settings={
                "sampled_parameters": list(names), "days": args.days,
                "dose_days": [dose.day for dose in effective_doses],
                "dose_amounts": [dose.amount for dose in effective_doses],
                "samples": args.samples, "fold_range": args.fold_range,
                "seed": args.seed,
            },
            generated_filenames=("sensitivity_samples.csv", "sensitivity_ranking.csv"),
        )
        top = ranking.groupby("outcome", sort=False).first().reset_index()
        print(json.dumps({"samples": len(samples), "parameters": len(names),
                          "top_effects": top[["outcome", "parameter", "spearman_rho"]].to_dict("records"),
                          "output": str(output), "method": "LHS Spearman screening"}, indent=2))
    elif args.command == "spatial-sensitivity":
        from .sensitivity import DEFAULT_SPATIAL_PARAMETERS
        names = (tuple(name.strip() for name in args.parameters.split(",") if name.strip())
                 or DEFAULT_SPATIAL_PARAMETERS)
        base_parameters, parameter_provenance = _load_parameter_file(
            args.spatial_parameters, SpatialParameters
        )
        systemic_parameters, systemic_provenance = _load_parameter_file(
            args.systemic_parameters, TwinParameters
        )
        initial, initial_provenance = _load_patient_spatial_initial(
            args.initial_zones, args.patient_id
        )
        dose_days = tuple(
            float(day) for day in args.dose_days.split(",") if day.strip()
        )
        if args.systemic_parameters and not dose_days:
            raise ValueError("--systemic-parameters requires --dose-days")
        systemic_trajectory = None
        if dose_days:
            doses = tuple(Dose(day, args.dose) for day in dose_days)
            systemic_trajectory = simulate(
                days=args.days, doses=doses, parameters=systemic_parameters,
            )
        samples, ranking = global_spatial_sensitivity(
            parameters=base_parameters, initial=initial, names=names, days=args.days,
            checkpoint_occupancy=args.checkpoint_occupancy,
            systemic_trajectory=systemic_trajectory,
            n=args.samples, fold_range=args.fold_range, seed=args.seed,
        )
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        samples.to_csv(output / "spatial_sensitivity_samples.csv", index=False)
        ranking.to_csv(output / "spatial_sensitivity_ranking.csv", index=False)
        write_sensitivity_audit(
            output, "spatial_sensitivity_provenance.json", "spatial_sensitivity",
            inputs={
                "parameter_source": parameter_provenance,
                "systemic_parameter_source": systemic_provenance,
                "initial_state_source": initial_provenance,
            },
            settings={
                "sampled_parameters": list(names), "days": args.days,
                "checkpoint_mode": (
                    "dose_driven_systemic_trajectory" if systemic_trajectory is not None
                    else "fixed_occupancy"
                ),
                "checkpoint_occupancy": (
                    None if systemic_trajectory is not None else args.checkpoint_occupancy
                ),
                "dose_days": list(dose_days), "dose": args.dose,
                "peak_pd1_occupancy": (
                    float(systemic_trajectory.pd1_occupancy.max())
                    if systemic_trajectory is not None else None
                ),
                "samples": args.samples, "fold_range": args.fold_range,
                "seed": args.seed,
            },
            generated_filenames=(
                "spatial_sensitivity_samples.csv",
                "spatial_sensitivity_ranking.csv",
            ),
        )
        top = ranking.groupby("outcome", sort=False).first().reset_index()
        print(json.dumps({
            "samples": len(samples), "parameters": len(names),
            "checkpoint_mode": (
                "dose_driven_systemic_trajectory" if systemic_trajectory is not None
                else "fixed_occupancy"
            ),
            "peak_pd1_occupancy": (
                float(systemic_trajectory.pd1_occupancy.max())
                if systemic_trajectory is not None else args.checkpoint_occupancy
            ),
            "top_effects": top[["outcome", "parameter", "spearman_rho"]].to_dict("records"),
            "output": str(output), "method": "paired six-zone LHS Spearman screening",
        }, indent=2))
    elif args.command == "qc":
        import pandas as pd
        baseline_values = tuple(value.strip() for value in args.baseline_values.split(",") if value.strip())
        report = cohort_qc(load_cells(
            args.cells, (args.time_column,) if args.time_column else (),
        ), pd.read_csv(args.labels),
                           args.time_column, baseline_values or None,
                           args.feature_tissue_pattern,
                           args.min_signature_gene_coverage)
        if args.manifest:
            report["provenance"] = verify_manifest(args.manifest)
            report["passed"] = report["passed"] and report["provenance"]["valid"]
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        if not report["passed"]:
            raise SystemExit(2)
    elif args.command == "verify-sensitivity":
        input_paths = {}
        for key, value in (
            ("parameter_source", args.base_parameters or args.spatial_parameters),
            ("systemic_parameter_source", args.systemic_parameters),
            ("initial_state_source", args.initial_zones),
        ):
            if value:
                input_paths[key] = value
        report = verify_sensitivity_audit(args.audit, input_paths)
        print(json.dumps({
            "valid": report["valid"], "stage": report["stage"],
            "verified_inputs": report["verified_inputs"],
            "verified_outputs": report["verified_outputs"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-calibration":
        input_paths = {"observations": args.observations}
        if args.initial_zones:
            input_paths["initial_state_source"] = args.initial_zones
        if args.base_mapping:
            input_paths["base_mapping"] = args.base_mapping
        report = verify_calibration_audit(args.audit, input_paths)
        print(json.dumps({
            "valid": report["valid"], "stage": report["stage"],
            "verified_inputs": report["verified_inputs"],
            "verified_outputs": report["verified_outputs"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-augmentation":
        input_paths = {"features": args.features}
        for key, value in (
            ("spatial_initial_zones", args.spatial_initial_zones),
            ("systemic_parameters", args.systemic_parameters),
            ("spatial_parameters", args.spatial_parameters),
            ("feature_state_mapping", args.feature_state_mapping),
            ("systemic_diagnostics", args.systemic_diagnostics),
            ("spatial_diagnostics", args.spatial_diagnostics),
        ):
            if value:
                input_paths[key] = value
        report = verify_augmentation_audit(
            args.audit, args.output, input_paths,
        )
        print(json.dumps({
            "valid": report["valid"],
            "verified_inputs": report["verified_inputs"],
            "output_sha256": report["output_sha256"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-features":
        report = verify_feature_audit(args.audit, args.cells, args.output)
        print(json.dumps({
            "valid": report["valid"],
            "input_sha256": report["input_sha256"],
            "output_sha256": report["output_sha256"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-spatial-features":
        report = verify_spatial_feature_audit(args.audit, args.spots, args.output)
        print(json.dumps({
            "valid": report["valid"],
            "input_sha256": report["input_sha256"],
            "output_sha256": report["output_sha256"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-uncertainty":
        report = verify_uncertainty_audit(args.audit, {
            "parameters": args.parameters, "diagnostics": args.diagnostics,
        })
        print(json.dumps({
            "valid": report["valid"],
            "verified_inputs": report["verified_inputs"],
            "verified_outputs": report["verified_outputs"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    elif args.command == "verify-spatial-uncertainty":
        input_paths = {
            "parameters": args.parameters, "diagnostics": args.diagnostics,
        }
        if args.systemic_parameters:
            input_paths["systemic_parameters"] = args.systemic_parameters
        if args.initial_zones:
            input_paths["initial_zones"] = args.initial_zones
        report = verify_uncertainty_audit(args.audit, input_paths)
        if report["stage"] != "spatial_calibration_uncertainty_propagation":
            raise ValueError("expected spatial uncertainty audit")
        print(json.dumps({
            "valid": report["valid"],
            "verified_inputs": report["verified_inputs"],
            "verified_outputs": report["verified_outputs"],
            "audit_sha256": report["audit_sha256"],
        }, indent=2))
    else:
        report = verify_twin()
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        if not report["passed"]:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
