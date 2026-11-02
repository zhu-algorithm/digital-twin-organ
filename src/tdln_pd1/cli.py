from __future__ import annotations

import argparse
from dataclasses import asdict, fields
import hashlib
import json
from pathlib import Path
import numpy as np
from .data import load_cells, load_labels
from .features import aggregate_cells
from .model import evaluate_external, leave_one_site_out, predict_locked, train_evaluate
from .simulator import Dose, TwinParameters, simulate
from .assimilation import calibrate, save_calibration
from .hybrid import augment_with_twin
from .spatial_dynamics import (SpatialParameters, simulate_spatial, simulate_spatial_coupled,
                               simulate_spatial_counterfactual)
from .spatial_data import aggregate_spatial_spots
from .spatial_assimilation import calibrate_spatial, save_spatial_calibration
from .sensitivity import global_sensitivity, global_spatial_sensitivity
from .sensitivity_audit import verify_sensitivity_audit, write_sensitivity_audit
from .provenance import verify_manifest
from .qc import cohort_qc
from .verification import verify_twin
from .validation_audit import (verify_prediction_audit, write_evaluation_audit,
                               write_prediction_audit)
from .calibration_audit import verify_calibration_audit, write_calibration_audit


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
    return selected, {
        "file": source.name,
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "patient_id": selected_id,
        "rows": int(len(selected)),
    }


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
        for name, entry in payload["sources"].items():
            digest = entry.get("sha256") if isinstance(entry, dict) else None
            if (name not in {"systemic", "spatial"} or not isinstance(digest, str)
                    or len(digest) != 64
                    or any(character not in "0123456789abcdef" for character in digest.lower())):
                raise ValueError("invalid parameter source fingerprint in feature sidecar")
        provenance["parameter_provenance"] = payload
        provenance["parameter_sidecar_sha256"] = hashlib.sha256(sidecar.read_bytes()).hexdigest()
    else:
        provenance["parameter_provenance"] = None
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
    aug = sub.add_parser("augment-twin")
    aug.add_argument("--features", required=True)
    aug.add_argument("--days", type=float, default=42.0)
    aug.add_argument("--dose-days", default="0,14,28")
    aug.add_argument("--dose", type=float, default=1.0)
    aug.add_argument("--spatial-initial-zones",
                     help="patient-zone-species states from spatial-featurize")
    aug.add_argument("--systemic-parameters", help="fitted TwinParameters JSON")
    aug.add_argument("--spatial-parameters", help="fitted SpatialParameters JSON")
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
    args = parser.parse_args()
    if args.command == "featurize":
        baseline_values = tuple(value.strip() for value in args.baseline_values.split(",") if value.strip())
        aggregate_cells(load_cells(args.cells, (args.time_column,) if args.time_column else ()),
                        args.time_column, baseline_values or None,
                        args.feature_tissue_pattern).to_csv(args.output, index=False)
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
        augmented = augment_with_twin(pd.read_csv(args.features), args.days, dose_days,
                                      args.dose, spatial_states,
                                      systemic_parameters=systemic_parameters,
                                      spatial_parameters=spatial_parameters)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        augmented.to_csv(args.output, index=False)
        _write_parameter_provenance(args.output, {
            "systemic": systemic_provenance, "spatial": spatial_provenance,
        })
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
                           args.time_column, baseline_values or None, args.feature_tissue_pattern)
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
        report = verify_calibration_audit(args.audit, input_paths)
        print(json.dumps({
            "valid": report["valid"], "stage": report["stage"],
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
