from __future__ import annotations

import argparse
from dataclasses import fields
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
from .sensitivity import global_sensitivity
from .provenance import verify_manifest
from .qc import cohort_qc
from .verification import verify_twin


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
    spatial_cal = sub.add_parser("calibrate-spatial")
    spatial_cal.add_argument("--observations", required=True)
    spatial_cal.add_argument("--fit", required=True,
                             help="comma-separated SpatialParameters fields")
    spatial_cal.add_argument("--initial-zones")
    spatial_cal.add_argument("--patient-id")
    spatial_cal.add_argument("--dose-days", default="")
    spatial_cal.add_argument("--dose", type=float, default=1.0)
    spatial_cal.add_argument("--max-nfev", type=int, default=200)
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
    sens.add_argument("--days", type=float, default=42.0)
    sens.add_argument("--dose-days", default="0,14,28")
    sens.add_argument("--dose", type=float, default=1.0)
    sens.add_argument("--samples", type=int, default=128)
    sens.add_argument("--fold-range", type=float, default=4.0)
    sens.add_argument("--seed", type=int, default=42)
    sens.add_argument("--output", required=True)
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
    args = parser.parse_args()
    if args.command == "featurize":
        baseline_values = tuple(value.strip() for value in args.baseline_values.split(",") if value.strip())
        aggregate_cells(load_cells(args.cells), args.time_column, baseline_values or None,
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
        print(json.dumps({"n_patients": len(predictions), "output": args.output,
                          "model_refit": False}, indent=2))
    elif args.command == "evaluate":
        import pandas as pd
        metrics = evaluate_external(pd.read_csv(args.predictions), load_labels(args.labels), args.bootstrap)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(metrics, indent=2) + "\n")
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
        doses = tuple(Dose(float(day), args.dose) for day in args.dose_days.split(",") if day.strip())
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
        )
        save_calibration(args.output, fitted, diagnostics, compared)
        print(json.dumps(diagnostics, indent=2))
    elif args.command == "calibrate-spatial":
        import pandas as pd
        observations = pd.read_csv(args.observations)
        initial = None
        if args.initial_zones:
            states = pd.read_csv(args.initial_zones)
            if "patient_id" not in states:
                raise ValueError("initial zone states need patient_id")
            patients = states.patient_id.astype(str).unique().tolist()
            patient_id = args.patient_id or (patients[0] if len(patients) == 1 else None)
            if patient_id is None:
                raise ValueError("--patient-id is required when initial states contain multiple patients")
            initial = states[states.patient_id.astype(str).eq(str(patient_id))]
            if initial.empty:
                raise ValueError(f"patient not found in initial zone states: {patient_id}")
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
        )
        save_spatial_calibration(args.output, fitted, diagnostics, compared)
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
        print(json.dumps({"patients": int(states.patient_id.nunique()),
                          "zones_observed": int(states.zone.nunique()),
                          "species": int(states.species.nunique()),
                          "output": args.output, "labels_used": False}, indent=2))
    elif args.command == "sensitivity":
        from .sensitivity import DEFAULT_PARAMETERS
        names = tuple(name.strip() for name in args.parameters.split(",") if name.strip()) or DEFAULT_PARAMETERS
        doses = tuple(Dose(float(day), args.dose) for day in args.dose_days.split(",") if day.strip())
        samples, ranking = global_sensitivity(names=names, days=args.days, doses=doses,
                                              n=args.samples, fold_range=args.fold_range, seed=args.seed)
        output = Path(args.output)
        output.mkdir(parents=True, exist_ok=True)
        samples.to_csv(output / "sensitivity_samples.csv", index=False)
        ranking.to_csv(output / "sensitivity_ranking.csv", index=False)
        top = ranking.groupby("outcome", sort=False).first().reset_index()
        print(json.dumps({"samples": len(samples), "parameters": len(names),
                          "top_effects": top[["outcome", "parameter", "spearman_rho"]].to_dict("records"),
                          "output": str(output), "method": "LHS Spearman screening"}, indent=2))
    elif args.command == "qc":
        import pandas as pd
        baseline_values = tuple(value.strip() for value in args.baseline_values.split(",") if value.strip())
        report = cohort_qc(load_cells(args.cells), pd.read_csv(args.labels),
                           args.time_column, baseline_values or None, args.feature_tissue_pattern)
        if args.manifest:
            report["provenance"] = verify_manifest(args.manifest)
            report["passed"] = report["passed"] and report["provenance"]["valid"]
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        if not report["passed"]:
            raise SystemExit(2)
    else:
        report = verify_twin()
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        if not report["passed"]:
            raise SystemExit(2)


if __name__ == "__main__":
    main()
