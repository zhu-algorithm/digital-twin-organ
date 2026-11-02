from __future__ import annotations

import argparse
import json
from pathlib import Path
from .data import load_cells, load_labels
from .features import aggregate_cells
from .model import evaluate_external, predict_locked, train_evaluate
from .simulator import Dose, simulate
from .assimilation import calibrate, save_calibration
from .hybrid import augment_with_twin
from .spatial_dynamics import simulate_spatial


def main() -> None:
    parser = argparse.ArgumentParser(description="TDLN PD-1 response baseline")
    sub = parser.add_subparsers(dest="command", required=True)
    featurize = sub.add_parser("featurize")
    featurize.add_argument("--cells", required=True)
    featurize.add_argument("--output", required=True)
    train = sub.add_parser("train")
    train.add_argument("--features", required=True)
    train.add_argument("--labels", required=True)
    train.add_argument("--output", required=True)
    train.add_argument("--folds", type=int, default=5)
    predict = sub.add_parser("predict")
    predict.add_argument("--model", required=True)
    predict.add_argument("--features", required=True)
    predict.add_argument("--output", required=True)
    evaluate = sub.add_parser("evaluate")
    evaluate.add_argument("--predictions", required=True)
    evaluate.add_argument("--labels", required=True)
    evaluate.add_argument("--output", required=True)
    evaluate.add_argument("--bootstrap", type=int, default=500)
    sim = sub.add_parser("simulate")
    sim.add_argument("--days", type=float, default=42.0)
    sim.add_argument("--dose-days", default="0,14,28")
    sim.add_argument("--dose", type=float, default=1.0)
    sim.add_argument("--output", required=True)
    cal = sub.add_parser("calibrate")
    cal.add_argument("--observations", required=True)
    cal.add_argument("--fit", required=True, help="comma-separated TwinParameters fields")
    cal.add_argument("--dose-days", default="")
    cal.add_argument("--dose", type=float, default=1.0)
    cal.add_argument("--output", required=True)
    aug = sub.add_parser("augment-twin")
    aug.add_argument("--features", required=True)
    aug.add_argument("--days", type=float, default=42.0)
    aug.add_argument("--dose-days", default="0,14,28")
    aug.add_argument("--dose", type=float, default=1.0)
    aug.add_argument("--output", required=True)
    spatial = sub.add_parser("simulate-spatial")
    spatial.add_argument("--days", type=float, default=14.0)
    spatial.add_argument("--step", type=float, default=0.1)
    spatial.add_argument("--transport-only", action="store_true")
    spatial.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "featurize":
        aggregate_cells(load_cells(args.cells)).to_csv(args.output, index=False)
    elif args.command == "train":
        import pandas as pd
        metrics = train_evaluate(pd.read_csv(args.features), load_labels(args.labels), Path(args.output), args.folds)
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
    elif args.command == "simulate":
        doses = tuple(Dose(float(day), args.dose) for day in args.dose_days.split(",") if day.strip())
        trajectory = simulate(days=args.days, doses=doses)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        trajectory.to_csv(args.output, index=False)
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
    elif args.command == "augment-twin":
        import pandas as pd
        dose_days = tuple(float(day) for day in args.dose_days.split(",") if day.strip())
        augmented = augment_with_twin(pd.read_csv(args.features), args.days, dose_days, args.dose)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        augmented.to_csv(args.output, index=False)
        print(json.dumps({"n_patients": len(augmented), "n_features": len(augmented.columns) - 1,
                          "output": args.output, "labels_used": False}, indent=2))
    else:
        trajectory = simulate_spatial(args.days, args.step, reactions=not args.transport_only)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        trajectory.to_csv(args.output, index=False)
        efferent = trajectory[(trajectory.species == "effector") & (trajectory.zone == "efferent_lymph")]
        print(json.dumps({"rows": len(trajectory), "zones": int(trajectory.zone.nunique()),
                          "species": int(trajectory.species.nunique()),
                          "final_efferent_effector": float(efferent.amount.iloc[-1]),
                          "output": args.output}, indent=2))


if __name__ == "__main__":
    main()
