"""Fit liver-tumour ↔ TDLN bridge parameters from longitudinal observations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.calibration_audit import write_calibration_audit
from tdln_pd1.liver_bridge_assimilation import (calibrate_liver_bridge,
                                                 save_liver_bridge_calibration)
from tdln_pd1.provenance import sha256_file
from tdln_pd1.simulator import Dose


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--fit",
        default="bridge.antigen_drainage_gain,bridge.dc_migration_gain",
    )
    parser.add_argument("--dose-days", default="")
    parser.add_argument("--dose", type=float, default=1.0)
    parser.add_argument("--step", type=float, default=0.5)
    parser.add_argument("--portal-flow-multiplier", type=float)
    parser.add_argument("--arterial-flow-multiplier", type=float)
    parser.add_argument("--max-nfev", type=int, default=100)
    args = parser.parse_args()
    fit = tuple(name.strip() for name in args.fit.split(",") if name.strip())
    doses = tuple(Dose(float(day), args.dose)
                  for day in args.dose_days.split(",") if day.strip())
    liver, coupling, bridge, diagnostics, compared = calibrate_liver_bridge(
        pd.read_csv(args.observations), fit_parameters=fit, doses=doses,
        step=args.step, max_nfev=args.max_nfev,
        portal_flow_multiplier=args.portal_flow_multiplier,
        arterial_flow_multiplier=args.arterial_flow_multiplier)
    save_liver_bridge_calibration(
        args.output, liver, coupling, bridge, diagnostics, compared)
    write_calibration_audit(
        args.output, "liver_tdln_bridge_provenance.json",
        "liver_tdln_bridge_calibration",
        inputs={"observations": {
            "file": args.observations.name, "sha256": sha256_file(args.observations),
        }},
        settings={"fit_parameters": list(fit), "dose_days": [d.day for d in doses],
                  "dose_amount": args.dose, "step": args.step,
                  "portal_flow_multiplier": args.portal_flow_multiplier,
                  "arterial_flow_multiplier": args.arterial_flow_multiplier,
                  "max_nfev": args.max_nfev},
        generated_filenames=("fitted_liver_tdln_bridge.json",
                             "liver_tdln_bridge_diagnostics.json",
                             "liver_tdln_bridge_residuals.csv"),
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

