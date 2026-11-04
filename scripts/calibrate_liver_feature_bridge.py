"""Calibrate single-cell features to liver–TDLN bridge multipliers."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.calibration_audit import write_calibration_audit
from tdln_pd1.liver_mapping_assimilation import (
    calibrate_feature_bridge_mapping, save_feature_bridge_mapping_calibration)
from tdln_pd1.provenance import sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fit", required=True)
    parser.add_argument("--max-nfev", type=int, default=200)
    args = parser.parse_args()
    fit = tuple(value.strip() for value in args.fit.split(",") if value.strip())
    mapping, diagnostics, residuals = calibrate_feature_bridge_mapping(
        pd.read_csv(args.observations), fit_parameters=fit, max_nfev=args.max_nfev)
    save_feature_bridge_mapping_calibration(
        args.output, mapping, diagnostics, residuals)
    write_calibration_audit(
        args.output, "liver_feature_bridge_mapping_provenance.json",
        "liver_feature_bridge_mapping_calibration",
        inputs={"observations": {"file": args.observations.name,
                                  "sha256": sha256_file(args.observations)}},
        settings={"fit_parameters": list(fit), "max_nfev": args.max_nfev,
                  "labels_used": False},
        generated_filenames=("fitted_liver_feature_bridge_mapping.json",
                             "liver_feature_bridge_mapping_diagnostics.json",
                             "liver_feature_bridge_mapping_residuals.csv"),
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
