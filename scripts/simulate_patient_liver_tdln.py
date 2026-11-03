"""Run the liver-tumour–TDLN closed loop from one label-free feature row."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_bridge_features import (simulate_patient_liver_tdln,
                                             verify_patient_liver_tdln_audit)
from tdln_pd1.provenance import sha256_file
from tdln_pd1.simulator import Dose


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--days", type=float, default=28.0)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--dose-days", default="")
    parser.add_argument("--dose", type=float, default=1.0)
    args = parser.parse_args()
    features = pd.read_csv(args.features)
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    matches = features[features.patient_id.astype(str).eq(str(args.patient_id))]
    if len(matches) != 1:
        raise ValueError("patient-id must select exactly one feature row")
    doses = tuple(Dose(float(day), args.dose)
                  for day in args.dose_days.split(",") if day.strip())
    systemic, spatial, audit = simulate_patient_liver_tdln(
        matches.iloc[0], days=args.days, step=args.step, doses=doses)
    args.output.mkdir(parents=True, exist_ok=True)
    systemic_path = args.output / "patient_liver_tdln_systemic.csv"
    spatial_path = args.output / "patient_liver_tdln_spatial.csv"
    systemic.to_csv(systemic_path, index=False)
    spatial.to_csv(spatial_path, index=False)
    audit.update({
        "patient_id": str(args.patient_id),
        "input": {"file": args.features.name, "sha256": sha256_file(args.features)},
        "settings": {"days": args.days, "step": args.step,
                     "dose_days": [dose.day for dose in doses], "dose": args.dose},
        "outputs": {
            systemic_path.name: sha256_file(systemic_path),
            spatial_path.name: sha256_file(spatial_path),
        },
    })
    (args.output / "patient_liver_tdln_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n")
    verify_patient_liver_tdln_audit(args.output, args.features)
    print(json.dumps({"patient_id": str(args.patient_id),
                      "final_liver_tumour_fraction": float(systemic.liver_tumour_fraction.iloc[-1]),
                      "peak_pd1_occupancy": float(systemic.pd1_occupancy.max()),
                      "labels_used": False, "research_only": True}, indent=2))


if __name__ == "__main__":
    main()

