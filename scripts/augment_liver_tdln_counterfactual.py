"""Append matched liver–TDLN anti-PD-1 contrasts to a patient feature table."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_bridge_features import (
    COUNTERFACTUAL_FEATURES, FeatureBridgeMapping,
    augment_with_liver_tdln_counterfactuals,
    verify_liver_counterfactual_augmentation,
)
from tdln_pd1.provenance import sha256_file
from tdln_pd1.simulator import Dose


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--days", type=float, default=28.0)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--dose-days", required=True)
    parser.add_argument("--dose", type=float, default=1.0)
    args = parser.parse_args()
    doses = tuple(Dose(float(day), args.dose)
                  for day in args.dose_days.split(",") if day.strip())
    mapping = FeatureBridgeMapping()
    augmented = augment_with_liver_tdln_counterfactuals(
        pd.read_csv(args.features), doses=doses, days=args.days,
        step=args.step, mapping=mapping)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(args.output, index=False)
    audit_path = Path(f"{args.output}.liver-counterfactual-audit.json")
    audit = {
        "schema_version": 1,
        "stage": "label_free_liver_tdln_counterfactual_augmentation",
        "labels_used": False,
        "research_only": True,
        "input": {"file": args.features.name, "sha256": sha256_file(args.features)},
        "output": {"file": args.output.name, "sha256": sha256_file(args.output)},
        "settings": {"days": args.days, "step": args.step,
                     "dose_days": [dose.day for dose in doses],
                     "dose_amounts": [dose.amount for dose in doses],
                     "paired_control": True},
        "mapping": asdict(mapping),
        "mechanistic_feature_columns": list(COUNTERFACTUAL_FEATURES.values()),
        "limitation": "Model-derived matched contrasts; not clinical response labels.",
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    verify_liver_counterfactual_augmentation(
        audit_path, args.features, args.output)
    print(json.dumps({"n_patients": len(augmented),
                      "n_mechanistic_features": len(COUNTERFACTUAL_FEATURES),
                      "labels_used": False, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
