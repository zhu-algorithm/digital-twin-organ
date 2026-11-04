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
    FeatureBridgeMapping,
    augment_with_liver_tdln_counterfactuals,
    liver_bridge_parameter_draws,
    load_liver_bridge_parameters,
    validate_perfusion_manifest,
    validate_treatment_manifest,
    verify_liver_counterfactual_augmentation,
)
from tdln_pd1.provenance import sha256_file
from tdln_pd1.liver_tdln_coupling import LiverTDLNCouplingParameters
from tdln_pd1.liver_tdln_spatial import LiverTDLNSpatialBridge
from tdln_pd1.simulator import Dose


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--days", type=float, default=28.0)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--dose-days")
    parser.add_argument("--dose", type=float, default=1.0)
    parser.add_argument("--bridge-parameters", type=Path)
    parser.add_argument("--bridge-diagnostics", type=Path)
    parser.add_argument("--uncertainty-samples", type=int, default=0)
    parser.add_argument("--uncertainty-seed", type=int, default=42)
    parser.add_argument("--treatment-manifest", type=Path)
    parser.add_argument("--perfusion-manifest", type=Path)
    args = parser.parse_args()
    if bool(args.dose_days) == bool(args.treatment_manifest):
        raise ValueError("provide exactly one of --dose-days or --treatment-manifest")
    features = pd.read_csv(args.features)
    doses = tuple(Dose(float(day), args.dose)
                  for day in (args.dose_days or "").split(",") if day.strip())
    patient_doses = treatment_contract = None
    if args.treatment_manifest:
        patient_doses, treatment_contract = validate_treatment_manifest(
            pd.read_csv(args.treatment_manifest),
            features.patient_id.astype(str).tolist(), args.days)
    patient_flows = perfusion_contract = None
    if args.perfusion_manifest:
        patient_flows, perfusion_contract = validate_perfusion_manifest(
            pd.read_csv(args.perfusion_manifest),
            features.patient_id.astype(str).tolist())
    mapping = FeatureBridgeMapping()
    if args.bridge_parameters:
        base_coupling, base_bridge = load_liver_bridge_parameters(args.bridge_parameters)
    else:
        base_coupling, base_bridge = None, None
    if args.uncertainty_samples < 0 or args.uncertainty_samples == 1:
        raise ValueError("--uncertainty-samples must be 0 or at least 2")
    if bool(args.bridge_diagnostics) != bool(args.uncertainty_samples):
        raise ValueError("--bridge-diagnostics and --uncertainty-samples require each other")
    if args.bridge_diagnostics and not args.bridge_parameters:
        raise ValueError("--bridge-diagnostics requires --bridge-parameters")
    parameter_draws = None
    if args.bridge_diagnostics:
        parameter_draws = liver_bridge_parameter_draws(
            base_coupling, base_bridge,
            json.loads(args.bridge_diagnostics.read_text()),
            args.uncertainty_samples, args.uncertainty_seed)
    augmented = augment_with_liver_tdln_counterfactuals(
        features, doses=doses, days=args.days,
        step=args.step, mapping=mapping, base_coupling=base_coupling,
        base_bridge=base_bridge, parameter_draws=parameter_draws,
        patient_doses=patient_doses, patient_flows=patient_flows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(args.output, index=False)
    mechanistic_columns = [
        column for column in augmented.columns if column.startswith("liver_cf_")
    ]
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
                     "paired_control": True,
                     "uncertainty_samples": args.uncertainty_samples,
                     "uncertainty_seed": args.uncertainty_seed},
        "mapping": asdict(mapping),
        "base_parameters": {
            "coupling": asdict(base_coupling or LiverTDLNCouplingParameters()),
            "spatial_bridge": asdict(base_bridge or LiverTDLNSpatialBridge()),
        },
        "bridge_parameters": ({"file": args.bridge_parameters.name,
                                "sha256": sha256_file(args.bridge_parameters)}
                               if args.bridge_parameters else None),
        "bridge_diagnostics": ({"file": args.bridge_diagnostics.name,
                                 "sha256": sha256_file(args.bridge_diagnostics)}
                                if args.bridge_diagnostics else None),
        "treatment_manifest": ({"file": args.treatment_manifest.name,
                                 "sha256": sha256_file(args.treatment_manifest)}
                                if args.treatment_manifest else None),
        "treatment_contract": treatment_contract,
        "perfusion_manifest": ({"file": args.perfusion_manifest.name,
                                 "sha256": sha256_file(args.perfusion_manifest)}
                                if args.perfusion_manifest else None),
        "perfusion_contract": perfusion_contract,
        "mechanistic_feature_columns": mechanistic_columns,
        "limitation": "Model-derived matched contrasts; not clinical response labels.",
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    verify_liver_counterfactual_augmentation(
        audit_path, args.features, args.output, args.bridge_parameters,
        args.bridge_diagnostics, args.treatment_manifest, args.perfusion_manifest)
    print(json.dumps({"n_patients": len(augmented),
                      "n_mechanistic_features": len(mechanistic_columns),
                      "labels_used": False, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
