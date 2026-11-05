"""Run the liver-tumour–TDLN closed loop from one label-free feature row."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_bridge_features import (simulate_patient_liver_tdln,
                                             simulate_patient_liver_tdln_counterfactual,
                                             load_feature_state_mapping,
                                             load_liver_bridge_parameters,
                                             load_spatial_parameters,
                                             load_systemic_parameters,
                                             verify_patient_liver_tdln_audit)
from tdln_pd1.hybrid import FeatureStateMapping
from tdln_pd1.liver_mapping_assimilation import load_feature_bridge_mapping
from tdln_pd1.liver_tdln_coupling import LiverTDLNCouplingParameters
from tdln_pd1.liver_tdln_spatial import LiverTDLNSpatialBridge
from tdln_pd1.liver_tumor import LiverTumourParameters
from tdln_pd1.provenance import sha256_file
from tdln_pd1.simulator import Dose, TwinParameters
from tdln_pd1.spatial_dynamics import SpatialParameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--days", type=float, default=28.0)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--dose-days", default="")
    parser.add_argument("--dose", type=float, default=1.0)
    parser.add_argument("--systemic-parameters", type=Path)
    parser.add_argument("--spatial-parameters", type=Path)
    parser.add_argument("--liver-bridge-parameters", type=Path)
    parser.add_argument("--feature-bridge-mapping", type=Path)
    parser.add_argument("--feature-state-mapping", type=Path)
    parser.add_argument("--portal-flow-multiplier", type=float)
    parser.add_argument("--arterial-flow-multiplier", type=float)
    parser.add_argument("--paired-counterfactual", action="store_true")
    args = parser.parse_args()
    features = pd.read_csv(args.features)
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    matches = features[features.patient_id.astype(str).eq(str(args.patient_id))]
    if len(matches) != 1:
        raise ValueError("patient-id must select exactly one feature row")
    doses = tuple(Dose(float(day), args.dose)
                  for day in args.dose_days.split(",") if day.strip())
    systemic_parameters = (load_systemic_parameters(args.systemic_parameters)
                           if args.systemic_parameters else None)
    spatial_parameters = (load_spatial_parameters(args.spatial_parameters)
                          if args.spatial_parameters else None)
    feature_state_mapping = (load_feature_state_mapping(args.feature_state_mapping)
                             if args.feature_state_mapping else None)
    feature_bridge_mapping = (load_feature_bridge_mapping(args.feature_bridge_mapping)
                              if args.feature_bridge_mapping else None)
    if args.liver_bridge_parameters:
        liver_parameters, coupling_parameters, bridge_parameters = (
            load_liver_bridge_parameters(args.liver_bridge_parameters))
    else:
        liver_parameters = coupling_parameters = bridge_parameters = None
    simulation_kwargs = {
        "base_immune": systemic_parameters,
        "base_spatial": spatial_parameters,
        "base_liver": liver_parameters,
        "base_coupling": coupling_parameters,
        "base_bridge": bridge_parameters,
        "mapping": feature_bridge_mapping,
        "feature_state_mapping": feature_state_mapping,
    }
    if args.paired_counterfactual:
        systemic, spatial, effect, audit = simulate_patient_liver_tdln_counterfactual(
            matches.iloc[0], days=args.days, step=args.step, doses=doses,
            **simulation_kwargs,
            portal_flow_multiplier=args.portal_flow_multiplier,
            arterial_flow_multiplier=args.arterial_flow_multiplier)
    else:
        systemic, spatial, audit = simulate_patient_liver_tdln(
            matches.iloc[0], days=args.days, step=args.step, doses=doses,
            **simulation_kwargs,
            portal_flow_multiplier=args.portal_flow_multiplier,
            arterial_flow_multiplier=args.arterial_flow_multiplier)
        effect = None
    args.output.mkdir(parents=True, exist_ok=True)
    systemic_path = args.output / "patient_liver_tdln_systemic.csv"
    spatial_path = args.output / "patient_liver_tdln_spatial.csv"
    systemic.to_csv(systemic_path, index=False)
    spatial.to_csv(spatial_path, index=False)
    if effect is not None:
        effect_path = args.output / "patient_liver_tdln_counterfactual.json"
        effect_path.write_text(json.dumps(effect, indent=2, sort_keys=True) + "\n")
    audit.update({
        "patient_id": str(args.patient_id),
        "input": {"file": args.features.name, "sha256": sha256_file(args.features)},
        "systemic_parameters": (
            {"file": args.systemic_parameters.name,
             "sha256": sha256_file(args.systemic_parameters)}
            if args.systemic_parameters else None),
        "spatial_parameters": (
            {"file": args.spatial_parameters.name,
             "sha256": sha256_file(args.spatial_parameters)}
            if args.spatial_parameters else None),
        "liver_bridge_parameters": (
            {"file": args.liver_bridge_parameters.name,
             "sha256": sha256_file(args.liver_bridge_parameters)}
            if args.liver_bridge_parameters else None),
        "feature_bridge_mapping": (
            {"file": args.feature_bridge_mapping.name,
             "sha256": sha256_file(args.feature_bridge_mapping)}
            if args.feature_bridge_mapping else None),
        "feature_state_mapping": (
            {"file": args.feature_state_mapping.name,
             "sha256": sha256_file(args.feature_state_mapping)}
            if args.feature_state_mapping else None),
        "base_parameters": {
            "immune": asdict(systemic_parameters or TwinParameters()),
            "spatial": asdict(spatial_parameters or SpatialParameters()),
            "liver_tumour": asdict(liver_parameters or LiverTumourParameters()),
            "coupling": asdict(coupling_parameters or LiverTDLNCouplingParameters()),
            "spatial_bridge": asdict(bridge_parameters or LiverTDLNSpatialBridge()),
            "feature_state_mapping": asdict(
                feature_state_mapping or FeatureStateMapping()),
        },
        "settings": {"days": args.days, "step": args.step,
                     "dose_days": [dose.day for dose in doses], "dose": args.dose,
                     "portal_flow_multiplier": args.portal_flow_multiplier,
                     "arterial_flow_multiplier": args.arterial_flow_multiplier},
        "outputs": {
            systemic_path.name: sha256_file(systemic_path),
            spatial_path.name: sha256_file(spatial_path),
        },
    })
    audit["settings"]["paired_counterfactual"] = args.paired_counterfactual
    if effect is not None:
        audit["outputs"][effect_path.name] = sha256_file(effect_path)
    (args.output / "patient_liver_tdln_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n")
    verify_patient_liver_tdln_audit(
        args.output, args.features, args.systemic_parameters,
        args.spatial_parameters, args.liver_bridge_parameters,
        args.feature_bridge_mapping, args.feature_state_mapping)
    print(json.dumps({"patient_id": str(args.patient_id),
                      "final_liver_tumour_fraction": float(systemic.liver_tumour_fraction.iloc[-1]),
                      "peak_pd1_occupancy": float(systemic.pd1_occupancy.max()),
                      "paired_counterfactual": args.paired_counterfactual,
                      "mechanistic_relative_tumour_reduction": (
                          effect["final_liver_tumour_fraction_relative_reduction"]
                          if effect is not None else None),
                      "labels_used": False, "research_only": True}, indent=2))


if __name__ == "__main__":
    main()

