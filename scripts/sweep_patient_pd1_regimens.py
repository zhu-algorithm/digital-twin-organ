"""Screen research-only anti-PD-1 schedules in one patient digital twin."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_bridge_features import (load_feature_state_mapping,
                                             load_liver_bridge_parameters,
                                             load_spatial_parameters,
                                             load_systemic_parameters)
from tdln_pd1.liver_mapping_assimilation import load_feature_bridge_mapping
from tdln_pd1.regimen_sweep import (sweep_patient_regimens,
                                    verify_regimen_sweep_audit,
                                    write_regimen_sweep_audit)


def _grid(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(",") if item.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dose-amounts", default="0.5,1,2")
    parser.add_argument("--intervals", default="7,14")
    parser.add_argument("--days", type=float, default=28.0)
    parser.add_argument("--step", type=float, default=0.25)
    parser.add_argument("--portal-flow-multiplier", type=float)
    parser.add_argument("--arterial-flow-multiplier", type=float)
    parser.add_argument("--systemic-parameters", type=Path)
    parser.add_argument("--spatial-parameters", type=Path)
    parser.add_argument("--liver-bridge-parameters", type=Path)
    parser.add_argument("--feature-bridge-mapping", type=Path)
    parser.add_argument("--feature-state-mapping", type=Path)
    args = parser.parse_args()
    features = pd.read_csv(args.features)
    if "patient_id" not in features:
        raise ValueError("features need patient_id")
    selected = features[features.patient_id.astype(str).eq(str(args.patient_id))]
    if len(selected) != 1:
        raise ValueError("patient-id must select exactly one feature row")
    immune = load_systemic_parameters(args.systemic_parameters) if args.systemic_parameters else None
    spatial = load_spatial_parameters(args.spatial_parameters) if args.spatial_parameters else None
    state_mapping = load_feature_state_mapping(args.feature_state_mapping) if args.feature_state_mapping else None
    mapping = load_feature_bridge_mapping(args.feature_bridge_mapping) if args.feature_bridge_mapping else None
    if args.liver_bridge_parameters:
        liver, coupling, bridge = load_liver_bridge_parameters(args.liver_bridge_parameters)
    else:
        liver = coupling = bridge = None
    settings = {
        "dose_amounts": list(_grid(args.dose_amounts)),
        "intervals": list(_grid(args.intervals)),
        "days": args.days, "step": args.step,
        "portal_flow_multiplier": args.portal_flow_multiplier,
        "arterial_flow_multiplier": args.arterial_flow_multiplier,
    }
    result = sweep_patient_regimens(
        selected.iloc[0], dose_amounts=settings["dose_amounts"],
        intervals=settings["intervals"], days=args.days, step=args.step,
        mapping=mapping, base_coupling=coupling, base_bridge=bridge,
        base_spatial=spatial, base_liver=liver, base_immune=immune,
        feature_state_mapping=state_mapping,
        portal_flow_multiplier=args.portal_flow_multiplier,
        arterial_flow_multiplier=args.arterial_flow_multiplier,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    audit = Path(f"{args.output}.audit.json")
    write_regimen_sweep_audit(
        audit, features_path=args.features, output_path=args.output,
        patient_id=args.patient_id, settings=settings, mapping=mapping,
        base_coupling=coupling, base_bridge=bridge, base_spatial=spatial,
        base_liver=liver, base_immune=immune, feature_state_mapping=state_mapping,
        source_paths={
            "systemic_parameters": args.systemic_parameters,
            "spatial_parameters": args.spatial_parameters,
            "liver_bridge_parameters": args.liver_bridge_parameters,
            "feature_bridge_mapping": args.feature_bridge_mapping,
            "feature_state_mapping": args.feature_state_mapping,
        },
    )
    verify_regimen_sweep_audit(
        audit, args.features, args.output,
        source_paths={
            "systemic_parameters": args.systemic_parameters,
            "spatial_parameters": args.spatial_parameters,
            "liver_bridge_parameters": args.liver_bridge_parameters,
            "feature_bridge_mapping": args.feature_bridge_mapping,
            "feature_state_mapping": args.feature_state_mapping,
        },
    )
    print(json.dumps({"regimens": len(result),
                      "pareto_regimens": int(result.pareto_front.sum()),
                      "labels_used": False, "research_only": True}, indent=2))


if __name__ == "__main__":
    main()
