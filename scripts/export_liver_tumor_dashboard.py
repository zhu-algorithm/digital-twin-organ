"""Export deterministic trajectories for the liver-tumour virtual organ."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_tumor import parameter_manifest, simulate_liver_tumour

SCHEMA_VERSION = 1
COLUMNS = ("normal_cells", "tumour_rim", "tumour_core", "necrotic_tumour",
           "functional_vessels", "oxygen_normal", "oxygen_rim", "oxygen_core",
           "perfusion", "portal_flow", "arterial_flow", "tumour_fraction")


def _records(frame):
    days = np.arange(0.0, 120.01, 2.0)
    return [{"day": float(day), **{name: round(float(np.interp(day, frame.day, frame[name])), 6)
                                    for name in COLUMNS}} for day in days]


def build_dashboard_data() -> dict:
    scenarios = {
        "untreated": simulate_liver_tumour(),
        "antiangiogenic": simulate_liver_tumour(angiogenesis_inhibition=0.75),
        "low_perfusion": simulate_liver_tumour(flow_multiplier=0.65),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "model_version": "liver-tumour/0.1.0",
        "research_only": True,
        "time_unit": "day",
        "state_units": "normalised fraction",
        "scenarios": {name: _records(frame) for name, frame in scenarios.items()},
        "parameter_manifest": parameter_manifest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("docs/data/liver_tumor_dashboard.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(build_dashboard_data(), ensure_ascii=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()

