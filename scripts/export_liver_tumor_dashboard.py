"""Export deterministic trajectories for the liver-tumour virtual organ."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.liver_tumor import parameter_manifest, simulate_liver_tumour
from tdln_pd1.liver_tdln_coupling import (LiverTDLNCouplingParameters,
                                          simulate_liver_tdln)
from tdln_pd1.simulator import Dose

SCHEMA_VERSION = 1
COLUMNS = ("normal_cells", "tumour_rim", "tumour_core", "necrotic_tumour",
           "functional_vessels", "oxygen_normal", "oxygen_rim", "oxygen_core",
           "perfusion", "portal_flow", "arterial_flow", "tumour_fraction")
IMMUNE_COLUMNS = ("tumour_antigen", "mature_dc_tdln", "tpex_tdln",
                  "effector_tdln", "effector_tumour", "pd1_occupancy",
                  "tumour_cells", "liver_tumour_fraction")


def _records(frame):
    days = np.arange(0.0, 120.01, 2.0)
    return [{"day": float(day), **{name: round(float(np.interp(day, frame.day, frame[name])), 6)
                                    for name in COLUMNS}} for day in days]


def _immune_records(frame):
    days = np.arange(0.0, 84.01, 2.0)
    # Dashboard values are deliberately quantised to their display resolution.
    # This avoids presenting solver/BLAS-scale digits as biological precision
    # and keeps the checked-in export byte-stable across the qualified matrix.
    count_columns = {"mature_dc_tdln", "tpex_tdln", "effector_tdln", "effector_tumour"}
    records = []
    for day in days:
        row = {"day": float(day)}
        for name in IMMUNE_COLUMNS:
            value = float(np.interp(day, frame.day, frame[name]))
            if name == "tumour_cells":
                value = round(value, -2)
            elif name in count_columns:
                value = round(value, -1)
            else:
                value = round(value, 4)
            row[name] = value
        records.append(row)
    return records


def build_dashboard_data() -> dict:
    scenarios = {
        "untreated": simulate_liver_tumour(),
        "antiangiogenic": simulate_liver_tumour(angiogenesis_inhibition=0.75),
        "low_perfusion": simulate_liver_tumour(flow_multiplier=0.65),
    }
    coupled = {
        "control": simulate_liver_tdln(days=84),
        "anti_pd1": simulate_liver_tdln(
            days=84, doses=tuple(Dose(day, 1.0) for day in (0, 21, 42, 63))),
    }
    return {
        "schema_version": SCHEMA_VERSION,
        "model_version": "liver-tumour/0.1.0",
        "research_only": True,
        "time_unit": "day",
        "state_units": "normalised fraction",
        "scenarios": {name: _records(frame) for name, frame in scenarios.items()},
        "immune_coupling": {name: _immune_records(frame) for name, frame in coupled.items()},
        "anti_pd1_dose_days": [0, 21, 42, 63],
        "parameter_manifest": parameter_manifest(),
        "coupling_parameter_manifest": asdict(LiverTDLNCouplingParameters()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("docs/data/liver_tumor_dashboard.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(build_dashboard_data(), ensure_ascii=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()

