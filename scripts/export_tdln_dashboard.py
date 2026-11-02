"""Export deterministic, versioned trajectories for the TDLN web dashboard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

# Keep the exporter runnable from a fresh checkout before editable installation.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tdln_pd1.simulator import Dose, parameter_manifest, simulate
from tdln_pd1.spatial import ZONES
from tdln_pd1.spatial_dynamics import SPECIES, simulate_spatial


SCHEMA_VERSION = 1


def _sample(frame, days: np.ndarray, columns: tuple[str, ...]) -> list[dict[str, float]]:
    return [
        {"day": float(day), **{name: round(float(np.interp(day, frame.day, frame[name])), 6)
                                for name in columns}}
        for day in days
    ]


def build_dashboard_data() -> dict:
    days = np.arange(0.0, 42.01, 2.0)
    columns = ("tumour_cells", "tpex_tdln", "effector_tdln", "effector_tumour", "pd1_occupancy")
    control = simulate(days=42)
    treated = simulate(days=42, doses=tuple(Dose(day, 1.0) for day in (0, 14, 28)))

    spatial = simulate_spatial(days=14, step=0.25)
    spatial_days = np.arange(0.0, 14.01, 1.0)
    spatial_records = []
    for day in spatial_days:
        values = {}
        for species in SPECIES:
            species_frame = spatial[spatial.species == species]
            values[species] = {
                zone.name: round(float(np.interp(
                    day,
                    species_frame[species_frame.zone == zone.name].day,
                    species_frame[species_frame.zone == zone.name].amount,
                )), 6)
                for zone in ZONES
            }
        spatial_records.append({"day": float(day), "values": values})

    return {
        "schema_version": SCHEMA_VERSION,
        "model_version": "tdln-pd1/0.1.0",
        "research_only": True,
        "time_unit": "day",
        "zones": [{"id": zone.name, "role": zone.role} for zone in ZONES],
        "species": list(SPECIES),
        "scenarios": {
            "control": _sample(control, days, columns),
            "anti_pd1": _sample(treated, days, columns),
        },
        "spatial": spatial_records,
        "dose_schedule_days": [0, 14, 28],
        "parameter_manifest": parameter_manifest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("docs/data/tdln_dashboard.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(build_dashboard_data(), ensure_ascii=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
