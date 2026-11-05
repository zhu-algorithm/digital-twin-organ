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
from tdln_pd1.liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                         liver_tdln_spatial_counterfactual)
from tdln_pd1.simulator import Dose
from tdln_pd1.spatial import ZONES

SCHEMA_VERSION = 3
COLUMNS = ("normal_cells", "tumour_rim", "tumour_core", "necrotic_tumour",
           "functional_vessels", "oxygen_normal", "oxygen_rim", "oxygen_core",
           "perfusion", "portal_flow", "arterial_flow", "arterial_buffer_flow",
           "oxygenated_inflow",
           "arterial_flow_fraction", "tumour_fraction")
IMMUNE_COLUMNS = ("tumour_antigen", "mature_dc_tdln", "tpex_tdln",
                  "effector_tdln", "effector_tumour", "pd1_occupancy",
                  "tumour_cells", "liver_tumour_fraction",
                  "liver_normal_cells", "liver_tumour_rim",
                  "liver_tumour_core", "liver_necrotic_tumour",
                  "liver_functional_vessels", "liver_oxygen_normal",
                  "liver_oxygen_rim", "liver_oxygen_core",
                  "liver_perfusion", "liver_portal_flow",
                  "liver_arterial_flow", "liver_arterial_buffer_flow",
                  "liver_oxygenated_inflow")


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
            elif name.startswith("liver_"):
                value = round(value, 6)
            else:
                value = round(value, 4)
            row[name] = value
        records.append(row)
    return records


def _spatial_records(frame):
    display = lambda value, digits: float(format(float(value), f".{digits}g"))
    records = {}
    for arm in ("control", "anti_pd1"):
        arm_frame = frame[frame.arm.eq(arm)]
        arm_records = []
        for day in np.arange(0.0, 28.01, 2.0):
            values = {}
            for species in ("antigen", "mature_dc", "tpex", "effector"):
                values[species] = {}
                for zone in (item.name for item in ZONES):
                    series = arm_frame[arm_frame.species.eq(species) & arm_frame.zone.eq(zone)]
                    value = float(np.interp(day, series.day, series.amount))
                    values[species][zone] = display(value, 3 if species == "antigen" else 2)
            source_rows = arm_frame[arm_frame.day.eq(day)]
            antigen_source = float(source_rows.afferent_antigen_source.iloc[0])
            dc_source = float(source_rows.afferent_dc_source.iloc[0])
            arm_records.append({"day": float(day),
                                "afferent_antigen_source": display(antigen_source, 3),
                                "afferent_dc_source": display(dc_source, 3),
                                "values": values})
        records[arm] = arm_records
    return records


def build_dashboard_data() -> dict:
    scenarios = {
        "untreated": simulate_liver_tumour(),
        "antiangiogenic": simulate_liver_tumour(angiogenesis_inhibition=0.75),
        "low_perfusion": simulate_liver_tumour(flow_multiplier=0.65),
        "portal_restriction": simulate_liver_tumour(portal_flow_multiplier=0.55),
        "arterial_restriction": simulate_liver_tumour(arterial_flow_multiplier=0.05),
    }
    coupled = {
        "control": simulate_liver_tdln(days=84),
        "anti_pd1": simulate_liver_tdln(
            days=84, doses=tuple(Dose(day, 1.0) for day in (0, 21, 42, 63))),
    }
    spatial = liver_tdln_spatial_counterfactual(
        doses=(Dose(0, 1.0), Dose(7, 1.0), Dose(14, 1.0), Dose(21, 1.0)),
        days=28, step=0.25,
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "model_version": "liver-tumour/0.3.0",
        "research_only": True,
        "time_unit": "day",
        "state_units": "normalised fraction",
        "scenarios": {name: _records(frame) for name, frame in scenarios.items()},
        "immune_coupling": {name: _immune_records(frame) for name, frame in coupled.items()},
        "tdln_zones": [{"id": zone.name, "role": zone.role} for zone in ZONES],
        "liver_driven_tdln_spatial": _spatial_records(spatial),
        "anti_pd1_dose_days": [0, 21, 42, 63],
        "parameter_manifest": parameter_manifest(),
        "coupling_parameter_manifest": asdict(LiverTDLNCouplingParameters()),
        "spatial_bridge_parameter_manifest": asdict(LiverTDLNSpatialBridge()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("docs/data/liver_tumor_dashboard.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(build_dashboard_data(), ensure_ascii=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()

