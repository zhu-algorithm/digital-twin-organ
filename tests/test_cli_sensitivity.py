from dataclasses import asdict
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

from tdln_pd1.cli import _load_patient_spatial_initial, main
from tdln_pd1.simulator import TwinParameters
from tdln_pd1.spatial_dynamics import SPECIES, ZONE_NAMES, SpatialParameters
from tdln_pd1.spatial_feature_audit import write_spatial_feature_audit


def _complete_states(path: Path) -> None:
    rows = []
    for patient_id, scale in (("p1", 1.0), ("p2", 2.0)):
        for species in SPECIES:
            for zone in ZONE_NAMES:
                rows.append({
                    "patient_id": patient_id, "species": species, "zone": zone,
                    "amount": scale if (species, zone) == ("antigen", "afferent_lymph") else 0.0,
                })
    pd.DataFrame(rows).to_csv(path, index=False)


def test_patient_spatial_initial_selection_is_hashed_and_fail_closed(tmp_path):
    path = tmp_path / "states.csv"
    _complete_states(path)
    selected, provenance = _load_patient_spatial_initial(str(path), "p2")
    assert set(selected.patient_id) == {"p2"}
    assert provenance["patient_id"] == "p2"
    assert provenance["rows"] == len(SPECIES) * len(ZONE_NAMES)
    assert len(provenance["sha256"]) == 64
    with pytest.raises(ValueError, match="patient-id is required"):
        _load_patient_spatial_initial(str(path), None)
    with pytest.raises(ValueError, match="requires --initial-zones"):
        _load_patient_spatial_initial(None, "p1")


def test_patient_spatial_initial_binds_upstream_spot_derivation(tmp_path):
    spots = tmp_path / "spots.csv"
    spots.write_text("patient_id,zone\np1,paracortex\n")
    states = tmp_path / "states.csv"
    _complete_states(states)
    audit = write_spatial_feature_audit(spots, states, None)
    _, provenance = _load_patient_spatial_initial(str(states), "p1")
    assert provenance["derivation_audit"]["file"] == audit.name
    assert provenance["derivation_audit"]["stage"] == (
        "label_free_spatial_state_featurization"
    )
    audit.write_text("replaced\n")
    with pytest.raises(ValueError, match="invalid derivation audit"):
        _load_patient_spatial_initial(str(states), "p1")


def test_spatial_sensitivity_cli_binds_fitted_parameters_and_patient_state(
    tmp_path, monkeypatch,
):
    parameter_path = tmp_path / "spatial_parameters.json"
    parameter_path.write_text(json.dumps(asdict(SpatialParameters())))
    state_path = tmp_path / "states.csv"
    _complete_states(state_path)
    output = tmp_path / "screen"
    monkeypatch.setattr(sys, "argv", [
        "tdln-pd1", "spatial-sensitivity",
        "--spatial-parameters", str(parameter_path),
        "--initial-zones", str(state_path), "--patient-id", "p1",
        "--parameters", "antigen_half_saturation,checkpoint_growth_gain",
        "--days", "1", "--samples", "8", "--output", str(output),
    ])
    main()
    provenance = json.loads((output / "spatial_sensitivity_provenance.json").read_text())
    assert provenance["inputs"]["parameter_source"]["parameter_class"] == "SpatialParameters"
    assert provenance["inputs"]["initial_state_source"]["patient_id"] == "p1"
    assert provenance["settings"]["sampled_parameters"] == [
        "antigen_half_saturation", "checkpoint_growth_gain",
    ]
    assert set(provenance["outputs"]) == {
        "spatial_sensitivity_samples.csv", "spatial_sensitivity_ranking.csv",
    }
    assert (output / "spatial_sensitivity_samples.csv").is_file()
    assert (output / "spatial_sensitivity_ranking.csv").is_file()


def test_spatial_sensitivity_cli_uses_dose_driven_occupancy(tmp_path, monkeypatch):
    spatial_path = tmp_path / "spatial.json"
    spatial_path.write_text(json.dumps(asdict(SpatialParameters())))
    systemic_path = tmp_path / "systemic.json"
    systemic_path.write_text(json.dumps(asdict(TwinParameters())))
    output = tmp_path / "screen"
    monkeypatch.setattr(sys, "argv", [
        "tdln-pd1", "spatial-sensitivity",
        "--spatial-parameters", str(spatial_path),
        "--systemic-parameters", str(systemic_path),
        "--dose-days", "0", "--dose", "1",
        "--parameters", "checkpoint_growth_gain",
        "--days", "2", "--samples", "8", "--output", str(output),
    ])
    main()
    audit = json.loads((output / "spatial_sensitivity_provenance.json").read_text())
    assert audit["settings"]["checkpoint_mode"] == "dose_driven_systemic_trajectory"
    assert audit["settings"]["peak_pd1_occupancy"] > 0
    assert audit["inputs"]["systemic_parameter_source"]["parameter_class"] == "TwinParameters"


def test_systemic_sensitivity_cli_adapts_omitted_schedule_to_short_window(
    tmp_path, monkeypatch,
):
    output = tmp_path / "systemic_screen"
    monkeypatch.setattr(sys, "argv", [
        "tdln-pd1", "sensitivity", "--parameters", "tumour_growth",
        "--days", "3", "--samples", "8", "--output", str(output),
    ])
    main()
    audit = json.loads((output / "sensitivity_provenance.json").read_text())
    assert audit["settings"]["dose_days"] == [0.0]
    assert audit["settings"]["dose_amounts"] == [1.0]


def test_spatial_uncertainty_cli_requires_a_real_treated_schedule(tmp_path, monkeypatch):
    parameters = tmp_path / "spatial.json"
    parameters.write_text(json.dumps(asdict(SpatialParameters())))
    diagnostics = tmp_path / "diagnostics.json"
    diagnostics.write_text(json.dumps({
        "locally_identifiable": True,
        "local_log_parameter_covariance": {
            "parameter_order": ["antigen_decay"],
            "log_scale_covariance": [[0.01]],
        },
    }))
    monkeypatch.setattr(sys, "argv", [
        "tdln-pd1", "propagate-spatial-uncertainty",
        "--parameters", str(parameters), "--diagnostics", str(diagnostics),
        "--dose-days", "", "--samples", "2", "--output", str(tmp_path / "out"),
    ])
    with pytest.raises(ValueError, match="explicit positive dose schedule"):
        main()
