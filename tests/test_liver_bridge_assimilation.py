import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tdln_pd1.liver_bridge_assimilation import (calibrate_liver_bridge,
                                                 save_liver_bridge_calibration,
                                                 validate_bridge_observations)
from tdln_pd1.calibration_audit import (verify_calibration_audit,
                                        write_calibration_audit)
from tdln_pd1.provenance import sha256_file
from tdln_pd1.liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                         simulate_liver_tdln_spatial)
from tdln_pd1.liver_tumor import LiverTumourParameters


ROOT = Path(__file__).parents[1]


def _synthetic_observations(gain=0.8):
    _, spatial = simulate_liver_tdln_spatial(
        days=8, step=0.5, bridge=LiverTDLNSpatialBridge(antigen_drainage_gain=gain))
    rows = []
    for day in (2.0, 4.0, 6.0, 8.0):
        for species, zone in (("antigen", "afferent_lymph"),
                              ("mature_dc", "paracortex")):
            series = spatial[spatial.species.eq(species) & spatial.zone.eq(zone)]
            value = float(np.interp(day, series.day, series.amount))
            rows.append({"day": day, "target": f"{species}@{zone}",
                         "value": value, "sd": max(value * 0.01, 1e-4),
                         "source": "synthetic_unit_test"})
    return pd.DataFrame(rows)


def test_calibration_recovers_spatial_drainage_gain():
    fitted_liver, fitted_coupling, fitted_bridge, diagnostics, residuals = calibrate_liver_bridge(
        _synthetic_observations(), base_bridge=LiverTDLNSpatialBridge(0.2), max_nfev=60)
    assert fitted_bridge.antigen_drainage_gain == pytest.approx(0.8, rel=2e-3)
    assert diagnostics["success"] is True
    assert diagnostics["optimizer_initialization"] == (
        "weighted_least_squares_then_soft_l1_refinement"
    )
    assert diagnostics["weighted_rmse"] < 1e-3
    assert len(residuals) == 8
    assert fitted_coupling.immune_killing_gain > 0
    assert fitted_liver.arterial_buffer_fraction > 0


def test_calibration_recovers_independent_migratory_dc_gain():
    truth = LiverTDLNSpatialBridge(antigen_drainage_gain=0, dc_migration_gain=180)
    _, spatial = simulate_liver_tdln_spatial(days=8, step=0.5, bridge=truth)
    rows = []
    for day in (2.0, 4.0, 6.0, 8.0):
        series = spatial[
            spatial.species.eq("mature_dc")
            & spatial.zone.eq("afferent_lymph")]
        value = float(np.interp(day, series.day, series.amount))
        rows.append({"day": day, "target": "mature_dc@afferent_lymph",
                     "value": value, "sd": max(value * 0.01, 1e-4)})
    _, _, fitted, diagnostics, _ = calibrate_liver_bridge(
        pd.DataFrame(rows), fit_parameters=("bridge.dc_migration_gain",),
        base_bridge=LiverTDLNSpatialBridge(
            antigen_drainage_gain=0, dc_migration_gain=25), max_nfev=60)
    assert fitted.dc_migration_gain == pytest.approx(180, rel=2e-3)
    assert diagnostics["success"] is True


def test_bundled_demo_jointly_recovers_both_afferent_route_gains():
    observations = pd.read_csv(ROOT / "examples/liver_tdln_bridge_observations.csv")
    _, _, fitted, diagnostics, residuals = calibrate_liver_bridge(
        observations,
        fit_parameters=("bridge.antigen_drainage_gain", "bridge.dc_migration_gain"),
        max_nfev=60,
    )
    assert fitted.antigen_drainage_gain == pytest.approx(0.8, rel=2e-3)
    assert fitted.dc_migration_gain == pytest.approx(180, rel=2e-3)
    assert diagnostics["success"] is True
    assert diagnostics["locally_identifiable"] is True
    assert diagnostics["identifiability_warning"] is False
    assert diagnostics["jacobian_rank"] == 2
    assert len(residuals) == 12


def test_calibration_recovers_hepatic_arterial_buffer_fraction():
    truth = LiverTumourParameters(arterial_buffer_fraction=0.40)
    systemic, _ = simulate_liver_tdln_spatial(
        days=5, step=0.25, liver_parameters=truth, portal_flow_multiplier=0.4)
    observations = pd.DataFrame([
        {"day": day, "target": "liver_arterial_flow",
         "value": float(np.interp(day, systemic.day, systemic.liver_arterial_flow)),
         "sd": 0.001, "source": "synthetic_unit_test"}
        for day in (1, 2, 3, 4, 5)
    ])
    liver, _, _, diagnostics, _ = calibrate_liver_bridge(
        observations, fit_parameters=("liver.arterial_buffer_fraction",),
        base_liver=LiverTumourParameters(arterial_buffer_fraction=0.10),
        portal_flow_multiplier=0.4, max_nfev=60)
    assert liver.arterial_buffer_fraction == pytest.approx(0.40, rel=2e-3)
    assert diagnostics["success"] is True
    assert diagnostics["weighted_rmse"] < 1e-3


def test_buffer_calibration_requires_an_isolated_portal_perturbation():
    observations = _synthetic_observations()
    with pytest.raises(ValueError, match="isolated portal perturbation"):
        calibrate_liver_bridge(
            observations,
            fit_parameters=("liver.arterial_buffer_fraction",))
    with pytest.raises(ValueError, match="isolated portal perturbation"):
        calibrate_liver_bridge(
            observations,
            fit_parameters=("liver.arterial_buffer_fraction",),
            portal_flow_multiplier=0.5, arterial_flow_multiplier=1.2)


def test_observation_contract_rejects_unknown_duplicates_and_bad_uncertainty():
    observations = _synthetic_observations()
    with pytest.raises(ValueError, match="unknown"):
        validate_bridge_observations(observations.assign(target="unknown"))
    with pytest.raises(ValueError, match="unique"):
        validate_bridge_observations(pd.concat([observations, observations.iloc[[0]]]))
    with pytest.raises(ValueError, match="strictly positive"):
        validate_bridge_observations(observations.assign(sd=0))
    with pytest.raises(ValueError, match="fit parameters"):
        calibrate_liver_bridge(observations, fit_parameters=("bridge.unknown",))


def test_saved_calibration_is_explicitly_research_only(tmp_path):
    liver, coupling, bridge, diagnostics, residuals = calibrate_liver_bridge(
        _synthetic_observations(), base_bridge=LiverTDLNSpatialBridge(0.4), max_nfev=30)
    save_liver_bridge_calibration(tmp_path, liver, coupling, bridge, diagnostics, residuals)
    payload = json.loads((tmp_path / "fitted_liver_tdln_bridge.json").read_text())
    assert payload["research_only"] is True
    assert "liver_tumour" in payload
    assert (tmp_path / "liver_tdln_bridge_diagnostics.json").is_file()
    assert (tmp_path / "liver_tdln_bridge_residuals.csv").is_file()


def test_bridge_calibration_outputs_are_cryptographically_bound(tmp_path):
    observations_path = tmp_path / "observations.csv"
    observations = _synthetic_observations()
    observations.to_csv(observations_path, index=False)
    liver, coupling, bridge, diagnostics, residuals = calibrate_liver_bridge(
        observations, base_bridge=LiverTDLNSpatialBridge(0.4), max_nfev=30)
    save_liver_bridge_calibration(tmp_path, liver, coupling, bridge, diagnostics, residuals)
    audit = write_calibration_audit(
        tmp_path, "liver_tdln_bridge_provenance.json",
        "liver_tdln_bridge_calibration",
        inputs={"observations": {"file": observations_path.name,
                                  "sha256": sha256_file(observations_path)}},
        settings={"fit_parameters": ["bridge.antigen_drainage_gain"]},
        generated_filenames=("fitted_liver_tdln_bridge.json",
                             "liver_tdln_bridge_diagnostics.json",
                             "liver_tdln_bridge_residuals.csv"),
    )
    verified = verify_calibration_audit(
        audit, input_paths={"observations": observations_path})
    assert verified["valid"] is True
    residuals_path = tmp_path / "liver_tdln_bridge_residuals.csv"
    residuals_path.write_text(residuals_path.read_text() + "\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_calibration_audit(audit, input_paths={"observations": observations_path})
