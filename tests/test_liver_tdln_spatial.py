import numpy as np
import pytest

from tdln_pd1.liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                         liver_tdln_spatial_counterfactual,
                                         simulate_liver_tdln_spatial)
from tdln_pd1.simulator import Dose, TwinParameters
from tdln_pd1.liver_tdln_coupling import simulate_liver_tdln


def _series(frame, species, zone):
    return frame[frame.species.eq(species) & frame.zone.eq(zone)].amount


def test_liver_antigen_and_migratory_dc_enter_the_spatial_tdln():
    systemic, spatial = simulate_liver_tdln_spatial(days=10, step=0.5)
    assert systemic.tumour_antigen.max() > 0
    assert spatial.afferent_antigen_source.max() > 0
    assert spatial.afferent_dc_source.max() > 0
    assert _series(spatial, "antigen", "subcapsular_sinus").max() > 0
    assert _series(spatial, "mature_dc", "afferent_lymph").max() > 0
    assert _series(spatial, "mature_dc", "paracortex").max() > 0


def test_zero_drainage_removes_liver_driven_antigen_and_dc_signal():
    _, connected = simulate_liver_tdln_spatial(days=8, step=0.5)
    _, disconnected = simulate_liver_tdln_spatial(
        days=8, step=0.5, bridge=LiverTDLNSpatialBridge(
            antigen_drainage_gain=0, dc_migration_gain=0))
    assert _series(disconnected, "antigen", "subcapsular_sinus").max() == pytest.approx(0)
    assert _series(disconnected, "mature_dc", "paracortex").max() == pytest.approx(0)
    assert _series(connected, "mature_dc", "paracortex").max() > 0


def test_soluble_antigen_and_migratory_dc_routes_are_independently_executable():
    _, dc_only = simulate_liver_tdln_spatial(
        days=8, step=0.5, bridge=LiverTDLNSpatialBridge(
            antigen_drainage_gain=0, dc_migration_gain=100))
    _, antigen_only = simulate_liver_tdln_spatial(
        days=8, step=0.5, bridge=LiverTDLNSpatialBridge(
            antigen_drainage_gain=0.45, dc_migration_gain=0))
    assert _series(dc_only, "antigen", "subcapsular_sinus").max() == pytest.approx(0)
    assert _series(dc_only, "mature_dc", "paracortex").max() > 0
    assert _series(antigen_only, "antigen", "subcapsular_sinus").max() > 0
    assert _series(antigen_only, "mature_dc", "paracortex").max() > 0


def test_dc_source_uses_explicit_systemic_migration_flux_and_perfusion():
    immune = TwinParameters(dc_migration=0.2)
    systemic, spatial = simulate_liver_tdln_spatial(
        days=3, step=0.5, immune_parameters=immune)
    source_rows = spatial[spatial.species.eq("antigen") &
                          spatial.zone.eq("afferent_lymph")]
    expected_flux = immune.dc_migration * np.interp(
        source_rows.day, systemic.day, systemic.tumour_antigen)
    assert source_rows.systemic_dc_migration_flux.to_numpy() == pytest.approx(expected_flux)
    expected_source = (LiverTDLNSpatialBridge().dc_migration_gain
                       * expected_flux
                       * source_rows.liver_perfusion_delivery_factor)
    assert source_rows.afferent_dc_source.to_numpy() == pytest.approx(expected_source)


def test_lower_portal_flow_reduces_both_tdln_delivery_routes():
    _, nominal = simulate_liver_tdln_spatial(days=3, step=0.5)
    _, low_flow = simulate_liver_tdln_spatial(
        days=3, step=0.5, portal_flow_multiplier=0.5)
    nominal_source = nominal[nominal.species.eq("antigen") &
                             nominal.zone.eq("afferent_lymph")].iloc[0]
    low_source = low_flow[low_flow.species.eq("antigen") &
                          low_flow.zone.eq("afferent_lymph")].iloc[0]
    assert low_source.liver_perfusion_delivery_factor < nominal_source.liver_perfusion_delivery_factor
    assert low_source.afferent_antigen_source < nominal_source.afferent_antigen_source
    assert low_source.afferent_dc_source < nominal_source.afferent_dc_source


def test_efferent_effectors_feed_back_to_reduce_viable_liver_tumour():
    connected_systemic, _ = simulate_liver_tdln_spatial(days=28, step=0.5)
    disconnected_systemic, _ = simulate_liver_tdln_spatial(
        days=28, step=0.5,
        bridge=LiverTDLNSpatialBridge(efferent_effector_killing_gain=0.0))
    assert connected_systemic.spatial_efferent_effector.iloc[-1] > 0
    assert connected_systemic.spatial_effector_killing_flux.max() > 0
    assert connected_systemic.systemic_effector_killing_flux.max() > 0
    assert connected_systemic.total_effector_killing_flux.to_numpy() == pytest.approx(
        connected_systemic.systemic_effector_killing_flux
        + connected_systemic.spatial_effector_killing_flux)
    assert connected_systemic.spatial_effector_killing_fraction.between(0, 1).all()
    assert connected_systemic.tumour_cells.iloc[-1] < disconnected_systemic.tumour_cells.iloc[-1]


def test_zero_spatial_feedback_recovers_systemic_closed_loop():
    systemic, _ = simulate_liver_tdln_spatial(
        days=8, step=0.5,
        bridge=LiverTDLNSpatialBridge(efferent_effector_killing_gain=0.0))
    baseline = simulate_liver_tdln(days=8, step=0.5)
    assert systemic.day.to_numpy() == pytest.approx(baseline.day)
    for column in ("liver_tumour_rim", "liver_tumour_core", "tumour_antigen",
                   "effector_tumour", "pd1_occupancy"):
        assert systemic[column].to_numpy() == pytest.approx(
            baseline[column], rel=2e-6, abs=1e-8)


def test_anti_pd1_increases_efferent_effector_output_in_matched_closed_loop():
    paired = liver_tdln_spatial_counterfactual(
        doses=(Dose(0, 1.0), Dose(7, 1.0)), days=14, step=0.5)
    endpoint = paired[paired.species.eq("effector") & paired.zone.eq("efferent_lymph")]
    values = endpoint.groupby("arm").amount.last()
    assert values["anti_pd1"] > values["control"]
    assert paired[paired.arm.eq("anti_pd1")].pd1_occupancy.max() > 0


def test_spatial_bridge_fails_closed_for_invalid_gain():
    with pytest.raises(ValueError, match="bridge parameters"):
        simulate_liver_tdln_spatial(
            days=2, bridge=LiverTDLNSpatialBridge(antigen_drainage_gain=np.nan))
