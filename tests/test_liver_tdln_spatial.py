import numpy as np
import pytest

from tdln_pd1.liver_tdln_spatial import (LiverTDLNSpatialBridge,
                                         liver_tdln_spatial_counterfactual,
                                         simulate_liver_tdln_spatial)
from tdln_pd1.simulator import Dose


def _series(frame, species, zone):
    return frame[frame.species.eq(species) & frame.zone.eq(zone)].amount


def test_liver_antigen_enters_sinus_and_drives_paracortical_dc():
    systemic, spatial = simulate_liver_tdln_spatial(days=10, step=0.5)
    assert systemic.tumour_antigen.max() > 0
    assert spatial.afferent_antigen_source.max() > 0
    assert _series(spatial, "antigen", "subcapsular_sinus").max() > 0
    assert _series(spatial, "mature_dc", "paracortex").max() > 0


def test_zero_drainage_removes_liver_driven_antigen_and_dc_signal():
    _, connected = simulate_liver_tdln_spatial(days=8, step=0.5)
    _, disconnected = simulate_liver_tdln_spatial(
        days=8, step=0.5, bridge=LiverTDLNSpatialBridge(antigen_drainage_gain=0))
    assert _series(disconnected, "antigen", "subcapsular_sinus").max() == pytest.approx(0)
    assert _series(disconnected, "mature_dc", "paracortex").max() == pytest.approx(0)
    assert _series(connected, "mature_dc", "paracortex").max() > 0


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
