import pandas as pd
import pytest

from tdln_pd1.spatial_data import aggregate_spatial_spots
from tdln_pd1.spatial_dynamics import simulate_spatial


def _spots():
    return pd.DataFrame({
        "spot_id": ["s1", "s2", "s3"], "patient_id": ["p1"] * 3,
        "zone": ["paracortex", "paracortex", "afferent_lymph"],
        "area": [2.0, 1.0, 0.5], "antigen": [0.0, 0.0, 2.0],
        "mature_dc": [1.0, 2.0, 4.0], "tpex": [3.0, 5.0, 0.0],
        "effector": [7.0, 11.0, 0.0],
    })


def test_spatial_spots_aggregate_to_auditable_zone_states():
    states = aggregate_spatial_spots(_spots(), "area")
    amount = states.set_index(["species", "zone"]).amount
    assert amount[("antigen", "afferent_lymph")] == 1.0
    assert amount[("tpex", "paracortex")] == 11.0
    assert amount[("effector", "paracortex")] == 25.0
    trajectory = simulate_spatial(days=1, step=1, initial=states)
    day_zero = trajectory[trajectory.day.eq(0)].set_index(["species", "zone"]).amount
    assert day_zero[("tpex", "paracortex")] == 11.0
    assert day_zero[("effector", "paracortex")] == 25.0


def test_spatial_ingestion_fails_closed_on_anatomy_and_abundance_errors():
    invalid = _spots()
    invalid.loc[0, "zone"] = "unknown"
    with pytest.raises(ValueError, match="unknown TDLN zones"):
        aggregate_spatial_spots(invalid)
    negative = _spots()
    negative.loc[0, "tpex"] = -1
    with pytest.raises(ValueError, match="finite and non-negative"):
        aggregate_spatial_spots(negative)
    duplicate = _spots()
    duplicate.loc[1, "spot_id"] = "s1"
    with pytest.raises(ValueError, match="unique within patient"):
        aggregate_spatial_spots(duplicate)
