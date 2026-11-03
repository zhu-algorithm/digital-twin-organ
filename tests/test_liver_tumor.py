import numpy as np
import pytest

from tdln_pd1.liver_tumor import LiverTumourInitialState, simulate_liver_tumour


def test_liver_tumour_simulation_is_finite_nonnegative_and_grows():
    frame = simulate_liver_tumour(days=60)
    assert np.isfinite(frame.to_numpy()).all()
    assert (frame.drop(columns="day") >= 0).all().all()
    assert frame.total_tumour.iloc[-1] > frame.total_tumour.iloc[0]
    assert frame.normal_cells.iloc[-1] < frame.normal_cells.iloc[0]


def test_reduced_blood_flow_reduces_oxygen_and_increases_necrosis():
    normal = simulate_liver_tumour(days=90)
    low = simulate_liver_tumour(days=90, flow_multiplier=0.55)
    assert low.oxygen_core.iloc[-1] < normal.oxygen_core.iloc[-1]
    assert low.necrotic_tumour.iloc[-1] > normal.necrotic_tumour.iloc[-1]


def test_zero_tumour_does_not_spontaneously_create_tumour():
    state = LiverTumourInitialState(tumour_rim=0, tumour_core=0, necrotic_tumour=0)
    frame = simulate_liver_tumour(days=30, initial=state)
    assert frame.total_tumour.max() == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("kwargs", ({"days": 0}, {"step": -1}, {"flow_multiplier": 0},
                                    {"angiogenesis_inhibition": 1.2}))
def test_invalid_simulation_contracts_fail(kwargs):
    with pytest.raises(ValueError):
        simulate_liver_tumour(**kwargs)
