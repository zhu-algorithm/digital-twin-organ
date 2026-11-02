import pandas as pd
import pytest

from tdln_pd1.sensitivity import global_sensitivity, global_spatial_sensitivity
from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.spatial_dynamics import SpatialParameters


def test_sensitivity_is_reproducible_and_detects_dominant_growth_effect():
    kwargs = dict(names=("tumour_growth", "killing_rate"), days=10, n=32, seed=9, fold_range=3)
    samples_a, ranking_a = global_sensitivity(**kwargs)
    samples_b, ranking_b = global_sensitivity(**kwargs)
    pd.testing.assert_frame_equal(samples_a, samples_b)
    pd.testing.assert_frame_equal(ranking_a, ranking_b)
    tumour = ranking_a[ranking_a.outcome == "tumour_cells"].set_index("parameter")
    assert tumour.loc["tumour_growth", "spearman_rho"] > 0
    assert tumour.loc["tumour_growth", "abs_rho"] > tumour.loc["killing_rate", "abs_rho"]


def test_sensitivity_contract_rejects_invalid_configuration():
    with pytest.raises(ValueError):
        global_sensitivity(names=("not_a_parameter",), n=8)
    with pytest.raises(ValueError):
        global_sensitivity(names=("tumour_growth",), n=4)
    with pytest.raises(ValueError, match="before the simulation endpoint"):
        global_sensitivity(
            names=("tumour_growth",), days=10, n=8,
            doses=(Dose(0, 1), Dose(14, 1)),
        )


def test_default_sensitivity_schedule_adapts_without_hiding_explicit_events():
    samples, _ = global_sensitivity(
        names=("tumour_growth",), days=10, n=8, seed=3,
    )
    assert len(samples) == 8
    assert samples.peak_pd1_occupancy.gt(0).all()


def test_spatial_sensitivity_is_reproducible_and_preserves_paired_pd1_contrast():
    kwargs = dict(
        names=("antigen_half_saturation", "checkpoint_growth_gain"),
        days=4, checkpoint_occupancy=0.8, n=12, seed=11, fold_range=2,
    )
    samples_a, ranking_a = global_spatial_sensitivity(**kwargs)
    samples_b, ranking_b = global_spatial_sensitivity(**kwargs)
    pd.testing.assert_frame_equal(samples_a, samples_b)
    pd.testing.assert_frame_equal(ranking_a, ranking_b)
    assert samples_a.pd1_efferent_effector_delta.gt(0).all()
    assert samples_a.pd1_paracortex_tpex_delta.gt(0).all()
    assert set(ranking_a.outcome) == {
        "control_efferent_effector", "occupied_efferent_effector",
        "pd1_efferent_effector_delta", "control_paracortex_tpex",
        "occupied_paracortex_tpex", "pd1_paracortex_tpex_delta",
    }


def test_spatial_capacity_is_available_to_uncertainty_screening():
    samples, ranking = global_spatial_sensitivity(
        names=("tpex_carrying_capacity",), days=4, n=8, seed=2,
    )
    assert samples.tpex_carrying_capacity.nunique() == 8
    assert set(ranking.parameter) == {"tpex_carrying_capacity"}


def test_spatial_sensitivity_rejects_invalid_log_domain_and_occupancy():
    with pytest.raises(ValueError, match="strictly positive"):
        global_spatial_sensitivity(
            parameters=SpatialParameters(checkpoint_growth_gain=0),
            names=("checkpoint_growth_gain",), n=8,
        )
    with pytest.raises(ValueError, match="checkpoint_occupancy"):
        global_spatial_sensitivity(names=("tpex_growth",), n=8, checkpoint_occupancy=0)


def test_spatial_sensitivity_accepts_time_varying_systemic_occupancy():
    systemic = simulate(days=4, doses=(Dose(0, 1.0),))
    samples, ranking = global_spatial_sensitivity(
        names=("checkpoint_growth_gain",), days=4, n=8,
        systemic_trajectory=systemic,
    )
    assert samples.pd1_efferent_effector_delta.gt(0).all()
    assert set(ranking.checkpoint_mode) == {"systemic_trajectory"}
    assert ranking.checkpoint_occupancy.isna().all()
    control = ranking[ranking.outcome == "control_efferent_effector"].iloc[0]
    assert bool(control.estimable) is False
    assert control.spearman_rho == 0
