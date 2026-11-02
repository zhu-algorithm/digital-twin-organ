import pandas as pd
import pytest

from tdln_pd1.sensitivity import global_sensitivity


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
