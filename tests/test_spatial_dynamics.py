import numpy as np

from tdln_pd1.spatial_dynamics import SPECIES, simulate_spatial, spatial_summary


def test_transport_only_conserves_each_species():
    trajectory = simulate_spatial(days=5, reactions=False)
    totals = trajectory.groupby(["day", "species"]).amount.sum().unstack()
    for species in SPECIES:
        assert np.allclose(totals[species], totals[species].iloc[0], rtol=1e-7, atol=1e-7)


def test_antigen_reaches_sinus_and_dc_reaches_paracortex():
    trajectory = simulate_spatial(days=4)
    summary = spatial_summary(trajectory)
    assert summary[("antigen", "subcapsular_sinus")].max() > 0
    assert summary[("mature_dc", "paracortex")].max() > 0


def test_effectors_exit_through_efferent_lymph_and_states_are_nonnegative():
    trajectory = simulate_spatial(days=6)
    efferent = trajectory[(trajectory.species == "effector") & (trajectory.zone == "efferent_lymph")]
    assert efferent.amount.iloc[-1] > 0
    assert trajectory.amount.ge(0).all()
