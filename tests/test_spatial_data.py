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
    assert len(states) == 24
    observed = states[["zone", "zone_observed"]].drop_duplicates().set_index("zone")
    assert observed.zone_observed.sum() == 2
    assert not observed.loc["medulla", "zone_observed"]
    support = states[["zone", "sampling_units_observed", "sampling_units_total",
                      "sampling_unit_coverage"]].drop_duplicates().set_index("zone")
    assert support.loc["paracortex"].to_dict() == {
        "sampling_units_observed": 1.0, "sampling_units_total": 1.0,
        "sampling_unit_coverage": 1.0,
    }
    assert support.loc["medulla", "sampling_unit_coverage"] == 0.0
    amount = states.set_index(["species", "zone"]).amount
    assert amount[("antigen", "afferent_lymph")] == 1.0
    assert amount[("tpex", "paracortex")] == 11.0
    assert amount[("effector", "paracortex")] == 25.0
    trajectory = simulate_spatial(days=1, step=1, initial=states)
    day_zero = trajectory[trajectory.day.eq(0)].set_index(["species", "zone"]).amount
    assert day_zero[("tpex", "paracortex")] == 11.0
    assert day_zero[("effector", "paracortex")] == 25.0


def test_spatial_simulator_rejects_an_implicit_partial_state_grid():
    partial = pd.DataFrame({
        "species": ["tpex", "effector"], "zone": ["paracortex"] * 2,
        "amount": [10.0, 5.0],
    })
    with pytest.raises(ValueError, match="complete species-zone grid"):
        simulate_spatial(days=1, initial=partial)


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


def test_spatial_temporal_input_requires_explicit_baseline_and_excludes_post_treatment():
    baseline = _spots().assign(timepoint="baseline")
    post = _spots().assign(
        spot_id=lambda x: x.spot_id + "_post", timepoint="post", tpex=1000000.0,
    )
    longitudinal = pd.concat([baseline, post], ignore_index=True)
    with pytest.raises(ValueError, match="temporal columns detected"):
        aggregate_spatial_spots(longitudinal, "area")
    states = aggregate_spatial_spots(
        longitudinal, "area", "timepoint", ("baseline",)
    )
    amount = states.set_index(["species", "zone"]).amount
    assert amount[("tpex", "paracortex")] == 11.0


def test_spatial_baseline_selection_must_retain_every_patient():
    first = _spots().assign(timepoint="baseline")
    second = _spots().assign(patient_id="p2", timepoint="post")
    longitudinal = pd.concat([first, second], ignore_index=True)
    with pytest.raises(ValueError, match="removed all cells for patients"):
        aggregate_spatial_spots(
            longitudinal, "area", "timepoint", ("baseline",)
        )


def test_spatial_samples_require_an_explicit_estimand_and_equal_weighting():
    small = _spots().iloc[[0]].assign(section_id="small", spot_id="same")
    large = pd.concat([_spots().iloc[[0]].assign(section_id="large", spot_id=f"l{i}")
                       for i in range(4)], ignore_index=True)
    sections = pd.concat([small, large], ignore_index=True)
    with pytest.raises(ValueError, match="sample identifiers detected"):
        aggregate_spatial_spots(sections, "area")
    equal = aggregate_spatial_spots(
        sections, "area", sample_column="section_id", sample_aggregation="equal_samples"
    ).set_index(["species", "zone"]).amount
    summed = aggregate_spatial_spots(
        sections, "area", sample_column="section_id", sample_aggregation="sum_samples"
    ).set_index(["species", "zone"]).amount
    assert equal[("tpex", "paracortex")] == 15.0
    assert summed[("tpex", "paracortex")] == 30.0


def test_spot_ids_may_repeat_between_declared_sections_only():
    sections = pd.concat([
        _spots().iloc[[0]].assign(section_id="a", spot_id="s1"),
        _spots().iloc[[0]].assign(section_id="b", spot_id="s1"),
    ], ignore_index=True)
    states = aggregate_spatial_spots(
        sections, "area", sample_column="section_id", sample_aggregation="equal_samples"
    )
    assert states.set_index(["species", "zone"]).amount[("tpex", "paracortex")] == 6.0


def test_spatial_sample_column_must_be_a_distinct_sampling_identifier():
    with pytest.raises(ValueError, match="sample column must be distinct"):
        aggregate_spatial_spots(
            _spots(), "area", sample_column="zone", sample_aggregation="equal_samples"
        )


def test_observed_sample_mean_does_not_treat_an_unsampled_zone_as_zero():
    first = _spots().iloc[[0]].assign(section_id="a", spot_id="a1")
    second = _spots().iloc[[2]].assign(section_id="b", spot_id="b1")
    partial_sections = pd.concat([first, second], ignore_index=True)
    zero_assumption = aggregate_spatial_spots(
        partial_sections, "area", sample_column="section_id",
        sample_aggregation="equal_samples",
    ).set_index(["species", "zone"]).amount
    observed_states = aggregate_spatial_spots(
        partial_sections, "area", sample_column="section_id",
        sample_aggregation="equal_observed_samples",
    )
    observed_only = observed_states.set_index(["species", "zone"]).amount
    assert zero_assumption[("tpex", "paracortex")] == 3.0
    assert observed_only[("tpex", "paracortex")] == 6.0
    assert observed_only[("mature_dc", "afferent_lymph")] == 2.0
    support = observed_states[["zone", "sampling_units_observed",
                               "sampling_units_total", "sampling_unit_coverage"]]
    support = support.drop_duplicates().set_index("zone")
    assert support.loc["paracortex", "sampling_units_observed"] == 1
    assert support.loc["paracortex", "sampling_units_total"] == 2
    assert support.loc["paracortex", "sampling_unit_coverage"] == 0.5
    assert support.loc["medulla", "sampling_unit_coverage"] == 0.0
