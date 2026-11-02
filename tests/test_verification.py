import pytest

from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.verification import verify_twin


def test_formal_verification_report_passes_all_prespecified_checks():
    report = verify_twin()
    assert report["schema_version"] == 1
    assert report["research_only"] is True
    assert report["passed"] is True
    assert len(report["checks"]) == 9
    assert all(check["passed"] for check in report["checks"])
    guards = next(check for check in report["checks"]
                  if check["name"] == "input_parameter_domain_guards")
    assert set(guards["rejected_invalid_layers"]) == {"systemic", "spatial"}
    spatial_dose = next(check for check in report["checks"]
                        if check["name"] == "spatial_zero_dose_and_dose_direction")
    assert spatial_dose["maximum_zero_dose_scaled_relative_error"] <= 1e-6
    assert [row["dose"] for row in spatial_dose["series"]] == [0, 0.25, 0.5, 1, 2]
    assert all(
        later["final_efferent_effector"] >= earlier["final_efferent_effector"]
        for earlier, later in zip(spatial_dose["series"], spatial_dose["series"][1:])
    )
    capacity = next(check for check in report["checks"]
                    if check["name"] == "spatial_tpex_niche_capacity")
    assert capacity["maximum_capacity_ratio"] <= 1.0002
    assert capacity["final_capacity_ratio"] == pytest.approx(1.0, rel=2e-4)


def test_simulator_rejects_nonfinite_inputs_and_invalid_tolerances():
    with pytest.raises(ValueError, match="finite"):
        simulate(days=float("nan"))
    with pytest.raises(ValueError, match="tolerances"):
        simulate(days=1, rtol=0)
    with pytest.raises(ValueError, match="doses"):
        simulate(days=1, doses=(Dose(float("nan"), 1),))
