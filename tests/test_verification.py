import pytest

from tdln_pd1.simulator import Dose, simulate
from tdln_pd1.verification import verify_twin


def test_formal_verification_report_passes_all_prespecified_checks():
    report = verify_twin()
    assert report["schema_version"] == 1
    assert report["research_only"] is True
    assert report["passed"] is True
    assert len(report["checks"]) == 7
    assert all(check["passed"] for check in report["checks"])
    guards = next(check for check in report["checks"]
                  if check["name"] == "input_parameter_domain_guards")
    assert set(guards["rejected_invalid_layers"]) == {"systemic", "spatial"}


def test_simulator_rejects_nonfinite_inputs_and_invalid_tolerances():
    with pytest.raises(ValueError, match="finite"):
        simulate(days=float("nan"))
    with pytest.raises(ValueError, match="tolerances"):
        simulate(days=1, rtol=0)
    with pytest.raises(ValueError, match="doses"):
        simulate(days=1, doses=(Dose(float("nan"), 1),))
