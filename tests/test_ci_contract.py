from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_ci_exercises_calibrated_patient_spatial_sensitivity_and_integrity():
    workflow = (ROOT / ".github/workflows/tdln-ci.yml").read_text()
    required = (
        "tdln-pd1 spatial-sensitivity",
        "--spatial-parameters /tmp/spatial_calibration/fitted_spatial_parameters.json",
        "--systemic-parameters /tmp/calibration/fitted_parameters.json",
        "--initial-zones /tmp/spatial_initial.csv",
        "--patient-id example_patient",
        "--dose-days 0",
        "verify-sensitivity --audit /tmp/sensitivity/sensitivity_provenance.json",
        "verify-sensitivity --audit /tmp/spatial_sensitivity/spatial_sensitivity_provenance.json",
    )
    for marker in required:
        assert marker in workflow
