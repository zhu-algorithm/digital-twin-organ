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
        "verify-sensitivity --audit /tmp/spatial_sensitivity/spatial_sensitivity_provenance.json --spatial-parameters /tmp/spatial_calibration/fitted_spatial_parameters.json --systemic-parameters /tmp/calibration/fitted_parameters.json --initial-zones /tmp/spatial_initial.csv",
        "verify-calibration --audit /tmp/calibration/calibration_provenance.json --observations examples/longitudinal_observations.csv",
        "verify-calibration --audit /tmp/spatial_calibration/spatial_calibration_provenance.json --observations examples/spatial_longitudinal_observations.csv",
        "verify-augmentation --audit /tmp/hybrid.csv.twin-audit.json",
        "verify-features --audit /tmp/features.csv.feature-audit.json --cells examples/cells.csv",
        "verify-spatial-features --audit /tmp/spatial_cohort_initial.csv.spatial-audit.json --spots examples/spatial_spots_cohort.csv",
        "propagate-uncertainty --parameters /tmp/calibration/fitted_parameters.json",
        "verify-uncertainty --audit /tmp/uncertainty/uncertainty_provenance.json",
        "propagate-spatial-uncertainty --parameters /tmp/spatial_calibration/fitted_spatial_parameters.json",
        "verify-spatial-uncertainty --audit /tmp/spatial_uncertainty/uncertainty_provenance.json",
    )
    for marker in required:
        assert marker in workflow
