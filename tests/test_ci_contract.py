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
        "--paired-counterfactual",
        "patient_liver_tdln_counterfactual.json",
        "verify-sensitivity --audit /tmp/sensitivity/sensitivity_provenance.json",
        "verify-sensitivity --audit /tmp/spatial_sensitivity/spatial_sensitivity_provenance.json --spatial-parameters /tmp/spatial_calibration/fitted_spatial_parameters.json --systemic-parameters /tmp/calibration/fitted_parameters.json --initial-zones /tmp/spatial_initial.csv",
        "verify-calibration --audit /tmp/calibration/calibration_provenance.json --observations examples/longitudinal_observations.csv",
        "verify-calibration --audit /tmp/spatial_calibration/spatial_calibration_provenance.json --observations examples/spatial_longitudinal_observations.csv",
        "calibrate_liver_tdln_bridge.py --observations examples/liver_tdln_bridge_observations.csv",
        "--fit bridge.antigen_drainage_gain --output /tmp/liver_bridge_calibration",
        "verify-calibration --audit /tmp/liver_bridge_calibration/liver_tdln_bridge_provenance.json --observations examples/liver_tdln_bridge_observations.csv",
        "verify-augmentation --audit /tmp/hybrid.csv.twin-audit.json",
        "--systemic-diagnostics /tmp/calibration/calibration_diagnostics.json",
        "--spatial-diagnostics /tmp/spatial_calibration/spatial_calibration_diagnostics.json",
        "--uncertainty-samples 4 --uncertainty-seed 42",
        "calibrate-feature-state-mapping --observations examples/feature_state_mapping_observations.csv",
        "mapping_calibration/mapping_calibration_provenance.json",
        "--feature-state-mapping /tmp/mapping_calibration/fitted_feature_state_mapping.json",
        "--mapping-diagnostics /tmp/mapping_calibration/mapping_calibration_diagnostics.json",
        "evaluate --predictions /tmp/predictions.csv --model /tmp/model/model.joblib --features /tmp/external_smoke.csv",
        "verify-features --audit /tmp/features.csv.feature-audit.json --cells examples/cells.csv",
        "simulate_patient_liver_tdln.py --features /tmp/features.csv --patient-id p0",
        "patient_liver_tdln/patient_liver_tdln_audit.json",
        "featurize --cells examples/cells.csv --cell-type-vocabulary 'CD8 T,DC,Treg'",
        "verify-spatial-features --audit /tmp/spatial_cohort_initial.csv.spatial-audit.json --spots examples/spatial_spots_cohort.csv",
        "spatial-featurize --spots examples/spatial_spots_cohort.csv --weight-column area --time-column sampling_day --baseline-values 0",
        "propagate-uncertainty --parameters /tmp/calibration/fitted_parameters.json",
        "verify-uncertainty --audit /tmp/uncertainty/uncertainty_provenance.json",
        "propagate-spatial-uncertainty --parameters /tmp/spatial_calibration/fitted_spatial_parameters.json",
        "verify-spatial-uncertainty --audit /tmp/spatial_uncertainty/uncertainty_provenance.json",
        "verify-evaluation --audit /tmp/external_metrics.json.audit.json --predictions /tmp/predictions.csv --model /tmp/model/model.joblib --features /tmp/external_smoke.csv --labels /tmp/external_smoke_labels.csv --metrics /tmp/external_metrics.json",
    )
    for marker in required:
        assert marker in workflow


def test_ci_uses_exact_supported_numerical_constraints():
    workflow = (ROOT / ".github/workflows/tdln-ci.yml").read_text()
    assert "pip install -c constraints/py310-312.txt -e '.[test]'" in workflow
    constraints = (ROOT / "constraints/py310-312.txt").read_text()
    for package in ("numpy", "pandas", "scipy", "scikit-learn", "joblib", "pytest"):
        lines = [line for line in constraints.splitlines() if line.startswith(f"{package}==")]
        assert len(lines) == 1
