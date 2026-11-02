# TDLN PD-1 response baseline

This branch adds a runnable, leakage-aware baseline for converting TDLN single-cell data into patient-level immune features and predicting binary anti-PD-1 response. It preserves the original organ visualisation on `main`.

## What was actually found

Repository history and all branches were checked. At reconstruction time the repository had only `main`; no TDLN/PD-1 training implementation or matching commit existed. The existing project contains organ visualisation, CELLxGENE-derived expression summaries, and attribution rows for immune-checkpoint datasets. Those rows are provenance, not a labelled TDLN response cohort. This implementation is consequently a new baseline, not a recovered validated clinical model.

## Pipeline

1. Load cell-level CSV or optional `.h5ad` (`obs.patient_id`, `obs.cell_type`).
2. Compute patient-level cell fractions and cytotoxic, exhaustion, interferon, Treg, TCF7-memory, and dendritic-activation signatures.
3. Join pre-specified binary patient labels.
4. Evaluate regularised logistic regression with stratified, patient-grouped cross-validation.
5. Write cross-validated predictions, metrics, and a fitted model artifact.

The mechanistic digital-twin layer additionally simulates the tumour–blood–TDLN axis, antigen transport, DC migration, Tpex self-renewal/differentiation, effector export, exhaustion/reinvigoration, anti-PD-1 pharmacokinetics, TDLN exposure, and receptor occupancy. See [`docs/TDLN_MODEL_CARD.md`](docs/TDLN_MODEL_CARD.md) for evidence, verification and the validation ladder.

The zone-resolved layer is a directed reaction–transport network across afferent lymph, subcapsular sinus, paracortex, follicle, medulla and efferent lymph. It transports antigen, mature DCs, Tpex and effectors, with antigen capture in the sinus and priming/differentiation in the paracortex. `--transport-only` disables reactions for mass-conservation audits.

The six-zone network can consume the systemic model's time-varying PD-1 receptor occupancy. This makes treated and untreated spatial counterfactuals distinct while preserving exact zero-occupancy equivalence and transport-only mass conservation. Coupling gains remain calibration parameters rather than fitted clinical effects.

Grouping is performed at the patient level; cells from one patient never become independent train/test observations.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
tdln-pd1 featurize --cells examples/cells.csv --output examples/features.csv
tdln-pd1 qc --cells examples/cells.csv --labels examples/labels.csv --manifest examples/cohort_manifest.json --output artifacts/qc.json
tdln-pd1 train --features examples/features.csv --labels examples/labels.csv --output artifacts
tdln-pd1 predict --model artifacts/model.joblib --features external_features.csv --output external_predictions.csv
tdln-pd1 evaluate --predictions external_predictions.csv --labels external_labels.csv --output external_metrics.json
tdln-pd1 site-cv --features features.csv --labels labels_with_center.csv --output artifacts/site_cv
tdln-pd1 simulate --days 42 --dose-days 0,14,28 --dose 1 --output artifacts/tdln_trajectory.csv
tdln-pd1 calibrate --observations examples/longitudinal_observations.csv --fit tumour_growth --output artifacts/calibration
tdln-pd1 augment-twin --features examples/features.csv --days 42 --dose-days 0,14,28 --output artifacts/hybrid_features.csv
tdln-pd1 simulate-spatial --days 14 --output artifacts/spatial_tdln.csv
tdln-pd1 sensitivity --days 42 --samples 128 --output artifacts/sensitivity
pytest
```

The interactive six-zone virtual organ is published from [`docs/tdln.html`](docs/tdln.html). Its trajectories are not hand-authored: `python scripts/export_tdln_dashboard.py` deterministically exports the mechanistic and spatial simulations to a versioned local JSON file. CI regenerates the file and rejects stale dashboard data.

Install `.[single-cell,test]` to read `.h5ad` directly.

`qc` is the mandatory pre-training gate for imported cohorts. It checks cryptographic provenance, licensing/source metadata, patient matching, duplicate labels, binary outcomes, class balance, cell counts, signature-gene coverage, paired TDLN–tumour samples and clonotype assignment. A failed gate exits non-zero so CI and workflow systems cannot silently continue.

## Interpretation and limitations

- The included demo is synthetic and proves software execution only.
- No clinical-performance claim is made. External cohorts, protocol harmonisation, class definitions, calibration, confidence intervals, and prospective validation are required before translational use.
- Small cohorts should report uncertainty and avoid feature/model selection outside nested cross-validation.
- Batch, cancer type, treatment regimen, sampling time, tissue site, and centre are likely confounders and should be added as cohort metadata when available.
- This software is research-only and is not a medical device.

## Patient/cohort data assimilation

Calibration observations use long format: `day,variable,value,sd`. Variables must be explicit simulator states or documented derived endpoints. Positive parameters are fitted in log space with uncertainty-weighted robust residuals. The output includes fitted parameters, residuals, weighted RMSE, an identifiability warning, and an AIC-like diagnostic. Parameter-uncertainty ensembles propagate log-normal uncertainty into tumour and immune outcome intervals. The bundled observation file is synthetic and exists only to verify recovery of a known parameter.

`sensitivity` uses seeded Latin-hypercube sampling over log-scaled parameter ranges and ranks monotonic influence on final tumour cells, tumour effectors, TDLN Tpex and peak PD-1 occupancy using Spearman correlations. These are screening statistics, not Sobol variance-decomposition indices or causal effects. The full sampled design is saved so every ranking is reproducible and inspectable.

## Mechanism–ML fusion

`augment-twin` maps each patient's measured dendritic, cytotoxic, exhaustion, TCF7-memory and immune-cell-fraction features to bounded initial states, runs matched untreated and anti-PD-1 counterfactuals, and appends mechanistic endpoints. These include receptor occupancy, Tpex/effector states and treated-minus-control tumour/effector effects. The mapping is explicit in `hybrid.py`; response labels are not accepted or used during augmentation. The resulting hybrid table can be supplied to the existing `train` command.

## Locked external evaluation

Training writes a versioned model artifact with its exact feature schema. `predict` applies that artifact without refitting and fails on missing features. `evaluate` reports ROC AUC, average precision, Brier score, log loss, balanced accuracy, sensitivity, specificity, calibration intercept/slope and patient-level bootstrap intervals. A result is external validation only when the evaluation cohort was never used for model fitting, hyperparameter selection, threshold choice or feature design; running the command on the bundled synthetic cohort does not make it external evidence.

`site-cv` performs leave-one-center-out evaluation: each centre is predicted only by a model trained on all other centres. Sites whose training or test partition lacks one response class are reported as skipped rather than silently producing invalid AUC estimates. This is a stronger internal transportability check than random folds, but it is still not an independently collected external validation.
