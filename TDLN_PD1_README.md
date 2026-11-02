# TDLN PD-1 response baseline

This branch adds a runnable, leakage-aware baseline for converting TDLN single-cell data into patient-level immune features and predicting binary anti-PD-1 response. It preserves the original organ visualisation on `main`.

## What was actually found

Repository history and all branches were checked. At reconstruction time the repository had only `main`; no TDLN/PD-1 training implementation or matching commit existed. The existing project contains organ visualisation, CELLxGENE-derived expression summaries, and attribution rows for immune-checkpoint datasets. Those rows are provenance, not a labelled TDLN response cohort. This implementation is consequently a new baseline, not a recovered validated clinical model.

## Pipeline

1. Load cell-level CSV or optional `.h5ad` (`obs.patient_id`, `obs.cell_type`), requiring an explicit baseline slice whenever temporal columns exist.
2. Compute TDLN-restricted patient-level cell fractions and cytotoxic, exhaustion, interferon, Treg, TCF7-memory, and dendritic-activation signatures while retaining paired tissues for TCR sharing.
3. Join pre-specified binary patient labels.
4. Select logistic-regression regularization only inside nested, stratified patient-grouped cross-validation and generate one out-of-fold probability per patient.
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
tdln-pd1 calibrate --observations examples/longitudinal_observations.csv --fit tumour_growth --output artifacts/calibration
tdln-pd1 calibrate-spatial --observations examples/spatial_longitudinal_observations.csv --fit antigen_decay --output artifacts/spatial_calibration
tdln-pd1 spatial-featurize --spots examples/spatial_spots_cohort.csv --weight-column area --output artifacts/spatial_cohort_initial.csv
tdln-pd1 augment-twin --features examples/features.csv --spatial-initial-zones artifacts/spatial_cohort_initial.csv --systemic-parameters artifacts/calibration/fitted_parameters.json --spatial-parameters artifacts/spatial_calibration/fitted_spatial_parameters.json --days 42 --dose-days 0,14,28 --output artifacts/hybrid_features.csv
tdln-pd1 train --features artifacts/hybrid_features.csv --labels examples/labels.csv --output artifacts --ensemble-members 50
tdln-pd1 predict --model artifacts/model.joblib --features external_hybrid_features.csv --output external_predictions.csv
tdln-pd1 evaluate --predictions external_predictions.csv --labels external_labels.csv --output external_metrics.json
tdln-pd1 site-cv --features artifacts/hybrid_features.csv --labels labels_with_center.csv --output artifacts/site_cv
tdln-pd1 simulate --days 42 --dose-days 0,14,28 --dose 1 --parameters artifacts/calibration/fitted_parameters.json --output artifacts/tdln_trajectory.csv
tdln-pd1 spatial-featurize --spots examples/spatial_spots.csv --weight-column area --output artifacts/spatial_initial.csv
tdln-pd1 simulate-spatial --initial-zones artifacts/spatial_initial.csv --patient-id example_patient --systemic-parameters artifacts/calibration/fitted_parameters.json --spatial-parameters artifacts/spatial_calibration/fitted_spatial_parameters.json --days 14 --dose-days 0 --dose 1 --paired-counterfactual --output artifacts/spatial_patient.csv
tdln-pd1 simulate-spatial --days 14 --output artifacts/spatial_tdln.csv
tdln-pd1 sensitivity --days 42 --samples 128 --output artifacts/sensitivity
tdln-pd1 verify --output artifacts/verification.json
pytest
```

For longitudinal paired-tissue data, use the same temporal and tissue selection for both gates, for example `tdln-pd1 qc ... --time-column timepoint --baseline-values baseline --feature-tissue-pattern 'tdln|lymph'` followed by the same selection in `tdln-pd1 featurize`. Expression features default to the TDLN/lymph-node compartment; `--feature-tissue-pattern all` is an explicit opt-out. Omitting the temporal selection when a recognized time column exists is an error, preventing post-treatment information from leaking into a baseline model. QC also fails if any patient lacks cells in the chosen feature compartment.

`verify` is a fail-closed numerical qualification gate. It checks finite/non-negative systemic states, scaled zero-dose numerical equivalence, monotone default dose/occupancy direction, endpoint stability under tighter solver tolerances, per-species conservation in the reaction-free spatial network, and the prespecified checkpoint-coupling direction. It writes every measured error and tolerance to JSON and exits non-zero if any invariant fails. Zero-dose schedules restart the adaptive solver at event boundaries, so equivalence uses a documented state-scaled relative tolerance rather than a dimensionally invalid near-zero absolute tolerance. These are software and numerical checks, not biological validation.

The interactive six-zone virtual organ is published from [`docs/tdln.html`](docs/tdln.html). Its trajectories are not hand-authored: `python scripts/export_tdln_dashboard.py` deterministically exports the mechanistic and spatial simulations to a versioned local JSON file. CI regenerates the file and rejects stale dashboard data.

`spatial-featurize` is the explicit bridge from segmented spatial transcriptomics, imaging or deconvolution output into the executable six-zone model. Each input spot must carry a patient identifier, a curated anatomical `zone`, and non-negative `antigen`, `mature_dc`, `tpex` and `effector` abundance estimates; an optional area or cell-count weight converts densities into amounts. Unknown anatomy, duplicate spot identifiers, negative values and non-finite values fail closed. The command does not guess zones from coordinates, because coordinates without a tissue mask and landmarks do not identify lymph-node microanatomy. `simulate-spatial --initial-zones` then runs a selected patient's complete zone-resolved measured state rather than the population default. Abundances must be calibrated to consistent model units before patient interpretation.

Adding `--dose-days` couples the systemic anti-PD-1 pharmacokinetic/receptor-occupancy trajectory into that measured spatial state. `--paired-counterfactual` runs control and treated arms from the exact same patient initialization and writes an explicit `arm` and time-varying `pd1_occupancy` on every spatial row. A zero-amount schedule is numerically equivalent to control, while the default positive schedule increases efferent effectors in the prespecified verification scenario. Dosing is rejected in transport-only mode because receptor occupancy acts on immune reactions, not conservative transport.

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

`calibrate-spatial` accepts longitudinal `day,species,zone,amount,sd` observations and fits only explicitly named positive `SpatialParameters` in log space. It can use a selected patient's measured zone initialization and an optional anti-PD-1 dose schedule, in which case the fixed systemic receptor-occupancy trajectory is coupled during every spatial fit evaluation. Reports include weighted residuals, Jacobian rank and condition number, observation-to-parameter sufficiency, an AIC-like diagnostic and a fail-visible identifiability warning. The bundled spatial observations are synthetic and recover a known antigen-decay value; they are software verification, not biological calibration.

Fitted parameter JSON is executable rather than archival: `simulate --parameters`, `simulate-spatial --systemic-parameters/--spatial-parameters`, and `augment-twin --systemic-parameters/--spatial-parameters` load it with exact dataclass-schema and finite-positive-value checks. Whenever a parameter file is used, the command writes `<output>.parameters.json` containing its SHA-256 fingerprint and parameter class. Hybrid output also carries `qc_*_parameters_supplied` audit fields, which the model feature-role gate excludes from fitting. This closes the calibration-to-prediction chain and makes the actual parameterization independently checkable.

`sensitivity` uses seeded Latin-hypercube sampling over log-scaled parameter ranges and ranks monotonic influence on final tumour cells, tumour effectors, TDLN Tpex and peak PD-1 occupancy using Spearman correlations. These are screening statistics, not Sobol variance-decomposition indices or causal effects. The full sampled design is saved so every ranking is reproducible and inspectable.

## Mechanism–ML fusion

`augment-twin` maps each patient's measured dendritic, cytotoxic, exhaustion, TCF7-memory and immune-cell-fraction features to bounded initial states, runs matched untreated and anti-PD-1 counterfactuals, and appends mechanistic endpoints. These include receptor occupancy, Tpex/effector states and treated-minus-control tumour/effector effects. The mapping is explicit in `hybrid.py`; response labels are not accepted or used during augmentation. The resulting hybrid table can be supplied to the existing `train` command.

The same patient row now initializes the six-zone spatial twin: dendritic activation scales afferent DC input, TCF7-memory scales paracortical Tpex, cytotoxic/CD8 features scale paracortical effectors, and an optional `tumour_antigen_proxy` scales afferent antigen. When that proxy is absent the population default is retained rather than inferred from an unrelated expression signature. Treated/control efferent-effector endpoints and paracortical Tpex are appended to the hybrid feature table.

When `--spatial-initial-zones` is supplied, `augment-twin` replaces that feature-derived spatial initialization with each patient's measured/deconvolved 24-state table. The spatial and single-cell patient sets must match exactly; missing or extra patients, duplicate state coordinates, invalid anatomy, and any outcome column fail closed. Mechanistic spatial counterfactual endpoints then enter the same patient-level hybrid table used by `train`. A `qc_spatial_initialization_measured` audit flag records the source but is excluded from fitting by the feature-role gate.

## Locked external evaluation

Training writes a versioned model artifact with its exact feature schema. `predict` applies that artifact without refitting and fails on missing features. `evaluate` reports ROC AUC, average precision, Brier score, log loss, balanced accuracy, sensitivity, specificity, calibration intercept/slope and patient-level bootstrap intervals. A result is external validation only when the evaluation cohort was never used for model fitting, hyperparameter selection, threshold choice or feature design; running the command on the bundled synthetic cohort does not make it external evidence.

Training, prediction and evaluation normalize patient identifiers and fail closed on missing, duplicate or mismatched patient sets; no inner join may silently discard a patient. Model artifact version 6 stores SHA-256 hashes of normalized training identifiers plus a whole-cohort fingerprint. Locked prediction emits a per-patient overlap flag, and `evaluate` refuses the `external_locked_model` label if the artifact lacks this audit or any evaluation patient overlaps training. The hashes prevent accidental reuse; they are pseudonymous governance metadata, not anonymization against dictionary attacks on predictable identifiers.

Internal evaluation uses an outer patient-grouped fold for performance estimation and an inner patient-grouped fold to choose `C` from a fixed grid using log loss. If an outer training fold is too small to support two-class inner validation, it records an explicit fixed-`C` fallback instead of presenting it as tuned. The report includes the selected value and train/test size for every outer fold, a prespecified 0.1–0.9 decision curve, and an artifact-level final selection performed only after out-of-fold evaluation is complete.

Before any split, a label-blind role audit removes technical cell-count/QC columns such as `n_cells` and records why they were excluded. Non-numeric metadata must remain outside the model matrix. Outcome-like names (`response`, `outcome`, `RECIST`, survival labels, and related tokens) are rejected, as is any numeric feature that exactly copies or inverts the response vector. The included and excluded schema is persisted in both metrics and the locked artifact.

Artifact version 6 stores a deterministic patient-bootstrap refit ensemble. Locked prediction reports the ensemble median, 2.5th/97.5th percentiles and interval width for every patient. The final `C` is held fixed across bootstrap members so the interval measures patient-resampling instability rather than mixing in an undocumented tuning procedure. This is model-resampling uncertainty only: it does not cover assay error, population shift, causal ambiguity or mechanistic structural uncertainty.

When training is invoked through the CLI, artifact version 6 also embeds the exact hybrid-feature CSV hash and, when present, its parameter-provenance sidecar and sidecar hash. This binds the fitted predictor to the calibrated systemic/spatial parameter fingerprints used to generate its mechanistic covariates rather than relying on a detached workflow log.

The artifact also freezes a label-free training applicability profile for every model feature (median, robust scale, finite range and missingness). Locked predictions report each patient's maximum robust distance, fraction of features outside the observed training range, missing-feature fraction and a prespecified domain-shift flag. External evaluation fails closed when this audit is absent and reports in-domain coverage. These diagnostics expose unsupported extrapolation; they do not demonstrate biological transportability or justify discarding difficult patients after outcomes are known.

`site-cv` performs leave-one-center-out evaluation: each centre is predicted only by a model trained on all other centres. Sites whose training or test partition lacks one response class are reported as skipped rather than silently producing invalid AUC estimates. This is a stronger internal transportability check than random folds, but it is still not an independently collected external validation.
