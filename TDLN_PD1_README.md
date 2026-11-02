# TDLN PD-1 response baseline

This branch adds a runnable, leakage-aware baseline for converting TDLN single-cell data into patient-level immune features and predicting binary anti-PD-1 response. It preserves the original organ visualisation on `main`.

## What was actually found

Repository history and all branches were checked. At reconstruction time the repository had only `main`; no TDLN/PD-1 training implementation or matching commit existed. The existing project contains organ visualisation, CELLxGENE-derived expression summaries, and attribution rows for immune-checkpoint datasets. Those rows are provenance, not a labelled TDLN response cohort. This implementation is consequently a new baseline, not a recovered validated clinical model.

## Pipeline

1. Load cell-level CSV or optional `.h5ad` (`obs.patient_id`, `obs.cell_type`), preserving tissue, clonotype, recognized temporal fields and any explicitly requested custom time field; both formats require an explicit baseline slice whenever temporal columns exist.
2. Compute TDLN-restricted patient-level cell fractions and cytotoxic, exhaustion, interferon, Treg, TCF7-memory, and dendritic-activation signatures while retaining paired tissues for TCR sharing. TCR features use one fail-closed assignment filter and report unique-clone sharing, TDLN cell-weighted sharing and set Jaccard overlap; missing/unassigned placeholders can never count as shared lineages.
3. Join pre-specified binary patient labels.
4. Select logistic-regression regularization only inside nested, stratified patient-grouped cross-validation and generate one out-of-fold probability per patient.
5. Write cross-validated predictions, metrics, and a fitted model artifact.

The mechanistic digital-twin layer additionally simulates the tumour–blood–TDLN axis, antigen transport, DC migration, Tpex self-renewal/differentiation, effector export, exhaustion/reinvigoration, anti-PD-1 pharmacokinetics, TDLN exposure, and receptor occupancy. See [`docs/TDLN_MODEL_CARD.md`](docs/TDLN_MODEL_CARD.md) for evidence, verification and the validation ladder.

The zone-resolved layer is a directed reaction–transport network across afferent lymph, subcapsular sinus, paracortex, follicle, medulla and efferent lymph. It transports antigen, mature DCs, Tpex and effectors, with antigen capture in the sinus and priming/differentiation in the paracortex. `--transport-only` disables reactions for mass-conservation audits.

The six-zone network can consume the systemic model's time-varying PD-1 receptor occupancy. This makes treated and untreated spatial counterfactuals distinct while preserving exact zero-occupancy equivalence and transport-only mass conservation. Coupling gains remain calibration parameters rather than fitted clinical effects.

Systemic-to-spatial coupling is fail-closed: every supplied timeline row must contain a numeric, finite, unique time and a complete finite receptor occupancy within `[0,1]`. Invalid rows are never dropped, and callable occupancy sources are rejected rather than clipped when they leave the physical range.

Grouping is performed at the patient level; cells from one patient never become independent train/test observations.

Single-cell expression columns preserve explicit missing values for downstream fold-internal imputation, while supplied non-numeric or infinite signature-gene values fail before aggregation. Negative finite values remain valid because some normalized/scaled assays legitimately produce them. Cell-type fraction names use deterministic normalization, and distinct labels that collapse to the same feature name are rejected rather than overwriting one another; cohort-specific aliases should be harmonized upstream under a recorded ontology map.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[test]'
tdln-pd1 featurize --cells examples/cells.csv --output examples/features.csv
tdln-pd1 qc --cells examples/cells.csv --labels examples/labels.csv --manifest examples/cohort_manifest.json --output artifacts/qc.json
tdln-pd1 calibrate --observations examples/longitudinal_observations.csv --fit tumour_growth --profile-points 7 --output artifacts/calibration
tdln-pd1 calibrate-spatial --observations examples/spatial_longitudinal_observations.csv --fit antigen_decay --profile-points 7 --output artifacts/spatial_calibration
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

For longitudinal paired-tissue data, use the same temporal and tissue selection for both gates, for example `tdln-pd1 qc ... --time-column timepoint --baseline-values baseline --feature-tissue-pattern 'tdln|lymph'` followed by the same selection in `tdln-pd1 featurize`. Expression features default to the TDLN/lymph-node compartment; `--feature-tissue-pattern all` is an explicit opt-out. Omitting the temporal selection when a recognized time column exists is an error, preventing post-treatment information from leaking into a baseline model. AnnData ingestion preserves these `obs` fields before applying exactly the same gate, including a custom field named by `--time-column`; it may not silently discard time metadata. Missing or empty patient/cell-type identifiers fail during loading in both formats. QC also fails if any patient lacks cells in the chosen feature compartment.

`verify` is a fail-closed numerical qualification gate. It checks finite/non-negative systemic states, scaled zero-dose numerical equivalence, monotone default dose/occupancy direction, endpoint stability under tighter solver tolerances, per-species conservation in the reaction-free spatial network, and the prespecified checkpoint-coupling direction. It writes every measured error and tolerance to JSON and exits non-zero if any invariant fails. Zero-dose schedules restart the adaptive solver at event boundaries, so equivalence uses a documented state-scaled relative tolerance rather than a dimensionally invalid near-zero absolute tolerance. These are software and numerical checks, not biological validation.

Both ODE layers inspect raw solver states before writing output. Only negative drift within a documented state-scaled `10 × (atol + rtol × scale)` numerical budget is clipped to zero; larger negative excursions, non-finite states, or material systemic receptor occupancy outside `[0,1]` abort the simulation. Output clipping therefore cannot conceal solver failure under extreme parameterizations.

The interactive six-zone virtual organ is published from [`docs/tdln.html`](docs/tdln.html). Its trajectories are not hand-authored: `python scripts/export_tdln_dashboard.py` deterministically exports the mechanistic and spatial simulations to a versioned local JSON file. CI regenerates the file and rejects stale dashboard data.

`spatial-featurize` is the explicit bridge from segmented spatial transcriptomics, imaging or deconvolution output into the executable six-zone model. Each input spot must carry a patient identifier, a curated anatomical `zone`, and non-negative `antigen`, `mature_dc`, `tpex` and `effector` abundance estimates; an optional area or cell-count weight converts densities into amounts. Unknown anatomy, duplicate spot identifiers, negative values and non-finite values fail closed. The command does not guess zones from coordinates, because coordinates without a tissue mask and landmarks do not identify lymph-node microanatomy. It emits an explicit 24-state grid for every patient: zones with source spots carry `zone_observed=true`, while absent zones are explicit zero-valued modelling assumptions with `zone_observed=false`. The simulator rejects partial long tables instead of silently creating their missing states. `simulate-spatial --initial-zones` then runs this auditable zone-resolved state rather than the population default. Abundances must be calibrated to consistent model units before patient interpretation.

Adding `--dose-days` couples the systemic anti-PD-1 pharmacokinetic/receptor-occupancy trajectory into that measured spatial state. `--paired-counterfactual` runs control and treated arms from the exact same patient initialization and writes an explicit `arm` and time-varying `pd1_occupancy` on every spatial row. A zero-amount schedule is numerically equivalent to control, while the default positive schedule increases efferent effectors in the prespecified verification scenario. Dosing is rejected in transport-only mode because receptor occupancy acts on immune reactions, not conservative transport.

Systemic trajectories use right-continuous bolus semantics: the row at a dosing boundary is the post-dose state and `administered_dose` records the coalesced amount. A schedule at exactly the terminal simulation time fails instead of being accepted and silently ignored.

Install `.[single-cell,test]` to read `.h5ad` directly.

`qc` is the mandatory pre-training gate for imported cohorts. It checks cryptographic provenance, licensing/source metadata, patient matching, duplicate labels, binary outcomes, class balance, cell counts, signature-gene coverage, paired TDLN–tumour samples and clonotype assignment. A failed gate exits non-zero so CI and workflow systems cannot silently continue.

Patient and outcome contracts are shared by training, site validation and locked external evaluation: identifiers must be present, non-empty and unique, while every response must be explicitly observed and equal to integer 0 or 1. Missing outcomes are never treated as a third class or silently removed.

## Interpretation and limitations

- The included demo is synthetic and proves software execution only.
- No clinical-performance claim is made. External cohorts, protocol harmonisation, class definitions, calibration, confidence intervals, and prospective validation are required before translational use.
- Small cohorts should report uncertainty and avoid feature/model selection outside nested cross-validation.
- Batch, cancer type, treatment regimen, sampling time, tissue site, and centre are likely confounders and should be added as cohort metadata when available.
- This software is research-only and is not a medical device.

## Patient/cohort data assimilation

Calibration observations use long format: `day,variable,value,sd`. Variables must be explicit simulator states or documented derived endpoints. Positive parameters are fitted in log space with uncertainty-weighted robust residuals. The output includes fitted parameters, residuals, weighted RMSE, Jacobian rank and condition number, parameters touching the search boundary, approximate local log-scale intervals when estimable, an identifiability warning, and an AIC-like diagnostic. Optional `--profile-points` first establishes a separate Gaussian-error weighted-SSE optimum, then fixes each target parameter and re-optimizes all nuisance parameters with that same objective. The report distinguishes a resolved asymptotic 95% profile boundary from support that merely reaches the edge of the requested fold range; the Gaussian profile optimum is recorded separately from the primary robust estimate. Parameter-uncertainty ensembles propagate log-normal uncertainty into tumour and immune outcome intervals. The bundled observation file is synthetic and exists only to verify recovery of a known parameter.

Calibration tables must be non-empty and contain one unambiguous observation per model coordinate (`day,variable` systemically and `day,species,zone` spatially). Repeated assays must be aggregated upstream under a declared measurement model rather than silently receiving extra optimizer weight. Log-space fit parameter names must be unique and their supplied base values finite and strictly positive; fixed parameters may still be zero when physically meaningful.

`calibrate-spatial` accepts longitudinal `day,species,zone,amount,sd` observations and fits only explicitly named positive `SpatialParameters` in log space. It can use a selected patient's measured zone initialization and an optional anti-PD-1 dose schedule, in which case the fixed systemic receptor-occupancy trajectory is coupled during every spatial fit evaluation. Reports use the same rank, condition, search-boundary and approximate local-interval diagnostics as the systemic layer, alongside observation-to-parameter sufficiency and an AIC-like score. The bundled spatial observations are synthetic and recover a known antigen-decay value; they are software verification, not biological calibration.

Fitted parameter JSON is executable rather than archival: `simulate --parameters`, `simulate-spatial --systemic-parameters/--spatial-parameters`, and `augment-twin --systemic-parameters/--spatial-parameters` load it with exact dataclass-schema and finite-positive-value checks. Whenever a parameter file is used, the command writes `<output>.parameters.json` containing its SHA-256 fingerprint and parameter class. Hybrid output also carries `qc_*_parameters_supplied` audit fields, which the model feature-role gate excludes from fitting. This closes the calibration-to-prediction chain and makes the actual parameterization independently checkable.

`sensitivity` uses seeded Latin-hypercube sampling over log-scaled parameter ranges and ranks monotonic influence on final tumour cells, tumour effectors, TDLN Tpex and peak PD-1 occupancy using Spearman correlations. These are screening statistics, not Sobol variance-decomposition indices or causal effects. The full sampled design is saved so every ranking is reproducible and inspectable.

## Mechanism–ML fusion

`augment-twin` maps each patient's measured dendritic, cytotoxic, exhaustion, TCF7-memory and immune-cell-fraction features to bounded initial states, runs matched untreated and anti-PD-1 counterfactuals, and appends mechanistic endpoints. These include receptor occupancy, Tpex/effector states and treated-minus-control tumour/effector effects. The mapping is explicit in `hybrid.py`; response labels are not accepted or used during augmentation. The resulting hybrid table can be supplied to the existing `train` command.

The label-free boundary is enforced before simulation on both inputs: patient feature tables and measured spatial-state tables reject response/responder, outcome, label, RECIST, PFS, progression-free and overall-survival-like column names. Outcomes therefore cannot be carried through the mechanistic transform as inert-looking columns even if a downstream caller bypasses the standard training command.

The same patient row now initializes the six-zone spatial twin: dendritic activation scales afferent DC input, TCF7-memory scales paracortical Tpex, cytotoxic/CD8 features scale paracortical effectors, and an optional `tumour_antigen_proxy` scales afferent antigen. When that proxy is absent the population default is retained rather than inferred from an unrelated expression signature. Treated/control efferent-effector endpoints and paracortical Tpex are appended to the hybrid feature table.

Optional absent/NaN mechanistic inputs use their documented defaults. A supplied measurement is never silently repaired: immune signatures must be finite and may be signed so centered/scaled assays remain compatible; their state multipliers are explicitly bounded. Antigen proxies must be finite and non-negative, cell fractions must lie in `[0,1]`, and non-numeric values fail before simulation. Assay-specific signature scaling still requires cohort calibration before biological interpretation.

When `--spatial-initial-zones` is supplied, `augment-twin` replaces that feature-derived spatial initialization with each patient's measured/deconvolved 24-state table. The spatial and single-cell patient sets must match exactly; missing or extra patients, incomplete or duplicate state coordinates, invalid anatomy, and any outcome column fail closed. Mechanistic spatial counterfactual endpoints then enter the same patient-level hybrid table used by `train`. `qc_spatial_initialization_measured` and `qc_spatial_zone_coverage` record the source and fraction of anatomically observed zones; both are excluded from fitting by the feature-role gate.

## Locked external evaluation

Training writes a versioned model artifact with its exact feature schema and decision threshold. `predict` applies that artifact without refitting, emits the locked threshold on every row and fails on missing features. `evaluate` requires one consistent valid threshold and uses it for all threshold-dependent metrics rather than silently reverting to 0.5. It reports ROC AUC, average precision, Brier score, log loss, balanced accuracy, sensitivity, specificity and calibration intercept/slope. Bootstrap intervals cover every reported performance metric except the two calibration coefficients, whose unstable small-sample resampling distribution is not presented as a reliable interval. A result is external validation only when the evaluation cohort was never used for model fitting, hyperparameter selection, threshold choice or feature design; running the command on the bundled synthetic cohort does not make it external evidence.

Training, prediction and evaluation normalize patient identifiers and fail closed on missing, duplicate or mismatched patient sets; no inner join may silently discard a patient. Model artifact version 7 stores SHA-256 hashes of normalized training identifiers plus a whole-cohort fingerprint. Locked prediction emits a per-patient overlap flag, and `evaluate` refuses the `external_locked_model` label if the artifact lacks this audit or any evaluation patient overlaps training. The hashes prevent accidental reuse; they are pseudonymous governance metadata, not anonymization against dictionary attacks on predictable identifiers.

The locked model matrix has the same measurement contract in training, site validation and external prediction: biological and mechanistic features must be numeric, and supplied values must be finite. Missing values remain explicit and are imputed by the training-fold or fitted-artifact median; positive or negative infinity and string placeholders fail before model execution rather than being silently coerced.

Patient identifiers alone are insufficient because a development cohort can be renamed. Version 7 therefore also hashes every locked-schema training feature vector independently of patient ID and stores a sorted whole-feature-cohort fingerprint. Prediction reports exact training-row reuse, and external evaluation fails if any row matches, even when all patient names differ. This detects exact reuse, not near-duplicate samples, shared biospecimens, batch leakage or related patients; those still require cohort governance and provenance review. Feature hashes are integrity metadata, not anonymization, and model artifacts remain governed cohort data.

Internal evaluation uses an outer patient-grouped fold for performance estimation and an inner patient-grouped fold to choose `C` from a fixed grid using log loss. If an outer training fold is too small to support two-class inner validation, it records an explicit fixed-`C` fallback instead of presenting it as tuned. The report includes the selected value and train/test size for every outer fold, a prespecified 0.1–0.9 decision curve, and an artifact-level final selection performed only after out-of-fold evaluation is complete.

Before any split, a label-blind role audit removes technical cell-count/QC columns such as `n_cells` and records why they were excluded. Non-numeric metadata must remain outside the model matrix. Outcome-like names (`response`, `outcome`, `RECIST`, survival labels, and related tokens) are rejected, as is any numeric feature that exactly copies or inverts the response vector. The included and excluded schema is persisted in both metrics and the locked artifact.

Artifact version 7 stores a deterministic patient-bootstrap refit ensemble. Locked prediction reports the ensemble median, 2.5th/97.5th percentiles and interval width for every patient. The final `C` is held fixed across bootstrap members so the interval measures patient-resampling instability rather than mixing in an undocumented tuning procedure. This is model-resampling uncertainty only: it does not cover assay error, population shift, causal ambiguity or mechanistic structural uncertainty.

When training is invoked through the CLI, artifact version 7 also embeds the exact hybrid-feature CSV hash and, when present, its parameter-provenance sidecar and sidecar hash. This binds the fitted predictor to the calibrated systemic/spatial parameter fingerprints used to generate its mechanistic covariates rather than relying on a detached workflow log.

The artifact also freezes a label-free training applicability profile for every model feature (median, robust scale, finite range and missingness). Locked predictions report each patient's maximum robust distance, fraction of features outside the observed training range, missing-feature fraction and a prespecified domain-shift flag. External evaluation fails closed when this audit is absent and reports in-domain coverage. These diagnostics expose unsupported extrapolation; they do not demonstrate biological transportability or justify discarding difficult patients after outcomes are known.

External metrics always use the complete submitted cohort as the primary analysis, including domain-shifted patients. A separate, prespecified in-domain sensitivity analysis is emitted only when at least four in-domain patients remain and each response class has at least two; otherwise its status is `not_evaluable`. The report records its locked selection rule, coverage and explicitly secondary role. This prevents applicability filtering from silently improving the headline result after outcomes are inspected.

When labels contain `center`, the pooled primary result is accompanied by a descriptive centre-stratified table. Every centre remains listed with sample size, responder count and prevalence; full metrics are calculated only when both response classes occur, while single-class centres are explicitly `not_evaluable`. The range of AUC and Brier score across evaluable centres exposes heterogeneity but is not a replacement primary endpoint or evidence of transportability.

`site-cv` performs leave-one-center-out evaluation: each centre is predicted only by a model trained on all other centres. Within every training partition, regularization is selected by a second stratified patient-grouped cross-validation scored by log loss; the held-out centre never participates in that choice. Candidate grouped splits are inspected before fitting and the fold count is reduced until every training and validation partition contains both response classes. If no valid two-class split exists, a fixed-`C` fallback is recorded rather than described as tuned. The same validated-split rule protects outer and final-selection folds. A held-out centre containing only one response class is still predicted in full and retained in pooled proper-score and threshold analyses; only its centre-specific discrimination metrics are marked `not_evaluable`. A centre is skipped only when the complementary training partition has one class and therefore cannot fit the classifier. This avoids conditioning cohort inclusion on held-out outcomes. It is a stronger internal transportability check than random folds, but it is still not an independently collected external validation.

The CLI creates a cryptographic chain of custody for locked evaluation. `predict` writes `<predictions>.audit.json`, binding the SHA-256 of the versioned model artifact and input feature table to the exact prediction CSV. `evaluate` requires and verifies that sidecar before reading outcomes, embeds the model/feature audit hashes in the metrics, and writes `<metrics>.audit.json` binding the prediction audit, raw label file, analysis settings and exact metric report. Replacing predictions after inference therefore fails closed. These hashes establish file integrity and reproducibility only: they cannot prove that a cohort is genuinely independent or that the analysis was prospectively registered.
