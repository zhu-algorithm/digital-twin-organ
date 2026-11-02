# TDLN digital-twin model card

## Intended use

Research simulation of the tumour–blood–tumour-draining-lymph-node axis, generation of falsifiable hypotheses, and evaluation of analysis pipelines. It is **not** validated for diagnosis, drug selection, dose selection, or individual clinical decisions.

## Architecture

The executable core currently contains 12 non-negative state variables spanning tumour burden and antigen, migratory mature dendritic cells, naive CD8 T cells, TCF1-like progenitor-exhausted cells (Tpex), TDLN and tumour effectors, exhausted tumour T cells, Tregs, plasma/TDLN antibody exposure, and PD-1 receptor occupancy. Bolus dosing is integrated piecewise so discontinuities are explicit rather than hidden inside the ODE solver. Output trajectories are right-continuous at dosing boundaries, retain the post-dose state, and expose `administered_dose`; a dose at the terminal time is rejected because it has no simulated follow-up interval.

The spatial ontology separates afferent lymph, subcapsular sinus, paracortex, follicle, medulla, and efferent lymph. Single-cell annotations can be mapped to these zones; unresolved labels remain unresolved instead of being guessed.

The executable spatial model uses directed conservative transport between these zones. Antigen and migratory DCs enter through afferent lymph; antigen capture occurs in the subcapsular sinus; antigen-bearing DCs reach the paracortex; Tpex expand and differentiate there; effectors pass through medulla into efferent lymph. The implementation exposes a reaction-free mode in which every species must be conserved, separating numerical transport verification from biological source/sink assumptions.

Segmented spatial-omics or imaging tables can initialize all 24 species-by-zone states through an explicit patient/zone/species bridge. Anatomy must be curated upstream; coordinates alone are never treated as sufficient evidence for a zone. Spot areas or cell counts may weight non-negative deconvolved abundances, and invalid zones, duplicate spots or non-finite/negative estimates are rejected. Aggregation always emits the complete grid with a per-zone observation flag: missing anatomical zones become explicit zero-valued assumptions, never indistinguishable from measured zeros, and their coverage is carried into hybrid QC. Direct partial state tables are rejected. This supplies a reproducible data interface, not an assay-to-cell-count calibration; input units still require modality-specific validation.

The same measured 24-state initialization can be run as matched untreated and anti-PD-1 arms. Only the systemic drug schedule differs: its time-varying receptor occupancy drives the spatial reaction terms, and both arm and occupancy are retained in the output. Zero-dose equivalence and the prespecified positive-dose effector direction are tested; these are implementation checks, not evidence of an individual treatment effect.

The systemic and spatial layers are coupled through the simulated TDLN PD-1 receptor occupancy trajectory. Occupancy increases the hypothesised Tpex expansion and differentiation reaction rates but never alters conservative transport terms. The coupling gains are explicit parameters, default to zero effect at zero occupancy, and require calibration; the directional response is a software invariant, not clinical validation.

## Evidence anchors

- Four-compartment QSP precedent (blood, peripheral tissue, tumour, TDLN) and explicit antigen/APC/T-cell/checkpoint modules: [Jusko et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11138134/) and [Milberg et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6677731/).
- Antibody exposure gradients across primary tumour and TDLN, calibrated against immuno-PET: [Li & Hung, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9036243/).
- TDLN PD-1+TCF-1+ reservoir contributing to tumour control after checkpoint treatment: [Nature Communications, 2025](https://www.nature.com/articles/s41467-025-58510-1).
- Human paired scRNA/TCR evidence that lymph-node Tpex expand after anti-PD-1 and share clonotypes with intratumoural exhausted CD8 T cells: [Nature Communications, 2026](https://www.nature.com/articles/s41467-026-70751-2).
- Anatomically based antigen redistribution and lymph-flow modelling after draining-node removal: [Distant lymph nodes compensate, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10541098/).

These sources motivate model structure. They do not validate the default numerical parameters.

## Verification completed

- State schema and finite/non-negative trajectories.
- Fail-closed finite/non-negative parameter and initial-state domains before numerical integration.
- Explicit antibody transfer to TDLN and increased receptor occupancy.
- Higher tumour effector abundance after anti-PD-1 exposure in the default scenario.
- Dose coalescing invariance for simultaneous boluses.
- Lymph-node microanatomy zone coverage and deterministic cell-zone mapping.
- Existing single-cell feature and prediction tests remain passing.
- Synthetic-truth parameter recovery for longitudinal data assimilation.
- Synthetic-truth recovery and Jacobian identifiability diagnostics for longitudinal spatial calibration.
- Deterministic uncertainty propagation under a fixed random seed.
- Label-free, patient-aligned mechanistic counterfactual feature generation.
- Per-species mass conservation in the reaction-free spatial transport network.
- Directional antigen, DC and effector movement across six TDLN zones.
- Validated spatial-omics aggregation and exact day-zero initialization of all zone/species states.
- Matched spatial treatment/control counterfactuals with explicit receptor-occupancy provenance and zero-dose equivalence.
- TCR clone diversity, dominance and paired TDLN–tumour lineage sharing by unique-clone fraction, TDLN cell-weighted fraction and Jaccard overlap, with a common invalid-assignment filter.
- Locked-model prediction with feature-schema enforcement, calibration metrics and patient-level bootstrap intervals.
- Reproducible global parameter screening with Latin-hypercube sampling and outcome-specific rank correlations.
- Leave-one-centre-out transportability evaluation retains every predictable held-out patient. Single-class held-out centres remain in pooled proper-score and threshold analyses while their centre-specific discrimination is explicitly non-evaluable; only centres whose complementary training partition cannot fit a two-class model are skipped.
- Fold-internal regularization selection for every evaluable held-out centre, with explicit fixed-value fallback when inner tuning is impossible.
- Cross-version-safe grouped splitting that rejects any candidate train/validation partition lacking both response classes and reduces fold count before fallback.
- Fail-closed cohort QC and SHA-256 provenance verification before training.
- Reproducible browser dashboard export: six-zone spatial states and treated/control trajectories are generated directly from the executable model and checked byte-for-structure in CI.
- A machine-readable, fail-closed verification report covering zero-dose limits, default dose direction, tightened-tolerance convergence, finite/non-negative states, spatial mass conservation and coupled checkpoint direction.

## Required validation ladder

1. **Unit verification:** dimensions, solver convergence, mass/flow balance, limit cases.
2. **Module calibration:** antibody PK/occupancy, antigen/DC migration, Tpex proliferation and egress, tumour killing.
3. **Multimodal data assimilation:** matched TDLN/tumour scRNA-seq, TCR clonotypes, spatial transcriptomics, imaging and clinical PK.
4. **Retrospective external validation:** locked parameters and endpoints on independent centres.
5. **Prospective validation:** preregistered protocol and uncertainty-aware predictions.

Until steps 2–5 are completed with real data, default trajectories are demonstrations of model behaviour, not estimates of a patient outcome.

The prediction report labels internal cross-validation separately from locked external evaluation. Discrimination alone is insufficient: Brier score, log loss, calibration intercept/slope, sensitivity/specificity and bootstrap intervals are emitted together. Threshold-dependent metrics use the single decision threshold carried by the locked artifact and every prediction row; a missing, variable or invalid threshold fails closed. Percentile intervals cover ROC AUC, average precision, Brier score, log loss, balanced accuracy, sensitivity and specificity; calibration intercept and slope remain point diagnostics because their small-sample resampling distributions can be pathologically unstable. When a `center` column is present, uncertainty intervals resample entire centres with replacement and retain all patients from every sampled centre; otherwise they use patient-level resampling. The chosen bootstrap unit, group count, requested draws and seed are persisted in the report. Centre bootstrap represents clustering more honestly but remains unreliable with very few centres. The software cannot prove that a supplied cohort is independent; cohort provenance and analysis preregistration remain governance requirements.

Regularization strength is selected inside each outer patient-grouped fold using an inner patient-grouped grid search scored by log loss. Every patient receives exactly one outer-fold probability, and the per-fold selection record is persisted. A prespecified decision curve reports model, treat-all and treat-none net benefit from thresholds 0.1–0.9. These calculations expose potential utility and overfitting; they do not choose a clinical threshold or establish clinical benefit.

Feature roles are audited before fitting. Technical cell-count/QC columns are excluded and recorded, non-numeric metadata cannot enter the matrix, outcome-like column names are rejected, and exact copies or inversions of the response vector cause a hard failure. This catches direct leakage but cannot prove that every apparently biological feature was measured before outcome ascertainment; provenance and protocol review remain required.

Patient-set integrity is fail-closed. Training and evaluation require exact, unique, non-empty patient alignment and a complete binary response for every patient rather than silently dropping non-matches or missing outcomes. Versioned artifacts contain SHA-256 training-identifier hashes and a cohort fingerprint; locked external evaluation rejects any detected training overlap or legacy artifact without an overlap audit. These hashes reduce accidental leakage but are not anonymous when source identifiers are predictable, so artifacts must remain under the cohort's data-governance controls.

Artifact version 7 freezes the label-free training support of each model feature, patient-independent hashes of every training feature vector, and the exact training-file hash plus any calibrated-parameter sidecar and its hash. Locked evaluation rejects both patient-ID overlap and exact feature-row reuse after renaming. Every prediction also reports robust distance, range extrapolation, missing-feature fraction and a prespecified applicability flag. These controls expose exact reuse and unsupported extrapolation but cannot detect near-duplicates, shared biospecimens or every form of cohort leakage.

Command-line locked inference additionally emits a SHA-256 chain-of-custody sidecar binding the model artifact and feature file to the exact prediction CSV. Command-line external evaluation refuses predictions without a valid sidecar, records the verified model and feature hashes in the metric report, and emits a second sidecar binding predictions, labels, bootstrap settings and metrics. This makes post-inference file substitution detectable. It is not a digital signature, trusted timestamp, proof of independent collection or substitute for prospective registration.

Applicability never changes the primary external-cohort denominator: headline metrics include every submitted patient. The evaluator may additionally report a locked-rule in-domain sensitivity analysis when at least four patients and two patients per class remain, together with its coverage and secondary role. It is marked non-evaluable otherwise. This is a robustness view, not permission to exclude domain-shifted cases from the primary result.

For multicentre external cohorts, pooled metrics remain primary and retain every patient. A descriptive centre table reports sample size, prevalence and metrics for each two-class centre; single-class centres stay visible as non-evaluable rather than disappearing. Across-centre AUC and Brier ranges expose gross heterogeneity but are unstable with small centres and do not establish geographic transportability.

Locked patient predictions include a 95% interval from deterministic patient-level bootstrap refits. The interval holds the fully selected regularization value fixed and therefore quantifies instability under resampling of the development patients. It is not a clinical coverage guarantee and omits distribution shift, measurement error, mechanistic misspecification and uncertainty in how single-cell features map into twin states.

## Identifiability policy

Calibration never silently fits every parameter. The caller must name the fitted parameters, measurement error (`sd`) weights every residual, and positive parameters are optimized in log space. Both layers report local Jacobian rank and condition number, search-boundary hits, and approximate log-scale intervals when rank and residual degrees of freedom permit. The identifiability warning activates for sparse observations, rank deficiency, severe conditioning or a boundary solution. An optional one-parameter profile establishes a separate Gaussian-error weighted-SSE optimum, fixes each target over a log-spaced range and re-optimizes all remaining fitted parameters under the same objective. It reports whether each asymptotic 95% support boundary was actually crossed; an unresolved scan edge is never presented as a confidence limit, and the Gaussian optimum is not conflated with the primary robust estimate. Profile likelihood is stronger than a local Hessian but still relies on the observation model and an asymptotic chi-square reference; posterior geometry and independent perturbation data remain required for a defensible calibrated twin.

The same policy now applies to the six-zone model. Longitudinal species-by-zone observations can constrain explicitly selected transport/reaction parameters, optionally from a measured spatial initial state and under a fixed systemic PD-1 occupancy trajectory. Spatial reports add numerical Jacobian rank and condition number to the observation-count warning. These local diagnostics expose obvious non-identifiability but do not replace profile likelihood, posterior predictive checking or independent perturbation experiments.

Systemic and spatial fitted-parameter JSON files can be loaded directly into simulation, spatial counterfactual and hybrid-feature generation commands. Loading requires an exact parameter schema and finite positive values. Each output receives a sidecar with the source file's SHA-256 and parameter class, while numeric `qc_` audit fields record calibrated-parameter use without entering model fitting. A successful calibration report therefore cannot be silently disconnected from the downstream prediction run.

Global sensitivity screening samples configurable positive parameter ranges on a logarithmic scale and writes both the design and the outcome ranking. A low rank correlation does not establish irrelevance: non-monotonic effects, interactions and saturation can hide influence. Parameters that dominate an outcome are candidates for focused measurement and profile-likelihood analysis rather than automatic tuning.

## Mechanism–ML boundary

Single-cell features initialize selected immune states through bounded, documented rules. For every patient, the same initial state is simulated under treated and untreated schedules, producing counterfactual deltas. Outcome labels are joined only after this transformation. This ordering prevents the most direct form of target leakage, but it does not remove cohort confounding or guarantee transportability; nested model selection and locked external validation are still required.

The mapping distinguishes absence from invalidity. Missing optional features use declared population defaults, whereas supplied non-numeric, non-finite or negative signature/proxy values and fractions outside `[0,1]` fail closed instead of being clipped or silently replaced. This protects numerical provenance but does not calibrate assay units.

Patient-specific initialization now covers both scales. Dendritic activation, TCF7-memory, cytotoxic/CD8 abundance and an optional tumour-antigen proxy initialize afferent and paracortical spatial states through bounded multipliers. Antigen is not silently inferred from transcriptomic interferon activity: without an explicit proxy the population antigen default is used. Spatial treated-minus-control endpoints are generated before outcomes are joined, preserving the same label-free boundary as the systemic features.

If matched spatial data are supplied, their patient-zone-species states replace the feature-derived spatial initialization before counterfactual generation. Exact patient-set equality is required across single-cell features and spatial states, outcome-like columns are rejected, and the measured-source audit field is excluded from model fitting. This closes the executable spatial-omics-to-response-prediction path without allowing labels into the mechanistic transform.

Baseline and tissue boundaries are explicit. CSV and AnnData ingestion preserve the same tissue, clonotype and temporal metadata, including a custom AnnData `obs` field explicitly named as the time column. Recognized temporal columns force a declared baseline value and every patient must retain selected cells; post-treatment and baseline cells are never silently pooled, and missing patient/cell-type identifiers fail during loading. When paired TDLN/tumour data are present, expression signatures and cell fractions default to the TDLN compartment while TCR sharing uses both tissues. This prevents input-format choice, tumour expression or on-treatment biology from silently changing a baseline TDLN predictor.
