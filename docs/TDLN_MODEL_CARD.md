# TDLN digital-twin model card

## Intended use

Research simulation of the tumour–blood–tumour-draining-lymph-node axis, generation of falsifiable hypotheses, and evaluation of analysis pipelines. It is **not** validated for diagnosis, drug selection, dose selection, or individual clinical decisions.

## Architecture

The executable core currently contains 12 non-negative state variables spanning tumour burden and antigen, migratory mature dendritic cells, naive CD8 T cells, TCF1-like progenitor-exhausted cells (Tpex), TDLN and tumour effectors, exhausted tumour T cells, Tregs, plasma/TDLN antibody exposure, and PD-1 receptor occupancy. Bolus dosing is integrated piecewise so discontinuities are explicit rather than hidden inside the ODE solver.

The spatial ontology separates afferent lymph, subcapsular sinus, paracortex, follicle, medulla, and efferent lymph. Single-cell annotations can be mapped to these zones; unresolved labels remain unresolved instead of being guessed.

The executable spatial model uses directed conservative transport between these zones. Antigen and migratory DCs enter through afferent lymph; antigen capture occurs in the subcapsular sinus; antigen-bearing DCs reach the paracortex; Tpex expand and differentiate there; effectors pass through medulla into efferent lymph. The implementation exposes a reaction-free mode in which every species must be conserved, separating numerical transport verification from biological source/sink assumptions.

## Evidence anchors

- Four-compartment QSP precedent (blood, peripheral tissue, tumour, TDLN) and explicit antigen/APC/T-cell/checkpoint modules: [Jusko et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11138134/) and [Milberg et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6677731/).
- Antibody exposure gradients across primary tumour and TDLN, calibrated against immuno-PET: [Li & Hung, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9036243/).
- TDLN PD-1+TCF-1+ reservoir contributing to tumour control after checkpoint treatment: [Nature Communications, 2025](https://www.nature.com/articles/s41467-025-58510-1).
- Human paired scRNA/TCR evidence that lymph-node Tpex expand after anti-PD-1 and share clonotypes with intratumoural exhausted CD8 T cells: [Nature Communications, 2026](https://www.nature.com/articles/s41467-026-70751-2).
- Anatomically based antigen redistribution and lymph-flow modelling after draining-node removal: [Distant lymph nodes compensate, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10541098/).

These sources motivate model structure. They do not validate the default numerical parameters.

## Verification completed

- State schema and finite/non-negative trajectories.
- Explicit antibody transfer to TDLN and increased receptor occupancy.
- Higher tumour effector abundance after anti-PD-1 exposure in the default scenario.
- Dose coalescing invariance for simultaneous boluses.
- Lymph-node microanatomy zone coverage and deterministic cell-zone mapping.
- Existing single-cell feature and prediction tests remain passing.
- Synthetic-truth parameter recovery for longitudinal data assimilation.
- Deterministic uncertainty propagation under a fixed random seed.
- Label-free, patient-aligned mechanistic counterfactual feature generation.
- Per-species mass conservation in the reaction-free spatial transport network.
- Directional antigen, DC and effector movement across six TDLN zones.
- TCR clone diversity, dominance and paired TDLN–tumour lineage-sharing features.
- Locked-model prediction with feature-schema enforcement, calibration metrics and patient-level bootstrap intervals.
- Reproducible global parameter screening with Latin-hypercube sampling and outcome-specific rank correlations.
- Leave-one-centre-out transportability evaluation with explicit reporting of non-evaluable sites.

## Required validation ladder

1. **Unit verification:** dimensions, solver convergence, mass/flow balance, limit cases.
2. **Module calibration:** antibody PK/occupancy, antigen/DC migration, Tpex proliferation and egress, tumour killing.
3. **Multimodal data assimilation:** matched TDLN/tumour scRNA-seq, TCR clonotypes, spatial transcriptomics, imaging and clinical PK.
4. **Retrospective external validation:** locked parameters and endpoints on independent centres.
5. **Prospective validation:** preregistered protocol and uncertainty-aware predictions.

Until steps 2–5 are completed with real data, default trajectories are demonstrations of model behaviour, not estimates of a patient outcome.

The prediction report labels internal cross-validation separately from locked external evaluation. Discrimination alone is insufficient: Brier score, log loss, calibration intercept/slope, sensitivity/specificity and bootstrap intervals are emitted together. The software cannot prove that a supplied cohort is independent; cohort provenance and analysis preregistration remain governance requirements.

## Identifiability policy

Calibration never silently fits every parameter. The caller must name the fitted parameters, measurement error (`sd`) weights every residual, positive parameters are optimized in log space, and the report raises an identifiability warning when there are fewer than five observations per fitted parameter. This warning is heuristic rather than proof of identifiability; profile likelihood, posterior geometry and independent perturbation data remain required for a defensible calibrated twin.

Global sensitivity screening samples configurable positive parameter ranges on a logarithmic scale and writes both the design and the outcome ranking. A low rank correlation does not establish irrelevance: non-monotonic effects, interactions and saturation can hide influence. Parameters that dominate an outcome are candidates for focused measurement and profile-likelihood analysis rather than automatic tuning.

## Mechanism–ML boundary

Single-cell features initialize selected immune states through bounded, documented rules. For every patient, the same initial state is simulated under treated and untreated schedules, producing counterfactual deltas. Outcome labels are joined only after this transformation. This ordering prevents the most direct form of target leakage, but it does not remove cohort confounding or guarantee transportability; nested model selection and locked external validation are still required.
