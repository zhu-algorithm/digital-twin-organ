# TDLN digital-twin model card

## Intended use

Research simulation of the tumour–blood–tumour-draining-lymph-node axis, generation of falsifiable hypotheses, and evaluation of analysis pipelines. It is **not** validated for diagnosis, drug selection, dose selection, or individual clinical decisions.

## Architecture

The executable core currently contains 12 non-negative state variables spanning tumour burden and antigen, migratory mature dendritic cells, naive CD8 T cells, TCF1-like progenitor-exhausted cells (Tpex), TDLN and tumour effectors, exhausted tumour T cells, Tregs, plasma/TDLN antibody exposure, and PD-1 receptor occupancy. Bolus dosing is integrated piecewise so discontinuities are explicit rather than hidden inside the ODE solver.

The spatial ontology separates afferent lymph, subcapsular sinus, paracortex, follicle, medulla, and efferent lymph. Single-cell annotations can be mapped to these zones; unresolved labels remain unresolved instead of being guessed.

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

## Required validation ladder

1. **Unit verification:** dimensions, solver convergence, mass/flow balance, limit cases.
2. **Module calibration:** antibody PK/occupancy, antigen/DC migration, Tpex proliferation and egress, tumour killing.
3. **Multimodal data assimilation:** matched TDLN/tumour scRNA-seq, TCR clonotypes, spatial transcriptomics, imaging and clinical PK.
4. **Retrospective external validation:** locked parameters and endpoints on independent centres.
5. **Prospective validation:** preregistered protocol and uncertainty-aware predictions.

Until steps 2–5 are completed with real data, default trajectories are demonstrations of model behaviour, not estimates of a patient outcome.
