# Data contract

Raw patient data is deliberately not committed. A cell table must contain one row per cell, `patient_id`, `cell_type`, and any available signature-gene columns. Labels must contain one row per patient with `patient_id,response`, where response is binary 0/1 and is defined before analysis (for example RECIST CR/PR vs SD/PD).

Optional paired scRNA/TCR input columns are `clonotype_id` and `tissue`. When present, the pipeline adds clone richness, Simpson concentration, normalized entropy, top-clone fraction and the fraction of TDLN clonotypes also observed in tumour. `tissue` values are matched case-insensitively on `tdln`/`lymph` and `tumour`/`tumor`; projects should harmonize these labels before analysis.

An optional patient-level `tumour_antigen_proxy` may be supplied to `augment-twin` after featurization when it comes from a documented assay such as neoantigen burden or lymph antigen measurement. It must be a finite non-negative relative scale. The pipeline does not derive this quantity from interferon-expression signatures; when absent it preserves the population default. Hybrid output includes systemic endpoints plus treated/control efferent-effector and treated paracortical-Tpex spatial endpoints. These simulated columns are model-derived covariates, not measured biomarkers.

The repository's existing CELLxGENE-derived organ files are useful provenance/context but do **not** contain a verified paired TDLN cohort with anti-PD-1 outcomes. They are therefore not silently reused as supervised labels.

Every imported cohort should ship a JSON provenance manifest with `cohort_id`, `source`, `license`, `modalities`, `outcome_definition`, and a list of files containing relative `path`, semantic `role`, and SHA-256. `tdln-pd1 qc` verifies the hashes and fails closed on missing files, altered content, missing/duplicate patients, invalid response labels or a single-class matched cohort. Its report also records cell counts, response balance, signature coverage, paired TDLN–tumour availability and clonotype assignment.
