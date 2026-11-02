# Data contract

Raw patient data is deliberately not committed. A cell table must contain one row per cell, `patient_id`, `cell_type`, and any available signature-gene columns. Labels must contain one row per patient with `patient_id,response`, where response is binary 0/1 and is defined before analysis (for example RECIST CR/PR vs SD/PD).

Optional paired scRNA/TCR input columns are `clonotype_id` and `tissue`. When present, the pipeline adds clone richness, Simpson concentration, normalized entropy, top-clone fraction and the fraction of TDLN clonotypes also observed in tumour. `tissue` values are matched case-insensitively on `tdln`/`lymph` and `tumour`/`tumor`; projects should harmonize these labels before analysis.

The repository's existing CELLxGENE-derived organ files are useful provenance/context but do **not** contain a verified paired TDLN cohort with anti-PD-1 outcomes. They are therefore not silently reused as supervised labels.
