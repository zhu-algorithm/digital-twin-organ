# Data contract

Raw patient data is deliberately not committed. A cell table must contain one row per cell, `patient_id`, `cell_type`, and any available signature-gene columns. Labels must contain one row per patient with `patient_id,response`, where response is binary 0/1 and is defined before analysis (for example RECIST CR/PR vs SD/PD).

The repository's existing CELLxGENE-derived organ files are useful provenance/context but do **not** contain a verified paired TDLN cohort with anti-PD-1 outcomes. They are therefore not silently reused as supervised labels.
