# Liver–TDLN bridge calibration example

`liver_tdln_bridge_observations.csv` is generated from the simulator with
`antigen_drainage_gain = 0.8`. It exists only to exercise parameter recovery,
file contracts and provenance verification in CI. It is not measured biology,
not a patient cohort and not evidence that the default mechanism is valid.

Every row repeats this status in its `source` column so that extracting the CSV
without this note does not silently turn synthetic values into purported data.
Replacing it for scientific work requires measured values, measurement-error
estimates and cohort-level provenance appropriate to the assay.
