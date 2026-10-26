# Validation and input provenance

Date: September 30, 2026. The accompanying records describe execution on existing
local inputs generated from pinned PhysioNet VitalDB v1.0.0 records. They do not
constitute a new download-to-result rerun.

## Executed checks

- Python 3.12.11: all 13 tests passed; lint and compilation passed. Tests include
  patient deletion against direct GLS, signal-processing and event-definition
  invariants, and rejection of mismatched component estimates, counts and P values.
- R 4.4.1: all 33 R source files parsed. The whole-patient exclusion test passed.
  An independent synthetic-data test verified WLS coefficients, standard errors
  and the four-interval Wald statistic against normal equations and a
  surgery-score HC1 sandwich, with unequal interval counts and multiple events
  per surgery.
- Both workflow dry runs passed (9 source-to-result and 31 analysis stages).
  The source structure check and focused component numerical validator passed.
- The alpha joint-test script was executed on the local event response, interval
  and subwindow files. All four alpha coefficients and the cohort counts matched
  the saved component-analysis JSON and CSV. Fitted design values, responses and
  event weights were checked directly; each event/pair has total weight one.

## Aggregate results checked

| Alpha joint test | Rows | Events/pairs | Surgeries | Patients | Chi-square (4 df) | P |
|---|---:|---:|---:|---:|---:|---:|
| Real event | 1359 | 414 | 323 | 322 | 27.2425948889 | 0.0000177554464 |
| Real minus reference | 938 | 302 | 251 | 251 | 9.1008972197 | 0.0586264186 |

The Table SVI displays are <0.001 and 0.059. Point estimates and patient-bootstrap
CIs concern the 0–1-minute interval: -6.36% (-8.91% to -3.72%) and -5.76%
(-10.45% to -1.01%). These interval-specific CIs and the global four-interval
Wald test address different hypotheses. Existing bootstrap output was reused;
the 2000 bootstrap replicates were not recalculated for the focused checks.

The recorded patient-influence and full sensitivity-script executions on the
same date yielded 894 patients, 904 surgeries and 20,725 windows after exclusion;
overall MAP P=2.2505996875e-29 and nonlinearity P=0.0003567427. Random-slope
reliability was 0.288137152 (quartiles 0.182443069 and 0.450208444); split-half
Spearman correlations were 0.078775772 and 0.090473419. These match Tables SII/SIII
at reported precision. Those model fits were not repeated solely for packaging.

On October 1, 2026, the 0–1-minute alpha percentage displays were checked
against the saved CSV and JSON, using 100 × (exp(log-scale value) − 1).
Point estimates and both CI limits matched at the displayed precision. This
was a display check of existing aggregate output, without model or bootstrap
re-execution.

## Boundaries

`VALIDATION.json` records input/output SHA-256 hashes, aggregate results and the
checks completed for this repository snapshot. Existing R binaries reported
being built under R 4.4.3; execution completed successfully with R 4.4.1 and the
recorded package versions, including sandwich 3.1-1.

No fresh PhysioNet download, waveform re-extraction, primary MAP–RA model refit,
full bootstrap rerun or complete source-to-result execution was performed for
these focused checks. No patient-, surgery-, window- or event-level data, excluded
identifiers or model RDS objects are redistributed. The manuscript package's
Word/PDF and manifest checks are recorded separately in its preparation folder.
