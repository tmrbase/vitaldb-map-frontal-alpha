# Implementation and manuscript result map

## Patient-influence path (Methods S5; Table SIII)

1. `scripts/generate_patient_influence.py` reads the regenerated
   `outputs/maintenance_ra_full_v1/full_windows.parquet`.
2. `src/vitaldb_eeg/patient_influence.py` builds the original diagnostic design and
   fits its random-intercept CAR(1) model. Function AST fingerprints are recorded
   in `docs/influence_implementation_provenance.json`.
3. For each patient, all surgery contributions to the GLS information matrix and
   score are removed. Covariance parameters, design and knots remain fixed. The
   ranking score is the maximum absolute change in the logit(RA) contrast against
   MAP 65 over MAP 50–90 mmHg in 1-mmHg increments. This is a conditional deletion
   diagnostic, not a full covariance refit per deleted patient.
4. The patient ranking is written **locally only** to
   `outputs/maintenance_ra_model_defense_v1/subject_deletion_influence.csv`.
5. `03_model_sensitivities.R` validates that this ranking contains exactly the
   primary cohort's patients, excludes `ceiling(0.01 * N)` patients and all their
   surgeries, and rebuilds the design using the original primary spline knots.
   It refits the primary random-intercept/random-MAP-slope CAR(1) models in R.
6. `r_model_sensitivity_results.json` and `r_model_sensitivity_comparison.csv`
   contain `top_one_percent_influence_removed`. The Table SIII target is
   **894 patients, 904 surgeries and 20,725 windows**. The manuscript label is
   “Top 1% influential patients removed.”

The Python ranking retains the diagnostic basis to reproduce the
selection actually used. The primary model and exclusion refit use `nlme` and
`splines::ns` in R. No custom Python model replaces the primary inferential analysis.

## Existing random-slope diagnostics (Methods S5; Table SII)

`03_model_sensitivities.R` retains the original calculations:

- Reliability is `max(0, min(1, 1 - conditional_variance / fitted_slope_variance))`.
  Conditional variances use fitted random-effect and CAR(1) covariance parameters.
- Windows are ordered by time within surgery and split by count at `floor(n/2)`,
  requiring at least four windows in each half.
- Slopes are predicted from each half using the full-model fixed effects and
  covariance parameters and original whole-surgery MAP centering. Models are not
  separately fitted in each half.
- Correlations across surgeries use Spearman's method. The restricted diagnostic
  requires an observed MAP range of at least 10 mmHg in both halves.

The reported targets are median reliability 0.288 (quartiles 0.182 and 0.450),
split-half correlation 0.079, and restricted correlation 0.090. Generated individual
diagnostic CSV files contain surgery-level data and must remain local.

## Event sensitivity matrix (Methods S4/S7; Table SVI)

| Table SVI rows | Estimates and patient-bootstrap CIs | Four-interval joint test |
|---|---|---|
| Time-varying covariates | `21_event_confounding_sensitivity.R` | Same script, surgery-cluster HC1 Wald |
| Strict measured-zero volatile | `revision_sensitivities/strict_volatile_event_sensitivity.R` | Same script |
| 80% completeness | `05_event_analysis.R 0.80` | Same script |
| logit(RA) | `05_event_analysis.R 0.50` | Same script |
| Absolute alpha, actual non-alpha, delta | `20_alpha_nonalpha_decomposition.R` | Alpha Global P: script 06 using the same shared `R/event_components.R` preparation |

The two unadjusted four-interval alpha P values printed in Table SVI are <0.001
and 0.059. Script 06 consumes exactly script 20's prepared data, design,
alpha response and inverse-interval-count weights. It uses surgery-cluster HC1
covariance (`cadjust=TRUE`) and verifies all four coefficients and cohort counts
against script 20's saved bootstrap output. Eligibility is common to all four
component outcomes; the real-event design includes baseline log alpha, log
non-alpha and log delta plus drug, end-tidal carbon dioxide and temperature covariates.
No Holm-adjusted or delta P values are reported. Blank P cells remain unreported.

## Reporting and filenames

`revision_sensitivities/build_figures_and_tables.py` generates Figure 1 using
circles for all surgeries and squares for range-spanning surgeries, explicitly
identifies the same fitted model, and exports a 600-dpi PNG. Output filenames containing `BJA` are retained as stable workflow identifiers;
the displayed figure content corresponds to the accompanying manuscript.

The CR1S covariance estimator is
`clubSandwich::vcovCR(type="CR1S")`.

## Validation boundary and data policy

`scripts/run_analysis.py` generates the patient ranking before running the
exclusion/refit and includes the spectral joint tests. The full command is
documented in the root README. `scripts/validate_release.py --diagnostics`
checks the focused outputs against fixed published-precision targets.

Source records, patient/surgery/window/event data, excluded identifiers, model RDS
objects and generated ranking files are not included in Git or the distributed
code archive. Input reuse and the actual executed checks are
recorded separately in `docs/VALIDATION.md`. A focused rerun on existing local
inputs must not be described as a fresh download-to-result validation.
