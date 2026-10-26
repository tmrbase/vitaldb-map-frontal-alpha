# Statistical analysis order and report mapping

Run the complete workflow from the repository root:

```powershell
uv run python scripts/run_source_to_result.py --execute --workers 6
```

`scripts/run_analysis.py` invokes the following analyses in dependency order. Numeric filename prefixes are stable script identifiers, not a claim that omitted numbers are missing workflow steps.

| Script | Role in the final report |
|---|---|
| `01_main_model.R` | Random-intercept-only sensitivity in Table SIII |
| `18_random_slope_primary.R` | Authoritative B0/B1 models; Results, Table II, and Table SII |
| `19_random_slope_patient_bootstrap.R` | Patient-bootstrap uncertainty and random-slope diagnostics in Table SII |
| `24_build_primary_curve.R` | Surgery-equal standardised curve for Figure 1 |
| `generate_patient_influence.py` | Original patient-deletion ranking for Methods S5; local input for Table SIII |
| `03_model_sensitivities.R` | Table SII reliability/split-half diagnostics; Table SIII covariance, outcome-scale, sampling, and top-1%-patient exclusion/refit |
| `04_secondary_eeg.R` | Secondary EEG outcomes in Table SIII |
| `20_alpha_nonalpha_decomposition.R` | Maintenance and event alpha/non-alpha results in Tables SIII and SVI |
| `05_event_analysis.R` at 50% | Primary event-centred trajectory in Tables SV and SVI |
| `05_event_analysis.R` at 80% | Interval-completeness sensitivity in Table SVI |
| `06_event_spectral.R` | Same-model alpha WLS joint tests in Table SVI; shared preparation and coefficient/cohort checks against script 20 |
| `07_additional_covariates.R` | Cardiac and oxygenation covariate sensitivity in Table SIII |
| `08_hrv_sensitivity.R` | Maintenance HRV sensitivity in Table SIII |
| `09_event_map_explanation.R` | Contemporaneous event MAP coefficient in Table SV |
| `10_table1_missingness.R` | Table I, cohort/retention summaries in Table SI, and Figure S1 |
| `11_event_hrv_sensitivity.R` | Event HRV sensitivity in Table SVI |
| `21_event_confounding_sensitivity.R` | Time-varying covariate and reference-epoch analyses in Tables SV and SVI |
| `22_selection_and_pk_audit.R` | Included/excluded cohort comparison and PK agreement in Tables SI and SII |
| `23_residual_diagnostics.R` | Residual summaries and plots in Table SII and Figure S2 |
| `26_window_observability_audit.R` | Window-retention audit in Table SI and Figure S3 |
| `common_support_audit.py` | Surgery-level observed MAP support for Table SIV |
| `overlap_restricted_contrasts.R` | Support-restricted contrasts in Table SIV |
| `run_hinge65_analysis.R` | Secondary 65-mm Hg piecewise-linear models in Table SIV |
| `check_hinge65_diagnostics.R` | Hinge-model convergence and fit diagnostics |
| `hinge65_cluster_robust.R` | Patient-cluster CR1S inference for Table SIV |
| `primary_map_contrast_patient_CR1S.R` | Patient-cluster CR1S CIs for Figure 1 and Table II |
| `analyse_event_pretrend.R` | Adjusted pre-event RA intervals in Figure S4 and Table SV |
| `run_pro_event_supplement.R` | Additional event-design diagnostics in Tables SV and SVI |
| `strict_volatile_event_sensitivity.R` | Measured-zero volatile event sensitivity in Table SVI |
| `negative_control_lead_analysis.py` | Lead-to-baseline negative control in Table SV |
| `build_figures_and_tables.py` | Figure 1, Figure S4, source data, and aggregate supporting tables |

The driver uses only regenerated files under `outputs/`. Actual event non-alpha power is calculated within each eligible 30-s subwindow in the shared `R/event_components.R` used by scripts 20 and 06 before interval aggregation.

See `docs/METHODS_CODE_MAP.md` for exact patient-influence, diagnostic, and Table SVI mappings.
