# Purpose: refit the maintenance model with available cardiac and oxygenation covariates.
source("R/common.R")

main_objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
windows <- fread(windows_path)
covariates <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_additional_covariates_v1",
  "additional_covariates.csv"
))
audit <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_additional_covariates_v1",
  "case_covariate_audit.csv"
))
windows <- merge(
  windows, covariates, by = c("caseid", "time_seconds"),
  all.x = TRUE, sort = FALSE
)
fixed_knots <- main_objects$bundle$knots

analyze_bundle <- function(bundle, label) {
  message(label)
  fits <- fit_standard_model_set(bundle)
  summary <- summarize_standard_model_set(bundle, fits, label)
  contrast <- map_contrast_matrix(bundle)
  summary$map_grid_mmhg <- as.list(50:90)
  summary$map_contrast_vs_65_logit_ra <- as.list(as.numeric(
    contrast %*% fixef(fits$nonlinear)[bundle$nonlinear_names]
  ))
  summary
}

build_subset <- function(data) {
  build_main_design(data, minimum_windows = 8L, fixed_knots = fixed_knots)
}

no_n2o <- audit[!logical_true(n2o_any_positive_native)]$caseid
strict_n2o <- audit[
  !logical_true(n2o_any_positive_native) &
    n2o_observed_windows >= 0.8 * scheduled_primary_windows
]$caseid

results <- list(
  status = "R reanalysis of post-hoc N2O, CO/CI, HR and SpO2 sensitivities",
  causal_interpretation = FALSE,
  n2o = list(
    operational_no_positive = analyze_bundle(
      build_subset(windows[caseid %in% no_n2o]),
      "n2o_operational_no_positive"
    ),
    strict_measured_zero = analyze_bundle(
      build_subset(windows[caseid %in% strict_n2o]),
      "n2o_strict_measured_zero"
    )
  )
)

ci_windows <- windows[logical_true(ci_qc)]
ci_reference <- build_subset(ci_windows)
ci_adjusted <- augment_bundle(ci_reference, list(c("ci", "ci")))
ci_within <- ci_adjusted$data$ci - ave(
  ci_adjusted$data$ci, ci_adjusted$data$caseid, FUN = mean
)
map_within <- ci_adjusted$data$art_map - ave(
  ci_adjusted$data$art_map, ci_adjusted$data$caseid, FUN = mean
)
results$cardiac_index <- list(
  within_ci_map_correlation = cor(ci_within, map_within),
  same_subset_unadjusted = analyze_bundle(ci_reference, "ci_subset_unadjusted"),
  ci_adjusted = analyze_bundle(ci_adjusted, "ci_adjusted")
)

ev_reference <- build_subset(ci_windows[hemodynamic_device == "EV1000"])
results$cardiac_index$ev1000_ci_adjusted <- analyze_bundle(
  augment_bundle(ev_reference, list(c("ci", "ci"))),
  "ev1000_ci_adjusted"
)

co_windows <- windows[logical_true(co_qc)]
co_reference <- build_subset(co_windows)
results$cardiac_index$cardiac_output_same_subset_unadjusted <- analyze_bundle(
  co_reference, "co_subset_unadjusted"
)
results$cardiac_index$cardiac_output_adjusted <- analyze_bundle(
  augment_bundle(co_reference, list(c("co", "co"))),
  "co_adjusted"
)

hs_windows <- windows[logical_true(hr_qc) & logical_true(spo2_qc)]
hs_reference <- build_subset(hs_windows)
results$hr_spo2 <- list(
  same_subset_unadjusted = analyze_bundle(
    hs_reference, "hr_spo2_subset_unadjusted"
  ),
  hr_spo2_adjusted = analyze_bundle(
    augment_bundle(hs_reference, list(c("hr", "hr"), c("spo2", "spo2"))),
    "hr_spo2_adjusted"
  )
)
normoxic_reference <- build_subset(hs_windows[spo2 >= 94])
results$hr_spo2$normoxic_hr_adjusted <- analyze_bundle(
  augment_bundle(normoxic_reference, list(c("hr", "hr"))),
  "normoxic_hr_adjusted"
)

write_json_pretty(results, file.path(output_root, "r_additional_covariate_results.json"))
print(lapply(results[c("n2o", "cardiac_index", "hr_spo2")], function(x) {
  lapply(x, function(y) if (is.list(y) && !is.null(y$case_equal_G)) y$case_equal_G else y)
}))
