# Purpose: refit the maintenance model with available ECG-derived HRV measures.
source("R/common.R")

main_objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
windows <- fread(windows_path)
hrv <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_hrv_sensitivity_v1", "hrv_windows.csv"
))
hrv_main <- hrv[kind == "main", .(
  caseid, time_seconds, ecg_track,
  hrv300_eligible, hrv300_log_rmssd_ms, hrv300_log_sdnn_ms,
  hrv300_qrs_hr_bpm, hrv300_valid_rr_fraction, hrv300_ecg_coverage,
  hrv60_eligible, hrv60_log_rmssd_ms, hrv60_qrs_hr_bpm,
  hrv60_valid_rr_fraction, hrv60_ecg_coverage
)]
windows <- merge(
  windows, hrv_main, by = c("caseid", "time_seconds"),
  all.x = TRUE, sort = FALSE
)
additional <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_additional_covariates_v1",
  "additional_covariates.csv"
))[, .(caseid, time_seconds, hr, hr_qc)]
windows <- merge(
  windows, additional, by = c("caseid", "time_seconds"),
  all.x = TRUE, sort = FALSE
)
fixed_knots <- main_objects$bundle$knots

fit_payload <- function(bundle, label) {
  message(label)
  fits <- fit_standard_model_set(bundle)
  list(
    summary = summarize_standard_model_set(bundle, fits, label),
    fits = fits
  )
}

coefficient_summary <- function(fit, name) {
  estimate <- unname(fixef(fit)[name])
  standard_error <- sqrt(vcov(fit)[name, name])
  list(
    estimate = estimate,
    standard_error = standard_error,
    normal_95ci_lower = estimate - 1.959963984540054 * standard_error,
    normal_95ci_upper = estimate + 1.959963984540054 * standard_error
  )
}

long_data <- windows[
  logical_true(hrv300_eligible) & logical_true(hr_qc) &
    is.finite(hrv300_log_rmssd_ms)
]
long_reference <- build_main_design(
  long_data, minimum_windows = 8L, fixed_knots = fixed_knots
)
hr_bundle <- augment_bundle(long_reference, list(c("hr", "hr")))
hrv_bundle <- augment_bundle(
  hr_bundle, list(c("hrv300_log_rmssd_ms", "log_rmssd_300"))
)
hrv_sdnn_bundle <- augment_bundle(
  hrv_bundle, list(c("hrv300_log_sdnn_ms", "log_sdnn_300"))
)

unadjusted <- fit_payload(long_reference, "hrv_same_subset_unadjusted")
hr_adjusted <- fit_payload(hr_bundle, "heart_rate_adjusted")
hrv_adjusted <- fit_payload(hrv_bundle, "heart_rate_plus_log_rmssd_300s")
sdnn_adjusted <- fit_payload(
  hrv_sdnn_bundle, "heart_rate_plus_log_rmssd_and_log_sdnn_300s"
)

pred_hr <- fixed_prediction(
  hr_adjusted$fits$baseline, hr_bundle, hr_bundle$baseline_names
)
pred_hrv <- fixed_prediction(
  hrv_adjusted$fits$baseline, hrv_bundle, hrv_bundle$baseline_names
)
hrv_block <- case_equal_within_reduction(
  hrv_bundle$data, hrv_bundle$y, pred_hr, pred_hrv
)

short_data <- windows[
  logical_true(hrv60_eligible) & logical_true(hr_qc) &
    is.finite(hrv60_log_rmssd_ms)
]
short_reference <- build_main_design(
  short_data, minimum_windows = 8L, fixed_knots = fixed_knots
)
short_bundle <- augment_bundle(
  short_reference,
  list(c("hr", "hr"), c("hrv60_log_rmssd_ms", "log_rmssd_60"))
)
short_adjusted <- fit_payload(short_bundle, "heart_rate_plus_log_rmssd_60s")

hr_difference <- abs(
  long_reference$data$hrv300_qrs_hr_bpm - long_reference$data$hr
)
results <- list(
  analysis_version = "R_RA_MAP_HRV_SENSITIVITY_v1",
  status = "R reanalysis of post-hoc HRV sensitivity",
  causal_interpretation = FALSE,
  hrv_interpretation = "autonomic-state marker; not a validated direct pain measure",
  qrs_vs_monitor_hr_absolute_difference_bpm = list(
    median = median(hr_difference),
    q90 = unname(quantile(hr_difference, 0.90)),
    q95 = unname(quantile(hr_difference, 0.95))
  ),
  same_subset_unadjusted = unadjusted$summary,
  heart_rate_adjusted = hr_adjusted$summary,
  heart_rate_plus_log_rmssd_300s_adjusted = c(
    hrv_adjusted$summary,
    list(
      hrv_block_residual_reduction_beyond_hr = list(
        case_equal_q0 = hrv_block$q0,
        case_equal_q1 = hrv_block$q1,
        case_equal_G = hrv_block$G
      ),
      log_rmssd_within_coefficient_in_map_model = coefficient_summary(
        hrv_adjusted$fits$nonlinear, "log_rmssd_300_within"
      ),
      log_rmssd_between_coefficient_in_map_model = coefficient_summary(
        hrv_adjusted$fits$nonlinear, "log_rmssd_300_between"
      )
    )
  ),
  heart_rate_plus_log_rmssd_and_log_sdnn_300s_adjusted = sdnn_adjusted$summary,
  heart_rate_plus_log_rmssd_60s_adjusted = c(
    short_adjusted$summary,
    list(
      log_rmssd_within_coefficient_in_map_model = coefficient_summary(
        short_adjusted$fits$nonlinear, "log_rmssd_60_within"
      )
    )
  )
)

write_json_pretty(results, file.path(output_root, "r_hrv_sensitivity_results.json"))
print(list(
  unadjusted_G = unadjusted$summary$case_equal_G,
  hr_adjusted_G = hr_adjusted$summary$case_equal_G,
  hrv_adjusted_G = hrv_adjusted$summary$case_equal_G,
  hrv_beyond_hr_G = hrv_block$G
))
