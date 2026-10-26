# Table SVI alpha joint tests: exactly the WLS model used for script-20 CIs.
# Run after script 20. Existing local outputs may be used for a focused check;
# this script does not rerun maintenance mixed models or bootstrap intervals.
source("R/common.R")
source("R/event_components.R")
source("R/component_joint_test.R")

responses <- load_component_responses(file.path(repo_root, "outputs", "maintenance_ra_full_v1"))
prepared <- list(
  real_event = prepare_real(responses),
  paired_real_vs_pseudo = prepare_paired(responses)
)
reference <- jsonlite::read_json(file.path(
  output_root, "r_alpha_nonalpha_decomposition_results.json"
), simplifyVector = TRUE)$low_map_event
reference_names <- c(
  real_event = "real_event_adjusted", paired_real_vs_pseudo = "paired_real_minus_pseudo"
)
alpha <- lapply(names(prepared), function(name) {
  result <- component_alpha_joint_test(prepared[[name]])
  saved <- reference[[reference_names[[name]]]]
  stopifnot(identical(names(result$cohort), names(saved$cohort)))
  if (!all(unlist(result$cohort) == unlist(saved$cohort))) {
    stop("Analysis cohort differs from script-20 bootstrap output: ", name)
  }
  saved_alpha <- saved$estimates[saved$estimates$component == "log_alpha", ]
  stopifnot(nrow(saved_alpha) == 4L, !anyDuplicated(saved_alpha$interval))
  saved_beta <- saved_alpha$estimate_log_scale[match(interval_levels, saved_alpha$interval)]
  current_beta <- vapply(result$interval_changes, function(x) x$estimate, numeric(1))
  difference <- max(abs(current_beta - saved_beta))
  # Allow small serialization/platform rounding without changing the model.
  if (!is.finite(difference) || difference > 1e-8) {
    stop("Alpha estimates differ from script-20 bootstrap output: ", name)
  }
  result$consistency$script20_cohort_match <- TRUE
  result$consistency$maximum_difference_from_saved_script20_estimates <- difference
  result
})
names(alpha) <- names(prepared)
results <- list(
  implementation = list(
    language = R.version.string,
    estimator = "WLS; inverse number of available intervals per event/pair",
    preparation = "R/event_components.R shared with 20_alpha_nonalpha_decomposition.R",
    covariance = "sandwich::vcovCL, surgery-cluster HC1, cadjust=TRUE",
    test = "Unadjusted joint chi-square Wald test of all four interval coefficients",
    ci_source = "Unchanged script-20 patient-bootstrap intervals; no CI calculated here"
  ),
  event_spectral = list(outcomes = list(absolute_alpha = alpha))
)
write_json_pretty(results, file.path(output_root, "r_event_spectral_results.json"))
print(alpha)
