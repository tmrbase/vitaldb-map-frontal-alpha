# Purpose: repeat the event-centred analysis in the measured-zero volatile cohort.
source("R/common.R")

suppressPackageStartupMessages({
  library(sandwich)
})

arguments <- commandArgs(trailingOnly = TRUE)
destination <- if (length(arguments)) arguments[[1L]] else file.path(output_root, "strict_volatile_event_sensitivity.json")
input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
responses <- fread(file.path(input_root, "event_responses.csv"))
intervals <- fread(file.path(input_root, "event_intervals.csv"))
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")
fraction <- 0.50

baseline <- intervals[
  interval == "baseline_m5_m2",
  .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * fraction))
]
responses <- merge(responses, baseline, by = "event_id", all.x = TRUE, sort = FALSE)
responses <- responses[logical_true(strict_measured_zero_volatile)]

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

fit_clustered <- function(frame, formula) {
  fit <- lm(formula, data = frame)
  covariance <- vcovCL(fit, cluster = frame$case_cluster, type = "HC1", cadjust = TRUE)
  names_interval <- grep("^interval", names(coef(fit)), value = TRUE)
  beta <- coef(fit)[names_interval]
  covariance_interval <- covariance[names_interval, names_interval, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% MASS::ginv(covariance_interval) %*% beta)
  estimates <- setNames(lapply(seq_along(names_interval), function(j) {
    estimate <- unname(beta[j])
    standard_error <- sqrt(covariance_interval[j, j])
    list(
      estimate = estimate,
      standard_error = standard_error,
      normal_95ci_lower = estimate - qnorm(0.975) * standard_error,
      normal_95ci_upper = estimate + qnorm(0.975) * standard_error
    )
  }), interval_levels)
  list(
    estimates = estimates,
    global_test_all_four_zero = list(
      chi_square = statistic,
      degrees_of_freedom = length(names_interval),
      p_value = pchisq(statistic, length(names_interval), lower.tail = FALSE)
    )
  )
}

prepare_real <- function(scale) {
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  covariates <- c(
    baseline_outcome, "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  frame <- copy(responses[event_type == "real"])
  eligible <- frame$primary_ra_subwindows >= ceiling(frame$planned_subwindows * fraction) &
    logical_true(frame$baseline_fraction_eligible) & finite_rows(frame, c(outcome, covariates))
  frame <- frame[eligible]
  frame[, response := if (scale == "raw") 100 * (get(outcome) - get(baseline_outcome)) else get(outcome) - get(baseline_outcome)]
  frame[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  frame[, case_cluster := as.character(caseid)]
  centered <- character()
  for (column in covariates) {
    name <- paste0("c_", column)
    frame[, (name) := get(column) - mean(get(column))]
    centered <- c(centered, name)
  }
  formula <- as.formula(sprintf("response ~ 0 + interval + %s", paste(centered, collapse = " + ")))
  list(data = frame, formula = formula)
}

prepare_paired <- function(scale) {
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  frame <- copy(responses)
  eligible <- frame$primary_ra_subwindows >= ceiling(frame$planned_subwindows * fraction) &
    logical_true(frame$baseline_fraction_eligible) & finite_rows(frame, c(outcome, baseline_outcome))
  frame <- frame[eligible]
  frame[, change := if (scale == "raw") 100 * (get(outcome) - get(baseline_outcome)) else get(outcome) - get(baseline_outcome)]
  keys <- c("caseid", "subjectid", "event_number", "interval")
  real <- frame[event_type == "real", c(keys, "change"), with = FALSE]
  setnames(real, "change", "real_change")
  pseudo <- frame[event_type == "pseudo", c(keys, "change"), with = FALSE]
  setnames(pseudo, "change", "pseudo_change")
  paired <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  paired[, response := real_change - pseudo_change]
  paired[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  paired[, case_cluster := as.character(caseid)]
  list(data = paired, formula = response ~ 0 + interval)
}

analyse <- function(scale) {
  real <- prepare_real(scale)
  paired <- prepare_paired(scale)
  list(
    real_event_adjusted = c(
      list(rows = nrow(real$data), events = uniqueN(real$data$event_id), surgeries = uniqueN(real$data$caseid), subjects = uniqueN(real$data$subjectid)),
      fit_clustered(real$data, real$formula)
    ),
    paired_real_vs_pseudo = c(
      list(rows = nrow(paired$data), pairs = uniqueN(paste(paired$data$caseid, paired$data$event_number)), surgeries = uniqueN(paired$data$caseid), subjects = uniqueN(paired$data$subjectid)),
      fit_clustered(paired$data, paired$formula)
    )
  )
}

results <- list(
  cohort_definition = "Events from surgeries with strict_measured_zero_volatile=TRUE; interval completeness >=50%",
  strict_volatile_event_rows_available = nrow(responses),
  raw_ra_percentage_points = analyse("raw"),
  logit_ra = analyse("logit")
)
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
write_json_pretty(results, destination)
print(results)
