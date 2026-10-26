# Purpose: estimate event-centred RA trajectories at the requested interval-completeness threshold.
source("R/common.R")

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
responses <- fread(file.path(input_root, "event_responses.csv"))
intervals <- fread(file.path(input_root, "event_intervals.csv"))
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")
arguments <- commandArgs(trailingOnly = TRUE)
strict_fraction <- if (length(arguments)) as.numeric(arguments[[1L]]) else 0.50
if (!is.finite(strict_fraction) || strict_fraction <= 0 || strict_fraction > 1) {
  stop("The interval-availability fraction must lie in (0, 1]")
}

baseline <- intervals[
  interval == "baseline_m5_m2",
  .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * strict_fraction))
]
responses <- merge(responses, baseline, by = "event_id", all.x = TRUE, sort = FALSE)

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

robust_fit <- function(frame, formula) {
  fit <- lm(formula, data = frame)
  covariance <- vcovCL(
    fit, cluster = frame$case_cluster, type = "HC1", cadjust = TRUE
  )
  list(fit = fit, covariance = covariance)
}

wald_zero <- function(fit_object, coefficient_names) {
  beta <- coef(fit_object$fit)[coefficient_names]
  covariance <- fit_object$covariance[coefficient_names, coefficient_names, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% MASS::ginv(covariance) %*% beta)
  list(
    chi_square = statistic,
    degrees_of_freedom = length(coefficient_names),
    p_value = pchisq(statistic, length(coefficient_names), lower.tail = FALSE)
  )
}

coefficient_summary <- function(fit_object, coefficient_names) {
  setNames(lapply(coefficient_names, function(name) {
    estimate <- unname(coef(fit_object$fit)[name])
    standard_error <- sqrt(fit_object$covariance[name, name])
    list(
      estimate = estimate,
      standard_error = standard_error,
      normal_95ci_lower = estimate - 1.959963984540054 * standard_error,
      normal_95ci_upper = estimate + 1.959963984540054 * standard_error
    )
  }), interval_levels)
}

prepare_real <- function(scale = c("raw", "logit")) {
  scale <- match.arg(scale)
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  covariates <- c(
    baseline_outcome, "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  frame <- copy(responses[event_type == "real"])
  eligible <- frame$primary_ra_subwindows >= ceiling(frame$planned_subwindows * strict_fraction) &
    logical_true(frame$baseline_fraction_eligible) & finite_rows(frame, c(outcome, covariates))
  frame <- frame[eligible]
  frame[, response := if (scale == "raw") {
    100 * (get(outcome) - get(baseline_outcome))
  } else {
    get(outcome) - get(baseline_outcome)
  }]
  frame[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  frame[, case_cluster := as.character(caseid)]
  frame[, event_cluster := as.character(event_id)]
  centered <- character()
  for (column in covariates) {
    name <- paste0("c_", column)
    frame[, (name) := get(column) - mean(get(column))]
    centered <- c(centered, name)
  }
  formula <- as.formula(
    sprintf("response ~ 0 + interval + %s", paste(centered, collapse = " + "))
  )
  list(data = frame, formula = formula)
}

prepare_paired <- function(scale = c("raw", "logit")) {
  scale <- match.arg(scale)
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  frame <- copy(responses)
  eligible <- frame$primary_ra_subwindows >= ceiling(frame$planned_subwindows * strict_fraction) &
    logical_true(frame$baseline_fraction_eligible) & finite_rows(frame, c(outcome, baseline_outcome))
  frame <- frame[eligible]
  frame[, change := if (scale == "raw") {
    100 * (get(outcome) - get(baseline_outcome))
  } else {
    get(outcome) - get(baseline_outcome)
  }]
  keys <- c("caseid", "subjectid", "event_number", "interval")
  real <- frame[event_type == "real", c(keys, "change"), with = FALSE]
  setnames(real, "change", "real_change")
  pseudo <- frame[event_type == "pseudo", c(keys, "change"), with = FALSE]
  setnames(pseudo, "change", "pseudo_change")
  paired <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  paired[, response := real_change - pseudo_change]
  paired[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  paired[, case_cluster := as.character(caseid)]
  paired[, event_cluster := paste(caseid, "pair", event_number, sep = ":")]
  list(data = paired, formula = response ~ 0 + interval)
}

point_analysis <- function(prepared) {
  fit <- robust_fit(prepared$data, prepared$formula)
  interval_names <- grep("^interval", names(coef(fit$fit)), value = TRUE)
  list(
    fit = fit,
    interval_names = interval_names,
    summary = list(
      formula = paste(deparse(prepared$formula), collapse = " "),
      interval_estimates = coefficient_summary(fit, interval_names),
      global_test_all_four_zero = wald_zero(fit, interval_names)
    )
  )
}

bootstrap_coefficients <- function(prepared, point, repetitions, seed) {
  set.seed(seed)
  model_frame <- model.frame(prepared$formula, data = prepared$data)
  x <- model.matrix(prepared$formula, data = model_frame)
  y <- model.response(model_frame)
  interval_positions <- match(point$interval_names, colnames(x))
  subjects <- unique(prepared$data$subjectid)
  subject_rows <- split(seq_len(nrow(prepared$data)), prepared$data$subjectid)
  values <- matrix(NA_real_, nrow = repetitions, ncol = length(interval_positions))
  for (b in seq_len(repetitions)) {
    sampled <- sample(subjects, length(subjects), replace = TRUE)
    rows <- unlist(subject_rows[as.character(sampled)], use.names = FALSE)
    fit <- lm.fit(x[rows, , drop = FALSE], y[rows])
    values[b, ] <- fit$coefficients[interval_positions]
  }
  colnames(values) <- interval_levels
  lapply(seq_along(interval_levels), function(j) {
    as.numeric(quantile(values[, j], c(0.025, 0.975), na.rm = TRUE))
  }) |> setNames(interval_levels)
}

analyze_scale <- function(scale, seed) {
  message(sprintf("Event analysis on %s scale", scale))
  real <- prepare_real(scale)
  paired <- prepare_paired(scale)
  real_point <- point_analysis(real)
  paired_point <- point_analysis(paired)
  list(
    real_event_adjusted = c(
      list(
        rows = nrow(real$data),
        events = uniqueN(real$data$event_cluster),
        surgeries = uniqueN(real$data$caseid),
        subjects = uniqueN(real$data$subjectid),
        estimation_engine = "ordinary least squares with surgery-cluster HC1 covariance",
        patient_bootstrap_percentile_95ci = bootstrap_coefficients(
          real, real_point, 2000L, seed
        )
      ),
      real_point$summary
    ),
    paired_real_vs_pseudo = c(
      list(
        rows = nrow(paired$data),
        pairs = uniqueN(paired$data$event_cluster),
        surgeries = uniqueN(paired$data$caseid),
        subjects = uniqueN(paired$data$subjectid),
        estimation_engine = "ordinary least squares with surgery-cluster HC1 covariance",
        patient_bootstrap_percentile_95ci = bootstrap_coefficients(
          paired, paired_point, 2000L, seed + 10000L
        )
      ),
      paired_point$summary
    )
  )
}

results <- list(
  implementation = list(
    language = R.version.string,
    point_estimator = "stats::lm",
    covariance = "sandwich::vcovCL, surgery-cluster HC1",
    bootstrap = "patient resampling, 2000 repetitions"
  ),
  causal_interpretation = FALSE,
  strict_fraction = strict_fraction,
  raw_ra_percentage_points = analyze_scale("raw", 20260825L),
  logit_ra = analyze_scale("logit", 20270825L)
)

output_name <- if (abs(strict_fraction - 0.50) < 1e-12) {
  "r_event_analysis_results.json"
} else {
  sprintf("r_event_analysis_results_fraction_%.2f.json", strict_fraction)
}
write_json_pretty(results, file.path(output_root, output_name))
print(results)
