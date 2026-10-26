# Purpose: estimate adjusted RA changes in the two intervals immediately before event t=0.
source("R/common.R")

suppressPackageStartupMessages({
  library(sandwich)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L) stop("Usage: analyse_event_pretrend.R <pretrend_csv> <output_json>")
pretrend_path <- arguments[[1L]]
destination <- arguments[[2L]]
input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
pretrend <- fread(pretrend_path)
intervals <- fread(file.path(input_root, "event_intervals.csv"))
interval_levels <- c("pretrend_m2_m1", "pretrend_m1_0")
fraction <- 0.50

baseline_columns <- c(
  "event_id", "relative_alpha", "logit_relative_alpha", "propofol_ce",
  "remifentanil_ce", "etco2", "body_temperature", "art_map",
  "planned_subwindows", "primary_ra_subwindows"
)
baseline <- intervals[interval == "baseline_m5_m2", ..baseline_columns]
setnames(
  baseline,
  setdiff(baseline_columns, "event_id"),
  paste0("baseline_", setdiff(baseline_columns, "event_id"))
)
pretrend <- merge(pretrend, baseline, by = "event_id", all.x = TRUE, sort = FALSE)
pretrend[, interval_available := primary_ra_subwindows >= ceiling(planned_subwindows * fraction)]
pretrend[, baseline_available := baseline_primary_ra_subwindows >= ceiling(baseline_planned_subwindows * fraction)]

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

fit_clustered <- function(frame, formula) {
  fit <- lm(formula, data = frame)
  covariance <- vcovCL(fit, cluster = frame$case_cluster, type = "HC1", cadjust = TRUE)
  names_interval <- grep("^interval", names(coef(fit)), value = TRUE)
  estimates <- setNames(lapply(seq_along(names_interval), function(j) {
    estimate <- unname(coef(fit)[names_interval[j]])
    standard_error <- sqrt(covariance[names_interval[j], names_interval[j]])
    list(
      estimate = estimate,
      standard_error = standard_error,
      normal_95ci_lower = estimate - qnorm(0.975) * standard_error,
      normal_95ci_upper = estimate + qnorm(0.975) * standard_error
    )
  }), interval_levels)
  contrast <- c(-1, 1)
  contrast_estimate <- as.numeric(crossprod(contrast, coef(fit)[names_interval]))
  contrast_se <- sqrt(as.numeric(t(contrast) %*% covariance[names_interval, names_interval, drop = FALSE] %*% contrast))
  list(
    estimates = estimates,
    immediate_pretrend_change_last_minus_first = list(
      estimate = contrast_estimate,
      standard_error = contrast_se,
      normal_95ci_lower = contrast_estimate - qnorm(0.975) * contrast_se,
      normal_95ci_upper = contrast_estimate + qnorm(0.975) * contrast_se,
      p_value = 2 * pnorm(abs(contrast_estimate / contrast_se), lower.tail = FALSE)
    )
  )
}

prepare_real <- function(frame, scale) {
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  covariates <- c(
    baseline_outcome, "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  current <- copy(frame[event_type == "real"])
  current <- current[logical_true(interval_available) & logical_true(baseline_available) & finite_rows(current, c(outcome, covariates))]
  current[, response := if (scale == "raw") 100 * (get(outcome) - get(baseline_outcome)) else get(outcome) - get(baseline_outcome)]
  current[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  current[, case_cluster := as.character(caseid)]
  centered <- character()
  for (column in covariates) {
    name <- paste0("c_", column)
    current[, (name) := get(column) - mean(get(column))]
    centered <- c(centered, name)
  }
  formula <- as.formula(sprintf("response ~ 0 + interval + %s", paste(centered, collapse = " + ")))
  list(data = current, formula = formula)
}

prepare_paired <- function(frame, scale) {
  baseline_outcome <- if (scale == "raw") "baseline_relative_alpha" else "baseline_logit_relative_alpha"
  outcome <- if (scale == "raw") "relative_alpha" else "logit_relative_alpha"
  current <- copy(frame)
  current <- current[logical_true(interval_available) & logical_true(baseline_available) & finite_rows(current, c(outcome, baseline_outcome))]
  current[, change := if (scale == "raw") 100 * (get(outcome) - get(baseline_outcome)) else get(outcome) - get(baseline_outcome)]
  keys <- c("caseid", "subjectid", "event_number", "interval")
  real <- current[event_type == "real", c(keys, "change"), with = FALSE]
  setnames(real, "change", "real_change")
  pseudo <- current[event_type == "pseudo", c(keys, "change"), with = FALSE]
  setnames(pseudo, "change", "pseudo_change")
  paired <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  paired[, response := real_change - pseudo_change]
  paired[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  paired[, case_cluster := as.character(caseid)]
  list(data = paired, formula = response ~ 0 + interval)
}

analyse_frame <- function(frame, scale) {
  real <- prepare_real(frame, scale)
  paired <- prepare_paired(frame, scale)
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

analyse_cohort <- function(frame) {
  list(
    raw_ra_percentage_points = analyse_frame(frame, "raw"),
    logit_ra = analyse_frame(frame, "logit")
  )
}

results <- list(
  interval_definition = list(pretrend_m2_m1 = c(-120, -60), pretrend_m1_0 = c(-60, 0)),
  baseline_definition = c(-300, -120),
  interpretation = "Negative values are lower than the -5 to -2 min baseline. The last-minus-first contrast tests change across the immediate two-minute pre-event period.",
  broad_no_positive_recorded_volatile = analyse_cohort(pretrend),
  strict_measured_zero_volatile = analyse_cohort(pretrend[logical_true(strict_measured_zero_volatile)])
)
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
write_json_pretty(results, destination)
print(results)
