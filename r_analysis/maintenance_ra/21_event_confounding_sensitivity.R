# Purpose: adjust event and reference-epoch comparisons for time-varying covariate changes.
source("R/common.R")

suppressPackageStartupMessages({
  library(sandwich)
})

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
responses <- fread(file.path(input_root, "event_responses.csv"))
intervals <- fread(file.path(input_root, "event_intervals.csv"))
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")
interval_labels <- c("0-1 min", "1-3 min", "3-5 min", "5-10 min")

baseline <- intervals[
  interval == "baseline_m5_m2",
  .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * 0.50))
]
responses <- merge(responses, baseline, by = "event_id", all.x = TRUE, sort = FALSE)

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

center_columns <- function(frame, columns) {
  result <- character()
  for (column in columns) {
    name <- paste0("c_", column)
    frame[, (name) := get(column) - mean(get(column))]
    result <- c(result, name)
  }
  result
}

fit_clustered <- function(frame, formula) {
  fit <- lm(formula, data = frame)
  covariance <- sandwich::vcovCL(
    fit, cluster = frame$caseid, type = "HC1", cadjust = TRUE
  )
  list(fit = fit, covariance = covariance)
}

bootstrap_intervals <- function(frame, formula, coefficient_names, repetitions, seed) {
  set.seed(seed)
  model_frame <- model.frame(formula, data = frame)
  x <- model.matrix(formula, data = model_frame)
  y <- model.response(model_frame)
  positions <- match(coefficient_names, colnames(x))
  subjects <- unique(frame$subjectid)
  rows_by_subject <- split(seq_len(nrow(frame)), frame$subjectid)
  values <- matrix(NA_real_, nrow = repetitions, ncol = length(positions))
  for (b in seq_len(repetitions)) {
    sampled <- sample(subjects, length(subjects), replace = TRUE)
    rows <- unlist(rows_by_subject[as.character(sampled)], use.names = FALSE)
    current <- lm.fit(x[rows, , drop = FALSE], y[rows])
    values[b, ] <- current$coefficients[positions]
  }
  t(apply(values, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
}

summarize_fit <- function(frame, formula, label, seed) {
  fit_object <- fit_clustered(frame, formula)
  coefficient_names <- grep("^interval", names(coef(fit_object$fit)), value = TRUE)
  estimates <- unname(coef(fit_object$fit)[coefficient_names])
  bootstrap_ci <- bootstrap_intervals(frame, formula, coefficient_names, 2000L, seed)
  covariance <- fit_object$covariance[coefficient_names, coefficient_names, drop = FALSE]
  statistic <- as.numeric(t(estimates) %*% MASS::ginv(covariance) %*% estimates)
  data.table(
    Analysis = label,
    Interval = interval_labels,
    Estimate_RA_percentage_points = estimates,
    Lower_95ci = bootstrap_ci[, 1],
    Upper_95ci = bootstrap_ci[, 2],
    Global_P = pchisq(statistic, length(estimates), lower.tail = FALSE),
    Rows = nrow(frame),
    Events_or_pairs = uniqueN(frame$event_cluster),
    Surgeries = uniqueN(frame$caseid),
    Patients = uniqueN(frame$subjectid)
  )
}

prepare_real <- function() {
  variables <- c(
    "relative_alpha", "baseline_relative_alpha",
    "propofol_ce", "baseline_propofol_ce",
    "remifentanil_ce", "baseline_remifentanil_ce",
    "etco2", "baseline_etco2",
    "body_temperature", "baseline_body_temperature"
  )
  frame <- copy(responses[event_type == "real"])
  frame <- frame[
    primary_ra_subwindows >= ceiling(planned_subwindows * 0.50) &
      logical_true(baseline_fraction_eligible) & finite_rows(frame, variables)
  ]
  frame[, response := 100 * (relative_alpha - baseline_relative_alpha)]
  frame[, `:=`(
    d_propofol_ce = propofol_ce - baseline_propofol_ce,
    d_remifentanil_ce = remifentanil_ce - baseline_remifentanil_ce,
    d_etco2 = etco2 - baseline_etco2,
    d_temperature = body_temperature - baseline_body_temperature,
    interval = factor(interval, levels = interval_levels, ordered = TRUE),
    event_cluster = as.character(event_id)
  )]
  covariates <- c(
    "baseline_relative_alpha", "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature",
    "d_propofol_ce", "d_remifentanil_ce", "d_etco2", "d_temperature"
  )
  centered <- center_columns(frame, covariates)
  formula <- as.formula(sprintf("response ~ 0 + interval + %s", paste(centered, collapse = " + ")))
  list(data = frame, formula = formula)
}

prepare_paired <- function() {
  variables <- c(
    "relative_alpha", "baseline_relative_alpha", "baseline_art_map",
    "propofol_ce", "baseline_propofol_ce",
    "remifentanil_ce", "baseline_remifentanil_ce",
    "etco2", "baseline_etco2",
    "body_temperature", "baseline_body_temperature"
  )
  frame <- copy(responses)
  frame <- frame[
    primary_ra_subwindows >= ceiling(planned_subwindows * 0.50) &
      logical_true(baseline_fraction_eligible) & finite_rows(frame, variables)
  ]
  frame[, `:=`(
    response = 100 * (relative_alpha - baseline_relative_alpha),
    d_propofol_ce = propofol_ce - baseline_propofol_ce,
    d_remifentanil_ce = remifentanil_ce - baseline_remifentanil_ce,
    d_etco2 = etco2 - baseline_etco2,
    d_temperature = body_temperature - baseline_body_temperature
  )]
  keys <- c("caseid", "subjectid", "event_number", "interval")
  keep <- c(
    keys, "event_time_seconds", "response", "baseline_relative_alpha", "baseline_art_map",
    "baseline_propofol_ce", "baseline_remifentanil_ce", "baseline_etco2",
    "baseline_body_temperature", "d_propofol_ce", "d_remifentanil_ce",
    "d_etco2", "d_temperature"
  )
  real <- frame[event_type == "real", ..keep]
  pseudo <- frame[event_type == "pseudo", ..keep]
  setnames(real, setdiff(keep, keys), paste0(setdiff(keep, keys), "_real"))
  setnames(pseudo, setdiff(keep, keys), paste0(setdiff(keep, keys), "_pseudo"))
  paired <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  paired[, `:=`(
    response = response_real - response_pseudo,
    interval = factor(interval, levels = interval_levels, ordered = TRUE),
    event_cluster = paste(caseid, "pair", event_number, sep = ":"),
    time_difference_minutes = (event_time_seconds_real - event_time_seconds_pseudo) / 60
  )]
  base_names <- c(
    "baseline_relative_alpha", "baseline_art_map", "baseline_propofol_ce",
    "baseline_remifentanil_ce", "baseline_etco2", "baseline_body_temperature",
    "d_propofol_ce", "d_remifentanil_ce", "d_etco2", "d_temperature"
  )
  difference_names <- paste0("diff_", base_names)
  for (j in seq_along(base_names)) {
    name <- base_names[j]
    paired[, (difference_names[j]) := get(paste0(name, "_real")) - get(paste0(name, "_pseudo"))]
  }
  covariates <- c(difference_names, "time_difference_minutes")
  centered <- center_columns(paired, covariates)
  formula <- as.formula(sprintf("response ~ 0 + interval + %s", paste(centered, collapse = " + ")))
  list(data = paired, formula = formula)
}

real <- prepare_real()
paired <- prepare_paired()
result <- rbindlist(list(
  summarize_fit(real$data, real$formula, "Real event with time-varying covariate adjustment", 202608271L),
  summarize_fit(paired$data, paired$formula, "Real-minus-reference with covariate-change adjustment", 202608272L)
))
fwrite(result, file.path(output_root, "r_event_timevarying_confounding_estimates.csv"))

diagnostic_variables <- c(
  "event_time_seconds", "baseline_relative_alpha", "baseline_art_map",
  "baseline_propofol_ce", "baseline_remifentanil_ce", "baseline_etco2",
  "baseline_body_temperature"
)
diagnostic_labels <- c(
  "Event time, min", "Baseline RA, proportion", "Baseline MAP, mmHg",
  "Baseline propofol Ce, microgram/mL", "Baseline remifentanil Ce, ng/mL",
  "Baseline EtCO2, mmHg", "Baseline temperature, degrees C"
)
pair_rows <- paired$data[interval == interval_levels[1]]
diagnostics <- rbindlist(lapply(seq_along(diagnostic_variables), function(j) {
  variable <- diagnostic_variables[j]
  real_value <- pair_rows[[paste0(variable, "_real")]]
  pseudo_value <- pair_rows[[paste0(variable, "_pseudo")]]
  if (variable == "event_time_seconds") {
    real_value <- real_value / 60
    pseudo_value <- pseudo_value / 60
  }
  pooled_sd <- sqrt((var(real_value) + var(pseudo_value)) / 2)
  data.table(
    Characteristic = diagnostic_labels[j],
    Real = sprintf("%.3f [%.3f, %.3f]", median(real_value), quantile(real_value, 0.25), quantile(real_value, 0.75)),
    Pseudo = sprintf("%.3f [%.3f, %.3f]", median(pseudo_value), quantile(pseudo_value, 0.25), quantile(pseudo_value, 0.75)),
    Standardized_difference = ifelse(pooled_sd > 0, (mean(real_value) - mean(pseudo_value)) / pooled_sd, NA_real_),
    Pairs = nrow(pair_rows),
    Surgeries = uniqueN(pair_rows$caseid),
    Patients = uniqueN(pair_rows$subjectid)
  )
}))
fwrite(diagnostics, file.path(output_root, "r_event_pseudo_match_diagnostics.csv"))

write_json(list(
  analysis = "Sensitivity analysis additionally adjusting for concurrent changes in propofol and remifentanil effect-site concentrations, EtCO2, and temperature; the paired model also adjusts baseline and event-time differences.",
  results = result,
  matching_diagnostics = diagnostics
), file.path(output_root, "r_event_timevarying_confounding_results.json"))

message("Event confounding sensitivity analysis complete")
