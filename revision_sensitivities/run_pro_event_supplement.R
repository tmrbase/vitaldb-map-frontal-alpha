# Purpose: assemble the additional event-centred diagnostics reported in Tables S5 and S6.
suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
  stop(paste(
    "Usage: Rscript run_pro_event_supplement.R PROJECT_ROOT",
    "PRETREND_INTERVALS PRETREND_SUBWINDOWS OUTPUT_DIR"
  ))
}
project_root <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
pretrend_path <- normalizePath(arguments[[2L]], winslash = "/", mustWork = TRUE)
pretrend_subwindow_path <- normalizePath(arguments[[3L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(arguments[[4L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

setwd(file.path(project_root, "r_analysis", "maintenance_ra"))
source("R/common.R")

intervals <- fread(file.path(
  project_root, "outputs", "maintenance_ra_full_v1", "event_intervals.csv"
))
pretrend <- fread(pretrend_path)
pretrend_subwindows <- fread(pretrend_subwindow_path)
pre_levels <- c("pretrend_m2_m1", "pretrend_m1_0")

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

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
pretrend[, interval_available :=
  primary_ra_subwindows >= ceiling(planned_subwindows * 0.50)]
pretrend[, baseline_available :=
  baseline_primary_ra_subwindows >= ceiling(baseline_planned_subwindows * 0.50)]

fit_robust <- function(frame, formula, coefficient_names) {
  fit <- lm(formula, data = frame)
  covariance <- vcovCL(fit, cluster = frame$case_cluster, type = "HC1", cadjust = TRUE)
  beta <- coef(fit)[coefficient_names]
  local_covariance <- covariance[coefficient_names, coefficient_names, drop = FALSE]
  standard_error <- sqrt(diag(local_covariance))
  wald <- as.numeric(t(beta) %*% solve(local_covariance, beta))
  list(
    fit = fit,
    covariance = covariance,
    estimates = data.table(
      interval = pre_levels,
      estimate_percentage_points = as.numeric(beta),
      standard_error = standard_error,
      lower_95ci = as.numeric(beta) - 1.959963984540054 * standard_error,
      upper_95ci = as.numeric(beta) + 1.959963984540054 * standard_error
    ),
    global = data.table(
      statistic = wald,
      degrees_of_freedom = length(coefficient_names),
      p_value = pchisq(wald, length(coefficient_names), lower.tail = FALSE)
    )
  )
}

prepare_real_pretrend <- function(frame) {
  covariates <- c(
    "baseline_relative_alpha", "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  current <- copy(frame[event_type == "real"])
  current <- current[
    interval_available & baseline_available &
      finite_rows(current, c("relative_alpha", covariates))
  ]
  current[, response := 100 * (relative_alpha - baseline_relative_alpha)]
  current[, interval := factor(interval, levels = pre_levels, ordered = TRUE)]
  current[, case_cluster := as.character(caseid)]
  centered <- character()
  for (column in covariates) {
    name <- paste0("c_", column)
    current[, (name) := get(column) - mean(get(column))]
    centered <- c(centered, name)
  }
  formula <- as.formula(sprintf(
    "response ~ 0 + interval + %s", paste(centered, collapse = " + ")
  ))
  list(data = current, formula = formula)
}

prepare_paired_pretrend <- function(frame) {
  current <- copy(frame[
    interval_available & baseline_available &
      finite_rows(frame, c("relative_alpha", "baseline_relative_alpha"))
  ])
  current[, change := 100 * (relative_alpha - baseline_relative_alpha)]
  keys <- c("caseid", "subjectid", "event_number", "interval")
  real <- current[event_type == "real", c(keys, "change"), with = FALSE]
  pseudo <- current[event_type == "pseudo", c(keys, "change"), with = FALSE]
  setnames(real, "change", "real_change")
  setnames(pseudo, "change", "pseudo_change")
  paired <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  paired[, response := real_change - pseudo_change]
  paired[, interval := factor(interval, levels = pre_levels, ordered = TRUE)]
  paired[, case_cluster := as.character(caseid)]
  list(data = paired, formula = response ~ 0 + interval)
}

pretrend_results <- list()
for (cohort in c("broad", "strict_measured_zero_volatile")) {
  frame <- if (cohort == "broad") pretrend else pretrend[
    logical_true(strict_measured_zero_volatile)
  ]
  real <- prepare_real_pretrend(frame)
  paired <- prepare_paired_pretrend(frame)
  coefficient_names <- paste0("interval", pre_levels)
  real_fit <- fit_robust(real$data, real$formula, coefficient_names)
  paired_fit <- fit_robust(paired$data, paired$formula, coefficient_names)
  real_fit$estimates[, `:=`(
    cohort = cohort, analysis = "real_event_adjusted",
    rows = nrow(real$data), events_or_pairs = uniqueN(real$data$event_id),
    surgeries = uniqueN(real$data$caseid)
  )]
  paired_fit$estimates[, `:=`(
    cohort = cohort, analysis = "real_minus_pseudo",
    rows = nrow(paired$data),
    events_or_pairs = uniqueN(paste(paired$data$caseid, paired$data$event_number)),
    surgeries = uniqueN(paired$data$caseid)
  )]
  real_fit$global[, `:=`(
    cohort = cohort, analysis = "real_event_adjusted",
    rows = nrow(real$data), events_or_pairs = uniqueN(real$data$event_id),
    surgeries = uniqueN(real$data$caseid)
  )]
  paired_fit$global[, `:=`(
    cohort = cohort, analysis = "real_minus_pseudo",
    rows = nrow(paired$data),
    events_or_pairs = uniqueN(paste(paired$data$caseid, paired$data$event_number)),
    surgeries = uniqueN(paired$data$caseid)
  )]
  pretrend_results[[cohort]] <- list(real = real_fit, paired = paired_fit)
}

pretrend_estimates <- rbindlist(lapply(pretrend_results, function(item) {
  rbindlist(list(item$real$estimates, item$paired$estimates), fill = TRUE)
}), fill = TRUE)
pretrend_global <- rbindlist(lapply(pretrend_results, function(item) {
  rbindlist(list(item$real$global, item$paired$global), fill = TRUE)
}), fill = TRUE)
setcolorder(pretrend_estimates, c(
  "cohort", "analysis", "interval", "estimate_percentage_points",
  "lower_95ci", "upper_95ci", "standard_error", "rows",
  "events_or_pairs", "surgeries"
))
setcolorder(pretrend_global, c(
  "cohort", "analysis", "statistic", "degrees_of_freedom", "p_value",
  "rows", "events_or_pairs", "surgeries"
))

# Secondary event-centred model: event fixed effects, event-time fixed effects,
# within-event contemporaneous MAP, and within-event measured time-varying covariates.
event_levels <- c(
  "baseline_m5_m2", "pretrend_m2_m1", "pretrend_m1_0",
  "response_0_1", "response_1_3", "response_3_5", "response_5_10"
)
event_frame <- rbindlist(list(
  intervals[event_type == "real" & interval %chin% event_levels],
  pretrend[event_type == "real" & interval %chin% event_levels,
    names(intervals), with = FALSE]
), use.names = TRUE, fill = TRUE)
event_frame <- unique(event_frame, by = c("event_id", "interval"))
event_frame[, available :=
  primary_ra_subwindows >= ceiling(planned_subwindows * 0.50)]
required <- c(
  "logit_relative_alpha", "art_map", "propofol_ce", "remifentanil_ce",
  "etco2", "body_temperature"
)
event_frame <- event_frame[
  available & finite_rows(event_frame, required)
]
event_counts <- event_frame[, .(intervals = uniqueN(interval)), by = event_id]
event_frame <- event_frame[event_id %in% event_counts[intervals >= 4, event_id]]
event_frame[, interval := factor(interval, levels = event_levels, ordered = TRUE)]

dummy <- model.matrix(~ 0 + interval, data = event_frame)
colnames(dummy) <- sub("^interval", "time_", colnames(dummy))
continuous_names <- c(
  "logit_relative_alpha", "art_map", "propofol_ce", "remifentanil_ce",
  "etco2", "body_temperature"
)
continuous <- as.data.table(event_frame[, ..continuous_names])
continuous[, event_id := event_frame$event_id]
for (name in continuous_names) {
  continuous[, (name) := get(name) - mean(get(name)), by = event_id]
}
dummy_frame <- as.data.table(dummy)
dummy_frame[, event_id := event_frame$event_id]
dummy_names <- setdiff(names(dummy_frame), "event_id")
for (name in dummy_names) {
  dummy_frame[, (name) := get(name) - mean(get(name)), by = event_id]
}
analysis_frame <- cbind(
  data.table(
    y = continuous$logit_relative_alpha,
    case_cluster = as.character(event_frame$caseid),
    event_id = event_frame$event_id
  ),
  dummy_frame[, setdiff(dummy_names, "time_baseline_m5_m2"), with = FALSE],
  data.table(
    map_within_event_10 = continuous$art_map / 10,
    propofol_within_event = continuous$propofol_ce,
    remifentanil_within_event = continuous$remifentanil_ce,
    etco2_within_event = continuous$etco2,
    temperature_within_event = continuous$body_temperature
  )
)
time_names <- grep("^time_", names(analysis_frame), value = TRUE)
covariate_names <- c(
  "propofol_within_event", "remifentanil_within_event",
  "etco2_within_event", "temperature_within_event"
)
base_formula <- as.formula(sprintf(
  "y ~ 0 + %s", paste(c(time_names, covariate_names), collapse = " + ")
))
map_formula <- as.formula(sprintf(
  "y ~ 0 + %s", paste(c(time_names, "map_within_event_10", covariate_names), collapse = " + ")
))
base_fit <- lm(base_formula, data = analysis_frame)
map_fit <- lm(map_formula, data = analysis_frame)
map_covariance <- vcovCL(
  map_fit, cluster = analysis_frame$case_cluster, type = "HC1", cadjust = TRUE
)
map_beta <- coef(map_fit)["map_within_event_10"]
map_se <- sqrt(map_covariance["map_within_event_10", "map_within_event_10"])

other_matrix <- model.matrix(base_fit)
map_vector <- analysis_frame$map_within_event_10
map_residual <- lm.fit(other_matrix, map_vector)$residuals
map_r2 <- 1 - sum(map_residual^2) / sum(map_vector^2)
map_vif <- 1 / (1 - map_r2)

base_time <- coef(base_fit)[time_names]
map_time <- coef(map_fit)[time_names]
time_attenuation <- data.table(
  interval = sub("^time_", "", time_names),
  coefficient_without_map = as.numeric(base_time),
  coefficient_with_map = as.numeric(map_time),
  change_after_map_adjustment = as.numeric(map_time - base_time)
)
event_model_summary <- data.table(
  analysis = "event_fixed_effects_plus_event_time_and_contemporaneous_map",
  rows = nrow(analysis_frame),
  events = uniqueN(analysis_frame$event_id),
  surgeries = uniqueN(event_frame$caseid),
  subjects = uniqueN(event_frame$subjectid),
  map_beta_per_10_mmhg_higher = map_beta,
  map_standard_error = map_se,
  map_lower_95ci = map_beta - 1.959963984540054 * map_se,
  map_upper_95ci = map_beta + 1.959963984540054 * map_se,
  map_p_value = 2 * pnorm(abs(map_beta / map_se), lower.tail = FALSE),
  map_vif_conditional_on_event_time_and_covariates = map_vif,
  base_AIC = AIC(base_fit),
  map_AIC = AIC(map_fit)
)

# Time-alignment audit based on the saved 30-s subwindows and the extraction code.
pre0 <- pretrend_subwindows[
  event_type == "real" & relative_center_seconds >= -60 &
    relative_center_seconds < 0
]
alignment_by_event <- pre0[, .(
  minimum_30s_median_map = min(art_map, na.rm = TRUE),
  any_30s_median_map_below65 = any(art_map < 65, na.rm = TRUE),
  analysed_subwindows = .N
), by = event_id]
alignment_audit <- list(
  event_time_definition = paste(
    "left edge of the first 10-s median MAP bin in a run of at least six",
    "consecutive bins below 65 mm Hg"
  ),
  pretrend_windows = "[-120,-90), [-90,-60), [-60,-30), and [-30,0) seconds",
  eeg_window_crosses_event_time = FALSE,
  filter_scope = paste(
    "zero-phase filtering is applied separately inside each 30-s subwindow;",
    "no EEG samples at or after event time enter a pre-event window"
  ),
  preceding_isolated_low_bins_possible = TRUE,
  real_events_with_saved_pre0_subwindows = nrow(alignment_by_event),
  events_with_any_pre0_30s_median_map_below65 = sum(
    alignment_by_event$any_30s_median_map_below65
  ),
  percent_with_any_pre0_30s_median_map_below65 = 100 * mean(
    alignment_by_event$any_30s_median_map_below65
  ),
  conclusion = paste(
    "Use 'the final minute before the sustained low-MAP event start point (t=0)';",
    "do not claim 'before MAP first crossed 65 mm Hg'."
  )
)

fwrite(pretrend_estimates, file.path(output_dir, "pretrend_interval_estimates.csv"))
fwrite(pretrend_global, file.path(output_dir, "pretrend_global_tests.csv"))
fwrite(event_model_summary, file.path(output_dir, "event_centred_map_model.csv"))
fwrite(time_attenuation, file.path(output_dir, "event_time_coefficients_with_without_map.csv"))
write_json_pretty(alignment_audit, file.path(output_dir, "event_time_alignment_audit.json"))
write_json_pretty(
  list(
    interpretation = paste(
      "Associational model only. Adjustment for contemporaneous MAP",
      "must not be interpreted as mediation or a causal effect."
    ),
    pretrend_global_tests = pretrend_global,
    event_centred_model = event_model_summary,
    event_time_coefficients = time_attenuation,
    time_alignment_audit = alignment_audit
  ),
  file.path(output_dir, "pro_event_supplement_results.json")
)

print(pretrend_global)
print(event_model_summary)
print(time_attenuation)
print(alignment_audit)
