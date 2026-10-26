# Purpose: add ECG-derived HRV measures to the event-centred sensitivity model.
source("R/common.R")

suppressPackageStartupMessages(library(sandwich))

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
responses <- fread(file.path(input_root, "event_responses.csv"))
intervals <- fread(file.path(input_root, "event_intervals.csv"))
hrv <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_hrv_sensitivity_v1", "hrv_windows.csv"
))[kind == "event"]
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")
baseline_covariates <- c(
  "baseline_logit_relative_alpha", "baseline_propofol_ce",
  "baseline_remifentanil_ce", "baseline_etco2", "baseline_body_temperature"
)

baseline_eligibility <- intervals[
  interval == "baseline_m5_m2",
  .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * 0.5))
]
responses <- merge(responses, baseline_eligibility, by = "event_id", all.x = TRUE, sort = FALSE)

baseline_hrv <- hrv[interval == "baseline_m5_m2", .(
  event_id = target_id,
  baseline_log_rmssd = hrv_event_log_rmssd_ms,
  baseline_qrs_hr = hrv_event_qrs_hr_bpm,
  baseline_hrv_eligible = hrv_event_eligible
)]

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

fit_clustered <- function(frame, formula) {
  fit <- lm(formula, data = frame)
  covariance <- vcovCL(fit, cluster = frame$case_cluster, type = "HC1", cadjust = TRUE)
  list(fit = fit, covariance = covariance)
}

model_summary <- function(object) {
  names <- grep("^interval", names(coef(object$fit)), value = TRUE)
  beta <- coef(object$fit)[names]
  covariance <- object$covariance[names, names, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% MASS::ginv(covariance) %*% beta)
  list(
    rows = nobs(object$fit),
    global_test = list(
      chi_square = statistic, degrees_of_freedom = length(names),
      p_value = pchisq(statistic, length(names), lower.tail = FALSE)
    ),
    intervals = setNames(lapply(names, function(name) {
      estimate <- unname(coef(object$fit)[name])
      standard_error <- sqrt(object$covariance[name, name])
      list(
        estimate = estimate, standard_error = standard_error,
        normal_95ci_lower = estimate - 1.959963984540054 * standard_error,
        normal_95ci_upper = estimate + 1.959963984540054 * standard_error
      )
    }), interval_levels)
  )
}

real <- copy(responses[event_type == "real"])
eligible <- real$primary_ra_subwindows >= ceiling(real$planned_subwindows * 0.5) &
  logical_true(real$baseline_fraction_eligible) &
  finite_rows(real, c("logit_relative_alpha", baseline_covariates))
real <- real[eligible]
real[, change_logit_ra := logit_relative_alpha - baseline_logit_relative_alpha]
real[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
real[, case_cluster := as.character(caseid)]
for (column in baseline_covariates) {
  real[, (paste0("c_", column)) := get(column) - mean(get(column))]
}
real <- merge(real, baseline_hrv, by = "event_id", all = FALSE, sort = FALSE)
real <- real[
  logical_true(baseline_hrv_eligible) & is.finite(baseline_log_rmssd) & is.finite(baseline_qrs_hr)
]
real[, c_baseline_log_rmssd := baseline_log_rmssd - mean(baseline_log_rmssd)]
real[, c_baseline_qrs_hr := baseline_qrs_hr - mean(baseline_qrs_hr)]
base_formula <- as.formula(sprintf(
  "change_logit_ra ~ 0 + interval + %s",
  paste(paste0("c_", baseline_covariates), collapse = " + ")
))
hrv_formula <- update(base_formula, . ~ . + c_baseline_qrs_hr + c_baseline_log_rmssd)
real_base <- fit_clustered(real, base_formula)
real_hrv <- fit_clustered(real, hrv_formula)

event_frame <- copy(responses)
eligible <- event_frame$primary_ra_subwindows >= ceiling(event_frame$planned_subwindows * 0.5) &
  logical_true(event_frame$baseline_fraction_eligible) &
  finite_rows(event_frame, c("logit_relative_alpha", "baseline_logit_relative_alpha"))
event_frame <- event_frame[eligible]
event_frame[, change := logit_relative_alpha - baseline_logit_relative_alpha]
keys <- c("caseid", "subjectid", "event_number", "interval")
real_change <- event_frame[event_type == "real", c(keys, "change"), with = FALSE]
pseudo_change <- event_frame[event_type == "pseudo", c(keys, "change"), with = FALSE]
setnames(real_change, "change", "real_change")
setnames(pseudo_change, "change", "pseudo_change")
paired <- merge(real_change, pseudo_change, by = keys, all = FALSE, sort = FALSE)
paired[, response := real_change - pseudo_change]

metadata <- unique(responses[, .(event_id, caseid, subjectid, event_number, event_type)])
metadata <- merge(metadata, baseline_hrv, by = "event_id", all = FALSE, sort = FALSE)
metadata <- metadata[
  logical_true(baseline_hrv_eligible) & is.finite(baseline_log_rmssd) & is.finite(baseline_qrs_hr)
]
real_meta <- metadata[event_type == "real", .(
  caseid, event_number,
  real_baseline_log_rmssd = baseline_log_rmssd,
  real_baseline_qrs_hr = baseline_qrs_hr
)]
pseudo_meta <- metadata[event_type == "pseudo", .(
  caseid, event_number,
  pseudo_baseline_log_rmssd = baseline_log_rmssd,
  pseudo_baseline_qrs_hr = baseline_qrs_hr
)]
differences <- merge(real_meta, pseudo_meta, by = c("caseid", "event_number"), all = FALSE)
differences[, baseline_log_rmssd_difference := real_baseline_log_rmssd - pseudo_baseline_log_rmssd]
differences[, baseline_qrs_hr_difference := real_baseline_qrs_hr - pseudo_baseline_qrs_hr]
paired <- merge(
  paired,
  differences[, .(caseid, event_number, baseline_log_rmssd_difference, baseline_qrs_hr_difference)],
  by = c("caseid", "event_number"), all = FALSE, sort = FALSE
)
paired[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
paired[, case_cluster := as.character(caseid)]
paired[, event_cluster := paste(caseid, event_number, sep = ":")]
paired[, c_baseline_log_rmssd_difference := baseline_log_rmssd_difference - mean(baseline_log_rmssd_difference)]
paired[, c_baseline_qrs_hr_difference := baseline_qrs_hr_difference - mean(baseline_qrs_hr_difference)]
paired_base <- fit_clustered(paired, response ~ 0 + interval)
paired_hrv <- fit_clustered(
  paired,
  response ~ 0 + interval + c_baseline_qrs_hr_difference + c_baseline_log_rmssd_difference
)

hrv_baseline <- baseline_hrv[, .(
  event_id, baseline_log_rmssd, baseline_qrs_hr, baseline_hrv_eligible
)]
hrv_response <- copy(hrv[interval %in% interval_levels])
setnames(hrv_response, "target_id", "event_id")
hrv_response <- merge(hrv_response, hrv_baseline, by = "event_id", all.x = TRUE, sort = FALSE)
hrv_response <- hrv_response[
  logical_true(hrv_event_eligible) & logical_true(baseline_hrv_eligible) &
    is.finite(hrv_event_log_rmssd_ms) & is.finite(baseline_log_rmssd)
]
hrv_response[, change_log_rmssd := hrv_event_log_rmssd_ms - baseline_log_rmssd]
hrv_response[, change_qrs_hr := hrv_event_qrs_hr_bpm - baseline_qrs_hr]
event_keys <- c("caseid", "event_number", "interval")
real_response <- hrv_response[event_type == "real", c(event_keys, "change_log_rmssd", "change_qrs_hr"), with = FALSE]
pseudo_response <- hrv_response[event_type == "pseudo", c(event_keys, "change_log_rmssd", "change_qrs_hr"), with = FALSE]
setnames(real_response, c("change_log_rmssd", "change_qrs_hr"), c("real_hrv", "real_hr"))
setnames(pseudo_response, c("change_log_rmssd", "change_qrs_hr"), c("pseudo_hrv", "pseudo_hr"))
parallel <- merge(real_response, pseudo_response, by = event_keys, all = FALSE, sort = FALSE)
parallel[, real_minus_pseudo_hrv := real_hrv - pseudo_hrv]
parallel[, real_minus_pseudo_hr := real_hr - pseudo_hr]
parallel[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
parallel[, case_cluster := as.character(caseid)]

hrv_coef <- unname(coef(real_hrv$fit)["c_baseline_log_rmssd"])
hrv_se <- sqrt(real_hrv$covariance["c_baseline_log_rmssd", "c_baseline_log_rmssd"])
results <- list(
  status = "R reanalysis of post-hoc event HRV sensitivity",
  causal_interpretation = FALSE,
  real_event_baseline_adjustment = list(
    events = uniqueN(real$event_id), surgeries = uniqueN(real$caseid), subjects = uniqueN(real$subjectid),
    same_subset_standard_adjustment = model_summary(real_base),
    plus_baseline_hr_and_hrv = model_summary(real_hrv),
    baseline_hrv_coefficient = list(
      estimate = hrv_coef, standard_error = hrv_se,
      p_value = 2 * pnorm(abs(hrv_coef / hrv_se), lower.tail = FALSE)
    )
  ),
  paired_real_vs_pseudo_baseline_adjustment = list(
    pairs = uniqueN(paired$event_cluster), surgeries = uniqueN(paired$caseid),
    same_subset_unadjusted = model_summary(paired_base),
    plus_real_minus_pseudo_baseline_hr_and_hrv = model_summary(paired_hrv)
  ),
  post_event_hrv_parallel_response = list(
    pairs = uniqueN(parallel[, .(caseid, event_number)]), surgeries = uniqueN(parallel$caseid),
    log_rmssd_real_minus_pseudo_change = model_summary(
      fit_clustered(parallel, real_minus_pseudo_hrv ~ 0 + interval)
    ),
    qrs_hr_real_minus_pseudo_change_bpm = model_summary(
      fit_clustered(parallel, real_minus_pseudo_hr ~ 0 + interval)
    )
  )
)

write_json_pretty(results, file.path(output_root, "r_event_hrv_sensitivity_results.json"))
print(results)
