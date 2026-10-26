# Purpose: fit the reported secondary EEG outcome models summarised in Table S3.
source("R/common.R")

main_objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
bundle <- main_objects$bundle
data <- bundle$data

if (any(data$alpha_power <= 0 | data$delta_power <= 0 | data$relative_alpha <= 0)) {
  stop("Secondary spectral outcomes contain non-positive values")
}

outcomes <- list(
  log_absolute_alpha = log(data$alpha_power),
  log_delta = log(data$delta_power),
  log_total_1_20 = log(data$alpha_power / data$relative_alpha),
  log_adr = data$log_adr
)
labels <- c(
  log_absolute_alpha = "Log absolute alpha power",
  log_delta = "Log delta power",
  log_total_1_20 = "Log total 1-20 Hz power",
  log_adr = "Log alpha/delta ratio"
)

case_metrics <- function(y, pred0, pred1) {
  frame <- data.table(caseid = data$caseid, subjectid = data$subjectid)
  frame[, r0 := y - pred0]
  frame[, r1 := y - pred1]
  frame[, r0 := r0 - mean(r0), by = caseid]
  frame[, r1 := r1 - mean(r1), by = caseid]
  frame[, .(q0 = mean(r0^2), q1 = mean(r1^2), windows = .N), by = .(subjectid, caseid)]
}

fixed_model_patient_bootstrap <- function(metrics, seed, repetitions = 2000L) {
  set.seed(seed)
  subjects <- unique(metrics$subjectid)
  subject_index <- match(metrics$subjectid, subjects)
  draws <- numeric(repetitions)
  for (b in seq_len(repetitions)) {
    weights_subject <- as.vector(rmultinom(1L, length(subjects), rep(1 / length(subjects), length(subjects))))
    weights <- weights_subject[subject_index]
    draws[b] <- 1 - weighted.mean(metrics$q1, weights) / weighted.mean(metrics$q0, weights)
  }
  c(
    lower = unname(quantile(draws, 0.025)),
    upper = unname(quantile(draws, 0.975)),
    bootstrap_standard_error = sd(draws)
  )
}

results <- list()
curves <- list()
selected_map <- c(60, 65, 75, 85)
for (i in seq_along(outcomes)) {
  outcome_name <- names(outcomes)[i]
  y <- outcomes[[i]]
  message(sprintf("[%d/%d] %s", i, length(outcomes), outcome_name))
  fit0 <- fit_car1_random_map_slope_lme(bundle, bundle$baseline_names, outcome = y)
  fit1 <- fit_car1_random_map_slope_lme(bundle, bundle$nonlinear_names, outcome = y)
  pred0 <- fixed_prediction(fit0, bundle, bundle$baseline_names)
  pred1 <- fixed_prediction(fit1, bundle, bundle$nonlinear_names)
  metrics <- case_metrics(y, pred0, pred1)
  q0 <- mean(metrics$q0)
  q1 <- mean(metrics$q1)
  interval <- fixed_model_patient_bootstrap(metrics, 2026082403L + i - 1L)
  statistic <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))), 0)
  degrees <- length(bundle$nonlinear_names) - length(bundle$baseline_names)
  contrast <- map_contrast_matrix(bundle, grid = selected_map, reference = 65)
  estimate <- as.numeric(contrast %*% fixef(fit1)[bundle$nonlinear_names])
  covariance <- vcov(fit1)[bundle$nonlinear_names, bundle$nonlinear_names]
  standard_error <- sqrt(pmax(rowSums((contrast %*% covariance) * contrast), 0))
  results[[outcome_name]] <- list(
    outcome = outcome_name,
    label = labels[[outcome_name]],
    subjects = uniqueN(data$subjectid),
    surgeries = uniqueN(data$caseid),
    windows = nrow(data),
    likelihood_ratio = statistic,
    df = degrees,
    raw_p = pchisq(statistic, degrees, lower.tail = FALSE),
    q0 = q0,
    q1 = q1,
    residual_reduction_fraction = 1 - q1 / q0,
    residual_reduction_lower = interval[["lower"]],
    residual_reduction_upper = interval[["upper"]],
    bootstrap_standard_error = interval[["bootstrap_standard_error"]]
  )
  curves[[outcome_name]] <- data.table(
    outcome = outcome_name,
    map_mmhg = selected_map,
    log_contrast_vs_65 = estimate,
    lower = estimate - 1.96 * standard_error,
    upper = estimate + 1.96 * standard_error,
    ratio_vs_65 = exp(estimate),
    ratio_lower = exp(estimate - 1.96 * standard_error),
    ratio_upper = exp(estimate + 1.96 * standard_error)
  )
}

raw_p <- vapply(results, function(x) x$raw_p, numeric(1))
holm_p <- p.adjust(raw_p, method = "holm")
for (name in names(results)) results[[name]]$holm_p <- unname(holm_p[[name]])

result_table <- rbindlist(lapply(results, as.data.table), fill = TRUE)
curve_table <- rbindlist(curves)
fwrite(result_table, file.path(output_root, "r_secondary_eeg_global_results.csv"))
fwrite(curve_table, file.path(output_root, "r_secondary_eeg_selected_curves.csv"))
write_json_pretty(results, file.path(output_root, "r_secondary_eeg_results.json"))
print(result_table)
print(curve_table)
