# Purpose: analyse absolute alpha and actual non-alpha power across maintenance and events.
source("R/common.R")

message("[1/6] Loading maintenance-window and event data")
main_objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
bundle <- main_objects$bundle
data <- bundle$data

data[, log_alpha := log(alpha_power)]
data[, log_nonalpha := log_alpha - logit_relative_alpha]
if (any(!is.finite(data$log_alpha)) || any(!is.finite(data$log_nonalpha))) {
  stop("Non-finite alpha/non-alpha power in the primary cohort")
}
identity_error <- max(abs(data$logit_relative_alpha - data$log_alpha + data$log_nonalpha))

case_metrics <- function(y, pred0, pred1) {
  frame <- data.table(
    subjectid = data$subjectid,
    caseid = data$caseid,
    r0 = as.numeric(y) - as.numeric(pred0),
    r1 = as.numeric(y) - as.numeric(pred1)
  )
  frame[, r0 := r0 - mean(r0), by = caseid]
  frame[, r1 := r1 - mean(r1), by = caseid]
  frame[, .(q0 = mean(r0^2), q1 = mean(r1^2)), by = .(subjectid, caseid)]
}

patient_bootstrap_g <- function(metrics, repetitions, seed) {
  set.seed(seed)
  subjects <- unique(metrics$subjectid)
  subject_index <- match(metrics$subjectid, subjects)
  draws <- numeric(repetitions)
  for (b in seq_len(repetitions)) {
    weights_subject <- as.vector(rmultinom(
      1L, length(subjects), rep(1 / length(subjects), length(subjects))
    ))
    weights <- weights_subject[subject_index]
    draws[b] <- 1 - weighted.mean(metrics$q1, weights) /
      weighted.mean(metrics$q0, weights)
  }
  list(
    repetitions = repetitions,
    percentile_95ci = as.list(as.numeric(quantile(draws, c(0.025, 0.975)))),
    standard_error = sd(draws)
  )
}

analyze_maintenance <- function(y, label, seed) {
  message(sprintf("[2/6] Maintenance mixed models: %s", label))
  fit0 <- fit_car1_random_map_slope_lme(bundle, bundle$baseline_names, outcome = y)
  fit1 <- fit_car1_random_map_slope_lme(bundle, bundle$nonlinear_names, outcome = y)
  pred0 <- fixed_prediction(fit0, bundle, bundle$baseline_names)
  pred1 <- fixed_prediction(fit1, bundle, bundle$nonlinear_names)
  metrics <- case_metrics(y, pred0, pred1)
  q0 <- mean(metrics$q0)
  q1 <- mean(metrics$q1)
  statistic <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))), 0)
  degrees <- length(bundle$nonlinear_names) - length(bundle$baseline_names)
  grid <- c(60, 65, 75, 85)
  contrast <- map_contrast_matrix(bundle, grid = grid, reference = 65)
  estimate <- as.numeric(contrast %*% fixef(fit1)[bundle$nonlinear_names])
  covariance <- vcov(fit1)[bundle$nonlinear_names, bundle$nonlinear_names]
  standard_error <- sqrt(pmax(rowSums((contrast %*% covariance) * contrast), 0))
  list(
    summary = list(
      outcome = label,
      subjects = uniqueN(data$subjectid),
      surgeries = uniqueN(data$caseid),
      windows = nrow(data),
      likelihood_ratio = statistic,
      degrees_of_freedom = degrees,
      p_value = pchisq(statistic, degrees, lower.tail = FALSE),
      q0 = q0,
      q1 = q1,
      residual_reduction_fraction = 1 - q1 / q0,
      bootstrap = patient_bootstrap_g(metrics, 2000L, seed)
    ),
    curve = data.table(
      outcome = label,
      map_mmhg = grid,
      log_contrast_vs_65 = estimate,
      lower_95ci = estimate - 1.959963984540054 * standard_error,
      upper_95ci = estimate + 1.959963984540054 * standard_error,
      power_ratio_vs_65 = exp(estimate),
      ratio_lower_95ci = exp(estimate - 1.959963984540054 * standard_error),
      ratio_upper_95ci = exp(estimate + 1.959963984540054 * standard_error)
    )
  )
}

maintenance_alpha <- analyze_maintenance(data$log_alpha, "log_alpha", 2026082701L)
maintenance_nonalpha <- analyze_maintenance(data$log_nonalpha, "log_nonalpha", 2026082702L)

message("[3/6] Exact case-equal decomposition of within-surgery RA variation")
variance_frame <- data[, .(
  subjectid,
  caseid,
  r = logit_relative_alpha,
  a = log_alpha,
  n = log_nonalpha
)]
variance_frame[, `:=`(
  r = r - mean(r),
  a = a - mean(a),
  n = n - mean(n)
), by = caseid]
variance_by_case <- variance_frame[, .(
  variance_logit_ra = mean(r^2),
  alpha_covariance_component = mean(r * a),
  nonalpha_covariance_component = mean(r * (-n))
), by = .(subjectid, caseid)]
variance_by_case[, identity_error :=
  variance_logit_ra - alpha_covariance_component - nonalpha_covariance_component]

variance_point <- variance_by_case[, .(
  variance_logit_ra = mean(variance_logit_ra),
  alpha_covariance_component = mean(alpha_covariance_component),
  nonalpha_covariance_component = mean(nonalpha_covariance_component)
)]
variance_point[, `:=`(
  alpha_fraction = alpha_covariance_component / variance_logit_ra,
  nonalpha_fraction = nonalpha_covariance_component / variance_logit_ra
)]

set.seed(2026082703L)
subjects <- unique(variance_by_case$subjectid)
subject_index <- match(variance_by_case$subjectid, subjects)
variance_boot <- matrix(NA_real_, nrow = 2000L, ncol = 2L)
colnames(variance_boot) <- c("alpha_fraction", "nonalpha_fraction")
for (b in seq_len(nrow(variance_boot))) {
  weights_subject <- as.vector(rmultinom(
    1L, length(subjects), rep(1 / length(subjects), length(subjects))
  ))
  weights <- weights_subject[subject_index]
  total <- weighted.mean(variance_by_case$variance_logit_ra, weights)
  variance_boot[b, 1L] <- weighted.mean(
    variance_by_case$alpha_covariance_component, weights
  ) / total
  variance_boot[b, 2L] <- weighted.mean(
    variance_by_case$nonalpha_covariance_component, weights
  ) / total
}

variance_decomposition <- list(
  definition = paste(
    "logit(RA)=log(alpha power)-log(non-alpha 1-20 Hz power);",
    "case-equal covariance components sum exactly to within-surgery logit(RA) variance"
  ),
  point = as.list(variance_point[1]),
  patient_bootstrap_2000_percentile_95ci = list(
    alpha_fraction = as.list(as.numeric(quantile(variance_boot[, 1L], c(0.025, 0.975)))),
    nonalpha_fraction = as.list(as.numeric(quantile(variance_boot[, 2L], c(0.025, 0.975))))
  ),
  maximum_case_identity_error = max(abs(variance_by_case$identity_error))
)

message("[4/6] Preparing descriptive low-MAP event spectral outcomes")
source("R/event_components.R")
responses <- load_component_responses(file.path(repo_root, "outputs", "maintenance_ra_full_v1"))

analyze_event <- function(prepared, seed, repetitions = 2000L) {
  interval_positions <- match(paste0("interval", interval_levels), colnames(prepared$x))
  point_beta <- multi_wls(
    prepared$x, prepared$y, event_weights(prepared$event)
  )[interval_positions, , drop = FALSE]
  rownames(point_beta) <- interval_levels
  colnames(point_beta) <- c("logit_ra", "log_alpha", "log_nonalpha", "log_delta")

  set.seed(seed)
  subjects <- unique(prepared$data$subjectid)
  subject_rows <- split(seq_len(nrow(prepared$data)), prepared$data$subjectid)
  draws <- array(
    NA_real_, dim = c(repetitions, length(interval_levels), 4L),
    dimnames = list(NULL, interval_levels, colnames(point_beta))
  )
  for (b in seq_len(repetitions)) {
    sampled <- sample(subjects, length(subjects), replace = TRUE)
    rows <- unlist(subject_rows[as.character(sampled)], use.names = FALSE)
    draw_id <- rep(seq_along(sampled), lengths(subject_rows[as.character(sampled)]))
    boot_event <- paste(prepared$event[rows], draw_id, sep = "#")
    draws[b, , ] <- multi_wls(
      prepared$x[rows, , drop = FALSE], prepared$y[rows, , drop = FALSE],
      event_weights(boot_event)
    )[interval_positions, , drop = FALSE]
  }

  table <- rbindlist(lapply(seq_along(interval_levels), function(j) {
    rbindlist(lapply(seq_len(4L), function(k) {
      estimate <- point_beta[j, k]
      ci <- quantile(draws[, j, k], c(0.025, 0.975), na.rm = TRUE)
      data.table(
        interval = interval_levels[j],
        component = colnames(point_beta)[k],
        estimate_log_scale = estimate,
        lower_95ci = ci[[1L]],
        upper_95ci = ci[[2L]],
        percent_power_change = if (k == 1L) NA_real_ else 100 * (exp(estimate) - 1),
        percent_power_change_lower = if (k == 1L) NA_real_ else 100 * (exp(ci[[1L]]) - 1),
        percent_power_change_upper = if (k == 1L) NA_real_ else 100 * (exp(ci[[2L]]) - 1)
      )
    }))
  }))

  list(
    cohort = list(
      rows = nrow(prepared$data),
      events_or_pairs = uniqueN(prepared$event),
      surgeries = uniqueN(prepared$data$caseid),
      subjects = uniqueN(prepared$data$subjectid)
    ),
    estimates = table,
    bootstrap_repetitions = repetitions
  )
}

message("[5/6] Running patient-bootstrap event decomposition")
event_real <- analyze_event(prepare_real(responses), 2026082704L)
event_paired <- analyze_event(prepare_paired(responses), 2026082705L)

message("[6/6] Writing decomposition outputs")
flatten_maintenance <- function(object) {
  interval <- unlist(object$bootstrap$percentile_95ci)
  data.table(
    outcome = object$outcome,
    subjects = object$subjects,
    surgeries = object$surgeries,
    windows = object$windows,
    likelihood_ratio = object$likelihood_ratio,
    degrees_of_freedom = object$degrees_of_freedom,
    p_value = object$p_value,
    q0 = object$q0,
    q1 = object$q1,
    residual_reduction_fraction = object$residual_reduction_fraction,
    residual_reduction_lower_95ci = interval[[1L]],
    residual_reduction_upper_95ci = interval[[2L]],
    bootstrap_standard_error = object$bootstrap$standard_error
  )
}
maintenance_table <- rbindlist(list(
  flatten_maintenance(maintenance_alpha$summary),
  flatten_maintenance(maintenance_nonalpha$summary)
), fill = TRUE)
maintenance_curves <- rbindlist(list(
  maintenance_alpha$curve,
  maintenance_nonalpha$curve
))
event_table <- rbindlist(list(
  copy(event_real$estimates)[, analysis := "real_event_adjusted"],
  copy(event_paired$estimates)[, analysis := "real_minus_pseudo"]
), fill = TRUE)
setcolorder(event_table, c("analysis", setdiff(names(event_table), "analysis")))
fwrite(maintenance_table, file.path(output_root, "r_alpha_nonalpha_maintenance_results.csv"))
fwrite(maintenance_curves, file.path(output_root, "r_alpha_nonalpha_maintenance_curves.csv"))
fwrite(event_table, file.path(output_root, "r_alpha_nonalpha_event_results.csv"))
fwrite(variance_by_case, file.path(output_root, "r_alpha_nonalpha_variance_by_case.csv"))

results <- list(
  interpretation_scope = paste(
    "Algebraic and associational decomposition only;",
    "this analysis does not establish a causal effect of MAP or alpha power"
  ),
  identity = list(
    equation = "logit(RA) = log(alpha power) - log(non-alpha 1-20 Hz power)",
    maximum_window_identity_error = identity_error
  ),
  maintenance = list(
    log_alpha = maintenance_alpha$summary,
    log_nonalpha = maintenance_nonalpha$summary,
    curves = maintenance_curves
  ),
  within_surgery_variance_decomposition = variance_decomposition,
  low_map_event = list(
    nonalpha_definition = paste(
      "Mean of subwindow-level log non-alpha power, calculated as",
      "log(alpha power)-logit(RA) within each eligible 30-s subwindow"
    ),
    real_event_adjusted = event_real,
    paired_real_minus_pseudo = event_paired
  )
)
write_json_pretty(results, file.path(output_root, "r_alpha_nonalpha_decomposition_results.json"))

real_01 <- event_table[analysis == "real_event_adjusted" & interval == "response_0_1"]
paired_01 <- event_table[analysis == "real_minus_pseudo" & interval == "response_0_1"]
summary_lines <- c(
  "# RAのα・非αパワー分解解析",
  "",
  "各30秒subwindowでNを1–20 Hzの非αパワーとしてlog(N)=log(A)−logit(RA)を計算し、その後にevent interval内で平均した。",
  "interval-levelのRA、αパワー、非αパワーは別々の記述的outcomeであり、厳密な寄与率分解は行わない。",
  "",
  "## 持続低MAPイベント後0–1分",
  "",
  sprintf(
    "実イベントの調整解析では、logit(RA)変化は%.4f、αパワー変化は%.2f%%（95%%信頼区間%.2f–%.2f%%）、δパワー変化は%.2f%%（%.2f–%.2f%%）、非αパワー変化は%.2f%%（%.2f–%.2f%%）であった。",
    real_01[component == "logit_ra", estimate_log_scale],
    real_01[component == "log_alpha", percent_power_change],
    real_01[component == "log_alpha", percent_power_change_lower],
    real_01[component == "log_alpha", percent_power_change_upper],
    real_01[component == "log_delta", percent_power_change],
    real_01[component == "log_delta", percent_power_change_lower],
    real_01[component == "log_delta", percent_power_change_upper],
    real_01[component == "log_nonalpha", percent_power_change],
    real_01[component == "log_nonalpha", percent_power_change_lower],
    real_01[component == "log_nonalpha", percent_power_change_upper]
  ),
  sprintf(
    "偽イベントとの差の差では、logit(RA)変化は%.4f、αパワー変化は%.2f%%（95%%信頼区間%.2f–%.2f%%）、δパワー変化は%.2f%%（%.2f–%.2f%%）、非αパワー変化は%.2f%%（%.2f–%.2f%%）であった。",
    paired_01[component == "logit_ra", estimate_log_scale],
    paired_01[component == "log_alpha", percent_power_change],
    paired_01[component == "log_alpha", percent_power_change_lower],
    paired_01[component == "log_alpha", percent_power_change_upper],
    paired_01[component == "log_delta", percent_power_change],
    paired_01[component == "log_delta", percent_power_change_lower],
    paired_01[component == "log_delta", percent_power_change_upper],
    paired_01[component == "log_nonalpha", percent_power_change],
    paired_01[component == "log_nonalpha", percent_power_change_lower],
    paired_01[component == "log_nonalpha", percent_power_change_upper]
  ),
  "最初の1分のRA低下は主として絶対αパワー低下と同時にみられ、δパワーまたは実測非αパワーの増加を示す明確な証拠はなかった。",
  "",
  "## 解釈上の限界",
  "",
  "この結果は集団平均の低MAPイベント応答を分解したものであり、個々のRA低下が常にαパワー低下だけから生じることを意味しない。",
  "また、MAP低下がαパワー低下を因果的に生じさせたことを証明するものではない。"
)
writeLines(
  summary_lines,
  file.path(output_root, "RA_ALPHA_NONALPHA_DECOMPOSITION_JA_20260827.md"),
  useBytes = TRUE
)

print(maintenance_table)
print(maintenance_curves)
print(event_table[interval == "response_0_1"])
