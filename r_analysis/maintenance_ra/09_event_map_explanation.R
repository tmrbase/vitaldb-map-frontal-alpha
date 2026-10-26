# Purpose: estimate the descriptive contemporaneous MAP term in the event-centred RA model.
source("R/common.R")

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")
responses <- fread(file.path(input_root, "event_responses.csv"))
intervals <- fread(file.path(input_root, "event_intervals.csv"))
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")
time_varying <- c("propofol_ce", "remifentanil_ce", "etco2", "body_temperature")

attach_eligibility <- function(fraction) {
  baseline <- intervals[
    interval == "baseline_m5_m2",
    .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * fraction))
  ]
  merge(responses, baseline, by = "event_id", all.x = TRUE, sort = FALSE)
}

finite_rows <- function(frame, columns) {
  rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
}

prepare_component <- function(
  component = c("logit_ra", "log_alpha", "log_nonalpha", "log_delta"),
  fraction = 0.50,
  suppression_free = FALSE
) {
  component <- match.arg(component)
  frame <- attach_eligibility(fraction)
  if (!all(c("log_nonalpha_power", "baseline_log_nonalpha_power") %in% names(frame))) {
    stop("Event inputs must contain subwindow-derived actual non-alpha power")
  }
  mapping <- list(
    logit_ra = c("logit_relative_alpha", "baseline_logit_relative_alpha"),
    log_alpha = c("log_alpha_power", "baseline_log_alpha_power"),
    log_nonalpha = c("log_nonalpha_power", "baseline_log_nonalpha_power"),
    log_delta = c("log_delta_power", "baseline_log_delta_power")
  )[[component]]
  needed <- c(
    mapping, "art_map", "baseline_art_map", time_varying,
    paste0("baseline_", time_varying)
  )
  eligible <- frame$event_type == "real" &
    frame$primary_ra_subwindows >= ceiling(frame$planned_subwindows * fraction) &
    logical_true(frame$baseline_fraction_eligible) & finite_rows(frame, needed)
  if (suppression_free) {
    eligible <- eligible & !logical_true(frame$suppression_present) &
      !logical_true(frame$baseline_suppression_present)
  }
  frame <- frame[eligible]
  frame[, outcome := get(mapping[1])]
  frame[, baseline_outcome := get(mapping[2])]
  frame[, delta_outcome := outcome - baseline_outcome]
  frame[, delta_art_map := art_map - baseline_art_map]
  for (name in time_varying) {
    frame[, (paste0("delta_", name)) := get(name) - get(paste0("baseline_", name))]
  }
  frame[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  frame[, event_cluster := as.character(event_id)]
  frame[, case_cluster := as.character(caseid)]
  frame[, subject_cluster := as.character(subjectid)]
  setorder(frame, subjectid, caseid, event_number, interval)
  frame
}

center <- function(x) as.numeric(x) - mean(as.numeric(x))

build_design <- function(frame) {
  interval_codes <- as.integer(frame$interval)
  dummies <- sapply(seq_along(interval_levels), function(j) as.numeric(interval_codes == j))
  colnames(dummies) <- paste0("interval_", interval_levels)
  baseline_columns <- c(
    "baseline_outcome", "baseline_art_map", "baseline_propofol_ce",
    "baseline_remifentanil_ce", "baseline_etco2", "baseline_body_temperature",
    paste0("delta_", time_varying)
  )
  centered <- lapply(baseline_columns, function(name) center(frame[[name]]))
  usable <- vapply(centered, sd, numeric(1)) > 1e-12
  centered_matrix <- do.call(cbind, centered[usable])
  colnames(centered_matrix) <- paste0("centered_", baseline_columns[usable])
  x0 <- cbind(dummies, centered_matrix)
  delta_map10 <- center(frame$delta_art_map) / 10
  map_terms <- dummies * delta_map10
  colnames(map_terms) <- paste0("delta_map_per10_x_", interval_levels)
  list(x0 = x0, x1 = cbind(x0, map_terms), map_positions = seq.int(ncol(x0) + 1L, ncol(x0) + 4L))
}

event_weights <- function(event) {
  counts <- table(event)
  1 / as.numeric(counts[as.character(event)])
}

weighted_fit <- function(y, x, weights, clusters) {
  # Closed-form WLS keeps the 2,000 patient-bootstrap refits tractable. The
  # coefficient solution is the standard weighted least-squares estimator;
  # the covariance is the surgery-cluster HC1 sandwich used for the joint test.
  root <- sqrt(weights)
  beta <- qr.solve(x * root, y * root, tol = 1e-10)
  residual <- as.numeric(y - x %*% beta)
  bread <- MASS::ginv(crossprod(x, x * weights))
  score_rows <- x * as.numeric(weights * residual)
  cluster_scores <- rowsum(score_rows, clusters, reorder = FALSE)
  correction <- (nrow(cluster_scores) / (nrow(cluster_scores) - 1)) *
    ((length(y) - 1) / (length(y) - ncol(x)))
  covariance <- bread %*% crossprod(cluster_scores) %*% bread * correction
  list(beta = beta, covariance = covariance, fitted = as.numeric(x %*% beta))
}

event_equal_mse <- function(residual, event) {
  mean(data.table(event = event, value = residual^2)[, mean(value), by = event]$V1)
}

within_event_mse <- function(residual, event) {
  frame <- data.table(event = event, residual = residual)
  frame <- frame[, if (.N >= 2L) .(value = (residual - mean(residual))^2) else NULL, by = event]
  mean(frame[, mean(value), by = event]$V1)
}

metrics <- function(y, pred0, pred1, event, interval) {
  r0 <- y - pred0
  r1 <- y - pred1
  q0 <- event_equal_mse(r0, event)
  q1 <- event_equal_mse(r1, event)
  w0 <- within_event_mse(r0, event)
  w1 <- within_event_mse(r1, event)
  by_interval <- setNames(lapply(interval_levels, function(name) {
    keep <- as.character(interval) == name
    iq0 <- mean(r0[keep]^2)
    iq1 <- mean(r1[keep]^2)
    list(rows = sum(keep), q0 = iq0, q1 = iq1, residual_reduction_fraction = 1 - iq1 / iq0)
  }), interval_levels)
  list(
    event_equal_q0 = q0,
    event_equal_q1 = q1,
    event_equal_residual_reduction_fraction = 1 - q1 / q0,
    within_event_q0 = w0,
    within_event_q1 = w1,
    within_event_residual_reduction_fraction = 1 - w1 / w0,
    window_equal_residual_reduction_fraction = 1 - mean(r1^2) / mean(r0^2),
    by_interval = by_interval
  )
}

analyze_component <- function(frame, repetitions, seed, label) {
  design <- build_design(frame)
  y <- frame$delta_outcome
  weights <- event_weights(frame$event_cluster)
  fit0 <- weighted_fit(y, design$x0, weights, frame$case_cluster)
  fit1 <- weighted_fit(y, design$x1, weights, frame$case_cluster)
  apparent <- metrics(
    y, fit0$fitted, fit1$fitted, frame$event_cluster, frame$interval
  )
  positions <- design$map_positions
  beta <- fit1$beta[positions]
  covariance <- fit1$covariance[positions, positions, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% MASS::ginv(covariance) %*% beta)

  set.seed(seed)
  subjects <- unique(frame$subject_cluster)
  subject_rows <- split(seq_len(nrow(frame)), frame$subject_cluster)
  train <- numeric(repetitions)
  test <- numeric(repetitions)
  for (b in seq_len(repetitions)) {
    sampled <- sample(subjects, length(subjects), replace = TRUE)
    rows <- unlist(subject_rows[sampled], use.names = FALSE)
    draw_id <- rep(seq_along(sampled), lengths(subject_rows[sampled]))
    boot_event <- paste(frame$event_cluster[rows], draw_id, sep = "#")
    boot_weights <- event_weights(boot_event)
    b0 <- weighted_fit(
      y[rows], design$x0[rows, , drop = FALSE], boot_weights,
      paste(frame$case_cluster[rows], draw_id, sep = "#")
    )$beta
    b1 <- weighted_fit(
      y[rows], design$x1[rows, , drop = FALSE], boot_weights,
      paste(frame$case_cluster[rows], draw_id, sep = "#")
    )$beta
    training <- metrics(
      y[rows], as.numeric(design$x0[rows, , drop = FALSE] %*% b0),
      as.numeric(design$x1[rows, , drop = FALSE] %*% b1),
      boot_event, frame$interval[rows]
    )
    testing <- metrics(
      y, as.numeric(design$x0 %*% b0), as.numeric(design$x1 %*% b1),
      frame$event_cluster, frame$interval
    )
    train[b] <- training$event_equal_residual_reduction_fraction
    test[b] <- testing$event_equal_residual_reduction_fraction
  }
  apparent_g <- apparent$event_equal_residual_reduction_fraction
  training_quantiles <- quantile(train, c(0.025, 0.975))
  list(
    label = label,
    rows = nrow(frame),
    events = uniqueN(frame$event_cluster),
    surgeries = uniqueN(frame$caseid),
    subjects = uniqueN(frame$subjectid),
    map_global_test = list(
      chi_square = statistic, degrees_of_freedom = 4L,
      p_value = pchisq(statistic, 4L, lower.tail = FALSE)
    ),
    residual_explanation = apparent,
    patient_bootstrap = list(
      requested_repetitions = repetitions,
      apparent_event_equal_fraction = apparent_g,
      percentile_95ci_event_equal_fraction = as.list(as.numeric(training_quantiles)),
      basic_95ci_event_equal_fraction = as.list(c(
        2 * apparent_g - training_quantiles[2],
        2 * apparent_g - training_quantiles[1]
      )),
      mean_optimism = mean(train - test),
      optimism_corrected_event_equal_fraction = apparent_g - mean(train - test)
    )
  )
}

results <- list(
  status = "R reanalysis of event-response MAP additional explanation",
  causal_interpretation = FALSE,
  components = list()
)
for (i in seq_along(c("logit_ra", "log_alpha", "log_nonalpha", "log_delta"))) {
  name <- c("logit_ra", "log_alpha", "log_nonalpha", "log_delta")[i]
  message(sprintf("[%d/4] %s", i, name))
  results$components[[name]] <- analyze_component(
    prepare_component(name), 2000L, 2026082501L + i - 1L, name
  )
}
results$sensitivities <- list(
  logit_ra_80pct = analyze_component(
    prepare_component("logit_ra", fraction = 0.80),
    500L, 2026182501L, "logit_ra_80pct"
  ),
  logit_ra_suppression_free = analyze_component(
    prepare_component("logit_ra", suppression_free = TRUE),
    500L, 2026282501L, "logit_ra_suppression_free"
  )
)

write_json_pretty(results, file.path(output_root, "r_event_map_explanation_results.json"))
print(lapply(results$components, function(x) list(
  G = x$residual_explanation$event_equal_residual_reduction_fraction,
  p = x$map_global_test$p_value,
  optimism_corrected_G = x$patient_bootstrap$optimism_corrected_event_equal_fraction
)))
