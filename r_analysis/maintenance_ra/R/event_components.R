# Shared, unchanged script-20 event eligibility, design, responses and WLS.
# Both bootstrap inference (20) and alpha joint tests (06) use these helpers.
interval_levels <- c("response_0_1", "response_1_3", "response_3_5", "response_5_10")

load_component_responses <- function(input_root) {
  responses <- fread(file.path(input_root, "event_responses.csv"))
  intervals <- fread(file.path(input_root, "event_intervals.csv"))
  subwindows <- fread(file.path(input_root, "event_subwindows.csv"))
  baseline <- intervals[
    interval == "baseline_m5_m2",
    .(event_id, baseline_fraction_eligible = primary_ra_subwindows >= ceiling(planned_subwindows * 0.50))
  ]
  responses <- merge(responses, baseline, by = "event_id", all.x = TRUE, sort = FALSE)

  # Calculate actual non-alpha power at the 30-s subwindow level and only then
  # average it within each event interval. This deliberately does not impose an
  # exact interval-level identity involving logit(mean(RA)).
  subwindows[, log_nonalpha_power := log_alpha_power - logit_relative_alpha]
  actual_nonalpha <- subwindows[
    logical_true(primary_ra_eligible) & is.finite(log_nonalpha_power),
    .(
      log_nonalpha_power = mean(log_nonalpha_power),
      actual_nonalpha_subwindows = .N
    ),
    by = .(event_id, interval)
  ]
  count_audit <- merge(
    intervals[, .(event_id, interval, primary_ra_subwindows)],
    actual_nonalpha,
    by = c("event_id", "interval"),
    all.x = TRUE,
    sort = FALSE
  )
  count_audit[is.na(actual_nonalpha_subwindows), actual_nonalpha_subwindows := 0L]
  if (any(count_audit$primary_ra_subwindows != count_audit$actual_nonalpha_subwindows)) {
    stop("Actual non-alpha subwindow counts do not match primary RA subwindow counts")
  }
  # Current source-to-result exports already contain these interval aggregates.
  # Remove them before merging the independently recomputed values so data.table
  # does not rename the columns to .x/.y and hide the canonical column names.
  existing_nonalpha_columns <- intersect(
    c("log_nonalpha_power", "baseline_log_nonalpha_power"),
    names(responses)
  )
  if (length(existing_nonalpha_columns)) {
    responses[, (existing_nonalpha_columns) := NULL]
  }
  response_nonalpha <- actual_nonalpha[interval %in% interval_levels]
  responses <- merge(
    responses,
    response_nonalpha,
    by = c("event_id", "interval"),
    all.x = TRUE,
    sort = FALSE
  )
  baseline_nonalpha <- actual_nonalpha[
    interval == "baseline_m5_m2",
    .(event_id, baseline_log_nonalpha_power = log_nonalpha_power)
  ]
  responses <- merge(responses, baseline_nonalpha, by = "event_id", all.x = TRUE, sort = FALSE)
  responses[, `:=`(
    change_logit_ra = logit_relative_alpha - baseline_logit_relative_alpha,
    change_log_alpha = log_alpha_power - baseline_log_alpha_power,
    change_log_nonalpha = log_nonalpha_power - baseline_log_nonalpha_power,
    change_log_delta = log_delta_power - baseline_log_delta_power
  )]

  finite_rows <- function(frame, columns) {
    rowSums(!is.finite(as.matrix(frame[, columns, with = FALSE]))) == 0
  }

  required <- c(
    "change_logit_ra", "change_log_alpha", "change_log_nonalpha", "change_log_delta",
    "baseline_log_alpha_power", "baseline_log_nonalpha_power", "baseline_log_delta_power",
    "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  eligible <- responses$primary_ra_subwindows >= ceiling(responses$planned_subwindows * 0.50) &
    logical_true(responses$baseline_fraction_eligible) & finite_rows(responses, required)
  responses <- responses[eligible]

  responses
}

prepare_real <- function(responses) {
  frame <- copy(responses[event_type == "real"])
  frame[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  common_covariates <- c(
    "baseline_log_alpha_power", "baseline_log_nonalpha_power", "baseline_log_delta_power",
    "baseline_propofol_ce", "baseline_remifentanil_ce",
    "baseline_etco2", "baseline_body_temperature"
  )
  centered <- character()
  for (column in common_covariates) {
    name <- paste0("c_", column)
    frame[, (name) := get(column) - mean(get(column))]
    centered <- c(centered, name)
  }
  x <- model.matrix(
    as.formula(sprintf("~ 0 + interval + %s", paste(centered, collapse = " + "))),
    data = frame
  )
  list(
    data = frame,
    x = x,
    y = as.matrix(frame[, .(
      change_logit_ra, change_log_alpha, change_log_nonalpha, change_log_delta
    )]),
    event = as.character(frame$event_id)
  )
}

prepare_paired <- function(responses) {
  keys <- c("caseid", "subjectid", "event_number", "interval")
  changes <- c(
    "change_logit_ra", "change_log_alpha", "change_log_nonalpha", "change_log_delta"
  )
  real <- responses[event_type == "real", c(keys, changes), with = FALSE]
  pseudo <- responses[event_type == "pseudo", c(keys, changes), with = FALSE]
  setnames(real, changes, paste0("real_", changes))
  setnames(pseudo, changes, paste0("pseudo_", changes))
  frame <- merge(real, pseudo, by = keys, all = FALSE, sort = FALSE)
  for (name in changes) {
    frame[, (name) := get(paste0("real_", name)) - get(paste0("pseudo_", name))]
  }
  frame[, interval := factor(interval, levels = interval_levels, ordered = TRUE)]
  x <- model.matrix(~ 0 + interval, data = frame)
  list(
    data = frame,
    x = x,
    y = as.matrix(frame[, ..changes]),
    event = paste(frame$caseid, frame$event_number, sep = ":")
  )
}

event_weights <- function(event) {
  counts <- table(event)
  1 / as.numeric(counts[event])
}

multi_wls <- function(x, y, weights) {
  # Standard multivariate WLS, evaluated directly because the same design is
  # refitted 2,000 times inside the patient bootstrap.
  root <- sqrt(weights)
  qr.solve(x * root, y * root, tol = 1e-10)
}
