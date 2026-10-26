# The design and response are supplied by script 20's shared preparation.
# No formula reconstruction, additional filtering or weight substitution occurs.
component_alpha_joint_test <- function(prepared) {
  x <- prepared$x
  y <- prepared$y[, "change_log_alpha"]
  weights <- event_weights(prepared$event)
  stopifnot(
    nrow(x) == nrow(prepared$data), length(y) == nrow(x),
    length(weights) == nrow(x), all(is.finite(x)), all(is.finite(y)),
    all(is.finite(weights)), all(weights > 0),
    !anyNA(prepared$data$caseid), !anyNA(prepared$data$subjectid),
    !anyNA(prepared$event)
  )
  fit <- lm(y ~ 0 + x, weights = weights, na.action = na.fail, x = TRUE, y = TRUE)
  stopifnot(
    fit$rank == ncol(x), nobs(fit) == nrow(prepared$data),
    identical(dim(fit$x), dim(x)), identical(as.vector(fit$x), as.vector(x)),
    identical(unname(fit$y), unname(y)),
    identical(unname(fit$weights), unname(weights))
  )
  beta <- multi_wls(x, prepared$y, weights)[, "change_log_alpha"]
  maximum_difference <- max(abs(unname(coef(fit)) - unname(beta)))
  if (maximum_difference > 1e-10) stop("Alpha lm and script-20 WLS coefficients differ")
  event_totals <- tapply(fit$weights, prepared$event, sum)
  stopifnot(all(abs(event_totals - 1) < 1e-12))
  covariance <- sandwich::vcovCL(
    fit, cluster = prepared$data$caseid, type = "HC1", cadjust = TRUE
  )
  positions <- match(paste0("interval", interval_levels), colnames(x))
  stopifnot(length(positions) == 4L, !anyNA(positions))
  b <- unname(coef(fit)[positions])
  v <- covariance[positions, positions, drop = FALSE]
  if (qr(v)$rank != 4L) stop("Four-interval covariance is singular")
  statistic <- as.numeric(crossprod(b, solve(v, b)))
  list(
    cohort = list(
      rows = nrow(prepared$data), events_or_pairs = data.table::uniqueN(prepared$event),
      surgeries = data.table::uniqueN(prepared$data$caseid),
      subjects = data.table::uniqueN(prepared$data$subjectid)
    ),
    interval_changes = setNames(lapply(seq_along(b), function(j) list(
      estimate = b[j], standard_error = sqrt(v[j, j])
    )), interval_levels),
    global_test_all_four_changes_zero = list(
      chi_square = statistic, degrees_of_freedom = 4L,
      p_value = pchisq(statistic, 4L, lower.tail = FALSE)
    ),
    consistency = list(
      identical_design_response_weights = TRUE,
      maximum_coefficient_difference_from_multi_wls = maximum_difference,
      maximum_event_total_weight_error = max(abs(event_totals - 1))
    )
  )
}
