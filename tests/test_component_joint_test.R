# Run from r_analysis/maintenance_ra with Rscript --vanilla ../../tests/test_component_joint_test.R
suppressPackageStartupMessages(library(data.table))
source("R/event_components.R")
source("R/component_joint_test.R")
set.seed(7319)
frame <- rbindlist(lapply(seq_len(72L), function(event) {
  k <- 1L + event %% 4L
  data.table(event = event, interval = factor(interval_levels[seq_len(k)], levels = interval_levels))
}))
# Two events per surgery and two surgeries per patient; unequal interval counts.
frame[, caseid := ceiling(event / 2)]
frame[, subjectid := ceiling(caseid / 2)]
frame[, covariate := rnorm(.N)]
x <- model.matrix(~ 0 + interval + covariate, data = frame)
y <- as.numeric(x %*% c(-0.15, -0.09, 0.04, 0.08, 0.2)) +
  rnorm(max(frame$caseid), sd = 0.3)[frame$caseid] + rnorm(nrow(frame), sd = 0.2)
prepared <- list(data = frame, x = x, y = cbind(change_log_alpha = y), event = as.character(frame$event))
result <- component_alpha_joint_test(prepared)
w <- event_weights(prepared$event)
# Independent normal equations and surgery-score sandwich, without vcovCL.
bread <- solve(crossprod(x, w * x))
beta <- as.numeric(bread %*% crossprod(x, w * y))
scores <- rowsum(x * as.numeric(w * (y - x %*% beta)), group = frame$caseid)
n <- nrow(x)
k <- ncol(x)
g <- nrow(scores)
covariance <- bread %*% crossprod(scores) %*% bread * ((n - 1) / (n - k)) * (g / (g - 1))
manual_chi <- as.numeric(crossprod(beta[1:4], solve(covariance[1:4, 1:4], beta[1:4])))
stopifnot(
  max(abs(vapply(result$interval_changes, function(z) z$estimate, numeric(1)) - beta[1:4])) < 1e-12,
  max(abs(vapply(result$interval_changes, function(z) z$standard_error, numeric(1)) - sqrt(diag(covariance)[1:4]))) < 1e-12,
  abs(result$global_test_all_four_changes_zero$chi_square - manual_chi) < 1e-10,
  result$cohort$events_or_pairs == 72L, result$cohort$surgeries == 36L,
  result$cohort$subjects == 18L
)
# This fixture would detect accidentally reverting to unweighted OLS.
stopifnot(max(abs(coef(lm(y ~ 0 + x)) - beta)) > 1e-3)
cat("WLS coefficients and surgery-cluster HC1 match independent normal-equation/sandwich calculation.\n")
