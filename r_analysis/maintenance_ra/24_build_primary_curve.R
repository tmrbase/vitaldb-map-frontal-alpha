# Purpose: standardise the fitted MAP–RA curve with equal surgery weights for Figure 1.
source("R/common.R")

objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))

case_equal_weights <- function(caseid) {
  counts <- table(caseid)
  1 / (length(counts) * as.numeric(counts[as.character(caseid)]))
}

standardized_ra_curve <- function(objects, grid = seq(50, 110, by = 0.5), reference = 65) {
  bundle <- objects$bundle
  fit <- objects$nonlinear
  fixed_names <- bundle$nonlinear_names
  beta <- fixef(fit)[fixed_names]
  covariance <- vcov(fit)[fixed_names, fixed_names, drop = FALSE]
  x_observed <- as.matrix(bundle$design[, fixed_names, with = FALSE])
  map_names <- grep("^map_basis[0-9]+_within$", fixed_names, value = TRUE)
  observed_basis <- natural_spline_basis(bundle$data$art_map, bundle$knots$map)
  case_mean_basis <- sapply(seq_len(ncol(observed_basis)), function(j) {
    ave(observed_basis[, j], bundle$data$caseid, FUN = mean)
  })
  weights <- case_equal_weights(bundle$data$caseid)

  counterfactual <- function(map_value) {
    design <- x_observed
    basis <- natural_spline_basis(map_value, bundle$knots$map)[1, ]
    for (j in seq_along(map_names)) {
      name <- map_names[j]
      raw_within <- basis[j] - case_mean_basis[, j]
      design[, name] <- (
        raw_within - bundle$scaling[[name]]$mean
      ) / bundle$scaling[[name]]$standard_deviation
    }
    eta <- as.numeric(design %*% beta)
    probability <- plogis(eta)
    estimate <- sum(weights * probability)
    gradient <- colSums(design * as.numeric(weights * probability * (1 - probability)))
    list(estimate = estimate, gradient = gradient)
  }

  ref <- counterfactual(reference)
  rbindlist(lapply(grid, function(value) {
    current <- counterfactual(value)
    gradient <- current$gradient - ref$gradient
    standard_error <- sqrt(max(as.numeric(t(gradient) %*% covariance %*% gradient), 0))
    estimate <- 100 * (current$estimate - ref$estimate)
    data.table(
      map_mmhg = value,
      adjusted_ra_difference_pp = estimate,
      lower_95ci = estimate - 1.959963984540054 * 100 * standard_error,
      upper_95ci = estimate + 1.959963984540054 * 100 * standard_error,
      standardized_ra_at_map_percent = 100 * current$estimate,
      reference_ra_percent = 100 * ref$estimate
    )
  }))
}

curve <- standardized_ra_curve(objects)
fwrite(curve, file.path(output_root, "standardized_map_ra_curve.csv"))
cat("Wrote", nrow(curve), "standardised MAP-RA curve points\n")
