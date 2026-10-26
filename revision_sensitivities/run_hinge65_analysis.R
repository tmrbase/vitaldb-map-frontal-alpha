# Purpose: fit the secondary continuous piecewise-linear models with a knot at 65 mm Hg.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript run_hinge65_analysis.R PROJECT_ROOT OUTPUT_DIR")
}
project_root <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(arguments[[2L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

setwd(file.path(project_root, "r_analysis", "maintenance_ra"))
source("R/common.R")

threshold <- 65
main_objects <- readRDS(file.path(
  project_root,
  "outputs", "maintenance_ra_r_reanalysis_v1",
  "r_random_slope_primary_objects.rds"
))
full_bundle <- main_objects$bundle
windows <- fread(file.path(
  project_root, "outputs", "maintenance_ra_full_v1", "full_windows.csv"
))

inverse_logit <- function(x) 1 / (1 + exp(-x))

build_segmented_bundle <- function(bundle, threshold = 65) {
  result <- bundle
  map <- as.numeric(result$data$art_map)
  hinge <- pmax(threshold - map, 0)
  map_mean <- ave(map, result$data$caseid, FUN = mean)
  hinge_mean <- ave(hinge, result$data$caseid, FUN = mean)

  result$design$map_mean_10 <- (map_mean - threshold) / 10
  result$design$hinge65_mean_10 <- hinge_mean / 10
  result$design$map_within_10_raw <- (map - map_mean) / 10
  result$design$hinge65_within_10_raw <- (hinge - hinge_mean) / 10

  existing <- names(bundle$design)
  adjustment_names <- existing[
    !startsWith(existing, "map_basis") & existing != "map_linear_within_10"
  ]
  result$segmented_between_names <- c(
    adjustment_names, "map_mean_10", "hinge65_mean_10"
  )
  result$segmented_linear_names <- c(
    result$segmented_between_names, "map_within_10_raw"
  )
  result$segmented_hinge_names <- c(
    result$segmented_linear_names, "hinge65_within_10_raw"
  )
  result
}

linear_combination <- function(fit, weights, label) {
  coefficients <- fixef(fit)
  covariance <- vcov(fit)
  vector <- setNames(rep(0, length(coefficients)), names(coefficients))
  vector[names(weights)] <- weights
  estimate <- sum(vector * coefficients)
  standard_error <- sqrt(as.numeric(t(vector) %*% covariance %*% vector))
  data.table(
    estimand = label,
    estimate_log_scale = estimate,
    standard_error = standard_error,
    lower_95ci_log_scale = estimate - 1.959963984540054 * standard_error,
    upper_95ci_log_scale = estimate + 1.959963984540054 * standard_error,
    p_value = 2 * pnorm(abs(estimate / standard_error), lower.tail = FALSE)
  )
}

standardised_levels <- function(fit, segmented, targets, outcome) {
  fixed_names <- segmented$segmented_hinge_names
  coefficients <- fixef(fit)[fixed_names]
  covariance <- vcov(fit)[fixed_names, fixed_names, drop = FALSE]
  base <- as.matrix(segmented$design[, fixed_names, with = FALSE])
  map_mean <- ave(
    as.numeric(segmented$data$art_map), segmented$data$caseid, FUN = mean
  )
  observed_hinge <- pmax(threshold - as.numeric(segmented$data$art_map), 0)
  hinge_mean <- ave(observed_hinge, segmented$data$caseid, FUN = mean)

  one_level <- function(target) {
    design <- base
    design[, "map_within_10_raw"] <- (target - map_mean) / 10
    design[, "hinge65_within_10_raw"] <- (
      pmax(threshold - target, 0) - hinge_mean
    ) / 10
    eta <- as.numeric(design %*% coefficients)
    if (outcome == "logit_ra") {
      response <- inverse_logit(eta)
      derivative <- response * (1 - response)
    } else {
      response <- exp(eta)
      derivative <- response
    }
    row_gradients <- design * derivative
    gradient_frame <- as.data.table(row_gradients)
    gradient_frame[, caseid := segmented$data$caseid]
    gradient <- as.numeric(
      gradient_frame[, lapply(.SD, mean), by = caseid][
        , lapply(.SD, mean), .SDcols = setdiff(names(gradient_frame), "caseid")
      ]
    )
    case_response <- data.table(
      caseid = segmented$data$caseid,
      response = response
    )[, .(response = mean(response)), by = caseid]
    list(
      target = target,
      estimate = mean(case_response$response),
      gradient = gradient,
      design_names = fixed_names,
      covariance = covariance
    )
  }

  levels <- lapply(targets, one_level)
  estimates <- rbindlist(lapply(levels, function(item) {
    standard_error <- sqrt(as.numeric(
      t(item$gradient) %*% item$covariance %*% item$gradient
    ))
    data.table(
      map_mmhg = item$target,
      estimate_response = item$estimate,
      standard_error_response = standard_error,
      lower_95ci_response = item$estimate - 1.959963984540054 * standard_error,
      upper_95ci_response = item$estimate + 1.959963984540054 * standard_error
    )
  }))

  contrasts <- rbindlist(lapply(list(c(60, 65), c(65, 70), c(60, 70)), function(pair) {
    target <- levels[[match(pair[1], targets)]]
    reference <- levels[[match(pair[2], targets)]]
    gradient <- target$gradient - reference$gradient
    difference <- target$estimate - reference$estimate
    standard_error <- sqrt(as.numeric(
      t(gradient) %*% target$covariance %*% gradient
    ))
    data.table(
      target_map_mmhg = pair[1],
      reference_map_mmhg = pair[2],
      contrast_response = difference,
      standard_error_response = standard_error,
      lower_95ci_response = difference - 1.959963984540054 * standard_error,
      upper_95ci_response = difference + 1.959963984540054 * standard_error
    )
  }))
  list(levels = estimates, contrasts = contrasts)
}

analyse_bundle <- function(bundle, cohort_label, outcome_label, outcome) {
  segmented <- build_segmented_bundle(bundle, threshold)
  message(sprintf("[%s/%s] fitting segmented linear model", cohort_label, outcome_label))
  linear_fit <- fit_car1_random_map_slope_lme(
    segmented, segmented$segmented_linear_names, outcome = outcome
  )
  message(sprintf("[%s/%s] fitting 65-mm Hg hinge model", cohort_label, outcome_label))
  hinge_fit <- fit_car1_random_map_slope_lme(
    segmented, segmented$segmented_hinge_names, outcome = outcome
  )
  message(sprintf("[%s/%s] fitting restricted cubic spline comparator", cohort_label, outcome_label))
  rcs_fit <- if (
    cohort_label == "broad_all_windows" && outcome_label == "logit_ra"
  ) main_objects$nonlinear else fit_car1_random_map_slope_lme(
    bundle, bundle$nonlinear_names, outcome = outcome
  )

  statistic <- max(
    2 * (as.numeric(logLik(hinge_fit)) - as.numeric(logLik(linear_fit))), 0
  )
  coefficient_table <- summary(hinge_fit)$tTable
  hinge_name <- "hinge65_within_10_raw"
  map_name <- "map_within_10_raw"

  slope_table <- rbindlist(list(
    linear_combination(
      hinge_fit,
      setNames(c(-0.5, 0), c(map_name, hinge_name)),
      "5-mm Hg MAP decrease at MAP >=65"
    ),
    linear_combination(
      hinge_fit,
      setNames(c(-0.5, 0.5), c(map_name, hinge_name)),
      "5-mm Hg MAP decrease at MAP <65"
    ),
    linear_combination(
      hinge_fit,
      setNames(0.5, hinge_name),
      "difference in 5-mm Hg decrease slopes: below minus above 65"
    )
  ))
  slope_table[, `:=`(cohort = cohort_label, outcome = outcome_label)]
  setcolorder(slope_table, c("cohort", "outcome", "estimand"))

  map_terms <- c(map_name, hinge_name, "map_mean_10", "hinge65_mean_10")
  map_covariance <- vcov(hinge_fit)[map_terms, map_terms, drop = FALSE]
  map_correlation <- cov2cor(map_covariance)
  off_diagonal <- abs(map_correlation[row(map_correlation) != col(map_correlation)])

  response <- standardised_levels(
    hinge_fit, segmented, c(60, 65, 70), outcome_label
  )
  response$levels[, `:=`(cohort = cohort_label, outcome = outcome_label)]
  response$contrasts[, `:=`(cohort = cohort_label, outcome = outcome_label)]

  summary_row <- data.table(
    cohort = cohort_label,
    outcome = outcome_label,
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data),
    hinge_beta_per_10_mmhg_hinge = fixef(hinge_fit)[hinge_name],
    hinge_standard_error = coefficient_table[hinge_name, "Std.Error"],
    hinge_lower_95ci = fixef(hinge_fit)[hinge_name] -
      1.959963984540054 * coefficient_table[hinge_name, "Std.Error"],
    hinge_upper_95ci = fixef(hinge_fit)[hinge_name] +
      1.959963984540054 * coefficient_table[hinge_name, "Std.Error"],
    hinge_wald_p_value = coefficient_table[hinge_name, "p-value"],
    hinge_lrt_statistic = statistic,
    hinge_lrt_df = 1,
    hinge_lrt_p_value = pchisq(statistic, 1, lower.tail = FALSE),
    linear_AIC = AIC(linear_fit),
    hinge_AIC = AIC(hinge_fit),
    rcs_AIC = AIC(rcs_fit),
    linear_BIC = BIC(linear_fit),
    hinge_BIC = BIC(hinge_fit),
    rcs_BIC = BIC(rcs_fit),
    maximum_absolute_correlation_among_map_terms = max(off_diagonal),
    hinge_iterations = if (is.null(hinge_fit$numIter)) NA_integer_ else as.integer(hinge_fit$numIter)
  )
  list(
    summary = summary_row,
    slopes = slope_table,
    levels = response$levels,
    contrasts = response$contrasts,
    fit = hinge_fit,
    segmented = segmented
  )
}

case_ranges <- full_bundle$data[, .(
  minimum_map = min(art_map), maximum_map = max(art_map)
), by = caseid]
crossing_ids <- case_ranges[
  minimum_map < threshold & maximum_map >= threshold, caseid
]

strict_bundle <- build_main_design(
  windows,
  strict_measured_zero = TRUE,
  minimum_windows = 8L,
  fixed_knots = full_bundle$knots
)
overlap_bundle <- build_main_design(
  windows[caseid %in% crossing_ids],
  minimum_windows = 8L,
  fixed_knots = full_bundle$knots
)

analyses <- list(
  analyse_bundle(
    full_bundle, "broad_all_windows", "logit_ra", full_bundle$y
  ),
  analyse_bundle(
    strict_bundle, "strict_measured_zero_volatile", "logit_ra", strict_bundle$y
  ),
  analyse_bundle(
    overlap_bundle, "surgeries_crossing_65", "logit_ra", overlap_bundle$y
  ),
  analyse_bundle(
    full_bundle, "broad_all_windows", "log_alpha_power",
    log(as.numeric(full_bundle$data$alpha_power))
  )
)

summary_table <- rbindlist(lapply(analyses, `[[`, "summary"), fill = TRUE)
slope_table <- rbindlist(lapply(analyses, `[[`, "slopes"), fill = TRUE)
level_table <- rbindlist(lapply(analyses, `[[`, "levels"), fill = TRUE)
contrast_table <- rbindlist(lapply(analyses, `[[`, "contrasts"), fill = TRUE)

if (any(!is.finite(unlist(summary_table[, .(
  hinge_beta_per_10_mmhg_hinge, hinge_standard_error,
  hinge_lrt_statistic, hinge_AIC, rcs_AIC
)])))) stop("Non-finite hinge-model result")

fwrite(summary_table, file.path(output_dir, "hinge65_model_summary.csv"))
fwrite(slope_table, file.path(output_dir, "hinge65_log_scale_slopes.csv"))
fwrite(level_table, file.path(output_dir, "hinge65_standardised_levels.csv"))
fwrite(contrast_table, file.path(output_dir, "hinge65_standardised_contrasts.csv"))
saveRDS(
  lapply(analyses, function(item) list(fit = item$fit, segmented = item$segmented)),
  file.path(output_dir, "hinge65_model_objects.rds"),
  compress = "xz"
)
write_json_pretty(
  list(
    threshold_mmhg = threshold,
    interpretation = paste(
      "The within-surgery hinge coefficient tests a change in slope below the",
      "operational threshold; it is not a causal or biological-breakpoint test."
    ),
    crossing_definition = "observed surgery-level MAP minimum <65 and maximum >=65 mm Hg",
    crossing_surgeries = length(crossing_ids),
    model_summary = summary_table,
    log_scale_slopes = slope_table,
    standardised_levels = level_table,
    standardised_contrasts = contrast_table
  ),
  file.path(output_dir, "hinge65_results.json")
)

print(summary_table)
print(slope_table)
print(contrast_table[outcome == "logit_ra"])
