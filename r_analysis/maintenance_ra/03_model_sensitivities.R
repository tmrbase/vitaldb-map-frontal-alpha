# Purpose: Table SII random-slope diagnostics and Table SIII sensitivities.
source("R/common.R")
source("R/influence.R")

suppressPackageStartupMessages({
  library(clubSandwich)
})

main_objects <- readRDS(file.path(output_root, "r_random_slope_primary_objects.rds"))
random_intercept_results <- fromJSON(
  file.path(output_root, "r_main_model_results.json"),
  simplifyVector = FALSE
)
windows <- fread(windows_path)
reference <- main_objects$bundle

fit_independent_lme <- function(bundle, fixed_names, outcome = NULL) {
  fit_data <- cbind(
    data.table(
      y = if (is.null(outcome)) bundle$y else as.numeric(outcome),
      caseid = factor(bundle$data$caseid),
      subjectid = factor(bundle$data$subjectid),
      time5 = as.numeric(bundle$data$time_seconds) / 300,
      random_map10 = (
        as.numeric(bundle$data$art_map) -
          ave(as.numeric(bundle$data$art_map), bundle$data$caseid, FUN = mean)
      ) / 10
    ),
    bundle$design[, fixed_names, with = FALSE]
  )
  fixed_formula <- as.formula(
    sprintf("y ~ -1 + %s", paste(fixed_names, collapse = " + "))
  )
  lme(
    fixed = fixed_formula,
    random = ~1 + random_map10 | caseid,
    data = fit_data,
    method = "ML",
    na.action = na.fail,
    control = lmeControl(
      opt = "optim", msMaxIter = 300, msMaxEval = 1000,
      tolerance = 1e-7, niterEM = 25, returnObject = TRUE
    )
  )
}

fit_model_set <- function(bundle, outcome = NULL, correlation = TRUE) {
  fitter <- if (correlation) fit_car1_random_map_slope_lme else fit_independent_lme
  list(
    baseline = fitter(bundle, bundle$baseline_names, outcome = outcome),
    linear = fitter(bundle, bundle$linear_names, outcome = outcome),
    nonlinear = fitter(bundle, bundle$nonlinear_names, outcome = outcome)
  )
}

summarize_set <- function(bundle, fits, label, outcome = NULL) {
  y <- if (is.null(outcome)) bundle$y else as.numeric(outcome)
  pred0 <- fixed_prediction(fits$baseline, bundle, bundle$baseline_names)
  pred1 <- fixed_prediction(fits$nonlinear, bundle, bundle$nonlinear_names)
  reduction <- case_equal_within_reduction(bundle$data, y, pred0, pred1)
  random_map10 <- (
    as.numeric(bundle$data$art_map) -
      ave(as.numeric(bundle$data$art_map), bundle$data$caseid, FUN = mean)
  ) / 10
  r2_0 <- mixed_r2_random_slope(fits$baseline, pred0, random_map10)
  r2_1 <- mixed_r2_random_slope(fits$nonlinear, pred1, random_map10)
  overall <- max(2 * (as.numeric(logLik(fits$nonlinear)) - as.numeric(logLik(fits$baseline))), 0)
  nonlinear <- max(2 * (as.numeric(logLik(fits$nonlinear)) - as.numeric(logLik(fits$linear))), 0)
  contrast <- map_contrast_matrix(bundle)
  list(
    label = label,
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data),
    case_equal_G = reduction$G,
    window_weighted_G = window_weighted_within_reduction(bundle$data, y, pred0, pred1),
    baseline_marginal_r2 = r2_0$marginal_r2,
    map_marginal_r2 = r2_1$marginal_r2,
    delta_marginal_r2 = r2_1$marginal_r2 - r2_0$marginal_r2,
    overall_map_lrt = overall,
    overall_map_df = length(bundle$nonlinear_names) - length(bundle$baseline_names),
    overall_map_p = pchisq(
      overall,
      length(bundle$nonlinear_names) - length(bundle$baseline_names),
      lower.tail = FALSE
    ),
    nonlinearity_lrt = nonlinear,
    nonlinearity_df = length(bundle$nonlinear_names) - length(bundle$linear_names),
    nonlinearity_p = pchisq(
      nonlinear,
      length(bundle$nonlinear_names) - length(bundle$linear_names),
      lower.tail = FALSE
    ),
    contrast_50_vs_65 = as.numeric(contrast[1, ] %*% fixef(fits$nonlinear)[bundle$nonlinear_names]),
    contrast_90_vs_65 = as.numeric(contrast[nrow(contrast), ] %*% fixef(fits$nonlinear)[bundle$nonlinear_names]),
    baseline = fit_components(fits$baseline),
    nonlinear = fit_components(fits$nonlinear)
  )
}

equal_spaced_windows <- function(data, count = 8L) {
  copy(data)[order(caseid, time_seconds), {
    positions <- unique(round(seq(1, .N, length.out = count)))
    if (length(positions) == count) .SD[positions] else .SD[0]
  }, by = caseid][, caseid := as.character(caseid)]
}

single_surgery_per_subject <- function(data) {
  first_cases <- unique(data[, .(subjectid, caseid)])[order(subjectid, caseid), .SD[1], by = subjectid]$caseid
  copy(data[caseid %in% first_cases])
}

message("[1/7] Independent-residual sensitivity")
independent_fits <- fit_model_set(reference, correlation = FALSE)
independent_summary <- summarize_set(
  reference, independent_fits, "independent_residual"
)

message("[2/7] Raw RA and arcsine-square-root outcome sensitivities")
ra <- pmin(pmax(reference$data$relative_alpha, 0), 1)
raw_fits <- fit_model_set(reference, outcome = ra)
arcsin_values <- asin(sqrt(ra))
arcsin_fits <- fit_model_set(reference, outcome = arcsin_values)
raw_summary <- summarize_set(reference, raw_fits, "raw_ra", ra)
arcsin_summary <- summarize_set(
  reference, arcsin_fits, "arcsin_sqrt_ra", arcsin_values
)

message("[3/7] Window-count sensitivities")
fixed_knots <- reference$knots
equal8_bundle <- build_main_design(
  equal_spaced_windows(reference$data, 8L), minimum_windows = 8L,
  fixed_knots = fixed_knots
)
minimum12_bundle <- build_main_design(
  windows, minimum_windows = 12L, fixed_knots = fixed_knots
)
minimum20_bundle <- build_main_design(
  windows, minimum_windows = 20L, fixed_knots = fixed_knots
)
strict_volatile_bundle <- build_main_design(
  windows,
  strict_measured_zero = TRUE,
  minimum_windows = 8L,
  fixed_knots = fixed_knots
)
window_summaries <- list()
window_bundles <- list(
  equal_spaced_8 = equal8_bundle,
  minimum_12 = minimum12_bundle,
  minimum_20 = minimum20_bundle,
  strict_measured_zero_volatile = strict_volatile_bundle
)
for (key in names(window_bundles)) {
  local_bundle <- window_bundles[[key]]
  window_summaries[[key]] <- summarize_set(
    local_bundle, fit_model_set(local_bundle), key
  )
}

message("[4/7] One-surgery-per-patient sensitivity")
single_bundle <- build_main_design(
  single_surgery_per_subject(reference$data), minimum_windows = 8L,
  fixed_knots = fixed_knots
)
single_summary <- summarize_set(
  single_bundle, fit_model_set(single_bundle), "single_surgery_per_subject"
)

message("[5/7] Top-one-percent influence-trimmed sensitivity")
influence_path <- file.path(
  repo_root, "outputs", "maintenance_ra_model_defense_v1",
  "subject_deletion_influence.csv"
)
influence <- fread(influence_path)
excluded_subjects <- select_influential_patients(influence, reference$data$subjectid)
trimmed_bundle <- build_main_design(
  windows[!subjectid %in% excluded_subjects], minimum_windows = 8L,
  fixed_knots = fixed_knots
)
trimmed_summary <- summarize_set(
  trimmed_bundle, fit_model_set(trimmed_bundle), "top_one_percent_influence_removed"
)
trimmed_summary$excluded_subjectids <- as.list(excluded_subjects)

message("[6/7] Random MAP-slope sensitivity")
random_data <- cbind(
  data.table(
    y = reference$y,
    caseid = factor(reference$data$caseid),
    time5 = reference$data$time_seconds / 300,
    random_map10 = (
      reference$data$art_map - ave(reference$data$art_map, reference$data$caseid, FUN = mean)
    ) / 10
  ),
  reference$design[, reference$nonlinear_names, with = FALSE]
)
random_formula <- as.formula(sprintf(
  "y ~ -1 + %s", paste(reference$nonlinear_names, collapse = " + ")
))
random_slope_fit <- lme(
  fixed = random_formula,
  random = ~1 + random_map10 | caseid,
  correlation = corCAR1(value = 0.5, form = ~time5 | caseid),
  data = random_data,
  method = "ML",
  na.action = na.fail,
  control = lmeControl(
    opt = "optim", msMaxIter = 500, msMaxEval = 2000,
    tolerance = 1e-7, niterEM = 50, returnObject = TRUE
  )
)
random_covariance <- as.matrix(getVarCov(random_slope_fit, type = "random.effects"))
random_slope_summary <- list(
  log_likelihood = as.numeric(logLik(random_slope_fit)),
  AIC = AIC(random_slope_fit),
  BIC = BIC(random_slope_fit),
  random_intercept_sd = sqrt(random_covariance[1, 1]),
  random_map_slope_sd_per_10_mmhg = sqrt(random_covariance[2, 2]),
  intercept_slope_correlation = random_covariance[1, 2] /
    sqrt(random_covariance[1, 1] * random_covariance[2, 2]),
  residual_variance = sigma(random_slope_fit)^2,
  correlation_at_five_minutes = as.numeric(
    coef(random_slope_fit$modelStruct$corStruct, unconstrained = FALSE)
  ),
  comparison_random_intercept_AIC = as.numeric(random_intercept_results$nonlinear$AIC),
  regular_lrt_not_used = TRUE
)

# Quantify whether surgery-specific MAP slopes are estimated precisely and
# reproduce across the first and second halves of the same surgery.  The
# conditional covariance below uses the fitted random-effect covariance and
# CAR(1) residual covariance; reliability is 1 minus the posterior-to-prior
# variance ratio for the random slope.
random_slope_blup <- function(indices) {
  local <- random_data[indices]
  z <- cbind(1, local$random_map10)
  distance <- abs(outer(local$time5, local$time5, "-"))
  rho <- random_slope_summary$correlation_at_five_minutes
  residual_covariance <- random_slope_summary$residual_variance * rho^distance
  marginal_covariance <- z %*% random_covariance %*% t(z) + residual_covariance
  diag(marginal_covariance) <- diag(marginal_covariance) + 1e-8
  fixed_design <- as.matrix(local[, reference$nonlinear_names, with = FALSE])
  residual <- local$y - as.numeric(
    fixed_design %*% fixef(random_slope_fit)[reference$nonlinear_names]
  )
  inverse_residual <- solve(marginal_covariance, residual)
  inverse_random <- solve(marginal_covariance, z)
  estimate <- random_covariance %*% t(z) %*% inverse_residual
  posterior_covariance <- random_covariance -
    random_covariance %*% t(z) %*% inverse_random %*% random_covariance
  reliability <- max(
    0,
    min(1, 1 - posterior_covariance[2, 2] / random_covariance[2, 2])
  )
  list(
    slope = as.numeric(estimate[2]),
    reliability = as.numeric(reliability),
    map_range = diff(range(local$random_map10 * 10))
  )
}

case_indices <- split(seq_len(nrow(random_data)), random_data$caseid)
individual_diagnostics <- rbindlist(lapply(names(case_indices), function(case_id) {
  indices <- case_indices[[case_id]]
  estimate <- random_slope_blup(indices)
  data.table(
    caseid = as.character(case_id),
    windows = length(indices),
    random_slope = estimate$slope,
    reliability = estimate$reliability,
    map_range = estimate$map_range
  )
}))

split_diagnostics <- rbindlist(lapply(names(case_indices), function(case_id) {
  indices <- case_indices[[case_id]]
  ordered <- indices[order(random_data$time5[indices])]
  split_point <- floor(length(ordered) / 2)
  early <- ordered[seq_len(split_point)]
  late <- ordered[split_point + seq_len(length(ordered) - split_point)]
  if (min(length(early), length(late)) < 4L) return(NULL)
  early_estimate <- random_slope_blup(early)
  late_estimate <- random_slope_blup(late)
  data.table(
    caseid = as.character(case_id),
    early_slope = early_estimate$slope,
    late_slope = late_estimate$slope,
    early_reliability = early_estimate$reliability,
    late_reliability = late_estimate$reliability,
    early_map_range = early_estimate$map_range,
    late_map_range = late_estimate$map_range
  )
}))

temporal_supported <- split_diagnostics[
  early_map_range >= 10 & late_map_range >= 10
]
random_slope_summary$individual_estimation <- list(
  median_reliability = median(individual_diagnostics$reliability),
  reliability_q25 = as.numeric(quantile(individual_diagnostics$reliability, 0.25)),
  reliability_q75 = as.numeric(quantile(individual_diagnostics$reliability, 0.75)),
  surgeries_reliability_ge_0_50 = sum(individual_diagnostics$reliability >= 0.50),
  surgeries_reliability_ge_0_70 = sum(individual_diagnostics$reliability >= 0.70)
)
random_slope_summary$temporal_reproducibility <- list(
  all_n = nrow(split_diagnostics),
  all_spearman = cor(
    split_diagnostics$early_slope,
    split_diagnostics$late_slope,
    method = "spearman"
  ),
  map_range_ge_10_n = nrow(temporal_supported),
  map_range_ge_10_spearman = cor(
    temporal_supported$early_slope,
    temporal_supported$late_slope,
    method = "spearman"
  )
)

message("[7/7] Patient-cluster robust block tests and serialization")
map_terms <- grep("^map_basis.*_within$", reference$nonlinear_names, value = TRUE)
nonlinear_terms <- setdiff(map_terms, "map_basis1_within")
map_positions <- match(map_terms, names(fixef(main_objects$nonlinear)))
nonlinear_positions <- match(nonlinear_terms, names(fixef(main_objects$nonlinear)))
robust_overall <- tryCatch(
  Wald_test(
    main_objects$nonlinear,
    constraints = constrain_zero(map_positions),
    vcov = "CR2", cluster = reference$data$subjectid, test = "HTZ"
  ),
  error = function(e) list(error = conditionMessage(e))
)
robust_nonlinear <- tryCatch(
  Wald_test(
    main_objects$nonlinear,
    constraints = constrain_zero(nonlinear_positions),
    vcov = "CR2", cluster = reference$data$subjectid, test = "HTZ"
  ),
  error = function(e) list(error = conditionMessage(e))
)

primary_summary <- summarize_set(
  reference,
  list(
    baseline = main_objects$baseline,
    linear = main_objects$linear,
    nonlinear = main_objects$nonlinear
  ),
  "primary_car1_logit_ra"
)
all_summaries <- c(
  list(
    primary = primary_summary,
    independent_residual = independent_summary,
    raw_ra = raw_summary,
    arcsin_sqrt_ra = arcsin_summary
  ),
  window_summaries,
  list(
    single_surgery_per_subject = single_summary,
    top_one_percent_influence_removed = trimmed_summary
  )
)

comparison <- rbindlist(lapply(names(all_summaries), function(key) {
  x <- all_summaries[[key]]
  data.table(
    analysis = key,
    subjects = x$subjects,
    surgeries = x$surgeries,
    windows = x$windows,
    case_equal_G = x$case_equal_G,
    window_weighted_G = x$window_weighted_G,
    delta_marginal_r2 = x$delta_marginal_r2,
    overall_map_p = x$overall_map_p,
    nonlinearity_p = x$nonlinearity_p,
    contrast_50_vs_65 = x$contrast_50_vs_65,
    contrast_90_vs_65 = x$contrast_90_vs_65
  )
}))

robust_to_list <- function(x) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x$error)) return(x)
  as.list(as.data.frame(x)[1, ])
}

results <- list(
  status = "R reanalysis of post-hoc model-defense sensitivities",
  summaries = all_summaries,
  random_map_slope = random_slope_summary,
  patient_cluster_robust = list(
    overall_map = robust_to_list(robust_overall),
    nonlinearity = robust_to_list(robust_nonlinear)
  )
)
fwrite(comparison, file.path(output_root, "r_model_sensitivity_comparison.csv"))
fwrite(
  individual_diagnostics,
  file.path(output_root, "r_random_slope_individual_diagnostics.csv")
)
fwrite(
  split_diagnostics,
  file.path(output_root, "r_random_slope_temporal_diagnostics.csv")
)
write_json_pretty(results, file.path(output_root, "r_model_sensitivity_results.json"))
saveRDS(
  list(
    independent = independent_fits,
    raw = raw_fits,
    arcsin = arcsin_fits,
    random_slope = random_slope_fit
  ),
  file.path(output_root, "r_model_sensitivity_objects.rds"),
  compress = "xz"
)
print(comparison)
