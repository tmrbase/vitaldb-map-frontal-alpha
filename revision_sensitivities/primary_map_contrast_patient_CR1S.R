# Purpose: calculate patient-cluster CR1S confidence intervals for the principal MAP contrasts.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(clubSandwich)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript primary_map_contrast_patient_CR1S.R REPOSITORY_ROOT OUTPUT_CSV")
}

repository_root <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
output_csv <- normalizePath(arguments[[2L]], winslash = "/", mustWork = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

source(file.path(
  repository_root,
  "r_analysis", "maintenance_ra", "R", "common.R"
))

objects <- readRDS(file.path(
  repository_root,
  "outputs", "maintenance_ra_r_reanalysis_v1",
  "r_random_slope_primary_objects.rds"
))
bundle <- objects$bundle
fit <- objects$nonlinear

patient_cluster_cr1s_covariance <- function(fit, bundle) {
  fixed_names <- bundle$nonlinear_names
  clusters <- uniqueN(bundle$data$subjectid)
  observations <- nrow(bundle$data)
  parameters <- length(fixed_names)
  cr1s_multiplier <- clusters / (clusters - 1) *
    (observations - 1) / (observations - parameters)
  robust <- as.matrix(clubSandwich::vcovCR(
    fit,
    cluster = bundle$data$subjectid,
    type = "CR1S"
  ))[fixed_names, fixed_names, drop = FALSE]
  model_covariance <- vcov(fit)[fixed_names, fixed_names, drop = FALSE]

  list(
    robust = robust,
    model = model_covariance,
    clusters = clusters,
    observations = observations,
    parameters = parameters,
    cr1s_multiplier = cr1s_multiplier,
    engine = "clubSandwich::vcovCR(type='CR1S')"
  )
}

inverse_logit <- function(x) 1 / (1 + exp(-x))

map_basis_observed <- natural_spline_basis(bundle$data$art_map, bundle$knots$map)
case_basis_means <- data.table(caseid = bundle$data$caseid)
for (j in seq_len(ncol(map_basis_observed))) {
  case_basis_means[[sprintf("basis%d", j)]] <- map_basis_observed[, j]
}
case_basis_means <- case_basis_means[, lapply(.SD, mean), by = caseid]

design_at <- function(target) {
  design <- copy(bundle$design[, bundle$nonlinear_names, with = FALSE])
  target_basis <- natural_spline_basis(target, bundle$knots$map)[1, ]
  case_index <- match(bundle$data$caseid, case_basis_means$caseid)
  for (j in seq_along(target_basis)) {
    name <- sprintf("map_basis%d_within", j)
    raw_within <- target_basis[j] -
      case_basis_means[[sprintf("basis%d", j)]][case_index]
    design[[name]] <- (
      raw_within - bundle$scaling[[name]]$mean
    ) / bundle$scaling[[name]]$standard_deviation
  }
  as.matrix(design)
}

covariance <- patient_cluster_cr1s_covariance(fit, bundle)
beta <- fixef(fit)[bundle$nonlinear_names]

standardised_contrast <- function(target, reference = 65) {
  target_design <- design_at(target)
  reference_design <- design_at(reference)
  eta_target <- as.numeric(target_design %*% beta)
  eta_reference <- as.numeric(reference_design %*% beta)
  p_target <- inverse_logit(eta_target)
  p_reference <- inverse_logit(eta_reference)
  case <- bundle$data$caseid

  case_values <- data.table(
    caseid = case,
    target = p_target,
    reference = p_reference
  )[, .(
    target = mean(target),
    reference = mean(reference)
  ), by = caseid]

  gradient_rows <-
    target_design * (p_target * (1 - p_target)) -
    reference_design * (p_reference * (1 - p_reference))
  gradient_table <- as.data.table(gradient_rows)
  gradient_table[, caseid := case]
  gradient_case <- gradient_table[, lapply(.SD, mean), by = caseid]
  gradient <- colMeans(as.matrix(gradient_case[, -"caseid"]))

  estimate <- mean(case_values$target - case_values$reference)
  model_se <- sqrt(as.numeric(
    t(gradient) %*% covariance$model %*% gradient
  ))
  robust_se <- sqrt(as.numeric(
    t(gradient) %*% covariance$robust %*% gradient
  ))

  data.table(
    target_map_mmhg = target,
    reference_map_mmhg = reference,
    surgeries = uniqueN(bundle$data$caseid),
    patient_clusters = covariance$clusters,
    windows = nrow(bundle$data),
    contrast_percentage_points = 100 * estimate,
    model_based_standard_error_percentage_points = 100 * model_se,
    model_based_lower_95ci_percentage_points =
      100 * (estimate - 1.959963984540054 * model_se),
    model_based_upper_95ci_percentage_points =
      100 * (estimate + 1.959963984540054 * model_se),
    patient_cluster_CR1S_standard_error_percentage_points = 100 * robust_se,
    patient_cluster_CR1S_lower_95ci_percentage_points =
      100 * (estimate - 1.959963984540054 * robust_se),
    patient_cluster_CR1S_upper_95ci_percentage_points =
      100 * (estimate + 1.959963984540054 * robust_se),
    CR1S_multiplier = covariance$cr1s_multiplier,
    covariance_engine = covariance$engine
  )
}

result <- rbindlist(lapply(c(60, 70, 80), standardised_contrast))
fwrite(result, output_csv)
print(result)
