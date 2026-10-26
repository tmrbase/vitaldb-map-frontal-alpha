# Purpose: calculate patient-cluster CR1S inference for the 65-mm Hg hinge model.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(clubSandwich)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript hinge65_cluster_robust.R MODEL_RDS OUTPUT_CSV")
}
objects <- readRDS(arguments[[1L]])

cluster_robust_covariance <- function(fit, segmented) {
  fixed_names <- segmented$segmented_hinge_names
  clusters <- uniqueN(segmented$data$subjectid)
  observations <- nrow(segmented$data)
  parameters <- length(fixed_names)
  cr1s_multiplier <- clusters / (clusters - 1) *
    (observations - 1) / (observations - parameters)
  robust <- as.matrix(clubSandwich::vcovCR(
    fit,
    cluster = segmented$data$subjectid,
    type = "CR1S"
  ))[fixed_names, fixed_names, drop = FALSE]
  model_covariance <- vcov(fit)[fixed_names, fixed_names, drop = FALSE]
  list(
    robust = robust,
    model_covariance = model_covariance,
    clusters = clusters,
    observations = observations,
    parameters = parameters,
    cr1s_multiplier = cr1s_multiplier,
    engine = "clubSandwich::vcovCR(type='CR1S')"
  )
}

labels <- c(
  "broad_all_windows_logit_ra",
  "strict_measured_zero_volatile_logit_ra",
  "surgeries_crossing_65_logit_ra",
  "broad_all_windows_log_alpha_power"
)
rows <- lapply(seq_along(objects), function(index) {
  fit <- objects[[index]]$fit
  segmented <- objects[[index]]$segmented
  covariance <- cluster_robust_covariance(fit, segmented)
  name <- "hinge65_within_10_raw"
  map_name <- "map_within_10_raw"
  beta <- fixef(fit)[name]
  model_se <- sqrt(vcov(fit)[name, name])
  robust_se <- sqrt(covariance$robust[name, name])
  coefficients <- fixef(fit)
  above_vector <- setNames(rep(0, length(coefficients)), names(coefficients))
  below_vector <- above_vector
  above_vector[map_name] <- -0.5
  below_vector[map_name] <- -0.5
  below_vector[name] <- 0.5
  above_estimate <- sum(above_vector * coefficients)
  below_estimate <- sum(below_vector * coefficients)
  above_robust_se <- sqrt(as.numeric(
    t(above_vector) %*% covariance$robust %*% above_vector
  ))
  below_robust_se <- sqrt(as.numeric(
    t(below_vector) %*% covariance$robust %*% below_vector
  ))
  data.table(
    analysis = labels[index],
    patient_clusters = covariance$clusters,
    observations = covariance$observations,
    fixed_effect_parameters = covariance$parameters,
    hinge_beta = beta,
    model_based_standard_error = model_se,
    patient_cluster_CR1S_standard_error = robust_se,
    patient_cluster_CR1S_lower_95ci = beta - 1.959963984540054 * robust_se,
    patient_cluster_CR1S_upper_95ci = beta + 1.959963984540054 * robust_se,
    patient_cluster_CR1S_p_value = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE),
    five_mmhg_decrease_above65_estimate = above_estimate,
    five_mmhg_decrease_above65_CR1S_se = above_robust_se,
    five_mmhg_decrease_above65_CR1S_lower =
      above_estimate - 1.959963984540054 * above_robust_se,
    five_mmhg_decrease_above65_CR1S_upper =
      above_estimate + 1.959963984540054 * above_robust_se,
    five_mmhg_decrease_below65_estimate = below_estimate,
    five_mmhg_decrease_below65_CR1S_se = below_robust_se,
    five_mmhg_decrease_below65_CR1S_lower =
      below_estimate - 1.959963984540054 * below_robust_se,
    five_mmhg_decrease_below65_CR1S_upper =
      below_estimate + 1.959963984540054 * below_robust_se,
    CR1S_multiplier = covariance$cr1s_multiplier,
    covariance_engine = covariance$engine
  )
})
result <- rbindlist(rows)
fwrite(result, arguments[[2L]])
print(result)
