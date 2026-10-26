# Purpose: re-standardise MAP contrasts within surgeries spanning each target and 65 mm Hg.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript overlap_restricted_contrasts.R REPOSITORY_ROOT OUTPUT_DIR")
}
repository_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
source(file.path(repository_root, "r_analysis", "maintenance_ra", "R", "common.R"))
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

objects <- readRDS(
  file.path(repository_root, "outputs", "maintenance_ra_r_reanalysis_v1", "r_random_slope_primary_objects.rds")
)
bundle <- objects$bundle
fit <- objects$nonlinear
beta <- fixef(fit)[bundle$nonlinear_names]
cov_beta <- vcov(fit)[bundle$nonlinear_names, bundle$nonlinear_names, drop = FALSE]

inverse_logit <- function(x) 1 / (1 + exp(-x))

map_basis_observed <- natural_spline_basis(bundle$data$art_map, bundle$knots$map)
case_basis_means <- data.table(caseid = bundle$data$caseid)
for (j in seq_len(ncol(map_basis_observed))) {
  case_basis_means[[sprintf("basis%d", j)]] <- map_basis_observed[, j]
}
case_basis_means <- case_basis_means[, lapply(.SD, mean), by = caseid]

case_ranges <- bundle$data[, .(
  min_map = min(art_map),
  max_map = max(art_map),
  windows = .N
), by = caseid]

values_by_case <- split(bundle$data$art_map, bundle$data$caseid)

design_at <- function(target) {
  design <- copy(bundle$design[, bundle$nonlinear_names, with = FALSE])
  target_basis <- natural_spline_basis(target, bundle$knots$map)[1, ]
  case_index <- match(bundle$data$caseid, case_basis_means$caseid)
  for (j in seq_along(target_basis)) {
    name <- sprintf("map_basis%d_within", j)
    raw_within <- target_basis[j] - case_basis_means[[sprintf("basis%d", j)]][case_index]
    design[[name]] <- (
      raw_within - bundle$scaling[[name]]$mean
    ) / bundle$scaling[[name]]$standard_deviation
  }
  as.matrix(design)
}

standardised_contrast <- function(target, case_ids, support_definition) {
  target_design <- design_at(target)
  reference_design <- design_at(65)
  keep <- bundle$data$caseid %in% case_ids
  target_design <- target_design[keep, , drop = FALSE]
  reference_design <- reference_design[keep, , drop = FALSE]
  local_case <- bundle$data$caseid[keep]

  eta_target <- as.numeric(target_design %*% beta)
  eta_reference <- as.numeric(reference_design %*% beta)
  p_target <- inverse_logit(eta_target)
  p_reference <- inverse_logit(eta_reference)

  rows <- data.table(
    caseid = local_case,
    p_target = p_target,
    p_reference = p_reference
  )[, .(
    p_target = mean(p_target),
    p_reference = mean(p_reference)
  ), by = caseid]

  target_weight <- p_target * (1 - p_target)
  reference_weight <- p_reference * (1 - p_reference)
  grad_rows <- (
    target_design * target_weight -
      reference_design * reference_weight
  )
  grad_dt <- as.data.table(grad_rows)
  grad_dt[, caseid := local_case]
  grad_case <- grad_dt[, lapply(.SD, mean), by = caseid]
  gradient <- colMeans(as.matrix(grad_case[, -"caseid"]))
  standard_error <- sqrt(as.numeric(t(gradient) %*% cov_beta %*% gradient))

  estimate <- mean(rows$p_target - rows$p_reference)
  data.table(
    support_definition = support_definition,
    target_map_mmhg = target,
    reference_map_mmhg = 65,
    surgeries = length(unique(case_ids)),
    windows_used_for_standardisation = sum(keep),
    target_standardised_ra_percent = 100 * mean(rows$p_target),
    reference_standardised_ra_percent = 100 * mean(rows$p_reference),
    contrast_percentage_points = 100 * estimate,
    lower_95ci_percentage_points = 100 * (estimate - 1.96 * standard_error),
    upper_95ci_percentage_points = 100 * (estimate + 1.96 * standard_error),
    ci_method = "delta method using fixed-effect covariance; pointwise"
  )
}

results <- list()
for (target in c(55, 60, 65, 70, 80, 90, 100)) {
  lower <- min(target, 65)
  upper <- max(target, 65)
  span_ids <- case_ranges[
    min_map <= lower & max_map >= upper,
    caseid
  ]
  near_ids <- as.numeric(names(values_by_case)[vapply(
    values_by_case,
    function(values) {
      any(abs(values - target) <= 2.5) && any(abs(values - 65) <= 2.5)
    },
    logical(1)
  )])
  all_ids <- unique(bundle$data$caseid)

  results[[length(results) + 1L]] <- standardised_contrast(
    target,
    all_ids,
    "all surgeries (model-based standardisation)"
  )
  if (length(span_ids)) {
    results[[length(results) + 1L]] <- standardised_contrast(
      target,
      span_ids,
      "observed range spans target and 65 mm Hg"
    )
  }
  if (length(near_ids)) {
    results[[length(results) + 1L]] <- standardised_contrast(
      target,
      near_ids,
      "observed windows within +/-2.5 mm Hg of target and reference"
    )
  }
}

result <- rbindlist(results, use.names = TRUE, fill = TRUE)
fwrite(result, file.path(output_dir, "overlap_restricted_map_contrasts.csv"))
print(result[target_map_mmhg %in% c(60, 70, 80)])
