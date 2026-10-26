# Purpose: obtain patient-bootstrap uncertainty for the primary random-slope estimands.
source("R/common.R")

args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args) >= 1L) as.integer(args[[1L]]) else 2000L
workers <- if (length(args) >= 2L) as.integer(args[[2L]]) else 24L
batch_size <- if (length(args) >= 3L) as.integer(args[[3L]]) else 100L
seed <- 20260827L

object_path <- file.path(output_root, "r_random_slope_primary_objects.rds")
if (!file.exists(object_path)) stop("Run 18_random_slope_primary.R first")
main <- readRDS(object_path)
bundle <- main$bundle

resample_subject_bundle <- function(bundle, replicate_id, seed) {
  set.seed(seed + as.integer(replicate_id))
  subjects <- unique(bundle$data$subjectid)
  sampled <- sample(subjects, size = length(subjects), replace = TRUE)
  indices <- integer()
  resampled_groups <- integer()
  next_group <- 0L
  subject_vector <- bundle$data$subjectid
  case_vector <- bundle$data$caseid
  for (subject in sampled) {
    subject_indices <- which(subject_vector == subject)
    subject_cases <- case_vector[subject_indices]
    for (case in unique(subject_cases)) {
      local <- subject_indices[subject_cases == case]
      indices <- c(indices, local)
      next_group <- next_group + 1L
      resampled_groups <- c(resampled_groups, rep(next_group, length(local)))
    }
  }
  out <- bundle
  out$data <- copy(bundle$data[indices])
  out$data[, caseid := resampled_groups]
  out$design <- copy(bundle$design[indices])
  out$y <- bundle$y[indices]
  out
}

bootstrap_one_random_slope <- function(replicate_id, bundle, seed) {
  tryCatch({
    sampled <- resample_subject_bundle(bundle, replicate_id, seed)
    fit0 <- fit_car1_random_map_slope_lme(
      sampled,
      sampled$baseline_names,
      rho_starts = 0.5
    )
    fit1 <- fit_car1_random_map_slope_lme(
      sampled,
      sampled$nonlinear_names,
      rho_starts = 0.5
    )

    train0 <- fixed_prediction(fit0, sampled, sampled$baseline_names)
    train1 <- fixed_prediction(fit1, sampled, sampled$nonlinear_names)
    training <- case_equal_within_reduction(
      sampled$data,
      sampled$y,
      train0,
      train1
    )$G

    original0 <- as.numeric(
      as.matrix(bundle$design[, bundle$baseline_names, with = FALSE]) %*%
        fixef(fit0)[bundle$baseline_names]
    )
    original1 <- as.numeric(
      as.matrix(bundle$design[, bundle$nonlinear_names, with = FALSE]) %*%
        fixef(fit1)[bundle$nonlinear_names]
    )
    test <- case_equal_within_reduction(
      bundle$data,
      bundle$y,
      original0,
      original1
    )$G

    list(
      replicate = replicate_id,
      success = TRUE,
      training_G = training,
      original_test_G = test,
      optimism = training - test
    )
  }, error = function(e) {
    list(
      replicate = replicate_id,
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

bootstrap_root <- file.path(output_root, "random_slope_patient_bootstrap_batches")
dir.create(bootstrap_root, recursive = TRUE, showWarnings = FALSE)

project_dir <- getwd()
cluster <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(cluster, "project_dir", envir = environment())
parallel::clusterEvalQ(cluster, {
  setwd(project_dir)
  source("R/common.R")
  NULL
})
parallel::clusterExport(
  cluster,
  c("bundle", "seed", "resample_subject_bundle", "bootstrap_one_random_slope"),
  envir = environment()
)

for (start in seq.int(1L, repetitions, by = batch_size)) {
  stop_at <- min(start + batch_size - 1L, repetitions)
  batch_path <- file.path(
    bootstrap_root,
    sprintf("bootstrap_%04d_%04d.rds", start, stop_at)
  )
  if (file.exists(batch_path)) {
    message(sprintf("[resume] %d-%d already exists", start, stop_at))
    next
  }
  ids <- seq.int(start, stop_at)
  message(sprintf("[random-slope bootstrap] %d-%d", start, stop_at))
  batch <- parallel::parLapply(
    cluster,
    ids,
    bootstrap_one_random_slope,
    bundle = bundle,
    seed = seed
  )
  saveRDS(batch, batch_path, compress = "xz")
  success <- sum(vapply(batch, function(x) isTRUE(x$success), logical(1)))
  message(sprintf("[saved] %d/%d successful", success, length(batch)))
}

batch_files <- sort(list.files(bootstrap_root, pattern = "\\.rds$", full.names = TRUE))
records <- unlist(lapply(batch_files, readRDS), recursive = FALSE)
records <- records[vapply(records, function(x) x$replicate <= repetitions, logical(1))]
records <- records[order(vapply(records, `[[`, integer(1), "replicate"))]
successful <- records[vapply(records, function(x) isTRUE(x$success), logical(1))]

training <- vapply(successful, `[[`, numeric(1), "training_G")
test <- vapply(successful, `[[`, numeric(1), "original_test_G")
optimism <- training - test

pred0 <- fixed_prediction(main$baseline, bundle, bundle$baseline_names)
pred1 <- fixed_prediction(main$nonlinear, bundle, bundle$nonlinear_names)
apparent <- case_equal_within_reduction(bundle$data, bundle$y, pred0, pred1)$G

summary <- list(
  requested_repetitions = repetitions,
  successful_repetitions = length(successful),
  failed_repetitions = repetitions - length(successful),
  seed = seed,
  workers = workers,
  apparent_G = apparent,
  bootstrap_percentile_interval = as.list(
    as.numeric(quantile(training, c(0.025, 0.975)))
  ),
  mean_optimism = mean(optimism),
  optimism_corrected_G = apparent - mean(optimism),
  training_mean = mean(training),
  original_test_mean = mean(test)
)

write_json_pretty(
  summary,
  file.path(output_root, "r_random_slope_patient_bootstrap_summary.json")
)
fwrite(
  rbindlist(lapply(records, as.data.table), fill = TRUE),
  file.path(output_root, "r_random_slope_patient_bootstrap_records.csv")
)
print(summary)
