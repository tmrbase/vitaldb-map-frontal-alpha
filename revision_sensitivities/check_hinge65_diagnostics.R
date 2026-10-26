# Purpose: write convergence and fit diagnostics for the 65-mm Hg hinge models.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript check_hinge65_diagnostics.R MODEL_RDS OUTPUT_CSV")
}
objects <- readRDS(arguments[[1L]])

rows <- lapply(seq_along(objects), function(index) {
  fit <- objects[[index]]$fit
  bundle <- objects[[index]]$segmented
  residual <- as.numeric(residuals(fit, type = "normalized"))
  centered <- residual - mean(residual)
  variance <- mean(centered^2)
  excess_kurtosis <- mean(centered^4) / variance^2 - 3
  acf_table <- ACF(fit, resType = "normalized", maxLag = 8)
  acf_nonzero <- acf_table$ACF[acf_table$lag > 0]
  ranges <- bundle$data[, .(
    minimum_map = min(art_map), maximum_map = max(art_map)
  ), by = caseid]
  data.table(
    analysis_index = index,
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data),
    windows_below65 = sum(bundle$data$art_map < 65),
    percent_windows_below65 = 100 * mean(bundle$data$art_map < 65),
    surgeries_with_below65 = uniqueN(bundle$data[art_map < 65, caseid]),
    surgeries_crossing65 = nrow(ranges[minimum_map < 65 & maximum_map >= 65]),
    normalised_residual_mean = mean(residual),
    normalised_residual_sd = sd(residual),
    normalised_residual_excess_kurtosis = excess_kurtosis,
    maximum_absolute_residual_acf_lags1_8 = max(abs(acf_nonzero)),
    minimum_map = min(bundle$data$art_map),
    maximum_map = max(bundle$data$art_map)
  )
})
result <- rbindlist(rows)
fwrite(result, arguments[[2L]])
print(result)
