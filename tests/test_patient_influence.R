# Run from repository root with Rscript tests/test_patient_influence.R.
suppressPackageStartupMessages(library(data.table))
source("r_analysis/maintenance_ra/R/influence.R")
ranking <- data.table(subjectid = 1:101, maximum_absolute_curve_change_logit_ra = 1:101)
# 101 patients require rounding up to two; repeated surgeries do not change N.
patients <- c(1:101, 101, 101)
selected <- select_influential_patients(ranking, patients)
stopifnot(identical(selected, c(101L, 100L)))
surgeries <- data.table(subjectid = patients, caseid = seq_along(patients))
retained <- surgeries[!subjectid %in% selected]
stopifnot(nrow(retained) == 99L, uniqueN(retained$subjectid) == 99L)
# A stale, duplicate, or nonfinite ranking must fail before any refitting.
fails <- function(x) inherits(try(select_influential_patients(x, patients), silent = TRUE), "try-error")
stopifnot(fails(ranking[-1]), fails(rbind(ranking, ranking[1])))
ranking[1, maximum_absolute_curve_change_logit_ra := NA_real_]
stopifnot(fails(ranking))
cat("Patient exclusion checks passed\n")
