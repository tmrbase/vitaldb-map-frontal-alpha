# Patient-level ranking for Table SIII. Validate the upstream cohort before
# excluding all surgeries belonging to the top ceiling(1% * N) patients.
select_influential_patients <- function(influence, subjectids) {
  required <- c("subjectid", "maximum_absolute_curve_change_logit_ra")
  if (!all(required %in% names(influence))) {
    stop("Influence input lacks the required patient/ranking columns")
  }
  patients <- unique(as.character(subjectids))
  if (anyNA(patients) || !length(patients) || anyDuplicated(influence$subjectid) ||
      !setequal(as.character(influence$subjectid), patients)) {
    stop("Influence ranking must contain exactly one row per primary-cohort patient")
  }
  score <- influence$maximum_absolute_curve_change_logit_ra
  if (any(!is.finite(score)) || any(score < 0)) {
    stop("Influence scores must be finite and nonnegative")
  }
  ranked <- influence[order(-score)]$subjectid
  head(ranked, max(1L, ceiling(0.01 * length(patients))))
}
