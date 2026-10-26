# Purpose: create Table 1 inputs and cohort/missingness summaries for Table S1.
source("R/common.R")

windows <- fread(windows_path)
cases <- fread(file.path(
  repo_root, "outputs", "maintenance_ra_full_v1", "eligible_cases.csv"
))

continuous_summary <- function(x) {
  values <- as.numeric(x[is.finite(as.numeric(x))])
  if (!length(values)) return("NA")
  q <- quantile(values, c(0.25, 0.5, 0.75))
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

binary_summary <- function(x) {
  values <- as.numeric(x[is.finite(as.numeric(x))])
  if (!length(values)) return("NA")
  count <- sum(values > 0)
  sprintf("%d (%.1f%%)", count, 100 * count / length(values))
}

cases[, sex_male := as.numeric(toupper(as.character(sex)) == "M")]
cases[, asa_ge3 := as.numeric(asa >= 3)]
cases[, emergency_derived := as.numeric(emop > 0)]
cases[, thoracic := as.numeric(grepl("thoracic", department, ignore.case = TRUE))]
cases[, operation_minutes := (opend - opstart) / 60]

complete <- windows[logical_true(primary_complete_window)]
broad_ids <- complete[, .N, by = caseid][N >= 8]$caseid
strict_ids <- complete[logical_true(strict_measured_zero_volatile), .N, by = caseid][N >= 8]$caseid
cohorts <- list(
  all_eligible = cases,
  broad_analytic = cases[caseid %in% broad_ids],
  strict_measured_zero = cases[caseid %in% strict_ids],
  excluded_from_broad = cases[!caseid %in% broad_ids]
)
selection_flow <- list(
  eligible_surgeries = nrow(cases),
  eligible_patients = uniqueN(cases$subjectid),
  scheduled_windows = nrow(windows),
  complete_windows = sum(logical_true(windows$primary_complete_window)),
  excluded_surgeries = nrow(cases) - length(broad_ids),
  analytic_surgeries = length(broad_ids),
  analytic_patients = uniqueN(cases[caseid %in% broad_ids]$subjectid),
  analytic_windows = sum(logical_true(windows$primary_complete_window) & windows$caseid %in% broad_ids)
)

draw_cohort_flow <- function() {
  par(mar = rep(0.15, 4), xpd = NA, family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  box_text <- function(left, bottom, right, top, label, fill = "white", cex = 0.88) {
    rect(left, bottom, right, top, col = fill, border = "black", lwd = 1.2)
    text((left + right) / 2, (bottom + top) / 2, label, cex = cex)
  }
  text(0.5, 0.96, "Selection of surgeries and analysis windows", cex = 1.35, font = 2)
  box_text(
    0.12, 0.75, 0.88, 0.84,
    sprintf(
      "Eligible surgeries with no positive end-tidal volatile anaesthetic concentration\n%s surgeries; %s patients",
      format(selection_flow$eligible_surgeries, big.mark = ","),
      format(selection_flow$eligible_patients, big.mark = ",")
    )
  )
  arrows(0.5, 0.73, 0.5, 0.66, length = 0.08, lwd = 1.2)
  box_text(
    0.29, 0.56, 0.71, 0.64,
    sprintf(
      "Scheduled maintenance-phase windows\n%s windows",
      format(selection_flow$scheduled_windows, big.mark = ",")
    )
  )
  arrows(0.5, 0.54, 0.5, 0.47, length = 0.08, lwd = 1.2)
  box_text(
    0.22, 0.37, 0.78, 0.45,
    sprintf(
      "Windows meeting all primary signal and covariate criteria\n%s windows",
      format(selection_flow$complete_windows, big.mark = ",")
    )
  )
  arrows(0.5, 0.35, 0.5, 0.25, length = 0.08, lwd = 1.2)
  segments(0.66, 0.37, 0.77, 0.29, lwd = 1.2)
  box_text(
    0.68, 0.21, 0.97, 0.29,
    sprintf(
      "Excluded: fewer than 8 complete windows\n%s surgeries",
      format(selection_flow$excluded_surgeries, big.mark = ",")
    ),
    fill = "#F5F5F5",
    cex = 0.82
  )
  box_text(
    0.29, 0.10, 0.71, 0.19,
    sprintf(
      "Primary analysis\n%s surgeries; %s patients; %s windows",
      format(selection_flow$analytic_surgeries, big.mark = ","),
      format(selection_flow$analytic_patients, big.mark = ","),
      format(selection_flow$analytic_windows, big.mark = ",")
    ),
    fill = "#DCEEF7"
  )
}

png(file.path(release_asset_root, "FigureS1_cohort_flow.png"), width = 2100, height = 1476, res = 200)
draw_cohort_flow()
dev.off()
tiff(file.path(release_asset_root, "FigureS1_cohort_flow.tiff"), width = 10.5, height = 7.38, units = "in", res = 300, compression = "lzw")
draw_cohort_flow()
dev.off()
pdf(file.path(release_asset_root, "FigureS1_cohort_flow.pdf"), width = 10.5, height = 7.38, useDingbats = FALSE)
draw_cohort_flow()
dev.off()

continuous <- c(
  age = "Age, years", bmi = "Body mass index, kg/m2",
  preop_hb = "Preoperative haemoglobin, g/dL",
  preop_alb = "Preoperative albumin, g/dL",
  preop_cr = "Preoperative creatinine, mg/dL",
  preop_gluc = "Preoperative glucose, mg/dL",
  intraop_ebl = "Estimated blood loss, mL",
  intraop_crystalloid = "Crystalloid, mL",
  operation_minutes = "Operation duration, min"
)
binary <- c(
  sex_male = "Male sex", asa_ge3 = "ASA physical status >=3",
  emergency_derived = "Emergency surgery", preop_htn = "Preoperative hypertension",
  preop_dm = "Preoperative diabetes mellitus", thoracic = "Thoracic surgery"
)

rows <- list(
  c(Variable = "Surgeries, n", vapply(cohorts, nrow, integer(1))),
  c(Variable = "Patients, n", vapply(cohorts, function(x) uniqueN(x$subjectid), integer(1)))
)
for (name in names(continuous)) {
  rows[[length(rows) + 1L]] <- c(
    Variable = unname(continuous[name]),
    vapply(cohorts, function(x) continuous_summary(x[[name]]), character(1))
  )
}
for (name in names(binary)) {
  rows[[length(rows) + 1L]] <- c(
    Variable = unname(binary[name]),
    vapply(cohorts, function(x) binary_summary(x[[name]]), character(1))
  )
}
table1 <- rbindlist(lapply(rows, as.list), fill = TRUE)
setnames(table1, c("Variable", names(cohorts)))
fwrite(table1, file.path(output_root, "r_table1.csv"))

analytic <- windows[
  logical_true(primary_complete_window) & caseid %in% broad_ids
]
analytic[, map_region := cut(
  art_map, breaks = c(-Inf, 65, 75, Inf), right = FALSE,
  labels = c("MAP <65", "MAP 65-<75", "MAP >=75")
)]
map_support <- analytic[, .(
  windows = .N,
  surgeries = uniqueN(caseid),
  patients = uniqueN(subjectid),
  surgeries_with_2plus_windows = sum(.N >= 2L),
  median_windows_per_contributing_surgery = median(.N)
), by = map_region]
# Recompute surgery-level support correctly rather than using the grouped row count above.
map_support <- rbindlist(lapply(levels(analytic$map_region), function(region) {
  local <- analytic[map_region == region]
  per_case <- local[, .N, by = caseid]
  data.table(
    map_region = region,
    windows = nrow(local), surgeries = uniqueN(local$caseid),
    patients = uniqueN(local$subjectid),
    surgeries_with_2plus_windows = sum(per_case$N >= 2L),
    median_windows_per_contributing_surgery = median(per_case$N)
  )
}))
fwrite(map_support, file.path(output_root, "r_map_region_support.csv"))

candidate <- windows[caseid %in% broad_ids]
eligible <- logical_true(candidate$primary_complete_window)
requirements <- list(
  eeg_signal_quality = logical_true(candidate$eeg_signal_eligible),
  primary_ra = logical_true(candidate$primary_ra_eligible),
  art_map_quality = logical_true(candidate$art_qc),
  broad_no_positive_volatile = logical_true(candidate$broad_no_positive_volatile),
  propofol_ce_finite = is.finite(candidate$propofol_ce_reconstructed),
  remifentanil_ce_finite = is.finite(candidate$remifentanil_ce_reconstructed),
  propofol_rate_coverage_80 = candidate$propofol_rate_observed_fraction >= 0.8,
  remifentanil_rate_coverage_80 = candidate$remifentanil_rate_observed_fraction >= 0.8,
  etco2_finite = is.finite(candidate$etco2),
  temperature_finite = is.finite(candidate$body_temperature)
)
failures <- lapply(requirements, function(x) list(
  failed_windows = sum(!x), failed_fraction = mean(!x)
))

results <- list(
  implementation = R.version.string,
  selection_flow = selection_flow,
  table1 = list(cohorts = lapply(cohorts, function(x) list(
    surgeries = nrow(x), patients = uniqueN(x$subjectid)
  ))),
  map_support = list(
    analytic_windows = nrow(analytic),
    analytic_surgeries = uniqueN(analytic$caseid),
    analytic_patients = uniqueN(analytic$subjectid),
    low_map_windows = sum(analytic$art_map < 65),
    low_map_window_fraction = mean(analytic$art_map < 65),
    surgeries_with_any_low_map = uniqueN(analytic[art_map < 65]$caseid)
  ),
  missingness = list(
    candidate_windows_in_primary_analysis_surgeries = nrow(candidate),
    eligible_windows = sum(eligible),
    ineligible_windows = sum(!eligible),
    eligible_fraction = mean(eligible),
    overlapping_requirement_failures = failures
  )
)

write_json_pretty(results, file.path(output_root, "r_table1_missingness_results.json"))
print(results[c("selection_flow", "map_support", "missingness")])
