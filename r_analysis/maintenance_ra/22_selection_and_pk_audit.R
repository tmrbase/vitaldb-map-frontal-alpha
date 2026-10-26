# Purpose: describe cohort selection and compare reconstructed with device-recorded drug concentrations.
source("R/common.R")

input_root <- file.path(repo_root, "outputs", "maintenance_ra_full_v1")

cases <- fread(file.path(input_root, "eligible_cases.csv"))
windows <- fread(file.path(input_root, "full_windows.csv"))
complete_counts <- windows[logical_true(primary_complete_window), .N, by = caseid]
analytic_caseids <- complete_counts[N >= 8L, caseid]
cases[, analytic_group := fifelse(caseid %in% analytic_caseids, "Included", "Excluded: <8 complete windows")]
cases[, `:=`(
  operation_minutes = (opend - opstart) / 60,
  thoracic_department = as.integer(grepl("thoracic", department, ignore.case = TRUE)),
  asa_ge3 = as.integer(asa >= 3),
  male = as.integer(toupper(as.character(sex)) == "M")
)]

median_iqr_text <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

continuous_row <- function(variable, label) {
  included <- cases[analytic_group == "Included"][[variable]]
  excluded <- cases[analytic_group != "Included"][[variable]]
  finite_in <- included[is.finite(included)]
  finite_ex <- excluded[is.finite(excluded)]
  pooled_sd <- sqrt((var(finite_in) + var(finite_ex)) / 2)
  data.table(
    Characteristic = label,
    Included = median_iqr_text(finite_in),
    Included_missing = sum(!is.finite(included)),
    Excluded = median_iqr_text(finite_ex),
    Excluded_missing = sum(!is.finite(excluded)),
    Standardized_mean_difference = ifelse(pooled_sd > 0, (mean(finite_in) - mean(finite_ex)) / pooled_sd, NA_real_)
  )
}

binary_row <- function(variable, label) {
  included <- cases[analytic_group == "Included"][[variable]]
  excluded <- cases[analytic_group != "Included"][[variable]]
  finite_in <- included[is.finite(included)]
  finite_ex <- excluded[is.finite(excluded)]
  p_in <- mean(finite_in > 0)
  p_ex <- mean(finite_ex > 0)
  pooled_sd <- sqrt((p_in * (1 - p_in) + p_ex * (1 - p_ex)) / 2)
  data.table(
    Characteristic = label,
    Included = sprintf("%d (%.1f%%)", sum(finite_in > 0), 100 * p_in),
    Included_missing = sum(!is.finite(included)),
    Excluded = sprintf("%d (%.1f%%)", sum(finite_ex > 0), 100 * p_ex),
    Excluded_missing = sum(!is.finite(excluded)),
    Standardized_mean_difference = ifelse(pooled_sd > 0, (p_in - p_ex) / pooled_sd, NA_real_)
  )
}

selection_table <- rbindlist(list(
  continuous_row("age", "Age, yr"),
  binary_row("male", "Male sex"),
  continuous_row("bmi", "Body mass index, kg/m2"),
  binary_row("asa_ge3", "ASA physical status >=3"),
  binary_row("emop", "Emergency surgery"),
  binary_row("preop_htn", "Preoperative hypertension"),
  binary_row("preop_dm", "Preoperative diabetes mellitus"),
  continuous_row("preop_hb", "Preoperative haemoglobin, g/dL"),
  continuous_row("preop_alb", "Preoperative albumin, g/dL"),
  continuous_row("preop_cr", "Preoperative creatinine, mg/dL"),
  continuous_row("preop_gluc", "Preoperative glucose, mg/dL"),
  binary_row("thoracic_department", "Thoracic surgery department"),
  continuous_row("operation_minutes", "Operation duration, min"),
  continuous_row("intraop_ebl", "Estimated blood loss, mL"),
  continuous_row("intraop_crystalloid", "Crystalloid administered, mL")
))
selection_table[, `:=`(
  Included_surgeries = sum(cases$analytic_group == "Included"),
  Excluded_surgeries = sum(cases$analytic_group != "Included")
)]
fwrite(selection_table, file.path(output_root, "r_selection_included_vs_excluded.csv"))

audit_drug <- function(label, reconstructed_name, device_name, unit) {
  frame <- windows[
    caseid %in% analytic_caseids & logical_true(primary_complete_window) &
      is.finite(get(reconstructed_name)) & is.finite(get(device_name)) &
      get(reconstructed_name) >= 0 & get(device_name) >= 0,
    .(
      caseid, subjectid,
      reconstructed = get(reconstructed_name),
      device = get(device_name)
    )
  ]
  difference <- frame$reconstructed - frame$device
  pearson <- cor(frame$reconstructed, frame$device, method = "pearson")
  spearman <- cor(frame$reconstructed, frame$device, method = "spearman")
  concordance <- 2 * cov(frame$reconstructed, frame$device) /
    (var(frame$reconstructed) + var(frame$device) + (mean(frame$reconstructed) - mean(frame$device))^2)
  bias <- mean(difference)
  data.table(
    Drug = label,
    Unit = unit,
    Patients = uniqueN(frame$subjectid),
    Surgeries = uniqueN(frame$caseid),
    Windows = nrow(frame),
    Pearson_r = pearson,
    Spearman_rho = spearman,
    Concordance_correlation = concordance,
    Median_difference = median(difference),
    Difference_Q1 = quantile(difference, 0.25),
    Difference_Q3 = quantile(difference, 0.75),
    Mean_difference = bias,
    Lower_95_agreement = bias - 1.96 * sd(difference),
    Upper_95_agreement = bias + 1.96 * sd(difference),
    Mean_absolute_error = mean(abs(difference)),
    Root_mean_squared_error = sqrt(mean(difference^2))
  )
}

pk_audit <- rbindlist(list(
  audit_drug("Propofol", "propofol_ce_reconstructed", "propofol_ce_device", "microgram/mL"),
  audit_drug("Remifentanil", "remifentanil_ce_reconstructed", "remifentanil_ce_device", "ng/mL")
))
fwrite(pk_audit, file.path(output_root, "r_pk_device_audit.csv"))

write_json(list(
  selection = list(
    candidate_surgeries = nrow(cases),
    included_surgeries = sum(cases$analytic_group == "Included"),
    excluded_surgeries = sum(cases$analytic_group != "Included")
  ),
  pk_audit = pk_audit
), file.path(output_root, "r_selection_and_pk_audit_results.json"))

message("Selection comparison and PK audit complete")
