# Purpose: fit the authoritative B0/B1 mixed models with surgery-specific MAP slopes.
source("R/common.R")

message("[1/6] Reading analysis windows")
windows <- fread(windows_path)
bundle <- build_main_design(windows, minimum_windows = 8L)
random_map10 <- (
  as.numeric(bundle$data$art_map) -
    ave(as.numeric(bundle$data$art_map), bundle$data$caseid, FUN = mean)
) / 10

message("[2/6] Fitting baseline model with surgery-specific MAP slope")
fit0 <- fit_car1_random_map_slope_lme(bundle, bundle$baseline_names)

message("[3/6] Fitting linear MAP model with surgery-specific MAP slope")
fit_linear <- fit_car1_random_map_slope_lme(bundle, bundle$linear_names)

message("[4/6] Fitting nonlinear MAP model with surgery-specific MAP slope")
fit1 <- fit_car1_random_map_slope_lme(bundle, bundle$nonlinear_names)

pred0 <- fixed_prediction(fit0, bundle, bundle$baseline_names)
pred_linear <- fixed_prediction(fit_linear, bundle, bundle$linear_names)
pred1 <- fixed_prediction(fit1, bundle, bundle$nonlinear_names)

reduction <- case_equal_within_reduction(bundle$data, bundle$y, pred0, pred1)
window_weighted <- window_weighted_within_reduction(bundle$data, bundle$y, pred0, pred1)

lrt_overall <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))), 0)
df_overall <- length(bundle$nonlinear_names) - length(bundle$baseline_names)
lrt_nonlinear <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit_linear))), 0)
df_nonlinear <- length(bundle$nonlinear_names) - length(bundle$linear_names)

contrast <- map_contrast_matrix(bundle, grid = 50:110, reference = 65)
curve <- as.numeric(contrast %*% fixef(fit1)[bundle$nonlinear_names])

r2_0 <- mixed_r2_random_slope(fit0, pred0, random_map10)
r2_1 <- mixed_r2_random_slope(fit1, pred1, random_map10)

results <- list(
  implementation = list(
    language = R.version.string,
    engine = "nlme::lme",
    method = "maximum likelihood",
    spline_engine = "splines::ns with two internal and two boundary knots",
    random_effect = "surgery random intercept and surgery-specific linear within-surgery MAP slope",
    residual_correlation = "nlme::corCAR1 with five-minute time unit"
  ),
  cohort = list(
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data)
  ),
  knots = bundle$knots,
  baseline = random_slope_components(fit0, random_map10),
  linear = random_slope_components(fit_linear, random_map10),
  nonlinear = random_slope_components(fit1, random_map10),
  convergence_diagnostics = list(
    baseline = lme_fit_diagnostics(fit0),
    linear = lme_fit_diagnostics(fit_linear),
    nonlinear = lme_fit_diagnostics(fit1)
  ),
  tests = list(
    overall_map = list(
      statistic = lrt_overall,
      df = df_overall,
      p_value = pchisq(lrt_overall, df_overall, lower.tail = FALSE)
    ),
    nonlinearity = list(
      statistic = lrt_nonlinear,
      df = df_nonlinear,
      p_value = pchisq(lrt_nonlinear, df_nonlinear, lower.tail = FALSE)
    )
  ),
  residual_reduction = list(
    case_equal_q0 = reduction$q0,
    case_equal_q1 = reduction$q1,
    case_equal_G = reduction$G,
    window_weighted_G = window_weighted
  ),
  mixed_r2 = list(
    baseline = r2_0,
    nonlinear = r2_1,
    delta_marginal_r2 = r2_1$marginal_r2 - r2_0$marginal_r2,
    delta_conditional_r2 = r2_1$conditional_r2 - r2_0$conditional_r2
  ),
  map_curve = list(
    grid_mmhg = as.list(50:110),
    contrast_logit_ra_vs_65 = as.list(curve)
  )
)

write_json_pretty(
  results,
  file.path(output_root, "r_random_slope_primary_results.json")
)
saveRDS(
  list(bundle = bundle, baseline = fit0, linear = fit_linear, nonlinear = fit1),
  file.path(output_root, "r_random_slope_primary_objects.rds"),
  compress = "xz"
)
fwrite(
  data.table(map_mmhg = 50:110, contrast_logit_ra_vs_65 = curve),
  file.path(output_root, "r_random_slope_primary_map_curve.csv")
)

print(results[c("cohort", "baseline", "linear", "nonlinear", "tests", "residual_reduction", "mixed_r2")])
