# Purpose: fit the reported random-intercept-only sensitivity models.
source("R/common.R")

message("[1/5] Reading analysis windows")
windows <- fread(windows_path)
bundle <- build_main_design(windows, minimum_windows = 8L)

message("[2/5] Fitting baseline CAR(1) random-intercept model")
fit0 <- fit_car1_lme(bundle, bundle$baseline_names)

message("[3/5] Fitting linear MAP model")
fit_linear <- fit_car1_lme(bundle, bundle$linear_names)

message("[4/5] Fitting nonlinear MAP model")
fit1 <- fit_car1_lme(bundle, bundle$nonlinear_names)

pred0 <- fixed_prediction(fit0, bundle, bundle$baseline_names)
pred_linear <- fixed_prediction(fit_linear, bundle, bundle$linear_names)
pred1 <- fixed_prediction(fit1, bundle, bundle$nonlinear_names)

residual_reduction <- case_equal_within_reduction(bundle$data, bundle$y, pred0, pred1)
window_weighted <- window_weighted_within_reduction(bundle$data, bundle$y, pred0, pred1)

lrt_overall <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))), 0)
df_overall <- length(bundle$nonlinear_names) - length(bundle$baseline_names)
lrt_nonlinear <- max(2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit_linear))), 0)
df_nonlinear <- length(bundle$nonlinear_names) - length(bundle$linear_names)

contrast <- map_contrast_matrix(bundle)
curve <- as.numeric(contrast %*% fixef(fit1)[bundle$nonlinear_names])

r2_0 <- mixed_r2(fit0, pred0)
r2_1 <- mixed_r2(fit1, pred1)

results <- list(
  implementation = list(
    language = R.version.string,
    engine = "nlme::lme",
    method = "maximum likelihood",
    random_effect = "surgery random intercept",
    residual_correlation = "nlme::corCAR1 with five-minute time unit"
  ),
  cohort = list(
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data),
    minimum_windows = min(bundle$data[, .N, by = caseid]$N),
    median_windows = median(bundle$data[, .N, by = caseid]$N)
  ),
  knots = bundle$knots,
  baseline = fit_components(fit0),
  linear = fit_components(fit_linear),
  nonlinear = fit_components(fit1),
  tests = list(
    overall_map = list(statistic = lrt_overall, df = df_overall, p_value = pchisq(lrt_overall, df_overall, lower.tail = FALSE)),
    nonlinearity = list(statistic = lrt_nonlinear, df = df_nonlinear, p_value = pchisq(lrt_nonlinear, df_nonlinear, lower.tail = FALSE))
  ),
  residual_reduction = list(
    case_equal_q0 = residual_reduction$q0,
    case_equal_q1 = residual_reduction$q1,
    case_equal_G = residual_reduction$G,
    window_weighted_G = window_weighted
  ),
  mixed_r2 = list(
    baseline = r2_0,
    nonlinear = r2_1,
    delta_marginal_r2 = r2_1$marginal_r2 - r2_0$marginal_r2,
    delta_conditional_r2 = r2_1$conditional_r2 - r2_0$conditional_r2
  ),
  map_curve = list(
    grid_mmhg = as.list(50:90),
    contrast_logit_ra_vs_65 = as.list(curve)
  )
)

message("[5/5] Writing R results and fitted objects")
write_json_pretty(results, file.path(output_root, "r_main_model_results.json"))
saveRDS(
  list(bundle = bundle, baseline = fit0, linear = fit_linear, nonlinear = fit1),
  file.path(output_root, "r_main_model_objects.rds"),
  compress = "xz"
)
fwrite(
  data.table(map_mmhg = 50:90, contrast_logit_ra_vs_65 = curve),
  file.path(output_root, "r_main_map_curve.csv")
)

print(results[c("cohort", "baseline", "linear", "nonlinear", "tests", "residual_reduction", "mixed_r2")])
