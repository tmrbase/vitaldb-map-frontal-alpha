options(stringsAsFactors = FALSE)

# Shared model, design-matrix, standardisation, diagnostics, and serialisation helpers.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(jsonlite)
})

repo_root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
output_root <- file.path(repo_root, "outputs", "maintenance_ra_r_reanalysis_v1")
release_asset_root <- Sys.getenv(
  "VITALDB_RELEASE_ASSETS",
  unset = file.path(repo_root, "outputs", "release_assets")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(release_asset_root, recursive = TRUE, showWarnings = FALSE)

windows_path <- file.path(
  repo_root,
  "outputs",
  "maintenance_ra_full_v1",
  "full_windows.csv"
)

natural_spline_basis <- function(x, knots) {
  x <- as.numeric(x)
  knots <- as.numeric(knots)
  stopifnot(length(knots) >= 3L, all(diff(knots) > 0))
  as.matrix(splines::ns(
    x,
    knots = knots[2L:(length(knots) - 1L)],
    Boundary.knots = knots[c(1L, length(knots))],
    intercept = FALSE
  ))
}

strict_knots <- function(x) {
  finite <- x[is.finite(x)]
  knots <- as.numeric(quantile(finite, probs = c(0.05, 0.35, 0.65, 0.95), type = 7))
  if (!all(diff(knots) > 0)) {
    unique_values <- sort(unique(finite))
    if (length(unique_values) < 4L) stop("Fewer than four unique values")
    knots <- as.numeric(quantile(unique_values, probs = c(0.05, 0.35, 0.65, 0.95), type = 7))
  }
  if (!all(diff(knots) > 0)) stop("Could not obtain strictly increasing knots")
  knots
}

standardize_column <- function(x) {
  center <- mean(x)
  scale <- sd(x)
  if (!is.finite(scale) || scale <= 1e-12) stop("Unusable design column")
  list(value = (x - center) / scale, mean = center, sd = scale)
}

logical_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(as.character(x)) %chin% c("true", "t", "1")
}

build_main_design <- function(
  windows,
  strict_measured_zero = FALSE,
  minimum_windows = 8L,
  fixed_knots = NULL
) {
  data <- copy(windows[logical_true(primary_complete_window)])
  if (strict_measured_zero) {
    data <- data[logical_true(strict_measured_zero_volatile)]
  }
  counts <- data[, .N, by = caseid]
  keep <- counts[N >= as.integer(minimum_windows), caseid]
  data <- data[caseid %in% keep]
  setorder(data, subjectid, caseid, time_seconds)
  if (!nrow(data)) stop("No complete surgery has enough windows")

  design <- data.table(intercept = rep(1, nrow(data)))
  scaling <- list()
  knots <- list()

  add_standardized <- function(name, raw) {
    item <- standardize_column(raw)
    design[[name]] <<- item$value
    scaling[[name]] <<- list(mean = item$mean, standard_deviation = item$sd)
  }

  spline_sources <- list(
    time = "minutes_since_opstart",
    propofol_ce = "propofol_ce_reconstructed",
    remifentanil_ce = "remifentanil_ce_reconstructed"
  )
  for (prefix in names(spline_sources)) {
    source <- spline_sources[[prefix]]
    local_knots <- if (is.null(fixed_knots)) strict_knots(data[[source]]) else fixed_knots[[prefix]]
    knots[[prefix]] <- as.numeric(local_knots)
    basis <- natural_spline_basis(data[[source]], local_knots)
    for (j in seq_len(ncol(basis))) {
      raw <- basis[, j]
      between <- ave(raw, data$caseid, FUN = mean)
      within <- raw - between
      add_standardized(sprintf("%s_basis%d_between", prefix, j), between)
      add_standardized(sprintf("%s_basis%d_within", prefix, j), within)
    }
  }

  linear_sources <- list(etco2 = "etco2", temperature = "body_temperature")
  for (prefix in names(linear_sources)) {
    source <- linear_sources[[prefix]]
    raw <- as.numeric(data[[source]])
    between <- ave(raw, data$caseid, FUN = mean)
    within <- raw - between
    add_standardized(sprintf("%s_between", prefix), between)
    add_standardized(sprintf("%s_within", prefix), within)
  }

  age_knots <- if (is.null(fixed_knots)) strict_knots(data$age) else fixed_knots$age
  knots$age <- as.numeric(age_knots)
  age_basis <- natural_spline_basis(data$age, age_knots)
  for (j in seq_len(ncol(age_basis))) {
    add_standardized(sprintf("age_basis%d", j), age_basis[, j])
  }

  map_knots <- if (is.null(fixed_knots)) strict_knots(data$art_map) else fixed_knots$map
  knots$map <- as.numeric(map_knots)
  map_basis <- natural_spline_basis(data$art_map, map_knots)
  for (j in seq_len(ncol(map_basis))) {
    raw <- map_basis[, j]
    between <- ave(raw, data$caseid, FUN = mean)
    within <- raw - between
    add_standardized(sprintf("map_basis%d_between", j), between)
    add_standardized(sprintf("map_basis%d_within", j), within)
  }
  raw_map_within <- (data$art_map - ave(data$art_map, data$caseid, FUN = mean)) / 10
  add_standardized("map_linear_within_10", raw_map_within)

  design$sex_male <- as.numeric(toupper(as.character(data$sex)) == "M")
  all_names <- names(design)
  common_names <- all_names[
    !startsWith(all_names, "map_basis") & all_names != "map_linear_within_10"
  ]
  map_between_names <- all_names[
    startsWith(all_names, "map_basis") & endsWith(all_names, "_between")
  ]
  map_within_names <- all_names[
    startsWith(all_names, "map_basis") & endsWith(all_names, "_within")
  ]
  baseline_names <- c(common_names, map_between_names)
  linear_names <- c(baseline_names, "map_linear_within_10")
  nonlinear_names <- c(baseline_names, map_within_names)

  list(
    data = data,
    y = as.numeric(data$logit_relative_alpha),
    design = design,
    baseline_names = baseline_names,
    linear_names = linear_names,
    nonlinear_names = nonlinear_names,
    knots = knots,
    scaling = scaling
  )
}

lme_fit_diagnostics <- function(fit) {
  fixed_covariance <- tryCatch(vcov(fit), error = function(e) e)
  random_covariance <- tryCatch(
    as.matrix(getVarCov(fit, type = "random.effects")),
    error = function(e) e
  )
  rho <- tryCatch(
    as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)),
    error = function(e) NA_real_
  )
  warnings <- attr(fit, "fit_warnings", exact = TRUE)
  if (is.null(warnings)) warnings <- character()
  optimizer_issue <- any(grepl(
    "convergence|iteration limit|false convergence|singular convergence",
    warnings,
    ignore.case = TRUE
  ))
  fixed_covariance_ok <- !inherits(fixed_covariance, "error") &&
    all(is.finite(fixed_covariance))
  random_covariance_finite <- !inherits(random_covariance, "error") &&
    all(is.finite(random_covariance))
  minimum_random_eigenvalue <- if (random_covariance_finite) {
    min(eigen(random_covariance, symmetric = TRUE, only.values = TRUE)$values)
  } else {
    NA_real_
  }
  random_covariance_ok <- random_covariance_finite &&
    minimum_random_eigenvalue >= -1e-10
  essential_ok <- is.finite(as.numeric(logLik(fit))) &&
    all(is.finite(fixef(fit))) &&
    fixed_covariance_ok && random_covariance_ok &&
    is.finite(sigma(fit)) && sigma(fit) > 0 &&
    is.finite(rho) && rho >= 0 && rho < 1 &&
    !optimizer_issue
  list(
    essential_checks_passed = essential_ok,
    optimiser = attr(fit, "optimizer", exact = TRUE) %||% NA_character_,
    optimiser_warnings = as.list(warnings),
    log_likelihood = as.numeric(logLik(fit)),
    residual_standard_deviation = sigma(fit),
    car1_phi = rho,
    fixed_effect_covariance_finite = fixed_covariance_ok,
    random_effect_covariance_positive_semidefinite = random_covariance_ok,
    minimum_random_effect_covariance_eigenvalue = minimum_random_eigenvalue,
    variance_parameter_covariance_available = !is.character(fit$apVar),
    variance_parameter_covariance_message = if (is.character(fit$apVar)) fit$apVar else NULL
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

fit_car1_lme <- function(
  bundle,
  fixed_names,
  outcome = NULL,
  rho_starts = c(0.1, 0.3, 0.7),
  optimizer = "optim"
) {
  fit_data <- cbind(
    data.table(
      y = if (is.null(outcome)) bundle$y else as.numeric(outcome),
      caseid = factor(bundle$data$caseid),
      subjectid = factor(bundle$data$subjectid),
      time5 = as.numeric(bundle$data$time_seconds) / 300
    ),
    bundle$design[, fixed_names, with = FALSE]
  )
  formula_fixed <- as.formula(sprintf("y ~ -1 + %s", paste(fixed_names, collapse = " + ")))
  candidates <- list()
  errors <- character()
  for (rho in rho_starts) {
    warning_messages <- character()
    fit <- tryCatch(
      withCallingHandlers(lme(
        fixed = formula_fixed,
        random = ~1 | caseid,
        correlation = corCAR1(value = rho, form = ~time5 | caseid),
        data = fit_data,
        method = "ML",
        na.action = na.fail,
        control = lmeControl(
          opt = optimizer,
          msMaxIter = 300,
          msMaxEval = 1000,
          tolerance = 1e-7,
          niterEM = 25,
          returnObject = TRUE
        )
      ), warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      errors <- c(errors, sprintf("rho=%s: %s", rho, conditionMessage(fit)))
    } else {
      attr(fit, "fit_warnings") <- unique(warning_messages)
      attr(fit, "optimizer") <- optimizer
      candidates[[length(candidates) + 1L]] <- fit
    }
  }
  if (!length(candidates)) stop(paste(errors, collapse = "\n"))
  diagnostics <- lapply(candidates, lme_fit_diagnostics)
  valid <- which(vapply(
    diagnostics,
    function(x) isTRUE(x$essential_checks_passed),
    logical(1)
  ))
  if (!length(valid)) stop("All nlme candidates failed essential convergence checks")
  likelihoods <- vapply(candidates[valid], function(x) as.numeric(logLik(x)), numeric(1))
  selected <- valid[[which.max(likelihoods)]]
  attr(candidates[[selected]], "fit_diagnostics") <- diagnostics[[selected]]
  candidates[[selected]]
}

fit_car1_random_map_slope_lme <- function(
  bundle,
  fixed_names,
  outcome = NULL,
  rho_starts = c(0.1, 0.3, 0.7),
  optimizer = "optim"
) {
  fit_data <- cbind(
    data.table(
      y = if (is.null(outcome)) bundle$y else as.numeric(outcome),
      caseid = factor(bundle$data$caseid),
      subjectid = factor(bundle$data$subjectid),
      time5 = as.numeric(bundle$data$time_seconds) / 300,
      random_map10 = (
        as.numeric(bundle$data$art_map) -
          ave(as.numeric(bundle$data$art_map), bundle$data$caseid, FUN = mean)
      ) / 10
    ),
    bundle$design[, fixed_names, with = FALSE]
  )
  formula_fixed <- as.formula(sprintf("y ~ -1 + %s", paste(fixed_names, collapse = " + ")))
  candidates <- list()
  errors <- character()
  for (rho in rho_starts) {
    warning_messages <- character()
    fit <- tryCatch(
      withCallingHandlers(lme(
        fixed = formula_fixed,
        random = ~1 + random_map10 | caseid,
        correlation = corCAR1(value = rho, form = ~time5 | caseid),
        data = fit_data,
        method = "ML",
        na.action = na.fail,
        control = lmeControl(
          opt = optimizer,
          msMaxIter = 500,
          msMaxEval = 2000,
          tolerance = 1e-7,
          niterEM = 50,
          returnObject = TRUE
        )
      ), warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      errors <- c(errors, sprintf("rho=%s: %s", rho, conditionMessage(fit)))
    } else {
      attr(fit, "fit_warnings") <- unique(warning_messages)
      attr(fit, "optimizer") <- optimizer
      candidates[[length(candidates) + 1L]] <- fit
    }
  }
  if (!length(candidates)) stop(paste(errors, collapse = "\n"))
  diagnostics <- lapply(candidates, lme_fit_diagnostics)
  valid <- which(vapply(
    diagnostics,
    function(x) isTRUE(x$essential_checks_passed),
    logical(1)
  ))
  if (!length(valid)) stop("All nlme candidates failed essential convergence checks")
  likelihoods <- vapply(candidates[valid], function(x) as.numeric(logLik(x)), numeric(1))
  selected <- valid[[which.max(likelihoods)]]
  attr(candidates[[selected]], "fit_diagnostics") <- diagnostics[[selected]]
  candidates[[selected]]
}

random_slope_components <- function(fit, random_map10) {
  covariance <- as.matrix(getVarCov(fit, type = "random.effects"))
  random_variance <-
    covariance[1, 1] +
    2 * covariance[1, 2] * random_map10 +
    covariance[2, 2] * random_map10^2
  list(
    random_intercept_variance = covariance[1, 1],
    random_slope_variance_per_10_mmhg = covariance[2, 2],
    random_slope_sd_per_10_mmhg = sqrt(covariance[2, 2]),
    intercept_slope_correlation = covariance[1, 2] /
      sqrt(covariance[1, 1] * covariance[2, 2]),
    mean_random_effect_variance = mean(random_variance),
    residual_variance = sigma(fit)^2,
    correlation_at_five_minutes = as.numeric(
      coef(fit$modelStruct$corStruct, unconstrained = FALSE)
    ),
    log_likelihood = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit)
  )
}

mixed_r2_random_slope <- function(fit, prediction, random_map10) {
  component <- random_slope_components(fit, random_map10)
  fixed_variance <- var(as.numeric(prediction)) * (length(prediction) - 1) / length(prediction)
  total <- fixed_variance + component$mean_random_effect_variance + component$residual_variance
  list(
    fixed_effect_variance = fixed_variance,
    mean_random_effect_variance = component$mean_random_effect_variance,
    residual_variance = component$residual_variance,
    marginal_r2 = fixed_variance / total,
    conditional_r2 = (fixed_variance + component$mean_random_effect_variance) / total
  )
}

fixed_prediction <- function(fit, bundle, fixed_names) {
  as.numeric(
    as.matrix(bundle$design[, fixed_names, with = FALSE]) %*%
      fixef(fit)[fixed_names]
  )
}

case_equal_within_reduction <- function(data, y, pred0, pred1) {
  frame <- data.table(
    caseid = data$caseid,
    r0 = as.numeric(y) - as.numeric(pred0),
    r1 = as.numeric(y) - as.numeric(pred1)
  )
  frame[, r0 := r0 - mean(r0), by = caseid]
  frame[, r1 := r1 - mean(r1), by = caseid]
  local <- frame[, .(q0 = mean(r0^2), q1 = mean(r1^2)), by = caseid]
  q0 <- mean(local$q0)
  q1 <- mean(local$q1)
  list(q0 = q0, q1 = q1, G = 1 - q1 / q0)
}

window_weighted_within_reduction <- function(data, y, pred0, pred1) {
  frame <- data.table(
    caseid = data$caseid,
    r0 = as.numeric(y) - as.numeric(pred0),
    r1 = as.numeric(y) - as.numeric(pred1)
  )
  frame[, r0 := r0 - mean(r0), by = caseid]
  frame[, r1 := r1 - mean(r1), by = caseid]
  1 - sum(frame$r1^2) / sum(frame$r0^2)
}

fit_components <- function(fit) {
  variance_table <- VarCorr(fit)
  list(
    log_likelihood = as.numeric(logLik(fit)),
    random_intercept_variance = as.numeric(variance_table[1, "Variance"]),
    residual_variance = sigma(fit)^2,
    correlation_at_five_minutes = as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)),
    AIC = AIC(fit),
    BIC = BIC(fit)
  )
}

mixed_r2 <- function(fit, prediction) {
  component <- fit_components(fit)
  fixed_variance <- var(as.numeric(prediction)) * (length(prediction) - 1) / length(prediction)
  total <- fixed_variance + component$random_intercept_variance + component$residual_variance
  list(
    fixed_effect_variance = fixed_variance,
    random_intercept_variance = component$random_intercept_variance,
    residual_variance = component$residual_variance,
    marginal_r2 = fixed_variance / total,
    conditional_r2 = (fixed_variance + component$random_intercept_variance) / total
  )
}

map_contrast_matrix <- function(bundle, grid = 50:90, reference = 65) {
  grid_basis <- natural_spline_basis(grid, bundle$knots$map)
  ref_basis <- natural_spline_basis(reference, bundle$knots$map)[1, ]
  out <- matrix(0, nrow = length(grid), ncol = length(bundle$nonlinear_names))
  colnames(out) <- bundle$nonlinear_names
  for (j in seq_len(ncol(grid_basis))) {
    name <- sprintf("map_basis%d_within", j)
    out[, name] <- (grid_basis[, j] - ref_basis[j]) /
      bundle$scaling[[name]]$standard_deviation
  }
  out
}

write_json_pretty <- function(x, path) {
  write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = 16, null = "null")
}

augment_bundle <- function(bundle, covariates) {
  result <- bundle
  for (item in covariates) {
    source <- item[[1L]]
    prefix <- item[[2L]]
    values <- as.numeric(result$data[[source]])
    between <- ave(values, result$data$caseid, FUN = mean)
    within <- values - between
    for (part in c("between", "within")) {
      raw <- if (part == "between") between else within
      standardized <- standardize_column(raw)
      name <- sprintf("%s_%s", prefix, part)
      result$design[[name]] <- standardized$value
      result$scaling[[name]] <- list(
        mean = standardized$mean,
        standard_deviation = standardized$sd
      )
      result$baseline_names <- c(result$baseline_names, name)
      result$linear_names <- c(result$linear_names, name)
      result$nonlinear_names <- c(result$nonlinear_names, name)
    }
  }
  result
}

fit_standard_model_set <- function(bundle, outcome = NULL) {
  list(
    baseline = fit_car1_random_map_slope_lme(bundle, bundle$baseline_names, outcome = outcome),
    linear = fit_car1_random_map_slope_lme(bundle, bundle$linear_names, outcome = outcome),
    nonlinear = fit_car1_random_map_slope_lme(bundle, bundle$nonlinear_names, outcome = outcome)
  )
}

summarize_standard_model_set <- function(bundle, fits, label, outcome = NULL) {
  y <- if (is.null(outcome)) bundle$y else as.numeric(outcome)
  pred0 <- fixed_prediction(fits$baseline, bundle, bundle$baseline_names)
  pred1 <- fixed_prediction(fits$nonlinear, bundle, bundle$nonlinear_names)
  reduction <- case_equal_within_reduction(bundle$data, y, pred0, pred1)
  random_map10 <- (
    as.numeric(bundle$data$art_map) -
      ave(as.numeric(bundle$data$art_map), bundle$data$caseid, FUN = mean)
  ) / 10
  overall <- max(2 * (as.numeric(logLik(fits$nonlinear)) - as.numeric(logLik(fits$baseline))), 0)
  nonlinear <- max(2 * (as.numeric(logLik(fits$nonlinear)) - as.numeric(logLik(fits$linear))), 0)
  list(
    label = label,
    subjects = uniqueN(bundle$data$subjectid),
    surgeries = uniqueN(bundle$data$caseid),
    windows = nrow(bundle$data),
    overall_map_likelihood_ratio = overall,
    overall_map_df = length(bundle$nonlinear_names) - length(bundle$baseline_names),
    overall_map_p_value = pchisq(
      overall,
      length(bundle$nonlinear_names) - length(bundle$baseline_names),
      lower.tail = FALSE
    ),
    nonlinearity_likelihood_ratio = nonlinear,
    nonlinearity_df = length(bundle$nonlinear_names) - length(bundle$linear_names),
    nonlinearity_p_value = pchisq(
      nonlinear,
      length(bundle$nonlinear_names) - length(bundle$linear_names),
      lower.tail = FALSE
    ),
    case_equal_q0 = reduction$q0,
    case_equal_q1 = reduction$q1,
    case_equal_G = reduction$G,
    window_weighted_G = window_weighted_within_reduction(
      bundle$data, y, pred0, pred1
    ),
    baseline = random_slope_components(fits$baseline, random_map10),
    nonlinear = random_slope_components(fits$nonlinear, random_map10)
  )
}
