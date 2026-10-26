"""Patient-deletion influence diagnostic for Methods S5.

Only ranks patients for Table SIII's sensitivity analysis. The primary and
post-exclusion mixed models use nlme in R. The original spline basis and GLS
diagnostic are retained to reproduce the ranking used in the reported result.
See docs/METHODS_CODE_MAP.md for provenance and validation boundaries.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.special import expit, logit

KNOT_PROBABILITIES = (0.05, 0.35, 0.65, 0.95)


@dataclass(frozen=True)
class Car1Fit:
    beta: np.ndarray
    covariance_beta: np.ndarray
    log_likelihood: float
    case_random_intercept_variance: float
    residual_variance: float
    correlation_at_five_minutes: float
    converged: bool
    optimizer_message: str
    fixed_effect_names: tuple[str, ...]

    def predict_fixed(self, design: np.ndarray) -> np.ndarray:
        return np.asarray(design, dtype=float) @ self.beta


def restricted_cubic_spline_basis(
    values: np.ndarray, knots: np.ndarray
) -> np.ndarray:
    """Harrell-style restricted cubic spline: linear plus K-2 nonlinear columns."""

    x = np.asarray(values, dtype=float)
    k = np.asarray(knots, dtype=float)
    if k.ndim != 1 or k.size < 3 or not np.all(np.diff(k) > 0):
        raise ValueError(f"Spline knots must be strictly increasing: {k}")
    span = float(k[-1] - k[0])
    z = (x - k[0]) / span
    z_knots = (k - k[0]) / span

    def positive_cube(source: np.ndarray, threshold: float) -> np.ndarray:
        return np.maximum(source - threshold, 0.0) ** 3

    columns = [z]
    penultimate = float(z_knots[-2])
    last = float(z_knots[-1])
    denominator = last - penultimate
    for current in z_knots[:-2]:
        nonlinear = (
            positive_cube(z, float(current))
            - positive_cube(z, penultimate) * (last - float(current)) / denominator
            + positive_cube(z, last) * (penultimate - float(current)) / denominator
        )
        columns.append(nonlinear)
    return np.column_stack(columns)


def _strict_knots(values: pd.Series) -> np.ndarray:
    finite = np.asarray(values[np.isfinite(values)], dtype=float)
    knots = np.quantile(finite, KNOT_PROBABILITIES)
    if not np.all(np.diff(knots) > 0):
        unique = np.unique(finite)
        if unique.size < 4:
            raise ValueError("Fewer than four unique values are available for a spline")
        knots = np.quantile(unique, KNOT_PROBABILITIES)
    if not np.all(np.diff(knots) > 0):
        raise ValueError(f"Could not obtain strictly increasing spline knots: {knots}")
    return knots


def _append_standardized(
    destination: dict[str, np.ndarray],
    manifest: dict[str, dict[str, float]],
    name: str,
    values: np.ndarray,
) -> None:
    raw = np.asarray(values, dtype=float)
    mean = float(np.mean(raw))
    standard_deviation = float(np.std(raw, ddof=1))
    if not np.isfinite(standard_deviation) or standard_deviation <= 1e-12:
        raise ValueError(f"Design column {name} has no usable variation")
    destination[name] = (raw - mean) / standard_deviation
    manifest[name] = {"mean": mean, "standard_deviation": standard_deviation}


def build_main_design(
    windows: pd.DataFrame,
    *,
    strict_measured_zero: bool = False,
    minimum_windows: int = 8,
    fixed_knots: dict[str, list[float]] | None = None,
) -> dict[str, Any]:
    data = windows.loc[windows["primary_complete_window"].astype(bool)].copy()
    if strict_measured_zero:
        data = data.loc[data["strict_measured_zero_volatile"].astype(bool)].copy()
    counts = data.groupby("caseid", observed=True).size()
    data = data.loc[data["caseid"].isin(counts[counts >= int(minimum_windows)].index)].copy()
    data = data.sort_values(["subjectid", "caseid", "time_seconds"]).reset_index(drop=True)
    if data.empty:
        raise ValueError("No complete surgery has enough windows for the main model")

    continuous: dict[str, np.ndarray] = {}
    scaling: dict[str, dict[str, float]] = {}
    knots: dict[str, list[float]] = {}
    case_groups = data["caseid"]

    for source, prefix in (
        ("minutes_since_opstart", "time"),
        ("propofol_ce_reconstructed", "propofol_ce"),
        ("remifentanil_ce_reconstructed", "remifentanil_ce"),
    ):
        source_knots = (
            np.asarray(fixed_knots[prefix], dtype=float)
            if fixed_knots is not None
            else _strict_knots(data[source])
        )
        knots[prefix] = source_knots.tolist()
        basis = restricted_cubic_spline_basis(data[source].to_numpy(), source_knots)
        for column in range(basis.shape[1]):
            raw = pd.Series(basis[:, column], index=data.index)
            between = raw.groupby(case_groups, observed=True).transform("mean").to_numpy()
            within = raw.to_numpy() - between
            _append_standardized(
                continuous, scaling, f"{prefix}_basis{column + 1}_between", between
            )
            _append_standardized(
                continuous, scaling, f"{prefix}_basis{column + 1}_within", within
            )

    for source, prefix in (("etco2", "etco2"), ("body_temperature", "temperature")):
        raw = data[source].astype(float)
        between = raw.groupby(case_groups, observed=True).transform("mean").to_numpy()
        within = raw.to_numpy() - between
        _append_standardized(continuous, scaling, f"{prefix}_between", between)
        _append_standardized(continuous, scaling, f"{prefix}_within", within)

    age_knots = (
        np.asarray(fixed_knots["age"], dtype=float)
        if fixed_knots is not None
        else _strict_knots(data["age"])
    )
    knots["age"] = age_knots.tolist()
    age_basis = restricted_cubic_spline_basis(data["age"].to_numpy(), age_knots)
    for column in range(age_basis.shape[1]):
        _append_standardized(
            continuous, scaling, f"age_basis{column + 1}", age_basis[:, column]
        )

    map_knots = (
        np.asarray(fixed_knots["map"], dtype=float)
        if fixed_knots is not None
        else _strict_knots(data["art_map"])
    )
    knots["map"] = map_knots.tolist()
    map_basis = restricted_cubic_spline_basis(data["art_map"].to_numpy(), map_knots)
    for column in range(map_basis.shape[1]):
        raw = pd.Series(map_basis[:, column], index=data.index)
        between = raw.groupby(case_groups, observed=True).transform("mean").to_numpy()
        within = raw.to_numpy() - between
        _append_standardized(
            continuous, scaling, f"map_basis{column + 1}_between", between
        )
        _append_standardized(
            continuous, scaling, f"map_basis{column + 1}_within", within
        )
    raw_map_within = (
        data["art_map"]
        - data.groupby("caseid", observed=True)["art_map"].transform("mean")
    ).to_numpy() / 10.0
    _append_standardized(continuous, scaling, "map_linear_within_10", raw_map_within)

    columns = {"intercept": np.ones(len(data), dtype=float)}
    columns.update(continuous)
    columns["sex_male"] = data["sex"].astype(str).str.upper().eq("M").astype(float).to_numpy()
    design = pd.DataFrame(columns, index=data.index)
    common_names = [
        name for name in design.columns if not name.startswith("map_basis")
        and name != "map_linear_within_10"
    ]
    map_between_names = [
        name for name in design.columns if name.startswith("map_basis") and name.endswith("_between")
    ]
    map_within_names = [
        name for name in design.columns if name.startswith("map_basis") and name.endswith("_within")
    ]
    baseline_names = common_names + map_between_names
    linear_names = baseline_names + ["map_linear_within_10"]
    nonlinear_names = baseline_names + map_within_names
    return {
        "data": data,
        "y": data["logit_relative_alpha"].to_numpy(dtype=float),
        "baseline": design[baseline_names].to_numpy(dtype=float),
        "linear": design[linear_names].to_numpy(dtype=float),
        "nonlinear": design[nonlinear_names].to_numpy(dtype=float),
        "baseline_names": tuple(baseline_names),
        "linear_names": tuple(linear_names),
        "nonlinear_names": tuple(nonlinear_names),
        "knots": knots,
        "scaling": scaling,
    }


def _group_indices(groups: np.ndarray) -> list[np.ndarray]:
    return [np.flatnonzero(groups == value) for value in pd.unique(groups)]


def _fit_for_covariance(
    theta: np.ndarray,
    y: np.ndarray,
    design: np.ndarray,
    times: np.ndarray,
    grouped_indices: list[np.ndarray],
    time_unit_seconds: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    tau_squared = float(np.exp(theta[0]))
    residual_variance = float(np.exp(theta[1]))
    rho = float(expit(theta[2]) * 0.999)
    information = np.zeros((design.shape[1], design.shape[1]), dtype=float)
    score = np.zeros(design.shape[1], dtype=float)
    blocks: list[tuple[np.ndarray, tuple[np.ndarray, bool], float]] = []
    for indices in grouped_indices:
        local_time = times[indices]
        distances = np.abs(local_time[:, None] - local_time[None, :]) / float(time_unit_seconds)
        correlation = np.power(rho, distances)
        covariance = residual_variance * correlation + tau_squared
        covariance.flat[:: len(indices) + 1] += 1e-8
        factor = cho_factor(covariance, lower=True, check_finite=False)
        local_design = design[indices]
        inverse_design = cho_solve(factor, local_design, check_finite=False)
        inverse_y = cho_solve(factor, y[indices], check_finite=False)
        information += local_design.T @ inverse_design
        score += local_design.T @ inverse_y
        log_determinant = 2.0 * float(np.log(np.diag(factor[0])).sum())
        blocks.append((indices, factor, log_determinant))
    covariance_beta = np.linalg.pinv(information, rcond=1e-10)
    beta = covariance_beta @ score
    quadratic = 0.0
    log_determinant = 0.0
    for indices, factor, block_log_determinant in blocks:
        residual = y[indices] - design[indices] @ beta
        quadratic += float(residual @ cho_solve(factor, residual, check_finite=False))
        log_determinant += block_log_determinant
    log_likelihood = -0.5 * (
        len(y) * math.log(2.0 * math.pi) + log_determinant + quadratic
    )
    return log_likelihood, beta, covariance_beta


def fit_car1_random_intercept(
    y: np.ndarray,
    design: np.ndarray,
    *,
    groups: np.ndarray,
    times: np.ndarray,
    fixed_effect_names: tuple[str, ...],
    time_unit_seconds: float = 300.0,
    initial_theta: np.ndarray | None = None,
) -> Car1Fit:
    y = np.asarray(y, dtype=float)
    design = np.asarray(design, dtype=float)
    groups = np.asarray(groups)
    times = np.asarray(times, dtype=float)
    grouped_indices = _group_indices(groups)
    if initial_theta is None:
        frame = pd.DataFrame({"group": groups, "y": y})
        case_means = frame.groupby("group", observed=True)["y"].mean()
        between = max(float(case_means.var(ddof=1)), 1e-3)
        centered = y - frame.groupby("group", observed=True)["y"].transform("mean").to_numpy()
        within = max(float(np.var(centered, ddof=1)), 1e-3)
        initial_theta = np.asarray([math.log(between), math.log(within), logit(0.3 / 0.999)])

    def objective(theta: np.ndarray) -> float:
        try:
            value = _fit_for_covariance(
                theta, y, design, times, grouped_indices, time_unit_seconds
            )[0]
            return -value if np.isfinite(value) else 1e100
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            return 1e100

    starts = [np.asarray(initial_theta, dtype=float)]
    for rho in (0.1, 0.7):
        candidate = np.asarray(initial_theta, dtype=float).copy()
        candidate[2] = logit(rho / 0.999)
        starts.append(candidate)
    optimization = None
    attempts: list[str] = []
    for start_number, start in enumerate(starts, start=1):
        candidate = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=[(-12.0, 5.0), (-12.0, 5.0), (-7.0, 7.0)],
            options={"maxiter": 300, "ftol": 1e-10},
        )
        attempts.append(f"start{start_number}:{candidate.message}")
        if optimization is None or candidate.fun < optimization.fun:
            optimization = candidate
        if candidate.success:
            optimization = candidate
            break
    assert optimization is not None
    log_likelihood, beta, covariance_beta = _fit_for_covariance(
        optimization.x, y, design, times, grouped_indices, time_unit_seconds
    )
    return Car1Fit(
        beta=beta,
        covariance_beta=covariance_beta,
        log_likelihood=float(log_likelihood),
        case_random_intercept_variance=float(np.exp(optimization.x[0])),
        residual_variance=float(np.exp(optimization.x[1])),
        correlation_at_five_minutes=float(expit(optimization.x[2]) * 0.999),
        converged=bool(optimization.success),
        optimizer_message="; ".join(attempts),
        fixed_effect_names=fixed_effect_names,
    )


def map_within_contrast_matrix(
    bundle: dict[str, Any],
    grid_mmhg: np.ndarray,
    *,
    reference_mmhg: float,
) -> np.ndarray:
    """Build same-surgery MAP contrasts with between-surgery MAP held fixed."""

    grid = np.asarray(grid_mmhg, dtype=float)
    knots = np.asarray(bundle["knots"]["map"], dtype=float)
    grid_basis = restricted_cubic_spline_basis(grid, knots)
    reference_basis = restricted_cubic_spline_basis(
        np.asarray([float(reference_mmhg)]), knots
    )[0]
    contrast = np.zeros((len(grid), len(bundle["nonlinear_names"])), dtype=float)
    for column in range(grid_basis.shape[1]):
        name = f"map_basis{column + 1}_within"
        position = bundle["nonlinear_names"].index(name)
        standard_deviation = float(bundle["scaling"][name]["standard_deviation"])
        contrast[:, position] = (
            grid_basis[:, column] - reference_basis[column]
        ) / standard_deviation
    return contrast


def _subject_influence(
    bundle: dict[str, Any], fit: Car1Fit, output: Path
) -> tuple[dict[str, Any], pd.DataFrame]:
    data = bundle["data"]
    design = np.asarray(bundle["nonlinear"], dtype=float)
    y = np.asarray(bundle["y"], dtype=float)
    information_by_subject: dict[int, np.ndarray] = {}
    score_by_subject: dict[int, np.ndarray] = {}
    windows_by_subject: dict[int, int] = {}
    surgeries_by_subject: dict[int, int] = {}
    for caseid, local in data.groupby("caseid", observed=True, sort=False):
        indices = local.index.to_numpy(dtype=int)
        times = local["time_seconds"].to_numpy(dtype=float)
        distances = np.abs(times[:, None] - times[None, :]) / 300.0
        covariance = fit.residual_variance * np.power(
            fit.correlation_at_five_minutes, distances
        ) + fit.case_random_intercept_variance
        covariance.flat[:: len(indices) + 1] += 1e-8
        factor = cho_factor(covariance, lower=True, check_finite=False)
        local_design = design[indices]
        local_information = local_design.T @ cho_solve(
            factor, local_design, check_finite=False
        )
        local_score = local_design.T @ cho_solve(
            factor, y[indices], check_finite=False
        )
        subject = int(local["subjectid"].iloc[0])
        information_by_subject[subject] = information_by_subject.get(
            subject, np.zeros_like(local_information)
        ) + local_information
        score_by_subject[subject] = score_by_subject.get(
            subject, np.zeros_like(local_score)
        ) + local_score
        windows_by_subject[subject] = windows_by_subject.get(subject, 0) + len(local)
        surgeries_by_subject[subject] = surgeries_by_subject.get(subject, 0) + 1

    total_information = sum(information_by_subject.values())
    total_score = sum(score_by_subject.values())
    inverse_information = np.linalg.pinv(total_information, rcond=1e-10)
    full_beta = inverse_information @ total_score
    grid = np.arange(50.0, 91.0, 1.0)
    contrast = map_within_contrast_matrix(bundle, grid, reference_mmhg=65.0)
    full_curve = contrast @ full_beta
    rows: list[dict[str, Any]] = []
    for subject in sorted(information_by_subject):
        reduced_information = total_information - information_by_subject[subject]
        reduced_score = total_score - score_by_subject[subject]
        reduced_beta = np.linalg.pinv(reduced_information, rcond=1e-10) @ reduced_score
        reduced_curve = contrast @ reduced_beta
        rows.append(
            {
                "subjectid": subject,
                "surgeries": surgeries_by_subject[subject],
                "windows": windows_by_subject[subject],
                "leverage_trace": float(
                    np.trace(inverse_information @ information_by_subject[subject])
                ),
                "maximum_absolute_curve_change_logit_ra": float(
                    np.max(np.abs(reduced_curve - full_curve))
                ),
                "curve_change_at_50_mmhg": float(reduced_curve[0] - full_curve[0]),
                "curve_change_at_90_mmhg": float(reduced_curve[-1] - full_curve[-1]),
            }
        )
    influence = pd.DataFrame(rows).sort_values(
        "maximum_absolute_curve_change_logit_ra", ascending=False
    )
    influence.to_csv(output / "subject_deletion_influence.csv", index=False)
    top_count = max(1, int(math.ceil(0.01 * len(influence))))
    return (
        {
            "subjects": int(len(influence)),
            "top_one_percent_count": top_count,
            "maximum_curve_change": float(
                influence["maximum_absolute_curve_change_logit_ra"].max()
            ),
            "median_curve_change": float(
                influence["maximum_absolute_curve_change_logit_ra"].median()
            ),
            "maximum_leverage_trace": float(influence["leverage_trace"].max()),
            "top_subjects": influence.head(10).to_dict(orient="records"),
        },
        influence,
    )
