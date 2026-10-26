"""Build reported Figure 1, Supplementary Figure S4, and aggregate source tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLOURS = {
    "navy": "#16324F",
    "blue": "#2F6F9F",
    "red": "#B33A3A",
    "grey": "#6B7280",
    "light_grey": "#D8DDE3",
    "green": "#3B7A57",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Figure 1, Figure S4, and aggregate reporting source tables."
    )
    parser.add_argument(
        "--primary-curve", required=True, help="Standardised primary MAP–RA curve CSV."
    )
    parser.add_argument(
        "--analysis", required=True, help="Directory containing regenerated analyses."
    )
    parser.add_argument(
        "--pro-analysis", required=True, help="Directory containing threshold/event diagnostics."
    )
    parser.add_argument(
        "--pretrend-intervals", required=True, help="Regenerated pre-event interval CSV."
    )
    parser.add_argument("--event-intervals", required=True, help="Regenerated event interval CSV.")
    parser.add_argument(
        "--event-results-80", required=True, help="80% interval-completeness sensitivity JSON."
    )
    parser.add_argument("--output", required=True, help="Directory for figures and source tables.")
    return parser.parse_args()


def read_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def style_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=8.5, colors="#27313A")
    axis.grid(axis="y", color="#E7EAEE", linewidth=0.7, zorder=0)


def build_figure1(primary_curve: Path, analysis: Path, pro_analysis: Path, output: Path) -> None:
    curve = pd.read_csv(primary_curve)
    contrasts = pd.read_csv(analysis / "overlap_restricted_map_contrasts.csv")
    hinge_contrasts = pd.read_csv(pro_analysis / "hinge65_standardised_contrasts.csv")
    display_curve = curve.loc[curve["map_mmhg"].between(55, 105)].copy()
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.7), gridspec_kw={"height_ratios": [1.15, 1]})
    ax = axes[0]
    style_axis(ax)
    ax.axvspan(55, 61, color=COLOURS["light_grey"], alpha=0.65, zorder=0)
    ax.fill_between(
        display_curve["map_mmhg"],
        display_curve["lower_95ci"],
        display_curve["upper_95ci"],
        color=COLOURS["blue"],
        alpha=0.20,
        linewidth=0,
    )
    ax.plot(
        display_curve["map_mmhg"],
        display_curve["adjusted_ra_difference_pp"],
        color=COLOURS["navy"],
        linewidth=2.0,
    )
    ax.axhline(0, color="#20252A", linewidth=0.8)
    ax.axvline(65, color=COLOURS["grey"], linestyle="--", linewidth=0.9)
    ax.text(
        65.7,
        0.97,
        "65 mm Hg: pragmatic\nreference threshold",
        transform=ax.get_xaxis_transform(),
        va="top",
        fontsize=7.4,
        color="#4B5563",
    )
    hinge_ra = hinge_contrasts.loc[
        hinge_contrasts["cohort"].eq("broad_all_windows")
        & hinge_contrasts["outcome"].eq("logit_ra")
    ].set_index(["target_map_mmhg", "reference_map_mmhg"])
    below = 100 * float(hinge_ra.loc[(60, 65), "contrast_response"])
    above = 100 * float(hinge_ra.loc[(65, 70), "contrast_response"])
    ax.text(
        0.98,
        0.05,
        "Targeted secondary hinge model\n5-mm Hg decrease\n65→60: "
        + f"{below:.2f}"
        + " pp\n70→65: "
        + f"{above:.2f}"
        + " pp",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#27313A",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#C9CED4",
            "alpha": 0.94,
        },
    )
    ax.text(
        0.012,
        0.97,
        "Below observed\n5th percentile",
        transform=ax.transAxes,
        va="top",
        fontsize=7.5,
        color="#4B5563",
    )
    ax.set_xlim(55, 105)
    ax.set_xlabel("Mean arterial pressure (mm Hg)", fontsize=9.5)
    ax.set_ylabel("Adjusted RA difference\n(percentage points)", fontsize=9.5)
    ax.text(-0.08, 1.03, "A", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = axes[1]
    style_axis(ax)
    targets = [55, 60, 70, 80, 90, 100]
    definitions = [
        ("all surgeries (model-based standardisation)", "All surgeries", COLOURS["blue"], -0.13),
        (
            "observed range spans target and 65 mm Hg",
            "Range-spanning surgeries (same model)",
            COLOURS["green"],
            0.13,
        ),
    ]
    for definition, label, colour, offset in definitions:
        frame = contrasts.loc[
            contrasts["support_definition"].eq(definition)
            & contrasts["target_map_mmhg"].isin(targets)
        ].copy()
        frame = frame.set_index("target_map_mmhg").loc[targets].reset_index()
        x = np.arange(len(targets)) + offset
        ax.errorbar(
            x,
            frame["contrast_percentage_points"],
            yerr=np.vstack(
                [
                    frame["contrast_percentage_points"] - frame["lower_95ci_percentage_points"],
                    frame["upper_95ci_percentage_points"] - frame["contrast_percentage_points"],
                ]
            ),
            fmt="o" if label == "All surgeries" else "s",
            color=colour,
            ecolor=colour,
            capsize=2.5,
            markersize=4.5,
            linewidth=1.1,
            label=label,
            zorder=3,
        )
        if "Range-spanning" in label:
            for xi, upper_ci, n in zip(
                x,
                frame["upper_95ci_percentage_points"],
                frame["surgeries"],
                strict=True,
            ):
                vertical_offset = 11 if -0.3 < upper_ci < 0 else 5
                ax.annotate(
                    f"n={int(n)}",
                    (xi, upper_ci),
                    xytext=(0, vertical_offset),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6.7,
                    color=colour,
                )
    ax.axhline(0, color="#20252A", linewidth=0.8)
    ax.set_xticks(np.arange(len(targets)), [str(item) for item in targets])
    ax.set_xlabel("Target MAP compared with 65 mm Hg", fontsize=9.5)
    ax.set_ylabel("Adjusted RA contrast\n(percentage points)", fontsize=9.5)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.text(-0.08, 1.03, "B", transform=ax.transAxes, fontsize=12, fontweight="bold")
    fig.tight_layout(h_pad=2.2)
    fig.savefig(output / "Figure_1_MAP_RA_support_BJA.pdf", bbox_inches="tight")
    fig.savefig(output / "Figure_1_MAP_RA_support_BJA.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    contrasts.to_csv(output / "Figure_1_source_data.csv", index=False)


def clustered_mean_ci(
    frame: pd.DataFrame,
    *,
    value_column: str,
    multiplier: float = 1.0,
    repetitions: int = 2000,
    seed: int = 20260829,
) -> tuple[float, float, float]:
    frame = frame.loc[np.isfinite(pd.to_numeric(frame[value_column], errors="coerce"))].copy()
    frame[value_column] = multiplier * pd.to_numeric(frame[value_column], errors="coerce")
    estimate = float(frame[value_column].mean())
    grouped = {key: values[value_column].to_numpy() for key, values in frame.groupby("subjectid")}
    subjects = np.asarray(list(grouped))
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=float)
    for index in range(repetitions):
        sampled = rng.choice(subjects, size=len(subjects), replace=True)
        values[index] = np.concatenate([grouped[item] for item in sampled]).mean()
    lower, upper = np.quantile(values, [0.025, 0.975])
    return estimate, float(lower), float(upper)


def baseline_centered_summaries(
    combined: pd.DataFrame,
    *,
    value_column: str,
    multiplier: float,
    interval_order: list[str],
    time_map: dict[str, float],
    seed: int,
) -> pd.DataFrame:
    rows: list[dict] = []
    for event_type in ["real", "pseudo"]:
        local = combined.loc[combined["event_type"].eq(event_type)].copy()
        baseline = local.loc[
            local["interval"].eq("baseline_m5_m2"),
            ["event_id", value_column],
        ].rename(columns={value_column: "baseline_value"})
        for index, interval in enumerate(interval_order):
            frame = local.loc[
                local["interval"].eq(interval),
                ["event_id", "subjectid", "caseid", value_column],
            ].merge(baseline, on="event_id", how="inner", validate="one_to_one")
            frame["change_from_baseline"] = multiplier * (
                pd.to_numeric(frame[value_column], errors="coerce")
                - pd.to_numeric(frame["baseline_value"], errors="coerce")
            )
            estimate, lower, upper = clustered_mean_ci(
                frame,
                value_column="change_from_baseline",
                seed=seed + index + (0 if event_type == "real" else 100),
            )
            rows.append(
                {
                    "analysis": "Real sustained low-MAP event"
                    if event_type == "real"
                    else "Time-nearest reference epoch",
                    "interval": interval,
                    "time_min": time_map[interval],
                    "estimate": estimate,
                    "lower": lower,
                    "upper": upper,
                    "events": frame["event_id"].nunique(),
                    "surgeries": frame["caseid"].nunique(),
                }
            )
    return pd.DataFrame(rows)


def build_figure2(
    analysis: Path,
    pretrend_intervals: Path,
    event_intervals: Path,
    output: Path,
) -> None:
    pre = pd.read_csv(pretrend_intervals)
    event = pd.read_csv(event_intervals)
    event = event.loc[
        event["interval"].isin(
            ["baseline_m5_m2", "response_0_1", "response_1_3", "response_3_5", "response_5_10"]
        )
    ].copy()
    combined = pd.concat([event, pre], ignore_index=True)
    combined["available"] = combined["primary_ra_subwindows"] >= np.ceil(
        combined["planned_subwindows"] * 0.5
    )
    combined = combined.loc[combined["available"]].copy()
    interval_order = [
        "baseline_m5_m2",
        "pretrend_m2_m1",
        "pretrend_m1_0",
        "response_0_1",
        "response_1_3",
        "response_3_5",
        "response_5_10",
    ]
    time_map = {
        "baseline_m5_m2": -3.5,
        "pretrend_m2_m1": -1.5,
        "pretrend_m1_0": -0.5,
        "response_0_1": 0.5,
        "response_1_3": 2.0,
        "response_3_5": 4.0,
        "response_5_10": 7.5,
    }
    map_rows = []
    for event_type in ["real", "pseudo"]:
        for index, interval in enumerate(interval_order):
            frame = combined.loc[
                combined["event_type"].eq(event_type) & combined["interval"].eq(interval)
            ]
            estimate, lower, upper = clustered_mean_ci(
                frame,
                value_column="art_map",
                seed=20260829 + index + (0 if event_type == "real" else 100),
            )
            map_rows.append(
                {
                    "analysis": "Real sustained low-MAP event"
                    if event_type == "real"
                    else "Time-nearest reference epoch",
                    "interval": interval,
                    "time_min": time_map[interval],
                    "estimate": estimate,
                    "lower": lower,
                    "upper": upper,
                    "events": frame["event_id"].nunique(),
                    "surgeries": frame["caseid"].nunique(),
                }
            )
    map_data = pd.DataFrame(map_rows)
    map_change_data = baseline_centered_summaries(
        combined,
        value_column="art_map",
        multiplier=1.0,
        interval_order=interval_order,
        time_map=time_map,
        seed=20261029,
    )
    ra_change_data = baseline_centered_summaries(
        combined,
        value_column="relative_alpha",
        multiplier=100.0,
        interval_order=interval_order,
        time_map=time_map,
        seed=20261129,
    )

    fig, axes = plt.subplots(
        3, 1, figsize=(7.2, 9.4), gridspec_kw={"height_ratios": [1.0, 1.0, 1.0]}
    )
    ax = axes[0]
    style_axis(ax)
    for label, colour, marker in [
        ("Real sustained low-MAP event", COLOURS["red"], "o"),
        ("Time-nearest reference epoch", COLOURS["grey"], "s"),
    ]:
        frame = map_data.loc[map_data["analysis"].eq(label)]
        ax.plot(
            frame["time_min"],
            frame["estimate"],
            marker=marker,
            color=colour,
            linewidth=1.6,
            markersize=4.3,
            label=label,
        )
        ax.fill_between(
            frame["time_min"], frame["lower"], frame["upper"], color=colour, alpha=0.14, linewidth=0
        )
    ax.axhline(65, color="#20252A", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#20252A", linewidth=0.8)
    ax.axvspan(-2, 0, color=COLOURS["light_grey"], alpha=0.22, zorder=0)
    ax.set_ylabel("Observed MAP\n(mm Hg)", fontsize=9.5)
    ax.text(-0.08, 1.03, "A", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = axes[1]
    style_axis(ax)
    for label, colour, marker in [
        ("Real sustained low-MAP event", COLOURS["red"], "o"),
        ("Time-nearest reference epoch", COLOURS["grey"], "s"),
    ]:
        frame = map_change_data.loc[map_change_data["analysis"].eq(label)]
        ax.plot(
            frame["time_min"],
            frame["estimate"],
            marker=marker,
            color=colour,
            linewidth=1.6,
            markersize=4.3,
            label=label,
        )
        ax.fill_between(
            frame["time_min"], frame["lower"], frame["upper"], color=colour, alpha=0.14, linewidth=0
        )
    ax.axhline(0, color="#20252A", linewidth=0.8)
    ax.axvline(0, color="#20252A", linewidth=0.8)
    ax.axvspan(-2, 0, color=COLOURS["light_grey"], alpha=0.22, zorder=0)
    ax.set_ylabel("MAP change from each epoch's\nown baseline (mm Hg)", fontsize=9.2)
    real_last_map = map_change_data.loc[
        map_change_data["analysis"].eq("Real sustained low-MAP event")
        & map_change_data["interval"].eq("pretrend_m1_0")
    ].iloc[0]
    ax.annotate(
        f"{real_last_map['estimate']:.1f} mm Hg",
        (real_last_map["time_min"], real_last_map["estimate"]),
        xytext=(-7, -14),
        textcoords="offset points",
        ha="right",
        fontsize=7.2,
        color=COLOURS["red"],
    )
    ax.text(-0.08, 1.03, "B", transform=ax.transAxes, fontsize=12, fontweight="bold")

    ax = axes[2]
    style_axis(ax)
    for label, colour, marker in [
        ("Real sustained low-MAP event", COLOURS["red"], "o"),
        ("Time-nearest reference epoch", COLOURS["grey"], "s"),
    ]:
        frame = ra_change_data.loc[ra_change_data["analysis"].eq(label)]
        ax.plot(
            frame["time_min"],
            frame["estimate"],
            marker=marker,
            color=colour,
            linewidth=1.6,
            markersize=4.3,
            label=label,
        )
        ax.fill_between(
            frame["time_min"], frame["lower"], frame["upper"], color=colour, alpha=0.14, linewidth=0
        )
    ax.axhline(0, color="#20252A", linewidth=0.8)
    ax.axvline(0, color="#20252A", linewidth=0.8)
    ax.axvspan(-2, 0, color=COLOURS["light_grey"], alpha=0.22, zorder=0)
    ax.set_ylabel("RA change from each epoch's own\nbaseline (percentage points)", fontsize=9.2)
    real_last_ra = ra_change_data.loc[
        ra_change_data["analysis"].eq("Real sustained low-MAP event")
        & ra_change_data["interval"].eq("pretrend_m1_0")
    ].iloc[0]
    ax.annotate(
        f"{real_last_ra['estimate']:.2f} pp",
        (real_last_ra["time_min"], real_last_ra["estimate"]),
        xytext=(-7, -15),
        textcoords="offset points",
        ha="right",
        fontsize=7.2,
        color=COLOURS["red"],
    )
    ax.text(-0.08, 1.03, "C", transform=ax.transAxes, fontsize=12, fontweight="bold")
    axes[2].set_xlabel("Interval relative to the event start point (t=0)", fontsize=9.5)
    for ax in axes:
        ax.set_xticks(
            [-3.5, -1.5, -0.5, 0.5, 2, 4, 7.5],
            ["−5 to −2", "−2 to −1", "−1 to 0", "0 to 1", "1 to 3", "3 to 5", "5 to 10"],
        )
        ax.tick_params(axis="x", labelrotation=0, labelsize=7.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=8,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.995),
    )
    fig.tight_layout(h_pad=2.0, rect=[0, 0, 1, 0.965])
    fig.savefig(output / "FigureS4_event_trajectory_diagnostic_BJA.pdf", bbox_inches="tight")
    fig.savefig(
        output / "FigureS4_event_trajectory_diagnostic_BJA.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    map_data.assign(panel="A").to_csv(output / "FigureS4A_source_data.csv", index=False)
    map_change_data.assign(panel="B").to_csv(output / "FigureS4B_source_data.csv", index=False)
    ra_change_data.assign(panel="C").to_csv(output / "FigureS4C_source_data.csv", index=False)


def result_table(
    result: dict, analysis_label: str, cohort_label: str, interval_names: list[str]
) -> list[dict]:
    block = result["raw_ra_percentage_points"][analysis_label]
    estimates = block.get("interval_estimates", block.get("estimates"))
    global_test = block.get("global_test_all_four_zero", {})
    rows = []
    for interval in interval_names:
        item = estimates[interval]
        rows.append(
            {
                "cohort": cohort_label,
                "analysis": analysis_label,
                "interval": interval,
                "estimate_pp": item["estimate"],
                "lower_95ci_pp": item["normal_95ci_lower"],
                "upper_95ci_pp": item["normal_95ci_upper"],
                "rows": block["rows"],
                "events_or_pairs": block.get("events", block.get("pairs")),
                "surgeries": block["surgeries"],
                "subjects": block["subjects"],
                "global_p": global_test.get("p_value"),
            }
        )
    return rows


def format_estimate_ci(
    estimate: float, lower: float, upper: float, *, digits: int = 3, multiplier: float = 1.0
) -> str:
    return f"{multiplier * estimate:.{digits}f} ({multiplier * lower:.{digits}f} to {multiplier * upper:.{digits}f})"


def format_log_contrast_as_percent(
    estimate: float, lower: float, upper: float, *, digits: int = 2
) -> str:
    values = 100.0 * (np.exp(np.asarray([estimate, lower, upper], dtype=float)) - 1.0)
    return f"{values[0]:.{digits}f} ({values[1]:.{digits}f} to {values[2]:.{digits}f})"


def build_tables(analysis: Path, pro_analysis: Path, event_results_80: Path, output: Path) -> None:
    contrasts = pd.read_csv(analysis / "overlap_restricted_map_contrasts.csv")
    contrasts.to_csv(output / "common_support_map_contrasts.csv", index=False)

    pretrend = read_json(analysis / "event_pretrend_results.json")
    rows = []
    intervals = ["pretrend_m2_m1", "pretrend_m1_0"]
    for cohort in ["broad_no_positive_recorded_volatile", "strict_measured_zero_volatile"]:
        result = pretrend[cohort]
        for analysis_name in ["real_event_adjusted", "paired_real_vs_pseudo"]:
            rows.extend(result_table(result, analysis_name, cohort, intervals))
            contrast = result["raw_ra_percentage_points"][analysis_name][
                "immediate_pretrend_change_last_minus_first"
            ]
            base = result["raw_ra_percentage_points"][analysis_name]
            rows.append(
                {
                    "cohort": cohort,
                    "analysis": analysis_name,
                    "interval": "pretrend_m1_0 minus pretrend_m2_m1",
                    "estimate_pp": contrast["estimate"],
                    "lower_95ci_pp": contrast["normal_95ci_lower"],
                    "upper_95ci_pp": contrast["normal_95ci_upper"],
                    "rows": base["rows"],
                    "events_or_pairs": base.get("events", base.get("pairs")),
                    "surgeries": base["surgeries"],
                    "subjects": base["subjects"],
                    "global_p": contrast["p_value"],
                }
            )
    pd.DataFrame(rows).to_csv(output / "immediate_pre_event_ra.csv", index=False)

    strict = read_json(analysis / "strict_volatile_event_sensitivity.json")
    post = ["response_0_1", "response_1_3", "response_3_5", "response_5_10"]
    rows = []
    for analysis_name in ["real_event_adjusted", "paired_real_vs_pseudo"]:
        rows.extend(result_table(strict, analysis_name, "strict_measured_zero_volatile", post))
    pd.DataFrame(rows).to_csv(output / "strict_volatile_event_sensitivity.csv", index=False)

    availability = read_json(event_results_80)
    rows = []
    for analysis_name in ["real_event_adjusted", "paired_real_vs_pseudo"]:
        rows.extend(
            result_table(availability, analysis_name, "80_percent_interval_completeness", post)
        )
    pd.DataFrame(rows).to_csv(output / "event_80_percent_completeness.csv", index=False)

    summaries = pd.read_csv(pro_analysis / "hinge65_model_summary.csv")
    hinge_contrasts = pd.read_csv(pro_analysis / "hinge65_standardised_contrasts.csv")
    slopes = pd.read_csv(pro_analysis / "hinge65_log_scale_slopes.csv")
    robust = pd.read_csv(pro_analysis / "hinge65_patient_cluster_CR1S.csv")
    diagnostics = pd.read_csv(pro_analysis / "hinge65_diagnostics.csv")
    support = pd.read_csv(analysis / "case_level_map_support.csv")
    specifications = [
        ("Full cohort", "broad_all_windows", "logit_ra", 1, "broad_no_positive_recorded_volatile"),
        ("Surgeries crossing 65 mm Hg", "surgeries_crossing_65", "logit_ra", 3, None),
        (
            "Strict measured-zero volatile cohort",
            "strict_measured_zero_volatile",
            "logit_ra",
            2,
            "strict_measured_zero_volatile",
        ),
        (
            "Absolute alpha power sensitivity",
            "broad_all_windows",
            "log_alpha_power",
            4,
            "broad_no_positive_recorded_volatile",
        ),
    ]
    threshold_rows = []
    threshold_response_rows = []
    threshold_slope_rows = []
    alpha_conversion_rows = []
    for label, cohort, outcome, diagnostic_index, support_cohort in specifications:
        summary = summaries.loc[
            summaries["cohort"].eq(cohort) & summaries["outcome"].eq(outcome)
        ].iloc[0]
        diag = diagnostics.loc[diagnostics["analysis_index"].eq(diagnostic_index)].iloc[0]
        cr1s = robust.loc[robust["analysis"].eq(f"{cohort}_{outcome}")].iloc[0]
        difference = slopes.loc[
            slopes["cohort"].eq(cohort)
            & slopes["outcome"].eq(outcome)
            & slopes["estimand"].str.startswith("difference in 5-mm Hg")
        ].iloc[0]
        above_slope = slopes.loc[
            slopes["cohort"].eq(cohort)
            & slopes["outcome"].eq(outcome)
            & slopes["estimand"].eq("5-mm Hg MAP decrease at MAP >=65")
        ].iloc[0]
        below_slope = slopes.loc[
            slopes["cohort"].eq(cohort)
            & slopes["outcome"].eq(outcome)
            & slopes["estimand"].eq("5-mm Hg MAP decrease at MAP <65")
        ].iloc[0]
        local_contrasts = hinge_contrasts.loc[
            hinge_contrasts["cohort"].eq(cohort) & hinge_contrasts["outcome"].eq(outcome)
        ].set_index(["target_map_mmhg", "reference_map_mmhg"])
        scale = 100.0 if outcome == "logit_ra" else 1.0
        response_scale_label = (
            "RA percentage points" if outcome == "logit_ra" else "absolute alpha power"
        )
        slope_outcome_label = "logit(RA)" if outcome == "logit_ra" else "absolute alpha power"
        support_60 = ""
        support_70 = ""
        if support_cohort is not None:
            support_frame = support.loc[support["cohort"].eq(support_cohort)].set_index(
                "target_map_mmhg"
            )
            support_60 = int(support_frame.loc[60, "surgeries_range_spans_target_and_reference"])
            support_70 = int(support_frame.loc[70, "surgeries_range_spans_target_and_reference"])
        if outcome == "logit_ra":
            contrast_60_65 = format_estimate_ci(
                *local_contrasts.loc[
                    (60, 65), ["contrast_response", "lower_95ci_response", "upper_95ci_response"]
                ],
                digits=2,
                multiplier=scale,
            )
            contrast_65_70 = format_estimate_ci(
                *local_contrasts.loc[
                    (65, 70), ["contrast_response", "lower_95ci_response", "upper_95ci_response"]
                ],
                digits=2,
                multiplier=scale,
            )
            response_estimand = "standardised RA difference, percentage points"
            slope_scale = "logit(RA) change"
            slope_above = format_estimate_ci(
                above_slope["estimate_log_scale"],
                above_slope["lower_95ci_log_scale"],
                above_slope["upper_95ci_log_scale"],
                digits=3,
            )
            slope_below = format_estimate_ci(
                below_slope["estimate_log_scale"],
                below_slope["lower_95ci_log_scale"],
                below_slope["upper_95ci_log_scale"],
                digits=3,
            )
            slope_difference = format_estimate_ci(
                difference["estimate_log_scale"],
                difference["lower_95ci_log_scale"],
                difference["upper_95ci_log_scale"],
                digits=3,
            )
        else:
            contrast_60_65 = format_log_contrast_as_percent(
                below_slope["estimate_log_scale"],
                below_slope["lower_95ci_log_scale"],
                below_slope["upper_95ci_log_scale"],
            )
            contrast_65_70 = format_log_contrast_as_percent(
                above_slope["estimate_log_scale"],
                above_slope["lower_95ci_log_scale"],
                above_slope["upper_95ci_log_scale"],
            )
            response_estimand = "geometric-mean ratio, expressed as percentage change"
            slope_scale = "percentage change in geometric mean"
            slope_above = contrast_65_70
            slope_below = contrast_60_65
            slope_difference = format_log_contrast_as_percent(
                difference["estimate_log_scale"],
                difference["lower_95ci_log_scale"],
                difference["upper_95ci_log_scale"],
            )
            for estimand_name, slope_row in [
                ("MAP 70 to 65 mm Hg", above_slope),
                ("MAP 65 to 60 mm Hg", below_slope),
                ("below-minus-above slope ratio", difference),
            ]:
                alpha_conversion_rows.append(
                    {
                        "estimand": estimand_name,
                        "log_contrast": float(slope_row["estimate_log_scale"]),
                        "lower_95ci_log": float(slope_row["lower_95ci_log_scale"]),
                        "upper_95ci_log": float(slope_row["upper_95ci_log_scale"]),
                        "transformation": "100*(exp(log contrast)-1)",
                        "percent_change_95ci": format_log_contrast_as_percent(
                            slope_row["estimate_log_scale"],
                            slope_row["lower_95ci_log_scale"],
                            slope_row["upper_95ci_log_scale"],
                        ),
                    }
                )

        response_row = {
            "analysis": label,
            "outcome": response_scale_label,
            "estimand_scale": response_estimand,
            "surgeries": int(summary["surgeries"]),
            "windows": int(summary["windows"]),
            "windows_below_65": int(diag["windows_below65"]),
            "surgeries_crossing_65": int(diag["surgeries_crossing65"]),
            "surgeries_spanning_60_65": support_60,
            "surgeries_spanning_65_70": support_70,
            "MAP_65_to_60_estimate_95ci": contrast_60_65,
            "MAP_70_to_65_estimate_95ci": contrast_65_70,
        }
        slope_row = {
            "analysis": label,
            "outcome": slope_outcome_label,
            "slope_scale": slope_scale,
            "5mm_decrease_above_65_estimate_95ci": slope_above,
            "5mm_decrease_below_65_estimate_95ci": slope_below,
            "below_minus_above_estimate_95ci": slope_difference,
            "model_based_p_for_slope_change": float(summary["hinge_lrt_p_value"]),
            "patient_cluster_cr1s_p_for_slope_change": float(cr1s["patient_cluster_CR1S_p_value"]),
        }
        threshold_response_rows.append(response_row)
        threshold_slope_rows.append(slope_row)
        threshold_rows.append({**response_row, **slope_row})
    pd.DataFrame(threshold_rows).to_csv(output / "threshold_hinge_analysis.csv", index=False)
    pd.DataFrame(threshold_response_rows).to_csv(
        output / "threshold_support_and_contrasts.csv", index=False
    )
    pd.DataFrame(threshold_slope_rows).to_csv(
        output / "threshold_slopes_and_inference.csv",
        index=False,
        float_format="%.3f",
    )
    pd.DataFrame(alpha_conversion_rows).to_csv(
        output / "absolute_alpha_hinge_conversion_audit.csv", index=False
    )

    global_tests = pd.read_csv(pro_analysis / "pretrend_global_tests.csv")
    interval_estimates = pd.read_csv(pro_analysis / "pretrend_interval_estimates.csv")
    event_map = pd.read_csv(pro_analysis / "event_centred_map_model.csv").iloc[0]
    alignment = read_json(pro_analysis / "event_time_alignment_audit.json")
    event_rows = []
    for analysis_name, label in [
        ("real_event_adjusted", "Real-event intervals before t=0: joint 2-df test"),
        ("real_minus_pseudo", "Real-minus-reference intervals before t=0: joint 2-df test"),
    ]:
        item = global_tests.loc[
            global_tests["cohort"].eq("broad") & global_tests["analysis"].eq(analysis_name)
        ].iloc[0]
        event_rows.append(
            {
                "analysis": label,
                "scale_or_estimand": "Two before-t=0 interval coefficients jointly equal zero",
                "estimate_95ci": "Joint test",
                "p_value": float(item["p_value"]),
                "rows": int(item["rows"]),
                "events_or_pairs": int(item["events_or_pairs"]),
                "surgeries": int(item["surgeries"]),
            }
        )
    final_pre = interval_estimates.loc[
        interval_estimates["cohort"].eq("broad")
        & interval_estimates["analysis"].eq("real_minus_pseudo")
        & interval_estimates["interval"].eq("pretrend_m1_0")
    ].iloc[0]
    event_rows.append(
        {
            "analysis": "Final minute before t=0: real-minus-reference difference",
            "scale_or_estimand": "RA change, percentage points",
            "estimate_95ci": format_estimate_ci(
                final_pre["estimate_percentage_points"],
                final_pre["lower_95ci"],
                final_pre["upper_95ci"],
                digits=2,
            ),
            "p_value": "",
            "rows": int(final_pre["rows"]),
            "events_or_pairs": int(final_pre["events_or_pairs"]),
            "surgeries": int(final_pre["surgeries"]),
        }
    )
    event_rows.append(
        {
            "analysis": "Secondary event-centred contemporaneous MAP model",
            "scale_or_estimand": "Change in logit(RA) per 10-mm Hg higher MAP",
            "estimate_95ci": format_estimate_ci(
                event_map["map_beta_per_10_mmhg_higher"],
                event_map["map_lower_95ci"],
                event_map["map_upper_95ci"],
                digits=3,
            ),
            "p_value": float(event_map["map_p_value"]),
            "rows": int(event_map["rows"]),
            "events_or_pairs": int(event_map["events"]),
            "surgeries": int(event_map["surgeries"]),
        }
    )
    event_rows.append(
        {
            "analysis": "Before-t=0 MAP alignment audit",
            "scale_or_estimand": "Events with any 30-s median MAP <65 mm Hg before t=0",
            "estimate_95ci": f"{alignment['events_with_any_pre0_30s_median_map_below65']}/{alignment['real_events_with_saved_pre0_subwindows']} ({alignment['percent_with_any_pre0_30s_median_map_below65']:.1f}%)",
            "p_value": "",
            "rows": "",
            "events_or_pairs": int(alignment["real_events_with_saved_pre0_subwindows"]),
            "surgeries": "",
        }
    )
    pd.DataFrame(event_rows).to_csv(output / "event_centred_diagnostics.csv", index=False)


def main() -> None:
    args = parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    primary_curve = Path(args.primary_curve)
    analysis = Path(args.analysis)
    pro_analysis = Path(args.pro_analysis)
    build_figure1(primary_curve, analysis, pro_analysis, output)
    build_figure2(
        analysis,
        Path(args.pretrend_intervals),
        Path(args.event_intervals),
        output,
    )
    build_tables(analysis, pro_analysis, Path(args.event_results_80), output)


if __name__ == "__main__":
    main()
