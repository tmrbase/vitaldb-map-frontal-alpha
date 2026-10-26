"""Describe surgery-level MAP support for each standardised MAP contrast."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def as_true(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "t", "1"})


def support_table(
    frame: pd.DataFrame,
    *,
    reference: float = 65.0,
    half_width: float = 2.5,
) -> pd.DataFrame:
    grouped = frame.groupby("caseid", sort=True)["art_map"]
    ranges = grouped.agg(min_map="min", max_map="max", windows="size")

    rows: list[dict[str, float | int]] = []
    values_by_case = {int(caseid): values.to_numpy(dtype=float) for caseid, values in grouped}
    total = len(ranges)
    for target in range(50, 111):
        lower = min(float(target), reference)
        upper = max(float(target), reference)
        span = (ranges["min_map"] <= lower) & (ranges["max_map"] >= upper)
        neighbourhood = np.fromiter(
            (
                np.any(np.abs(values - float(target)) <= half_width)
                and np.any(np.abs(values - reference) <= half_width)
                for values in values_by_case.values()
            ),
            dtype=bool,
            count=total,
        )
        rows.append(
            {
                "target_map_mmhg": target,
                "reference_map_mmhg": reference,
                "surgeries_total": total,
                "surgeries_range_spans_target_and_reference": int(span.sum()),
                "range_span_percent": 100.0 * float(span.mean()),
                "surgeries_with_windows_near_target_and_reference": int(neighbourhood.sum()),
                "near_both_percent": 100.0 * float(neighbourhood.mean()),
                "neighbourhood_half_width_mmhg": half_width,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise surgery-level observed MAP support for reported contrasts."
    )
    parser.add_argument(
        "--windows", type=Path, required=True, help="Regenerated maintenance-window CSV."
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Directory for audit tables."
    )
    parser.add_argument(
        "--half-width",
        type=float,
        default=2.5,
        help="Support radius around each MAP target, mm Hg.",
    )
    args = parser.parse_args()
    if args.half_width <= 0:
        parser.error("--half-width must be positive")

    usecols = [
        "caseid",
        "subjectid",
        "art_map",
        "primary_complete_window",
        "strict_measured_zero_volatile",
    ]
    frame = pd.read_csv(args.windows, usecols=usecols)
    frame = frame.loc[
        as_true(frame["primary_complete_window"]) & np.isfinite(frame["art_map"])
    ].copy()

    counts = frame.groupby("caseid").size()
    eligible_cases = counts.index[counts >= 8]
    frame = frame.loc[frame["caseid"].isin(eligible_cases)].copy()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    broad = support_table(frame, half_width=args.half_width)
    strict = support_table(
        frame.loc[as_true(frame["strict_measured_zero_volatile"])],
        half_width=args.half_width,
    )
    broad.insert(0, "cohort", "broad_no_positive_recorded_volatile")
    strict.insert(0, "cohort", "strict_measured_zero_volatile")
    combined = pd.concat([broad, strict], ignore_index=True)
    combined.to_csv(output_dir / "case_level_map_support.csv", index=False)

    selected = combined.loc[combined["target_map_mmhg"].isin([55, 60, 65, 70, 80, 90, 100])].copy()
    selected.to_csv(output_dir / "case_level_map_support_selected.csv", index=False)

    summary = {
        "input": str(args.windows),
        "definition": {
            "range_span": (
                "observed surgery-level MAP range contains both the target and "
                "the 65-mm Hg reference"
            ),
            "near_both": (
                "the surgery contains at least one observed window within the stated "
                "half-width of both the target and 65-mm Hg reference"
            ),
        },
        "broad": {
            "patients": int(frame["subjectid"].nunique()),
            "surgeries": int(frame["caseid"].nunique()),
            "windows": len(frame),
        },
        "strict": {
            "patients": int(
                frame.loc[as_true(frame["strict_measured_zero_volatile"]), "subjectid"].nunique()
            ),
            "surgeries": int(
                frame.loc[as_true(frame["strict_measured_zero_volatile"]), "caseid"].nunique()
            ),
            "windows": int(as_true(frame["strict_measured_zero_volatile"]).sum()),
        },
        "selected_targets": selected.to_dict(orient="records"),
    }
    (output_dir / "case_level_map_support_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(selected.to_string(index=False))


if __name__ == "__main__":
    main()
