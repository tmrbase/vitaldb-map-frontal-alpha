"""Descriptive lead-to-baseline negative-control analysis for event epochs.

The script estimates the relative-alpha (RA) change from the lead interval
(-10 to -5 min) to the event-specific baseline (-5 to -2 min). Confidence
intervals resample patients, preserving all events within a patient.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def cluster_bootstrap_mean(
    frame: pd.DataFrame,
    value: str,
    cluster: str,
    repetitions: int,
    seed: int,
) -> tuple[float, float, float, int, int]:
    data = frame[[cluster, value]].dropna().copy()
    clusters = data[cluster].drop_duplicates().to_numpy()
    by_cluster = {
        key: group[value].to_numpy(dtype=float) for key, group in data.groupby(cluster, sort=False)
    }
    estimate = float(data[value].mean())
    rng = np.random.default_rng(seed)
    boot = np.empty(repetitions, dtype=float)
    for index in range(repetitions):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        values = np.concatenate([by_cluster[key] for key in sampled])
        boot[index] = values.mean()
    lower, upper = np.quantile(boot, [0.025, 0.975])
    return estimate, float(lower), float(upper), len(data), len(clusters)


def prepare_changes(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    keep = raw.loc[
        raw["interval"].isin(["lead_m10_m5", "baseline_m5_m2"])
        & raw["ra_interval_eligible"].fillna(False).astype(bool),
        [
            "caseid",
            "subjectid",
            "event_number",
            "event_id",
            "event_type",
            "interval",
            "relative_alpha",
        ],
    ].copy()
    wide = keep.pivot(
        index=["caseid", "subjectid", "event_number", "event_id", "event_type"],
        columns="interval",
        values="relative_alpha",
    ).reset_index()
    wide = wide.dropna(subset=["lead_m10_m5", "baseline_m5_m2"])
    wide["lead_to_baseline_ra_change_pp"] = 100.0 * (wide["baseline_m5_m2"] - wide["lead_m10_m5"])
    return wide


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Estimate the lead-to-baseline negative-control event contrast."
    )
    parser.add_argument(
        "--event-intervals", required=True, type=Path, help="Regenerated event-interval CSV."
    )
    parser.add_argument(
        "--output-csv", required=True, type=Path, help="Destination estimate table."
    )
    parser.add_argument(
        "--output-json", required=True, type=Path, help="Destination audit summary."
    )
    parser.add_argument(
        "--bootstrap-repetitions", type=int, default=5000, help="Patient-bootstrap repetitions."
    )
    parser.add_argument("--seed", type=int, default=20260829, help="Bootstrap random seed.")
    args = parser.parse_args()
    if args.bootstrap_repetitions < 1:
        parser.error("--bootstrap-repetitions must be at least 1")

    changes = prepare_changes(args.event_intervals)
    rows: list[dict[str, object]] = []
    for offset, event_type in enumerate(["real", "pseudo"]):
        subset = changes.loc[changes["event_type"] == event_type]
        estimate, lower, upper, events, patients = cluster_bootstrap_mean(
            subset,
            "lead_to_baseline_ra_change_pp",
            "subjectid",
            args.bootstrap_repetitions,
            args.seed + offset,
        )
        rows.append(
            {
                "contrast": f"{event_type}_lead_to_baseline",
                "estimate_pp": estimate,
                "ci95_lower_pp": lower,
                "ci95_upper_pp": upper,
                "events": events,
                "patients": patients,
            }
        )

    match_keys = ["caseid", "subjectid", "event_number"]
    paired = changes.loc[
        changes["event_type"] == "real", match_keys + ["lead_to_baseline_ra_change_pp"]
    ].merge(
        changes.loc[
            changes["event_type"] == "pseudo",
            match_keys + ["lead_to_baseline_ra_change_pp"],
        ],
        on=match_keys,
        how="inner",
        suffixes=("_real", "_pseudo"),
        validate="one_to_one",
    )
    paired["paired_real_minus_reference_pp"] = (
        paired["lead_to_baseline_ra_change_pp_real"]
        - paired["lead_to_baseline_ra_change_pp_pseudo"]
    )
    estimate, lower, upper, events, patients = cluster_bootstrap_mean(
        paired,
        "paired_real_minus_reference_pp",
        "subjectid",
        args.bootstrap_repetitions,
        args.seed + 2,
    )
    rows.append(
        {
            "contrast": "paired_real_minus_reference_lead_to_baseline",
            "estimate_pp": estimate,
            "ci95_lower_pp": lower,
            "ci95_upper_pp": upper,
            "events": events,
            "patients": patients,
        }
    )

    result = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_csv, index=False)
    args.output_json.write_text(
        json.dumps(
            {
                "analysis": "lead-to-baseline RA negative control",
                "intervals": {"lead": "-10 to -5 min", "baseline": "-5 to -2 min"},
                "unit": "RA percentage points",
                "confidence_interval": "patient-level percentile bootstrap",
                "bootstrap_repetitions": args.bootstrap_repetitions,
                "seed": args.seed,
                "results": result.to_dict(orient="records"),
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
