"""Generate the local patient-ranking input required by Table SIII."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from vitaldb_eeg.patient_influence import (
    _subject_influence,
    build_main_design,
    fit_car1_random_intercept,
)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--windows", type=Path,
        default=root / "outputs/maintenance_ra_full_v1/full_windows.parquet",
    )
    parser.add_argument(
        "--output", type=Path,
        default=root / "outputs/maintenance_ra_model_defense_v1",
    )
    args = parser.parse_args()
    windows = (
        pd.read_parquet(args.windows) if args.windows.suffix == ".parquet"
        else pd.read_csv(args.windows)
    )
    reference = build_main_design(windows, minimum_windows=8)
    print("Fitting the original random-intercept CAR(1) influence diagnostic", flush=True)
    fit = fit_car1_random_intercept(
        reference["y"], reference["nonlinear"],
        groups=reference["data"]["caseid"].to_numpy(),
        times=reference["data"]["time_seconds"].to_numpy(dtype=float),
        fixed_effect_names=reference["nonlinear_names"],
    )
    if not fit.converged:
        raise RuntimeError(f"Influence diagnostic did not converge: {fit.optimizer_message}")
    args.output.mkdir(parents=True, exist_ok=True)
    summary, _ = _subject_influence(reference, fit, args.output)
    provenance = {
        "input": str(args.windows.resolve()),
        "input_sha256": hashlib.sha256(args.windows.read_bytes()).hexdigest(),
        "cohort": {
            "patients": int(reference["data"]["subjectid"].nunique()),
            "surgeries": int(reference["data"]["caseid"].nunique()),
            "windows": len(reference["data"]),
        },
        "covariance_held_fixed_during_deletion": True,
        "reference_mmhg": 65, "grid_mmhg": [50, 90, 1],
        "converged": fit.converged,
        "covariance": {
            "random_intercept_variance": fit.case_random_intercept_variance,
            "residual_variance": fit.residual_variance,
            "correlation_at_five_minutes": fit.correlation_at_five_minutes,
        },
        "summary": summary,
        "data_policy": "Local derived patient-level output; do not redistribute.",
    }
    (args.output / "influence_generation.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(f"Ranked {summary['subjects']} patients; exclude {summary['top_one_percent_count']}.")


if __name__ == "__main__":
    main()
