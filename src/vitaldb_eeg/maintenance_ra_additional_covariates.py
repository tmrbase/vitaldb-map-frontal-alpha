"""Extract additional haemodynamic and oxygenation covariates by analysis window."""

from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import ensure_output_directories
from .data import import_vitaldb, vital_path
from .maintenance_ra_extraction import _native_records, _track_cadence, _window_native

DEVICES = {
    "EV1000": ("EV1000/CO", "EV1000/CI"),
    "Vigileo": ("Vigileo/CO", "Vigileo/CI"),
    "Vigilance": ("Vigilance/CO", "Vigilance/CI"),
    "CardioQ": ("CardioQ/CO", "CardioQ/CI"),
}
DEVICE_PRIORITY = tuple(DEVICES)
N2O_TRACK = "Primus/FEN2O"
HR_TRACK = "Solar8000/HR"
SPO2_TRACK = "Solar8000/PLETH_SPO2"


def _output_root(cfg: dict[str, Any]) -> Path:
    return ensure_output_directories(cfg).output_root / "maintenance_ra_additional_covariates_v1"


def _window_median_with_coverage(
    times: np.ndarray,
    values: np.ndarray,
    *,
    center: float,
    seconds: float,
    cadence: float,
    valid_range: tuple[float, float],
    coverage_min: float = 0.80,
) -> tuple[float, bool, int, int]:
    _, local_values = _window_native(times, values, center, seconds)
    low, high = valid_range
    valid = np.isfinite(local_values) & (local_values >= low) & (local_values <= high)
    expected = max(round(seconds / cadence), 1) if np.isfinite(cadence) and cadence > 0 else 0
    required = max(math.ceil(expected * coverage_min), 1) if expected else 0
    eligible = bool(expected and int(valid.sum()) >= required)
    median = float(np.median(local_values[valid])) if eligible else float("nan")
    return median, eligible, int(valid.sum()), int(expected)


def _simple_window_median(
    times: np.ndarray,
    values: np.ndarray,
    *,
    center: float,
    seconds: float,
    valid_range: tuple[float, float],
) -> tuple[float, bool]:
    _, local = _window_native(times, values, center, seconds)
    low, high = valid_range
    valid = local[np.isfinite(local) & (local >= low) & (local <= high)]
    return (float(np.median(valid)), True) if valid.size else (float("nan"), False)


def _extract_case(
    cfg: dict[str, Any], caseid: int, centers: list[float]
) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = ensure_output_directories(cfg)
    requested = [N2O_TRACK, HR_TRACK, SPO2_TRACK]
    requested.extend(track for pair in DEVICES.values() for track in pair)
    vitaldb = import_vitaldb()
    vital = vitaldb.VitalFile(str(vital_path(paths.vital_root, caseid)), track_names=requested)
    native = {name: _native_records(vital, name) for name in requested}
    maintenance_start = float(min(centers) - 15.0)
    maintenance_stop = float(max(centers) + 15.0)

    n2o_times, n2o_values = native[N2O_TRACK]
    native_n2o = n2o_values[
        (n2o_times >= maintenance_start)
        & (n2o_times <= maintenance_stop)
        & np.isfinite(n2o_values)
        & (n2o_values >= 0.0)
        & (n2o_values <= 100.0)
    ]
    n2o_any_positive_native = bool((native_n2o >= 1.0).any())

    hr_times, hr_values = native[HR_TRACK]
    spo2_times, spo2_values = native[SPO2_TRACK]
    hr_cadence = _track_cadence(hr_times, start=maintenance_start, stop=maintenance_stop)
    spo2_cadence = _track_cadence(spo2_times, start=maintenance_start, stop=maintenance_stop)

    per_device: dict[str, dict[str, Any]] = {}
    for device, (co_name, ci_name) in DEVICES.items():
        co_times, co_values = native[co_name]
        ci_times, ci_values = native[ci_name]
        per_device[device] = {
            "co_times": co_times,
            "co_values": co_values,
            "ci_times": ci_times,
            "ci_values": ci_values,
            "ci_valid_windows": 0,
        }
        for center in centers:
            _, valid = _simple_window_median(
                ci_times,
                ci_values,
                center=center,
                seconds=30.0,
                valid_range=(0.2, 10.0),
            )
            per_device[device]["ci_valid_windows"] += int(valid)

    selected_device = max(
        DEVICE_PRIORITY,
        key=lambda device: (
            int(per_device[device]["ci_valid_windows"]),
            -DEVICE_PRIORITY.index(device),
        ),
    )
    if int(per_device[selected_device]["ci_valid_windows"]) == 0:
        selected_device = "none"

    rows: list[dict[str, Any]] = []
    for center in centers:
        n2o, n2o_observed = _simple_window_median(
            n2o_times,
            n2o_values,
            center=center,
            seconds=30.0,
            valid_range=(0.0, 100.0),
        )
        hr, hr_qc, hr_valid, hr_expected = _window_median_with_coverage(
            hr_times,
            hr_values,
            center=center,
            seconds=30.0,
            cadence=hr_cadence,
            valid_range=(30.0, 220.0),
        )
        spo2, spo2_qc, spo2_valid, spo2_expected = _window_median_with_coverage(
            spo2_times,
            spo2_values,
            center=center,
            seconds=30.0,
            cadence=spo2_cadence,
            valid_range=(70.0, 100.0),
        )
        if selected_device == "none":
            co = float("nan")
            ci = float("nan")
            co_qc = False
            ci_qc = False
        else:
            selected = per_device[selected_device]
            co, co_qc = _simple_window_median(
                selected["co_times"],
                selected["co_values"],
                center=center,
                seconds=30.0,
                valid_range=(0.5, 20.0),
            )
            ci, ci_qc = _simple_window_median(
                selected["ci_times"],
                selected["ci_values"],
                center=center,
                seconds=30.0,
                valid_range=(0.2, 10.0),
            )
        rows.append(
            {
                "caseid": int(caseid),
                "time_seconds": float(center),
                "n2o_expired": n2o,
                "n2o_observed": bool(n2o_observed),
                "n2o_positive": bool(n2o_observed and n2o >= 1.0),
                # Persist the native-cadence surgery-level flag in every shard so
                # that a resumed extraction produces the same N2O cohort.
                "n2o_any_positive_native": n2o_any_positive_native,
                "hr": hr,
                "hr_qc": bool(hr_qc),
                "hr_valid_samples": hr_valid,
                "hr_expected_samples": hr_expected,
                "spo2": spo2,
                "spo2_qc": bool(spo2_qc),
                "spo2_valid_samples": spo2_valid,
                "spo2_expected_samples": spo2_expected,
                "hemodynamic_device": selected_device,
                "co": co,
                "co_qc": bool(co_qc),
                "ci": ci,
                "ci_qc": bool(ci_qc),
            }
        )
    frame = pd.DataFrame(rows)
    audit = {
        "caseid": int(caseid),
        "scheduled_primary_windows": len(frame),
        "n2o_track_has_records": bool(n2o_times.size),
        "n2o_observed_windows": int(frame["n2o_observed"].sum()),
        "n2o_positive_windows": int(frame["n2o_positive"].sum()),
        "n2o_any_positive_native": n2o_any_positive_native,
        "hr_cadence_seconds": float(hr_cadence),
        "hr_valid_windows": int(frame["hr_qc"].sum()),
        "spo2_cadence_seconds": float(spo2_cadence),
        "spo2_valid_windows": int(frame["spo2_qc"].sum()),
        "hemodynamic_device": selected_device,
        "co_valid_windows": int(frame["co_qc"].sum()),
        "ci_valid_windows": int(frame["ci_qc"].sum()),
    }
    return frame, audit


def _extract_worker(
    cfg: dict[str, Any], caseid: int, centers: list[float], destination: str
) -> dict[str, Any]:
    target = Path(destination)
    if target.exists():
        frame = pd.read_parquet(target)
        return {
            "caseid": int(caseid),
            "status": "cached",
            "scheduled_primary_windows": len(frame),
        }
    temporary = target.with_suffix(f".partial-{os.getpid()}.parquet")
    try:
        frame, audit = _extract_case(cfg, caseid, centers)
        target.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(temporary, index=False)
        temporary.replace(target)
        return {"status": "ok", **audit}
    except Exception as error:  # noqa: BLE001 - every case failure is retained for audit
        if temporary.exists():
            temporary.unlink()
        return {
            "caseid": int(caseid),
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error),
        }


def extract_additional_covariates(
    cfg: dict[str, Any], *, workers: int | None = None, force: bool = False
) -> dict[str, Any]:
    paths = ensure_output_directories(cfg)
    source = paths.output_root / "maintenance_ra_full_v1" / "full_windows.parquet"
    windows = pd.read_parquet(source)
    complete = windows.loc[windows["primary_complete_window"].fillna(False)].copy()
    counts = complete.groupby("caseid", observed=True).size()
    caseids = counts.loc[counts >= 8].index.astype(int)
    complete = complete.loc[complete["caseid"].isin(caseids)]
    output = _output_root(cfg)
    shards = output / "case_shards"
    jobs: list[tuple[int, list[float], Path]] = []
    for caseid, local in complete.groupby("caseid", observed=True):
        target = shards / f"case_{int(caseid):04d}.parquet"
        if force and target.exists():
            target.unlink()
        jobs.append((int(caseid), local["time_seconds"].astype(float).tolist(), target))

    audits: list[dict[str, Any]] = []
    worker_count = max(1, int(workers or cfg["runtime"]["workers"]))
    with ProcessPoolExecutor(max_workers=min(worker_count, len(jobs))) as pool:
        futures = {
            pool.submit(_extract_worker, cfg, caseid, centers, str(target)): caseid
            for caseid, centers, target in jobs
        }
        for future in as_completed(futures):
            audits.append(future.result())

    audit_frame = pd.DataFrame(audits).sort_values("caseid")
    # Cached rows are re-audited from their shards to keep summaries complete.
    frames = [pd.read_parquet(target) for _, _, target in jobs if target.exists()]
    covariates = pd.concat(frames, ignore_index=True).sort_values(["caseid", "time_seconds"])
    case_summary = (
        covariates.groupby("caseid", observed=True)
        .agg(
            scheduled_primary_windows=("time_seconds", "size"),
            n2o_observed_windows=("n2o_observed", "sum"),
            n2o_positive_windows=("n2o_positive", "sum"),
            n2o_any_positive_native=("n2o_any_positive_native", "max"),
            hr_valid_windows=("hr_qc", "sum"),
            spo2_valid_windows=("spo2_qc", "sum"),
            co_valid_windows=("co_qc", "sum"),
            ci_valid_windows=("ci_qc", "sum"),
            hemodynamic_device=("hemodynamic_device", "first"),
        )
        .reset_index()
    )
    case_summary["n2o_any_positive_native"] = case_summary["n2o_any_positive_native"].astype(bool)

    output.mkdir(parents=True, exist_ok=True)
    covariates.to_parquet(output / "additional_covariates.parquet", index=False)
    covariates.to_csv(output / "additional_covariates.csv", index=False)
    case_summary.to_csv(output / "case_covariate_audit.csv", index=False)
    audit_frame.to_csv(output / "extraction_job_audit.csv", index=False)
    summary = {
        "status": "post hoc additional covariate extraction",
        "surgeries": int(covariates["caseid"].nunique()),
        "windows": len(covariates),
        "n2o_positive_surgeries": int(case_summary["n2o_any_positive_native"].sum()),
        "n2o_strict_measured_zero_surgeries": int(
            (
                ~case_summary["n2o_any_positive_native"]
                & case_summary["n2o_observed_windows"].ge(
                    0.8 * case_summary["scheduled_primary_windows"]
                )
            ).sum()
        ),
        "ci_8plus_surgeries": int(case_summary["ci_valid_windows"].ge(8).sum()),
        "hr_spo2_8plus_surgeries": int(
            (
                case_summary["hr_valid_windows"].ge(8) & case_summary["spo2_valid_windows"].ge(8)
            ).sum()
        ),
        "device_counts_ci_8plus": {
            str(device): int(count)
            for device, count in case_summary.loc[
                case_summary["ci_valid_windows"].ge(8), "hemodynamic_device"
            ]
            .value_counts()
            .items()
        },
    }
    (output / "extraction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary
