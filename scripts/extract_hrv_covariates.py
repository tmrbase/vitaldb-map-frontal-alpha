"""Extract window- and event-level ECG heart-rate-variability covariates."""

from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from wfdb.processing import xqrs_detect

from vitaldb_eeg.config import ensure_output_directories, load_config
from vitaldb_eeg.data import import_vitaldb, vital_path

ECG_TRACKS = ("SNUADC/ECG_II", "SNUADC/ECG_V5")
MAIN_WINDOW_SECONDS = 300.0
SHORT_WINDOW_SECONDS = 60.0
MINIMUM_ECG_COVERAGE = 0.90
MAX_INTERPOLATED_GAP_SECONDS = 0.20
MINIMUM_VALID_RR_FRACTION = 0.80
RR_RANGE_SECONDS = (0.30, 2.00)
LOCAL_RR_RATIO_RANGE = (0.80, 1.20)
EVENT_INTERVALS = {
    "baseline_m5_m2": (-300.0, -120.0),
    "response_0_1": (0.0, 60.0),
    "response_1_3": (60.0, 180.0),
    "response_3_5": (180.0, 300.0),
    "response_5_10": (300.0, 600.0),
}


def _output_root(cfg: dict[str, Any]) -> Path:
    return ensure_output_directories(cfg).output_root / "maintenance_ra_hrv_sensitivity_v1"


def _fill_short_gaps(signal: np.ndarray, maximum_samples: int) -> np.ndarray:
    output = np.asarray(signal, dtype=float).copy()
    finite = np.isfinite(output)
    if finite.all() or not finite.any():
        return output
    changes = np.diff(np.r_[False, ~finite, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    for start, stop in zip(starts, stops, strict=True):
        if (
            stop - start <= int(maximum_samples)
            and start > 0
            and stop < len(output)
            and np.isfinite(output[start - 1])
            and np.isfinite(output[stop])
        ):
            output[start:stop] = np.linspace(output[start - 1], output[stop], stop - start + 2)[
                1:-1
            ]
    return output


def _finite_segments(signal: np.ndarray, minimum_samples: int) -> list[tuple[int, int]]:
    finite = np.isfinite(signal)
    changes = np.diff(np.r_[False, finite, False].astype(np.int8))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return [
        (int(start), int(stop))
        for start, stop in zip(starts, stops, strict=True)
        if stop - start >= int(minimum_samples)
    ]


def _materialize(
    track: Any,
    *,
    vital_start: float,
    record_starts: np.ndarray,
    record_stops: np.ndarray,
    start: float,
    stop: float,
) -> np.ndarray:
    fs = float(track.srate)
    length = max(round((float(stop) - float(start)) * fs), 1)
    signal = np.full(length, np.nan, dtype=float)
    first = int(np.searchsorted(record_stops, float(start), side="right"))
    last = int(np.searchsorted(record_starts, float(stop), side="left"))
    for record in track.recs[first:last]:
        relative = float(record["dt"] - vital_start)
        values = np.asarray(record["val"])
        source_start = max(0, round((float(start) - relative) * fs))
        source_stop = min(len(values), round((float(stop) - relative) * fs))
        if source_stop <= source_start:
            continue
        destination_start = round((relative - float(start)) * fs) + source_start
        destination_stop = destination_start + source_stop - source_start
        if destination_start < 0:
            source_start -= destination_start
            destination_start = 0
        if destination_stop > length:
            source_stop -= destination_stop - length
            destination_stop = length
        if destination_stop <= destination_start:
            continue
        signal[destination_start:destination_stop] = values[source_start:source_stop].astype(
            float
        ) * float(track.gain) + float(track.offset)
    return signal


def _detect_rr(signal: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray, int]:
    filled = _fill_short_gaps(
        signal, maximum_samples=max(round(MAX_INTERPOLATED_GAP_SECONDS * fs), 1)
    )
    segments = _finite_segments(filled, minimum_samples=max(round(10.0 * fs), 1))
    rr_values: list[np.ndarray] = []
    rr_midpoints: list[np.ndarray] = []
    detected_beats = 0
    for start, stop in segments:
        local = filled[start:stop]
        try:
            peaks = np.asarray(xqrs_detect(sig=local, fs=fs, verbose=False), dtype=int)
        except Exception:  # noqa: BLE001, S112 - a failed ECG segment is skipped
            continue
        detected_beats += len(peaks)
        if len(peaks) < 3:
            continue
        absolute_peaks = peaks + int(start)
        rr = np.diff(absolute_peaks).astype(float) / fs
        midpoints = (absolute_peaks[1:] + absolute_peaks[:-1]) / (2.0 * fs)
        rr_values.append(rr)
        rr_midpoints.append(midpoints)
    if not rr_values:
        return np.asarray([], dtype=float), np.asarray([], dtype=float), detected_beats
    return np.concatenate(rr_values), np.concatenate(rr_midpoints), detected_beats


def _rr_metrics(
    rr: np.ndarray,
    midpoints: np.ndarray,
    *,
    start_seconds: float,
    stop_seconds: float,
    coverage: float,
    detected_beats: int,
) -> dict[str, Any]:
    selected = (midpoints >= float(start_seconds)) & (midpoints < float(stop_seconds))
    local_rr = np.asarray(rr[selected], dtype=float)
    duration = float(stop_seconds - start_seconds)
    minimum_nn = max(30, math.floor(duration * 0.40))
    if local_rr.size:
        local_median = median_filter(
            local_rr, size=min(11, max(3, local_rr.size | 1)), mode="nearest"
        )
        valid = (
            np.isfinite(local_rr)
            & (local_rr >= RR_RANGE_SECONDS[0])
            & (local_rr <= RR_RANGE_SECONDS[1])
            & (local_rr >= LOCAL_RR_RATIO_RANGE[0] * local_median)
            & (local_rr <= LOCAL_RR_RATIO_RANGE[1] * local_median)
        )
    else:
        valid = np.asarray([], dtype=bool)
    valid_rr = local_rr[valid]
    adjacent = valid[1:] & valid[:-1] if valid.size >= 2 else np.asarray([], dtype=bool)
    successive = np.diff(local_rr)[adjacent] if adjacent.size else np.asarray([], dtype=float)
    valid_fraction = float(valid.mean()) if valid.size else 0.0
    eligible = bool(
        float(coverage) >= MINIMUM_ECG_COVERAGE
        and len(valid_rr) >= minimum_nn
        and len(successive) >= max(minimum_nn - 2, 1)
        and valid_fraction >= MINIMUM_VALID_RR_FRACTION
    )
    rmssd_ms = float(np.sqrt(np.mean(successive**2)) * 1000.0) if eligible else math.nan
    sdnn_ms = float(np.std(valid_rr, ddof=1) * 1000.0) if eligible else math.nan
    qrs_hr = float(60.0 / np.median(valid_rr)) if eligible else math.nan
    return {
        "eligible": eligible,
        "rmssd_ms": rmssd_ms,
        "log_rmssd_ms": float(np.log(rmssd_ms)) if eligible and rmssd_ms > 0 else math.nan,
        "sdnn_ms": sdnn_ms,
        "log_sdnn_ms": float(np.log(sdnn_ms)) if eligible and sdnn_ms > 0 else math.nan,
        "qrs_hr_bpm": qrs_hr,
        "valid_rr": len(valid_rr),
        "successive_rr_pairs": len(successive),
        "valid_rr_fraction": valid_fraction,
        "detected_beats": int(detected_beats),
        "minimum_required_nn": int(minimum_nn),
        "duration_seconds": duration,
        "ecg_coverage": float(coverage),
    }


def _measure_window(
    track: Any,
    *,
    vital_start: float,
    record_starts: np.ndarray,
    record_stops: np.ndarray,
    start: float,
    stop: float,
) -> dict[str, Any]:
    signal = _materialize(
        track,
        vital_start=vital_start,
        record_starts=record_starts,
        record_stops=record_stops,
        start=start,
        stop=stop,
    )
    coverage = float(np.isfinite(signal).mean())
    rr, midpoints, detected = _detect_rr(signal, float(track.srate))
    return _rr_metrics(
        rr,
        midpoints,
        start_seconds=0.0,
        stop_seconds=float(stop - start),
        coverage=coverage,
        detected_beats=detected,
    )


def _measure_main_window(
    track: Any,
    *,
    vital_start: float,
    record_starts: np.ndarray,
    record_stops: np.ndarray,
    center: float,
) -> dict[str, Any]:
    start = float(center) - MAIN_WINDOW_SECONDS
    stop = float(center)
    signal = _materialize(
        track,
        vital_start=vital_start,
        record_starts=record_starts,
        record_stops=record_stops,
        start=start,
        stop=stop,
    )
    finite = np.isfinite(signal)
    coverage_300 = float(finite.mean())
    short_start = max(len(signal) - round(SHORT_WINDOW_SECONDS * float(track.srate)), 0)
    coverage_60 = float(finite[short_start:].mean()) if short_start < len(signal) else 0.0
    rr, midpoints, detected = _detect_rr(signal, float(track.srate))
    long_metrics = _rr_metrics(
        rr,
        midpoints,
        start_seconds=0.0,
        stop_seconds=MAIN_WINDOW_SECONDS,
        coverage=coverage_300,
        detected_beats=detected,
    )
    short_metrics = _rr_metrics(
        rr,
        midpoints,
        start_seconds=MAIN_WINDOW_SECONDS - SHORT_WINDOW_SECONDS,
        stop_seconds=MAIN_WINDOW_SECONDS,
        coverage=coverage_60,
        detected_beats=detected,
    )
    output: dict[str, Any] = {}
    output.update({f"hrv300_{key}": value for key, value in long_metrics.items()})
    output.update({f"hrv60_{key}": value for key, value in short_metrics.items()})
    return output


def _extract_case(
    cfg: dict[str, Any],
    caseid: int,
    centers: list[float],
    event_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = ensure_output_directories(cfg)
    vitaldb = import_vitaldb()
    vital = vitaldb.VitalFile(
        str(vital_path(paths.vital_root, caseid)), track_names=list(ECG_TRACKS)
    )
    track_name = next(
        (name for name in ECG_TRACKS if name in vital.trks and len(vital.trks[name].recs)),
        None,
    )
    if track_name is None:
        return pd.DataFrame(), {
            "caseid": int(caseid),
            "status": "no_ecg_track",
            "main_windows": len(centers),
            "event_windows": int(len(event_rows) * len(EVENT_INTERVALS)),
        }
    track = vital.trks[track_name]
    fs = float(track.srate)
    record_starts = np.asarray(
        [float(record["dt"] - vital.dtstart) for record in track.recs], dtype=float
    )
    record_stops = record_starts + np.asarray(
        [len(np.asarray(record["val"])) / fs for record in track.recs], dtype=float
    )
    rows: list[dict[str, Any]] = []
    for center in centers:
        metrics = _measure_main_window(
            track,
            vital_start=float(vital.dtstart),
            record_starts=record_starts,
            record_stops=record_stops,
            center=float(center),
        )
        rows.append(
            {
                "caseid": int(caseid),
                "kind": "main",
                "target_id": f"{int(caseid)}:main:{float(center):.3f}",
                "event_type": "",
                "event_number": np.nan,
                "interval": "main",
                "time_seconds": float(center),
                "ecg_track": track_name,
                "ecg_sample_hz": fs,
                **metrics,
            }
        )
    for event in event_rows:
        onset = float(event["event_time_seconds"])
        for interval, (relative_start, relative_stop) in EVENT_INTERVALS.items():
            metrics = _measure_window(
                track,
                vital_start=float(vital.dtstart),
                record_starts=record_starts,
                record_stops=record_stops,
                start=onset + relative_start,
                stop=onset + relative_stop,
            )
            rows.append(
                {
                    "caseid": int(caseid),
                    "kind": "event",
                    "target_id": str(event["event_id"]),
                    "event_type": str(event["event_type"]),
                    "event_number": int(event["event_number"]),
                    "interval": interval,
                    "time_seconds": onset,
                    "ecg_track": track_name,
                    "ecg_sample_hz": fs,
                    **{f"hrv_event_{key}": value for key, value in metrics.items()},
                }
            )
    frame = pd.DataFrame(rows)
    main_mask = frame["kind"].eq("main")
    event_mask = frame["kind"].eq("event")
    return frame, {
        "caseid": int(caseid),
        "status": "ok",
        "ecg_track": track_name,
        "ecg_sample_hz": fs,
        "main_windows": len(centers),
        "main_hrv300_eligible": int(frame.loc[main_mask, "hrv300_eligible"].sum()),
        "main_hrv60_eligible": int(frame.loc[main_mask, "hrv60_eligible"].sum()),
        "event_windows": int(event_mask.sum()),
        "event_hrv_eligible": int(frame.loc[event_mask, "hrv_event_eligible"].sum())
        if "hrv_event_eligible" in frame.columns
        else 0,
    }


def _worker(
    cfg: dict[str, Any],
    caseid: int,
    centers: list[float],
    event_rows: list[dict[str, Any]],
    destination: str,
) -> dict[str, Any]:
    target = Path(destination)
    if target.exists():
        frame = pd.read_parquet(target)
        return {
            "caseid": int(caseid),
            "status": "cached",
            "rows": len(frame),
        }
    temporary = target.with_suffix(f".partial-{os.getpid()}.parquet")
    try:
        frame, audit = _extract_case(cfg, caseid, centers, event_rows)
        target.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(temporary, index=False)
        temporary.replace(target)
        return audit
    except Exception as error:  # noqa: BLE001 - every case failure is retained for audit
        if temporary.exists():
            temporary.unlink()
        return {
            "caseid": int(caseid),
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error),
        }


def extract_hrv(
    cfg: dict[str, Any],
    *,
    output: Path,
    workers: int,
    force: bool,
) -> dict[str, Any]:
    paths = ensure_output_directories(cfg)
    windows = pd.read_parquet(paths.output_root / "maintenance_ra_full_v1" / "full_windows.parquet")
    complete = windows.loc[windows["primary_complete_window"].fillna(False)].copy()
    counts = complete.groupby("caseid", observed=True).size()
    complete = complete.loc[complete["caseid"].isin(counts[counts >= 8].index)].copy()
    events = pd.read_parquet(
        paths.output_root / "maintenance_ra_full_v1" / "event_intervals.parquet"
    )
    event_keys = (
        events[
            ["caseid", "subjectid", "event_id", "event_type", "event_number", "event_time_seconds"]
        ]
        .drop_duplicates("event_id")
        .sort_values(["caseid", "event_type", "event_number"])
    )
    caseids = np.asarray(sorted(complete["caseid"].unique()), dtype=int)
    main_by_case = {
        int(caseid): local["time_seconds"].astype(float).tolist()
        for caseid, local in complete.loc[complete["caseid"].isin(caseids)].groupby(
            "caseid", observed=True
        )
    }
    event_by_case = {
        int(caseid): local.to_dict("records")
        for caseid, local in event_keys.loc[event_keys["caseid"].isin(caseids)].groupby(
            "caseid", observed=True
        )
    }
    shards = output / "case_shards"
    jobs: list[tuple[int, list[float], list[dict[str, Any]], Path]] = []
    for caseid in caseids:
        target = shards / f"case_{int(caseid):04d}.parquet"
        if force and target.exists():
            target.unlink()
        jobs.append(
            (
                int(caseid),
                main_by_case.get(int(caseid), []),
                event_by_case.get(int(caseid), []),
                target,
            )
        )
    audits: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max(1, min(int(workers), len(jobs)))) as pool:
        futures = {
            pool.submit(_worker, cfg, caseid, centers, event_rows, str(target)): caseid
            for caseid, centers, event_rows, target in jobs
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            audits.append(future.result())
            if completed % 25 == 0 or completed == len(jobs):
                print(f"[HRV extraction] {completed}/{len(jobs)} surgeries", flush=True)
    frames = [pd.read_parquet(target) for _, _, _, target in jobs if target.exists()]
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    output.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(output / "hrv_windows.parquet", index=False)
    combined.to_csv(output / "hrv_windows.csv", index=False)
    audit_frame = pd.DataFrame(audits).sort_values("caseid")
    audit_frame.to_csv(output / "extraction_job_audit.csv", index=False)
    main = combined.loc[combined["kind"].eq("main")].copy() if not combined.empty else combined
    event = combined.loc[combined["kind"].eq("event")].copy() if not combined.empty else combined
    summary = {
        "status": "post hoc ECG-derived HRV extraction",
        "surgeries": len(caseids),
        "surgeries_with_output": int(combined["caseid"].nunique()) if not combined.empty else 0,
        "lead_counts": {
            str(key): int(value) for key, value in main["ecg_track"].value_counts().items()
        }
        if not main.empty
        else {},
        "main_windows": len(main),
        "main_hrv300_eligible": int(main["hrv300_eligible"].sum()) if not main.empty else 0,
        "main_hrv300_eligible_fraction": float(main["hrv300_eligible"].mean())
        if not main.empty
        else math.nan,
        "main_hrv60_eligible": int(main["hrv60_eligible"].sum()) if not main.empty else 0,
        "main_hrv60_eligible_fraction": float(main["hrv60_eligible"].mean())
        if not main.empty
        else math.nan,
        "surgeries_with_8plus_hrv300_windows": int(
            main.loc[main["hrv300_eligible"].fillna(False)]
            .groupby("caseid", observed=True)
            .size()
            .ge(8)
            .sum()
        )
        if not main.empty
        else 0,
        "event_windows": len(event),
        "event_hrv_eligible": int(event["hrv_event_eligible"].sum()) if not event.empty else 0,
    }
    (output / "extraction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract heart-rate-variability covariates.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/analysis_public.yaml"),
        help="Analysis configuration file.",
    )
    parser.add_argument("--output", type=Path, help="Optional output CSV path.")
    parser.add_argument(
        "--workers", type=int, help="Parallel case workers; defaults to the configuration."
    )
    parser.add_argument("--force", action="store_true", help="Rebuild existing case-level shards.")
    args = parser.parse_args()
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be at least 1")
    cfg = load_config(args.config)
    output = args.output or _output_root(cfg)
    summary = extract_hrv(
        cfg,
        output=output,
        workers=int(args.workers or cfg["runtime"]["workers"]),
        force=args.force,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
