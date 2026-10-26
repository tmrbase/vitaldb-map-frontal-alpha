"""Construct event-centred EEG/MAP subwindows and interval summaries."""

from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm import tqdm

from .config import ensure_output_directories
from .data import import_vitaldb, normalize_track_name, vital_path
from .maintenance_ra_extraction import (
    _art_ten_second_series,
    _art_window_qc,
    _condition_fraction_native,
    _low_map_events,
    _native_median,
    _native_records,
    _positive_summary_native,
    _raw_eeg_qc,
    _series_window_fraction,
    _series_window_median,
    _track_cadence,
    _volatile_classification,
    _zoh_rate_grid,
)
from .pkpd import (
    effect_site_concentration_from_rate,
    marsh_propofol_model,
    minto_remifentanil_model,
)
from .signal import (
    build_welch_psd_grid,
    filter_valid_blocks,
    spectral_endpoints_at_centers_from_grid,
)

EVENT_INTERVALS = (
    ("lead_m10_m5", -600.0, -300.0),
    ("baseline_m5_m2", -300.0, -120.0),
    ("response_0_1", 0.0, 60.0),
    ("response_1_3", 60.0, 180.0),
    ("response_3_5", 180.0, 300.0),
    ("response_5_10", 300.0, 600.0),
)
RESPONSE_INTERVALS = tuple(name for name, _, _ in EVENT_INTERVALS if name.startswith("response"))


def actual_log_nonalpha(log_alpha_power: float, logit_relative_alpha: float) -> float:
    """Recover log non-alpha power from values measured in the same subwindow."""

    if not np.isfinite(log_alpha_power) or not np.isfinite(logit_relative_alpha):
        return float("nan")
    return float(log_alpha_power - logit_relative_alpha)


def _pseudo_event_times(
    centers: np.ndarray,
    medians: np.ndarray,
    real_events: list[float],
    *,
    maintenance_start: float,
    maintenance_stop: float,
    threshold: float = 65.0,
    exclusion_seconds: float = 900.0,
) -> list[float | None]:
    """Match each real event to the closest clean within-surgery reference epoch."""

    candidates: list[float] = []
    for candidate in centers - 5.0:
        candidate = float(candidate)
        if candidate < maintenance_start + 600.0 or candidate > maintenance_stop - 600.0:
            continue
        if any(abs(candidate - real) < exclusion_seconds for real in real_events):
            continue
        interval = (centers >= candidate - 300.0) & (centers < candidate + 600.0)
        local = medians[interval]
        if local.size < 90 or int(np.isfinite(local).sum()) < 72:
            continue
        if np.any(local[np.isfinite(local)] < float(threshold)):
            continue
        candidates.append(candidate)
    matches: list[float | None] = []
    used: list[float] = []
    for real in real_events:
        available = [
            value
            for value in candidates
            if all(abs(value - prior) >= exclusion_seconds for prior in used)
        ]
        if not available:
            matches.append(None)
            continue
        chosen = min(available, key=lambda value: (abs(value - real), value))
        used.append(chosen)
        matches.append(float(chosen))
    return matches


def _event_eligible(events: list[float], start: float, stop: float) -> list[float]:
    return [float(value) for value in events if value >= start + 600.0 and value <= stop - 600.0]


def _subwindow_endpoint(
    raw: np.ndarray,
    *,
    left: float,
    right: float,
    cfg: dict[str, Any],
) -> dict[str, Any]:
    spec = cfg["analysis"]
    signal_cfg = cfg["signal"]
    sample_hz = float(spec["sample_hz"])
    start = max(round(left * sample_hz), 0)
    stop = min(round(right * sample_hz), len(raw))
    expected = round((right - left) * sample_hz)
    local_raw = np.asarray(raw[start:stop], dtype=float)
    if len(local_raw) < expected:
        local_raw = np.pad(local_raw, (0, expected - len(local_raw)), constant_values=np.nan)
    eeg_qc = _raw_eeg_qc(local_raw, spec)
    filtered = filter_valid_blocks(
        local_raw,
        sample_hz=sample_hz,
        low_hz=float(signal_cfg["primary_band_hz"][0]),
        high_hz=float(signal_cfg["primary_band_hz"][1]),
        order=int(signal_cfg["filter_order"]),
        guard_seconds=float(signal_cfg["guard_seconds"]),
    )
    welch = signal_cfg["welch"]
    grid = build_welch_psd_grid(
        filtered.values,
        sample_hz=sample_hz,
        nperseg=int(welch["nperseg"]),
        noverlap=int(welch["noverlap"]),
        detrend=str(welch["detrend"]),
    )
    duration = round(right - left)
    center = duration // 2
    endpoint = spectral_endpoints_at_centers_from_grid(
        grid,
        [center],
        window_seconds=duration,
        alpha_band_hz=tuple(float(value) for value in cfg["maintenance_ra"]["alpha_band_hz"]),
        delta_band_hz=tuple(float(value) for value in cfg["maintenance_ra"]["delta_band_hz"]),
        relative_total_band_hz=tuple(
            float(value) for value in cfg["maintenance_ra"]["relative_total_band_hz"]
        ),
        minimum_valid_segments=int(spec["welch_minimum_valid_segments"]),
    )[center]
    return {**eeg_qc, **endpoint}


def _aggregate_interval(
    rows: pd.DataFrame,
    *,
    event_id: str,
    event_type: str,
    event_time: float,
    interval_name: str,
    relative_start: float,
    relative_stop: float,
    sr_times: np.ndarray,
    sr_values: np.ndarray,
    maintenance_start: float,
    maintenance_stop: float,
    spec: dict[str, Any],
) -> dict[str, Any]:
    planned = len(rows)
    required = math.ceil(planned * float(spec["event_interval_valid_fraction_min"]))
    primary = rows.loc[rows["primary_ra_eligible"].astype(bool)]
    signal = rows.loc[rows["eeg_signal_eligible"].astype(bool)]
    ra_eligible = len(primary) >= required
    signal_eligible = len(signal) >= required

    def mean_if(frame: pd.DataFrame, column: str, eligible: bool) -> float:
        values = pd.to_numeric(frame[column], errors="coerce")
        return float(values.mean()) if eligible and values.notna().any() else float("nan")

    relative_alpha = mean_if(primary, "relative_alpha", ra_eligible)
    logit_ra = (
        float(
            np.log(np.clip(relative_alpha, 1e-6, 1 - 1e-6) / np.clip(1 - relative_alpha, 1e-6, 1))
        )
        if np.isfinite(relative_alpha)
        else float("nan")
    )
    absolute_start = event_time + relative_start
    absolute_stop = event_time + relative_stop
    center = (absolute_start + absolute_stop) / 2.0
    duration = absolute_stop - absolute_start
    sr_coverage, sr_positive_fraction = _condition_fraction_native(
        sr_times,
        sr_values,
        center=center,
        seconds=duration,
        predicate=lambda values: values > float(spec["sr_presence_positive_threshold"]),
        maintenance_start=maintenance_start,
        maintenance_stop=maintenance_stop,
    )
    _, sr_positive_median = _positive_summary_native(
        sr_times,
        sr_values,
        center=center,
        seconds=duration,
        positive_threshold=float(spec["sr_presence_positive_threshold"]),
    )
    return {
        "event_id": event_id,
        "event_type": event_type,
        "event_time_seconds": event_time,
        "interval": interval_name,
        "relative_start_seconds": relative_start,
        "relative_stop_seconds": relative_stop,
        "planned_subwindows": planned,
        "required_subwindows": required,
        "primary_ra_subwindows": len(primary),
        "signal_eligible_subwindows": len(signal),
        "ra_interval_eligible": bool(ra_eligible),
        "signal_interval_eligible": bool(signal_eligible),
        "relative_alpha": relative_alpha,
        "logit_relative_alpha": logit_ra,
        "log_alpha_power": mean_if(primary, "log_alpha_power", ra_eligible),
        "log_nonalpha_power": mean_if(primary, "log_nonalpha_power", ra_eligible),
        "log_delta_power": mean_if(primary, "log_delta_power", ra_eligible),
        "log_adr": mean_if(primary, "log_adr", ra_eligible),
        "sr_coverage": sr_coverage,
        "sr_positive_fraction": sr_positive_fraction,
        "sr_positive_median": sr_positive_median,
        "suppression_interval_eligible": bool(
            signal_eligible and sr_coverage >= float(spec["sr_window_fraction_min"])
        ),
        "suppression_present": bool(
            signal_eligible
            and sr_coverage >= float(spec["sr_window_fraction_min"])
            and sr_positive_fraction >= float(spec["sr_presence_window_fraction_min"])
        ),
        "propofol_ce": mean_if(rows, "propofol_ce", True),
        "remifentanil_ce": mean_if(rows, "remifentanil_ce", True),
        "etco2": mean_if(rows, "etco2", True),
        "body_temperature": mean_if(rows, "body_temperature", True),
        "art_map": mean_if(rows, "art_map", True),
    }


def extract_event_case(
    cfg: dict[str, Any], case: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    paths = ensure_output_directories(cfg)
    spec = cfg["analysis"]
    tracks = cfg["tracks"]
    requested = [
        tracks["eeg1"],
        tracks["sqi"],
        tracks["sr"],
        tracks["art_map"],
        tracks["body_temperature"],
        tracks["ppf_rate"],
        tracks["ppf_ce"],
        tracks["rftn20_rate"],
        tracks["rftn50_rate"],
        tracks["rftn20_ce"],
        tracks["rftn50_ce"],
        *tracks["etco2"],
        *tracks["inhaled_expired"],
    ]
    vitaldb = import_vitaldb()
    vital = vitaldb.VitalFile(
        str(vital_path(paths.vital_root, int(case["caseid"]))), track_names=requested
    )
    native = {name: _native_records(vital, name) for name in requested if name != tracks["eeg1"]}
    loaded = {normalize_track_name(name): name for name in vital.trks}
    eeg_loaded = loaded.get(tracks["eeg1"])
    if eeg_loaded is None:
        raise RuntimeError("EEG track is unavailable")
    matrix = np.asarray(
        vital.to_numpy([eeg_loaded], interval=1.0 / float(spec["sample_hz"])), dtype=float
    )
    raw = matrix.reshape(-1) if matrix.ndim == 1 else matrix[:, -1]

    maintenance_start = max(
        float(case["anestart"]) + float(spec["minutes_after_anestart"]) * 60.0,
        float(case["opstart"]) + float(spec["minutes_after_opstart"]) * 60.0,
        float(spec["pk_burn_in_after_record_start_minutes"]) * 60.0,
    )
    maintenance_stop = float(case["opend"]) - float(spec["minutes_before_opend"]) * 60.0
    if maintenance_stop <= maintenance_start:
        raise RuntimeError("No maintenance interval remains")
    volatile_state, volatile_audit = _volatile_classification(
        native, start=maintenance_start, stop=maintenance_stop, spec=spec
    )
    if not volatile_audit["broad_no_positive_volatile"]:
        raise RuntimeError("Case is no longer in the broad no-positive-volatile cohort")

    ppf_rate, ppf_observed, _ = _zoh_rate_grid(
        *native[tracks["ppf_rate"]],
        stop_seconds=maintenance_stop,
        maximum_gap_seconds=float(spec["infusion_max_gap_seconds"]),
    )
    ppf_ce = effect_site_concentration_from_rate(
        ppf_rate * float(spec["propofol_concentration_mg_ml"]) / 60.0,
        model=marsh_propofol_model(
            float(case["weight"]), ke0_per_minute=float(spec["marsh_ke0_per_minute"])
        ),
    )
    remi_inputs: list[np.ndarray] = []
    remi_observed_masks: list[np.ndarray] = []
    for name, concentration in spec["remifentanil_concentrations_ug_ml"].items():
        times, values = native[name]
        if not times.size:
            continue
        rate, observed, _ = _zoh_rate_grid(
            times,
            values,
            stop_seconds=maintenance_stop,
            maximum_gap_seconds=float(spec["infusion_max_gap_seconds"]),
        )
        remi_inputs.append(rate * float(concentration) / 60.0)
        remi_observed_masks.append(observed)
    if not remi_inputs:
        raise RuntimeError("No remifentanil rate is available")
    remi_input = np.nansum(np.vstack(remi_inputs), axis=0)
    remi_input[np.all(~np.isfinite(np.vstack(remi_inputs)), axis=0)] = np.nan
    remi_observed = np.logical_and.reduce(remi_observed_masks)
    remi_ce = effect_site_concentration_from_rate(
        remi_input,
        model=minto_remifentanil_model(
            age_years=float(case["age"]),
            sex=str(case["sex"]),
            weight_kg=float(case["weight"]),
            height_cm=float(case["height"]),
        ),
    )

    art_times, art_values = native[tracks["art_map"]]
    art_cadence = _track_cadence(art_times, start=maintenance_start, stop=maintenance_stop)
    art_centers, art_medians = _art_ten_second_series(
        art_times,
        art_values,
        start=maintenance_start,
        stop=maintenance_stop,
        cadence=art_cadence,
        spec=spec,
    )
    real_events = _event_eligible(
        _low_map_events(
            art_centers,
            art_medians,
            threshold=float(spec["event_map_threshold_mmhg"]),
            minimum_low_seconds=float(spec["event_minimum_low_seconds"]),
            separation_seconds=float(spec["event_separation_seconds"]),
            maximum_events=int(spec["event_maximum_per_surgery"]),
        ),
        maintenance_start,
        maintenance_stop,
    )
    pseudo_events = _pseudo_event_times(
        art_centers,
        art_medians,
        real_events,
        maintenance_start=maintenance_start,
        maintenance_stop=maintenance_stop,
        threshold=float(spec["event_map_threshold_mmhg"]),
        exclusion_seconds=float(spec["event_pseudo_exclusion_seconds"]),
    )

    sqi_times, sqi_values = native[tracks["sqi"]]
    sr_times, sr_values = native[tracks["sr"]]
    etco2_name = max(
        tracks["etco2"],
        key=lambda name: int(
            np.sum((native[name][0] >= maintenance_start) & (native[name][0] <= maintenance_stop))
        ),
    )
    etco2_times, etco2_values = native[etco2_name]
    temperature_times, temperature_values = native[tracks["body_temperature"]]
    subwindow_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    event_pairs = [("real", value, index) for index, value in enumerate(real_events, start=1)]
    event_pairs += [
        ("pseudo", value, index)
        for index, value in enumerate(pseudo_events, start=1)
        if value is not None
    ]
    subwindow_seconds = float(spec["event_subwindow_seconds"])
    for event_type, event_time, event_number in event_pairs:
        event_id = f"{int(case['caseid'])}:{event_type}:{event_number}"
        for interval_name, relative_start, relative_stop in EVENT_INTERVALS:
            centers = np.arange(
                event_time + relative_start + subwindow_seconds / 2.0,
                event_time + relative_stop,
                subwindow_seconds,
            )
            current_rows: list[dict[str, Any]] = []
            for center in centers:
                left = float(center - subwindow_seconds / 2.0)
                right = float(center + subwindow_seconds / 2.0)
                endpoint = _subwindow_endpoint(raw, left=left, right=right, cfg=cfg)
                # Calculate actual non-alpha power before interval aggregation.
                # In general, mean(log(A) - logit(RA)) is not equal to
                # mean(log(A)) - logit(mean(RA)).
                log_nonalpha_power = actual_log_nonalpha(
                    endpoint["log_alpha_power"], endpoint["logit_relative_alpha"]
                )
                sqi_coverage, sqi_fraction = _condition_fraction_native(
                    sqi_times,
                    sqi_values,
                    center=float(center),
                    seconds=subwindow_seconds,
                    predicate=lambda values: values >= float(spec["sqi_threshold"]),
                    maintenance_start=maintenance_start,
                    maintenance_stop=maintenance_stop,
                )
                sr_coverage, sr_lt20_fraction = _condition_fraction_native(
                    sr_times,
                    sr_values,
                    center=float(center),
                    seconds=subwindow_seconds,
                    predicate=lambda values: values < float(spec["sr_threshold"]),
                    maintenance_start=maintenance_start,
                    maintenance_stop=maintenance_stop,
                )
                signal_eligible = bool(
                    endpoint["eeg_raw_qc"]
                    and endpoint["welch_valid_segments"]
                    >= int(spec["welch_minimum_valid_segments"])
                    and sqi_coverage >= float(spec["sqi_window_fraction_min"])
                    and sqi_fraction >= float(spec["sqi_window_fraction_min"])
                )
                primary_ra_eligible = bool(
                    signal_eligible
                    and sr_coverage >= float(spec["sr_window_fraction_min"])
                    and sr_lt20_fraction >= float(spec["sr_window_fraction_min"])
                )
                art_qc = _art_window_qc(
                    art_times,
                    art_values,
                    center=float(center),
                    seconds=subwindow_seconds,
                    cadence=art_cadence,
                    spec=spec,
                )
                row = {
                    "caseid": int(case["caseid"]),
                    "subjectid": int(case["subjectid"]),
                    "event_id": event_id,
                    "event_type": event_type,
                    "event_number": int(event_number),
                    "event_time_seconds": float(event_time),
                    "interval": interval_name,
                    "subwindow_center_seconds": float(center),
                    "relative_center_seconds": float(center - event_time),
                    **endpoint,
                    "log_nonalpha_power": log_nonalpha_power,
                    "sqi_coverage": sqi_coverage,
                    "sqi_ge80_fraction": sqi_fraction,
                    "sr_coverage": sr_coverage,
                    "sr_lt20_fraction": sr_lt20_fraction,
                    "eeg_signal_eligible": signal_eligible,
                    "primary_ra_eligible": primary_ra_eligible,
                    **art_qc,
                    "propofol_ce": _series_window_median(ppf_ce, float(center), subwindow_seconds),
                    "remifentanil_ce": _series_window_median(
                        remi_ce, float(center), subwindow_seconds
                    ),
                    "propofol_rate_observed_fraction": _series_window_fraction(
                        ppf_observed, float(center), subwindow_seconds
                    ),
                    "remifentanil_rate_observed_fraction": _series_window_fraction(
                        remi_observed, float(center), subwindow_seconds
                    ),
                    "etco2": _native_median(
                        etco2_times, etco2_values, float(center), subwindow_seconds
                    ),
                    "body_temperature": _native_median(
                        temperature_times,
                        temperature_values,
                        float(center),
                        subwindow_seconds,
                    ),
                    "volatile_state": volatile_state,
                    "strict_measured_zero_volatile": bool(
                        volatile_audit["strict_measured_zero_volatile"]
                    ),
                }
                current_rows.append(row)
                subwindow_rows.append(row)
            interval_rows.append(
                {
                    "caseid": int(case["caseid"]),
                    "subjectid": int(case["subjectid"]),
                    "event_number": int(event_number),
                    "strict_measured_zero_volatile": bool(
                        volatile_audit["strict_measured_zero_volatile"]
                    ),
                    **_aggregate_interval(
                        pd.DataFrame(current_rows),
                        event_id=event_id,
                        event_type=event_type,
                        event_time=float(event_time),
                        interval_name=interval_name,
                        relative_start=relative_start,
                        relative_stop=relative_stop,
                        sr_times=sr_times,
                        sr_values=sr_values,
                        maintenance_start=maintenance_start,
                        maintenance_stop=maintenance_stop,
                        spec=spec,
                    ),
                }
            )
    audit = {
        "caseid": int(case["caseid"]),
        "subjectid": int(case["subjectid"]),
        "status": "ok",
        "real_events_detected_before_boundary_gate": len(
            _low_map_events(
                art_centers,
                art_medians,
                threshold=float(spec["event_map_threshold_mmhg"]),
                minimum_low_seconds=float(spec["event_minimum_low_seconds"]),
                separation_seconds=float(spec["event_separation_seconds"]),
                maximum_events=int(spec["event_maximum_per_surgery"]),
            )
        ),
        "real_events_with_full_lead_and_response": len(real_events),
        "pseudo_events_matched": int(sum(value is not None for value in pseudo_events)),
        "volatile_state": volatile_state,
    }
    return pd.DataFrame(subwindow_rows), pd.DataFrame(interval_rows), audit


def _extract_event_one(
    cfg: dict[str, Any], case: dict[str, Any], shard_text: str, force: bool
) -> dict[str, Any]:
    shard = Path(shard_text)
    interval_shard = shard.with_name(shard.stem + "_intervals.parquet")
    audit_path = shard.with_suffix(".audit.json")
    if shard.exists() and interval_shard.exists() and audit_path.exists() and not force:
        return json.loads(audit_path.read_text(encoding="utf-8"))
    try:
        subwindows, intervals, audit = extract_event_case(cfg, case)
        shard.parent.mkdir(parents=True, exist_ok=True)
        for frame, target in ((subwindows, shard), (intervals, interval_shard)):
            temporary = target.with_suffix(f".partial-{os.getpid()}.parquet")
            frame.to_parquet(temporary, index=False)
            temporary.replace(target)
    except Exception as error:  # noqa: BLE001 - every case failure is recorded in the audit
        audit = {
            "caseid": int(case["caseid"]),
            "subjectid": int(case["subjectid"]),
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error),
        }
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    return audit


def _attach_baseline(intervals: pd.DataFrame) -> pd.DataFrame:
    baseline_columns = [
        "event_id",
        "ra_interval_eligible",
        "relative_alpha",
        "logit_relative_alpha",
        "log_alpha_power",
        "log_nonalpha_power",
        "log_delta_power",
        "log_adr",
        "suppression_interval_eligible",
        "suppression_present",
        "propofol_ce",
        "remifentanil_ce",
        "etco2",
        "body_temperature",
        "art_map",
    ]
    baseline = intervals.loc[intervals["interval"].eq("baseline_m5_m2"), baseline_columns].copy()
    baseline = baseline.rename(
        columns={name: f"baseline_{name}" for name in baseline_columns if name != "event_id"}
    )
    response = intervals.loc[intervals["interval"].isin(RESPONSE_INTERVALS)].copy()
    response = response.merge(baseline, on="event_id", how="left", validate="many_to_one")
    response["event_primary_complete"] = (
        response["event_type"].eq("real")
        & response["ra_interval_eligible"].astype(bool)
        & response["baseline_ra_interval_eligible"].fillna(False).astype(bool)
        & response[
            [
                "baseline_logit_relative_alpha",
                "baseline_propofol_ce",
                "baseline_remifentanil_ce",
                "baseline_etco2",
                "baseline_body_temperature",
            ]
        ]
        .notna()
        .all(axis=1)
    )
    return response


def _availability_at_fraction(
    intervals: pd.DataFrame,
    *,
    fraction: float,
) -> dict[str, Any]:
    """Summarize outcome availability without examining any MAP--EEG association."""

    frame = intervals.loc[intervals["event_type"].eq("real")].copy()
    required = np.ceil(
        pd.to_numeric(frame["planned_subwindows"], errors="coerce") * float(fraction)
    )
    frame["available"] = pd.to_numeric(frame["primary_ra_subwindows"], errors="coerce").ge(required)
    baseline = frame.loc[frame["interval"].eq("baseline_m5_m2"), ["event_id", "available"]].rename(
        columns={"available": "baseline_available"}
    )
    response = frame.loc[frame["interval"].isin(RESPONSE_INTERVALS)].merge(
        baseline, on="event_id", how="left", validate="many_to_one"
    )
    eligible = response.loc[response["available"] & response["baseline_available"].fillna(False)]
    return {
        "minimum_interval_fraction": float(fraction),
        "events_with_eligible_baseline_and_at_least_one_response": int(
            eligible["event_id"].nunique()
        ),
        "surgeries_with_eligible_baseline_and_at_least_one_response": int(
            eligible["caseid"].nunique()
        ),
        "events_with_all_four_responses": int(
            eligible.groupby("event_id", observed=True).size().eq(4).sum()
        ),
        "eligible_response_rows": len(eligible),
        "eligible_rows_by_interval": {
            str(key): int(value)
            for key, value in eligible.groupby("interval", observed=True).size().items()
        },
    }


def prepare_maintenance_ra_full_events(
    cfg: dict[str, Any], *, workers: int | None = None, force: bool = False
) -> dict[str, Any]:
    """Create the event and lag datasets for the full analysis cohort."""

    output = ensure_output_directories(cfg).output_root / "maintenance_ra_full_v1"
    spec = cfg["analysis"]
    cases = pd.read_parquet(output / "eligible_cases.parquet")
    audit = pd.read_csv(output / "case_audit.csv")
    event_caseids = audit.loc[
        audit["status"].eq("ok")
        & pd.to_numeric(audit["low_map_event_count"], errors="coerce").fillna(0).gt(0),
        "caseid",
    ].astype(int)
    cases = cases.loc[cases["caseid"].astype(int).isin(event_caseids)].copy()
    shards = output / "event_shards"
    jobs = [
        (row._asdict(), shards / f"case_{int(row.caseid):04d}_subwindows.parquet")
        for row in cases.itertuples(index=False)
    ]
    audits: list[dict[str, Any]] = []
    worker_count = max(1, int(workers or cfg["runtime"]["workers"]))
    if worker_count == 1:
        for case, shard in tqdm(jobs, desc="maintenance RA full events"):
            audits.append(_extract_event_one(cfg, case, str(shard), force))
    else:
        with ProcessPoolExecutor(max_workers=min(worker_count, len(jobs))) as pool:
            futures = {
                pool.submit(_extract_event_one, cfg, case, str(shard), force): int(case["caseid"])
                for case, shard in jobs
            }
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="maintenance RA full events"
            ):
                audits.append(future.result())
    audit_frame = pd.DataFrame(audits).sort_values(["subjectid", "caseid"])
    audit_frame.to_csv(output / "event_case_audit.csv", index=False)
    good_caseids = audit_frame.loc[audit_frame["status"].eq("ok"), "caseid"].astype(int)
    subwindow_paths = [shards / f"case_{caseid:04d}_subwindows.parquet" for caseid in good_caseids]
    interval_paths = [
        shards / f"case_{caseid:04d}_subwindows_intervals.parquet" for caseid in good_caseids
    ]
    subwindows = pd.concat(
        [pd.read_parquet(path) for path in subwindow_paths if path.exists()], ignore_index=True
    )
    intervals = pd.concat(
        [pd.read_parquet(path) for path in interval_paths if path.exists()], ignore_index=True
    )
    responses = _attach_baseline(intervals)
    subwindows.to_parquet(output / "event_subwindows.parquet", index=False)
    intervals.to_parquet(output / "event_intervals.parquet", index=False)
    responses.to_parquet(output / "event_responses.parquet", index=False)
    subwindows.to_csv(output / "event_subwindows.csv", index=False)
    intervals.to_csv(output / "event_intervals.csv", index=False)
    responses.to_csv(output / "event_responses.csv", index=False)
    real_intervals = intervals.loc[intervals["event_type"].eq("real")]
    real_responses = responses.loc[responses["event_type"].eq("real")]
    summary = {
        "status": "full-cohort event/lag dataset",
        "candidate_event_cases": len(cases),
        "cases_processed_ok": int(audit_frame["status"].eq("ok").sum()),
        "cases_failed": int(audit_frame["status"].ne("ok").sum()),
        "real_events_with_full_lead_and_response": int(real_intervals["event_id"].nunique()),
        "pseudo_events_matched": int(
            intervals.loc[intervals["event_type"].eq("pseudo"), "event_id"].nunique()
        ),
        "real_events_with_eligible_baseline": int(
            real_intervals.loc[
                real_intervals["interval"].eq("baseline_m5_m2")
                & real_intervals["ra_interval_eligible"].astype(bool),
                "event_id",
            ].nunique()
        ),
        "real_events_with_at_least_one_primary_response": int(
            real_responses.loc[
                real_responses["event_primary_complete"].astype(bool), "event_id"
            ].nunique()
        ),
        "real_events_with_all_four_primary_responses": int(
            real_responses.loc[real_responses["event_primary_complete"].astype(bool)]
            .groupby("event_id", observed=True)
            .size()
            .eq(4)
            .sum()
        ),
        "primary_complete_response_rows": int(real_responses["event_primary_complete"].sum()),
        "availability_sensitivity_80_percent": _availability_at_fraction(
            intervals,
            fraction=float(spec["event_interval_valid_fraction_sensitivity"]),
        ),
    }
    (output / "event_flow.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary
