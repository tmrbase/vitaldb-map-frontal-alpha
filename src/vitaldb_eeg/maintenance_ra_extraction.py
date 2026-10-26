"""Build quality-controlled maintenance windows and the primary analytic cohort."""

from __future__ import annotations

import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import ensure_output_directories
from .data import clinical_frame, import_vitaldb, normalize_track_name, read_track_index, vital_path
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


def _output_root(cfg: dict[str, Any]) -> Path:
    return ensure_output_directories(cfg).output_root / "maintenance_ra_prescreen_v1"


def _has_any(available: frozenset[str], candidates: list[str]) -> bool:
    return any(name in available for name in candidates)


def _metadata_eligible_cases(cfg: dict[str, Any]) -> pd.DataFrame:
    """Build the result-independent sampling frame for the full analysis cohort."""

    paths = ensure_output_directories(cfg)
    clinical = clinical_frame(paths.clinical_data)
    track_rows = {row.caseid: row for row in read_track_index(paths.track_index)}
    tracks = cfg["tracks"]
    required = {
        tracks["eeg1"],
        tracks["sqi"],
        tracks["sr"],
        tracks["art_map"],
        tracks["ppf_rate"],
    }
    remi_rate_tracks = [tracks["rftn20_rate"], tracks["rftn50_rate"]]

    def track_eligible(caseid: int) -> bool:
        row = track_rows.get(int(caseid))
        if row is None or not vital_path(paths.vital_root, caseid).exists():
            return False
        available = row.tracks
        return required.issubset(available) and _has_any(available, remi_rate_tracks)

    numeric_columns = ["age", "height", "weight", "anestart", "opstart", "opend"]
    values = clinical[numeric_columns].apply(pd.to_numeric, errors="coerce")
    metadata_ok = values.notna().all(axis=1)
    metadata_ok &= values["age"].ge(float(cfg["cohort"]["adult_age"]))
    metadata_ok &= values["height"].gt(0) & values["weight"].gt(0)
    metadata_ok &= clinical["ane_type"].eq(cfg["cohort"]["primary_ane_type"])
    metadata_ok &= (
        clinical["sex"].astype(str).str.upper().isin(["M", "F", "MALE", "FEMALE", "0", "1"])
    )
    metadata_ok &= clinical["caseid"].map(track_eligible)
    duration_ok = clinical.apply(
        lambda row: (
            len(_maintenance_centers(row.to_dict(), cfg["analysis"]))
            >= int(cfg["analysis"]["minimum_valid_windows"])
        ),
        axis=1,
    )
    return clinical.loc[metadata_ok & duration_ok].copy().sort_values(["subjectid", "caseid"])


def _prescreen_volatile_case(cfg: dict[str, Any], case: dict[str, Any]) -> dict[str, Any]:
    paths = ensure_output_directories(cfg)
    spec = cfg["analysis"]
    tracks = cfg["tracks"]
    centers = _maintenance_centers(case, spec)
    start = float(centers[0] - spec["window_seconds"] / 2.0)
    stop = float(centers[-1] + spec["window_seconds"] / 2.0)
    vitaldb = import_vitaldb()
    vital = vitaldb.VitalFile(
        str(vital_path(paths.vital_root, int(case["caseid"]))),
        track_names=list(tracks["inhaled_expired"]),
    )
    native = {name: _native_records(vital, name) for name in tracks["inhaled_expired"]}
    _, audit = _volatile_classification(native, start=start, stop=stop, spec=spec)
    return {
        "caseid": int(case["caseid"]),
        "subjectid": int(case["subjectid"]),
        "scheduled_windows": len(centers),
        **audit,
    }


def prepare_maintenance_ra_prescreen(
    cfg: dict[str, Any], *, workers: int | None = None, force: bool = False
) -> dict[str, Any]:
    """Screen every metadata-eligible case for expired volatile anaesthetic."""

    metadata_cases = _metadata_eligible_cases(cfg)
    spec = cfg["analysis"]
    output = _output_root(cfg)
    output.mkdir(parents=True, exist_ok=True)
    prescreen_path = output / "volatile_prescreen.csv"
    prescreen_manifest_path = output / "volatile_prescreen_manifest.json"
    prescreen_payload = {
        "caseids": metadata_cases["caseid"].astype(int).tolist(),
        "volatile_positive_thresholds": spec["volatile_positive_thresholds"],
        "volatile_zero_tolerances": spec["volatile_zero_tolerances"],
        "volatile_valid_range": spec["volatile_valid_range"],
        "volatile_boundary_tolerance_seconds": spec["volatile_boundary_tolerance_seconds"],
        "volatile_max_gap_seconds": spec["volatile_max_gap_seconds"],
        "maintenance_window": {
            key: spec[key]
            for key in (
                "pk_burn_in_after_record_start_minutes",
                "minutes_after_anestart",
                "minutes_after_opstart",
                "minutes_before_opend",
                "grid_seconds",
                "window_seconds",
                "minimum_valid_windows",
            )
        },
    }
    prescreen_fingerprint = hashlib.sha256(
        json.dumps(prescreen_payload, sort_keys=True).encode()
    ).hexdigest()
    cached_fingerprint = None
    if prescreen_manifest_path.exists():
        cached_fingerprint = json.loads(prescreen_manifest_path.read_text(encoding="utf-8")).get(
            "fingerprint"
        )
    if prescreen_path.exists() and cached_fingerprint == prescreen_fingerprint and not force:
        prescreen = pd.read_csv(prescreen_path)
    else:
        jobs = [row._asdict() for row in metadata_cases.itertuples(index=False)]
        worker_count = max(1, min(int(workers or cfg["runtime"]["workers"]), len(jobs)))
        if worker_count == 1:
            records = [_prescreen_volatile_case(cfg, case) for case in jobs]
        else:
            with ProcessPoolExecutor(max_workers=worker_count) as pool:
                records = list(pool.map(_prescreen_volatile_case, [cfg] * len(jobs), jobs))
        prescreen = pd.DataFrame(records)
        prescreen.to_csv(prescreen_path, index=False)
        prescreen_manifest_path.write_text(
            json.dumps(
                {"fingerprint": prescreen_fingerprint, **prescreen_payload},
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
    return {
        "status": "volatile prescreen complete",
        "metadata_eligible_cases": len(prescreen),
        "broad_no_positive_volatile_cases": int(
            prescreen["broad_no_positive_volatile"].astype(bool).sum()
        ),
        "strict_measured_zero_volatile_cases": int(
            prescreen["strict_measured_zero_volatile"].astype(bool).sum()
        ),
        "raw_eeg_or_map_used_for_selection": False,
    }


def _maintenance_centers(case: dict[str, Any], spec: dict[str, Any]) -> np.ndarray:
    start = max(
        60.0 * float(spec["pk_burn_in_after_record_start_minutes"]),
        float(case["anestart"]) + 60.0 * float(spec["minutes_after_anestart"]),
        float(case["opstart"]) + 60.0 * float(spec["minutes_after_opstart"]),
    )
    stop = float(case["opend"]) - 60.0 * float(spec["minutes_before_opend"])
    half = float(spec["window_seconds"]) / 2.0
    if stop - half < start + half:
        return np.asarray([], dtype=float)
    return np.arange(start + half, stop - half + 1e-9, float(spec["grid_seconds"]), dtype=float)


def _native_records(vital: Any, track_name: str) -> tuple[np.ndarray, np.ndarray]:
    loaded = {normalize_track_name(name): name for name in vital.trks}
    key = loaded.get(track_name)
    if key is None:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    recs = vital.trks[key].recs
    if not recs:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    times = np.asarray([float(rec["dt"]) - float(vital.dtstart) for rec in recs], dtype=float)
    values = np.asarray([float(rec["val"]) for rec in recs], dtype=float)
    order = np.argsort(times, kind="stable")
    return times[order], values[order]


def _track_cadence(times: np.ndarray, *, start: float, stop: float) -> float:
    selected = np.asarray(times[(times >= start) & (times <= stop)], dtype=float)
    differences = np.diff(selected)
    differences = differences[np.isfinite(differences) & (differences > 0)]
    return float(np.median(differences)) if differences.size else float("nan")


def _window_native(
    times: np.ndarray, values: np.ndarray, center: float, seconds: float
) -> tuple[np.ndarray, np.ndarray]:
    half = float(seconds) / 2.0
    selected = (times >= center - half) & (times < center + half)
    return times[selected], values[selected]


def _longest_true_seconds(mask: np.ndarray, sample_hz: float) -> float:
    values = np.asarray(mask, dtype=bool)
    if values.size == 0 or not values.any():
        return 0.0
    padded = np.pad(values.astype(np.int8), (1, 1))
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return float(np.max(stops - starts) / float(sample_hz))


def _raw_eeg_qc(raw_window: np.ndarray, spec: dict[str, Any]) -> dict[str, Any]:
    values = np.asarray(raw_window, dtype=float)
    finite = np.isfinite(values)
    coverage = float(finite.mean()) if values.size else 0.0
    usable = values[finite]
    if usable.size < 2:
        return {
            "eeg_raw_coverage": coverage,
            "eeg_flatline_fraction": float("nan"),
            "eeg_flatline_longest_seconds": float("nan"),
            "eeg_saturation_fraction": float("nan"),
            "eeg_saturation_longest_seconds": float("nan"),
            "eeg_raw_qc": False,
        }
    pair_finite = finite[1:] & finite[:-1]
    flatline = np.zeros(max(values.size - 1, 0), dtype=bool)
    flatline[pair_finite] = np.abs(np.diff(values)[pair_finite]) <= float(
        spec["flatline_difference_absolute_max"]
    )
    saturated = finite & (np.abs(values) >= float(spec["saturation_absolute_min"]))
    flat_fraction = float(flatline.sum() / max(pair_finite.sum(), 1))
    saturation_fraction = float(saturated.sum() / max(finite.sum(), 1))
    flat_longest = _longest_true_seconds(flatline, float(spec["sample_hz"]))
    saturation_longest = _longest_true_seconds(saturated, float(spec["sample_hz"]))
    eligible = (
        coverage >= float(spec["eeg_raw_coverage_min"])
        and flat_fraction < float(spec["flatline_fraction_max"])
        and flat_longest < float(spec["flatline_longest_seconds_max"])
        and saturation_fraction < float(spec["saturation_fraction_max"])
        and saturation_longest < float(spec["saturation_longest_seconds_max"])
    )
    return {
        "eeg_raw_coverage": coverage,
        "eeg_flatline_fraction": flat_fraction,
        "eeg_flatline_longest_seconds": flat_longest,
        "eeg_saturation_fraction": saturation_fraction,
        "eeg_saturation_longest_seconds": saturation_longest,
        "eeg_raw_qc": bool(eligible),
    }


def _condition_fraction_native(
    times: np.ndarray,
    values: np.ndarray,
    *,
    center: float,
    seconds: float,
    predicate: Any,
    maintenance_start: float,
    maintenance_stop: float,
) -> tuple[float, float]:
    _, local_values = _window_native(times, values, center, seconds)
    cadence = _track_cadence(times, start=maintenance_start, stop=maintenance_stop)
    expected = max(round(float(seconds) / cadence), 1) if np.isfinite(cadence) else 0
    finite = np.isfinite(local_values)
    coverage = float(finite.sum() / expected) if expected else 0.0
    fraction = float(np.mean(predicate(local_values[finite]))) if finite.any() else 0.0
    return min(coverage, 1.0), fraction


def _positive_summary_native(
    times: np.ndarray,
    values: np.ndarray,
    *,
    center: float,
    seconds: float,
    positive_threshold: float,
) -> tuple[float, float]:
    """Return positive-value fraction and positive-part median in one native window."""

    _, local_values = _window_native(times, values, center, seconds)
    finite = local_values[np.isfinite(local_values)]
    if not finite.size:
        return 0.0, float("nan")
    positive = finite[finite > float(positive_threshold)]
    fraction = float(positive.size / finite.size)
    median = float(np.median(positive)) if positive.size else float("nan")
    return fraction, median


def _art_window_qc(
    times: np.ndarray,
    values: np.ndarray,
    *,
    center: float,
    seconds: float,
    cadence: float,
    spec: dict[str, Any],
) -> dict[str, Any]:
    local_times, local_values = _window_native(times, values, center, seconds)
    low, high = (float(value) for value in spec["art_valid_range_mmhg"])
    valid = np.isfinite(local_values) & (local_values >= low) & (local_values <= high)
    expected = round(float(seconds) / cadence) if np.isfinite(cadence) and cadence > 0 else 0
    required = max(
        int(spec["art_minimum_samples"]), math.ceil(expected * float(spec["art_coverage_min"]))
    )
    valid_times = local_times[valid]
    maximum_gap = float(np.max(np.diff(valid_times))) if valid_times.size >= 2 else float("inf")
    eligible = valid.sum() >= required and maximum_gap <= float(spec["art_max_gap_seconds"])
    return {
        "art_map": float(np.median(local_values[valid])) if valid.any() else float("nan"),
        "art_valid_samples": int(valid.sum()),
        "art_expected_samples": int(expected),
        "art_required_samples": int(required),
        "art_max_gap_seconds": maximum_gap,
        "art_qc": bool(eligible),
    }


def _native_median(times: np.ndarray, values: np.ndarray, center: float, seconds: float) -> float:
    _, local = _window_native(times, values, center, seconds)
    local = local[np.isfinite(local)]
    return float(np.median(local)) if local.size else float("nan")


def _art_ten_second_series(
    times: np.ndarray,
    values: np.ndarray,
    *,
    start: float,
    stop: float,
    cadence: float,
    spec: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(math.floor(start / 10.0) * 10.0, stop + 10.0, 10.0)
    centers = edges[:-1] + 5.0
    medians = np.full(centers.shape, np.nan, dtype=float)
    low, high = (float(value) for value in spec["art_valid_range_mmhg"])
    expected = max(round(10.0 / cadence), 1) if np.isfinite(cadence) and cadence > 0 else 0
    required = max(math.ceil(expected * float(spec["art_coverage_min"])), 1)
    for index, (left, right) in enumerate(pairwise(edges)):
        selected = (times >= left) & (times < right)
        local = values[selected]
        valid = np.isfinite(local) & (local >= low) & (local <= high)
        if int(valid.sum()) >= required:
            medians[index] = float(np.median(local[valid]))
    return centers, medians


def _preceding_map_features(
    centers: np.ndarray,
    medians: np.ndarray,
    *,
    endpoint: float,
    threshold: float = 65.0,
) -> dict[str, float]:
    selected = (centers >= endpoint - 300.0) & (centers < endpoint)
    values = medians[selected]
    valid = np.isfinite(values)
    if int(valid.sum()) < 24:
        return {
            "preceding_map_10s_bins": int(valid.sum()),
            "preceding_map_arv": float("nan"),
            "preceding_map_below65_auc_mmhg_min": float("nan"),
        }
    clean = values[valid]
    # ARV is defined only across adjacent 10-second bins.  Do not bridge a
    # missing bin, because that would turn a 20+ second change into a nominal
    # 10-second change and overstate short-term variability.
    adjacent = valid[:-1] & valid[1:]
    differences = np.abs(np.diff(values))[adjacent]
    arv = float(np.mean(differences)) if differences.size else float("nan")
    area = float(np.sum(np.maximum(float(threshold) - clean, 0.0)) * 10.0 / 60.0)
    return {
        "preceding_map_10s_bins": int(valid.sum()),
        "preceding_map_arv": arv,
        "preceding_map_below65_auc_mmhg_min": area,
    }


def _low_map_events(
    centers: np.ndarray,
    medians: np.ndarray,
    *,
    threshold: float = 65.0,
    minimum_low_seconds: float = 60.0,
    separation_seconds: float = 900.0,
    maximum_events: int = 2,
) -> list[float]:
    low = np.isfinite(medians) & (medians < float(threshold))
    padded = np.pad(low.astype(np.int8), (1, 1))
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    candidates: list[float] = []
    for start_index, stop_index in zip(starts, stops, strict=True):
        bin_seconds = float(np.median(np.diff(centers))) if len(centers) > 1 else 10.0
        minimum_bins = math.ceil(float(minimum_low_seconds) / bin_seconds)
        if stop_index - start_index < minimum_bins:
            continue
        event_time = float(centers[start_index] - 5.0)
        baseline = (
            (centers >= event_time - 300.0) & (centers < event_time - 120.0) & np.isfinite(medians)
        )
        if int(baseline.sum()) < 14:
            continue
        if float(np.median(medians[baseline])) < float(threshold):
            continue
        candidates.append(event_time)
    selected: list[float] = []
    for event_time in candidates:
        if selected and event_time - selected[-1] < float(separation_seconds):
            continue
        selected.append(event_time)
        if len(selected) == int(maximum_events):
            break
    return selected


def _volatile_classification(
    native: dict[str, tuple[np.ndarray, np.ndarray]],
    *,
    start: float,
    stop: float,
    spec: dict[str, Any],
) -> tuple[str, dict[str, Any]]:
    present = []
    positive = []
    strict_zero = []
    audit: dict[str, Any] = {}
    low, high = (float(value) for value in spec["volatile_valid_range"])
    for name, threshold in spec["volatile_positive_thresholds"].items():
        times, values = native.get(name, (np.asarray([]), np.asarray([])))
        if times.size == 0:
            audit[f"volatile_{name}_state"] = "absent"
            continue
        present.append(name)
        selected = (
            (times >= start)
            & (times <= stop)
            & np.isfinite(values)
            & (values >= low)
            & (values <= high)
        )
        local_times = times[selected]
        local_values = values[selected]
        if local_times.size == 0:
            audit[f"volatile_{name}_state"] = "invalid"
            continue
        span_ok = local_times[0] <= start + float(
            spec["volatile_boundary_tolerance_seconds"]
        ) and local_times[-1] >= stop - float(spec["volatile_boundary_tolerance_seconds"])
        max_gap = float(np.max(np.diff(local_times))) if local_times.size >= 2 else float("inf")
        coverage_ok = span_ok and max_gap <= float(spec["volatile_max_gap_seconds"])
        is_positive = bool(np.any(local_values >= float(threshold)))
        is_zero = coverage_ok and bool(
            np.all(np.abs(local_values) <= float(spec["volatile_zero_tolerances"][name]))
        )
        positive.append(is_positive)
        strict_zero.append(is_zero)
        state = (
            "positive"
            if is_positive
            else "measured_zero"
            if is_zero
            else "present_inactive_or_incomplete"
        )
        audit[f"volatile_{name}_state"] = state
        audit[f"volatile_{name}_max"] = float(np.max(local_values))
        audit[f"volatile_{name}_max_gap_seconds"] = max_gap
    if any(positive):
        state = "positive"
    elif not present:
        state = "track_absent"
    elif strict_zero and all(strict_zero) and len(strict_zero) == len(present):
        state = "measured_zero"
    else:
        state = "indeterminate"
    audit["volatile_state"] = state
    audit["broad_no_positive_volatile"] = state in {"track_absent", "measured_zero"}
    audit["strict_measured_zero_volatile"] = state == "measured_zero"
    return state, audit


def _zoh_rate_grid(
    times: np.ndarray,
    values: np.ndarray,
    *,
    stop_seconds: float,
    maximum_gap_seconds: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    grid = np.arange(0, math.ceil(stop_seconds) + 1, dtype=float)
    valid = np.isfinite(times) & np.isfinite(values) & (values >= 0)
    times = times[valid]
    values = values[valid]
    if times.size == 0:
        return (
            np.full(grid.shape, np.nan),
            np.zeros(grid.shape, dtype=bool),
            {"rate_coverage": 0.0, "rate_max_gap_seconds": float("inf")},
        )
    order = np.argsort(times, kind="stable")
    times = times[order]
    values = values[order]
    positions = np.searchsorted(times, grid, side="right") - 1
    observed = positions >= 0
    ages = np.full(grid.shape, np.inf, dtype=float)
    ages[observed] = grid[observed] - times[positions[observed]]
    coverage = observed & (ages <= float(maximum_gap_seconds))
    rates = np.full(grid.shape, np.nan, dtype=float)
    # Pump rates are state variables. Carry the latest rate through recorder gaps
    # for PK integration, but retain a separate observed-coverage mask so windows
    # crossing a long gap cannot enter the primary analysis.
    rates[observed] = values[positions[observed]]
    differences = np.diff(times)
    return (
        rates,
        coverage,
        {
            "rate_coverage": float(np.mean(coverage)),
            "rate_max_gap_seconds": float(np.max(differences))
            if differences.size
            else float("inf"),
        },
    )


def _series_window_median(values: np.ndarray, center: float, seconds: float) -> float:
    start = max(math.floor(center - seconds / 2.0), 0)
    stop = min(math.ceil(center + seconds / 2.0), len(values))
    local = values[start:stop]
    local = local[np.isfinite(local)]
    return float(np.median(local)) if local.size else float("nan")


def _series_window_fraction(values: np.ndarray, center: float, seconds: float) -> float:
    start = max(math.floor(center - seconds / 2.0), 0)
    stop = min(math.ceil(center + seconds / 2.0), len(values))
    local = np.asarray(values[start:stop], dtype=bool)
    return float(np.mean(local)) if local.size else 0.0


def extract_case(cfg: dict[str, Any], case: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    paths = ensure_output_directories(cfg)
    spec = cfg["analysis"]
    tracks = cfg["tracks"]
    signal_cfg = cfg["signal"]
    centers = _maintenance_centers(case, spec)
    maintenance_start = float(centers[0] - spec["window_seconds"] / 2.0)
    maintenance_stop = float(centers[-1] + spec["window_seconds"] / 2.0)
    requested = [
        tracks["eeg1"],
        tracks["sqi"],
        tracks["sr"],
        tracks["art_map"],
        tracks["body_temperature"],
        tracks["ppf_rate"],
        tracks["ppf_vol"],
        tracks["ppf_ce"],
        tracks["rftn20_rate"],
        tracks["rftn20_vol"],
        tracks["rftn50_rate"],
        tracks["rftn50_vol"],
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
        raise RuntimeError("EEG track disappeared after metadata selection")
    matrix = np.asarray(
        vital.to_numpy([eeg_loaded], interval=1.0 / float(spec["sample_hz"])), dtype=float
    )
    raw = matrix.reshape(-1) if matrix.ndim == 1 else matrix[:, -1]

    filtered = filter_valid_blocks(
        raw,
        sample_hz=float(spec["sample_hz"]),
        low_hz=float(signal_cfg["primary_band_hz"][0]),
        high_hz=float(signal_cfg["primary_band_hz"][1]),
        order=int(signal_cfg["filter_order"]),
        guard_seconds=float(signal_cfg["guard_seconds"]),
    )
    welch = signal_cfg["welch"]
    grid = build_welch_psd_grid(
        filtered.values,
        sample_hz=float(spec["sample_hz"]),
        nperseg=int(welch["nperseg"]),
        noverlap=int(welch["noverlap"]),
        detrend=str(welch["detrend"]),
    )
    endpoints = spectral_endpoints_at_centers_from_grid(
        grid,
        centers,
        window_seconds=int(spec["window_seconds"]),
        alpha_band_hz=tuple(float(value) for value in cfg["maintenance_ra"]["alpha_band_hz"]),
        delta_band_hz=tuple(float(value) for value in cfg["maintenance_ra"]["delta_band_hz"]),
        relative_total_band_hz=tuple(
            float(value) for value in cfg["maintenance_ra"]["relative_total_band_hz"]
        ),
        minimum_valid_segments=int(spec["welch_minimum_valid_segments"]),
    )

    volatile_state, volatile_audit = _volatile_classification(
        native, start=maintenance_start, stop=maintenance_stop, spec=spec
    )
    ppf_times, ppf_values = native[tracks["ppf_rate"]]
    ppf_rate, ppf_rate_observed, ppf_rate_audit = _zoh_rate_grid(
        ppf_times,
        ppf_values,
        stop_seconds=maintenance_stop,
        maximum_gap_seconds=float(spec["infusion_max_gap_seconds"]),
    )
    ppf_input = ppf_rate * float(spec["propofol_concentration_mg_ml"]) / 60.0
    ppf_ce = effect_site_concentration_from_rate(
        ppf_input,
        model=marsh_propofol_model(
            float(case["weight"]), ke0_per_minute=float(spec["marsh_ke0_per_minute"])
        ),
    )
    remi_inputs = []
    remi_audits = []
    remi_display_sources: list[tuple[np.ndarray, np.ndarray]] = []
    for name, concentration in spec["remifentanil_concentrations_ug_ml"].items():
        times, values = native[name]
        if times.size == 0:
            continue
        rates, observed, audit = _zoh_rate_grid(
            times,
            values,
            stop_seconds=maintenance_stop,
            maximum_gap_seconds=float(spec["infusion_max_gap_seconds"]),
        )
        remi_inputs.append(rates * float(concentration) / 60.0)
        audit["coverage_mask"] = observed
        remi_audits.append(audit)
        display_name = tracks["rftn20_ce"] if name == tracks["rftn20_rate"] else tracks["rftn50_ce"]
        remi_display_sources.append(native[display_name])
    if not remi_inputs:
        raise RuntimeError("Remifentanil rate track disappeared after metadata selection")
    remi_input = np.nansum(np.vstack(remi_inputs), axis=0)
    all_missing = np.all(~np.isfinite(np.vstack(remi_inputs)), axis=0)
    remi_input[all_missing] = np.nan
    remi_rate_observed = np.logical_and.reduce(
        [np.asarray(item.pop("coverage_mask"), dtype=bool) for item in remi_audits]
    )
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
    art_10s_centers, art_10s_medians = _art_ten_second_series(
        art_times,
        art_values,
        start=maintenance_start,
        stop=maintenance_stop,
        cadence=art_cadence,
        spec=spec,
    )
    low_map_events = _low_map_events(art_10s_centers, art_10s_medians)
    sqi_times, sqi_values = native[tracks["sqi"]]
    sr_times, sr_values = native[tracks["sr"]]
    etco2_name = max(
        tracks["etco2"],
        key=lambda name: int(
            np.sum((native[name][0] >= maintenance_start) & (native[name][0] <= maintenance_stop))
        ),
    )
    etco2_times, etco2_values = native[etco2_name]
    bt_times, bt_values = native[tracks["body_temperature"]]
    half_samples = round(float(spec["window_seconds"]) * float(spec["sample_hz"]) / 2.0)
    rows: list[dict[str, Any]] = []
    for center in centers:
        midpoint = round(float(center) * float(spec["sample_hz"]))
        raw_window = raw[max(midpoint - half_samples, 0) : min(midpoint + half_samples, len(raw))]
        eeg_qc = _raw_eeg_qc(raw_window, spec)
        art_qc = _art_window_qc(
            art_times,
            art_values,
            center=float(center),
            seconds=float(spec["window_seconds"]),
            cadence=art_cadence,
            spec=spec,
        )
        sqi_coverage, sqi_fraction = _condition_fraction_native(
            sqi_times,
            sqi_values,
            center=float(center),
            seconds=float(spec["window_seconds"]),
            predicate=lambda values: values >= float(spec["sqi_threshold"]),
            maintenance_start=maintenance_start,
            maintenance_stop=maintenance_stop,
        )
        sr_coverage, sr_fraction = _condition_fraction_native(
            sr_times,
            sr_values,
            center=float(center),
            seconds=float(spec["window_seconds"]),
            predicate=lambda values: values < float(spec["sr_threshold"]),
            maintenance_start=maintenance_start,
            maintenance_stop=maintenance_stop,
        )
        sr_positive_fraction, sr_positive_median = _positive_summary_native(
            sr_times,
            sr_values,
            center=float(center),
            seconds=float(spec["window_seconds"]),
            positive_threshold=float(spec["sr_presence_positive_threshold"]),
        )
        endpoint = endpoints[round(center)]
        spectral_qc = int(endpoint["welch_valid_segments"]) >= int(
            spec["welch_minimum_valid_segments"]
        )
        signal_eligible = (
            bool(eeg_qc["eeg_raw_qc"])
            and spectral_qc
            and sqi_coverage >= float(spec["sqi_window_fraction_min"])
            and sqi_fraction >= float(spec["sqi_window_fraction_min"])
        )
        sr_observed = sr_coverage >= float(spec["sr_window_fraction_min"])
        primary_ra = (
            signal_eligible and sr_observed and sr_fraction >= float(spec["sr_window_fraction_min"])
        )
        suppression_present = bool(
            signal_eligible
            and sr_observed
            and sr_positive_fraction >= float(spec["sr_presence_window_fraction_min"])
        )
        broad = bool(volatile_audit["broad_no_positive_volatile"])
        ppf_value = _series_window_median(ppf_ce, float(center), float(spec["window_seconds"]))
        remi_value = _series_window_median(remi_ce, float(center), float(spec["window_seconds"]))
        ppf_display = _native_median(
            *native[tracks["ppf_ce"]], float(center), float(spec["window_seconds"])
        )
        remi_display_values = [
            _native_median(times, values, float(center), float(spec["window_seconds"]))
            for times, values in remi_display_sources
        ]
        remi_display = (
            float(np.nansum(remi_display_values))
            if any(np.isfinite(value) for value in remi_display_values)
            else float("nan")
        )
        ppf_observed_fraction = _series_window_fraction(
            ppf_rate_observed, float(center), float(spec["window_seconds"])
        )
        remi_observed_fraction = _series_window_fraction(
            remi_rate_observed, float(center), float(spec["window_seconds"])
        )
        etco2_value = _native_median(
            etco2_times, etco2_values, float(center), float(spec["window_seconds"])
        )
        temperature_value = _native_median(
            bt_times, bt_values, float(center), float(spec["window_seconds"])
        )
        complete = bool(
            primary_ra
            and art_qc["art_qc"]
            and broad
            and np.isfinite(ppf_value)
            and np.isfinite(remi_value)
            and np.isfinite(etco2_value)
            and np.isfinite(temperature_value)
            and ppf_observed_fraction >= 0.8
            and remi_observed_fraction >= 0.8
        )
        suppression_complete = bool(
            signal_eligible
            and sr_observed
            and art_qc["art_qc"]
            and broad
            and np.isfinite(ppf_value)
            and np.isfinite(remi_value)
            and np.isfinite(etco2_value)
            and np.isfinite(temperature_value)
            and ppf_observed_fraction >= 0.8
            and remi_observed_fraction >= 0.8
        )
        rows.append(
            {
                "caseid": int(case["caseid"]),
                "subjectid": int(case["subjectid"]),
                "time_seconds": float(center),
                "minutes_since_opstart": (float(center) - float(case["opstart"])) / 60.0,
                "operation_progress": (float(center) - float(case["opstart"]))
                / max(float(case["opend"]) - float(case["opstart"]), 1.0),
                "age": float(case["age"]),
                "sex": str(case["sex"]),
                "asa": case.get("asa"),
                "department": case.get("department"),
                "optype": case.get("optype"),
                "emergency": case.get("emop"),
                "relative_alpha": endpoint["relative_alpha"],
                "logit_relative_alpha": endpoint["logit_relative_alpha"],
                "sef95": endpoint["sef95"],
                "alpha_peak_frequency": endpoint["alpha_peak_frequency"],
                "alpha_peak_is_edge": bool(endpoint["alpha_peak_is_edge"]),
                "alpha_power": endpoint["alpha_power"],
                "delta_power": endpoint["delta_power"],
                "log_adr": endpoint["log_adr"],
                "welch_valid_segments": int(endpoint["welch_valid_segments"]),
                "welch_total_segments": int(endpoint["welch_total_segments"]),
                **eeg_qc,
                "sqi_coverage": sqi_coverage,
                "sqi_ge80_fraction": sqi_fraction,
                "sr_coverage": sr_coverage,
                "sr_lt20_fraction": sr_fraction,
                "sr_positive_fraction": sr_positive_fraction,
                "sr_positive_median": sr_positive_median,
                "suppression_present": suppression_present,
                "high_suppression_burden_window": bool(
                    signal_eligible
                    and sr_observed
                    and sr_fraction < float(spec["sr_window_fraction_min"])
                ),
                "suppressed_window": bool(
                    signal_eligible
                    and sr_observed
                    and sr_fraction < float(spec["sr_window_fraction_min"])
                ),
                "eeg_signal_eligible": bool(signal_eligible),
                "primary_ra_eligible": bool(primary_ra),
                **art_qc,
                "propofol_ce_reconstructed": ppf_value,
                "remifentanil_ce_reconstructed": remi_value,
                "propofol_ce_device": ppf_display,
                "remifentanil_ce_device": remi_display,
                "propofol_rate_observed_fraction": ppf_observed_fraction,
                "remifentanil_rate_observed_fraction": remi_observed_fraction,
                "etco2": etco2_value,
                "body_temperature": temperature_value,
                "volatile_state": volatile_state,
                "broad_no_positive_volatile": broad,
                "strict_measured_zero_volatile": bool(
                    volatile_audit["strict_measured_zero_volatile"]
                ),
                "primary_complete_window": complete,
                "suppression_complete_window": suppression_complete,
                **_preceding_map_features(art_10s_centers, art_10s_medians, endpoint=float(center)),
            }
        )
    audit = {
        "caseid": int(case["caseid"]),
        "subjectid": int(case["subjectid"]),
        "scheduled_windows": len(centers),
        "art_native_cadence_seconds": art_cadence,
        "etco2_source": etco2_name,
        **volatile_audit,
        **{f"ppf_{key}": value for key, value in ppf_rate_audit.items()},
        "remi_min_rate_coverage": float(min(item["rate_coverage"] for item in remi_audits)),
        "remi_max_rate_gap_seconds": float(
            max(item["rate_max_gap_seconds"] for item in remi_audits)
        ),
        "low_map_event_count": len(low_map_events),
        "low_map_event_times_seconds": low_map_events,
        "art_10s_bins": len(art_10s_medians),
        "art_10s_valid_bins": int(np.isfinite(art_10s_medians).sum()),
    }
    frame = pd.DataFrame(rows)
    audit["primary_complete_windows"] = int(frame["primary_complete_window"].sum())
    audit["primary_complete_8plus"] = bool(
        frame["primary_complete_window"].sum() >= int(spec["minimum_valid_windows"])
    )
    return frame, audit


def _extract_one(
    cfg: dict[str, Any], case: dict[str, Any], shard_text: str, force: bool
) -> dict[str, Any]:
    shard = Path(shard_text)
    audit_path = shard.with_suffix(".audit.json")
    if shard.exists() and audit_path.exists() and not force:
        return json.loads(audit_path.read_text(encoding="utf-8"))
    try:
        frame, audit = extract_case(cfg, case)
        audit["status"] = "ok"
        shard.parent.mkdir(parents=True, exist_ok=True)
        temporary = shard.with_suffix(f".partial-{os.getpid()}.parquet")
        frame.to_parquet(temporary, index=False)
        temporary.replace(shard)
    except Exception as exc:  # noqa: BLE001 - every selected case is retained in the audit
        audit = {
            "caseid": int(case["caseid"]),
            "subjectid": int(case["subjectid"]),
            "status": "error",
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    audit_path.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )
    return audit
