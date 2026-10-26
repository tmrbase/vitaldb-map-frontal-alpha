import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from vitaldb_eeg.maintenance_ra_events import (
    _aggregate_interval,
    _subwindow_endpoint,
    actual_log_nonalpha,
)
from vitaldb_eeg.maintenance_ra_extraction import _low_map_events
from vitaldb_eeg.signal import build_welch_psd_grid, filter_valid_blocks

ROOT = Path(__file__).resolve().parents[1]
TRACK_CATALOG_SHA256 = "29374aae9449c18f3d1ea360dbf7712614a70dcfe63eedf639bab27c8a1556c9"


def test_renv_lock_is_restorable_and_binary_compatible() -> None:
    project = ROOT / "r_analysis" / "maintenance_ra"
    lock = json.loads((project / "renv.lock").read_text(encoding="utf-8-sig"))
    assert lock["R"]["Version"] == "4.4.1"
    expected = {
        "data.table": "1.18.2.1",
        "nlme": "3.1-169",
        "jsonlite": "2.0.0",
        "clubSandwich": "0.6.2",
        "sandwich": "3.1-1",
        "lmtest": "0.9-40",
        "ggplot2": "4.0.3",
        "MASS": "7.3-65",
        "renv": "1.2.2",
    }
    assert {name: lock["Packages"][name]["Version"] for name in expected} == expected
    for record in lock["Packages"].values():
        for field in ("Depends", "Imports", "Suggests", "LinkingTo", "Enhances"):
            if field in record:
                assert isinstance(record[field], str), f"{record['Package']} {field} is not scalar"
    profile = (project / ".Rprofile").read_text(encoding="utf-8")
    assert "RENV_CONFIG_CACHE_ENABLED" in profile


def test_pinned_vitaldb_v100_track_catalog() -> None:
    catalog = ROOT / "source_catalog" / "vitaldb_v1.0.0_track_index.tsv.gz"
    assert hashlib.sha256(catalog.read_bytes()).hexdigest() == TRACK_CATALOG_SHA256
    with gzip.open(catalog, "rt", encoding="utf-8") as handle:
        rows = handle.readlines()
    assert len(rows) == 6_388
    assert rows[0].split("\t", 1)[0].endswith("0001.vital")
    assert rows[-1].split("\t", 1)[0].endswith("6388.vital")


def load_public_config() -> dict:
    return yaml.safe_load((ROOT / "configs" / "analysis_public.yaml").read_text(encoding="utf-8"))


def test_filter_guard_and_subwindow_future_leakage() -> None:
    sample_hz = 128
    time = np.arange(45 * sample_hz) / sample_hz
    base = np.sin(2 * np.pi * 10 * time)
    filtered = filter_valid_blocks(
        base[: 30 * sample_hz],
        sample_hz=sample_hz,
        low_hz=0.5,
        high_hz=30.0,
        order=4,
        guard_seconds=1.0,
    )
    assert np.isnan(filtered.values[:sample_hz]).all()
    assert np.isnan(filtered.values[-sample_hz:]).all()

    changed_future = base.copy()
    changed_future[30 * sample_hz :] = 1_000.0
    cfg = load_public_config()
    first = _subwindow_endpoint(base, left=0.0, right=30.0, cfg=cfg)
    second = _subwindow_endpoint(changed_future, left=0.0, right=30.0, cfg=cfg)
    for key in ("relative_alpha", "logit_relative_alpha", "log_alpha_power"):
        assert np.isclose(first[key], second[key], equal_nan=True)


def test_welch_segment_count_for_thirty_seconds() -> None:
    sample_hz = 128
    values = np.ones(30 * sample_hz)
    grid = build_welch_psd_grid(
        values,
        sample_hz=sample_hz,
        nperseg=256,
        noverlap=128,
        detrend="constant",
    )
    expected = (len(values) - 256) // (256 - 128) + 1
    assert expected == 29
    assert len(grid.cumulative_valid) - 1 == expected
    assert grid.cumulative_valid[-1] == expected


def test_event_time_is_start_of_first_ten_second_low_map_bin() -> None:
    centers = np.arange(5.0, 405.0, 10.0)
    medians = np.full_like(centers, 70.0)
    medians[(centers >= 305.0) & (centers < 365.0)] = 64.0
    events = _low_map_events(centers, medians, maximum_events=2)
    assert events == [300.0]


def event_rows(*, total: int = 4, eligible_count: int = 2) -> pd.DataFrame:
    eligible = [True] * eligible_count + [False] * (total - eligible_count)
    return pd.DataFrame(
        {
            "primary_ra_eligible": eligible,
            "eeg_signal_eligible": eligible,
            "relative_alpha": [0.20 + 0.02 * index for index in range(eligible_count)]
            + [np.nan] * (total - eligible_count),
            "log_alpha_power": [1.0 + 0.1 * index for index in range(eligible_count)]
            + [np.nan] * (total - eligible_count),
            "log_nonalpha_power": [2.0 + 0.1 * index for index in range(eligible_count)]
            + [np.nan] * (total - eligible_count),
            "log_delta_power": [3.0 + 0.1 * index for index in range(eligible_count)]
            + [np.nan] * (total - eligible_count),
            "log_adr": [-2.0] * eligible_count + [np.nan] * (total - eligible_count),
            "propofol_ce": [2.0] * total,
            "remifentanil_ce": [3.0] * total,
            "etco2": [35.0] * total,
            "body_temperature": [36.5] * total,
            "art_map": [64.0] * total,
        }
    )


def aggregate_at_fraction(fraction: float, *, total: int = 4, eligible_count: int = 2) -> dict:
    spec = {
        "event_interval_valid_fraction_min": fraction,
        "sr_presence_positive_threshold": 0.0,
        "sr_window_fraction_min": 0.8,
        "sr_presence_window_fraction_min": 0.2,
    }
    return _aggregate_interval(
        event_rows(total=total, eligible_count=eligible_count),
        event_id="1:real:1",
        event_type="real",
        event_time=300.0,
        interval_name="response_0_1",
        relative_start=0.0,
        relative_stop=60.0,
        sr_times=np.array([300.0, 360.0]),
        sr_values=np.array([0.0, 0.0]),
        maintenance_start=0.0,
        maintenance_stop=1_000.0,
        spec=spec,
    )


def test_event_interval_completeness_primary_and_sensitivity() -> None:
    primary = aggregate_at_fraction(0.50)
    sensitivity = aggregate_at_fraction(0.80)
    assert primary["required_subwindows"] == 2
    assert primary["ra_interval_eligible"] is True
    assert sensitivity["required_subwindows"] == 4
    assert sensitivity["ra_interval_eligible"] is False


def test_first_event_minute_two_subwindow_completeness() -> None:
    primary = aggregate_at_fraction(0.50, total=2, eligible_count=1)
    sensitivity = aggregate_at_fraction(0.80, total=2, eligible_count=1)
    assert primary["planned_subwindows"] == 2
    assert primary["required_subwindows"] == 1
    assert primary["ra_interval_eligible"] is True
    assert sensitivity["required_subwindows"] == 2
    assert sensitivity["ra_interval_eligible"] is False


def test_nonalpha_is_computed_before_interval_aggregation() -> None:
    log_alpha = np.array([1.0, 2.0])
    relative_alpha = np.array([0.20, 0.60])
    logit_ra = np.log(relative_alpha / (1.0 - relative_alpha))
    actual = np.array(
        [actual_log_nonalpha(alpha, logit) for alpha, logit in zip(log_alpha, logit_ra)]
    )
    aggregated_after_calculation = actual.mean()
    incorrectly_derived_after_aggregation = log_alpha.mean() - np.log(
        relative_alpha.mean() / (1.0 - relative_alpha.mean())
    )
    assert not np.isclose(aggregated_after_calculation, incorrectly_derived_after_aggregation)
