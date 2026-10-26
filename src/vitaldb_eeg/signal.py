"""Apply EEG filtering, Welch spectra, and signal-quality calculations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class FilterResult:
    values: np.ndarray
    valid_blocks: int
    filtered_samples: int
    guard_samples_removed: int


@dataclass(frozen=True)
class WelchPsdGrid:
    frequencies: np.ndarray
    cumulative_psd: np.ndarray
    cumulative_valid: np.ndarray
    sample_hz: float
    segment_seconds: float
    step_seconds: float


def contiguous_true_blocks(mask: np.ndarray) -> list[tuple[int, int]]:
    values = np.asarray(mask, dtype=bool)
    if values.size == 0:
        return []
    padded = np.pad(values.astype(np.int8), (1, 1))
    transitions = np.diff(padded)
    starts = np.flatnonzero(transitions == 1)
    stops = np.flatnonzero(transitions == -1)
    return [(int(start), int(stop)) for start, stop in zip(starts, stops, strict=True)]


def design_bandpass_sos(
    sample_hz: float,
    low_hz: float,
    high_hz: float,
    order: int,
) -> np.ndarray:
    nyquist = float(sample_hz) / 2.0
    if not 0 < low_hz < high_hz < nyquist:
        raise ValueError(f"Invalid band [{low_hz}, {high_hz}] for fs={sample_hz}")
    return signal.butter(
        int(order),
        [float(low_hz), float(high_hz)],
        btype="bandpass",
        fs=float(sample_hz),
        output="sos",
    )


def filter_valid_blocks(
    eeg: np.ndarray,
    *,
    sample_hz: float,
    low_hz: float,
    high_hz: float,
    order: int,
    guard_seconds: float,
) -> FilterResult:
    values = np.asarray(eeg, dtype=float).reshape(-1)
    output = np.full(values.shape, np.nan, dtype=float)
    sos = design_bandpass_sos(sample_hz, low_hz, high_hz, order)
    guard = round(float(guard_seconds) * float(sample_hz))
    min_block = max(3 * (2 * sos.shape[0] + 1), 2 * guard + 3)
    block_count = 0
    filtered_samples = 0
    guard_removed = 0
    for start, stop in contiguous_true_blocks(np.isfinite(values)):
        if stop - start < min_block:
            continue
        filtered = signal.sosfiltfilt(sos, values[start:stop])
        keep_start = start + guard
        keep_stop = stop - guard
        if keep_stop <= keep_start:
            continue
        output[keep_start:keep_stop] = filtered[guard : len(filtered) - guard]
        block_count += 1
        filtered_samples += keep_stop - keep_start
        guard_removed += 2 * guard
    return FilterResult(
        values=output,
        valid_blocks=block_count,
        filtered_samples=filtered_samples,
        guard_samples_removed=guard_removed,
    )


def sef95_from_psd(frequencies: np.ndarray, psd: np.ndarray) -> float:
    freq = np.asarray(frequencies, dtype=float)
    power = np.asarray(psd, dtype=float)
    valid = np.isfinite(freq) & np.isfinite(power) & (power >= 0)
    freq = freq[valid]
    power = power[valid]
    if freq.size < 2 or float(np.sum(power)) <= 0:
        return float("nan")
    cumulative = np.cumsum(power)
    target = 0.95 * cumulative[-1]
    index = int(np.searchsorted(cumulative, target, side="left"))
    if index <= 0:
        return float(freq[0])
    if index >= len(freq):
        return float(freq[-1])
    lower_power = cumulative[index - 1]
    upper_power = cumulative[index]
    if upper_power <= lower_power:
        return float(freq[index])
    fraction = (target - lower_power) / (upper_power - lower_power)
    return float(freq[index - 1] + fraction * (freq[index] - freq[index - 1]))


def build_welch_psd_grid(
    filtered_eeg: np.ndarray,
    *,
    sample_hz: float,
    nperseg: int,
    noverlap: int,
    detrend: str,
) -> WelchPsdGrid:
    """Compute aligned periodograms once for exact rolling Welch means."""
    values = np.asarray(filtered_eeg, dtype=float).reshape(-1)
    step = int(nperseg - noverlap)
    if step <= 0:
        raise ValueError("Welch step must be positive.")
    frequencies = np.fft.rfftfreq(nperseg, d=1.0 / sample_hz)
    if values.size < nperseg:
        return WelchPsdGrid(
            frequencies=frequencies,
            cumulative_psd=np.zeros((1, len(frequencies)), dtype=float),
            cumulative_valid=np.zeros(1, dtype=np.int64),
            sample_hz=float(sample_hz),
            segment_seconds=float(nperseg / sample_hz),
            step_seconds=float(step / sample_hz),
        )
    windows = np.lib.stride_tricks.sliding_window_view(values, nperseg)[::step]
    valid = np.isfinite(windows).all(axis=1)
    safe = np.where(valid[:, None], windows, 0.0)
    frequencies, periodograms = signal.periodogram(
        safe,
        fs=float(sample_hz),
        window="hann",
        detrend=detrend,
        scaling="density",
        axis=-1,
    )
    periodograms[~valid] = 0.0
    return WelchPsdGrid(
        frequencies=np.asarray(frequencies, dtype=float),
        cumulative_psd=np.vstack(
            [
                np.zeros((1, periodograms.shape[1]), dtype=float),
                np.cumsum(periodograms, axis=0, dtype=float),
            ]
        ),
        cumulative_valid=np.r_[0, np.cumsum(valid.astype(np.int64))],
        sample_hz=float(sample_hz),
        segment_seconds=float(nperseg / sample_hz),
        step_seconds=float(step / sample_hz),
    )


def spectral_endpoints_at_centers_from_grid(
    grid: WelchPsdGrid,
    centers_seconds: Iterable[int | float],
    *,
    window_seconds: int,
    sef_band_hz: tuple[float, float] = (0.5, 30.0),
    alpha_band_hz: tuple[float, float] = (8.0, 13.0),
    delta_band_hz: tuple[float, float] = (1.0, 4.0),
    relative_total_band_hz: tuple[float, float] = (1.0, 20.0),
    minimum_valid_segments: int | None = None,
) -> dict[int, dict[str, float]]:
    """Return transparent short-window spectral endpoints from one PSD grid.

    Power ratios are computed from the same mean Welch PSD used for SEF95 so
    endpoint comparisons do not differ because of window placement or PSD
    estimation. Natural-log transforms are returned only for positive powers.
    """
    if not np.isclose(grid.step_seconds, 1.0):
        raise ValueError("Current center mapping requires a 1-second Welch step.")
    centers = np.asarray([round(float(value)) for value in centers_seconds], dtype=int)
    empty = {
        "sef95": float("nan"),
        "alpha_peak_frequency": float("nan"),
        "alpha_peak_is_edge": False,
        "alpha_power": float("nan"),
        "delta_power": float("nan"),
        "relative_alpha": float("nan"),
        "logit_relative_alpha": float("nan"),
        "log_adr": float("nan"),
        "log_alpha_power": float("nan"),
        "log_delta_power": float("nan"),
        "welch_valid_segments": 0,
        "welch_total_segments": 0,
    }
    output = {int(center): dict(empty) for center in centers}
    if centers.size == 0:
        return output
    segment_count = round((float(window_seconds) - grid.segment_seconds) / grid.step_seconds) + 1
    starts = centers - int(window_seconds) // 2
    stops = starts + segment_count
    available = len(grid.cumulative_valid) - 1
    selected = np.flatnonzero((starts >= 0) & (stops <= available))
    if selected.size == 0:
        return output
    valid_counts = grid.cumulative_valid[stops[selected]] - grid.cumulative_valid[starts[selected]]
    required_segments = (
        int(segment_count) if minimum_valid_segments is None else int(minimum_valid_segments)
    )
    if required_segments < 1 or required_segments > segment_count:
        raise ValueError(
            f"minimum_valid_segments must be in [1, {segment_count}], got {required_segments}."
        )
    selected = selected[valid_counts >= required_segments]
    if selected.size == 0:
        return output
    valid_counts = grid.cumulative_valid[stops[selected]] - grid.cumulative_valid[starts[selected]]
    mean_psd = (
        grid.cumulative_psd[stops[selected]] - grid.cumulative_psd[starts[selected]]
    ) / valid_counts[:, None]

    def band_power(psd: np.ndarray, band_hz: tuple[float, float]) -> float:
        mask = (grid.frequencies >= band_hz[0]) & (grid.frequencies <= band_hz[1])
        if int(mask.sum()) < 2:
            return float("nan")
        return float(np.trapezoid(psd[mask], grid.frequencies[mask]))

    def alpha_peak(psd: np.ndarray) -> tuple[float, bool]:
        """Return the raw-PSD alpha peak used by Manquat et al.

        The primary definition is deliberately simple and literature-comparable:
        the frequency with maximum power from 8 through 14 Hz, inclusive.  A
        maximum at either band edge is retained for the primary analysis but is
        flagged for a prespecified interior-peak sensitivity analysis.
        """
        mask = (grid.frequencies >= 8.0) & (grid.frequencies <= 14.0)
        frequencies = grid.frequencies[mask]
        powers = np.asarray(psd[mask], dtype=float)
        valid = np.isfinite(frequencies) & np.isfinite(powers) & (powers >= 0)
        if not bool(valid.any()):
            return float("nan"), False
        frequencies = frequencies[valid]
        powers = powers[valid]
        index = int(np.argmax(powers))
        peak = float(frequencies[index])
        edge = bool(np.isclose(peak, 8.0, atol=1e-12) or np.isclose(peak, 14.0, atol=1e-12))
        return peak, edge

    sef_mask = (grid.frequencies >= sef_band_hz[0]) & (grid.frequencies <= sef_band_hz[1])
    for source_row, center_index in enumerate(selected):
        psd = mean_psd[source_row]
        alpha = band_power(psd, alpha_band_hz)
        delta = band_power(psd, delta_band_hz)
        total = band_power(psd, relative_total_band_hz)
        alpha_peak_frequency, alpha_peak_is_edge = alpha_peak(psd)
        relative_alpha = alpha / total if alpha > 0 and total > 0 else float("nan")
        if np.isfinite(relative_alpha) and 0 < relative_alpha < 1:
            logit_relative_alpha = float(np.log(relative_alpha / (1.0 - relative_alpha)))
        else:
            logit_relative_alpha = float("nan")
        output[int(centers[center_index])] = {
            "sef95": sef95_from_psd(grid.frequencies[sef_mask], psd[sef_mask]),
            "alpha_peak_frequency": alpha_peak_frequency,
            "alpha_peak_is_edge": alpha_peak_is_edge,
            "alpha_power": alpha,
            "delta_power": delta,
            "relative_alpha": float(relative_alpha),
            "logit_relative_alpha": logit_relative_alpha,
            "log_adr": float(np.log(alpha) - np.log(delta))
            if alpha > 0 and delta > 0
            else float("nan"),
            "log_alpha_power": float(np.log(alpha)) if alpha > 0 else float("nan"),
            "log_delta_power": float(np.log(delta)) if delta > 0 else float("nan"),
            "welch_valid_segments": int(valid_counts[source_row]),
            "welch_total_segments": int(segment_count),
        }
    return output
