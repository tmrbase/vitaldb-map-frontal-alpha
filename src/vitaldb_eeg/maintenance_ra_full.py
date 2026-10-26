"""Orchestrate maintenance-window and sustained-low-MAP-event extraction."""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

from .config import ensure_output_directories
from .maintenance_ra_extraction import (
    _extract_one,
    _metadata_eligible_cases,
)


def _output_root(cfg: dict[str, Any]) -> Path:
    return ensure_output_directories(cfg).output_root / "maintenance_ra_full_v1"


def build_full_manifest(cfg: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build the full result-independent cohort from the volatile prescreen."""

    paths = ensure_output_directories(cfg)
    prescreen_root = paths.output_root / "maintenance_ra_prescreen_v1"
    prescreen_path = prescreen_root / "volatile_prescreen.csv"
    prescreen_manifest_path = prescreen_root / "volatile_prescreen_manifest.json"
    if not prescreen_path.exists() or not prescreen_manifest_path.exists():
        raise FileNotFoundError(
            "The native-cadence volatile prescreen must be completed before the full run"
        )

    metadata = _metadata_eligible_cases(cfg)
    prescreen = pd.read_csv(prescreen_path)
    expected = set(metadata["caseid"].astype(int))
    observed = set(prescreen["caseid"].astype(int))
    if expected != observed:
        raise RuntimeError(
            "The cached volatile prescreen does not match the current metadata frame: "
            f"missing={len(expected - observed)}, unexpected={len(observed - expected)}"
        )
    broad_caseids = set(
        prescreen.loc[prescreen["broad_no_positive_volatile"].astype(bool), "caseid"].astype(int)
    )
    cases = metadata.loc[metadata["caseid"].astype(int).isin(broad_caseids)].copy()

    cases = cases.sort_values(["subjectid", "caseid"]).reset_index(drop=True)
    output = _output_root(cfg)
    output.mkdir(parents=True, exist_ok=True)
    cases.to_parquet(output / "eligible_cases.parquet", index=False)
    cases.to_csv(output / "eligible_cases.csv", index=False)
    manifest = {
        "status": "full broad operational cohort",
        "selection": "all metadata-eligible cases with no positive expired volatile signal",
        "subjects": int(cases["subjectid"].nunique()),
        "surgeries": int(cases["caseid"].nunique()),
        "prescreen_manifest": str(prescreen_manifest_path),
        "raw_outcome_used_for_selection": False,
        "raw_map_used_for_selection": False,
    }
    (output / "full_cohort_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return cases, manifest


def prepare_maintenance_ra_full(
    cfg: dict[str, Any],
    *,
    workers: int | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Extract the window-level feature set for all eligible cases."""

    cases, manifest = build_full_manifest(cfg)
    output = _output_root(cfg)
    shards = output / "case_shards"
    jobs = [
        (row._asdict(), shards / f"case_{int(row.caseid):04d}.parquet")
        for row in cases.itertuples(index=False)
    ]
    audits: list[dict[str, Any]] = []
    worker_count = max(1, int(workers or cfg["runtime"]["workers"]))
    if worker_count == 1:
        for case, shard in tqdm(jobs, desc="maintenance RA full"):
            audits.append(_extract_one(cfg, case, str(shard), force))
    else:
        with ProcessPoolExecutor(max_workers=min(worker_count, len(jobs))) as pool:
            futures = {
                pool.submit(_extract_one, cfg, case, str(shard), force): int(case["caseid"])
                for case, shard in jobs
            }
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="maintenance RA full"
            ):
                audits.append(future.result())

    audit_frame = pd.DataFrame(audits).sort_values(["subjectid", "caseid"])
    audit_frame.to_csv(output / "case_audit.csv", index=False)
    good_caseids = audit_frame.loc[audit_frame["status"].eq("ok"), "caseid"].astype(int)
    paths = [shards / f"case_{caseid:04d}.parquet" for caseid in good_caseids]
    frames = [pd.read_parquet(path) for path in paths if path.exists()]
    windows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not windows.empty:
        windows = windows.sort_values(["subjectid", "caseid", "time_seconds"])
        windows.to_parquet(output / "full_windows.parquet", index=False)
        windows.to_csv(output / "full_windows.csv", index=False)

    minimum = int(cfg["analysis"]["minimum_valid_windows"])
    summary = {
        **manifest,
        "cases_processed_ok": int(audit_frame["status"].eq("ok").sum()),
        "cases_failed": int(audit_frame["status"].ne("ok").sum()),
        "scheduled_windows": len(windows),
        "primary_complete_windows": (
            int(windows["primary_complete_window"].sum()) if not windows.empty else 0
        ),
        "surgeries_with_minimum_complete_windows": (
            int(windows.groupby("caseid")["primary_complete_window"].sum().ge(minimum).sum())
            if not windows.empty
            else 0
        ),
        "strict_measured_zero_surgeries": (
            int(windows.loc[windows["strict_measured_zero_volatile"], "caseid"].nunique())
            if not windows.empty
            else 0
        ),
    }
    (output / "full_extraction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary
