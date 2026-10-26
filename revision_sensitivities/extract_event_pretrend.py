"""Re-extract the two EEG/MAP intervals immediately before each event start."""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

from vitaldb_eeg.config import load_config

PRETREND_INTERVALS = (
    ("pretrend_m2_m1", -120.0, -60.0),
    ("pretrend_m1_0", -60.0, 0.0),
)


def extract_one(
    config_path: str,
    case: dict[str, Any],
    shard_text: str,
    force: bool,
) -> dict[str, Any]:
    from vitaldb_eeg import maintenance_ra_events as events

    cfg = load_config(config_path)
    events.EVENT_INTERVALS = PRETREND_INTERVALS
    shard = Path(shard_text)
    interval_shard = shard.with_name(shard.stem + "_intervals.parquet")
    audit_path = shard.with_suffix(".audit.json")
    if shard.exists() and interval_shard.exists() and audit_path.exists() and not force:
        return json.loads(audit_path.read_text(encoding="utf-8"))
    try:
        subwindows, intervals, audit = events.extract_event_case(cfg, case)
        shard.parent.mkdir(parents=True, exist_ok=True)
        for frame, target in ((subwindows, shard), (intervals, interval_shard)):
            temporary = target.with_suffix(f".partial-{os.getpid()}.parquet")
            frame.to_parquet(temporary, index=False)
            temporary.replace(target)
    except Exception as error:  # noqa: BLE001 - every failure is retained for audit
        audit = {
            "caseid": int(case["caseid"]),
            "subjectid": int(case["subjectid"]),
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error),
        }
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-extract the two intervals immediately before sustained low-MAP event t=0."
    )
    parser.add_argument("--config", required=True, help="Analysis configuration file.")
    parser.add_argument(
        "--source-root", required=True, help="Directory containing the event case audit."
    )
    parser.add_argument(
        "--output", required=True, help="Directory for the extracted pre-event intervals."
    )
    parser.add_argument("--workers", type=int, default=8, help="Parallel case workers.")
    parser.add_argument("--force", action="store_true", help="Re-extract existing case shards.")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    return args


def main() -> None:
    args = parse_args()
    source_root = Path(args.source_root)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    cases = pd.read_parquet(source_root / "eligible_cases.parquet")
    case_audit = pd.read_csv(source_root / "case_audit.csv")
    ids = case_audit.loc[
        case_audit["status"].eq("ok")
        & pd.to_numeric(case_audit["low_map_event_count"], errors="coerce").fillna(0).gt(0),
        "caseid",
    ].astype(int)
    cases = cases.loc[cases["caseid"].astype(int).isin(ids)].copy()
    shards = output / "shards"
    jobs = [
        (
            str(Path(args.config).resolve()),
            row._asdict(),
            str(shards / f"case_{int(row.caseid):04d}_pretrend.parquet"),
            bool(args.force),
        )
        for row in cases.itertuples(index=False)
    ]
    audits: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        futures = {pool.submit(extract_one, *job): int(job[1]["caseid"]) for job in jobs}
        for future in tqdm(as_completed(futures), total=len(futures), desc="pre-event -2 to 0 min"):
            audits.append(future.result())
    audit_frame = pd.DataFrame(audits).sort_values(["subjectid", "caseid"])
    audit_frame.to_csv(output / "pretrend_case_audit.csv", index=False)
    good_ids = audit_frame.loc[audit_frame["status"].eq("ok"), "caseid"].astype(int)
    subwindow_paths = [shards / f"case_{caseid:04d}_pretrend.parquet" for caseid in good_ids]
    interval_paths = [
        shards / f"case_{caseid:04d}_pretrend_intervals.parquet" for caseid in good_ids
    ]
    subwindows = pd.concat([pd.read_parquet(path) for path in subwindow_paths], ignore_index=True)
    intervals = pd.concat([pd.read_parquet(path) for path in interval_paths], ignore_index=True)
    subwindows.to_parquet(output / "event_pretrend_subwindows.parquet", index=False)
    intervals.to_parquet(output / "event_pretrend_intervals.parquet", index=False)
    subwindows.to_csv(output / "event_pretrend_subwindows.csv", index=False)
    intervals.to_csv(output / "event_pretrend_intervals.csv", index=False)
    summary = {
        "definition": [list(item) for item in PRETREND_INTERVALS],
        "candidate_cases": len(cases),
        "cases_processed_ok": int(audit_frame["status"].eq("ok").sum()),
        "cases_failed": int(audit_frame["status"].ne("ok").sum()),
        "real_events": int(intervals.loc[intervals["event_type"].eq("real"), "event_id"].nunique()),
        "pseudo_events": int(
            intervals.loc[intervals["event_type"].eq("pseudo"), "event_id"].nunique()
        ),
    }
    (output / "pretrend_extraction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
