"""Run extraction of additional physiological covariates for sensitivity models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from vitaldb_eeg.config import load_config
from vitaldb_eeg.maintenance_ra_additional_covariates import extract_additional_covariates


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract N2O, heart rate, SpO2, cardiac output and cardiac index."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/analysis_public.yaml"),
        help="Analysis configuration file.",
    )
    parser.add_argument(
        "--workers", type=int, help="Parallel case workers; defaults to the configuration."
    )
    parser.add_argument("--force", action="store_true", help="Rebuild existing case-level shards.")
    args = parser.parse_args()
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be at least 1")
    cfg = load_config(args.config)
    summary = extract_additional_covariates(cfg, workers=args.workers, force=args.force)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
