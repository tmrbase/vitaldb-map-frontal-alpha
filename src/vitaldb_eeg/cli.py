"""Command-line interface for cohort prescreening, window extraction, and events."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .config import load_config
from .maintenance_ra_events import prepare_maintenance_ra_full_events
from .maintenance_ra_extraction import prepare_maintenance_ra_prescreen
from .maintenance_ra_full import prepare_maintenance_ra_full


def _common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/analysis_public.yaml"),
        help="Analysis configuration file.",
    )
    parser.add_argument(
        "--workers", type=int, help="Parallel case workers; defaults to the configuration."
    )
    parser.add_argument("--force", action="store_true", help="Rebuild verified cached products.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare the datasets used in the VitalDB MAP-alpha analysis."
    )
    commands = parser.add_subparsers(dest="command", required=True)
    descriptions = {
        "prepare-volatile-prescreen": (
            "Identify surgeries without recorded positive expired volatile concentration."
        ),
        "prepare-maintenance-windows": (
            "Construct quality-controlled maintenance-phase EEG/MAP windows."
        ),
        "prepare-low-map-events": "Construct sustained low-MAP event and reference intervals.",
    }
    for name, description in descriptions.items():
        _common_options(commands.add_parser(name, help=description, description=description))
    return parser


def dispatch(args: argparse.Namespace, cfg: dict[str, Any]) -> dict[str, Any]:
    if args.command == "prepare-volatile-prescreen":
        return prepare_maintenance_ra_prescreen(cfg, workers=args.workers, force=args.force)
    if args.command == "prepare-maintenance-windows":
        return prepare_maintenance_ra_full(cfg, workers=args.workers, force=args.force)
    return prepare_maintenance_ra_full_events(cfg, workers=args.workers, force=args.force)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be at least 1")
    result = dispatch(args, load_config(args.config))
    print(json.dumps(result, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
