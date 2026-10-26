"""Rebuild the MAP–RA analysis from VitalDB source records to reported results."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def quoted(command: list[str]) -> str:
    return subprocess.list2cmdline(command)


def resolve_rscript(requested: str) -> str:
    located = shutil.which(requested)
    if located:
        return located
    if os.name == "nt" and requested.lower() in {"rscript", "rscript.exe"}:
        program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
        candidates = sorted((program_files / "R").glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    return requested


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Download VitalDB inputs and regenerate analysis datasets and results."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--execute", action="store_true", help="Execute the displayed commands.")
    mode.add_argument(
        "--dry-run", action="store_true", help="Display commands without executing them."
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Parallel download, extraction, and bootstrap workers."
    )
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or full path.")
    parser.add_argument("--force", action="store_true", help="Rebuild verified cached products.")
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download even source records that already pass v1.0.0 checksum validation.",
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Skip the PhysioNet v1.0.0 download stage and use the existing local data/vitaldb directory.",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")

    python = sys.executable
    config = root / "configs" / "analysis_public.yaml"
    full = root / "outputs" / "maintenance_ra_full_v1"
    scripts = root / "scripts"
    revision = root / "revision_sensitivities"
    force = ["--force"] if args.force else []
    force_download = ["--force"] if args.force_download else []
    commands: list[tuple[list[str], Path]] = []

    if not args.skip_fetch:
        commands.append(
            (
                [
                    python,
                    str(scripts / "fetch_vitaldb_inputs.py"),
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force_download,
                ],
                root,
            )
        )
    commands.extend(
        [
            (
                [
                    python,
                    "-m",
                    "vitaldb_eeg",
                    "prepare-volatile-prescreen",
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    "-m",
                    "vitaldb_eeg",
                    "prepare-maintenance-windows",
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    "-m",
                    "vitaldb_eeg",
                    "prepare-low-map-events",
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    str(revision / "extract_event_pretrend.py"),
                    "--config",
                    str(config),
                    "--source-root",
                    str(full),
                    "--output",
                    str(full),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    str(scripts / "extract_additional_covariates.py"),
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    str(scripts / "extract_hrv_covariates.py"),
                    "--config",
                    str(config),
                    "--workers",
                    str(args.workers),
                    *force,
                ],
                root,
            ),
            (
                [
                    python,
                    str(scripts / "run_analysis.py"),
                    "--execute",
                    "--workers",
                    str(args.workers),
                    "--rscript",
                    resolve_rscript(args.rscript),
                ],
                root,
            ),
            (
                [
                    python,
                    str(scripts / "validate_release.py"),
                    "--generated-results",
                ],
                root,
            ),
        ]
    )

    execute = args.execute and not args.dry_run
    for index, (command, cwd) in enumerate(commands, start=1):
        print(f"[{index}/{len(commands)}] cwd={cwd}")
        print(quoted(command), flush=True)
        if execute:
            subprocess.run(command, cwd=cwd, check=True)
    if not execute:
        print("Dry run only. Re-run with --execute to perform the source-to-result rebuild.")
    else:
        print("Source-to-result workflow completed successfully.")


if __name__ == "__main__":
    main()
