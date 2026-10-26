"""Run the reported analyses and build the final figures in dependency order."""

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
        candidates = sorted(
            (program_files / "R").glob("R-*/bin/Rscript.exe"),
            reverse=True,
        )
        if candidates:
            return str(candidates[0])
    return requested


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Run every statistical analysis and reporting step used in the manuscript."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--execute", action="store_true", help="Execute the displayed commands.")
    mode.add_argument(
        "--dry-run", action="store_true", help="Display commands without executing them."
    )
    parser.add_argument(
        "--workers", type=int, default=6, help="Parallel workers for the bootstrap."
    )
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or full path.")
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=root / "outputs" / "release_assets",
        help="Directory for regenerated figures and aggregate source tables.",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")

    rscript = resolve_rscript(args.rscript)
    r_work = root / "r_analysis" / "maintenance_ra"
    revision = root / "revision_sensitivities"
    full = root / "outputs" / "maintenance_ra_full_v1"
    analysis = root / "outputs" / "maintenance_ra_r_reanalysis_v1"
    model_objects = analysis / "hinge65_model_objects.rds"

    commands: list[tuple[list[str], Path]] = [
        ([rscript, "01_main_model.R"], r_work),
        ([rscript, "18_random_slope_primary.R"], r_work),
        (
            [rscript, "19_random_slope_patient_bootstrap.R", "2000", str(args.workers), "100"],
            r_work,
        ),
        ([rscript, "24_build_primary_curve.R"], r_work),
        (
            [sys.executable, str(root / "scripts" / "generate_patient_influence.py")],
            root,
        ),
        ([rscript, "03_model_sensitivities.R"], r_work),
        ([rscript, "04_secondary_eeg.R"], r_work),
        ([rscript, "20_alpha_nonalpha_decomposition.R"], r_work),
        ([rscript, "05_event_analysis.R", "0.50"], r_work),
        ([rscript, "05_event_analysis.R", "0.80"], r_work),
        ([rscript, "06_event_spectral.R"], r_work),
        ([rscript, "07_additional_covariates.R"], r_work),
        ([rscript, "08_hrv_sensitivity.R"], r_work),
        ([rscript, "09_event_map_explanation.R"], r_work),
        ([rscript, "10_table1_missingness.R"], r_work),
        ([rscript, "11_event_hrv_sensitivity.R"], r_work),
        ([rscript, "21_event_confounding_sensitivity.R"], r_work),
        ([rscript, "22_selection_and_pk_audit.R"], r_work),
        ([rscript, "23_residual_diagnostics.R"], r_work),
        ([rscript, "26_window_observability_audit.R"], r_work),
        (
            [
                sys.executable,
                str(revision / "common_support_audit.py"),
                "--windows",
                str(full / "full_windows.csv"),
                "--output-dir",
                str(analysis),
            ],
            root,
        ),
        (
            [rscript, str(revision / "overlap_restricted_contrasts.R"), str(root), str(analysis)],
            r_work,
        ),
        ([rscript, str(revision / "run_hinge65_analysis.R"), str(root), str(analysis)], r_work),
        (
            [
                rscript,
                str(revision / "check_hinge65_diagnostics.R"),
                str(model_objects),
                str(analysis / "hinge65_diagnostics.csv"),
            ],
            r_work,
        ),
        (
            [
                rscript,
                str(revision / "hinge65_cluster_robust.R"),
                str(model_objects),
                str(analysis / "hinge65_patient_cluster_CR1S.csv"),
            ],
            r_work,
        ),
        (
            [
                rscript,
                str(revision / "primary_map_contrast_patient_CR1S.R"),
                str(root),
                str(analysis / "primary_map_contrast_patient_CR1S.csv"),
            ],
            r_work,
        ),
        (
            [
                rscript,
                str(revision / "analyse_event_pretrend.R"),
                str(full / "event_pretrend_intervals.csv"),
                str(analysis / "event_pretrend_results.json"),
            ],
            r_work,
        ),
        (
            [
                rscript,
                str(revision / "run_pro_event_supplement.R"),
                str(root),
                str(full / "event_pretrend_intervals.csv"),
                str(full / "event_pretrend_subwindows.csv"),
                str(analysis),
            ],
            r_work,
        ),
        (
            [
                rscript,
                str(revision / "strict_volatile_event_sensitivity.R"),
                str(analysis / "strict_volatile_event_sensitivity.json"),
            ],
            r_work,
        ),
        (
            [
                sys.executable,
                str(revision / "negative_control_lead_analysis.py"),
                "--event-intervals",
                str(full / "event_intervals.csv"),
                "--output-csv",
                str(analysis / "event_lead_negative_control.csv"),
                "--output-json",
                str(analysis / "event_lead_negative_control.json"),
            ],
            root,
        ),
        (
            [
                sys.executable,
                str(revision / "build_figures_and_tables.py"),
                "--primary-curve",
                str(analysis / "standardized_map_ra_curve.csv"),
                "--analysis",
                str(analysis),
                "--pro-analysis",
                str(analysis),
                "--pretrend-intervals",
                str(full / "event_pretrend_intervals.csv"),
                "--event-intervals",
                str(full / "event_intervals.csv"),
                "--event-results-80",
                str(analysis / "r_event_analysis_results_fraction_0.80.json"),
                "--output",
                str(args.figure_output.resolve()),
            ],
            root,
        ),
    ]

    execute = args.execute and not args.dry_run
    execution_environment = os.environ.copy()
    execution_environment["VITALDB_RELEASE_ASSETS"] = str(args.figure_output.resolve())
    for index, (command, cwd) in enumerate(commands, start=1):
        print(f"[{index}/{len(commands)}] cwd={cwd}")
        print(quoted(command), flush=True)
        if execute:
            subprocess.run(command, cwd=cwd, check=True, env=execution_environment)
    if not execute:
        print("Dry run only. Re-run with --execute to perform the analyses.")
    else:
        print("Reported analysis workflow completed successfully.")


if __name__ == "__main__":
    main()
