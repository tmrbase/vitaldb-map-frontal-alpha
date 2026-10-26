"""Validate the retained source-to-result path and regenerated results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

PINNED_TRACK_CATALOG_SHA256 = (
    "29374aae9449c18f3d1ea360dbf7712614a70dcfe63eedf639bab27c8a1556c9"
)


def assert_close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(f"{label}: expected {expected}, observed {actual}")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def expected_results(root: Path) -> dict[str, object]:
    return json.loads((root / "expected_results.json").read_text(encoding="utf-8"))


def validate_authoritative_path(root: Path) -> None:
    driver = (root / "scripts" / "run_source_to_result.py").read_text(encoding="utf-8")
    required = (
        "fetch_vitaldb_inputs.py",
        "prepare-volatile-prescreen",
        "prepare-maintenance-windows",
        "prepare-low-map-events",
        "run_analysis.py",
    )
    for value in required:
        if value not in driver:
            raise AssertionError(f"Source-to-result driver is missing: {value}")

    fetcher = (root / "scripts" / "fetch_vitaldb_inputs.py").read_text(encoding="utf-8")
    for value in (
        'DATASET_VERSION = "1.0.0"',
        "https://physionet-open.s3.amazonaws.com/vitaldb/",
        "SHA256SUMS_SHA256",
        "ThreadPoolExecutor",
        "download_with_resume",
    ):
        if value not in fetcher:
            raise AssertionError(f"Pinned parallel v1.0.0 source retrieval is missing: {value}")
    for value in ("vitaldb.find_cases(", "VitalFile(caseid"):
        if value in fetcher:
            raise AssertionError(
                "A retrieval path outside the pinned PhysioNet v1.0.0 source was detected: "
                f"{value}"
            )

    catalog = root / "source_catalog" / "vitaldb_v1.0.0_track_index.tsv.gz"
    if not catalog.is_file():
        raise AssertionError("The pinned v1.0.0 track catalogue is missing")
    observed_catalog_digest = hashlib.sha256(catalog.read_bytes()).hexdigest()
    if observed_catalog_digest != PINNED_TRACK_CATALOG_SHA256:
        raise AssertionError(
            "The bundled v1.0.0 track catalogue does not match its pinned SHA-256 digest"
        )

    r_project = root / "r_analysis" / "maintenance_ra"
    r_lock = json.loads((r_project / "renv.lock").read_text(encoding="utf-8-sig"))
    if r_lock["R"]["Version"] != "4.4.1":
        raise AssertionError("The R lockfile is not pinned to R 4.4.1")
    for record in r_lock["Packages"].values():
        for field in ("Depends", "Imports", "Suggests", "LinkingTo", "Enhances"):
            if field in record and not isinstance(record[field], str):
                raise AssertionError(
                    f"The R lockfile is not restorable: {record['Package']} {field} is not scalar"
                )
    r_profile = (r_project / ".Rprofile").read_text(encoding="utf-8")
    if "RENV_CONFIG_CACHE_ENABLED" not in r_profile:
        raise AssertionError("The R project does not disable the mutable global renv cache")

    analysis_driver = (root / "scripts" / "run_analysis.py").read_text(encoding="utf-8")
    for required_step in ("generate_patient_influence.py", "03_model_sensitivities.R", "06_event_spectral.R"):
        if required_step not in analysis_driver:
            raise AssertionError(f"Reported workflow step is missing: {required_step}")
    if analysis_driver.index("generate_patient_influence.py") > analysis_driver.index("03_model_sensitivities.R"):
        raise AssertionError("Patient influence must be generated before the exclusion/refit")
    if "27_optimizer_sensitivity.R" in analysis_driver:
        raise AssertionError("The public workflow still invokes the unreported optimizer audit")
    if (root / "r_analysis" / "maintenance_ra" / "27_optimizer_sensitivity.R").exists():
        raise AssertionError("The unreported optimizer audit remains in the release repository")

    nonalpha = (r_project / "R" / "event_components.R").read_text(encoding="utf-8")
    for value in (
        "actual_nonalpha <- subwindows[",
        "average it within each event interval",
    ):
        if value not in nonalpha:
            raise AssertionError(f"Subwindow-level non-alpha implementation is missing: {value}")
    for script in ("20_alpha_nonalpha_decomposition.R", "06_event_spectral.R"):
        if 'source("R/event_components.R")' not in (r_project / script).read_text(encoding="utf-8"):
            raise AssertionError(f"Shared component preparation missing in {script}")
    if analysis_driver.index("20_alpha_nonalpha_decomposition.R") > analysis_driver.index("06_event_spectral.R"):
        raise AssertionError("Component estimates must precede the same-model joint tests")

    observability = (
        root / "r_analysis" / "maintenance_ra" / "26_window_observability_audit.R"
    ).read_text(encoding="utf-8")
    for unreported in ("IPOW", "selection_formula", "TableS20", "FigureS4_window_observability"):
        if unreported in observability:
            raise AssertionError(f"Unreported observability analysis remains: {unreported}")

    common = (root / "r_analysis" / "maintenance_ra" / "R" / "common.R").read_text(encoding="utf-8")
    if "splines::ns(" not in common:
        raise AssertionError("The primary design does not use splines::ns")

    robust_code = (
        root / "revision_sensitivities" / "primary_map_contrast_patient_CR1S.R"
    ).read_text(encoding="utf-8")
    for value in ("clubSandwich::vcovCR(", 'type = "CR1S"'):
        if value not in robust_code:
            raise AssertionError(f"Package-based cluster-robust implementation is missing: {value}")
    print("Source-to-result path: retained implementation checks passed")


def validate_generated_numbers(root: Path) -> None:
    expected = expected_results(root)
    analysis = root / "outputs" / "maintenance_ra_r_reanalysis_v1"
    result = json.loads(
        (analysis / "r_random_slope_primary_results.json").read_text(encoding="utf-8")
    )
    if result["cohort"] != expected["cohort"]:
        raise AssertionError(
            f"Generated cohort mismatch: expected {expected['cohort']}, observed {result['cohort']}"
        )

    selection = json.loads(
        (analysis / "r_table1_missingness_results.json").read_text(encoding="utf-8")
    )["selection_flow"]
    for key, target in expected["selection_flow"].items():
        if int(selection[key]) != int(target):
            raise AssertionError(
                f"Selection-flow mismatch for {key}: expected {target}, observed {selection[key]}"
            )
    for key, label in (
        ("case_equal_G", "G"),
        ("delta_marginal_r2", "delta marginal R2"),
        ("delta_conditional_r2", "delta conditional R2"),
    ):
        source = (
            result["residual_reduction"][key] if key == "case_equal_G" else result["mixed_r2"][key]
        )
        assert_close(source, float(expected["model"][key]), 1e-8, label)

    robust = read_csv_rows(analysis / "primary_map_contrast_patient_CR1S.csv")
    engines = {row["covariance_engine"] for row in robust}
    if engines != {"clubSandwich::vcovCR(type='CR1S')"}:
        raise AssertionError(f"Unexpected cluster-robust covariance engine: {engines}")
    for target, values in expected["patient_cluster_CR1S_contrasts_percentage_points"].items():
        row = next(item for item in robust if int(item["target_map_mmhg"]) == int(target))
        columns = (
            "contrast_percentage_points",
            "patient_cluster_CR1S_lower_95ci_percentage_points",
            "patient_cluster_CR1S_upper_95ci_percentage_points",
        )
        for column, expected_value in zip(columns, values, strict=True):
            assert_close(
                float(row[column]), float(expected_value), 1e-6, f"{target} vs 65 {column}"
            )

    negative = read_csv_rows(analysis / "event_lead_negative_control.csv")
    row = next(
        item
        for item in negative
        if item["contrast"] == "paired_real_minus_reference_lead_to_baseline"
    )
    assert_close(
        float(row["estimate_pp"]),
        float(expected["negative_control"]["estimate_percentage_points"]),
        1e-8,
        "negative control",
    )
    if int(row["events"]) != int(expected["negative_control"]["events"]) or int(
        row["patients"]
    ) != int(expected["negative_control"]["patients"]):
        raise AssertionError("Negative-control sample size mismatch")

    spectral = read_csv_rows(analysis / "r_alpha_nonalpha_event_results.csv")
    labels = {
        "real_event_adjusted": "real_event",
        "real_minus_pseudo": "real_minus_reference",
    }
    for analysis_name, expected_name in labels.items():
        row = next(
            item
            for item in spectral
            if item["analysis"] == analysis_name
            and item["interval"] == "response_0_1"
            and item["component"] == "log_nonalpha"
        )
        observed = (
            float(row["percent_power_change"]),
            float(row["percent_power_change_lower"]),
            float(row["percent_power_change_upper"]),
        )
        for index, (actual, target) in enumerate(
            zip(
                observed,
                expected["actual_nonalpha_first_minute_percent_change"][expected_name],
                strict=True,
            )
        ):
            assert_close(actual, float(target), 1e-6, f"{expected_name} non-alpha {index}")

    assets = root / "outputs" / "release_assets"
    for name in (
        "FigureS1_cohort_flow.pdf",
        "FigureS2_residual_diagnostics.pdf",
        "FigureS3_window_observability.pdf",
        "FigureS4_event_trajectory_diagnostic_BJA.pdf",
    ):
        if not (assets / name).is_file():
            raise AssertionError(f"Missing regenerated supplementary figure: {name}")
    figure_b = read_csv_rows(assets / "FigureS4B_source_data.csv")
    map_row = next(
        item
        for item in figure_b
        if item["analysis"] == "Real sustained low-MAP event"
        and item["interval"] == "pretrend_m1_0"
    )
    figure_c = read_csv_rows(assets / "FigureS4C_source_data.csv")
    ra_row = next(
        item
        for item in figure_c
        if item["analysis"] == "Real sustained low-MAP event"
        and item["interval"] == "pretrend_m1_0"
    )
    assert_close(
        float(map_row["estimate"]),
        float(expected["final_pre_index_minute"]["map_change_mmhg"]),
        1e-8,
        "final pre-index MAP",
    )
    assert_close(
        float(ra_row["estimate"]),
        float(expected["final_pre_index_minute"]["ra_change_percentage_points"]),
        1e-8,
        "final pre-index RA",
    )
    print("Generated numerical checks: all expected results verified")


def validate_diagnostics(root: Path) -> None:
    """Compare diagnostic outputs with fixed manuscript-precision targets."""
    expected = expected_results(root)["diagnostics"]
    analysis = root / "outputs" / "maintenance_ra_r_reanalysis_v1"
    result = json.loads((analysis / "r_model_sensitivity_results.json").read_text())
    slope = result["random_map_slope"]
    for key in ("median_reliability", "reliability_q25", "reliability_q75"):
        assert_close(slope["individual_estimation"][key], expected[key], 0.0005, key)
    for key in ("all_spearman", "map_range_ge_10_spearman"):
        assert_close(slope["temporal_reproducibility"][key], expected[key], 0.0005, key)
    trimmed = result["summaries"]["top_one_percent_influence_removed"]
    for key, value in expected["influence_trimmed"].items():
        if trimmed[key] != value:
            raise AssertionError(f"Influence-trimmed {key}: expected {value}, observed {trimmed[key]}")
    for key, bound in (("overall_map_p", "influence_trimmed_overall_p_below"),
                       ("nonlinearity_p", "influence_trimmed_nonlinearity_p_below")):
        if not 0 <= trimmed[key] < expected[bound]:
            raise AssertionError(f"Influence-trimmed {key} no longer matches Table SIII")
    ranking = read_csv_rows(root / "outputs/maintenance_ra_model_defense_v1/subject_deletion_influence.csv")
    ordered = sorted(ranking, key=lambda row: -float(row["maximum_absolute_curve_change_logit_ra"]))
    ids = [int(row["subjectid"]) for row in ordered[:math.ceil(0.01 * len(ordered))]]
    if [int(value) for value in trimmed["excluded_subjectids"]] != ids:
        raise AssertionError("The R refit did not use the generated top-one-percent patient ranking")
    validate_spectral_consistency(root)
    print("Diagnostic checks: Tables SII, SIII, and spectral Global P in SVI passed")


def validate_spectral_consistency(root: Path) -> None:
    """Check shared-model coefficients and cohorts, not just matching P values."""
    expected = expected_results(root)["spectral_alpha_same_model"]
    analysis = root / "outputs" / "maintenance_ra_r_reanalysis_v1"
    spectral = json.loads((analysis / "r_event_spectral_results.json").read_text())
    alpha = spectral["event_spectral"]["outcomes"]["absolute_alpha"]
    saved = json.loads((analysis / "r_alpha_nonalpha_decomposition_results.json").read_text())["low_map_event"]
    table = read_csv_rows(analysis / "r_alpha_nonalpha_event_results.csv")
    for name, saved_name, csv_name in (
        ("real_event", "real_event_adjusted", "real_event_adjusted"),
        ("paired_real_vs_pseudo", "paired_real_minus_pseudo", "real_minus_pseudo"),
    ):
        result, reference, target = alpha[name], saved[saved_name], expected[name]
        if result["cohort"] != reference["cohort"] or result["cohort"] != target["cohort"]:
            raise AssertionError(f"{name}: joint-test and script-20 cohort mismatch")
        audit = result["consistency"]
        if not audit["identical_design_response_weights"] or not audit["script20_cohort_match"]:
            raise AssertionError(f"{name}: shared-model consistency check failed")
        for field in ("maximum_coefficient_difference_from_multi_wls",
                      "maximum_difference_from_saved_script20_estimates",
                      "maximum_event_total_weight_error"):
            assert_close(audit[field], 0, 1e-8, f"{name} {field}")
        if len(result["interval_changes"]) != 4:
            raise AssertionError(f"{name}: expected four interval estimates")
        for interval, values in result["interval_changes"].items():
            rows = [r for r in table if r["analysis"] == csv_name and r["component"] == "log_alpha"
                    and r["interval"] == interval]
            refs = [r for r in reference["estimates"] if r["component"] == "log_alpha"
                    and r["interval"] == interval]
            if len(rows) != 1 or len(refs) != 1:
                raise AssertionError(f"{name} {interval}: missing or duplicate estimate")
            for record in (rows[0], refs[0]):
                assert_close(values["estimate"], float(record["estimate_log_scale"]), 1e-8,
                             f"{name} {interval} script-20 point estimate")
        test = result["global_test_all_four_changes_zero"]
        if test["degrees_of_freedom"] != 4:
            raise AssertionError("Table SVI spectral tests must jointly test four intervals")
        assert_close(test["p_value"], target["global_p"], 1e-8, f"{name} Global P")
    print("Table SVI alpha: script-20 coefficients, cohorts and same-model Global P verified")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the public workflow and, optionally, regenerated numerical results."
    )
    parser.add_argument(
        "--generated-results",
        action="store_true",
        help="Also compare locally regenerated outputs with expected_results.json.",
    )
    parser.add_argument(
        "--diagnostics", action="store_true",
        help="Check regenerated reliability, split-half, influence refit and spectral joint tests only.",
    )
    parser.add_argument(
        "--spectral-consistency", action="store_true",
        help="Check only Table SVI alpha WLS coefficients, cohorts and joint tests.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    validate_authoritative_path(root)
    if args.generated_results:
        validate_generated_numbers(root)
    if args.generated_results or args.diagnostics:
        validate_diagnostics(root)
    elif args.spectral_consistency:
        validate_spectral_consistency(root)


if __name__ == "__main__":
    main()
