"""A matching P value must not conceal a different component estimate/cohort."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release_validator", ROOT / "scripts/validate_release.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


@pytest.mark.parametrize("mismatch", [None, "estimate", "cohort", "p_value"])
def test_same_model_required_even_when_other_targets_match(tmp_path, mismatch):
    output = tmp_path / "outputs/maintenance_ra_r_reanalysis_v1"
    output.mkdir(parents=True)
    cohort = {"rows": 20, "events_or_pairs": 8, "surgeries": 6, "subjects": 5}
    expected, alpha, saved, csv_rows = {}, {}, {}, []
    for name, saved_name, csv_name in (
        ("real_event", "real_event_adjusted", "real_event_adjusted"),
        ("paired_real_vs_pseudo", "paired_real_minus_pseudo", "real_minus_pseudo"),
    ):
        changes = {f"response_{index}": {"estimate": 0.01 * index} for index in range(4)}
        estimates = [{"component": "log_alpha", "interval": key, "estimate_log_scale": value["estimate"]}
                     for key, value in changes.items()]
        expected[name] = {"cohort": dict(cohort), "global_p": 0.2}
        saved[saved_name] = {"cohort": dict(cohort), "estimates": estimates}
        csv_rows.extend(f"{csv_name},{r['component']},{r['interval']},{r['estimate_log_scale']}\n"
                        for r in estimates)
        alpha[name] = {
            "cohort": dict(cohort), "interval_changes": changes,
            "global_test_all_four_changes_zero": {"degrees_of_freedom": 4, "p_value": 0.2},
            "consistency": {
                "identical_design_response_weights": True, "script20_cohort_match": True,
                "maximum_coefficient_difference_from_multi_wls": 0,
                "maximum_difference_from_saved_script20_estimates": 0,
                "maximum_event_total_weight_error": 0,
            },
        }
    if mismatch == "estimate":
        alpha["real_event"]["interval_changes"]["response_0"]["estimate"] = 0.5
    elif mismatch == "cohort":
        alpha["real_event"]["cohort"]["rows"] += 1
    elif mismatch == "p_value":
        alpha["real_event"]["global_test_all_four_changes_zero"]["p_value"] = 0.019
    (tmp_path / "expected_results.json").write_text(json.dumps({"spectral_alpha_same_model": expected}))
    (output / "r_event_spectral_results.json").write_text(json.dumps({"event_spectral": {"outcomes": {"absolute_alpha": alpha}}}))
    (output / "r_alpha_nonalpha_decomposition_results.json").write_text(json.dumps({"low_map_event": saved}))
    (output / "r_alpha_nonalpha_event_results.csv").write_text(
        "analysis,component,interval,estimate_log_scale\n" + "".join(csv_rows))
    if mismatch:
        with pytest.raises(AssertionError):
            validator.validate_spectral_consistency(tmp_path)
    else:
        validator.validate_spectral_consistency(tmp_path)
