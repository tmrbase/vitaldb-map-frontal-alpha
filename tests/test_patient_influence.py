"""Check the deletion diagnostic against direct GLS on a synthetic cohort."""

import numpy as np
import pandas as pd
from scipy.linalg import block_diag

from vitaldb_eeg.patient_influence import (
    Car1Fit,
    _subject_influence,
    map_within_contrast_matrix,
)


def test_deletion_removes_every_surgery_of_a_patient_and_matches_direct_gls(tmp_path):
    rng = np.random.default_rng(20)
    caseid = np.repeat(np.arange(5), 4)
    subjectid = np.repeat([10, 10, 20, 30, 40], 4)
    times = np.tile([0.0, 300.0, 650.0, 1000.0], 5)
    x = np.column_stack([np.ones(20), rng.normal(size=(20, 3))])
    y = x @ np.array([0.1, 0.4, -0.2, 0.6]) + rng.normal(scale=0.2, size=20)
    names = ("intercept", "map_basis1_within", "map_basis2_within", "map_basis3_within")
    bundle = {
        "data": pd.DataFrame({"caseid": caseid, "subjectid": subjectid, "time_seconds": times}),
        "nonlinear": x, "y": y, "nonlinear_names": names,
        "knots": {"map": [50, 65, 80, 95]},
        "scaling": {name: {"standard_deviation": 1.0} for name in names[1:]},
    }
    fit = Car1Fit(
        np.zeros(4), np.eye(4), 0.0, 0.25, 0.4, 0.7, True, "synthetic", names
    )
    distance = np.abs(times[:4, None] - times[None, :4]) / 300
    covariance = block_diag(*[0.25 + 0.4 * 0.7**distance + 1e-8 * np.eye(4)] * 5)

    def direct_beta(keep):
        local_x = x[keep]
        inverse = np.linalg.inv(covariance[np.ix_(keep, keep)])
        return np.linalg.solve(local_x.T @ inverse @ local_x, local_x.T @ inverse @ y[keep])

    contrast = map_within_contrast_matrix(bundle, np.arange(50, 91), reference_mmhg=65)
    full = contrast @ direct_beta(np.ones(20, dtype=bool))
    summary, rows = _subject_influence(bundle, fit, tmp_path)
    assert summary["subjects"] == 4
    for row in rows.itertuples():
        expected = contrast @ direct_beta(subjectid != row.subjectid) - full
        assert np.isclose(row.maximum_absolute_curve_change_logit_ra, np.abs(expected).max())
    patient = rows.set_index("subjectid").loc[10]
    assert patient["surgeries"] == 2 and patient["windows"] == 8
    assert np.allclose(contrast[15], 0)  # Reference at MAP 65.
