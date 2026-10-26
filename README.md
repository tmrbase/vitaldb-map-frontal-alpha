# VitalDB MAP–frontal alpha analysis

Canonical repository: https://github.com/tmrbase/vitaldb-map-frontal-alpha

Code accompanying the **Minerva Anestesiologica manuscript, September 30, 2026**.

Source-to-result code for:

> *Mean arterial pressure and frontal relative alpha power during propofol–remifentanil-based anesthesia*

This release reconstructs the eligible cohort from the immutable PhysioNet
VitalDB v1.0.0 release, quality-controlled EEG/MAP windows, sustained low-MAP
event intervals, reported statistical analyses, figures, and aggregate source
tables. Primary statistical modeling is performed in R; Python is used for
version-matched public-data discovery and download, signal processing, dataset
construction, HRV extraction, the original patient-deletion influence diagnostic,
figure assembly, and validation. The influence diagnostic ranks patients for the
reported sensitivity; the model after exclusion is refitted in R.

See [the implementation map](docs/METHODS_CODE_MAP.md) for Methods S5,
Tables SII/SIII/SVI, the execution steps, and the distinction between interval CIs
and joint tests. [Validation records](docs/VALIDATION.md) describe the checks
actually performed and the use of previously generated local inputs.

## Reproduce the analysis

Requirements are Python 3.12, [uv](https://docs.astral.sh/uv/), R 4.4, at
least 45 GB of free storage (50 GB recommended), and network access to PhysioNet.

```powershell
uv sync --frozen --extra dev
Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv'); renv::restore(project='r_analysis/maintenance_ra', prompt=FALSE)"
uv run python scripts/run_source_to_result.py --execute --workers 8
```

The download stage fixes the source to [PhysioNet VitalDB
v1.0.0](https://physionet.org/content/vitaldb/1.0.0/) (DOI:
10.13026/czw8-9p62). It selects 1,577 candidate records from the bundled
version-matched track catalogue, downloads approximately 35.7 GiB in parallel,
resumes partial transfers, and verifies every file against the official
SHA-256 manifest. Bulk records are read from the public AWS S3 mirror documented
on the PhysioNet release page.

Without `--execute`, the driver only prints the ordered commands. Use
`--skip-fetch` when the source files created by
`scripts/fetch_vitaldb_inputs.py` are already present, and `--force` to rebuild
verified cached extraction products. Source records that pass their v1.0.0
checksum are retained even with `--force`; use the separate `--force-download`
option only when a fresh network transfer is specifically required.

Regenerated Figure 1, Supplementary Figures S1–S4, and their aggregate source tables are written to `outputs/release_assets/`.

## Fast checks

These checks do not download VitalDB records:

```powershell
uv run python -m pytest -q
uv run python -m ruff check scripts revision_sensitivities src tests
uv run python scripts/run_source_to_result.py --dry-run
uv run python scripts/validate_release.py
Rscript tests/test_patient_influence.R
```

## Repository layout

- `configs/analysis_public.yaml`: public analysis specification and quality rules.
- `source_catalog/`: compact public v1.0.0 header catalogue used to avoid a 95.4 GB full-archive download.
- `src/vitaldb_eeg/`: VitalDB access, signal processing, pharmacokinetic reconstruction, window construction, and event extraction.
- `scripts/`: data acquisition, covariate extraction, ordered execution, and validation entry points.
- `r_analysis/maintenance_ra/`: primary and reported sensitivity analyses.
- `revision_sensitivities/`: reported threshold-support and event-diagnostic analyses.
- `expected_results.json`: numerical targets used after a complete rebuild.
- `tests/`: invariant tests for signal processing, event timing, interval completeness, and non-alpha calculation.

`r_analysis/maintenance_ra/RUN_ORDER.md` maps every statistical script to the corresponding manuscript, table, or figure result. `REPRODUCIBILITY.md` documents the data flow, validation boundaries, and expected runtime.

## Data policy

No VitalDB waveform, clinical value, or derived case-, window-, subwindow-, or
event-level dataset is redistributed. The bundled catalogue contains only
public v1.0.0 filenames and track-header metadata needed to reconstruct the
result-independent download set. Local `data/` and `outputs/` directories are
excluded from Git.

In particular, `subject_deletion_influence.csv`, excluded patient identifiers,
model objects, and all generated patient-level diagnostics stay local. Only
their generation code and aggregate validation targets are distributed.

The software is released under the MIT License. The bundled public header
catalogue is derived from PhysioNet VitalDB v1.0.0 and remains subject to the
source dataset's Creative Commons Attribution 4.0 licence.
