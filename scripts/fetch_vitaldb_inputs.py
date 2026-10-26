"""Download the immutable PhysioNet VitalDB v1.0.0 source records.

This downloader fixes the DOI-versioned source URL, selects candidates from a
version-matched track catalogue, downloads records in parallel, resumes partial
files, and verifies every file against the official PhysioNet SHA-256 manifest.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

from vitaldb_eeg.config import load_config
from vitaldb_eeg.data import read_track_index
from vitaldb_eeg.maintenance_ra_extraction import _maintenance_centers

DATASET_VERSION = "1.0.0"
DATASET_DOI = "10.13026/czw8-9p62"
DATASET_PAGE = f"https://physionet.org/content/vitaldb/{DATASET_VERSION}/"
# PhysioNet documents this public S3 mirror for bulk v1.0.0 retrieval.
BASE_URL = f"https://physionet-open.s3.amazonaws.com/vitaldb/{DATASET_VERSION}"
SHA256SUMS_SHA256 = "874732d0548b5ba4de2efe851375c4f25f45488a2a19c771e2f9ab0288a5de1d"
CHUNK_BYTES = 1024 * 1024
MAX_ATTEMPTS = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(path: Path, expected_sha256: str, expected_bytes: int | None = None) -> bool:
    if not path.exists() or not path.is_file():
        return False
    if expected_bytes is not None and path.stat().st_size != expected_bytes:
        return False
    return sha256(path).lower() == expected_sha256.lower()


def download_with_resume(
    *,
    url: str,
    destination: Path,
    expected_sha256: str,
    expected_bytes: int | None,
    force: bool,
) -> str:
    """Download one file atomically, resuming an existing partial response."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    if not force and verify(destination, expected_sha256, expected_bytes):
        return "verified-existing"

    partial = destination.with_name(f"{destination.name}.partial")
    if force:
        partial.unlink(missing_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": "vitaldb-map-frontal-alpha-reproduction/2.0"}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", response.getcode())
                append = bool(offset and status == 206)
                mode = "ab" if append else "wb"
                with partial.open(mode) as handle:
                    shutil.copyfileobj(response, handle, length=CHUNK_BYTES)

            if verify(partial, expected_sha256, expected_bytes):
                os.replace(partial, destination)
                return "downloaded"
            last_error = RuntimeError(
                f"Checksum or size mismatch after downloading {destination.name}"
            )
            partial.unlink(missing_ok=True)
        except (OSError, urllib.error.URLError, RuntimeError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 416:
                partial.unlink(missing_ok=True)
        if attempt < MAX_ATTEMPTS:
            time.sleep(min(2**attempt, 10))
    raise RuntimeError(f"Failed to download {url} after {MAX_ATTEMPTS} attempts") from last_error


def official_hashes(path: Path) -> dict[str, str]:
    if sha256(path).lower() != SHA256SUMS_SHA256:
        raise RuntimeError("PhysioNet SHA256SUMS.txt does not match the pinned v1.0.0 digest")
    hashes: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split(maxsplit=1)
        hashes[relative.strip().replace("\\", "/")] = digest.lower()
    return hashes


def unpack_track_catalog(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", encoding="utf-8") as input_handle:
        destination.write_text(input_handle.read(), encoding="utf-8", newline="")


def discover_candidates(
    cfg: dict[str, Any], clinical: pd.DataFrame, track_index_path: Path
) -> tuple[list[int], pd.DataFrame]:
    required = {
        "BIS/EEG1_WAV",
        "BIS/SQI",
        "BIS/SR",
        "Solar8000/ART_MBP",
        "Orchestra/PPF20_RATE",
    }
    remifentanil = {"Orchestra/RFTN20_RATE", "Orchestra/RFTN50_RATE"}
    indexed = {row.caseid: row for row in read_track_index(track_index_path)}

    def track_eligible(caseid: int) -> bool:
        row = indexed.get(int(caseid))
        return bool(row and required.issubset(row.tracks) and remifentanil.intersection(row.tracks))

    numeric = clinical[["age", "height", "weight", "anestart", "opstart", "opend"]].apply(
        pd.to_numeric, errors="coerce"
    )
    eligible = numeric.notna().all(axis=1)
    eligible &= numeric["age"].ge(float(cfg["cohort"]["adult_age"]))
    eligible &= numeric["height"].gt(0) & numeric["weight"].gt(0)
    eligible &= clinical["ane_type"].eq(cfg["cohort"]["primary_ane_type"])
    eligible &= clinical["sex"].astype(str).str.upper().isin(
        ["M", "F", "MALE", "FEMALE", "0", "1"]
    )
    eligible &= clinical["caseid"].map(track_eligible)
    eligible &= clinical.apply(
        lambda row: (
            len(_maintenance_centers(row.to_dict(), cfg["analysis"]))
            >= int(cfg["analysis"]["minimum_valid_windows"])
        ),
        axis=1,
    )
    selected = clinical.loc[eligible].copy().sort_values("caseid")
    return selected["caseid"].astype(int).tolist(), selected


def main() -> None:
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=f"Download and verify PhysioNet VitalDB v{DATASET_VERSION} source records."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=repository_root / "configs" / "analysis_public.yaml",
        help="Public analysis configuration file.",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=repository_root / "data" / "vitaldb",
        help="Directory for downloaded records, source metadata, and the local manifest.",
    )
    parser.add_argument(
        "--workers", type=int, default=8, help="Concurrent PhysioNet downloads."
    )
    parser.add_argument("--force", action="store_true", help="Replace verified source records.")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")

    data_root = args.data_root.resolve()
    data_root.mkdir(parents=True, exist_ok=True)
    vital_root = data_root / "vital_files"
    vital_root.mkdir(parents=True, exist_ok=True)

    sums_path = data_root / "SHA256SUMS.txt"
    if args.force or not verify(sums_path, SHA256SUMS_SHA256):
        download_with_resume(
            url=f"{BASE_URL}/SHA256SUMS.txt",
            destination=sums_path,
            expected_sha256=SHA256SUMS_SHA256,
            expected_bytes=None,
            force=args.force,
        )
    hashes = official_hashes(sums_path)

    clinical_path = data_root / "clinical_data.csv"
    clinical_relative = "clinical_data.csv"
    download_with_resume(
        url=f"{BASE_URL}/{clinical_relative}",
        destination=clinical_path,
        expected_sha256=hashes[clinical_relative],
        expected_bytes=None,
        force=args.force,
    )

    bundled_catalog = (
        repository_root / "source_catalog" / "vitaldb_v1.0.0_track_index.tsv.gz"
    )
    if not bundled_catalog.exists():
        raise FileNotFoundError(f"Missing version-matched track catalogue: {bundled_catalog}")
    full_index = data_root / "vitaldb_v1.0.0_track_index.tsv"
    unpack_track_catalog(bundled_catalog, full_index)

    clinical = pd.read_csv(clinical_path, low_memory=False)
    caseids, selected = discover_candidates(load_config(args.config), clinical, full_index)
    if not caseids:
        raise ValueError("No v1.0.0 candidate case identifiers selected")

    all_index_rows = {row.caseid: row for row in read_track_index(full_index)}
    total_bytes = sum(all_index_rows[caseid].size_bytes for caseid in caseids)
    print(
        f"PhysioNet VitalDB v{DATASET_VERSION}: {len(caseids)} candidate files, "
        f"{total_bytes / 1024**3:.1f} GiB; {args.workers} parallel downloads",
        flush=True,
    )

    def fetch_case(caseid: int) -> tuple[int, str]:
        relative = f"vital_files/{caseid:04d}.vital"
        destination = vital_root / f"{caseid:04d}.vital"
        status = download_with_resume(
            url=f"{BASE_URL}/{relative}",
            destination=destination,
            expected_sha256=hashes[relative],
            expected_bytes=all_index_rows[caseid].size_bytes,
            force=args.force,
        )
        return caseid, status

    completed = 0
    completed_bytes = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch_case, caseid): caseid for caseid in caseids}
        for future in as_completed(futures):
            caseid, status = future.result()
            completed += 1
            completed_bytes += all_index_rows[caseid].size_bytes
            print(
                f"[{completed}/{len(caseids)}; {completed_bytes / 1024**3:.1f}/"
                f"{total_bytes / 1024**3:.1f} GiB] case {caseid:04d}: {status}",
                flush=True,
            )

    selected.to_csv(data_root / "candidate_clinical_data.csv", index=False)
    selected_ids = set(caseids)
    with full_index.open("r", encoding="utf-8") as input_handle, (
        vital_root / "trks.tsv"
    ).open("w", encoding="utf-8", newline="") as output_handle:
        for line in input_handle:
            caseid = int(Path(line.split("\t", 1)[0]).stem)
            if caseid in selected_ids:
                output_handle.write(line)

    manifest_rows = [
        {
            "caseid": caseid,
            "repository_relative_path": f"data/vitaldb/vital_files/{caseid:04d}.vital",
            "bytes": all_index_rows[caseid].size_bytes,
            "sha256": hashes[f"vital_files/{caseid:04d}.vital"],
            "source_version": DATASET_VERSION,
            "source_doi": DATASET_DOI,
        }
        for caseid in caseids
    ]
    pd.DataFrame(manifest_rows).to_csv(data_root / "source_download_manifest.csv", index=False)
    (data_root / "source_release.json").write_text(
        json.dumps(
            {
                "dataset": "VitalDB, a high-fidelity multi-parameter vital signs database in surgical patients",
                "version": DATASET_VERSION,
                "doi": DATASET_DOI,
                "release_page": DATASET_PAGE,
                "base_url": BASE_URL,
                "candidate_files": len(caseids),
                "candidate_bytes": total_bytes,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Prepared and verified {len(caseids)} v{DATASET_VERSION} cases under {data_root}")


if __name__ == "__main__":
    main()
