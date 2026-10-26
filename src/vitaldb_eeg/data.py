"""Access VitalDB catalogues, track metadata, and locally downloaded records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TrackIndexRow:
    caseid: int
    filename: str
    dtstart: float
    dtend: float
    tracks: frozenset[str]
    size_bytes: int


def normalize_track_name(name: str) -> str:
    return name.removeprefix("/")


def read_track_index(path: str | Path) -> list[TrackIndexRow]:
    rows: list[TrackIndexRow] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                raise ValueError(f"Invalid track index row {line_number}: {len(parts)} columns")
            filename, timestamp_text, track_text, size_text = parts
            metadata_and_tracks = track_text.split(",")
            dtstart = np.nan
            dtend = np.nan
            tracks: set[str] = set()
            for item in metadata_and_tracks:
                if item.startswith("#dtstart="):
                    dtstart = float(item.split("=", 1)[1])
                elif item.startswith("#dtend="):
                    dtend = float(item.split("=", 1)[1])
                elif item:
                    tracks.add(normalize_track_name(item))
            rows.append(
                TrackIndexRow(
                    caseid=int(Path(filename).stem),
                    filename=filename,
                    dtstart=float(dtstart if np.isfinite(dtstart) else timestamp_text),
                    dtend=float(dtend),
                    tracks=frozenset(tracks),
                    size_bytes=int(size_text),
                )
            )
    return rows


def clinical_frame(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    required = {
        "caseid",
        "subjectid",
        "age",
        "ane_type",
        "anestart",
        "opstart",
        "caseend",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"clinical_data.csv lacks required columns: {sorted(missing)}")
    frame["caseid"] = pd.to_numeric(frame["caseid"], errors="raise").astype(int)
    frame["subjectid"] = pd.to_numeric(frame["subjectid"], errors="raise").astype(int)
    return frame


def vital_path(vital_root: str | Path, caseid: int) -> Path:
    return Path(vital_root) / f"{int(caseid):04d}.vital"


def import_vitaldb():
    try:
        import vitaldb  # type: ignore
    except ImportError as exc:
        raise RuntimeError("Install dependencies with `uv sync --extra dev`.") from exc
    return vitaldb
