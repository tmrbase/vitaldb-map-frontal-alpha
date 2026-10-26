"""Load analysis configuration and resolve repository-local data paths."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ProjectPaths:
    workspace: Path
    data_root: Path
    output_root: Path

    @property
    def clinical_data(self) -> Path:
        return self.data_root / "clinical_data.csv"

    @property
    def track_index(self) -> Path:
        return self.data_root / "vital_files" / "trks.tsv"

    @property
    def vital_root(self) -> Path:
        return self.data_root / "vital_files"


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise TypeError(f"Configuration must be a mapping: {config_path}")
    cfg["_config_path"] = str(config_path)
    cfg["_workspace"] = str(config_path.parent.parent.resolve())
    return cfg


def project_paths(cfg: dict[str, Any]) -> ProjectPaths:
    workspace = Path(cfg["_workspace"])

    def resolve(value: str) -> Path:
        candidate = Path(value)
        return candidate.resolve() if candidate.is_absolute() else (workspace / candidate).resolve()

    paths = cfg["paths"]
    return ProjectPaths(
        workspace=workspace,
        data_root=resolve(paths["data_root"]),
        output_root=resolve(paths["output_root"]),
    )


def ensure_output_directories(cfg: dict[str, Any]) -> ProjectPaths:
    paths = project_paths(cfg)
    paths.output_root.mkdir(parents=True, exist_ok=True)
    return paths
