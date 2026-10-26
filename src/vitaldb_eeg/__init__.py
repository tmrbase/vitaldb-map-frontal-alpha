"""Reproducible VitalDB MAP–frontal-alpha analysis."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("vitaldb-map-ra-analysis")
except PackageNotFoundError:
    __version__ = "0.1.0"
