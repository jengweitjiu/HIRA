"""Shared fixtures for HIRA validation tests."""

from pathlib import Path

import pytest

# Project root: D:\HIRA
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"


@pytest.fixture
def results_dir():
    """Return the results directory path."""
    return RESULTS_DIR
