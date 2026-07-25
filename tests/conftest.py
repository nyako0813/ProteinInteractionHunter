"""Shared fixture paths and representative models."""

from pathlib import Path

import pytest


@pytest.fixture
def fixture_dir() -> Path:
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def valid_config_path(fixture_dir: Path) -> Path:
    return fixture_dir / "config.valid.yaml"
