"""Manifest hashing and timestamp provenance."""

from datetime import datetime
from pathlib import Path

import pytest
from pydantic import ValidationError

from protein_interaction_hunter.manifest import (
    build_input_file_manifest,
    build_run_manifest,
    sha256_file,
)
from protein_interaction_hunter.models.run import InputFileManifest


def test_sha256_and_input_manifest(valid_config_path: Path) -> None:
    manifest = build_input_file_manifest("config", valid_config_path, required=True)
    assert manifest.exists is True
    assert manifest.sha256 == sha256_file(valid_config_path)
    assert manifest.size_bytes == valid_config_path.stat().st_size
    assert manifest.modified_time is not None
    assert manifest.modified_time.utcoffset() is not None


def test_run_manifest_timestamp_is_timezone_aware(valid_config_path: Path) -> None:
    manifest = build_run_manifest(
        run_id="run-001",
        run_name="fixture",
        config_path=valid_config_path,
        config_snapshot_path=None,
        input_files=[],
        random_seed=0,
        command_line=["validate-inputs"],
    )
    assert manifest.started_at.utcoffset() is not None
    assert manifest.config_sha256 == sha256_file(valid_config_path)


def test_naive_input_timestamp_is_rejected() -> None:
    with pytest.raises(ValidationError, match="timezone-aware"):
        InputFileManifest(
            logical_name="x",
            path=Path("x"),
            exists=True,
            required=True,
            modified_time=datetime(2026, 7, 25),
        )
