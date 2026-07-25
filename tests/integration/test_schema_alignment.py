"""Checked-in schemas must match current Pydantic models."""

import subprocess
import sys
from pathlib import Path


def test_checked_in_json_schemas_have_no_drift() -> None:
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [sys.executable, "scripts/validate_schemas.py"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
