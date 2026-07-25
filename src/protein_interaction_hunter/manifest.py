"""Deterministic file fingerprints and run-manifest helpers."""

import hashlib
import importlib.metadata
import platform as platform_module
import sys
from datetime import UTC, datetime
from pathlib import Path

from protein_interaction_hunter import __version__
from protein_interaction_hunter.models.enums import RunStatus
from protein_interaction_hunter.models.run import InputFileManifest, RunManifest


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Hash a file without loading it entirely into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_input_file_manifest(
    logical_name: str,
    path: Path,
    *,
    required: bool,
) -> InputFileManifest:
    """Fingerprint a local input; missing optional files remain representable."""
    resolved = path.expanduser().resolve()
    if not resolved.exists():
        return InputFileManifest(
            logical_name=logical_name,
            path=resolved,
            exists=False,
            required=required,
        )
    stat = resolved.stat()
    return InputFileManifest(
        logical_name=logical_name,
        path=resolved,
        sha256=sha256_file(resolved),
        size_bytes=stat.st_size,
        modified_time=datetime.fromtimestamp(stat.st_mtime, tz=UTC),
        exists=True,
        required=required,
    )


def _dependency_versions() -> dict[str, str]:
    packages = ("pydantic", "PyYAML", "typer", "openpyxl")
    return {
        package: importlib.metadata.version(package)
        for package in packages
        if _distribution_exists(package)
    }


def _distribution_exists(package: str) -> bool:
    try:
        importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return False
    return True


def build_run_manifest(
    *,
    run_id: str,
    run_name: str,
    config_path: Path,
    config_snapshot_path: Path | None,
    input_files: list[InputFileManifest],
    random_seed: int,
    command_line: list[str] | None = None,
) -> RunManifest:
    """Build an initialized, timezone-aware local run manifest."""
    return RunManifest(
        run_id=run_id,
        run_name=run_name,
        started_at=datetime.now(UTC),
        status=RunStatus.INITIALIZED,
        command_line=list(command_line or sys.argv),
        working_directory=Path.cwd().resolve(),
        config_path=config_path.resolve(),
        config_sha256=sha256_file(config_path),
        config_snapshot_path=config_snapshot_path,
        input_files=input_files,
        python_version=platform_module.python_version(),
        platform=platform_module.platform(),
        package_version=__version__,
        dependency_versions=_dependency_versions(),
        git_commit=None,
        random_seed=random_seed,
    )
