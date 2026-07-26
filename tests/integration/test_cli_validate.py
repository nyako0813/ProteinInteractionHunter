"""Typer CLI smoke and validation tests."""

from pathlib import Path

from typer.testing import CliRunner

from protein_interaction_hunter.cli import app

runner = CliRunner()


def test_cli_version() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == "0.1.0"


def test_cli_validate_config(valid_config_path: Path) -> None:
    result = runner.invoke(app, ["validate-config", "--config", str(valid_config_path)])
    assert result.exit_code == 0
    assert "Configuration valid" in result.stdout
    assert "queries=1" in result.stdout


def test_cli_validate_inputs(valid_config_path: Path) -> None:
    result = runner.invoke(app, ["validate-inputs", "--config", str(valid_config_path)])
    assert result.exit_code == 0
    assert "proteins=13" in result.stdout
    assert "duplicates=1" in result.stdout
    assert "missing_coordinates=2" in result.stdout


def test_cli_inspect_fixture_has_counts_but_no_ranking(
    valid_config_path: Path,
) -> None:
    result = runner.invoke(app, ["inspect-fixture", "--config", str(valid_config_path)])
    assert result.exit_code == 0
    assert "protein_count: 13" in result.stdout
    assert "hypothetical_protein_count: 1" in result.stdout
    assert "score" not in result.stdout.lower()
    assert "rank" not in result.stdout.lower()


def test_cli_invalid_config_reports_field_location(fixture_dir: Path) -> None:
    result = runner.invoke(
        app,
        ["validate-config", "--config", str(fixture_dir / "config.invalid.yaml")],
    )
    assert result.exit_code == 1
    assert "minimum_length_aa" in result.stderr
