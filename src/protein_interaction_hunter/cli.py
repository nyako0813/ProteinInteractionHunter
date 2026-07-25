"""Validation-only MVP-0 command line interface."""

from pathlib import Path
from typing import Annotated

import typer

from protein_interaction_hunter import __version__
from protein_interaction_hunter.application.validation import validate_local_inputs
from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import ProteinInteractionHunterError

app = typer.Typer(
    name="protein-interaction-hunter",
    help="Validate ProteinInteractionHunter MVP-0 configuration and local fixtures.",
    no_args_is_help=True,
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option("--version", callback=_version_callback, is_eager=True),
    ] = None,
) -> None:
    """MVP-0 validation commands; biological analysis is not implemented."""


def _fail(exc: Exception) -> None:
    typer.echo(f"Error: {exc}", err=True)
    raise typer.Exit(code=1)


@app.command("validate-config")
def validate_config_command(
    config: Annotated[Path, typer.Option("--config", exists=True, dir_okay=False)],
) -> None:
    """Validate YAML structure without requiring optional inputs to exist."""
    try:
        loaded = load_config(config)
    except ProteinInteractionHunterError as exc:
        _fail(exc)
        return
    typer.echo(
        f"Configuration valid: run={loaded.project.run_name}, "
        f"queries={len(loaded.query.protein_ids)}, local_only={loaded.project.local_only}"
    )


@app.command("validate-inputs")
def validate_inputs_command(
    config: Annotated[Path, typer.Option("--config", exists=True, dir_okay=False)],
) -> None:
    """Validate required local inputs and identifier correspondence."""
    try:
        summary = validate_local_inputs(load_config(config))
    except ProteinInteractionHunterError as exc:
        _fail(exc)
        return
    typer.echo(
        "Inputs valid: "
        f"proteins={summary.protein_count}, coordinates={summary.gff_coordinate_count}, "
        f"identifier_matches={summary.identifier_match_count}, "
        f"duplicates={summary.duplicate_sequence_group_count}, "
        f"missing_coordinates={summary.missing_coordinate_count}"
    )


@app.command("inspect-fixture")
def inspect_fixture_command(
    config: Annotated[Path, typer.Option("--config", exists=True, dir_okay=False)],
) -> None:
    """Print fixture counts only; no scores or ranks are calculated."""
    try:
        summary = validate_local_inputs(load_config(config))
    except ProteinInteractionHunterError as exc:
        _fail(exc)
        return
    for label, value in (
        ("protein_count", summary.protein_count),
        ("query_count", summary.query_count),
        ("gff_coordinate_count", summary.gff_coordinate_count),
        ("annotation_count", summary.annotation_count),
        ("duplicate_sequence_group_count", summary.duplicate_sequence_group_count),
        ("missing_coordinate_count", summary.missing_coordinate_count),
        ("hypothetical_protein_count", summary.hypothetical_protein_count),
    ):
        typer.echo(f"{label}: {value}")
