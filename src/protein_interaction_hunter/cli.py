"""Local validation and MVP-1 evidence-prioritization command line interface."""

from pathlib import Path
from typing import Annotated

import typer

from protein_interaction_hunter import __version__
from protein_interaction_hunter.application.pipeline import InteractionCandidatePipeline
from protein_interaction_hunter.application.validation import validate_local_inputs
from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import ProteinInteractionHunterError

app = typer.Typer(
    name="protein-interaction-hunter",
    help="Validate local inputs and generate auditable MVP-1 evidence, scores, ranks, and tiers.",
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
    """Validate inputs and run the local-only evidence-prioritization pipeline."""


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
        f"missing_coordinates={summary.missing_coordinate_count}, "
        f"domains={summary.domain_annotation_count}, "
        f"domain_proteins={summary.domain_protein_count}, "
        f"unknown_domain_ids={summary.unknown_domain_id_count}, "
        f"query_domains={summary.query_domain_annotation_count}, "
        f"orthology_records={summary.orthology_annotation_count}, "
        f"orthology_proteins={summary.orthology_protein_count}, "
        f"unknown_orthology_ids={summary.unknown_orthology_id_count}, "
        f"profile_observations={summary.phylogenetic_profile_observation_count}, "
        f"profile_proteins={summary.phylogenetic_profile_protein_count}, "
        f"profile_species={summary.phylogenetic_profile_species_count}, "
        f"unknown_profile_ids={summary.unknown_phylogenetic_profile_id_count}"
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
        ("domain_annotation_count", summary.domain_annotation_count),
        ("domain_protein_count", summary.domain_protein_count),
        ("duplicate_sequence_group_count", summary.duplicate_sequence_group_count),
        ("missing_coordinate_count", summary.missing_coordinate_count),
        ("hypothetical_protein_count", summary.hypothetical_protein_count),
    ):
        typer.echo(f"{label}: {value}")


@app.command("generate-candidates")
def generate_candidates_command(
    config: Annotated[Path, typer.Option("--config", exists=True, dir_okay=False)],
) -> None:
    """Generate all query-protein pairs with optional evidence, scoring, and tiers."""
    try:
        result = InteractionCandidatePipeline().run(
            config,
            command_line=[
                "protein-interaction-hunter",
                "generate-candidates",
                "--config",
                str(config),
            ],
        )
    except ProteinInteractionHunterError as exc:
        _fail(exc)
        return
    summary = result.summary
    for label, value in (
        ("query_count", summary.query_count),
        ("protein_count", summary.protein_count),
        ("query_candidate_pair_count", summary.pair_count),
        ("included_count", summary.included_count),
        ("flagged_count", summary.flagged_count),
        ("excluded_count", summary.excluded_count),
        ("duplicate_group_count", summary.duplicate_group_count),
        ("fragment_candidate_count", summary.fragment_candidate_count),
        ("hypothetical_protein_count", summary.hypothetical_protein_count),
        ("missing_coordinate_count", summary.missing_coordinate_count),
        ("missing_annotation_count", summary.missing_annotation_count),
        ("ambiguous_mapping_count", summary.ambiguous_mapping_count),
        ("same_contig_pair_count", summary.same_contig_pair_count),
        ("different_contig_pair_count", summary.different_contig_pair_count),
        ("overlapping_pair_count", summary.overlapping_pair_count),
        ("missing_context_pair_count", summary.missing_context_pair_count),
        ("ambiguous_context_pair_count", summary.ambiguous_context_pair_count),
        ("neighborhood_pair_count", summary.neighborhood_pair_count),
        ("incomplete_context_pair_count", summary.incomplete_context_pair_count),
        ("output_path", summary.output_path),
    ):
        typer.echo(f"{label}: {value}")
