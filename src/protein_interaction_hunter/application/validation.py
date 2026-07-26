"""Local input validation and fixture summaries; no ranking is performed."""

from pydantic import Field

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.domains import LocalDomainTsvLoader
from protein_interaction_hunter.adapters.local.fasta import (
    LocalFastaLoader,
    duplicate_sequence_groups,
)
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.application.gene_context import build_coordinate_index
from protein_interaction_hunter.config import AppConfig
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.base import StrictModel


class InputValidationSummary(StrictModel):
    protein_count: int = Field(ge=0)
    query_count: int = Field(ge=0)
    gff_feature_count: int = Field(ge=0)
    gff_coordinate_count: int = Field(ge=0)
    annotation_count: int = Field(ge=0)
    domain_annotation_count: int = Field(ge=0)
    domain_protein_count: int = Field(ge=0)
    unknown_domain_id_count: int = Field(ge=0)
    query_domain_annotation_count: int = Field(ge=0)
    identifier_match_count: int = Field(ge=0)
    duplicate_sequence_group_count: int = Field(ge=0)
    missing_coordinate_count: int = Field(ge=0)
    hypothetical_protein_count: int = Field(ge=0)
    warnings: list[str] = Field(default_factory=list)


def validate_local_inputs(config: AppConfig) -> InputValidationSummary:
    """Validate required local files and summarize identifiers without analysis."""
    proteins = LocalFastaLoader().load(config.input.proteome_fasta)
    protein_ids = {protein.protein_id for protein in proteins}
    missing_queries = sorted(set(config.query.protein_ids) - protein_ids)
    if missing_queries:
        raise InputValidationError(
            "Query IDs not found in proteome FASTA: " + ", ".join(missing_queries)
        )

    document = LocalGff3Loader().load_document(config.input.genome_gff)
    coordinate_index = build_coordinate_index(proteins, document)
    matched_ids = protein_ids & set(coordinate_index.by_protein)
    missing_coordinates = protein_ids - set(coordinate_index.by_protein)

    annotations = []
    if config.input.annotation_table is not None:
        annotations = LocalAnnotationTsvLoader().load(config.input.annotation_table)

    domain_records = []
    if config.domains.local_table is not None:
        domain_records = LocalDomainTsvLoader().load(config.domains.local_table)
    domain_ids = {record.protein_id for record in domain_records}

    annotation_by_id = {record.protein_id: record for record in annotations}
    hypothetical_count = 0
    for protein in proteins:
        annotation = annotation_by_id.get(protein.protein_id)
        searchable = " ".join(
            part
            for part in (
                protein.description,
                annotation.product if annotation is not None else None,
            )
            if part
        ).lower()
        if "hypothetical" in searchable or "uncharacterized" in searchable:
            hypothetical_count += 1

    warnings: list[str] = []
    unknown_annotation_ids = sorted(set(annotation_by_id) - protein_ids)
    if unknown_annotation_ids:
        warnings.append("Annotation IDs absent from proteome: " + ", ".join(unknown_annotation_ids))
    unknown_domain_ids = sorted(domain_ids - protein_ids)
    if unknown_domain_ids:
        warnings.append("Domain IDs absent from proteome: " + ", ".join(unknown_domain_ids))

    return InputValidationSummary(
        protein_count=len(proteins),
        query_count=len(config.query.protein_ids),
        gff_feature_count=len(document.features),
        gff_coordinate_count=len(coordinate_index.by_protein),
        annotation_count=len(annotations),
        domain_annotation_count=len(domain_records),
        domain_protein_count=len(domain_ids & protein_ids),
        unknown_domain_id_count=len(domain_ids - protein_ids),
        query_domain_annotation_count=sum(
            record.protein_id in config.query.protein_ids for record in domain_records
        ),
        identifier_match_count=len(matched_ids),
        duplicate_sequence_group_count=len(duplicate_sequence_groups(proteins)),
        missing_coordinate_count=len(missing_coordinates),
        hypothetical_protein_count=hypothetical_count,
        warnings=warnings,
    )
