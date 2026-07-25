"""Protein, query, and candidate identity models."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceStatus,
    IdentifierMatchStatus,
)

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class ProteinRecord(StrictModel):
    protein_id: NonEmptyStr
    description: str = ""
    sequence: NonEmptyStr
    gene_id: NonEmptyStr | None = None
    locus_tag: NonEmptyStr | None = None
    aliases: list[NonEmptyStr] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


class QueryProtein(StrictModel):
    query_id: NonEmptyStr
    protein_id: NonEmptyStr
    resolution_method: NonEmptyStr = "exact_id"
    warnings: list[str] = Field(default_factory=list)


class CandidateProtein(StrictModel):
    query_id: NonEmptyStr
    protein_id: NonEmptyStr
    disposition: CandidateDisposition = CandidateDisposition.INCLUDED
    disposition_reasons: list[str] = Field(default_factory=list)
    duplicate_sequence_group: NonEmptyStr | None = None
    paralog_group: NonEmptyStr | None = None
    sequence_length: int = Field(ge=1)
    description: str = ""
    gene_id: NonEmptyStr | None = None
    locus_tag: NonEmptyStr | None = None
    old_locus_tag: NonEmptyStr | None = None
    contig: NonEmptyStr | None = None
    strand: str | None = Field(default=None, pattern=r"^[+\-?]$")
    has_coordinate: bool = False
    has_annotation: bool = False
    coordinate_status: EvidenceStatus = EvidenceStatus.MISSING
    annotation_status: EvidenceStatus = EvidenceStatus.MISSING
    same_contig_as_query: bool | None = None
    is_duplicate_sequence: bool = False
    is_fragment_candidate: bool = False
    fragment_reasons: list[str] = Field(default_factory=list)
    is_hypothetical: bool = False
    identifier_match_status: IdentifierMatchStatus = IdentifierMatchStatus.NO_MATCH
    original_identifiers: dict[str, list[str]] = Field(default_factory=dict)
    normalized_identifiers: dict[str, list[str]] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
