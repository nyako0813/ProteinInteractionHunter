"""Protein, query, and candidate identity models."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import CandidateDisposition

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
    warnings: list[str] = Field(default_factory=list)
