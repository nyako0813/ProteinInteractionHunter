"""Genome coordinate models."""

from typing import Annotated

from pydantic import Field, StringConstraints, model_validator

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import IdentifierMatchStatus

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class GeneCoordinate(StrictModel):
    seqid: NonEmptyStr
    feature_type: NonEmptyStr
    source: NonEmptyStr | None = None
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    strand: str | None = Field(default=None, pattern=r"^[+\-?]$")
    feature_id: NonEmptyStr | None = None
    parent_id: NonEmptyStr | None = None
    parent_ids: list[str] = Field(default_factory=list)
    protein_id: NonEmptyStr | None = None
    locus_tag: NonEmptyStr | None = None
    old_locus_tag: NonEmptyStr | None = None
    attributes: dict[str, list[str]] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_interval(self) -> "GeneCoordinate":
        if self.end < self.start:
            raise ValueError("end must be greater than or equal to start")
        return self


class SequenceRegion(StrictModel):
    seqid: NonEmptyStr
    start: int = Field(ge=1)
    end: int = Field(ge=1)

    @model_validator(mode="after")
    def validate_interval(self) -> "SequenceRegion":
        if self.end < self.start:
            raise ValueError("sequence-region end must be greater than or equal to start")
        return self


class GffDocument(StrictModel):
    features: list[GeneCoordinate]
    sequence_regions: dict[str, SequenceRegion] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)


class NormalizedFeature(StrictModel):
    representative_id: NonEmptyStr
    protein_id: NonEmptyStr | None = None
    gene_id: NonEmptyStr | None = None
    locus_tag: NonEmptyStr | None = None
    old_locus_tag: NonEmptyStr | None = None
    seqid: NonEmptyStr
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    strand: str = Field(default="?", pattern=r"^[+\-?]$")
    feature_type: NonEmptyStr
    coordinate_source: NonEmptyStr | None = None
    identifier_status: IdentifierMatchStatus
    parent_identifiers: list[str] = Field(default_factory=list)
    source_feature_ids: list[str] = Field(default_factory=list)
    is_gene: bool = False
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_interval(self) -> "NormalizedFeature":
        if self.end < self.start:
            raise ValueError("end must be greater than or equal to start")
        return self
