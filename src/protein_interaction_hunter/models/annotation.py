"""Normalized local annotation record."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import EvidenceStatus

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class AnnotationRecord(StrictModel):
    protein_id: NonEmptyStr
    gene_name: str | None = None
    locus_tag: str | None = None
    product: str | None = None
    functional_category: str | None = None
    localization_annotation: str | None = None
    transmembrane_annotation: str | None = None
    annotation_source: str | None = None
    annotation_confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    status: EvidenceStatus
    warnings: list[str] = Field(default_factory=list)
