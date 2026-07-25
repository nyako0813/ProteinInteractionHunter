"""Normalized local domain annotation records."""

from typing import Annotated

from pydantic import Field, StringConstraints, model_validator

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import EvidenceStatus

NonEmptyStr = Annotated[
    str,
    StringConstraints(strip_whitespace=True, min_length=1),
]


class DomainAnnotationRecord(StrictModel):
    protein_id: NonEmptyStr
    source: NonEmptyStr
    accession: NonEmptyStr
    name: NonEmptyStr | None = None
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    architecture_index: int = Field(ge=0)
    status: EvidenceStatus = EvidenceStatus.AVAILABLE
    warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_coordinates(self) -> "DomainAnnotationRecord":
        if self.end < self.start:
            raise ValueError("domain end must be greater than or equal to start")
        return self