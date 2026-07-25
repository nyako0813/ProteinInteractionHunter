"""Identifier normalization and auditable resolution records."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import IdentifierMatchStatus

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class IdentifierAlias(StrictModel):
    original: NonEmptyStr
    normalized: NonEmptyStr
    kind: NonEmptyStr
    source: NonEmptyStr


class IdentifierResolution(StrictModel):
    input_identifier: NonEmptyStr
    normalized_identifier: NonEmptyStr
    status: IdentifierMatchStatus
    canonical_protein_id: NonEmptyStr | None = None
    candidate_protein_ids: list[NonEmptyStr] = Field(default_factory=list)
    matched_aliases: list[IdentifierAlias] = Field(default_factory=list)
