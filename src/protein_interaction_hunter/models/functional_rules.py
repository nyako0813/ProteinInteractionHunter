"""Validated functional-complementarity rule models."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import PredictedRelationshipType

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class FunctionalRoleRule(StrictModel):
    role_id: NonEmptyStr
    include_terms: list[NonEmptyStr] = Field(min_length=1)
    exclude_terms: list[NonEmptyStr] = Field(default_factory=list)


class FunctionalPairRule(StrictModel):
    rule_id: NonEmptyStr
    query_role: NonEmptyStr
    candidate_role: NonEmptyStr
    relationship_hint: PredictedRelationshipType
    support_terms: list[NonEmptyStr] = Field(default_factory=list)
    conflicting_terms: list[NonEmptyStr] = Field(default_factory=list)


class FunctionalComplementarityRuleset(StrictModel):
    ruleset_version: NonEmptyStr
    roles: list[FunctionalRoleRule] = Field(min_length=1)
    pair_rules: list[FunctionalPairRule] = Field(min_length=1)