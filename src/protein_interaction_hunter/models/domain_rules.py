"""Validated domain-pair rule models."""

from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel

NonEmptyStr = Annotated[
    str,
    StringConstraints(strip_whitespace=True, min_length=1),
]


class DomainRoleRule(StrictModel):
    role_id: NonEmptyStr
    accessions: list[NonEmptyStr] = Field(min_length=1)


class DomainPairRule(StrictModel):
    rule_id: NonEmptyStr
    query_role: NonEmptyStr
    candidate_role: NonEmptyStr
    support_terms: list[NonEmptyStr] = Field(default_factory=list)
    conflicting_terms: list[NonEmptyStr] = Field(default_factory=list)
    allow_shared_accession: bool = False


class DomainPairRuleset(StrictModel):
    ruleset_version: NonEmptyStr
    roles: list[DomainRoleRule] = Field(min_length=1)
    pair_rules: list[DomainPairRule] = Field(min_length=1)