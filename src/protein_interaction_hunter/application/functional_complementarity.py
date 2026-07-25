"""Deterministic rule-based functional-complementarity evaluation."""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path

from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    FunctionalEvidence,
)
from protein_interaction_hunter.models.functional_rules import (
    FunctionalComplementarityRuleset,
    FunctionalRoleRule,
)

FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION = (
    "mvp1d-functional-complementarity-v1"
)


def normalize_annotation_text(value: str | None) -> str:
    """Normalize annotation text for deterministic term matching."""

    normalized = unicodedata.normalize("NFKC", value or "").casefold()
    normalized = re.sub(r"[^0-9a-z]+", " ", normalized)
    return " ".join(normalized.split())


def _term_is_present(text: str, term: str) -> bool:
    normalized_term = normalize_annotation_text(term)
    if not text or not normalized_term:
        return False
    return f" {normalized_term} " in f" {text} "


def _matched_terms(
    text: str,
    terms: list[str],
) -> list[str]:
    return sorted(
        {
            term
            for term in terms
            if _term_is_present(text, term)
        }
    )


def _role_matches(
    text: str,
    role: FunctionalRoleRule,
) -> tuple[bool, list[str], list[str]]:
    included = _matched_terms(text, role.include_terms)
    excluded = _matched_terms(text, role.exclude_terms)

    return bool(included) and not excluded, included, excluded


def evaluate_functional_complementarity(
    query_annotation_text: str | None,
    candidate_annotation_text: str | None,
    ruleset: FunctionalComplementarityRuleset,
    ruleset_path: Path | None = None,
) -> list[FunctionalEvidence]:
    """Evaluate functional pair rules without scoring or ranking."""

    query_text = normalize_annotation_text(query_annotation_text)
    candidate_text = normalize_annotation_text(candidate_annotation_text)
    serialized_path = str(ruleset_path.resolve()) if ruleset_path else None

    provenance = [
        EvidenceProvenance(
            source_name="functional_complementarity_rules",
            source_version=ruleset.ruleset_version,
            method="deterministic_normalized_term_matching",
            metadata={
                "engine_version": FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION,
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]

    if not query_text or not candidate_text:
        missing_entities: list[str] = []
        if not query_text:
            missing_entities.append("query_annotation")
        if not candidate_text:
            missing_entities.append("candidate_annotation")

        return [
            FunctionalEvidence(
                status=EvidenceStatus.MISSING,
                origin=EvidenceOrigin.ANNOTATION,
                calculation_rule_version=(
                    FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION
                ),
                query_annotation_text=query_text or None,
                candidate_annotation_text=candidate_text or None,
                ruleset_path=serialized_path,
                matched=None,
                conflicting_terms=missing_entities,
                provenance=provenance,
                warnings=[
                    "missing_functional_annotation:"
                    + ",".join(missing_entities)
                ],
            )
        ]

    query_roles: dict[str, tuple[list[str], list[str]]] = {}
    candidate_roles: dict[str, tuple[list[str], list[str]]] = {}

    for role in ruleset.roles:
        query_matches, query_included, query_excluded = _role_matches(
            query_text,
            role,
        )
        if query_matches:
            query_roles[role.role_id] = (
                query_included,
                query_excluded,
            )

        candidate_matches, candidate_included, candidate_excluded = (
            _role_matches(candidate_text, role)
        )
        if candidate_matches:
            candidate_roles[role.role_id] = (
                candidate_included,
                candidate_excluded,
            )

    evidence_records: list[FunctionalEvidence] = []

    for pair_rule in ruleset.pair_rules:
        if pair_rule.query_role not in query_roles:
            continue
        if pair_rule.candidate_role not in candidate_roles:
            continue

        query_terms, query_conflicts = query_roles[
            pair_rule.query_role
        ]
        candidate_terms, candidate_conflicts = candidate_roles[
            pair_rule.candidate_role
        ]

        evidence_records.append(
            FunctionalEvidence(
                status=EvidenceStatus.AVAILABLE,
                origin=EvidenceOrigin.ANNOTATION,
                calculation_rule_version=(
                    FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION
                ),
                query_role=pair_rule.query_role,
                candidate_role=pair_rule.candidate_role,
                relationship_hint=pair_rule.relationship_hint,
                rule_id=pair_rule.rule_id,
                query_matched_terms=query_terms,
                candidate_matched_terms=candidate_terms,
                support_terms=sorted(
                    set(
                        pair_rule.support_terms
                        + [f"query:{term}" for term in query_terms]
                        + [
                            f"candidate:{term}"
                            for term in candidate_terms
                        ]
                    )
                ),
                conflicting_terms=sorted(
                    set(
                        pair_rule.conflicting_terms
                        + query_conflicts
                        + candidate_conflicts
                    )
                ),
                query_annotation_text=query_text,
                candidate_annotation_text=candidate_text,
                ruleset_path=serialized_path,
                matched=True,
                provenance=provenance,
            )
        )

    if evidence_records:
        return evidence_records

    return [
        FunctionalEvidence(
            status=EvidenceStatus.AVAILABLE,
            origin=EvidenceOrigin.ANNOTATION,
            calculation_rule_version=(
                FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION
            ),
            query_annotation_text=query_text,
            candidate_annotation_text=candidate_text,
            ruleset_path=serialized_path,
            matched=False,
            support_terms=[],
            conflicting_terms=["no_matching_pair_rule"],
            provenance=provenance,
        )
    ]