from protein_interaction_hunter.models.enums import (
    EvidenceStatus,
    PredictedRelationshipType,
)
from protein_interaction_hunter.models.evidence import FunctionalEvidence


def test_functional_evidence_defaults_to_not_run() -> None:
    evidence = FunctionalEvidence()

    assert evidence.status is EvidenceStatus.NOT_RUN
    assert evidence.matched is None
    assert evidence.relationship_hint is None
    assert evidence.support_terms == []
    assert evidence.conflicting_terms == []


def test_functional_evidence_records_rule_match_without_score() -> None:
    evidence = FunctionalEvidence(
        status=EvidenceStatus.AVAILABLE,
        calculation_rule_version="mvp1d-functional-complementarity-v1",
        query_role="enzyme",
        candidate_role="accessory_factor",
        relationship_hint=PredictedRelationshipType.ACCESSORY_FACTOR,
        rule_id="enzyme-accessory-factor-v1",
        query_matched_terms=["enzyme"],
        candidate_matched_terms=["accessory"],
        support_terms=["query:enzyme", "candidate:accessory"],
        query_annotation_text="query enzyme",
        candidate_annotation_text="nearby accessory candidate",
        ruleset_path="rules/functional_complementarity.v1.yaml",
        matched=True,
    )

    assert evidence.matched is True
    assert evidence.relationship_hint is PredictedRelationshipType.ACCESSORY_FACTOR
    assert evidence.quality is None


def test_functional_evidence_list_defaults_are_independent() -> None:
    first = FunctionalEvidence()
    second = FunctionalEvidence()

    first.support_terms.append("enzyme")

    assert second.support_terms == []