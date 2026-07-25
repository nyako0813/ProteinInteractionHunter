"""Core model semantics."""

import json

import pytest
from pydantic import ValidationError

from protein_interaction_hunter.models import (
    CandidateDisposition,
    CandidateEvidenceBundle,
    CandidateScore,
    EvidenceOrigin,
    EvidenceStatus,
    PredictedRelationshipType,
    ProteinRecord,
)
from protein_interaction_hunter.models.evidence import DomainEvidence
from protein_interaction_hunter.outputs.jsonl import deterministic_json


def make_bundle() -> CandidateEvidenceBundle:
    return CandidateEvidenceBundle(
        run_id="run-001",
        query_id="QUERY_001",
        candidate_id="NEAR_001",
        candidate_disposition=CandidateDisposition.INCLUDED,
        predicted_relationship_type=PredictedRelationshipType.INSUFFICIENT_EVIDENCE,
        engine_statuses={"gene_context": EvidenceStatus.NOT_RUN},
    )


def test_score_range_is_validated() -> None:
    with pytest.raises(ValidationError):
        CandidateScore(physical_interaction_score=1.01)
    with pytest.raises(ValidationError):
        CandidateScore(contradiction_penalty=-0.01)


def test_none_and_zero_scores_remain_distinct() -> None:
    score = CandidateScore(physical_interaction_score=0.0)
    assert score.physical_interaction_score == 0.0
    assert score.functional_association_score is None
    assert score.total_ranking_score is None


def test_mutable_defaults_are_not_shared() -> None:
    left = make_bundle()
    right = make_bundle().model_copy(update={"candidate_id": "OTHER_001"}, deep=True)
    left.warnings.append("left only")
    assert right.warnings == []


def test_empty_identifier_is_rejected() -> None:
    with pytest.raises(ValidationError):
        ProteinRecord(protein_id="", sequence="MSTK")


def test_enums_serialize_to_stable_values() -> None:
    bundle = make_bundle()
    payload = json.loads(deterministic_json(bundle))
    assert payload["predicted_relationship_type"] == "insufficient_evidence"
    assert payload["candidate_disposition"] == "included"
    assert payload["engine_statuses"]["gene_context"] == "not_run"


def test_exact_and_ortholog_transferred_origins_are_distinct() -> None:
    exact = DomainEvidence(
        protein_id="P1",
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.EXACT_PROTEIN,
    )
    transferred = exact.model_copy(update={"origin": EvidenceOrigin.ORTHOLOG_TRANSFERRED})
    assert exact.origin != transferred.origin


def test_deterministic_serialization() -> None:
    bundle = make_bundle()
    assert deterministic_json(bundle) == deterministic_json(bundle.model_copy(deep=True))
