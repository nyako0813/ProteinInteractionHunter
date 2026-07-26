from typing import Any

from protein_interaction_hunter.application.fusion import (
    FUSION_ENGINE_VERSION,
    build_fusion_index,
    evaluate_fusion_pair,
)
from protein_interaction_hunter.models.enums import EvidenceStatus
from protein_interaction_hunter.models.evidence import FusionObservation


def observation(
    *,
    query: str = "Q",
    candidate: str = "C",
    fusion: str = "F1",
    organism: str = "Org1",
    query_start: int = 1,
    query_end: int = 10,
    candidate_start: int = 20,
    candidate_end: int = 30,
    query_coverage: float | None = 0.8,
    candidate_coverage: float | None = 0.8,
    record_id: str = "r1",
) -> FusionObservation:
    return FusionObservation(
        query_protein_id=query,
        candidate_protein_id=candidate,
        fusion_protein_id=fusion,
        reference_organism=organism,
        query_component_start=query_start,
        query_component_end=query_end,
        candidate_component_start=candidate_start,
        candidate_component_end=candidate_end,
        fusion_protein_length=max(query_end, candidate_end),
        query_component_coverage=query_coverage,
        candidate_component_coverage=candidate_coverage,
        source="test",
        source_record_id=record_id,
    )


def evaluate(
    records: list[FusionObservation],
    *,
    query: str = "Q",
    candidate: str = "C",
    minimum_records: int = 1,
    coverage: float = 0.6,
    overlap: float = 0.2,
) -> Any:
    return evaluate_fusion_pair(
        query,
        candidate,
        build_fusion_index(records),
        minimum_supporting_records=minimum_records,
        minimum_component_coverage=coverage,
        maximum_component_overlap_fraction=overlap,
    )


def test_separate_domain_fusion_supports_pair() -> None:
    evidence = evaluate([observation()])
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.supporting_record_count == 1
    assert evidence.qualifying_record_count == 1
    assert evidence.pair_supported is True
    assert evidence.minimum_component_overlap_fraction == 0.0
    assert evidence.calculation_rule_version == FUSION_ENGINE_VERSION
    assert evidence.quality is None


def test_component_order_on_fusion_protein_does_not_matter() -> None:
    record = observation(query_start=20, query_end=30, candidate_start=1, candidate_end=10)
    assert evaluate([record]).pair_supported is True


def test_reversed_pair_evaluation_preserves_support_and_orients_coverage() -> None:
    record = observation(query_coverage=0.9, candidate_coverage=0.7)
    forward = evaluate([record])
    reverse = evaluate([record], query="C", candidate="Q")
    assert forward.pair_supported == reverse.pair_supported is True
    assert reverse.query_protein_id == "C"
    assert reverse.candidate_protein_id == "Q"
    assert reverse.best_query_component_coverage == 0.7
    assert reverse.best_candidate_component_coverage == 0.9


def test_coverage_threshold_is_inclusive() -> None:
    evidence = evaluate([observation(query_coverage=0.6, candidate_coverage=0.6)])
    assert evidence.pair_supported is True


def test_low_coverage_is_not_supported() -> None:
    evidence = evaluate([observation(query_coverage=0.59)])
    assert evidence.pair_supported is False
    assert "low_query_component_coverage" in evidence.conflicting_terms


def test_overlap_threshold_is_inclusive() -> None:
    record = observation(candidate_start=9, candidate_end=18)
    evidence = evaluate([record])
    assert evidence.minimum_component_overlap_fraction == 0.2
    assert evidence.pair_supported is True


def test_excessive_overlap_is_not_supported() -> None:
    record = observation(candidate_start=8, candidate_end=17)
    evidence = evaluate([record])
    assert evidence.pair_supported is False
    assert "excessive_component_overlap" in evidence.conflicting_terms
    assert "invalid_component_separation" in evidence.conflicting_terms


def test_supporting_record_threshold_is_enforced() -> None:
    evidence = evaluate([observation()], minimum_records=2)
    assert evidence.pair_supported is False
    assert "insufficient_supporting_records" in evidence.conflicting_terms


def test_multiple_fusions_and_organisms_are_aggregated() -> None:
    records = [
        observation(fusion="F2", organism="Org2", record_id="r2"),
        observation(fusion="F1", organism="Org1", record_id="r1"),
    ]
    evidence = evaluate(records, minimum_records=2)
    assert evidence.qualifying_record_count == 2
    assert evidence.fusion_protein_ids == ["F1", "F2"]
    assert evidence.reference_organisms == ["Org1", "Org2"]
    assert evidence.source_record_ids == ["r1", "r2"]
    assert "multiple_supporting_fusions" in evidence.support_terms
    assert "multiple_reference_organisms" in evidence.support_terms


def test_no_pair_record_is_available_and_false() -> None:
    evidence = evaluate([observation()], query="Q", candidate="X")
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.supporting_record_count == 0
    assert evidence.pair_supported is False
    assert evidence.conflicting_terms == ["no_fusion_record"]


def test_missing_coverage_is_unknown_not_false() -> None:
    evidence = evaluate([observation(candidate_coverage=None)])
    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.qualifying_record_count == 0
    assert evidence.pair_supported is None
    assert "missing_component_coverage" in evidence.conflicting_terms
    assert evidence.warnings == ["fusion_records_with_missing_component_coverage:1"]


def test_evaluation_is_independent_of_record_order() -> None:
    records = [
        observation(fusion="F2", organism="Org2", record_id="r2"),
        observation(fusion="F1", organism="Org1", record_id="r1"),
    ]
    assert evaluate(records, minimum_records=2) == evaluate(
        list(reversed(records)), minimum_records=2
    )
