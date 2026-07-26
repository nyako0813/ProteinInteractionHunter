from typing import Any

from protein_interaction_hunter.application.known_interactions import (
    KNOWN_INTERACTIONS_ENGINE_VERSION,
    build_known_interaction_index,
    evaluate_known_interaction_pair,
)
from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import KnownInteractionObservation


def _observation(record_id: str, **updates: object) -> KnownInteractionObservation:
    data: dict[str, object] = {
        "status": EvidenceStatus.AVAILABLE,
        "origin": EvidenceOrigin.EXACT_PAIR,
        "protein_a_id": "QUERY",
        "protein_b_id": "CANDIDATE",
        "interaction_type": "direct",
        "reference_organism": "Test organism",
        "detection_method": "Y2H",
        "normalized_detection_method": "yeast_two_hybrid",
        "publication_id": "PMID:1",
        "confidence": 0.9,
        "is_direct": True,
        "is_physical": True,
        "is_biological": True,
        "source": "DB1",
        "source_record_id": record_id,
    }
    data.update(updates)
    return KnownInteractionObservation.model_validate(data)


def _evaluate(
    records: list[KnownInteractionObservation],
    **updates: object,
) -> Any:
    options: Any = {
        "minimum_supporting_records": 1,
        "minimum_direct_records": 1,
        "accepted_interaction_types": [
            "physical",
            "direct",
            "genetic",
            "functional_association",
        ],
        "accepted_evidence_methods": [],
        "excluded_evidence_methods": ["database inference"],
        "minimum_confidence": 0.5,
    }
    options.update(updates)
    return evaluate_known_interaction_pair(
        "QUERY", "CANDIDATE", build_known_interaction_index(records), **options
    )


def test_multiple_direct_physical_records_and_independence_counts() -> None:
    records = [
        _observation("R1"),
        _observation(
            "R2",
            protein_a_id="CANDIDATE",
            protein_b_id="QUERY",
            interaction_type="physical",
            detection_method="SPR",
            normalized_detection_method="surface_plasmon_resonance",
            source="DB2",
            publication_id="PMID:1",
            is_direct=True,
        ),
    ]
    evidence = _evaluate(records)
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.pair_supported is True
    assert evidence.direct_record_count == 2
    assert evidence.physical_record_count == 2
    assert evidence.independent_publication_count == 1
    assert evidence.independent_source_count == 2
    assert evidence.calculation_rule_version == KNOWN_INTERACTIONS_ENGINE_VERSION


def test_pair_orientation_is_canonical() -> None:
    record = _observation("R1", protein_a_id="CANDIDATE", protein_b_id="QUERY")
    assert _evaluate([record]).pair_supported is True


def test_functional_only_is_separate_from_direct_and_physical() -> None:
    record = _observation(
        "R1",
        interaction_type="functional_association",
        is_direct=False,
        is_physical=False,
        is_biological=True,
    )
    evidence = _evaluate([record])
    assert evidence.pair_supported is False
    assert evidence.direct_interaction_supported is False
    assert evidence.physical_interaction_supported is False
    assert evidence.functional_association_supported is True
    assert "functional_association_only" in evidence.conflicting_terms


def test_missing_confidence_is_missing_not_false() -> None:
    evidence = _evaluate([_observation("R1", confidence=None)])
    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.pair_supported is None
    assert evidence.direct_interaction_supported is None
    assert "missing_confidence" in evidence.conflicting_terms


def test_no_record_is_available_false() -> None:
    evidence = _evaluate([])
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.pair_supported is False
    assert evidence.supporting_record_count == 0
    assert evidence.direct_interaction_supported is False


def test_excluded_method_is_known_false() -> None:
    evidence = _evaluate([_observation("R1", normalized_detection_method="database_inference")])
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.pair_supported is False
    assert "excluded_detection_method" in evidence.conflicting_terms


def test_low_confidence_is_known_false() -> None:
    evidence = _evaluate([_observation("R1", confidence=0.2)])
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.pair_supported is False
    assert "low_confidence_record" in evidence.conflicting_terms


def test_predicted_only_is_not_accepted_by_default_policy() -> None:
    evidence = _evaluate(
        [_observation("R1", interaction_type="predicted", is_direct=False, is_physical=False)]
    )
    assert evidence.pair_supported is False
    assert "unsupported_interaction_type" in evidence.conflicting_terms


def test_method_alone_never_makes_record_direct() -> None:
    record = _observation(
        "R1",
        interaction_type="co_complex",
        is_direct=None,
        is_physical=True,
        normalized_detection_method="yeast_two_hybrid",
    )
    evidence = _evaluate(
        [record],
        minimum_direct_records=0,
        accepted_interaction_types=["co_complex"],
    )
    assert evidence.pair_supported is True
    assert evidence.direct_record_count == 0
    assert evidence.direct_interaction_supported is False
    assert evidence.physical_interaction_supported is True


def test_missing_required_accepted_method_is_unknown() -> None:
    record = _observation("R1", detection_method=None, normalized_detection_method=None)
    evidence = _evaluate([record], accepted_evidence_methods=["Y2H"])
    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.pair_supported is None


def test_minimum_direct_zero_allows_functional_pair_support() -> None:
    record = _observation("R1", interaction_type="genetic", is_direct=False, is_physical=False)
    evidence = _evaluate([record], minimum_direct_records=0)
    assert evidence.pair_supported is True
    assert evidence.functional_association_supported is True
