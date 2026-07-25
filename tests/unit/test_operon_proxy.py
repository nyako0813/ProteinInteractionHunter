import pytest

from protein_interaction_hunter.application.operon_proxy import (
    OPERON_PROXY_RULE_VERSION,
    calculate_operon_proxy,
)
from protein_interaction_hunter.models.enums import (
    CoordinatePosition,
    EvidenceStatus,
    OperonProxyStatus,
    RelativePosition,
    TranscriptionalOrder,
)
from protein_interaction_hunter.models.evidence import GenomeContextEvidence


def context(
    *,
    status: EvidenceStatus = EvidenceStatus.AVAILABLE,
    same_contig: bool | None = True,
    query_strand: str | None = "+",
    candidate_strand: str | None = "+",
    distance: int | None = 29,
    overlap: int | None = 0,
    intervening: int | None = 0,
    relative_position: RelativePosition | None = RelativePosition.DOWNSTREAM,
    coordinate_position: CoordinatePosition | None = CoordinatePosition.RIGHT_OF_QUERY,
) -> GenomeContextEvidence:
    return GenomeContextEvidence(
        status=status,
        same_contig=same_contig,
        same_seqid=same_contig,
        query_strand=query_strand,
        candidate_strand=candidate_strand,
        distance_bp=distance,
        edge_to_edge_distance_bp=distance,
        overlap_bp=overlap,
        intervening_gene_count=intervening,
        relative_position=relative_position,
        coordinate_position=coordinate_position,
    )


def test_supported_same_strand_adjacent_short_gap() -> None:
    evidence = calculate_operon_proxy(context(), 200)

    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.proxy_status is OperonProxyStatus.SUPPORTED
    assert evidence.same_contig is True
    assert evidence.same_strand is True
    assert evidence.is_adjacent is True
    assert evidence.passes_distance_threshold is True
    assert evidence.transcriptional_order is TranscriptionalOrder.QUERY_THEN_CANDIDATE
    assert evidence.support is None
    assert evidence.calculation_rule_version == OPERON_PROXY_RULE_VERSION
    assert evidence.conflicting_conditions == []


def test_reverse_strand_transcriptional_order() -> None:
    evidence = calculate_operon_proxy(
        context(
            query_strand="-",
            candidate_strand="-",
            relative_position=RelativePosition.UPSTREAM,
            coordinate_position=CoordinatePosition.RIGHT_OF_QUERY,
        ),
        200,
    )

    assert evidence.proxy_status is OperonProxyStatus.SUPPORTED
    assert evidence.transcriptional_order is TranscriptionalOrder.CANDIDATE_THEN_QUERY


def test_opposite_strand_is_not_supported() -> None:
    evidence = calculate_operon_proxy(context(candidate_strand="-"), 200)

    assert evidence.proxy_status is OperonProxyStatus.NOT_SUPPORTED
    assert evidence.same_strand is False
    assert "opposite_strand" in evidence.conflicting_conditions


def test_intervening_gene_is_partial_support() -> None:
    evidence = calculate_operon_proxy(context(intervening=1), 200)

    assert evidence.proxy_status is OperonProxyStatus.PARTIAL_SUPPORT
    assert evidence.is_adjacent is False
    assert "intervening_gene_present" in evidence.conflicting_conditions


def test_excessive_distance_is_partial_support() -> None:
    evidence = calculate_operon_proxy(context(distance=201), 200)

    assert evidence.proxy_status is OperonProxyStatus.PARTIAL_SUPPORT
    assert evidence.passes_distance_threshold is False
    assert "exceeds_intergenic_threshold" in evidence.conflicting_conditions


def test_overlap_is_not_supported() -> None:
    evidence = calculate_operon_proxy(
        context(
            distance=0,
            overlap=21,
            relative_position=RelativePosition.OVERLAPPING,
            coordinate_position=CoordinatePosition.OVERLAPPING,
        ),
        200,
    )

    assert evidence.proxy_status is OperonProxyStatus.NOT_SUPPORTED
    assert evidence.transcriptional_order is TranscriptionalOrder.OVERLAPPING
    assert "overlapping_features" in evidence.conflicting_conditions


def test_different_contig_is_not_applicable() -> None:
    evidence = calculate_operon_proxy(
        context(
            status=EvidenceStatus.NOT_APPLICABLE,
            same_contig=False,
            query_strand="+",
            candidate_strand="-",
            distance=None,
            overlap=None,
            intervening=None,
            relative_position=RelativePosition.DIFFERENT_CONTIG,
            coordinate_position=CoordinatePosition.DIFFERENT_CONTIG,
        ),
        200,
    )

    assert evidence.status is EvidenceStatus.NOT_APPLICABLE
    assert evidence.proxy_status is OperonProxyStatus.NOT_APPLICABLE
    assert evidence.transcriptional_order is TranscriptionalOrder.DIFFERENT_CONTIG
    assert evidence.intergenic_distance_bp is None


def test_missing_context_remains_missing() -> None:
    evidence = calculate_operon_proxy(
        context(
            status=EvidenceStatus.MISSING,
            same_contig=None,
            query_strand=None,
            candidate_strand=None,
            distance=None,
            overlap=None,
            intervening=None,
            relative_position=None,
            coordinate_position=None,
        ),
        200,
    )

    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.proxy_status is OperonProxyStatus.UNKNOWN
    assert "missing_gene_context" in evidence.conflicting_conditions


def test_same_feature_is_not_applicable() -> None:
    evidence = calculate_operon_proxy(
        context(
            distance=0,
            overlap=171,
            intervening=0,
            relative_position=RelativePosition.SAME_FEATURE,
            coordinate_position=CoordinatePosition.SAME_FEATURE,
        ),
        200,
    )

    assert evidence.status is EvidenceStatus.NOT_APPLICABLE
    assert evidence.proxy_status is OperonProxyStatus.NOT_APPLICABLE
    assert evidence.transcriptional_order is TranscriptionalOrder.SAME_FEATURE
    assert "same_feature" in evidence.conflicting_conditions


def test_unknown_strand_is_unknown() -> None:
    evidence = calculate_operon_proxy(
        context(query_strand="?", candidate_strand="+"),
        200,
    )

    assert evidence.proxy_status is OperonProxyStatus.UNKNOWN
    assert evidence.same_strand is None
    assert "unknown_strand" in evidence.conflicting_conditions


def test_negative_threshold_is_rejected() -> None:
    with pytest.raises(ValueError, match="must be non-negative"):
        calculate_operon_proxy(context(), -1)