from protein_interaction_hunter.models.enums import (
    EvidenceStatus,
    OperonProxyStatus,
    TranscriptionalOrder,
)
from protein_interaction_hunter.models.evidence import OperonEvidence


def test_operon_evidence_defaults_to_not_run() -> None:
    evidence = OperonEvidence()

    assert evidence.status is EvidenceStatus.NOT_RUN
    assert evidence.proxy_status is None
    assert evidence.support is None
    assert evidence.supporting_conditions == []
    assert evidence.conflicting_conditions == []


def test_operon_evidence_observations_do_not_require_score() -> None:
    evidence = OperonEvidence(
        status=EvidenceStatus.AVAILABLE,
        calculation_rule_version="mvp1c-operon-proxy-v1",
        same_contig=True,
        same_strand=True,
        is_adjacent=True,
        intergenic_distance_bp=29,
        overlap_bp=0,
        intervening_gene_count=0,
        transcriptional_order=TranscriptionalOrder.QUERY_THEN_CANDIDATE,
        maximum_intergenic_distance_bp=200,
        passes_distance_threshold=True,
        proxy_status=OperonProxyStatus.SUPPORTED,
        proxy_rule_id="same-strand-adjacent-short-gap-v1",
        supporting_conditions=[
            "same_contig",
            "same_strand",
            "no_intervening_gene",
            "within_intergenic_threshold",
        ],
    )

    assert evidence.proxy_status is OperonProxyStatus.SUPPORTED
    assert evidence.support is None
    assert evidence.intergenic_distance_bp == 29


def test_operon_evidence_mutable_defaults_are_independent() -> None:
    first = OperonEvidence()
    second = OperonEvidence()

    first.supporting_conditions.append("same_contig")

    assert second.supporting_conditions == []
