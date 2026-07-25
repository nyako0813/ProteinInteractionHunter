"""Deterministic operon-proxy interpretation of observed gene context.

This module does not claim that two genes form an operon. It records whether
simple coordinate-derived conditions are compatible with an operon-like
arrangement. The result is not propagated to ranking, scores, or evidence tiers.
"""

from __future__ import annotations

from protein_interaction_hunter.models.enums import (
    CoordinatePosition,
    EvidenceOrigin,
    EvidenceStatus,
    OperonProxyStatus,
    RelativePosition,
    TranscriptionalOrder,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    GenomeContextEvidence,
    OperonEvidence,
)

OPERON_PROXY_RULE_VERSION = "mvp1c-operon-proxy-v1"
OPERON_PROXY_RULE_ID = "same-contig-same-strand-adjacency-distance-v1"


def _transcriptional_order(context: GenomeContextEvidence) -> TranscriptionalOrder:
    """Describe pair order in the shared transcription direction."""

    if context.same_contig is False:
        return TranscriptionalOrder.DIFFERENT_CONTIG

    if context.coordinate_position is CoordinatePosition.SAME_FEATURE:
        return TranscriptionalOrder.SAME_FEATURE

    if context.coordinate_position is CoordinatePosition.OVERLAPPING:
        return TranscriptionalOrder.OVERLAPPING

    if context.relative_position is RelativePosition.DOWNSTREAM:
        return TranscriptionalOrder.QUERY_THEN_CANDIDATE

    if context.relative_position is RelativePosition.UPSTREAM:
        return TranscriptionalOrder.CANDIDATE_THEN_QUERY

    return TranscriptionalOrder.UNKNOWN


def _base_provenance(maximum_intergenic_distance_bp: int) -> list[EvidenceProvenance]:
    return [
        EvidenceProvenance(
            source_name="gene_context_operon_proxy",
            source_version=OPERON_PROXY_RULE_VERSION,
            method="deterministic_coordinate_rule_evaluation",
            metadata={
                "maximum_intergenic_distance_bp": maximum_intergenic_distance_bp,
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]


def calculate_operon_proxy(
    context: GenomeContextEvidence,
    maximum_intergenic_distance_bp: int,
) -> OperonEvidence:
    """Calculate an auditable operon proxy from gene-context observations.

    Rules:
    - Missing or failed gene context remains missing/failed.
    - Different contigs are not applicable.
    - The same feature is not treated as an operon pair.
    - A supported proxy requires:
      same contig, known same strand, no intervening gene, no overlap,
      and edge-to-edge distance within the configured threshold.
    - Overlap, opposite/unknown strand, intervening genes, and excessive
      distance are recorded transparently as conflicting conditions.
    - No numerical score is produced.
    """

    if maximum_intergenic_distance_bp < 0:
        raise ValueError("maximum_intergenic_distance_bp must be non-negative")

    provenance = _base_provenance(maximum_intergenic_distance_bp)

    if context.status is EvidenceStatus.NOT_RUN:
        return OperonEvidence(
            status=EvidenceStatus.NOT_RUN,
            origin=EvidenceOrigin.INFERRED,
            calculation_rule_version=OPERON_PROXY_RULE_VERSION,
            maximum_intergenic_distance_bp=maximum_intergenic_distance_bp,
            proxy_status=OperonProxyStatus.UNKNOWN,
            proxy_rule_id=OPERON_PROXY_RULE_ID,
            provenance=provenance,
            warnings=["gene_context_not_run"],
        )

    if context.status is EvidenceStatus.FAILED:
        return OperonEvidence(
            status=EvidenceStatus.FAILED,
            origin=EvidenceOrigin.INFERRED,
            calculation_rule_version=OPERON_PROXY_RULE_VERSION,
            maximum_intergenic_distance_bp=maximum_intergenic_distance_bp,
            proxy_status=OperonProxyStatus.UNKNOWN,
            proxy_rule_id=OPERON_PROXY_RULE_ID,
            provenance=provenance,
            warnings=sorted(set(context.warnings + ["gene_context_failed"])),
        )

    if context.status is EvidenceStatus.MISSING:
        return OperonEvidence(
            status=EvidenceStatus.MISSING,
            origin=EvidenceOrigin.INFERRED,
            calculation_rule_version=OPERON_PROXY_RULE_VERSION,
            same_contig=context.same_contig,
            same_strand=None,
            intergenic_distance_bp=None,
            overlap_bp=context.overlap_bp,
            intervening_gene_count=context.intervening_gene_count,
            transcriptional_order=_transcriptional_order(context),
            maximum_intergenic_distance_bp=maximum_intergenic_distance_bp,
            passes_distance_threshold=None,
            proxy_status=OperonProxyStatus.UNKNOWN,
            proxy_rule_id=OPERON_PROXY_RULE_ID,
            conflicting_conditions=["missing_gene_context"],
            provenance=provenance,
            warnings=sorted(set(context.warnings)),
        )

    if context.same_contig is False or context.status is EvidenceStatus.NOT_APPLICABLE:
        return OperonEvidence(
            status=EvidenceStatus.NOT_APPLICABLE,
            origin=EvidenceOrigin.INFERRED,
            calculation_rule_version=OPERON_PROXY_RULE_VERSION,
            same_contig=False,
            same_strand=None,
            is_adjacent=None,
            intergenic_distance_bp=None,
            overlap_bp=None,
            intervening_gene_count=None,
            transcriptional_order=TranscriptionalOrder.DIFFERENT_CONTIG,
            maximum_intergenic_distance_bp=maximum_intergenic_distance_bp,
            passes_distance_threshold=None,
            proxy_status=OperonProxyStatus.NOT_APPLICABLE,
            proxy_rule_id=OPERON_PROXY_RULE_ID,
            conflicting_conditions=["different_contig"],
            provenance=provenance,
            warnings=sorted(set(context.warnings)),
        )

    same_strand: bool | None
    if context.query_strand not in {"+", "-"} or context.candidate_strand not in {"+", "-"}:
        same_strand = None
    else:
        same_strand = context.query_strand == context.candidate_strand

    intervening_gene_count = context.intervening_gene_count
    is_adjacent = intervening_gene_count == 0 if intervening_gene_count is not None else None

    overlap_bp = context.overlap_bp
    is_overlapping = overlap_bp is not None and overlap_bp > 0

    distance = context.edge_to_edge_distance_bp
    if distance is None:
        distance = context.distance_bp

    passes_distance_threshold = (
        distance <= maximum_intergenic_distance_bp if distance is not None else None
    )

    supporting_conditions: list[str] = ["same_contig"]
    conflicting_conditions: list[str] = []

    if same_strand is True:
        supporting_conditions.append("same_strand")
    elif same_strand is False:
        conflicting_conditions.append("opposite_strand")
    else:
        conflicting_conditions.append("unknown_strand")

    if is_adjacent is True:
        supporting_conditions.append("no_intervening_gene")
    elif is_adjacent is False:
        conflicting_conditions.append("intervening_gene_present")
    else:
        conflicting_conditions.append("unknown_intervening_gene_count")

    if is_overlapping:
        conflicting_conditions.append("overlapping_features")
    elif overlap_bp == 0:
        supporting_conditions.append("non_overlapping")
    else:
        conflicting_conditions.append("unknown_overlap")

    if passes_distance_threshold is True:
        supporting_conditions.append("within_intergenic_threshold")
    elif passes_distance_threshold is False:
        conflicting_conditions.append("exceeds_intergenic_threshold")
    else:
        conflicting_conditions.append("unknown_intergenic_distance")

    if context.coordinate_position is CoordinatePosition.SAME_FEATURE:
        proxy_status = OperonProxyStatus.NOT_APPLICABLE
        conflicting_conditions.append("same_feature")
    elif (
        same_strand is True
        and is_adjacent is True
        and not is_overlapping
        and passes_distance_threshold is True
    ):
        proxy_status = OperonProxyStatus.SUPPORTED
    elif same_strand is False or is_overlapping:
        proxy_status = OperonProxyStatus.NOT_SUPPORTED
    elif (
        same_strand is None
        or is_adjacent is None
        or passes_distance_threshold is None
        or overlap_bp is None
    ):
        proxy_status = OperonProxyStatus.UNKNOWN
    elif supporting_conditions and conflicting_conditions:
        proxy_status = OperonProxyStatus.PARTIAL_SUPPORT
    else:
        proxy_status = OperonProxyStatus.NOT_SUPPORTED

    return OperonEvidence(
        status=(
            EvidenceStatus.NOT_APPLICABLE
            if proxy_status is OperonProxyStatus.NOT_APPLICABLE
            else EvidenceStatus.AVAILABLE
        ),
        origin=EvidenceOrigin.INFERRED,
        calculation_rule_version=OPERON_PROXY_RULE_VERSION,
        same_contig=True,
        same_strand=same_strand,
        is_adjacent=is_adjacent,
        intergenic_distance_bp=distance,
        overlap_bp=overlap_bp,
        intervening_gene_count=intervening_gene_count,
        transcriptional_order=_transcriptional_order(context),
        maximum_intergenic_distance_bp=maximum_intergenic_distance_bp,
        passes_distance_threshold=passes_distance_threshold,
        proxy_status=proxy_status,
        proxy_rule_id=OPERON_PROXY_RULE_ID,
        supporting_conditions=sorted(set(supporting_conditions)),
        conflicting_conditions=sorted(set(conflicting_conditions)),
        support=None,
        provenance=provenance,
        warnings=sorted(set(context.warnings)),
    )
