"""Deterministic local gene-fusion pair evaluation."""

from __future__ import annotations

from collections import defaultdict

from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    FusionEvidence,
    FusionObservation,
)

FUSION_ENGINE_VERSION = "mvp1i-gene-fusion-v1"
FusionIndex = dict[tuple[str, str], list[FusionObservation]]


def canonical_fusion_pair(protein_a: str, protein_b: str) -> tuple[str, str]:
    return (protein_a, protein_b) if protein_a <= protein_b else (protein_b, protein_a)


def _record_regions(record: FusionObservation) -> tuple[tuple[int, int], tuple[int, int]]:
    pair = canonical_fusion_pair(record.query_protein_id, record.candidate_protein_id)
    if record.query_protein_id == pair[0]:
        return (
            (record.query_component_start, record.query_component_end),
            (record.candidate_component_start, record.candidate_component_end),
        )
    return (
        (record.candidate_component_start, record.candidate_component_end),
        (record.query_component_start, record.query_component_end),
    )


def build_fusion_index(records: list[FusionObservation]) -> FusionIndex:
    index: dict[tuple[str, str], list[FusionObservation]] = defaultdict(list)
    for record in records:
        index[canonical_fusion_pair(record.query_protein_id, record.candidate_protein_id)].append(
            record
        )
    return {
        pair: sorted(
            pair_records,
            key=lambda record: (
                record.reference_organism,
                record.fusion_protein_id,
                _record_regions(record),
                record.source_record_id or "",
                record.source or "",
            ),
        )
        for pair, pair_records in sorted(index.items())
    }


def _oriented_coverages(
    record: FusionObservation,
    query_protein_id: str,
) -> tuple[float | None, float | None]:
    if record.query_protein_id == query_protein_id:
        return record.query_component_coverage, record.candidate_component_coverage
    return record.candidate_component_coverage, record.query_component_coverage


def evaluate_fusion_pair(
    query_protein_id: str,
    candidate_protein_id: str,
    fusion_index: FusionIndex,
    *,
    minimum_supporting_records: int,
    minimum_component_coverage: float,
    maximum_component_overlap_fraction: float,
) -> FusionEvidence:
    records = fusion_index.get(canonical_fusion_pair(query_protein_id, candidate_protein_id), [])
    provenance = [
        EvidenceProvenance(
            source_name="local_fusion_table",
            source_version=FUSION_ENGINE_VERSION,
            method="canonical_pair_component_coverage_and_overlap_evaluation",
            metadata={
                "coordinate_system": "1-based inclusive",
                "minimum_supporting_records": minimum_supporting_records,
                "minimum_component_coverage": minimum_component_coverage,
                "maximum_component_overlap_fraction": (maximum_component_overlap_fraction),
                "ranking_propagation": False,
                "score_propagation": False,
                "relationship_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]
    if not records:
        return FusionEvidence(
            status=EvidenceStatus.AVAILABLE,
            origin=EvidenceOrigin.INFERRED,
            query_protein_id=query_protein_id,
            candidate_protein_id=candidate_protein_id,
            supporting_record_count=0,
            qualifying_record_count=0,
            pair_supported=False,
            conflicting_terms=["no_fusion_record"],
            calculation_rule_version=FUSION_ENGINE_VERSION,
            source="local_table",
            provenance=provenance,
        )

    qualifying: list[FusionObservation] = []
    undetermined: list[FusionObservation] = []
    conflicting_terms: set[str] = set()
    oriented_coverages = [_oriented_coverages(record, query_protein_id) for record in records]
    for record, (query_coverage, candidate_coverage) in zip(
        records, oriented_coverages, strict=True
    ):
        overlap = record.component_overlap_fraction
        assert overlap is not None
        known_failure = False
        if query_coverage is not None and query_coverage < minimum_component_coverage:
            conflicting_terms.add("low_query_component_coverage")
            known_failure = True
        if candidate_coverage is not None and candidate_coverage < minimum_component_coverage:
            conflicting_terms.add("low_candidate_component_coverage")
            known_failure = True
        if overlap > maximum_component_overlap_fraction:
            conflicting_terms.update(
                {"excessive_component_overlap", "invalid_component_separation"}
            )
            known_failure = True
        missing_coverage = query_coverage is None or candidate_coverage is None
        if missing_coverage:
            conflicting_terms.add("missing_component_coverage")
        if not known_failure and missing_coverage:
            undetermined.append(record)
        elif not known_failure and query_coverage is not None and candidate_coverage is not None:
            qualifying.append(record)

    qualifying_fusion_ids = sorted({record.fusion_protein_id for record in qualifying})
    qualifying_organisms = sorted({record.reference_organism for record in qualifying})
    potential_count = len(qualifying) + len(undetermined)
    if len(qualifying) >= minimum_supporting_records:
        status = EvidenceStatus.AVAILABLE
        pair_supported: bool | None = bool(qualifying_fusion_ids and qualifying_organisms)
    elif potential_count >= minimum_supporting_records:
        status = EvidenceStatus.MISSING
        pair_supported = None
    else:
        status = EvidenceStatus.AVAILABLE
        pair_supported = False

    support_terms: set[str] = set()
    if qualifying:
        support_terms.update(
            {
                "qualifying_fusion_record",
                "separate_fusion_components",
                "sufficient_query_component_coverage",
                "sufficient_candidate_component_coverage",
            }
        )
    if len(qualifying_fusion_ids) > 1:
        support_terms.add("multiple_supporting_fusions")
    if len(qualifying_organisms) > 1:
        support_terms.add("multiple_reference_organisms")
    if len(qualifying) < minimum_supporting_records:
        conflicting_terms.add("insufficient_supporting_records")

    query_coverages = [value[0] for value in oriented_coverages if value[0] is not None]
    candidate_coverages = [value[1] for value in oriented_coverages if value[1] is not None]
    overlaps = [
        record.component_overlap_fraction
        for record in records
        if record.component_overlap_fraction is not None
    ]
    sources = sorted({record.source for record in records if record.source})
    source_record_ids = sorted(
        {record.source_record_id for record in records if record.source_record_id}
    )
    warnings = (
        [f"fusion_records_with_missing_component_coverage:{len(undetermined)}"]
        if undetermined
        else []
    )
    return FusionEvidence(
        status=status,
        origin=EvidenceOrigin.INFERRED,
        query_protein_id=query_protein_id,
        candidate_protein_id=candidate_protein_id,
        supporting_record_count=len(records),
        qualifying_record_count=len(qualifying),
        reference_organisms=sorted({record.reference_organism for record in records}),
        fusion_protein_ids=sorted({record.fusion_protein_id for record in records}),
        best_query_component_coverage=max(query_coverages) if query_coverages else None,
        best_candidate_component_coverage=(
            max(candidate_coverages) if candidate_coverages else None
        ),
        minimum_component_overlap_fraction=min(overlaps) if overlaps else None,
        pair_supported=pair_supported,
        support_terms=sorted(support_terms),
        conflicting_terms=sorted(conflicting_terms),
        calculation_rule_version=FUSION_ENGINE_VERSION,
        source="|".join(sources) or "local_table",
        source_record_ids=source_record_ids,
        warnings=warnings,
        provenance=provenance,
    )
