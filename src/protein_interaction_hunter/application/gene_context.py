"""Deterministic, observation-only gene-context evidence for MVP-1B."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any

from protein_interaction_hunter.application.identifiers import IdentifierIndex
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import (
    ContextCompleteness,
    CoordinatePosition,
    EvidenceOrigin,
    EvidenceStatus,
    IdentifierMatchStatus,
    RelativePosition,
    StrandRelationship,
)
from protein_interaction_hunter.models.evidence import EvidenceProvenance, GenomeContextEvidence
from protein_interaction_hunter.models.genome import (
    GeneCoordinate,
    GffDocument,
    NormalizedFeature,
    SequenceRegion,
)
from protein_interaction_hunter.models.protein import ProteinRecord

GENE_CONTEXT_RULE_VERSION = "mvp1b-gene-context-v1"


@dataclass(frozen=True)
class CoordinateIndex:
    by_protein: dict[str, NormalizedFeature]
    ambiguous_proteins: frozenset[str]
    features_by_contig: dict[str, tuple[NormalizedFeature, ...]]
    sequence_regions: dict[str, SequenceRegion]
    warnings: tuple[str, ...]


def _coordinate_key(record: GeneCoordinate) -> tuple[str, int, int, str]:
    return (record.seqid, record.start, record.end, record.strand or "?")


def _choose_representative(
    protein_id: str, records: list[GeneCoordinate], status: IdentifierMatchStatus
) -> NormalizedFeature | None:
    contigs = {record.seqid for record in records}
    if len(contigs) > 1:
        raise InputValidationError(
            f"Contradictory replicons for protein {protein_id}: {', '.join(sorted(contigs))}"
        )
    cds = {
        key
        for record in records
        if record.feature_type.casefold() == "cds"
        for key in [_coordinate_key(record)]
    }
    all_coordinates = {_coordinate_key(record) for record in records}
    selected_key: tuple[str, int, int, str] | None
    warning_list = [warning for record in records for warning in record.warnings]
    if len(cds) == 1:
        selected_key = next(iter(cds))
        if len(all_coordinates) > 1:
            warning_list.append("gene_cds_extent_difference:cds_selected")
    elif len(cds) > 1 or len(all_coordinates) > 1:
        return None
    else:
        selected_key = next(iter(all_coordinates))
    selected_records = [record for record in records if _coordinate_key(record) == selected_key]
    selected = min(
        selected_records,
        key=lambda record: (
            0 if record.feature_type.casefold() == "cds" else 1,
            record.feature_id or "",
        ),
    )
    parents = sorted({parent for record in records for parent in record.parent_ids})
    feature_ids = sorted({record.feature_id for record in records if record.feature_id})
    locus_tag = selected.locus_tag or next(
        (record.locus_tag for record in records if record.locus_tag), None
    )
    gene_id = selected.parent_id or next(
        (
            record.feature_id
            for record in records
            if record.feature_type.casefold() == "gene" and record.feature_id
        ),
        None,
    )
    return NormalizedFeature(
        representative_id=f"protein:{protein_id}",
        protein_id=protein_id,
        gene_id=gene_id,
        locus_tag=locus_tag,
        old_locus_tag=selected.old_locus_tag
        or next((record.old_locus_tag for record in records if record.old_locus_tag), None),
        seqid=selected.seqid,
        start=selected.start,
        end=selected.end,
        strand=selected.strand or "?",
        feature_type=selected.feature_type,
        coordinate_source=selected.source,
        identifier_status=status,
        parent_identifiers=parents,
        source_feature_ids=feature_ids,
        is_gene=True,
        warnings=sorted(set(warning_list)),
    )


def _unmapped_representatives(
    records: list[GeneCoordinate], mapped_record_ids: set[int]
) -> list[NormalizedFeature]:
    """Collapse unmapped gene/child records without double-counting gene/CDS/RNA units."""
    by_unit: dict[str, list[GeneCoordinate]] = defaultdict(list)
    known_gene_ids = {
        record.feature_id
        for record in records
        if record.feature_type.casefold() == "gene" and record.feature_id
    }
    for record in records:
        if id(record) in mapped_record_ids:
            continue
        parent = next((value for value in record.parent_ids if value in known_gene_ids), None)
        unit = (
            parent
            or record.feature_id
            or (f"anonymous:{record.seqid}:{record.start}:{record.end}:{record.feature_type}")
        )
        by_unit[unit].append(record)
    representatives: list[NormalizedFeature] = []
    for unit, values in sorted(by_unit.items()):
        genes = [record for record in values if record.feature_type.casefold() == "gene"]
        selected = min(
            genes or values,
            key=lambda record: (
                record.start,
                record.end,
                record.feature_type,
                record.feature_id or "",
            ),
        )
        representatives.append(
            NormalizedFeature(
                representative_id=f"gff:{unit}",
                gene_id=unit if unit in known_gene_ids else selected.parent_id,
                locus_tag=selected.locus_tag,
                seqid=selected.seqid,
                old_locus_tag=selected.old_locus_tag,
                start=selected.start,
                end=selected.end,
                strand=selected.strand or "?",
                feature_type=selected.feature_type,
                coordinate_source=selected.source,
                identifier_status=IdentifierMatchStatus.NO_MATCH,
                parent_identifiers=sorted(
                    {parent for record in values for parent in record.parent_ids}
                ),
                source_feature_ids=sorted(
                    {record.feature_id for record in values if record.feature_id}
                ),
                is_gene=bool(genes)
                or any(record.feature_type.casefold() == "cds" for record in values),
                warnings=sorted(
                    {"unmapped_gff_feature"}
                    | {warning for record in values for warning in record.warnings}
                ),
            )
        )
    return representatives


def build_coordinate_index(proteins: list[ProteinRecord], document: GffDocument) -> CoordinateIndex:
    identifier_index = IdentifierIndex(proteins)
    identifier_index.add_gff_records(document.features)
    grouped: dict[str, list[GeneCoordinate]] = defaultdict(list)
    statuses: dict[str, IdentifierMatchStatus] = {}
    ambiguous: set[str] = set()
    mapped_record_ids: set[int] = set()
    for record in document.features:
        resolution = identifier_index.resolve_gff(record)
        if resolution is None:
            continue
        if resolution.status is IdentifierMatchStatus.AMBIGUOUS_MATCH:
            ambiguous.update(resolution.candidate_protein_ids)
            continue
        if resolution.canonical_protein_id:
            protein_id = resolution.canonical_protein_id
            grouped[protein_id].append(record)
            mapped_record_ids.add(id(record))
            if record.protein_id == protein_id:
                statuses[protein_id] = IdentifierMatchStatus.EXACT_MATCH
            else:
                statuses.setdefault(protein_id, IdentifierMatchStatus.UNIQUE_ALIAS_MATCH)
    by_protein: dict[str, NormalizedFeature] = {}
    for protein_id, records in sorted(grouped.items()):
        representative = _choose_representative(protein_id, records, statuses[protein_id])
        if representative is None:
            ambiguous.add(protein_id)
        else:
            by_protein[protein_id] = representative
    all_representatives = list(by_protein.values()) + _unmapped_representatives(
        document.features, mapped_record_ids
    )
    by_contig: dict[str, list[NormalizedFeature]] = defaultdict(list)
    for feature in all_representatives:
        by_contig[feature.seqid].append(feature)
    ordered = {
        seqid: tuple(
            sorted(values, key=lambda item: (item.start, item.end, item.representative_id))
        )
        for seqid, values in sorted(by_contig.items())
    }
    return CoordinateIndex(
        by_protein=by_protein,
        ambiguous_proteins=frozenset(ambiguous),
        features_by_contig=ordered,
        sequence_regions=document.sequence_regions,
        warnings=tuple(document.warnings),
    )


def _base_evidence(
    status: EvidenceStatus,
    warning: str,
    query: NormalizedFeature | None,
    candidate: NormalizedFeature | None,
) -> GenomeContextEvidence:
    return GenomeContextEvidence(
        status=status,
        origin=EvidenceOrigin.ANNOTATION,
        calculation_rule_version=GENE_CONTEXT_RULE_VERSION,
        same_contig=(query.seqid == candidate.seqid) if query and candidate else None,
        same_seqid=(query.seqid == candidate.seqid) if query and candidate else None,
        query_contig=query.seqid if query else None,
        candidate_contig=candidate.seqid if candidate else None,
        query_start=query.start if query else None,
        query_end=query.end if query else None,
        query_strand=query.strand if query else None,
        candidate_start=candidate.start if candidate else None,
        candidate_end=candidate.end if candidate else None,
        candidate_strand=candidate.strand if candidate else None,
        context_completeness=ContextCompleteness.UNKNOWN,
        provenance=[
            EvidenceProvenance(
                source_name="gff3_gene_context",
                source_version=GENE_CONTEXT_RULE_VERSION,
                method="deterministic_interval_and_representative_feature_analysis",
            )
        ],
        warnings=[warning],
    )


def _coordinate_position(
    query: NormalizedFeature, candidate: NormalizedFeature
) -> CoordinatePosition:
    if query.representative_id == candidate.representative_id:
        return CoordinatePosition.SAME_FEATURE
    if candidate.end < query.start:
        return CoordinatePosition.LEFT_OF_QUERY
    if candidate.start > query.end:
        return CoordinatePosition.RIGHT_OF_QUERY
    return CoordinatePosition.OVERLAPPING


def _relative_position(
    query: NormalizedFeature, coordinate_position: CoordinatePosition
) -> RelativePosition:
    if coordinate_position is CoordinatePosition.SAME_FEATURE:
        return RelativePosition.SAME_FEATURE
    if coordinate_position is CoordinatePosition.OVERLAPPING:
        return RelativePosition.OVERLAPPING
    if query.strand == "+":
        return (
            RelativePosition.UPSTREAM
            if coordinate_position is CoordinatePosition.LEFT_OF_QUERY
            else RelativePosition.DOWNSTREAM
        )
    if query.strand == "-":
        return (
            RelativePosition.DOWNSTREAM
            if coordinate_position is CoordinatePosition.LEFT_OF_QUERY
            else RelativePosition.UPSTREAM
        )
    return RelativePosition.UNKNOWN


def _strand_relationship(
    query: NormalizedFeature, candidate: NormalizedFeature
) -> StrandRelationship:
    if query.strand not in {"+", "-"} or candidate.strand not in {"+", "-"}:
        return StrandRelationship.UNKNOWN
    if query.strand == candidate.strand:
        return StrandRelationship.SAME_DIRECTION
    if max(query.start, candidate.start) <= min(query.end, candidate.end):
        return StrandRelationship.OPPOSITE_PARALLEL
    left, right = sorted((query, candidate), key=lambda item: (item.start, item.end))
    if left.strand == "+" and right.strand == "-":
        return StrandRelationship.CONVERGENT
    return StrandRelationship.DIVERGENT


def _completeness(
    query_index: int,
    candidate_index: int,
    feature_count: int,
    neighborhood_gene_count: int,
    region: SequenceRegion | None,
) -> ContextCompleteness:
    if region is None:
        return ContextCompleteness.UNKNOWN
    left_truncated = min(query_index, candidate_index) < neighborhood_gene_count
    right_truncated = (
        feature_count - 1 - max(query_index, candidate_index) < neighborhood_gene_count
    )
    if left_truncated and right_truncated:
        return ContextCompleteness.BOTH_TRUNCATED
    if left_truncated:
        return ContextCompleteness.LEFT_TRUNCATED
    if right_truncated:
        return ContextCompleteness.RIGHT_TRUNCATED
    return ContextCompleteness.COMPLETE


def calculate_gene_context(
    query_id: str,
    candidate_id: str,
    index: CoordinateIndex,
    neighborhood_gene_count: int,
) -> GenomeContextEvidence:
    query = index.by_protein.get(query_id)
    candidate = index.by_protein.get(candidate_id)
    if query_id in index.ambiguous_proteins or candidate_id in index.ambiguous_proteins:
        ambiguous_id = query_id if query_id in index.ambiguous_proteins else candidate_id
        return _base_evidence(
            EvidenceStatus.FAILED,
            f"ambiguous_coordinate_mapping:{ambiguous_id}",
            query,
            candidate,
        )
    if query is None or candidate is None:
        missing_id = query_id if query is None else candidate_id
        return _base_evidence(
            EvidenceStatus.MISSING,
            f"missing_coordinate:{missing_id}",
            query,
            candidate,
        )
    common: dict[str, Any] = {
        "origin": EvidenceOrigin.ANNOTATION,
        "calculation_rule_version": GENE_CONTEXT_RULE_VERSION,
        "query_contig": query.seqid,
        "candidate_contig": candidate.seqid,
        "query_start": query.start,
        "query_end": query.end,
        "query_strand": query.strand,
        "candidate_start": candidate.start,
        "candidate_end": candidate.end,
        "candidate_strand": candidate.strand,
        "provenance": [
            EvidenceProvenance(
                source_name="gff3_gene_context",
                source_version=GENE_CONTEXT_RULE_VERSION,
                method="deterministic_interval_and_representative_feature_analysis",
            )
        ],
    }
    if query.seqid != candidate.seqid:
        return GenomeContextEvidence(
            status=EvidenceStatus.NOT_APPLICABLE,
            same_contig=False,
            same_seqid=False,
            strand_relationship=StrandRelationship.DIFFERENT_CONTIG,
            strand_relation=StrandRelationship.DIFFERENT_CONTIG.value,
            relative_position=RelativePosition.DIFFERENT_CONTIG,
            coordinate_position=CoordinatePosition.DIFFERENT_CONTIG,
            context_completeness=ContextCompleteness.UNKNOWN,
            warnings=["different_contig:distance_not_applicable"],
            **common,
        )
    features = index.features_by_contig[query.seqid]
    positions = {feature.representative_id: position for position, feature in enumerate(features)}
    query_index = positions[query.representative_id]
    candidate_index = positions[candidate.representative_id]
    coordinate_position = _coordinate_position(query, candidate)
    overlap_bp = max(0, min(query.end, candidate.end) - max(query.start, candidate.start) + 1)
    distance_bp = (
        0 if overlap_bp else max(query.start, candidate.start) - min(query.end, candidate.end) - 1
    )
    if coordinate_position in {
        CoordinatePosition.SAME_FEATURE,
        CoordinatePosition.OVERLAPPING,
    }:
        between: list[NormalizedFeature] = []
    else:
        left, right = sorted((query, candidate), key=lambda item: (item.start, item.end))
        between = [
            feature
            for feature in features
            if feature.representative_id
            not in {query.representative_id, candidate.representative_id}
            and feature.end > left.end
            and feature.start < right.start
        ]
    region = index.sequence_regions.get(query.seqid)
    warnings = list(query.warnings) + list(candidate.warnings)
    boundary_flags: list[str] = []
    if region is None:
        warnings.append(f"missing_sequence_region:{query.seqid}")
        boundary_flags.append("unknown_contig_boundaries")
    strand_relationship = _strand_relationship(query, candidate)
    feature_delta = abs(candidate_index - query_index)
    return GenomeContextEvidence(
        status=EvidenceStatus.AVAILABLE,
        same_contig=True,
        same_seqid=True,
        distance_bp=distance_bp,
        edge_to_edge_distance_bp=distance_bp,
        overlap_bp=overlap_bp,
        relative_position=_relative_position(query, coordinate_position),
        coordinate_position=coordinate_position,
        strand_relationship=strand_relationship,
        strand_relation=strand_relationship.value,
        intervening_feature_count=len(between),
        intervening_gene_count=sum(feature.is_gene for feature in between),
        query_feature_index=query_index,
        candidate_feature_index=candidate_index,
        feature_index_delta=feature_delta,
        within_neighborhood_window=feature_delta <= neighborhood_gene_count,
        within_neighborhood_gene_count=feature_delta <= neighborhood_gene_count,
        query_left_edge_distance_bp=query.start - region.start if region else None,
        query_right_edge_distance_bp=region.end - query.end if region else None,
        candidate_left_edge_distance_bp=candidate.start - region.start if region else None,
        candidate_right_edge_distance_bp=region.end - candidate.end if region else None,
        query_distance_to_contig_left_edge=query.start - region.start if region else None,
        query_distance_to_contig_right_edge=region.end - query.end if region else None,
        candidate_distance_to_contig_left_edge=(candidate.start - region.start if region else None),
        candidate_distance_to_contig_right_edge=(region.end - candidate.end if region else None),
        context_completeness=_completeness(
            query_index,
            candidate_index,
            len(features),
            neighborhood_gene_count,
            region,
        ),
        boundary_flags=boundary_flags,
        warnings=sorted(set(warnings)),
        **common,
    )
