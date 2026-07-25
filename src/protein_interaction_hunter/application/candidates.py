"""MVP-1A query resolution, duplicate grouping, and candidate generation."""

import hashlib
from collections import defaultdict

from pydantic import Field

from protein_interaction_hunter.application.identifiers import IdentifierIndex, normalize_identifier
from protein_interaction_hunter.config import CandidateGenerationConfig
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceStatus,
    IdentifierMatchStatus,
)
from protein_interaction_hunter.models.genome import GeneCoordinate
from protein_interaction_hunter.models.protein import CandidateProtein, ProteinRecord, QueryProtein

_FRAGMENT_TERMS = ("fragment", "partial", "truncated", "incomplete")
_HYPOTHETICAL_TERMS = ("hypothetical", "unknown", "uncharacterized")


class CandidateGenerationResult(StrictModel):
    queries: list[QueryProtein]
    candidates: list[CandidateProtein]
    duplicate_groups: dict[str, list[str]] = Field(default_factory=dict)
    ambiguous_mapping_count: int = 0


def build_duplicate_groups(records: list[ProteinRecord]) -> dict[str, list[str]]:
    by_sequence: dict[str, list[str]] = defaultdict(list)
    for record in records:
        by_sequence[record.sequence].append(record.protein_id)
    groups: dict[str, list[str]] = {}
    for sequence, identifiers in sorted(by_sequence.items()):
        if len(identifiers) > 1:
            group_id = f"dup-{hashlib.sha256(sequence.encode('ascii')).hexdigest()[:16]}"
            groups[group_id] = sorted(identifiers)
    return dict(sorted(groups.items()))


def resolve_queries(query_ids: list[str], index: IdentifierIndex) -> list[QueryProtein]:
    resolved: list[QueryProtein] = []
    canonical_seen: set[str] = set()
    for query_id in query_ids:
        resolution = index.resolve(query_id)
        if resolution.status is IdentifierMatchStatus.AMBIGUOUS_MATCH:
            raise InputValidationError(
                f"Query identifier is ambiguous: {query_id} -> "
                + ", ".join(resolution.candidate_protein_ids)
            )
        if resolution.canonical_protein_id is None:
            raise InputValidationError(f"Query identifier not found: {query_id}")
        if resolution.canonical_protein_id in canonical_seen:
            raise InputValidationError(
                "Multiple query identifiers resolve to the same protein: "
                + resolution.canonical_protein_id
            )
        canonical_seen.add(resolution.canonical_protein_id)
        resolved.append(
            QueryProtein(
                query_id=query_id,
                protein_id=resolution.canonical_protein_id,
                resolution_method=resolution.status.value,
            )
        )
    return resolved


def _mapped_coordinates(
    records: list[GeneCoordinate], index: IdentifierIndex
) -> tuple[dict[str, GeneCoordinate], set[str]]:
    mapped: dict[str, list[tuple[int, GeneCoordinate]]] = defaultdict(list)
    ambiguous: set[str] = set()
    for record in records:
        resolution = index.resolve_gff(record)
        if resolution is None:
            continue
        if resolution.status is IdentifierMatchStatus.AMBIGUOUS_MATCH:
            ambiguous.update(resolution.candidate_protein_ids)
        elif resolution.canonical_protein_id:
            priority = 0 if record.protein_id == resolution.canonical_protein_id else 1
            priority += 0 if record.feature_type.casefold() == "cds" else 2
            mapped[resolution.canonical_protein_id].append((priority, record))
    selected: dict[str, GeneCoordinate] = {}
    for key, values in mapped.items():
        records_for_protein = [
            item[1] for item in sorted(values, key=lambda item: (item[0], item[1].start))
        ]
        primary = records_for_protein[0]
        old_locus_tag = primary.old_locus_tag or next(
            (record.old_locus_tag for record in records_for_protein if record.old_locus_tag),
            None,
        )
        locus_tag = primary.locus_tag or next(
            (record.locus_tag for record in records_for_protein if record.locus_tag), None
        )
        selected[key] = primary.model_copy(
            update={"old_locus_tag": old_locus_tag, "locus_tag": locus_tag}
        )
    return selected, ambiguous


def _mapped_annotations(
    records: list[AnnotationRecord], index: IdentifierIndex
) -> tuple[dict[str, AnnotationRecord], set[str]]:
    mapped: dict[str, AnnotationRecord] = {}
    ambiguous: set[str] = set()
    for record in records:
        resolution = index.resolve_annotation(record)
        if resolution.status is IdentifierMatchStatus.AMBIGUOUS_MATCH:
            ambiguous.update(resolution.candidate_protein_ids)
        elif resolution.canonical_protein_id:
            mapped.setdefault(resolution.canonical_protein_id, record)
    return mapped, ambiguous


def _fragment_reasons(
    protein: ProteinRecord, query: ProteinRecord, minimum_length: int
) -> list[str]:
    reasons: list[str] = []
    if len(protein.sequence) < minimum_length:
        reasons.append(f"length_below_minimum:{minimum_length}")
    description = protein.description.casefold()
    reasons.extend(f"description_keyword:{term}" for term in _FRAGMENT_TERMS if term in description)
    if len(protein.sequence) * 2 < len(query.sequence):
        reasons.append("less_than_half_query_length")
    return reasons


def _is_hypothetical(protein: ProteinRecord, annotation: AnnotationRecord | None) -> bool:
    searchable = " ".join(
        value
        for value in (
            protein.description,
            annotation.product if annotation else None,
            annotation.functional_category if annotation else None,
        )
        if value
    ).casefold()
    return any(term in searchable for term in _HYPOTHETICAL_TERMS)


def _apply_policies(
    *,
    protein: ProteinRecord,
    query: QueryProtein,
    duplicate_group: str | None,
    fragment_reasons: list[str],
    hypothetical: bool,
    has_coordinate: bool,
    policy: CandidateGenerationConfig,
) -> tuple[CandidateDisposition, list[str], list[str]]:
    disposition = CandidateDisposition.INCLUDED
    reasons: list[str] = []
    warnings: list[str] = []
    if protein.protein_id == query.protein_id:
        disposition = CandidateDisposition.EXCLUDED
        reasons.append("self_candidate")
    if duplicate_group:
        warnings.append("duplicate_sequence")
        if protein.protein_id != query.protein_id:
            if policy.duplicate_sequence_policy == "error":
                raise InputValidationError(
                    f"Duplicate sequence group rejected by policy: {duplicate_group}"
                )
            if policy.duplicate_sequence_policy == "exclude":
                disposition = CandidateDisposition.EXCLUDED
                reasons.append("duplicate_sequence_policy_exclude")
            elif disposition is not CandidateDisposition.EXCLUDED:
                disposition = CandidateDisposition.FLAGGED
                reasons.append("duplicate_sequence_policy_flag")
    if fragment_reasons:
        warnings.append("fragment_candidate")
        if policy.fragment_policy == "exclude":
            disposition = CandidateDisposition.EXCLUDED
            reasons.append("fragment_policy_exclude")
        elif policy.fragment_policy == "flag" and disposition is CandidateDisposition.INCLUDED:
            disposition = CandidateDisposition.FLAGGED
            reasons.append("fragment_policy_flag")
    if hypothetical:
        warnings.append("hypothetical_or_uncharacterized")
        if not policy.include_hypothetical_proteins:
            disposition = CandidateDisposition.EXCLUDED
            reasons.append("hypothetical_policy_exclude")
    if not has_coordinate:
        warnings.append("missing_coordinate")
        if not policy.include_missing_coordinates:
            disposition = CandidateDisposition.EXCLUDED
            reasons.append("missing_coordinate_policy_exclude")
    return disposition, reasons or ["eligible_candidate"], warnings


def generate_candidates(
    *,
    proteins: list[ProteinRecord],
    coordinates: list[GeneCoordinate],
    annotations: list[AnnotationRecord],
    query_ids: list[str],
    policy: CandidateGenerationConfig,
) -> CandidateGenerationResult:
    index = IdentifierIndex(proteins)
    index.add_gff_records(coordinates)
    index.add_annotations(annotations)
    queries = resolve_queries(query_ids, index)
    by_id = {record.protein_id: record for record in proteins}
    coordinate_map, coordinate_ambiguous = _mapped_coordinates(coordinates, index)
    annotation_map, annotation_ambiguous = _mapped_annotations(annotations, index)
    ambiguous = coordinate_ambiguous | annotation_ambiguous
    duplicate_groups = build_duplicate_groups(proteins)
    duplicate_by_id = {
        protein_id: group_id
        for group_id, identifiers in duplicate_groups.items()
        for protein_id in identifiers
    }
    candidates: list[CandidateProtein] = []
    for query in queries:
        query_record = by_id[query.protein_id]
        query_coordinate = coordinate_map.get(query.protein_id)
        for protein in sorted(proteins, key=lambda item: item.protein_id):
            coordinate = coordinate_map.get(protein.protein_id)
            annotation = annotation_map.get(protein.protein_id)
            annotation_available = (
                annotation is not None and annotation.status is EvidenceStatus.AVAILABLE
            )
            duplicate_group = duplicate_by_id.get(protein.protein_id)
            fragment_reasons = _fragment_reasons(protein, query_record, policy.minimum_length_aa)
            hypothetical = _is_hypothetical(protein, annotation)
            disposition, reasons, warnings = _apply_policies(
                protein=protein,
                query=query,
                duplicate_group=duplicate_group,
                fragment_reasons=fragment_reasons,
                hypothetical=hypothetical,
                has_coordinate=coordinate is not None,
                policy=policy,
            )
            if not annotation_available:
                warnings.append("missing_annotation")
            if protein.protein_id in ambiguous:
                warnings.append("ambiguous_identifier_mapping")
            match_status = IdentifierMatchStatus.EXACT_MATCH
            if protein.protein_id in ambiguous:
                match_status = IdentifierMatchStatus.AMBIGUOUS_MATCH
            elif coordinate is None and annotation is None:
                match_status = IdentifierMatchStatus.NO_MATCH
            elif (coordinate and coordinate.protein_id != protein.protein_id) or (
                annotation and annotation.protein_id != protein.protein_id
            ):
                match_status = IdentifierMatchStatus.UNIQUE_ALIAS_MATCH
            original = {
                "protein_id": [protein.protein_id],
                "gene_id": [protein.gene_id] if protein.gene_id else [],
                "locus_tag": [protein.locus_tag] if protein.locus_tag else [],
                "old_locus_tag": [coordinate.old_locus_tag]
                if coordinate and coordinate.old_locus_tag
                else [],
                "gff_id": [coordinate.feature_id] if coordinate and coordinate.feature_id else [],
                "gff_parent": [coordinate.parent_id] if coordinate and coordinate.parent_id else [],
                "gff_protein_id": [coordinate.protein_id]
                if coordinate and coordinate.protein_id
                else [],
                "annotation_protein_id": [annotation.protein_id] if annotation else [],
                "annotation_gene_id": [annotation.gene_name]
                if annotation and annotation.gene_name
                else [],
                "annotation_locus_tag": [annotation.locus_tag]
                if annotation and annotation.locus_tag
                else [],
            }
            candidates.append(
                CandidateProtein(
                    query_id=query.query_id,
                    protein_id=protein.protein_id,
                    disposition=disposition,
                    disposition_reasons=reasons,
                    duplicate_sequence_group=duplicate_group,
                    sequence_length=len(protein.sequence),
                    description=protein.description,
                    gene_id=protein.gene_id or (coordinate.parent_id if coordinate else None),
                    locus_tag=protein.locus_tag or (coordinate.locus_tag if coordinate else None),
                    old_locus_tag=coordinate.old_locus_tag if coordinate else None,
                    contig=coordinate.seqid if coordinate else None,
                    strand=coordinate.strand if coordinate else None,
                    has_coordinate=coordinate is not None,
                    has_annotation=annotation_available,
                    coordinate_status=(
                        EvidenceStatus.AVAILABLE if coordinate else EvidenceStatus.MISSING
                    ),
                    annotation_status=annotation.status if annotation else EvidenceStatus.MISSING,
                    same_contig_as_query=(
                        coordinate.seqid == query_coordinate.seqid
                        if coordinate and query_coordinate
                        else None
                    ),
                    is_duplicate_sequence=duplicate_group is not None,
                    is_fragment_candidate=bool(fragment_reasons),
                    fragment_reasons=fragment_reasons,
                    is_hypothetical=hypothetical,
                    identifier_match_status=match_status,
                    original_identifiers=original,
                    normalized_identifiers={
                        key: [normalize_identifier(value) for value in values]
                        for key, values in original.items()
                    },
                    warnings=warnings,
                )
            )
    return CandidateGenerationResult(
        queries=queries,
        candidates=candidates,
        duplicate_groups=duplicate_groups,
        ambiguous_mapping_count=len(ambiguous),
    )
