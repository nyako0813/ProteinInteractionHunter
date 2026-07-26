"""Deterministic local orthology pair evaluation."""

from __future__ import annotations

from collections import defaultdict

from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    OrthologRecord,
)

ORTHOLOGY_ENGINE_VERSION = "mvp1g-orthology-v1"


def build_orthology_index(
    records: list[OrthologRecord],
) -> dict[str, list[OrthologRecord]]:
    index: dict[str, list[OrthologRecord]] = defaultdict(list)

    for record in records:
        index[record.protein_id].append(record)

    return {
        protein_id: sorted(
            protein_records,
            key=lambda item: (
                item.reference_id,
                item.orthogroup or "",
                item.ortholog_id or "",
                item.evalue if item.evalue is not None else float("inf"),
            ),
        )
        for protein_id, protein_records in index.items()
    }


def evaluate_orthology_pair(
    query_protein_id: str,
    candidate_protein_id: str,
    orthology_index: dict[str, list[OrthologRecord]],
) -> list[OrthologRecord]:
    """Evaluate shared orthology support without scoring or ranking."""

    query_records = orthology_index.get(query_protein_id, [])
    candidate_records = orthology_index.get(candidate_protein_id, [])

    provenance = [
        EvidenceProvenance(
            source_name="local_orthology_table",
            source_version=ORTHOLOGY_ENGINE_VERSION,
            method="deterministic_reference_and_orthogroup_matching",
            metadata={
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]

    if not query_records or not candidate_records:
        missing: list[str] = []

        if not query_records:
            missing.append("query_orthology_annotation")
        if not candidate_records:
            missing.append("candidate_orthology_annotation")

        return [
            OrthologRecord(
                status=EvidenceStatus.MISSING,
                origin=EvidenceOrigin.ORTHOLOG_TRANSFERRED,
                calculation_rule_version=ORTHOLOGY_ENGINE_VERSION,
                protein_id=candidate_protein_id,
                reference_id="missing",
                paired_protein_id=query_protein_id,
                pair_supported=None,
                conflicting_terms=missing,
                provenance=provenance,
                warnings=["missing_orthology_annotation:" + ",".join(missing)],
            )
        ]

    evidence_records: list[OrthologRecord] = []

    for query_record in query_records:
        for candidate_record in candidate_records:
            if query_record.reference_id != candidate_record.reference_id:
                continue

            shared_orthogroup = (
                query_record.orthogroup is not None
                and candidate_record.orthogroup is not None
                and query_record.orthogroup == candidate_record.orthogroup
            )

            same_ortholog = (
                query_record.ortholog_id is not None
                and candidate_record.ortholog_id is not None
                and query_record.ortholog_id == candidate_record.ortholog_id
            )

            pair_supported = shared_orthogroup or same_ortholog
            paralog_ambiguity = query_record.paralog_ambiguity or candidate_record.paralog_ambiguity

            support_terms: list[str] = []
            conflicting_terms: list[str] = []

            if shared_orthogroup:
                support_terms.append("shared_orthogroup")
            if same_ortholog:
                support_terms.append("shared_ortholog_id")
            if not pair_supported:
                conflicting_terms.append("no_shared_orthology_support")
            if paralog_ambiguity:
                conflicting_terms.append("paralog_ambiguity")

            evidence_records.append(
                candidate_record.model_copy(
                    update={
                        "calculation_rule_version": (ORTHOLOGY_ENGINE_VERSION),
                        "paired_protein_id": query_protein_id,
                        "paired_reference_id": query_record.reference_id,
                        "paired_ortholog_id": query_record.ortholog_id,
                        "paired_orthogroup": query_record.orthogroup,
                        "shared_orthogroup": shared_orthogroup,
                        "pair_supported": pair_supported,
                        "paralog_ambiguity": paralog_ambiguity,
                        "support_terms": sorted(set(support_terms)),
                        "conflicting_terms": sorted(set(conflicting_terms)),
                        "provenance": provenance,
                    }
                )
            )

    if evidence_records:
        return evidence_records

    candidate_record = candidate_records[0]

    return [
        candidate_record.model_copy(
            update={
                "calculation_rule_version": ORTHOLOGY_ENGINE_VERSION,
                "paired_protein_id": query_protein_id,
                "shared_orthogroup": False,
                "pair_supported": False,
                "conflicting_terms": ["no_shared_reference"],
                "provenance": provenance,
            }
        )
    ]
