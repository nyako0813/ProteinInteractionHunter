"""Deterministic local domain-pair evaluation."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from protein_interaction_hunter.models.domain import DomainAnnotationRecord
from protein_interaction_hunter.models.domain_rules import DomainPairRuleset
from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import (
    DomainEvidence,
    EvidenceProvenance,
)

DOMAIN_PAIR_ENGINE_VERSION = "mvp1e-domain-pair-v1"


def build_domain_index(
    records: list[DomainAnnotationRecord],
) -> dict[str, list[DomainAnnotationRecord]]:
    index: dict[str, list[DomainAnnotationRecord]] = defaultdict(list)

    for record in records:
        index[record.protein_id].append(record)

    return {
        protein_id: sorted(
            protein_records,
            key=lambda item: (
                item.architecture_index,
                item.start,
                item.end,
                item.accession,
            ),
        )
        for protein_id, protein_records in index.items()
    }


def _role_accessions(
    ruleset: DomainPairRuleset,
) -> dict[str, set[str]]:
    return {
        role.role_id: set(role.accessions)
        for role in ruleset.roles
    }


def evaluate_domain_pairs(
    query_protein_id: str,
    candidate_protein_id: str,
    domain_index: dict[str, list[DomainAnnotationRecord]],
    ruleset: DomainPairRuleset,
    ruleset_path: Path | None = None,
) -> list[DomainEvidence]:
    """Evaluate domain-pair rules without scoring or ranking."""

    query_domains = domain_index.get(query_protein_id, [])
    candidate_domains = domain_index.get(candidate_protein_id, [])
    serialized_path = str(ruleset_path.resolve()) if ruleset_path else None

    provenance = [
        EvidenceProvenance(
            source_name="local_domain_pair_rules",
            source_version=ruleset.ruleset_version,
            method="deterministic_accession_pair_matching",
            metadata={
                "engine_version": DOMAIN_PAIR_ENGINE_VERSION,
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
                "shared_domain_physical_weight": 0.0,
            },
        )
    ]

    if not query_domains or not candidate_domains:
        missing: list[str] = []
        if not query_domains:
            missing.append("query_domain_annotation")
        if not candidate_domains:
            missing.append("candidate_domain_annotation")

        return [
            DomainEvidence(
                status=EvidenceStatus.MISSING,
                origin=EvidenceOrigin.ANNOTATION,
                calculation_rule_version=DOMAIN_PAIR_ENGINE_VERSION,
                protein_id=candidate_protein_id,
                paired_protein_id=query_protein_id,
                pair_matched=None,
                conflicting_terms=missing,
                ruleset_path=serialized_path,
                provenance=provenance,
                warnings=[
                    "missing_domain_annotation:" + ",".join(missing)
                ],
            )
        ]

    roles = _role_accessions(ruleset)
    evidence_records: list[DomainEvidence] = []

    for pair_rule in ruleset.pair_rules:
        orientations = [
            (pair_rule.query_role, pair_rule.candidate_role),
        ]
        if pair_rule.query_role != pair_rule.candidate_role:
            orientations.append(
                (pair_rule.candidate_role, pair_rule.query_role)
            )

        for query_role, candidate_role in orientations:
            query_accessions = roles.get(query_role, set())
            candidate_accessions = roles.get(
                candidate_role,
                set(),
            )

            for query_domain in query_domains:
                if query_domain.accession not in query_accessions:
                    continue

                for candidate_domain in candidate_domains:
                    if candidate_domain.accession not in candidate_accessions:
                        continue

                    is_shared = (
                        query_domain.accession
                        == candidate_domain.accession
                    )

                    if is_shared and not pair_rule.allow_shared_accession:
                        continue

                    evidence_records.append(
                        DomainEvidence(
                            status=EvidenceStatus.AVAILABLE,
                            origin=EvidenceOrigin.ANNOTATION,
                            calculation_rule_version=(
                                DOMAIN_PAIR_ENGINE_VERSION
                            ),
                            protein_id=candidate_protein_id,
                            source=candidate_domain.source,
                            accession=candidate_domain.accession,
                            name=candidate_domain.name,
                            start=candidate_domain.start,
                            end=candidate_domain.end,
                            architecture_index=(
                                candidate_domain.architecture_index
                            ),
                            role=candidate_role,
                            pair_rule_id=pair_rule.rule_id,
                            paired_protein_id=query_protein_id,
                            paired_accession=query_domain.accession,
                            is_shared=is_shared,
                            pair_matched=True,
                            support_terms=sorted(
                                set(
                                    pair_rule.support_terms
                                    + [
                                        f"query:{query_domain.accession}",
                                        (
                                            "candidate:"
                                            f"{candidate_domain.accession}"
                                        ),
                                    ]
                                )
                            ),
                            conflicting_terms=sorted(
                                set(pair_rule.conflicting_terms)
                            ),
                            ruleset_path=serialized_path,
                            provenance=provenance,
                        )
                    )

    if evidence_records:
        return evidence_records

    candidate_domain = candidate_domains[0]

    return [
        DomainEvidence(
            status=EvidenceStatus.AVAILABLE,
            origin=EvidenceOrigin.ANNOTATION,
            calculation_rule_version=DOMAIN_PAIR_ENGINE_VERSION,
            protein_id=candidate_protein_id,
            source=candidate_domain.source,
            accession=candidate_domain.accession,
            name=candidate_domain.name,
            start=candidate_domain.start,
            end=candidate_domain.end,
            architecture_index=candidate_domain.architecture_index,
            paired_protein_id=query_protein_id,
            is_shared=any(
                query.accession == candidate.accession
                for query in query_domains
                for candidate in candidate_domains
            ),
            pair_matched=False,
            conflicting_terms=["no_matching_domain_pair_rule"],
            ruleset_path=serialized_path,
            provenance=provenance,
        )
    ]