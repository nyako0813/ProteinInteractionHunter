"""Deterministic pairwise phylogenetic profile evaluation."""

from __future__ import annotations

from collections import defaultdict

from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    PhylogeneticProfileEvidence,
)
from protein_interaction_hunter.models.phylogenetic_profile import (
    PhylogeneticProfileObservation,
)

PHYLOGENETIC_PROFILE_ENGINE_VERSION = "mvp1h-phylogenetic-profile-v1"

PhylogeneticProfileIndex = dict[
    str,
    dict[str, PhylogeneticProfileObservation],
]


def build_phylogenetic_profile_index(
    records: list[PhylogeneticProfileObservation],
) -> PhylogeneticProfileIndex:
    index: dict[str, dict[str, PhylogeneticProfileObservation]] = defaultdict(dict)
    for record in records:
        index[record.protein_id][record.species_id] = record
    return {
        protein_id: dict(sorted(species_records.items()))
        for protein_id, species_records in sorted(index.items())
    }


def evaluate_phylogenetic_profile_pair(
    query_protein_id: str,
    candidate_protein_id: str,
    profile_index: PhylogeneticProfileIndex,
    *,
    minimum_shared_species: int,
    minimum_informative_species: int,
    minimum_profile_similarity: float,
) -> PhylogeneticProfileEvidence:
    query_profile = profile_index.get(query_protein_id)
    candidate_profile = profile_index.get(candidate_protein_id)
    provenance = [
        EvidenceProvenance(
            source_name="local_phylogenetic_profile_table",
            source_version=PHYLOGENETIC_PROFILE_ENGINE_VERSION,
            method="deterministic_presence_absence_concordance",
            metadata={
                "minimum_shared_species": minimum_shared_species,
                "minimum_informative_species": minimum_informative_species,
                "minimum_profile_similarity": minimum_profile_similarity,
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]

    if not query_profile or not candidate_profile:
        missing: list[str] = []
        if not query_profile:
            missing.append("query_phylogenetic_profile")
        if not candidate_profile:
            missing.append("candidate_phylogenetic_profile")
        return PhylogeneticProfileEvidence(
            status=EvidenceStatus.MISSING,
            origin=EvidenceOrigin.INFERRED,
            query_protein_id=query_protein_id,
            candidate_protein_id=candidate_protein_id,
            pair_supported=None,
            conflicting_terms=missing,
            calculation_rule_version=PHYLOGENETIC_PROFILE_ENGINE_VERSION,
            source="local_table",
            warnings=["missing_phylogenetic_profile:" + ",".join(missing)],
            provenance=provenance,
        )

    shared_presence = 0
    shared_absence = 0
    discordant = 0
    unknown = 0
    all_species = sorted(set(query_profile) | set(candidate_profile))

    for species_id in all_species:
        query_record = query_profile.get(species_id)
        candidate_record = candidate_profile.get(species_id)
        query_presence = query_record.presence if query_record else None
        candidate_presence = candidate_record.presence if candidate_record else None
        if query_presence is None or candidate_presence is None:
            unknown += 1
        elif query_presence and candidate_presence:
            shared_presence += 1
        elif not query_presence and not candidate_presence:
            shared_absence += 1
        else:
            discordant += 1

    informative = shared_presence + shared_absence + discordant
    similarity = (shared_presence + shared_absence) / informative if informative else None
    pair_supported = (
        informative >= minimum_informative_species
        and shared_presence >= minimum_shared_species
        and similarity is not None
        and similarity >= minimum_profile_similarity
    )

    support_terms: list[str] = []
    conflicting_terms: list[str] = []
    if shared_presence:
        support_terms.append(f"shared_presence:{shared_presence}")
    if pair_supported:
        support_terms.append("profile_support_thresholds_met")
    if informative < minimum_informative_species:
        conflicting_terms.append("insufficient_informative_species")
    if shared_presence < minimum_shared_species:
        conflicting_terms.append("insufficient_shared_presence")
    if similarity is None:
        conflicting_terms.append("profile_similarity_unavailable")
    elif similarity < minimum_profile_similarity:
        conflicting_terms.append("profile_similarity_below_threshold")
    if discordant:
        conflicting_terms.append(f"discordant_species:{discordant}")

    combined_records = [*query_profile.values(), *candidate_profile.values()]
    source_values = sorted({record.source for record in combined_records if record.source})
    record_ids = sorted(
        {record.source_record_id for record in combined_records if record.source_record_id}
    )
    warnings = [f"unknown_species_observations:{unknown}"] if unknown else []

    return PhylogeneticProfileEvidence(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.INFERRED,
        query_protein_id=query_protein_id,
        candidate_protein_id=candidate_protein_id,
        informative_species_count=informative,
        shared_presence_count=shared_presence,
        shared_absence_count=shared_absence,
        discordant_count=discordant,
        unknown_count=unknown,
        profile_similarity=similarity,
        pair_supported=pair_supported,
        support_terms=sorted(support_terms),
        conflicting_terms=sorted(conflicting_terms),
        calculation_rule_version=PHYLOGENETIC_PROFILE_ENGINE_VERSION,
        source="|".join(source_values) or "local_table",
        source_record_id="|".join(record_ids) or None,
        warnings=warnings,
        provenance=provenance,
    )
