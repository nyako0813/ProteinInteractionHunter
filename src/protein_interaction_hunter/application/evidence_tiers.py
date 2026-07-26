"""Deterministic evidence tier classification for MVP-1L."""

from __future__ import annotations

from collections.abc import Sequence

from protein_interaction_hunter.config import (
    EvidenceTiersConfig,
    EvidenceTierThresholdConfig,
)
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceOrigin,
    EvidenceStatus,
    EvidenceTier,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    EvidenceTierResult,
    FusionEvidence,
    IntegratedScore,
    KnownInteractionEvidence,
    OrthologRecord,
    PhylogeneticProfileEvidence,
)
from protein_interaction_hunter.models.protein import CandidateProtein

EVIDENCE_TIER_RULE_VERSION = "mvp1l-evidence-tiers-v1"
HIGH_SPECIFICITY_DEFINITION_VERSION = "mvp1l-high-specificity-v1"
TIER_ORDER = (
    EvidenceTier.TIER_1,
    EvidenceTier.TIER_2,
    EvidenceTier.TIER_3,
    EvidenceTier.TIER_4,
)
TIER_PRIORITY = {
    EvidenceTier.TIER_1: 1,
    EvidenceTier.TIER_2: 2,
    EvidenceTier.TIER_3: 3,
    EvidenceTier.TIER_4: 4,
    EvidenceTier.UNCLASSIFIED: 5,
}
CAP_ORDER = (
    "explicit_conflict",
    "predicted_only",
    "functional_association_only",
)


def _supported(records: Sequence[object], attribute: str) -> bool:
    return any(
        getattr(record, "status", None) is EvidenceStatus.AVAILABLE
        and getattr(record, attribute, None) is True
        for record in records
    )


def _requirements(
    tier: EvidenceTier,
    threshold: EvidenceTierThresholdConfig,
    score: IntegratedScore,
    high_specificity_count: int,
) -> tuple[list[str], list[str]]:
    checks = (
        (
            "minimum_score",
            score.output_score is not None and score.output_score >= threshold.minimum_score,
        ),
        (
            "minimum_categories",
            score.evidence_category_count >= threshold.minimum_categories,
        ),
        (
            "minimum_components",
            score.evidence_component_count >= threshold.minimum_components,
        ),
        (
            "minimum_available_weight",
            score.available_weight >= threshold.minimum_available_weight,
        ),
        (
            "maximum_negative_components",
            score.negative_component_count <= threshold.maximum_negative_components,
        ),
        (
            "minimum_high_specificity_components",
            high_specificity_count >= threshold.minimum_high_specificity_components,
        ),
        (
            "require_high_specificity_evidence",
            not threshold.require_high_specificity_evidence or high_specificity_count > 0,
        ),
    )
    prefix = tier.value
    satisfied = [f"{prefix}:{name}" for name, passed in checks if passed]
    failed = [f"{prefix}:{name}" for name, passed in checks if not passed]
    return satisfied, failed


def _cap(tier: EvidenceTier, cap: EvidenceTier) -> EvidenceTier:
    return cap if TIER_PRIORITY[cap] > TIER_PRIORITY[tier] else tier


def evaluate_evidence_tier(
    query_protein_id: str,
    candidate: CandidateProtein,
    score: IntegratedScore | None,
    config: EvidenceTiersConfig,
    *,
    known_interactions: KnownInteractionEvidence | None = None,
    fusion: FusionEvidence | None = None,
    orthology: Sequence[OrthologRecord] = (),
    phylogenetic_profile: PhylogeneticProfileEvidence | None = None,
) -> EvidenceTierResult:
    """Classify one pair without modifying its score, rank, evidence, or disposition."""
    known_records = [known_interactions] if known_interactions else []
    fusion_records = [fusion] if fusion else []
    profile_records = [phylogenetic_profile] if phylogenetic_profile else []
    direct = _supported(known_records, "direct_interaction_supported")
    physical = _supported(known_records, "physical_interaction_supported")
    functional = _supported(known_records, "functional_association_supported")
    fusion_supported = _supported(fusion_records, "pair_supported")
    orthology_supported = _supported(orthology, "pair_supported")
    profile_supported = _supported(profile_records, "pair_supported")

    high_specificity: list[str] = []
    if direct:
        high_specificity.append("direct_interaction_evidence")
    elif physical:
        high_specificity.append("physical_association_evidence")
    if fusion_supported:
        high_specificity.append("fusion_association_evidence")
    high_specificity = sorted(set(high_specificity))

    types = {
        interaction_type
        for record in known_records
        if record.status is EvidenceStatus.AVAILABLE
        for interaction_type in record.interaction_types
    }
    predicted_only = bool(types) and types <= {"predicted"} and not high_specificity
    functional_only = functional and not (direct or physical or fusion_supported)

    empty = {
        "formal_score": score.output_score if score else None,
        "rank": score.rank if score else None,
        "sufficient_evidence": score.sufficient_evidence if score else False,
        "evidence_category_count": score.evidence_category_count if score else 0,
        "evidence_component_count": score.evidence_component_count if score else 0,
        "available_weight": score.available_weight if score else 0.0,
        "positive_component_count": score.positive_component_count if score else 0,
        "neutral_component_count": score.neutral_component_count if score else 0,
        "negative_component_count": score.negative_component_count if score else 0,
    }
    negative_components = (
        sorted(
            component.component_name
            for component in score.component_scores
            if component.applied and component.direction == "negative"
        )
        if score
        else []
    )
    explicit_conflict = bool(negative_components)
    support_terms = sorted(
        set(score.support_terms if score else [])
        | {
            term
            for record in (*known_records, *fusion_records, *orthology, *profile_records)
            for term in getattr(record, "support_terms", [])
        }
    )
    conflicting_terms = sorted(
        set(score.conflicting_terms if score else [])
        | {
            term
            for record in (*known_records, *fusion_records, *orthology, *profile_records)
            for term in getattr(record, "conflicting_terms", [])
        }
    )
    provenance = [
        EvidenceProvenance(
            source_name="evidence_tier_engine",
            source_version=config.rule_version,
            method="quantitative_base_tier_then_semantic_caps",
            metadata={
                "high_specificity_definition_version": HIGH_SPECIFICITY_DEFINITION_VERSION,
                "scoring_rule_version_dependency": (
                    score.calculation_rule_version if score else None
                ),
                "score_recalculated": False,
                "rank_recalculated": False,
            },
        )
    ]

    ineligible: str | None = None
    status = EvidenceStatus.AVAILABLE
    if score is None:
        ineligible = "scoring_result_missing"
        status = EvidenceStatus.MISSING
    elif score.status is not EvidenceStatus.AVAILABLE:
        ineligible = f"scoring_status_{score.status.value}"
        status = EvidenceStatus.MISSING
    elif candidate.protein_id == query_protein_id:
        ineligible = "self_pair"
        status = EvidenceStatus.NOT_APPLICABLE
    elif candidate.disposition is CandidateDisposition.EXCLUDED:
        ineligible = "candidate_excluded"
        status = EvidenceStatus.NOT_APPLICABLE
    elif not score.sufficient_evidence:
        ineligible = "insufficient_evidence"
    elif score.output_score is None:
        ineligible = "formal_score_missing"
    elif score.rank is None:
        ineligible = "rank_missing"

    common = dict(
        status=status,
        origin=EvidenceOrigin.INFERRED,
        query_protein_id=query_protein_id,
        candidate_protein_id=candidate.protein_id,
        high_specificity_component_count=len(high_specificity),
        high_specificity_components=high_specificity,
        direct_interaction_supported=direct,
        physical_interaction_supported=physical,
        functional_association_supported=functional,
        fusion_supported=fusion_supported,
        orthology_supported=orthology_supported,
        phylogenetic_profile_supported=profile_supported,
        explicit_conflict_present=explicit_conflict,
        predicted_only=predicted_only,
        functional_association_only=functional_only,
        support_terms=support_terms,
        conflicting_terms=conflicting_terms,
        calculation_rule_version=config.rule_version,
        provenance=provenance,
        **empty,
    )
    if ineligible:
        return EvidenceTierResult.model_validate(
            {
                **common,
                "assigned_tier": EvidenceTier.UNCLASSIFIED,
                "base_tier": EvidenceTier.UNCLASSIFIED,
                "tier_eligible": False,
                "failed_requirements": [ineligible],
            }
        )

    assert score is not None
    base_tier = EvidenceTier.UNCLASSIFIED
    satisfied: list[str] = []
    failed: list[str] = []
    for tier in TIER_ORDER:
        threshold = getattr(config, tier.value)
        tier_satisfied, tier_failed = _requirements(tier, threshold, score, len(high_specificity))
        satisfied.extend(tier_satisfied)
        if not tier_failed:
            base_tier = tier
            break
        failed.extend(tier_failed)

    if base_tier is EvidenceTier.UNCLASSIFIED:
        return EvidenceTierResult.model_validate(
            {
                **common,
                "assigned_tier": base_tier,
                "base_tier": base_tier,
                "tier_eligible": True,
                "satisfied_requirements": satisfied,
                "failed_requirements": failed,
            }
        )

    cap_values = {
        "explicit_conflict": EvidenceTier(config.explicit_conflict_tier_cap),
        "predicted_only": EvidenceTier(config.predicted_only_tier_cap),
        "functional_association_only": EvidenceTier(config.functional_association_only_tier_cap),
    }
    cap_conditions = {
        "explicit_conflict": explicit_conflict,
        "predicted_only": predicted_only,
        "functional_association_only": functional_only,
    }
    assigned: EvidenceTier = base_tier
    applied_caps: list[str] = []
    for cap_name in CAP_ORDER:
        if not cap_conditions[cap_name]:
            continue
        if TIER_PRIORITY[cap_values[cap_name]] > TIER_PRIORITY[base_tier]:
            applied_caps.append(f"{cap_name}:{cap_values[cap_name].value}")
        assigned = _cap(assigned, cap_values[cap_name])

    return EvidenceTierResult.model_validate(
        {
            **common,
            "assigned_tier": assigned,
            "base_tier": base_tier,
            "tier_eligible": True,
            "applied_tier_caps": applied_caps,
            "satisfied_requirements": satisfied,
            "failed_requirements": failed,
        }
    )
