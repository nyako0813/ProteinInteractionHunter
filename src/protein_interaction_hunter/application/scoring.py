# Deterministic integrated evidence scoring for MVP-1K.

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence
from typing import cast

from protein_interaction_hunter.config import ScoringConfig
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceOrigin,
    EvidenceStatus,
    IdentifierMatchStatus,
    OperonProxyStatus,
)
from protein_interaction_hunter.models.evidence import (
    AppliedScoringPenalty,
    DomainEvidence,
    EvidenceProvenance,
    FunctionalEvidence,
    FusionEvidence,
    GenomeContextEvidence,
    IntegratedScore,
    KnownInteractionEvidence,
    LocalizationEvidence,
    OperonEvidence,
    OrthologRecord,
    PhylogeneticProfileEvidence,
    ScoreCategory,
    ScoreComponent,
    ScoringDirection,
)
from protein_interaction_hunter.models.protein import CandidateProtein
from protein_interaction_hunter.models.scoring import CandidateScore

SCORING_RULE_VERSION = "mvp1k-integrated-scoring-v1"
CALCULATION_PRECISION = 12
COMPONENT_ORDER = (
    "genome_context",
    "operon_proxy",
    "domain_pair",
    "functional_complementarity",
    "localization",
    "orthology",
    "phylogenetic_profile",
    "fusion",
    "known_interactions",
)
COMPONENT_CATEGORIES = {
    "genome_context": "genomic_context",
    "operon_proxy": "genomic_context",
    "domain_pair": "functional_annotation",
    "functional_complementarity": "functional_annotation",
    "localization": "cellular_compatibility",
    "orthology": "evolutionary",
    "phylogenetic_profile": "evolutionary",
    "fusion": "direct_interaction",
    "known_interactions": "direct_interaction",
}


def _status(records: Sequence[object]) -> EvidenceStatus:
    if not records:
        return EvidenceStatus.NOT_RUN
    statuses = [
        cast(EvidenceStatus, getattr(record, "status"))  # noqa: B009
        for record in records
    ]
    if EvidenceStatus.AVAILABLE in statuses:
        return EvidenceStatus.AVAILABLE
    return statuses[0]


def _terms(records: Sequence[object], attribute: str) -> list[str]:
    return sorted({str(value) for record in records for value in getattr(record, attribute, [])})


def _versions(records: Sequence[object]) -> str | None:
    values = sorted(
        {
            str(value)
            for record in records
            if (value := getattr(record, "calculation_rule_version", None))
        }
    )
    return "|".join(values) or None


def _safe_value(value: float | None) -> float | None:
    if value is None or not math.isfinite(value) or value < -1.0 or value > 1.0:
        return None
    return value


def _component(
    name: str,
    status: EvidenceStatus,
    raw_value: str | float | int | bool | None,
    normalized_value: float | None,
    config: ScoringConfig,
    *,
    support_terms: Sequence[str] = (),
    conflicting_terms: Sequence[str] = (),
    source_rule_version: str | None = None,
    warnings: Sequence[str] = (),
) -> ScoreComponent:
    weight = getattr(config.weights, name)
    normalized = _safe_value(normalized_value)
    exclusion_reason: str | None = None
    if status is not EvidenceStatus.AVAILABLE:
        exclusion_reason = f"evidence_{status.value}"
    elif normalized is None:
        exclusion_reason = "invalid_or_unavailable_value"
    elif weight == 0:
        exclusion_reason = "weight_zero"
    applied = exclusion_reason is None
    direction: ScoringDirection
    if normalized is None:
        direction = "unknown"
    elif normalized > 0:
        direction = "positive"
    elif normalized < 0:
        direction = "negative"
    else:
        direction = "neutral"
    return ScoreComponent(
        component_name=name,
        category_name=COMPONENT_CATEGORIES[name],
        evidence_status=status,
        raw_value=raw_value,
        normalized_value=normalized,
        configured_weight=weight,
        direction=direction,
        applied=applied,
        exclusion_reason=exclusion_reason,
        support_terms=sorted(set(support_terms)),
        conflicting_terms=sorted(set(conflicting_terms)),
        source_rule_version=source_rule_version,
        warnings=sorted(set(warnings)),
    )


def normalize_genome_context(
    evidence: GenomeContextEvidence | None, config: ScoringConfig
) -> ScoreComponent:
    if evidence is None:
        return _component("genome_context", EvidenceStatus.NOT_RUN, None, None, config)
    value: float | None = None
    raw: str | int | bool | None = evidence.distance_bp
    if evidence.status is EvidenceStatus.AVAILABLE:
        if evidence.coordinate_position and evidence.coordinate_position.value == "same_feature":
            value = 0.0
            raw = "same_feature"
        elif evidence.overlap_bp and evidence.overlap_bp > 0:
            value = 0.0
            raw = f"overlap:{evidence.overlap_bp}"
        elif evidence.same_contig is True and evidence.within_neighborhood_window is True:
            value = 0.8 if evidence.intervening_gene_count == 0 else 0.5
        elif evidence.same_contig is True:
            value = 0.1
    return _component(
        "genome_context",
        evidence.status,
        raw,
        value,
        config,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    )


def normalize_operon(evidence: OperonEvidence | None, config: ScoringConfig) -> ScoreComponent:
    if evidence is None:
        return _component("operon_proxy", EvidenceStatus.NOT_RUN, None, None, config)
    mapping = {
        OperonProxyStatus.SUPPORTED: 1.0,
        OperonProxyStatus.PARTIAL_SUPPORT: 0.3,
        OperonProxyStatus.NOT_SUPPORTED: 0.0,
    }
    value = mapping.get(evidence.proxy_status) if evidence.proxy_status is not None else None
    return _component(
        "operon_proxy",
        evidence.status,
        evidence.proxy_status.value if evidence.proxy_status else None,
        value,
        config,
        support_terms=evidence.supporting_conditions,
        conflicting_terms=evidence.conflicting_conditions,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    )


def normalize_domain(records: Sequence[DomainEvidence], config: ScoringConfig) -> ScoreComponent:
    status = _status(records)
    value = (
        None
        if status is not EvidenceStatus.AVAILABLE
        else (1.0 if any(record.pair_matched is True for record in records) else 0.0)
    )
    return _component(
        "domain_pair",
        status,
        any(record.pair_matched is True for record in records) if records else None,
        value,
        config,
        support_terms=_terms(records, "support_terms"),
        conflicting_terms=_terms(records, "conflicting_terms"),
        source_rule_version=_versions(records),
        warnings=_terms(records, "warnings"),
    )


def normalize_functional(
    records: Sequence[FunctionalEvidence], config: ScoringConfig
) -> ScoreComponent:
    status = _status(records)
    value = (
        None
        if status is not EvidenceStatus.AVAILABLE
        else (1.0 if any(record.matched is True for record in records) else 0.0)
    )
    return _component(
        "functional_complementarity",
        status,
        any(record.matched is True for record in records) if records else None,
        value,
        config,
        support_terms=_terms(records, "support_terms"),
        conflicting_terms=_terms(records, "conflicting_terms"),
        source_rule_version=_versions(records),
        warnings=_terms(records, "warnings"),
    )


def normalize_localization(
    evidence: LocalizationEvidence | None, config: ScoringConfig
) -> tuple[ScoreComponent, list[AppliedScoringPenalty]]:
    if evidence is None:
        return _component("localization", EvidenceStatus.NOT_RUN, None, None, config), []
    penalties: list[AppliedScoringPenalty] = []
    value: float | None = None
    if evidence.status is EvidenceStatus.AVAILABLE:
        if evidence.compatibility is True:
            value = 1.0
        elif evidence.compatibility is False:
            value = -min(1.0, config.penalties.contradictory_evidence)
            penalties.append(
                AppliedScoringPenalty(
                    penalty_name="contradictory_evidence",
                    component_name="localization",
                    configured_value=config.penalties.contradictory_evidence,
                    explanation="explicit_localization_incompatibility",
                )
            )
        elif evidence.compatibility is None:
            value = 0.0
    return _component(
        "localization",
        evidence.status,
        evidence.compatibility,
        value,
        config,
        support_terms=evidence.matched_terms,
        conflicting_terms=evidence.conflicting_terms,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    ), penalties


def normalize_orthology(
    records: Sequence[OrthologRecord], config: ScoringConfig
) -> tuple[ScoreComponent, list[AppliedScoringPenalty]]:
    status = _status(records)
    penalties: list[AppliedScoringPenalty] = []
    value: float | None = None
    if status is EvidenceStatus.AVAILABLE:
        if any(record.pair_supported is True for record in records):
            value = 1.0
            if any(record.paralog_ambiguity for record in records):
                reduction = min(1.0, config.penalties.ambiguous_mapping)
                value *= 1.0 - reduction
                penalties.append(
                    AppliedScoringPenalty(
                        penalty_name="ambiguous_mapping",
                        component_name="orthology",
                        configured_value=config.penalties.ambiguous_mapping,
                        explanation="paralog_ambiguity_reduces_positive_strength",
                    )
                )
        else:
            value = 0.0
    return _component(
        "orthology",
        status,
        any(record.pair_supported is True for record in records) if records else None,
        value,
        config,
        support_terms=_terms(records, "support_terms"),
        conflicting_terms=_terms(records, "conflicting_terms"),
        source_rule_version=_versions(records),
        warnings=_terms(records, "warnings"),
    ), penalties


def normalize_phylogenetic_profile(
    evidence: PhylogeneticProfileEvidence | None, config: ScoringConfig
) -> tuple[ScoreComponent, list[AppliedScoringPenalty]]:
    if evidence is None:
        return _component("phylogenetic_profile", EvidenceStatus.NOT_RUN, None, None, config), []
    penalties: list[AppliedScoringPenalty] = []
    value: float | None = None
    similarity = evidence.profile_similarity
    if (
        evidence.status is EvidenceStatus.AVAILABLE
        and "insufficient_informative_species" not in evidence.conflicting_terms
        and similarity is not None
        and evidence.shared_presence_count is not None
    ):
        concordant = evidence.shared_presence_count + (evidence.shared_absence_count or 0)
        presence_fraction = evidence.shared_presence_count / concordant if concordant else 0.0
        value = similarity * presence_fraction
        informative = evidence.informative_species_count or 0
        discordant = evidence.discordant_count or 0
        if informative >= 3 and discordant / informative >= 0.8:
            value = -min(1.0, config.penalties.contradictory_evidence)
            penalties.append(
                AppliedScoringPenalty(
                    penalty_name="contradictory_evidence",
                    component_name="phylogenetic_profile",
                    configured_value=config.penalties.contradictory_evidence,
                    explanation="strong_profile_discordance",
                )
            )
    return _component(
        "phylogenetic_profile",
        evidence.status,
        similarity,
        value,
        config,
        support_terms=evidence.support_terms,
        conflicting_terms=evidence.conflicting_terms,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    ), penalties


def normalize_fusion(
    evidence: FusionEvidence | None, config: ScoringConfig
) -> tuple[ScoreComponent, list[AppliedScoringPenalty]]:
    if evidence is None:
        return _component("fusion", EvidenceStatus.NOT_RUN, None, None, config), []
    penalties: list[AppliedScoringPenalty] = []
    value: float | None = None
    if evidence.status is EvidenceStatus.AVAILABLE:
        if evidence.pair_supported is True:
            value = 1.0
        elif "excessive_component_overlap" in evidence.conflicting_terms:
            value = -min(1.0, config.penalties.contradictory_evidence)
            penalties.append(
                AppliedScoringPenalty(
                    penalty_name="contradictory_evidence",
                    component_name="fusion",
                    configured_value=config.penalties.contradictory_evidence,
                    explanation="excessive_fusion_component_overlap",
                )
            )
        else:
            value = 0.0
    return _component(
        "fusion",
        evidence.status,
        evidence.pair_supported,
        value,
        config,
        support_terms=evidence.support_terms,
        conflicting_terms=evidence.conflicting_terms,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    ), penalties


def normalize_known_interactions(
    evidence: KnownInteractionEvidence | None, config: ScoringConfig
) -> ScoreComponent:
    if evidence is None:
        return _component("known_interactions", EvidenceStatus.NOT_RUN, None, None, config)
    value: float | None = None
    if evidence.status is EvidenceStatus.AVAILABLE:
        if evidence.direct_interaction_supported is True:
            value = 1.0
        elif evidence.physical_interaction_supported is True:
            value = 0.9
        elif evidence.functional_association_supported is True:
            value = 0.5
        elif evidence.pair_supported is True:
            value = 0.4
        elif "predicted" in evidence.interaction_types:
            value = 0.2
        else:
            value = 0.0
    return _component(
        "known_interactions",
        evidence.status,
        evidence.pair_supported,
        value,
        config,
        support_terms=evidence.support_terms,
        conflicting_terms=evidence.conflicting_terms,
        source_rule_version=evidence.calculation_rule_version,
        warnings=evidence.warnings,
    )


def score_pair(
    query_protein_id: str,
    candidate: CandidateProtein,
    config: ScoringConfig,
    *,
    genome_context: GenomeContextEvidence | None = None,
    operon: OperonEvidence | None = None,
    domains: Sequence[DomainEvidence] = (),
    functional: Sequence[FunctionalEvidence] = (),
    localization: LocalizationEvidence | None = None,
    orthology: Sequence[OrthologRecord] = (),
    phylogenetic_profile: PhylogeneticProfileEvidence | None = None,
    fusion: FusionEvidence | None = None,
    known_interactions: KnownInteractionEvidence | None = None,
) -> IntegratedScore:
    localization_component, localization_penalties = normalize_localization(localization, config)
    orthology_component, orthology_penalties = normalize_orthology(orthology, config)
    profile_component, profile_penalties = normalize_phylogenetic_profile(
        phylogenetic_profile, config
    )
    fusion_component, fusion_penalties = normalize_fusion(fusion, config)
    components = [
        normalize_genome_context(genome_context, config),
        normalize_operon(operon, config),
        normalize_domain(domains, config),
        normalize_functional(functional, config),
        localization_component,
        orthology_component,
        profile_component,
        fusion_component,
        normalize_known_interactions(known_interactions, config),
    ]
    penalties = localization_penalties + orthology_penalties + profile_penalties + fusion_penalties
    if (
        candidate.identifier_match_status is IdentifierMatchStatus.AMBIGUOUS_MATCH
        and config.penalties.ambiguous_mapping > 0
    ):
        factor = max(0.0, 1.0 - config.penalties.ambiguous_mapping)
        components = [
            item.model_copy(
                update={
                    "normalized_value": (item.normalized_value or 0.0) * factor,
                    "direction": (
                        "positive" if (item.normalized_value or 0.0) * factor > 0 else "neutral"
                    ),
                }
            )
            if item.applied and item.direction == "positive"
            else item
            for item in components
        ]
        penalties.append(
            AppliedScoringPenalty(
                penalty_name="ambiguous_mapping",
                component_name="all_positive_components",
                configured_value=config.penalties.ambiguous_mapping,
                explanation="candidate_identifier_mapping_ambiguity",
            )
        )
    by_category: dict[str, list[ScoreComponent]] = defaultdict(list)
    for component in components:
        if component.applied:
            by_category[component.category_name].append(component)
    caps = config.category_caps.model_dump()
    finalized: list[ScoreComponent] = []
    for component in components:
        if not component.applied:
            finalized.append(component)
            continue
        category_components = by_category[component.category_name]
        configured_total = sum(item.configured_weight for item in category_components)
        scale = min(1.0, caps[component.category_name] / configured_total)
        effective = round(component.configured_weight * scale, CALCULATION_PRECISION)
        contribution = round(
            (component.normalized_value or 0.0) * effective,
            CALCULATION_PRECISION,
        )
        finalized.append(
            component.model_copy(
                update={
                    "effective_weight": effective,
                    "weighted_contribution": contribution,
                }
            )
        )
    category_scores: list[ScoreCategory] = []
    for category in sorted(caps):
        values = [item for item in finalized if item.applied and item.category_name == category]
        weight = round(sum(item.effective_weight for item in values), CALCULATION_PRECISION)
        raw = round(
            sum(item.weighted_contribution or 0.0 for item in values),
            CALCULATION_PRECISION,
        )
        category_scores.append(
            ScoreCategory(
                category_name=category,
                raw_weighted_sum=raw,
                available_weight=weight,
                normalized_score=(round(raw / weight, CALCULATION_PRECISION) if weight else 0.0),
                configured_cap=caps[category],
            )
        )
    applied = [item for item in finalized if item.applied]
    available_weight = round(sum(item.effective_weight for item in applied), CALCULATION_PRECISION)
    raw_sum = round(
        sum(item.weighted_contribution or 0.0 for item in applied),
        CALCULATION_PRECISION,
    )
    normalized = (
        round(
            max(0.0, min(1.0, raw_sum / available_weight)),
            CALCULATION_PRECISION,
        )
        if available_weight
        else None
    )
    provisional = (
        round(normalized * config.output_scale, CALCULATION_PRECISION)
        if normalized is not None
        else None
    )
    categories = len({item.category_name for item in applied})
    sufficient = (
        available_weight >= config.minimum_evidence_weight
        and categories >= config.minimum_evidence_categories
    )
    output = provisional if sufficient else None
    status = EvidenceStatus.AVAILABLE if available_weight else EvidenceStatus.MISSING
    warnings: list[str] = []
    conflicts: list[str] = []
    support: list[str] = []
    if not available_weight:
        warnings.append("no_available_scoring_components")
        conflicts.append("no_available_evidence")
    if not sufficient:
        warnings.append("insufficient_evidence_for_formal_score")
        conflicts.append("insufficient_evidence")
    else:
        support.append("sufficient_evidence_for_formal_score")
    if any(item.direction == "positive" for item in applied):
        support.append("positive_evidence_components")
    if any(item.direction == "negative" for item in applied):
        conflicts.append("explicit_negative_evidence_component")
    provenance = [
        EvidenceProvenance(
            source_name="integrated_scoring_engine",
            source_version=config.rule_version,
            method="category_capped_missing_excluded_weighted_support_strength",
            metadata={
                "missing_policy": config.missing_policy,
                "output_scale": config.output_scale,
                "tie_precision": config.tie_precision,
                "evidence_tier_propagation": False,
                "relationship_propagation": False,
                "disposition_propagation": False,
            },
        )
    ]
    return IntegratedScore(
        status=status,
        origin=EvidenceOrigin.INFERRED,
        query_protein_id=query_protein_id,
        candidate_protein_id=candidate.protein_id,
        raw_weighted_sum=raw_sum,
        available_weight=available_weight,
        normalized_score=normalized,
        provisional_score=provisional,
        output_score=output,
        evidence_category_count=categories,
        evidence_component_count=len(applied),
        positive_component_count=sum(item.direction == "positive" for item in applied),
        neutral_component_count=sum(item.direction == "neutral" for item in applied),
        negative_component_count=sum(item.direction == "negative" for item in applied),
        sufficient_evidence=sufficient,
        component_scores=finalized,
        category_scores=category_scores,
        applied_weights={item.component_name: item.effective_weight for item in finalized},
        applied_penalties=penalties,
        support_terms=sorted(support),
        conflicting_terms=sorted(conflicts),
        calculation_rule_version=config.rule_version,
        warnings=sorted(warnings),
        provenance=provenance,
    )


def rank_scores(
    scores: dict[tuple[str, str], IntegratedScore],
    candidates: Sequence[CandidateProtein],
    tie_precision: int,
) -> dict[tuple[str, str], IntegratedScore]:
    candidate_by_pair = {(item.query_id, item.protein_id): item for item in candidates}
    ranked = dict(scores)
    by_query: dict[str, list[tuple[tuple[str, str], IntegratedScore]]] = defaultdict(list)
    for pair, score in scores.items():
        candidate = candidate_by_pair[pair]
        if (
            candidate.disposition is not CandidateDisposition.EXCLUDED
            and score.sufficient_evidence
            and score.output_score is not None
        ):
            by_query[pair[0]].append((pair, score))
    for items in by_query.values():
        items.sort(
            key=lambda item: (
                -round(item[1].output_score or 0.0, tie_precision),
                -item[1].evidence_category_count,
                -item[1].available_weight,
                item[0][1],
            )
        )
        dense_rank = 0
        previous: float | None = None
        for pair, score in items:
            tie_value = round(score.output_score or 0.0, tie_precision)
            if previous is None or tie_value != previous:
                dense_rank += 1
                previous = tie_value
            ranked[pair] = score.model_copy(update={"rank": dense_rank, "tied_rank": dense_rank})
    return ranked


def candidate_score_from_integrated(score: IntegratedScore) -> CandidateScore:
    categories = {item.category_name: item for item in score.category_scores}

    def positive(name: str) -> float | None:
        item = categories.get(name)
        return max(0.0, item.normalized_score) if item and item.available_weight > 0 else None

    negative_sum = sum(
        abs(item.weighted_contribution or 0.0)
        for item in score.component_scores
        if item.applied and item.direction == "negative"
    )
    completeness_denominator = sum(item.configured_cap for item in score.category_scores)
    return CandidateScore(
        physical_interaction_score=positive("direct_interaction"),
        functional_association_score=positive("functional_annotation"),
        gene_context_score=positive("genomic_context"),
        evolutionary_coupling_score=positive("evolutionary"),
        annotation_confidence_score=positive("cellular_compatibility"),
        contradiction_penalty=(
            negative_sum / score.available_weight if score.available_weight else None
        ),
        evidence_completeness=(
            score.available_weight / completeness_denominator if completeness_denominator else None
        ),
        total_ranking_score=score.normalized_score if score.sufficient_evidence else None,
        calculation_trace=[
            f"raw_weighted_sum={score.raw_weighted_sum}",
            f"available_weight={score.available_weight}",
            f"normalized_score={score.normalized_score}",
            f"rule_version={score.calculation_rule_version}",
        ],
        warnings=score.warnings,
    )
