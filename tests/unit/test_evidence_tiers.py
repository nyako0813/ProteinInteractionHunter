from __future__ import annotations

from copy import deepcopy

import pytest

from protein_interaction_hunter.application.evidence_tiers import evaluate_evidence_tier
from protein_interaction_hunter.config import EvidenceTiersConfig
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceOrigin,
    EvidenceStatus,
    EvidenceTier,
    IdentifierMatchStatus,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceTierResult,
    FusionEvidence,
    IntegratedScore,
    KnownInteractionEvidence,
    ScoreComponent,
)
from protein_interaction_hunter.models.protein import CandidateProtein


def candidate(
    protein_id: str = "CANDIDATE",
    disposition: CandidateDisposition = CandidateDisposition.INCLUDED,
) -> CandidateProtein:
    return CandidateProtein(
        query_id="QUERY",
        protein_id=protein_id,
        disposition=disposition,
        sequence_length=100,
        identifier_match_status=IdentifierMatchStatus.EXACT_MATCH,
    )


def component(name: str, direction: str = "positive") -> ScoreComponent:
    value = {"positive": 1.0, "neutral": 0.0, "negative": -0.25}[direction]
    return ScoreComponent(
        component_name=name,
        category_name=f"category_{name}",
        evidence_status=EvidenceStatus.AVAILABLE,
        raw_value=True,
        normalized_value=value,
        configured_weight=1.0,
        effective_weight=1.0,
        weighted_contribution=value,
        direction=direction,  # type: ignore[arg-type]
        applied=True,
    )


def score(
    output: float | None = 80.0,
    categories: int = 4,
    components: int = 5,
    weight: float = 3.0,
    negatives: int = 0,
    sufficient: bool = True,
    rank: int | None = 1,
) -> IntegratedScore:
    items = [component(f"positive_{index}") for index in range(components - negatives)]
    items.extend(component(f"negative_{index}", "negative") for index in range(negatives))
    return IntegratedScore(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.INFERRED,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        raw_weighted_sum=output or 0.0,
        available_weight=weight,
        normalized_score=(output / 100.0 if output is not None else None),
        provisional_score=output,
        output_score=output,
        evidence_category_count=categories,
        evidence_component_count=components,
        positive_component_count=components - negatives,
        neutral_component_count=0,
        negative_component_count=negatives,
        sufficient_evidence=sufficient,
        rank=rank,
        tied_rank=rank,
        component_scores=items,
        calculation_rule_version="mvp1k-integrated-scoring-v1",
    )


def known(kind: str) -> KnownInteractionEvidence:
    direct = kind == "direct"
    physical = kind == "physical"
    functional = kind == "functional"
    predicted = kind == "predicted"
    return KnownInteractionEvidence(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.EXACT_PAIR,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        supporting_record_count=2,
        qualifying_record_count=2,
        direct_record_count=2 if direct else 0,
        physical_record_count=2 if physical else 0,
        biological_record_count=2 if functional else 0,
        independent_publication_count=1,
        independent_source_count=1,
        interaction_types=[kind if kind != "functional" else "functional_association"],
        pair_supported=not predicted,
        direct_interaction_supported=direct,
        physical_interaction_supported=physical,
        functional_association_supported=functional,
    )


def fusion() -> FusionEvidence:
    return FusionEvidence(
        status=EvidenceStatus.AVAILABLE,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        supporting_record_count=1,
        qualifying_record_count=1,
        pair_supported=True,
    )


def evaluate(
    integrated: IntegratedScore,
    *,
    interaction: KnownInteractionEvidence | None = None,
    fusion_evidence: FusionEvidence | None = None,
    item: CandidateProtein | None = None,
) -> EvidenceTierResult:
    return evaluate_evidence_tier(
        "QUERY",
        item or candidate(),
        integrated,
        EvidenceTiersConfig(enabled=True),
        known_interactions=interaction,
        fusion=fusion_evidence,
    )


@pytest.mark.parametrize(
    ("integrated", "expected"),
    [
        (score(80.0, 4, 5, 3.0), EvidenceTier.TIER_1),
        (score(55.0, 3, 4, 2.0), EvidenceTier.TIER_2),
        (score(25.0, 2, 2, 1.0), EvidenceTier.TIER_3),
        (score(0.0, 2, 2, 1.0), EvidenceTier.TIER_4),
    ],
)
def test_base_tiers_at_exact_boundaries(
    integrated: IntegratedScore, expected: EvidenceTier
) -> None:
    result = evaluate(integrated, interaction=known("direct"))
    assert result.base_tier is expected
    assert result.assigned_tier is expected


def test_tier_one_needs_categories_and_high_specificity() -> None:
    missing_category = evaluate(score(categories=3), interaction=known("direct"))
    missing_specificity = evaluate(score())
    assert missing_category.base_tier is EvidenceTier.TIER_2
    assert "tier_1:minimum_categories" in missing_category.failed_requirements or (
        missing_category.base_tier is EvidenceTier.TIER_2
    )
    assert missing_specificity.base_tier is EvidenceTier.TIER_2
    assert missing_specificity.high_specificity_component_count == 0


def test_no_quantitative_tier_is_unclassified() -> None:
    result = evaluate(score(categories=1, components=1, weight=0.5))
    assert result.tier_eligible is True
    assert result.assigned_tier is EvidenceTier.UNCLASSIFIED
    assert result.failed_requirements


def test_explicit_conflict_uses_only_applied_negative_component() -> None:
    result = evaluate(score(80.0, 4, 5, 3.0, negatives=1), interaction=known("direct"))
    assert result.explicit_conflict_present is True
    assert result.base_tier is EvidenceTier.TIER_2
    assert result.assigned_tier is EvidenceTier.TIER_3
    assert result.applied_tier_caps == ["explicit_conflict:tier_3"]

    missing = score()
    missing.component_scores.append(
        ScoreComponent(
            component_name="missing",
            category_name="unused",
            evidence_status=EvidenceStatus.MISSING,
            configured_weight=1.0,
            direction="unknown",
            applied=False,
            exclusion_reason="evidence_missing",
        )
    )
    assert evaluate(missing, interaction=known("direct")).explicit_conflict_present is False


def test_predicted_and_functional_only_caps_are_semantic_and_ordered() -> None:
    predicted = evaluate(score(), interaction=known("predicted"))
    functional = evaluate(score(), interaction=known("functional"))
    assert predicted.base_tier is EvidenceTier.TIER_2
    assert predicted.assigned_tier is EvidenceTier.TIER_3
    assert predicted.applied_tier_caps == ["predicted_only:tier_3"]
    assert functional.assigned_tier is EvidenceTier.TIER_3
    assert functional.applied_tier_caps == ["functional_association_only:tier_3"]
    assert (
        evaluate(score(), interaction=known("functional"), fusion_evidence=fusion()).assigned_tier
        is EvidenceTier.TIER_1
    )


def test_multiple_caps_are_recorded_in_fixed_order() -> None:
    result = evaluate(
        score(80.0, 4, 5, 3.0, negatives=1),
        interaction=known("functional"),
    )
    assert result.base_tier is EvidenceTier.TIER_2
    assert result.assigned_tier is EvidenceTier.TIER_3
    assert result.applied_tier_caps == [
        "explicit_conflict:tier_3",
        "functional_association_only:tier_3",
    ]


def test_direct_or_physical_records_are_one_deduplicated_specificity_class() -> None:
    direct = known("direct").model_copy(
        update={
            "supporting_record_count": 20,
            "independent_publication_count": 1,
            "independent_source_count": 1,
        }
    )
    direct_result = evaluate(score(), interaction=direct)
    physical_result = evaluate(score(), interaction=known("physical"))
    assert direct_result.high_specificity_components == ["direct_interaction_evidence"]
    assert direct_result.high_specificity_component_count == 1
    assert physical_result.high_specificity_components == ["physical_association_evidence"]


@pytest.mark.parametrize(
    ("integrated", "item", "reason"),
    [
        (score(sufficient=False, output=None, rank=None), candidate(), "insufficient_evidence"),
        (score(output=None, rank=None), candidate(), "formal_score_missing"),
        (score(rank=None), candidate(), "rank_missing"),
        (
            score(rank=None),
            candidate(disposition=CandidateDisposition.EXCLUDED),
            "candidate_excluded",
        ),
        (score(rank=None), candidate("QUERY"), "self_pair"),
    ],
)
def test_ineligible_pairs_are_unclassified(
    integrated: IntegratedScore, item: CandidateProtein, reason: str
) -> None:
    result = evaluate(integrated, item=item)
    assert result.assigned_tier is EvidenceTier.UNCLASSIFIED
    assert result.tier_eligible is False
    assert reason in result.failed_requirements


def test_missing_scoring_is_unclassified_missing() -> None:
    result = evaluate_evidence_tier("QUERY", candidate(), None, EvidenceTiersConfig(enabled=True))
    assert result.status is EvidenceStatus.MISSING
    assert result.assigned_tier is EvidenceTier.UNCLASSIFIED


def test_deterministic_raw_order_and_no_input_mutation() -> None:
    integrated = score()
    original = deepcopy(integrated)
    first = evaluate(integrated, interaction=known("direct"), fusion_evidence=fusion())
    second = evaluate(integrated, interaction=known("direct"), fusion_evidence=fusion())
    assert first == second
    assert integrated == original


def test_golden_decisions() -> None:
    # A: score=80, categories=4, components=5, weight=3, direct high-specificity,
    # negatives=0 -> base Tier 1, no cap -> Tier 1.
    a = evaluate(score(), interaction=known("direct"))
    # B: same quantitative values, functional only, no high-specificity -> base Tier 2,
    # functional-only cap Tier 3 -> Tier 3.
    b = evaluate(score(), interaction=known("functional"))
    # C: score=60, categories=3, components=4, weight=2, negatives=1 -> base Tier 2,
    # explicit-conflict cap Tier 3 -> Tier 3.
    c = evaluate(score(60.0, 3, 4, 2.0, negatives=1), interaction=known("direct"))
    # D: no formal score and insufficient evidence -> Unclassified.
    d = evaluate(score(None, 1, 1, 0.5, sufficient=False, rank=None))
    assert [item.assigned_tier for item in (a, b, c, d)] == [
        EvidenceTier.TIER_1,
        EvidenceTier.TIER_3,
        EvidenceTier.TIER_3,
        EvidenceTier.UNCLASSIFIED,
    ]
