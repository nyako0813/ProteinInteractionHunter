from __future__ import annotations

from typing import Any

import pytest

from protein_interaction_hunter.application.scoring import (
    candidate_score_from_integrated,
    normalize_domain,
    normalize_functional,
    normalize_fusion,
    normalize_genome_context,
    normalize_known_interactions,
    normalize_operon,
    normalize_orthology,
    normalize_phylogenetic_profile,
    rank_scores,
    score_pair,
)
from protein_interaction_hunter.config import ScoringConfig
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    EvidenceOrigin,
    EvidenceStatus,
    IdentifierMatchStatus,
    OperonProxyStatus,
)
from protein_interaction_hunter.models.evidence import (
    DomainEvidence,
    FunctionalEvidence,
    FusionEvidence,
    GenomeContextEvidence,
    KnownInteractionEvidence,
    LocalizationEvidence,
    OperonEvidence,
    OrthologRecord,
    PhylogeneticProfileEvidence,
)
from protein_interaction_hunter.models.protein import CandidateProtein


def candidate(
    protein_id: str = "CANDIDATE",
    query_id: str = "QUERY",
    disposition: CandidateDisposition = CandidateDisposition.INCLUDED,
) -> CandidateProtein:
    return CandidateProtein(
        query_id=query_id,
        protein_id=protein_id,
        disposition=disposition,
        sequence_length=100,
        identifier_match_status=IdentifierMatchStatus.EXACT_MATCH,
    )


def config(**updates: object) -> ScoringConfig:
    data = ScoringConfig().model_dump()
    data["enabled"] = True
    data.update(updates)
    return ScoringConfig.model_validate(data)


def context(value: bool = True) -> GenomeContextEvidence:
    return GenomeContextEvidence(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.ANNOTATION,
        same_contig=True,
        within_neighborhood_window=value,
        intervening_gene_count=0 if value else 4,
        distance_bp=25 if value else 1000,
        overlap_bp=0,
        calculation_rule_version="context-v1",
    )


def operon(status: OperonProxyStatus = OperonProxyStatus.SUPPORTED) -> OperonEvidence:
    return OperonEvidence(
        status=EvidenceStatus.AVAILABLE,
        proxy_status=status,
        calculation_rule_version="operon-v1",
    )


def domain(matched: bool | None = True) -> DomainEvidence:
    return DomainEvidence(
        status=EvidenceStatus.AVAILABLE if matched is not None else EvidenceStatus.MISSING,
        protein_id="CANDIDATE",
        pair_matched=matched,
        calculation_rule_version="domain-v1",
    )


def functional(matched: bool | None = True) -> FunctionalEvidence:
    return FunctionalEvidence(
        status=EvidenceStatus.AVAILABLE if matched is not None else EvidenceStatus.MISSING,
        matched=matched,
        calculation_rule_version="functional-v1",
    )


def localization(compatibility: bool | None = True) -> LocalizationEvidence:
    return LocalizationEvidence(
        status=EvidenceStatus.AVAILABLE if compatibility is not None else EvidenceStatus.MISSING,
        protein_id="CANDIDATE",
        compatibility=compatibility,
        calculation_rule_version="localization-v1",
    )


def orthology(supported: bool = True, ambiguity: bool = False) -> OrthologRecord:
    return OrthologRecord(
        status=EvidenceStatus.AVAILABLE,
        protein_id="CANDIDATE",
        reference_id="REF",
        pair_supported=supported,
        paralog_ambiguity=ambiguity,
        calculation_rule_version="orthology-v1",
    )


def profile(
    similarity: float = 0.8,
    shared_presence: int = 4,
    shared_absence: int = 0,
    discordant: int = 1,
) -> PhylogeneticProfileEvidence:
    return PhylogeneticProfileEvidence(
        status=EvidenceStatus.AVAILABLE,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        informative_species_count=shared_presence + shared_absence + discordant,
        shared_presence_count=shared_presence,
        shared_absence_count=shared_absence,
        discordant_count=discordant,
        profile_similarity=similarity,
        pair_supported=similarity >= 0.8,
        calculation_rule_version="profile-v1",
    )


def fusion(supported: bool | None = True, conflict: bool = False) -> FusionEvidence:
    return FusionEvidence(
        status=EvidenceStatus.MISSING if supported is None else EvidenceStatus.AVAILABLE,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        supporting_record_count=1,
        qualifying_record_count=1 if supported else 0,
        pair_supported=supported,
        conflicting_terms=["excessive_component_overlap"] if conflict else [],
        calculation_rule_version="fusion-v1",
    )


def known(kind: str = "direct") -> KnownInteractionEvidence:
    return KnownInteractionEvidence(
        status=EvidenceStatus.AVAILABLE,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        supporting_record_count=1 if kind != "none" else 0,
        qualifying_record_count=1 if kind not in {"none", "predicted"} else 0,
        direct_record_count=1 if kind == "direct" else 0,
        physical_record_count=1 if kind in {"direct", "physical"} else 0,
        biological_record_count=1 if kind == "functional" else 0,
        independent_publication_count=1 if kind != "none" else 0,
        independent_source_count=1 if kind != "none" else 0,
        interaction_types=["predicted"] if kind == "predicted" else [],
        pair_supported=kind in {"direct", "physical"},
        direct_interaction_supported=kind == "direct",
        physical_interaction_supported=kind == "physical",
        functional_association_supported=kind == "functional",
        calculation_rule_version="known-v1",
    )


def test_component_normalization_table() -> None:
    cfg = config()
    assert normalize_genome_context(context(), cfg).normalized_value == 0.8
    assert normalize_operon(operon(), cfg).normalized_value == 1.0
    assert normalize_operon(operon(OperonProxyStatus.PARTIAL_SUPPORT), cfg).normalized_value == 0.3
    assert normalize_domain([domain(True)], cfg).normalized_value == 1.0
    assert normalize_domain([domain(False)], cfg).normalized_value == 0.0
    assert normalize_functional([functional(True)], cfg).normalized_value == 1.0
    assert normalize_functional([functional(False)], cfg).normalized_value == 0.0
    assert normalize_orthology([orthology(False)], cfg)[0].normalized_value == 0.0
    assert normalize_phylogenetic_profile(profile(), cfg)[0].normalized_value == 0.8
    assert normalize_fusion(fusion(False), cfg)[0].normalized_value == 0.0
    assert normalize_known_interactions(known("direct"), cfg).normalized_value == 1.0
    assert normalize_known_interactions(known("physical"), cfg).normalized_value == 0.9
    assert normalize_known_interactions(known("functional"), cfg).normalized_value == 0.5
    assert normalize_known_interactions(known("predicted"), cfg).normalized_value == 0.2
    assert normalize_known_interactions(known("none"), cfg).normalized_value == 0.0


def test_missing_not_run_and_weight_zero_are_excluded() -> None:
    cfg = config(weights={**ScoringConfig().weights.model_dump(), "fusion": 0.0})
    missing = normalize_domain([domain(None)], cfg)
    not_run = normalize_operon(None, cfg)
    zero = normalize_fusion(fusion(True), cfg)[0]
    assert (missing.applied, missing.exclusion_reason) == (False, "evidence_missing")
    assert (not_run.applied, not_run.exclusion_reason) == (False, "evidence_not_run")
    assert (zero.applied, zero.exclusion_reason) == (False, "weight_zero")


def test_explicit_negative_only_and_penalties_are_audited() -> None:
    cfg = config()
    score = score_pair(
        "QUERY",
        candidate(),
        cfg,
        domains=[domain(True)],
        localization=localization(False),
        fusion=fusion(False, conflict=True),
    )
    components = {item.component_name: item for item in score.component_scores}
    assert components["localization"].normalized_value == -0.25
    assert components["fusion"].normalized_value == -0.25
    assert components["domain_pair"].normalized_value == 1.0
    assert score.negative_component_count == 2
    assert [item.component_name for item in score.applied_penalties] == [
        "localization",
        "fusion",
    ]


def test_paralog_ambiguity_reduces_positive_strength() -> None:
    component, penalties = normalize_orthology([orthology(True, True)], config())
    assert component.normalized_value == pytest.approx(0.9)
    assert penalties[0].penalty_name == "ambiguous_mapping"


def test_profile_shared_absence_does_not_create_full_support() -> None:
    component, _ = normalize_phylogenetic_profile(
        profile(similarity=1.0, shared_presence=1, shared_absence=9, discordant=0),
        config(),
    )
    assert component.normalized_value == pytest.approx(0.1)


def test_strong_profile_discordance_is_explicit_negative() -> None:
    component, penalties = normalize_phylogenetic_profile(
        profile(similarity=0.0, shared_presence=0, shared_absence=0, discordant=5),
        config(),
    )
    assert component.normalized_value == -0.25
    assert penalties[0].explanation == "strong_profile_discordance"


def test_manual_category_cap_and_score_calculation() -> None:
    # Genomic: (0.8 * 0.75) + (1.0 * 0.75) = 1.35, cap weight 1.5.
    # Functional: (1.0 * 0.75) + (0.0 * 0.75) = 0.75, cap weight 1.5.
    # Cellular: 1.0 * 0.5 = 0.5. Evolutionary: 1.0 * 0.75 = 0.75.
    # Direct: two 1.5 weights scale to 1.0 each = 2.0, cap weight 2.0.
    # raw=5.35, available=6.25, normalized=0.856, output=85.6.
    score = score_pair(
        "QUERY",
        candidate(),
        config(),
        genome_context=context(),
        operon=operon(),
        domains=[domain(True)],
        functional=[functional(False)],
        localization=localization(True),
        orthology=[orthology(True)],
        fusion=fusion(True),
        known_interactions=known("direct"),
    )
    assert score.raw_weighted_sum == pytest.approx(5.35, abs=1e-12)
    assert score.available_weight == pytest.approx(6.25, abs=1e-12)
    assert score.normalized_score == pytest.approx(0.856, abs=1e-12)
    assert score.provisional_score == pytest.approx(85.6, abs=1e-12)
    assert score.output_score == pytest.approx(85.6, abs=1e-12)
    assert score.sufficient_evidence is True
    weights = {item.component_name: item.effective_weight for item in score.component_scores}
    assert weights["genome_context"] == pytest.approx(0.75)
    assert weights["operon_proxy"] == pytest.approx(0.75)
    assert weights["fusion"] == pytest.approx(1.0)
    assert weights["known_interactions"] == pytest.approx(1.0)


def test_no_available_evidence_and_insufficient_evidence() -> None:
    empty = score_pair("QUERY", candidate(), config())
    assert empty.status is EvidenceStatus.MISSING
    assert empty.normalized_score is None
    assert empty.output_score is None
    one = score_pair("QUERY", candidate(), config(), domains=[domain(True)])
    assert one.provisional_score == 100.0
    assert one.output_score is None
    assert one.sufficient_evidence is False


def test_output_scale_and_weight_change_are_config_driven() -> None:
    cfg = config(
        output_scale=10.0,
        minimum_evidence_categories=1,
        weights={**ScoringConfig().weights.model_dump(), "domain_pair": 2.0},
    )
    score = score_pair("QUERY", candidate(), cfg, domains=[domain(True)])
    assert score.output_score == 10.0
    assert score.available_weight == 1.5


def test_candidate_score_mirrors_formal_normalized_score() -> None:
    integrated = score_pair(
        "QUERY", candidate(), config(minimum_evidence_categories=1), domains=[domain(True)]
    )
    score = candidate_score_from_integrated(integrated)
    assert score.total_ranking_score == 1.0
    assert score.functional_association_score == 1.0
    assert score.evidence_completeness == pytest.approx(1.0 / 7.5)


def _rankable(
    protein_id: str,
    query_id: str,
    value: float,
    disposition: CandidateDisposition = CandidateDisposition.INCLUDED,
) -> Any:
    cfg = config(minimum_evidence_categories=1)
    record = known("direct" if value == 1.0 else "functional")
    score = score_pair(
        query_id, candidate(protein_id, query_id, disposition), cfg, known_interactions=record
    )
    if value not in {1.0, 0.5}:
        score = score.model_copy(
            update={
                "normalized_score": value,
                "provisional_score": value * 100,
                "output_score": value * 100,
            }
        )
    return score


def test_dense_ranking_query_partition_tie_and_exclusions() -> None:
    candidates = [
        candidate("B", "Q1"),
        candidate("A", "Q1"),
        candidate("C", "Q1"),
        candidate("SELF", "Q1", CandidateDisposition.EXCLUDED),
        candidate("D", "Q2"),
    ]
    scores = {
        ("Q1", "B"): _rankable("B", "Q1", 0.8),
        ("Q1", "A"): _rankable("A", "Q1", 0.8),
        ("Q1", "C"): _rankable("C", "Q1", 0.5),
        ("Q1", "SELF"): _rankable("SELF", "Q1", 1.0, CandidateDisposition.EXCLUDED),
        ("Q2", "D"): _rankable("D", "Q2", 0.5),
    }
    ranked = rank_scores(scores, candidates, 8)
    assert ranked[("Q1", "A")].rank == 1
    assert ranked[("Q1", "B")].rank == 1
    assert ranked[("Q1", "C")].rank == 2
    assert ranked[("Q1", "SELF")].rank is None
    assert ranked[("Q2", "D")].rank == 1


def test_insufficient_score_is_not_ranked() -> None:
    record = score_pair("Q", candidate("X", "Q"), config(), domains=[domain(True)])
    ranked = rank_scores({("Q", "X"): record}, [candidate("X", "Q")], 8)
    assert ranked[("Q", "X")].rank is None


def test_component_order_is_input_order_independent() -> None:
    forward = score_pair("QUERY", candidate(), config(), domains=[domain(True), domain(False)])
    reverse = score_pair("QUERY", candidate(), config(), domains=[domain(False), domain(True)])
    assert forward == reverse
    assert [item.component_name for item in forward.component_scores] == [
        "genome_context",
        "operon_proxy",
        "domain_pair",
        "functional_complementarity",
        "localization",
        "orthology",
        "phylogenetic_profile",
        "fusion",
        "known_interactions",
    ]


def test_malformed_numeric_value_is_excluded_defensively() -> None:
    malformed = PhylogeneticProfileEvidence.model_construct(
        status=EvidenceStatus.AVAILABLE,
        query_protein_id="QUERY",
        candidate_protein_id="CANDIDATE",
        informative_species_count=5,
        shared_presence_count=5,
        shared_absence_count=0,
        discordant_count=0,
        profile_similarity=float("nan"),
        conflicting_terms=[],
        warnings=[],
    )
    component, _ = normalize_phylogenetic_profile(malformed, config())
    assert component.applied is False
    assert component.normalized_value is None
    assert component.exclusion_reason == "invalid_or_unavailable_value"


def test_candidate_identifier_ambiguity_reduces_all_positive_components() -> None:
    ambiguous = candidate().model_copy(
        update={"identifier_match_status": IdentifierMatchStatus.AMBIGUOUS_MATCH}
    )
    score = score_pair(
        "QUERY",
        ambiguous,
        config(minimum_evidence_categories=1),
        domains=[domain(True)],
    )
    domain_component = next(
        item for item in score.component_scores if item.component_name == "domain_pair"
    )
    assert domain_component.normalized_value == pytest.approx(0.9)
    assert score.applied_penalties[0].explanation == "candidate_identifier_mapping_ambiguity"
