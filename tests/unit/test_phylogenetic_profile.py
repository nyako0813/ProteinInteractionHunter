from pathlib import Path

from protein_interaction_hunter.adapters.local.phylogenetic_profile import (
    LocalPhylogeneticProfileTsvLoader,
)
from protein_interaction_hunter.application.phylogenetic_profile import (
    PHYLOGENETIC_PROFILE_ENGINE_VERSION,
    build_phylogenetic_profile_index,
    evaluate_phylogenetic_profile_pair,
)
from protein_interaction_hunter.models.enums import EvidenceStatus
from protein_interaction_hunter.models.phylogenetic_profile import (
    PhylogeneticProfileObservation,
)


def evaluate(
    query: str,
    candidate: str,
    records: list[PhylogeneticProfileObservation],
):
    return evaluate_phylogenetic_profile_pair(
        query,
        candidate,
        build_phylogenetic_profile_index(records),
        minimum_shared_species=2,
        minimum_informative_species=3,
        minimum_profile_similarity=0.8,
    )


def fixture_records(fixture_dir: Path) -> list[PhylogeneticProfileObservation]:
    return LocalPhylogeneticProfileTsvLoader().load(
        fixture_dir / "synthetic_phylogenetic_profiles.tsv"
    )


def test_shared_profile_supports_pair(fixture_dir: Path) -> None:
    evidence = evaluate("QUERY_001", "NEAR_001", fixture_records(fixture_dir))
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.informative_species_count == 4
    assert evidence.shared_presence_count == 2
    assert evidence.shared_absence_count == 2
    assert evidence.discordant_count == 0
    assert evidence.unknown_count == 1
    assert evidence.profile_similarity == 1.0
    assert evidence.pair_supported is True
    assert evidence.calculation_rule_version == PHYLOGENETIC_PROFILE_ENGINE_VERSION
    assert evidence.quality is None


def test_discordant_profile_is_not_supported(fixture_dir: Path) -> None:
    evidence = evaluate("QUERY_001", "MEM_001", fixture_records(fixture_dir))
    assert evidence.informative_species_count == 4
    assert evidence.shared_presence_count == 1
    assert evidence.shared_absence_count == 1
    assert evidence.discordant_count == 2
    assert evidence.unknown_count == 1
    assert evidence.profile_similarity == 0.5
    assert evidence.pair_supported is False
    assert "profile_similarity_below_threshold" in evidence.conflicting_terms


def test_shared_absence_alone_does_not_support_pair() -> None:
    records = [
        PhylogeneticProfileObservation(protein_id=protein, species_id=species, presence=False)
        for protein in ("Q", "C")
        for species in ("S1", "S2", "S3")
    ]
    evidence = evaluate("Q", "C", records)
    assert evidence.profile_similarity == 1.0
    assert evidence.shared_absence_count == 3
    assert evidence.shared_presence_count == 0
    assert evidence.pair_supported is False
    assert "insufficient_shared_presence" in evidence.conflicting_terms


def test_missing_query_and_candidate_are_explicit(fixture_dir: Path) -> None:
    records = fixture_records(fixture_dir)
    query_missing = evaluate("UNKNOWN", "NEAR_001", records)
    candidate_missing = evaluate("QUERY_001", "UNKNOWN", records)
    assert query_missing.status is EvidenceStatus.MISSING
    assert query_missing.pair_supported is None
    assert query_missing.conflicting_terms == ["query_phylogenetic_profile"]
    assert candidate_missing.status is EvidenceStatus.MISSING
    assert candidate_missing.pair_supported is None
    assert candidate_missing.conflicting_terms == ["candidate_phylogenetic_profile"]


def test_insufficient_informative_species_is_not_supported() -> None:
    records = [
        PhylogeneticProfileObservation(protein_id=protein, species_id=species, presence=True)
        for protein in ("Q", "C")
        for species in ("S1", "S2")
    ]
    evidence = evaluate("Q", "C", records)
    assert evidence.profile_similarity == 1.0
    assert evidence.pair_supported is False
    assert "insufficient_informative_species" in evidence.conflicting_terms


def test_evaluation_is_independent_of_input_row_order(fixture_dir: Path) -> None:
    records = fixture_records(fixture_dir)
    forward = evaluate("QUERY_001", "NEAR_001", records)
    reverse = evaluate("QUERY_001", "NEAR_001", list(reversed(records)))
    assert forward == reverse
