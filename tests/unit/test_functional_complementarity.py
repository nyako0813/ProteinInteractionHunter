from pathlib import Path
from typing import Any

from protein_interaction_hunter.adapters.local.functional_rules import (
    LocalFunctionalRulesLoader,
)
from protein_interaction_hunter.application.functional_complementarity import (
    evaluate_functional_complementarity,
    normalize_annotation_text,
)
from protein_interaction_hunter.models.enums import (
    EvidenceStatus,
    PredictedRelationshipType,
)


def load_rules(fixture_dir: Path) -> Any:
    path = (
        fixture_dir
        / "rules"
        / "functional_complementarity.v1.yaml"
    )
    return LocalFunctionalRulesLoader().load(path), path


def test_annotation_normalization_is_deterministic() -> None:
    assert (
        normalize_annotation_text("  Query—ENZYME  ")
        == "query enzyme"
    )


def test_enzyme_accessory_rule_matches(fixture_dir: Path) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "Query enzyme",
        "nearby accessory candidate",
        rules,
        path,
    )

    assert len(evidence) == 1
    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].matched is True
    assert evidence[0].query_role == "enzyme"
    assert evidence[0].candidate_role == "accessory_factor"
    assert (
        evidence[0].relationship_hint
        is PredictedRelationshipType.ACCESSORY_FACTOR
    )
    assert evidence[0].quality is None


def test_membrane_pair_rule_matches(fixture_dir: Path) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "Query enzyme",
        "membrane protein candidate",
        rules,
        path,
    )

    assert evidence[0].matched is True
    assert evidence[0].candidate_role == "membrane_protein"
    assert (
        evidence[0].relationship_hint
        is PredictedRelationshipType.PATHWAY_ASSOCIATION
    )


def test_unmatched_annotations_are_available_not_negative(
    fixture_dir: Path,
) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "Query enzyme",
        "housekeeping protein",
        rules,
        path,
    )

    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].matched is False
    assert evidence[0].relationship_hint is None
    assert evidence[0].conflicting_terms == [
        "no_matching_pair_rule"
    ]


def test_missing_candidate_annotation_remains_missing(
    fixture_dir: Path,
) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "Query enzyme",
        None,
        rules,
        path,
    )

    assert evidence[0].status is EvidenceStatus.MISSING
    assert evidence[0].matched is None
    assert "candidate_annotation" in evidence[0].conflicting_terms


def test_excluded_role_term_prevents_match(
    fixture_dir: Path,
) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "inactive enzyme",
        "accessory factor",
        rules,
        path,
    )

    assert evidence[0].matched is False


def test_term_matching_does_not_match_inside_larger_word(
    fixture_dir: Path,
) -> None:
    rules, path = load_rules(fixture_dir)

    evidence = evaluate_functional_complementarity(
        "enzyme",
        "membranesome protein",
        rules,
        path,
    )

    assert evidence[0].matched is False