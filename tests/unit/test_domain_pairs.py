from pathlib import Path

from protein_interaction_hunter.adapters.local.domain_rules import (
    LocalDomainRulesLoader,
)
from protein_interaction_hunter.adapters.local.domains import (
    LocalDomainTsvLoader,
)
from protein_interaction_hunter.application.domain_pairs import (
    build_domain_index,
    evaluate_domain_pairs,
)
from protein_interaction_hunter.models.enums import EvidenceStatus


def load_inputs(fixture_dir: Path):
    domain_records = LocalDomainTsvLoader().load(
        fixture_dir / "synthetic_domains.tsv"
    )
    rules_path = (
        fixture_dir
        / "rules"
        / "domain_pairs.v1.yaml"
    )
    rules = LocalDomainRulesLoader().load(rules_path)

    return build_domain_index(domain_records), rules, rules_path


def test_domain_index_is_deterministic(
    fixture_dir: Path,
) -> None:
    index, _, _ = load_inputs(fixture_dir)

    assert index["QUERY_001"][0].accession == "PF00001"
    assert index["NEAR_001"][0].architecture_index == 0


def test_catalytic_accessory_pair_matches(
    fixture_dir: Path,
) -> None:
    index, rules, path = load_inputs(fixture_dir)

    evidence = evaluate_domain_pairs(
        "QUERY_001",
        "NEAR_001",
        index,
        rules,
        path,
    )

    assert len(evidence) == 1
    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].pair_matched is True
    assert evidence[0].pair_rule_id == (
        "catalytic-accessory-domain-v1"
    )
    assert evidence[0].paired_accession == "PF00001"
    assert evidence[0].accession == "PF00002"
    assert evidence[0].quality is None


def test_catalytic_membrane_pair_matches(
    fixture_dir: Path,
) -> None:
    index, rules, path = load_inputs(fixture_dir)

    evidence = evaluate_domain_pairs(
        "QUERY_001",
        "MEM_001",
        index,
        rules,
        path,
    )

    assert evidence[0].pair_matched is True
    assert evidence[0].pair_rule_id == (
        "catalytic-membrane-domain-v1"
    )


def test_unmatched_pair_remains_available(
    fixture_dir: Path,
) -> None:
    index, rules, path = load_inputs(fixture_dir)

    evidence = evaluate_domain_pairs(
        "QUERY_001",
        "HOUSE_001",
        index,
        rules,
        path,
    )

    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].pair_matched is False
    assert evidence[0].pair_rule_id is None
    assert evidence[0].conflicting_terms == [
        "no_matching_domain_pair_rule"
    ]


def test_missing_candidate_domain_is_missing(
    fixture_dir: Path,
) -> None:
    index, rules, path = load_inputs(fixture_dir)

    evidence = evaluate_domain_pairs(
        "QUERY_001",
        "UNKNOWN_001",
        index,
        rules,
        path,
    )

    assert evidence[0].status is EvidenceStatus.MISSING
    assert evidence[0].pair_matched is None
    assert "candidate_domain_annotation" in (
        evidence[0].conflicting_terms
    )


def test_missing_query_domain_is_missing(
    fixture_dir: Path,
) -> None:
    index, rules, path = load_inputs(fixture_dir)

    evidence = evaluate_domain_pairs(
        "UNKNOWN_QUERY",
        "NEAR_001",
        index,
        rules,
        path,
    )

    assert evidence[0].status is EvidenceStatus.MISSING
    assert "query_domain_annotation" in (
        evidence[0].conflicting_terms
    )