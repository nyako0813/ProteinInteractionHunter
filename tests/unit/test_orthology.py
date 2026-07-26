from pathlib import Path

from protein_interaction_hunter.adapters.local.orthology import (
    LocalOrthologyTsvLoader,
)
from protein_interaction_hunter.application.orthology import (
    ORTHOLOGY_ENGINE_VERSION,
    build_orthology_index,
    evaluate_orthology_pair,
)
from protein_interaction_hunter.models.enums import EvidenceStatus
from protein_interaction_hunter.models.evidence import OrthologRecord


def load_index(
    fixture_dir: Path,
) -> dict[str, list[OrthologRecord]]:
    records = LocalOrthologyTsvLoader().load(fixture_dir / "synthetic_orthology.tsv")
    return build_orthology_index(records)


def test_orthology_index_is_deterministic(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    assert index["QUERY_001"][0].reference_id == "ref_archaea"
    assert index["NEAR_001"][0].ortholog_id == "REF_NEAR_001"


def test_shared_orthogroup_supports_pair(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    evidence = evaluate_orthology_pair(
        "QUERY_001",
        "NEAR_001",
        index,
    )

    assert len(evidence) == 1
    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].pair_supported is True
    assert evidence[0].shared_orthogroup is True
    assert evidence[0].paired_orthogroup == "OG0001"
    assert evidence[0].support_terms == ["shared_orthogroup"]
    assert evidence[0].calculation_rule_version == ORTHOLOGY_ENGINE_VERSION
    assert evidence[0].quality is None


def test_different_orthogroup_is_not_supported(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    evidence = evaluate_orthology_pair(
        "QUERY_001",
        "MEM_001",
        index,
    )

    assert evidence[0].status is EvidenceStatus.AVAILABLE
    assert evidence[0].pair_supported is False
    assert evidence[0].shared_orthogroup is False
    assert evidence[0].conflicting_terms == ["no_shared_orthology_support"]


def test_paralog_ambiguity_is_preserved(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    evidence = evaluate_orthology_pair(
        "QUERY_001",
        "PARA_001",
        index,
    )

    assert evidence[0].pair_supported is True
    assert evidence[0].shared_orthogroup is True
    assert evidence[0].paralog_ambiguity is True
    assert "paralog_ambiguity" in evidence[0].conflicting_terms


def test_missing_candidate_orthology_is_missing(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    evidence = evaluate_orthology_pair(
        "QUERY_001",
        "UNKNOWN_001",
        index,
    )

    assert evidence[0].status is EvidenceStatus.MISSING
    assert evidence[0].pair_supported is None
    assert evidence[0].conflicting_terms == ["candidate_orthology_annotation"]


def test_missing_query_orthology_is_missing(
    fixture_dir: Path,
) -> None:
    index = load_index(fixture_dir)

    evidence = evaluate_orthology_pair(
        "UNKNOWN_QUERY",
        "NEAR_001",
        index,
    )

    assert evidence[0].status is EvidenceStatus.MISSING
    assert evidence[0].pair_supported is None
    assert evidence[0].conflicting_terms == ["query_orthology_annotation"]


def test_different_reference_is_not_supported() -> None:
    from protein_interaction_hunter.models.evidence import (
        OrthologRecord,
    )

    index = build_orthology_index(
        [
            OrthologRecord(
                status=EvidenceStatus.AVAILABLE,
                protein_id="QUERY_001",
                reference_id="REF_A",
                ortholog_id="Q1",
                orthogroup="OG1",
            ),
            OrthologRecord(
                status=EvidenceStatus.AVAILABLE,
                protein_id="CANDIDATE_001",
                reference_id="REF_B",
                ortholog_id="C1",
                orthogroup="OG1",
            ),
        ]
    )

    evidence = evaluate_orthology_pair(
        "QUERY_001",
        "CANDIDATE_001",
        index,
    )

    assert evidence[0].pair_supported is False
    assert evidence[0].shared_orthogroup is False
    assert evidence[0].conflicting_terms == ["no_shared_reference"]
