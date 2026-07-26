"""MVP-1A candidate generation and policy behavior."""

from pathlib import Path
from typing import Any

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.application.candidates import (
    build_duplicate_groups,
    generate_candidates,
)
from protein_interaction_hunter.config import CandidateGenerationConfig
from protein_interaction_hunter.models import (
    AnnotationRecord,
    CandidateDisposition,
    EvidenceStatus,
    GeneCoordinate,
    IdentifierMatchStatus,
    ProteinRecord,
)


def generate(fixture_dir: Path, query_ids: list[str] | None = None, **updates: object) -> Any:
    policy = CandidateGenerationConfig().model_copy(update=updates)
    return generate_candidates(
        proteins=LocalFastaLoader().load(fixture_dir / "synthetic_proteome.fasta"),
        coordinates=LocalGff3Loader().load(fixture_dir / "synthetic_genome.gff3"),
        annotations=LocalAnnotationTsvLoader().load(fixture_dir / "synthetic_annotations.tsv"),
        query_ids=query_ids or ["QUERY_001"],
        policy=policy,
    )


def by_id(result: Any, protein_id: str) -> Any:
    return next(candidate for candidate in result.candidates if candidate.protein_id == protein_id)


def test_self_exclusion_and_duplicate_at_other_locus(fixture_dir: Path) -> None:
    result = generate(fixture_dir)
    query = by_id(result, "QUERY_001")
    duplicate = by_id(result, "DUP_001")
    assert query.disposition is CandidateDisposition.EXCLUDED
    assert "self_candidate" in query.disposition_reasons
    assert duplicate.disposition is CandidateDisposition.FLAGGED
    assert query.old_locus_tag == "OLD0001"
    assert query.original_identifiers["gff_id"] == ["cds_query"]
    assert query.normalized_identifiers["old_locus_tag"] == ["old0001"]
    assert duplicate.duplicate_sequence_group == query.duplicate_sequence_group


def test_duplicate_group_id_and_candidate_order_are_deterministic(fixture_dir: Path) -> None:
    first = generate(fixture_dir)
    second = generate(fixture_dir)
    assert first.duplicate_groups == second.duplicate_groups
    assert [item.protein_id for item in first.candidates] == sorted(
        item.protein_id for item in first.candidates
    )
    group_id = next(iter(first.duplicate_groups))
    assert group_id.startswith("dup-")
    assert first.duplicate_groups[group_id] == ["DUP_001", "QUERY_001"]


def test_multiple_queries_generate_deterministic_pair_sets(fixture_dir: Path) -> None:
    result = generate(fixture_dir, query_ids=["QUERY_001", "NEAR_001"])
    assert len(result.candidates) == 26
    self_pairs = [item for item in result.candidates if item.query_id == item.protein_id]
    assert all(item.disposition is CandidateDisposition.EXCLUDED for item in self_pairs)


def test_duplicate_policy_flag_and_exclude(fixture_dir: Path) -> None:
    assert by_id(generate(fixture_dir), "DUP_001").disposition is CandidateDisposition.FLAGGED
    excluded = by_id(generate(fixture_dir, duplicate_sequence_policy="exclude"), "DUP_001")
    assert excluded.disposition is CandidateDisposition.EXCLUDED
    assert "duplicate_sequence_policy_exclude" in excluded.disposition_reasons


def test_fragment_minimum_description_flag_and_exclude(fixture_dir: Path) -> None:
    flagged = by_id(generate(fixture_dir, minimum_length_aa=30), "FRAG_001")
    assert flagged.is_fragment_candidate
    assert "length_below_minimum:30" in flagged.fragment_reasons
    assert "description_keyword:fragment" in flagged.fragment_reasons
    assert flagged.disposition is CandidateDisposition.FLAGGED
    excluded = by_id(generate(fixture_dir, fragment_policy="exclude"), "FRAG_001")
    assert excluded.disposition is CandidateDisposition.EXCLUDED


def test_hypothetical_and_missing_data_are_retained(fixture_dir: Path) -> None:
    result = generate(fixture_dir)
    hypothetical = by_id(result, "HYP_001")
    missing_coordinate = by_id(result, "FRAG_001")
    missing_annotation = by_id(result, "FAR_001")
    assert hypothetical.is_hypothetical
    assert hypothetical.disposition is not CandidateDisposition.EXCLUDED
    assert not missing_coordinate.has_coordinate
    assert "missing_coordinate" in missing_coordinate.warnings
    assert missing_coordinate.coordinate_status is EvidenceStatus.MISSING
    assert not missing_annotation.has_annotation
    assert "missing_annotation" in missing_annotation.warnings
    assert missing_annotation.annotation_status is EvidenceStatus.MISSING


def test_duplicate_group_does_not_merge_distinct_proteins() -> None:
    records = [
        ProteinRecord(protein_id="B", sequence="MSTK"),
        ProteinRecord(protein_id="A", sequence="MSTK"),
    ]
    groups = build_duplicate_groups(records)
    assert list(groups.values()) == [["A", "B"]]


def test_ambiguous_annotation_alias_does_not_override_exact_gff_mapping() -> None:
    proteins = [
        ProteinRecord(protein_id="Q", sequence="MSTKAA"),
        ProteinRecord(protein_id="C", sequence="MSTKCC"),
    ]
    coordinates = [
        GeneCoordinate(
            seqid="c",
            feature_type="gene",
            start=10,
            end=30,
            strand="+",
            feature_id="gene-q",
            locus_tag="LQ",
        ),
        GeneCoordinate(
            seqid="c",
            feature_type="CDS",
            start=10,
            end=30,
            strand="+",
            feature_id="cds-q",
            parent_id="gene-q",
            parent_ids=["gene-q"],
            protein_id="Q",
            locus_tag="LQ",
        ),
        GeneCoordinate(
            seqid="c",
            feature_type="gene",
            start=100,
            end=120,
            strand="+",
            feature_id="gene-c",
            locus_tag="LC",
        ),
        GeneCoordinate(
            seqid="c",
            feature_type="CDS",
            start=100,
            end=120,
            strand="+",
            feature_id="cds-c",
            parent_id="gene-c",
            parent_ids=["gene-c"],
            protein_id="C",
            locus_tag="LC",
        ),
    ]
    annotations = [
        AnnotationRecord(
            protein_id="Q",
            gene_name="gene-c",
            locus_tag="LQ",
            product="query",
            status=EvidenceStatus.AVAILABLE,
        ),
        AnnotationRecord(
            protein_id="C",
            gene_name="gene-q",
            locus_tag="LC",
            product="candidate",
            status=EvidenceStatus.AVAILABLE,
        ),
    ]

    result = generate_candidates(
        proteins=proteins,
        coordinates=coordinates,
        annotations=annotations,
        query_ids=["Q"],
        policy=CandidateGenerationConfig(),
    )

    assert result.ambiguous_mapping_count == 0
    assert all(
        candidate.identifier_match_status is IdentifierMatchStatus.EXACT_MATCH
        for candidate in result.candidates
    )
    assert all(
        "ambiguous_identifier_mapping" not in candidate.warnings for candidate in result.candidates
    )
