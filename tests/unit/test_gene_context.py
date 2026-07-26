"""MVP-1B coordinate normalization and gene-context semantics."""

from pathlib import Path
from typing import Any

import pytest

from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.application.gene_context import (
    GENE_CONTEXT_RULE_VERSION,
    build_coordinate_index,
    calculate_gene_context,
)
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models import (
    ContextCompleteness,
    CoordinatePosition,
    EvidenceStatus,
    GeneCoordinate,
    GffDocument,
    ProteinRecord,
    RelativePosition,
    SequenceRegion,
    StrandRelationship,
)


def coordinate(
    protein_id: str,
    start: int,
    end: int,
    strand: str = "+",
    seqid: str = "c",
    feature_type: str = "CDS",
    feature_id: str | None = None,
    parent_id: str | None = None,
) -> GeneCoordinate:
    return GeneCoordinate(
        seqid=seqid,
        source="test",
        feature_type=feature_type,
        start=start,
        end=end,
        strand=strand,
        feature_id=feature_id or f"cds_{protein_id}",
        parent_id=parent_id,
        parent_ids=[parent_id] if parent_id else [],
        protein_id=protein_id if feature_type == "CDS" else None,
    )


def index_for(
    records: list[GeneCoordinate],
    protein_ids: list[str],
    *,
    with_region: bool = True,
) -> Any:
    regions = {"c": SequenceRegion(seqid="c", start=1, end=1000)} if with_region else {}
    return build_coordinate_index(
        [ProteinRecord(protein_id=value, sequence="MSTK") for value in protein_ids],
        GffDocument(features=records, sequence_regions=regions),
    )


def fixture_context(fixture_dir: Path, candidate_id: str) -> Any:
    proteins = LocalFastaLoader().load(fixture_dir / "synthetic_proteome.fasta")
    document = LocalGff3Loader().load_document(fixture_dir / "synthetic_genome.gff3")
    return calculate_gene_context(
        "QUERY_001", candidate_id, build_coordinate_index(proteins, document), 5
    )


def test_one_based_closed_interval_distance_and_adjacency() -> None:
    index = index_for([coordinate("Q", 100, 200), coordinate("C", 201, 250)], ["Q", "C"])
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.distance_bp == 0
    assert evidence.edge_to_edge_distance_bp == 0
    assert evidence.overlap_bp == 0
    assert evidence.coordinate_position is CoordinatePosition.RIGHT_OF_QUERY


def test_overlap_is_inclusive_and_separate_from_distance(fixture_dir: Path) -> None:
    evidence = fixture_context(fixture_dir, "OVERLAP_001")
    assert evidence.distance_bp == 0
    assert evidence.edge_to_edge_distance_bp == 0
    assert evidence.overlap_bp == 21
    assert evidence.relative_position is RelativePosition.OVERLAPPING
    assert evidence.coordinate_position is CoordinatePosition.OVERLAPPING


def test_self_pair_has_explicit_same_feature_semantics(fixture_dir: Path) -> None:
    evidence = fixture_context(fixture_dir, "QUERY_001")
    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.relative_position is RelativePosition.SAME_FEATURE
    assert evidence.coordinate_position is CoordinatePosition.SAME_FEATURE
    assert evidence.feature_index_delta == 0


def test_different_contig_is_not_applicable_not_negative(fixture_dir: Path) -> None:
    evidence = fixture_context(fixture_dir, "CONTIG2_001")
    assert evidence.status is EvidenceStatus.NOT_APPLICABLE
    assert evidence.same_contig is False
    assert evidence.distance_bp is None
    assert evidence.edge_to_edge_distance_bp is None
    assert evidence.overlap_bp is None
    assert evidence.strand_relationship is StrandRelationship.DIFFERENT_CONTIG


def test_missing_coordinate_remains_missing(fixture_dir: Path) -> None:
    evidence = fixture_context(fixture_dir, "FRAG_001")
    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.query_contig == "contig1"
    assert evidence.query_start == 50
    assert evidence.candidate_contig is None
    assert evidence.candidate_start is None
    assert evidence.distance_bp is None
    assert "missing_coordinate:FRAG_001" in evidence.warnings
    assert evidence.edge_to_edge_distance_bp is None


@pytest.mark.parametrize(
    ("query_strand", "candidate_strand", "expected"),
    [
        ("+", "+", StrandRelationship.SAME_DIRECTION),
        ("+", "-", StrandRelationship.CONVERGENT),
        ("-", "+", StrandRelationship.DIVERGENT),
        ("-", "-", StrandRelationship.SAME_DIRECTION),
        ("?", "+", StrandRelationship.UNKNOWN),
    ],
)
def test_strand_relationships_for_non_overlapping_features(
    query_strand: str, candidate_strand: str, expected: StrandRelationship
) -> None:
    index = index_for(
        [
            coordinate("Q", 100, 200, query_strand),
            coordinate("C", 300, 400, candidate_strand),
        ],
        ["Q", "C"],
    )
    assert calculate_gene_context("Q", "C", index, 1).strand_relationship is expected


def test_opposite_strand_overlap_is_opposite_parallel() -> None:
    index = index_for(
        [coordinate("Q", 100, 250, "+"), coordinate("C", 200, 300, "-")],
        ["Q", "C"],
    )
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.strand_relationship is StrandRelationship.OPPOSITE_PARALLEL


def test_relative_position_uses_query_transcription_direction() -> None:
    index = index_for(
        [coordinate("Q", 100, 200, "-"), coordinate("C", 300, 400, "+")],
        ["Q", "C"],
    )
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.coordinate_position is CoordinatePosition.RIGHT_OF_QUERY
    assert evidence.relative_position is RelativePosition.UPSTREAM


def test_representative_units_prevent_gene_cds_double_counting() -> None:
    records = [
        coordinate("Q", 100, 200, feature_id="cds_q", parent_id="gene_q"),
        coordinate("Q", 90, 210, feature_type="gene", feature_id="gene_q"),
        coordinate("X", 220, 260, feature_id="cds_x", parent_id="gene_x"),
        coordinate("X", 215, 265, feature_type="gene", feature_id="gene_x"),
        coordinate("C", 300, 400, feature_id="cds_c", parent_id="gene_c"),
        coordinate("C", 290, 410, feature_type="gene", feature_id="gene_c"),
    ]
    index = index_for(records, ["Q", "X", "C"])
    evidence = calculate_gene_context("Q", "C", index, 3)
    assert evidence.intervening_feature_count == 1
    assert evidence.intervening_gene_count == 1
    assert evidence.feature_index_delta == 2


def test_unique_cds_coordinate_is_preferred_over_gene_extent() -> None:
    records = [
        coordinate("Q", 100, 200, feature_id="cds_q", parent_id="gene_q"),
        coordinate("Q", 80, 220, feature_type="gene", feature_id="gene_q"),
    ]
    index = index_for(records, ["Q"])
    normalized = index.by_protein["Q"]
    assert (normalized.start, normalized.end, normalized.feature_type) == (100, 200, "CDS")
    assert "gene_cds_extent_difference:cds_selected" in normalized.warnings


def test_multiple_distinct_cds_coordinates_are_ambiguous() -> None:
    records = [coordinate("Q", 100, 200), coordinate("Q", 300, 400)]
    index = index_for(records, ["Q"])
    evidence = calculate_gene_context("Q", "Q", index, 1)
    assert "Q" in index.ambiguous_proteins
    assert evidence.status is EvidenceStatus.FAILED
    assert evidence.distance_bp is None


def test_contradictory_replicons_are_fatal() -> None:
    records = [coordinate("Q", 100, 200), coordinate("Q", 100, 200, seqid="other")]
    with pytest.raises(InputValidationError, match="Contradictory replicons"):
        index_for(records, ["Q"])


def test_missing_sequence_region_gives_unknown_edges_and_completeness() -> None:
    index = index_for(
        [coordinate("Q", 100, 200), coordinate("C", 300, 400)],
        ["Q", "C"],
        with_region=False,
    )
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.query_left_edge_distance_bp is None
    assert evidence.context_completeness is ContextCompleteness.UNKNOWN
    assert "missing_sequence_region:c" in evidence.warnings


def test_contig_edge_distances_and_window_completeness(fixture_dir: Path) -> None:
    near = fixture_context(fixture_dir, "NEAR_001")
    assert near.query_left_edge_distance_bp == 49
    assert near.candidate_right_edge_distance_bp == 50000 - 420
    assert near.context_completeness is ContextCompleteness.LEFT_TRUNCATED
    assert near.within_neighborhood_window is True

    assert near.query_distance_to_contig_left_edge == 49
    assert near.within_neighborhood_gene_count is True


def test_rule_version_is_exposed_in_evidence(fixture_dir: Path) -> None:
    evidence = fixture_context(fixture_dir, "NEAR_001")
    assert evidence.calculation_rule_version == GENE_CONTEXT_RULE_VERSION
    assert evidence.provenance[0].source_version == GENE_CONTEXT_RULE_VERSION


def test_non_gene_feature_is_counted_only_in_all_feature_count() -> None:
    records = [
        coordinate("Q", 100, 200),
        GeneCoordinate(
            seqid="c",
            source="test",
            feature_type="tRNA",
            start=220,
            end=260,
            strand="+",
            feature_id="trna_1",
        ),
        coordinate("C", 300, 400),
    ]
    index = index_for(records, ["Q", "C"])
    evidence = calculate_gene_context("Q", "C", index, 3)
    assert evidence.intervening_feature_count == 1
    assert evidence.intervening_gene_count == 0


def test_same_contig_left_side_and_separated_distance() -> None:
    index = index_for(
        [coordinate("Q", 300, 400), coordinate("C", 100, 200)],
        ["Q", "C"],
    )
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.coordinate_position is CoordinatePosition.LEFT_OF_QUERY
    assert evidence.relative_position is RelativePosition.UPSTREAM
    assert evidence.edge_to_edge_distance_bp == 99


def test_neighborhood_window_includes_exact_index_boundary() -> None:
    records = [
        coordinate("Q", 100, 150),
        coordinate("X", 200, 250),
        coordinate("C", 300, 350),
    ]
    index = index_for(records, ["Q", "X", "C"])
    assert calculate_gene_context("Q", "C", index, 2).within_neighborhood_gene_count is True
    assert calculate_gene_context("Q", "C", index, 1).within_neighborhood_gene_count is False


def test_right_truncated_context_is_distinct() -> None:
    records = [
        coordinate("A", 100, 150),
        coordinate("Q", 300, 350),
        coordinate("C", 500, 550),
    ]
    index = index_for(records, ["A", "Q", "C"])
    evidence = calculate_gene_context("Q", "C", index, 1)
    assert evidence.context_completeness is ContextCompleteness.RIGHT_TRUNCATED


def test_representative_feature_order_is_deterministic() -> None:
    records = [
        coordinate("Q", 100, 150),
        coordinate("B", 200, 250),
        coordinate("A", 200, 250),
        coordinate("C", 300, 350),
    ]
    first = index_for(records, ["Q", "B", "A", "C"])
    second = index_for(list(reversed(records)), ["C", "A", "B", "Q"])
    first_ids = [item.representative_id for item in first.features_by_contig["c"]]
    second_ids = [item.representative_id for item in second.features_by_contig["c"]]
    assert first_ids == second_ids == ["protein:Q", "protein:A", "protein:B", "protein:C"]
