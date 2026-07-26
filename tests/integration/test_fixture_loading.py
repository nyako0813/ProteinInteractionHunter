"""Synthetic local-only fixture loading."""

from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.fasta import (
    LocalFastaLoader,
    duplicate_sequence_groups,
)
from protein_interaction_hunter.adapters.local.gff import (
    LocalGff3Loader,
    coordinates_by_protein,
)
from protein_interaction_hunter.application.validation import validate_local_inputs
from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import EvidenceStatus


def test_synthetic_fasta_loads_and_normalizes(fixture_dir: Path) -> None:
    records = LocalFastaLoader().load(fixture_dir / "synthetic_proteome.fasta")
    assert len(records) == 13
    assert records[0].protein_id == "QUERY_001"
    assert records[0].gene_id == "gene_query"
    assert records[0].locus_tag == "LT0001"
    assert all(record.sequence == record.sequence.upper() for record in records)
    assert duplicate_sequence_groups(records) == [["DUP_001", "QUERY_001"]]


def test_synthetic_gff_loads_identifiers_and_decodes_attributes(
    fixture_dir: Path,
) -> None:
    records = LocalGff3Loader().load(fixture_dir / "synthetic_genome.gff3")
    index = coordinates_by_protein(records)
    assert len(index) == 11
    assert index["QUERY_001"].parent_id == "gene_query"
    assert index["CONTIG2_001"].attributes["Note"] == ["URL decoded"]
    assert "FRAG_001" not in index


def test_synthetic_annotation_allows_missing_values(fixture_dir: Path) -> None:
    records = LocalAnnotationTsvLoader().load(fixture_dir / "synthetic_annotations.tsv")
    assert len(records) == 11
    hypothetical = next(record for record in records if record.protein_id == "HYP_001")
    assert hypothetical.annotation_confidence is None
    assert hypothetical.status is EvidenceStatus.AVAILABLE


def test_combined_input_summary_and_query_existence(valid_config_path: Path) -> None:
    summary = validate_local_inputs(load_config(valid_config_path))
    assert summary.protein_count == 13
    assert summary.query_count == 1
    assert summary.gff_coordinate_count == 11
    assert summary.annotation_count == 11
    assert summary.duplicate_sequence_group_count == 1
    assert summary.missing_coordinate_count == 2
    assert summary.hypothetical_protein_count == 1


def test_gff_document_retains_sequence_region(fixture_dir: Path) -> None:
    document = LocalGff3Loader().load_document(fixture_dir / "synthetic_genome.gff3")
    assert document.sequence_regions["contig1"].start == 1
    assert document.sequence_regions["contig1"].end == 50000
    assert document.warnings == []


def test_gff_reversed_coordinate_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "reversed.gff3"
    path.write_text(
        "##gff-version 3\nc\ttest\tCDS\t200\t100\t.\t+\t0\tID=x;protein_id=P1\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="Invalid GFF line"):
        LocalGff3Loader().load_document(path)


def test_conflicting_sequence_regions_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "conflicting.gff3"
    path.write_text(
        "##gff-version 3\n"
        "##sequence-region c 1 1000\n"
        "##sequence-region c 1 2000\n"
        "c\ttest\tgene\t10\t20\t.\t+\t.\tID=x\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="Conflicting ##sequence-region"):
        LocalGff3Loader().load_document(path)
