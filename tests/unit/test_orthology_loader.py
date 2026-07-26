from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.orthology import (
    LocalOrthologyTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)


def test_load_valid_orthology_table(fixture_dir: Path) -> None:
    records = LocalOrthologyTsvLoader().load(fixture_dir / "synthetic_orthology.tsv")

    assert len(records) == 4
    assert records[0].protein_id == "QUERY_001"
    assert records[0].ortholog_id == "REF_QUERY_001"
    assert records[0].identity == 0.84
    assert records[0].status is EvidenceStatus.AVAILABLE
    assert records[0].origin is EvidenceOrigin.ORTHOLOG_TRANSFERRED
    assert records[0].paralog_ambiguity is False
    assert records[3].paralog_ambiguity is True


def test_missing_orthology_table_is_rejected(
    tmp_path: Path,
) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalOrthologyTsvLoader().load(tmp_path / "missing.tsv")


def test_missing_required_column_is_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "missing-column.tsv"
    path.write_text(
        ("protein_id\treference_id\tortholog_id\nP1\tREF\tORTHO1\n"),
        encoding="utf-8",
    )

    with pytest.raises(
        InputValidationError,
        match="missing columns",
    ):
        LocalOrthologyTsvLoader().load(path)


def test_invalid_fraction_is_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "invalid-fraction.tsv"
    path.write_text(
        (
            "protein_id\treference_id\tortholog_id\t"
            "reference_organism\tidentity\tquery_coverage\t"
            "subject_coverage\tevalue\torthogroup\t"
            "relationship\tparalog_ambiguity\tsource\t"
            "source_record_id\n"
            "P1\tREF\tORTHO1\tOrganism\t1.2\t0.8\t"
            "0.9\t1e-10\tOG1\treciprocal_best_hit\t"
            "false\ttest\trow-1\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        InputValidationError,
        match="Invalid orthology annotation",
    ):
        LocalOrthologyTsvLoader().load(path)


def test_invalid_boolean_is_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "invalid-bool.tsv"
    path.write_text(
        (
            "protein_id\treference_id\tortholog_id\t"
            "reference_organism\tidentity\tquery_coverage\t"
            "subject_coverage\tevalue\torthogroup\t"
            "relationship\tparalog_ambiguity\tsource\t"
            "source_record_id\n"
            "P1\tREF\tORTHO1\tOrganism\t0.8\t0.8\t"
            "0.9\t1e-10\tOG1\treciprocal_best_hit\t"
            "maybe\ttest\trow-1\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        InputValidationError,
        match="Invalid orthology annotation",
    ):
        LocalOrthologyTsvLoader().load(path)


def test_duplicate_orthology_annotation_is_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "duplicate.tsv"
    header = (
        "protein_id\treference_id\tortholog_id\t"
        "reference_organism\tidentity\tquery_coverage\t"
        "subject_coverage\tevalue\torthogroup\t"
        "relationship\tparalog_ambiguity\tsource\t"
        "source_record_id\n"
    )
    row = (
        "P1\tREF\tORTHO1\tOrganism\t0.8\t0.8\t"
        "0.9\t1e-10\tOG1\treciprocal_best_hit\t"
        "false\ttest\trow-1\n"
    )
    path.write_text(
        header + row + row,
        encoding="utf-8",
    )

    with pytest.raises(
        InputValidationError,
        match="Duplicate orthology annotation",
    ):
        LocalOrthologyTsvLoader().load(path)
