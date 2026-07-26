import csv
from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.known_interactions import (
    LocalKnownInteractionsTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError

HEADER = [
    "protein_a_id",
    "protein_b_id",
    "interaction_type",
    "reference_organism",
    "detection_method",
    "publication_id",
    "confidence",
    "is_direct",
    "is_physical",
    "is_biological",
    "database_version",
    "protein_a_reference_id",
    "protein_b_reference_id",
    "source",
    "source_record_id",
    "notes",
    "identifier_mapping_status",
]


def _row(**updates: str) -> dict[str, str]:
    row = {
        "protein_a_id": "QUERY",
        "protein_b_id": "CANDIDATE",
        "interaction_type": "direct",
        "reference_organism": "Test organism",
        "detection_method": "Y2H",
        "publication_id": "PMID:1",
        "confidence": "0.8",
        "is_direct": "true",
        "is_physical": "yes",
        "is_biological": "1",
        "database_version": "v1",
        "protein_a_reference_id": "QA",
        "protein_b_reference_id": "CB",
        "source": "TestDB",
        "source_record_id": "R1",
        "notes": "fixture",
        "identifier_mapping_status": "mapped",
    }
    row.update(updates)
    return row


def _write(path: Path, rows: list[dict[str, str]], header: list[str] = HEADER) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter=chr(9), lineterminator=chr(10))
        writer.writeheader()
        writer.writerows({key: row[key] for key in header} for row in rows)
    return path


def test_loader_preserves_original_and_normalizes_method_and_booleans(tmp_path: Path) -> None:
    record = LocalKnownInteractionsTsvLoader().load(_write(tmp_path / "x.tsv", [_row()]))[0]
    assert record.detection_method == "Y2H"
    assert record.normalized_detection_method == "yeast_two_hybrid"
    assert record.is_direct is True
    assert record.is_physical is True
    assert record.is_biological is True
    assert record.confidence == 0.8


def test_loader_accepts_reverse_orientation_as_data(tmp_path: Path) -> None:
    record = LocalKnownInteractionsTsvLoader().load(
        _write(tmp_path / "x.tsv", [_row(protein_a_id="CANDIDATE", protein_b_id="QUERY")])
    )[0]
    assert (record.protein_a_id, record.protein_b_id) == ("CANDIDATE", "QUERY")


def test_loader_blank_optional_values_become_none(tmp_path: Path) -> None:
    record = LocalKnownInteractionsTsvLoader().load(
        _write(tmp_path / "x.tsv", [_row(confidence="", is_direct="", detection_method="")])
    )[0]
    assert record.confidence is None
    assert record.is_direct is None
    assert record.normalized_detection_method is None


def test_loader_rejects_missing_required_columns(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="missing columns"):
        LocalKnownInteractionsTsvLoader().load(
            _write(
                tmp_path / "x.tsv",
                [_row()],
                [column for column in HEADER if column != "source_record_id"],
            )
        )


def test_loader_rejects_same_protein_pair(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="must differ"):
        LocalKnownInteractionsTsvLoader().load(
            _write(tmp_path / "x.tsv", [_row(protein_b_id="QUERY")])
        )


def test_loader_rejects_duplicate_source_record_even_for_reverse_pair(tmp_path: Path) -> None:
    rows = [_row(), _row(protein_a_id="CANDIDATE", protein_b_id="QUERY")]
    with pytest.raises(InputValidationError, match="Duplicate known interaction"):
        LocalKnownInteractionsTsvLoader().load(_write(tmp_path / "x.tsv", rows))


def test_loader_rejects_unknown_interaction_type(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="interaction_type"):
        LocalKnownInteractionsTsvLoader().load(
            _write(tmp_path / "x.tsv", [_row(interaction_type="mystery")])
        )


@pytest.mark.parametrize("value", ["maybe", "unknown", "2"])
def test_loader_rejects_unknown_boolean(value: str, tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="boolean"):
        LocalKnownInteractionsTsvLoader().load(_write(tmp_path / "x.tsv", [_row(is_direct=value)]))


@pytest.mark.parametrize("value", ["-0.1", "1.1"])
def test_loader_rejects_confidence_outside_unit_interval(value: str, tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="confidence"):
        LocalKnownInteractionsTsvLoader().load(_write(tmp_path / "x.tsv", [_row(confidence=value)]))


def test_loader_rejects_empty_required_value(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="source"):
        LocalKnownInteractionsTsvLoader().load(_write(tmp_path / "x.tsv", [_row(source="")]))


def test_loader_rejects_invalid_mapping_status(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="identifier_mapping_status"):
        LocalKnownInteractionsTsvLoader().load(
            _write(tmp_path / "x.tsv", [_row(identifier_mapping_status="guess")])
        )


def test_loader_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalKnownInteractionsTsvLoader().load(tmp_path / "missing.tsv")
