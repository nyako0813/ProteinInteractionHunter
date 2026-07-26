from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.fusion import LocalFusionTsvLoader
from protein_interaction_hunter.exceptions import InputValidationError

COLUMNS = (
    "query_protein_id",
    "candidate_protein_id",
    "fusion_protein_id",
    "reference_organism",
    "query_component_start",
    "query_component_end",
    "candidate_component_start",
    "candidate_component_end",
    "fusion_protein_length",
    "query_component_coverage",
    "candidate_component_coverage",
    "source_record_id",
)
HEADER = chr(9).join(COLUMNS) + chr(10)


def row(*values: object) -> str:
    return chr(9).join(str(value) for value in values) + chr(10)


def write_row(tmp_path: Path, values: tuple[object, ...]) -> Path:
    path = tmp_path / "fusion.tsv"
    path.write_text(HEADER + row(*values), encoding="utf-8")
    return path


def test_loads_valid_fusion_observations(fixture_dir: Path) -> None:
    records = LocalFusionTsvLoader().load(fixture_dir / "synthetic_fusions.tsv")
    assert len(records) == 4
    assert records[0].fusion_protein_id == "FUSION_A"
    assert records[0].component_overlap_length == 0
    assert records[0].component_overlap_fraction == 0.0
    assert records[0].query_component_coverage == 0.8


def test_optional_numeric_fields_can_be_blank(fixture_dir: Path) -> None:
    records = LocalFusionTsvLoader().load(fixture_dir / "synthetic_fusions.tsv")
    record = records[2]
    assert record.query_component_identity is None
    assert record.evalue_query is None
    assert record.query_component_reference_id is None


def test_coordinate_start_must_be_positive(tmp_path: Path) -> None:
    path = write_row(tmp_path, ("Q", "C", "F", "Org", 0, 10, 20, 30, 40, 0.8, 0.8, "r1"))
    with pytest.raises(InputValidationError, match="Invalid fusion observation"):
        LocalFusionTsvLoader().load(path)


def test_component_cannot_exceed_fusion_length(tmp_path: Path) -> None:
    path = write_row(tmp_path, ("Q", "C", "F", "Org", 1, 10, 20, 41, 40, 0.8, 0.8, "r1"))
    with pytest.raises(InputValidationError, match="exceeds fusion protein length"):
        LocalFusionTsvLoader().load(path)


def test_end_before_start_is_rejected(tmp_path: Path) -> None:
    path = write_row(tmp_path, ("Q", "C", "F", "Org", 10, 1, 20, 30, 40, 0.8, 0.8, "r1"))
    with pytest.raises(InputValidationError, match="end must be >= start"):
        LocalFusionTsvLoader().load(path)


def test_invalid_coverage_is_rejected(tmp_path: Path) -> None:
    path = write_row(tmp_path, ("Q", "C", "F", "Org", 1, 10, 20, 30, 40, 1.1, 0.8, "r1"))
    with pytest.raises(InputValidationError, match="Invalid fusion observation"):
        LocalFusionTsvLoader().load(path)


def test_identical_query_and_candidate_are_rejected(tmp_path: Path) -> None:
    path = write_row(tmp_path, ("Q", "Q", "F", "Org", 1, 10, 20, 30, 40, 0.8, 0.8, "r1"))
    with pytest.raises(InputValidationError, match="protein IDs must differ"):
        LocalFusionTsvLoader().load(path)


def test_swapped_duplicate_record_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.tsv"
    path.write_text(
        HEADER
        + row("Q", "C", "F", "Org", 1, 10, 20, 30, 40, 0.8, 0.8, "r1")
        + row("C", "Q", "F", "Org", 20, 30, 1, 10, 40, 0.8, 0.8, "r2"),
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="Duplicate fusion observation"):
        LocalFusionTsvLoader().load(path)


def test_missing_required_column_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "missing.tsv"
    path.write_text("query_protein_id candidate_protein_id", encoding="utf-8")
    with pytest.raises(InputValidationError, match="missing columns"):
        LocalFusionTsvLoader().load(path)


def test_missing_file_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalFusionTsvLoader().load(tmp_path / "missing.tsv")
