from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.domains import (
    LocalDomainTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError


def test_load_valid_domain_table(fixture_dir: Path) -> None:
    records = LocalDomainTsvLoader().load(
        fixture_dir / "synthetic_domains.tsv"
    )

    assert len(records) == 5
    assert records[0].protein_id == "QUERY_001"
    assert records[0].accession == "PF00001"
    assert records[0].architecture_index == 0


def test_missing_domain_table_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalDomainTsvLoader().load(tmp_path / "missing.tsv")


def test_invalid_domain_coordinates_are_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "invalid.tsv"
    path.write_text(
        (
            "protein_id\tsource\taccession\tname\tstart\tend\t"
            "architecture_index\n"
            "P1\tPfam\tPF00001\tBad domain\t50\t10\t0\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(InputValidationError, match="Invalid domain"):
        LocalDomainTsvLoader().load(path)


def test_duplicate_domain_annotation_is_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "duplicate.tsv"
    row = "P1\tPfam\tPF00001\tDomain\t1\t20\t0\n"
    path.write_text(
        (
            "protein_id\tsource\taccession\tname\tstart\tend\t"
            "architecture_index\n"
            + row
            + row
        ),
        encoding="utf-8",
    )

    with pytest.raises(InputValidationError, match="Duplicate domain"):
        LocalDomainTsvLoader().load(path)