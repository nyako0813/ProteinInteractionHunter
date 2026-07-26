from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.phylogenetic_profile import (
    LocalPhylogeneticProfileTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError


def test_loads_profile_and_presence_representations(fixture_dir: Path) -> None:
    records = LocalPhylogeneticProfileTsvLoader().load(
        fixture_dir / "synthetic_phylogenetic_profiles.tsv"
    )
    assert len(records) == 15
    assert [record.presence for record in records[:5]] == [True, True, False, False, None]
    assert records[0].taxonomic_group == "archaea"
    assert records[0].source == "synthetic_profiles"


def test_all_allowed_presence_representations(tmp_path: Path) -> None:
    path = tmp_path / "presence.tsv"
    path.write_text(
        "protein_id\tspecies_id\tpresence\n"
        "P1\tS1\tTRUE\nP1\tS2\t1\nP1\tS3\tYes\n"
        "P1\tS4\tFALSE\nP1\tS5\t0\nP1\tS6\tNo\nP1\tS7\t\n",
        encoding="utf-8",
    )
    records = LocalPhylogeneticProfileTsvLoader().load(path)
    assert [record.presence for record in records] == [True, True, True, False, False, False, None]


def test_missing_required_column_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "missing.tsv"
    path.write_text("protein_id\tspecies_id\nP1\tS1\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="missing columns.*presence"):
        LocalPhylogeneticProfileTsvLoader().load(path)


def test_invalid_presence_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "invalid.tsv"
    path.write_text(
        "protein_id\tspecies_id\tpresence\nP1\tS1\tmaybe\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="Invalid phylogenetic profile"):
        LocalPhylogeneticProfileTsvLoader().load(path)


def test_duplicate_protein_species_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.tsv"
    path.write_text(
        "protein_id\tspecies_id\tpresence\nP1\tS1\ttrue\nP1\tS1\tfalse\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="Duplicate phylogenetic profile"):
        LocalPhylogeneticProfileTsvLoader().load(path)


def test_missing_file_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalPhylogeneticProfileTsvLoader().load(tmp_path / "missing.tsv")
