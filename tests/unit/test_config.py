"""Strict configuration validation."""

from pathlib import Path

import pytest
import yaml

from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import ConfigurationError


def test_valid_config_resolves_paths_relative_to_yaml(valid_config_path: Path) -> None:
    config = load_config(valid_config_path)
    assert config.input.proteome_fasta == (
        valid_config_path.parent / "synthetic_proteome.fasta"
    ).resolve()
    assert config.query.protein_ids == ["QUERY_001"]


def test_invalid_config_rejects_multiple_range_errors(fixture_dir: Path) -> None:
    with pytest.raises(ConfigurationError) as error:
        load_config(fixture_dir / "config.invalid.yaml")
    message = str(error.value)
    assert "minimum_length_aa" in message
    assert "automatic_structure_prediction" in message
    assert "workers" in message


def test_unknown_config_field_is_rejected(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    raw["project"]["unexpected"] = "forbidden"
    path = tmp_path / "unknown.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="unexpected"):
        load_config(path)


def test_empty_query_id_is_rejected(valid_config_path: Path, tmp_path: Path) -> None:
    raw = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    raw["query"]["protein_ids"] = [""]
    path = tmp_path / "empty-query.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="query protein IDs must not be empty"):
        load_config(path)


def test_automatic_structure_prediction_true_is_rejected(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    raw["structure_prediction_queue"]["automatic_structure_prediction"] = True
    path = tmp_path / "automatic.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="automatic_structure_prediction"):
        load_config(path)
