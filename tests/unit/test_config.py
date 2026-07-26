"""Strict configuration validation."""

from pathlib import Path

import pytest
import yaml

from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import ConfigurationError


def test_valid_config_resolves_paths_relative_to_yaml(valid_config_path: Path) -> None:
    config = load_config(valid_config_path)
    assert (
        config.input.proteome_fasta
        == (valid_config_path.parent / "synthetic_proteome.fasta").resolve()
    )
    assert config.query.protein_ids == ["QUERY_001"]


def test_invalid_config_rejects_multiple_range_errors(fixture_dir: Path) -> None:
    with pytest.raises(ConfigurationError) as error:
        load_config(fixture_dir / "config.invalid.yaml")
    message = str(error.value)
    assert "minimum_length_aa" in message
    assert "automatic_structure_prediction" in message
    assert "workers" in message


def test_unknown_config_field_is_rejected(valid_config_path: Path, tmp_path: Path) -> None:
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


def test_orthology_local_table_path_is_resolved(
    valid_config_path: Path,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    data["orthology"]["source"] = "local_table"
    data["orthology"]["local_table"] = "synthetic_orthology.tsv"

    config_path = valid_config_path.parent / "config.orthology-test.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False),
        encoding="utf-8",
    )

    try:
        config = load_config(config_path)
    finally:
        config_path.unlink(missing_ok=True)

    assert config.orthology.source == "local_table"
    assert (
        config.orthology.local_table
        == (valid_config_path.parent / "synthetic_orthology.tsv").resolve()
    )


def test_phylogenetic_profile_path_and_thresholds_are_resolved(
    valid_config_path: Path,
) -> None:
    config = load_config(valid_config_path)
    assert config.phylogenetic_profile.enabled is True
    assert config.phylogenetic_profile.minimum_shared_species == 2
    assert config.phylogenetic_profile.minimum_informative_species == 3
    assert config.phylogenetic_profile.minimum_profile_similarity == 0.8
    assert (
        config.phylogenetic_profile.local_table
        == (valid_config_path.parent / "synthetic_phylogenetic_profiles.tsv").resolve()
    )


def test_absent_phylogenetic_profile_section_defaults_to_disabled(
    valid_config_path: Path,
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    del data["phylogenetic_profile"]
    path = tmp_path / "without-profile.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    config = load_config(path)
    assert config.phylogenetic_profile.enabled is False
    assert config.phylogenetic_profile.local_table is None


def test_fusion_path_and_thresholds_are_resolved(valid_config_path: Path) -> None:
    config = load_config(valid_config_path)
    assert config.fusion.enabled is True
    assert config.fusion.minimum_supporting_records == 1
    assert config.fusion.minimum_component_coverage == 0.6
    assert config.fusion.maximum_component_overlap_fraction == 0.2
    assert (
        config.fusion.local_table == (valid_config_path.parent / "synthetic_fusions.tsv").resolve()
    )


def test_absent_fusion_section_defaults_to_disabled(
    valid_config_path: Path,
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    del data["fusion"]
    path = tmp_path / "without-fusion.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    config = load_config(path)
    assert config.fusion.enabled is False
    assert config.fusion.local_table is None


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("minimum_supporting_records", 0),
        ("minimum_component_coverage", 1.1),
        ("maximum_component_overlap_fraction", -0.1),
    ],
)
def test_invalid_fusion_threshold_is_rejected(
    valid_config_path: Path,
    tmp_path: Path,
    field: str,
    value: int | float,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    data["fusion"][field] = value
    path = tmp_path / f"invalid-fusion-{field}.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ConfigurationError, match=field):
        load_config(path)


def test_known_interactions_path_and_policy_are_resolved(valid_config_path: Path) -> None:
    config = load_config(valid_config_path)
    assert config.known_interactions.enabled is True
    assert config.known_interactions.source == "local_table"
    assert config.known_interactions.minimum_supporting_records == 1
    assert config.known_interactions.minimum_direct_records == 1
    assert config.known_interactions.minimum_confidence == 0.5
    assert config.known_interactions.excluded_evidence_methods == ["database_inference"]
    assert (
        config.known_interactions.local_table
        == (valid_config_path.parent / "synthetic_known_interactions.tsv").resolve()
    )


def test_absent_known_interactions_section_defaults_to_disabled(
    valid_config_path: Path,
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    del data["known_interactions"]
    path = tmp_path / "without-known.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    config = load_config(path)
    assert config.known_interactions.enabled is False
    assert config.known_interactions.local_table is None


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("minimum_supporting_records", 0),
        ("minimum_direct_records", -1),
        ("minimum_confidence", -0.1),
        ("minimum_confidence", 1.1),
    ],
)
def test_invalid_known_interactions_threshold_is_rejected(
    valid_config_path: Path,
    tmp_path: Path,
    field: str,
    value: int | float,
) -> None:
    data = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    data["known_interactions"][field] = value
    path = tmp_path / f"invalid-known-{field}.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ConfigurationError, match=field):
        load_config(path)
