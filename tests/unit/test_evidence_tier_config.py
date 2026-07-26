from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
import yaml

from protein_interaction_hunter.config import load_config
from protein_interaction_hunter.exceptions import ConfigurationError


def write_config(source: Path, target: Path, update: Callable[[dict[str, Any]], object]) -> Path:
    data = yaml.safe_load(source.read_text(encoding="utf-8"))
    update(data)
    target.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return target


def test_omitted_evidence_tier_section_defaults_disabled(
    valid_config_path: Path, tmp_path: Path
) -> None:
    path = write_config(
        valid_config_path,
        tmp_path / "omitted.yaml",
        lambda data: data.pop("evidence_tiers"),
    )
    config = load_config(path)
    assert config.evidence_tiers.enabled is False
    assert config.evidence_tiers.tier_1.minimum_score == 75.0


def test_valid_enabled_evidence_tier_config(valid_config_path: Path, tmp_path: Path) -> None:
    def enable(data: dict[str, Any]) -> None:
        data["scoring"]["enabled"] = True
        data["evidence_tiers"]["enabled"] = True

    config = load_config(write_config(valid_config_path, tmp_path / "enabled.yaml", enable))
    assert config.evidence_tiers.enabled is True
    assert config.evidence_tiers.tier_4.maximum_negative_components == 99


@pytest.mark.parametrize(
    ("field", "upper", "lower"),
    [
        ("minimum_score", 20.0, 30.0),
        ("minimum_categories", 1, 2),
        ("minimum_components", 1, 2),
        ("minimum_available_weight", 0.5, 1.0),
        ("minimum_high_specificity_components", 0, 1),
        ("maximum_negative_components", 2, 1),
    ],
)
def test_tier_requirement_order_is_validated(
    valid_config_path: Path,
    tmp_path: Path,
    field: str,
    upper: float,
    lower: float,
) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["tier_1"][field] = upper
        data["evidence_tiers"]["tier_2"][field] = lower

    with pytest.raises(ConfigurationError, match=field):
        load_config(write_config(valid_config_path, tmp_path / f"bad-{field}.yaml", invalidate))


def test_tier_one_high_specificity_gate_cannot_be_disabled(
    valid_config_path: Path, tmp_path: Path
) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["tier_1"]["require_high_specificity_evidence"] = False
        for name in ("tier_2", "tier_3", "tier_4"):
            data["evidence_tiers"][name]["require_high_specificity_evidence"] = False

    with pytest.raises(ConfigurationError, match="minimum_high_specificity_components"):
        load_config(write_config(valid_config_path, tmp_path / "tier1-high.yaml", invalidate))


def test_high_specificity_boolean_order_is_validated(
    valid_config_path: Path, tmp_path: Path
) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["tier_1"]["require_high_specificity_evidence"] = False
        data["evidence_tiers"]["tier_2"]["require_high_specificity_evidence"] = True

    with pytest.raises(ConfigurationError, match="require_high_specificity"):
        load_config(write_config(valid_config_path, tmp_path / "bad-high.yaml", invalidate))


def test_invalid_tier_cap_is_rejected(valid_config_path: Path, tmp_path: Path) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["explicit_conflict_tier_cap"] = "tier_9"

    with pytest.raises(ConfigurationError, match="explicit_conflict_tier_cap"):
        load_config(write_config(valid_config_path, tmp_path / "bad-cap.yaml", invalidate))


def test_enabled_tiers_require_scoring(valid_config_path: Path, tmp_path: Path) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["enabled"] = True

    with pytest.raises(ConfigurationError, match="requires scoring.enabled"):
        load_config(write_config(valid_config_path, tmp_path / "dependency.yaml", invalidate))


def test_enabled_tier_threshold_cannot_exceed_output_scale(
    valid_config_path: Path, tmp_path: Path
) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["scoring"]["enabled"] = True
        data["scoring"]["output_scale"] = 50.0
        data["evidence_tiers"]["enabled"] = True

    with pytest.raises(ConfigurationError, match="output_scale"):
        load_config(write_config(valid_config_path, tmp_path / "scale.yaml", invalidate))


def test_disabled_tiers_do_not_restrict_scoring_output_scale(
    valid_config_path: Path, tmp_path: Path
) -> None:
    def update(data: dict[str, Any]) -> None:
        data["scoring"]["enabled"] = True
        data["scoring"]["output_scale"] = 1.0

    config = load_config(write_config(valid_config_path, tmp_path / "disabled-scale.yaml", update))
    assert config.evidence_tiers.enabled is False
    assert config.scoring.output_scale == 1.0


def test_evidence_tier_unknown_field_is_rejected(valid_config_path: Path, tmp_path: Path) -> None:
    def invalidate(data: dict[str, Any]) -> None:
        data["evidence_tiers"]["tier_1"]["unexpected"] = True

    with pytest.raises(ConfigurationError, match="unexpected"):
        load_config(write_config(valid_config_path, tmp_path / "unknown.yaml", invalidate))
