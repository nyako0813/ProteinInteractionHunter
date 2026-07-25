from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.functional_rules import (
    LocalFunctionalRulesLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import PredictedRelationshipType


def test_load_valid_functional_rules(fixture_dir: Path) -> None:
    rules = LocalFunctionalRulesLoader().load(
        fixture_dir / "rules" / "functional_complementarity.v1.yaml"
    )

    assert rules.ruleset_version == "mvp1d-functional-complementarity-v1"
    assert len(rules.roles) == 3
    assert len(rules.pair_rules) == 2
    assert (
        rules.pair_rules[0].relationship_hint
        is PredictedRelationshipType.ACCESSORY_FACTOR
    )


def test_missing_functional_rules_file_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError, match="rules not found"):
        LocalFunctionalRulesLoader().load(tmp_path / "missing.yaml")


def test_invalid_functional_rules_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "invalid.yaml"
    path.write_text(
        "ruleset_version: test\nroles: []\npair_rules: []\n",
        encoding="utf-8",
    )

    with pytest.raises(InputValidationError, match="Invalid functional"):
        LocalFunctionalRulesLoader().load(path)