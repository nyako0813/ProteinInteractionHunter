from pathlib import Path

import pytest

from protein_interaction_hunter.adapters.local.domain_rules import (
    LocalDomainRulesLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError


def test_load_valid_domain_rules(fixture_dir: Path) -> None:
    rules = LocalDomainRulesLoader().load(
        fixture_dir / "rules" / "domain_pairs.v1.yaml"
    )

    assert rules.ruleset_version == "mvp1e-domain-pair-v1"
    assert len(rules.roles) == 4
    assert len(rules.pair_rules) == 2
    assert rules.pair_rules[0].allow_shared_accession is False


def test_missing_domain_rules_file_is_rejected(
    tmp_path: Path,
) -> None:
    with pytest.raises(InputValidationError, match="not found"):
        LocalDomainRulesLoader().load(
            tmp_path / "missing.yaml"
        )


def test_invalid_domain_rules_are_rejected(
    tmp_path: Path,
) -> None:
    path = tmp_path / "invalid.yaml"
    path.write_text(
        (
            "ruleset_version: test\n"
            "roles: []\n"
            "pair_rules: []\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        InputValidationError,
        match="Invalid domain pair",
    ):
        LocalDomainRulesLoader().load(path)