"""Load and validate domain-pair YAML rules."""

from pathlib import Path

import yaml
from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.domain_rules import (
    DomainPairRuleset,
)


class LocalDomainRulesLoader:
    def load(self, path: Path) -> DomainPairRuleset:
        rules_path = path.expanduser().resolve()

        if not rules_path.is_file():
            raise InputValidationError(
                f"Domain pair rules not found: {rules_path}"
            )

        try:
            payload = yaml.safe_load(
                rules_path.read_text(encoding="utf-8")
            )
        except yaml.YAMLError as exc:
            raise InputValidationError(
                f"Invalid domain pair YAML: {rules_path}: {exc}"
            ) from exc

        if not isinstance(payload, dict):
            raise InputValidationError(
                f"Domain pair rules must be a mapping: {rules_path}"
            )

        try:
            return DomainPairRuleset.model_validate(payload)
        except ValidationError as exc:
            raise InputValidationError(
                f"Invalid domain pair rules: {rules_path}: {exc}"
            ) from exc