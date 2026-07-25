"""Load and validate functional-complementarity YAML rules."""

from pathlib import Path

import yaml
from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.functional_rules import (
    FunctionalComplementarityRuleset,
)


class LocalFunctionalRulesLoader:
    def load(self, path: Path) -> FunctionalComplementarityRuleset:
        rules_path = path.expanduser().resolve()

        if not rules_path.is_file():
            raise InputValidationError(
                f"Functional complementarity rules not found: {rules_path}"
            )

        try:
            payload = yaml.safe_load(rules_path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            raise InputValidationError(
                f"Invalid functional complementarity YAML: {rules_path}: {exc}"
            ) from exc

        if not isinstance(payload, dict):
            raise InputValidationError(
                f"Functional complementarity rules must be a mapping: {rules_path}"
            )

        try:
            return FunctionalComplementarityRuleset.model_validate(payload)
        except ValidationError as exc:
            raise InputValidationError(
                f"Invalid functional complementarity rules: {rules_path}: {exc}"
            ) from exc