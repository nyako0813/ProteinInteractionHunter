"""Proteome input boundary."""

from pathlib import Path
from typing import Protocol

from protein_interaction_hunter.models.protein import ProteinRecord


class ProteomeRepository(Protocol):
    def load(self, path: Path) -> list[ProteinRecord]:
        """Load and validate normalized protein records."""
        ...
