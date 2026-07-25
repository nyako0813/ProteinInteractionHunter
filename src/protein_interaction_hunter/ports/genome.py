"""Genome input boundary."""

from pathlib import Path
from typing import Protocol

from protein_interaction_hunter.models.genome import GeneCoordinate


class GenomeRepository(Protocol):
    def load(self, path: Path) -> list[GeneCoordinate]:
        """Load normalized coordinates without neighborhood analysis."""
        ...
