"""Annotation input boundary."""

from pathlib import Path
from typing import Protocol

from protein_interaction_hunter.models.annotation import AnnotationRecord


class AnnotationProvider(Protocol):
    def load(self, path: Path) -> list[AnnotationRecord]:
        """Load local annotations without network access."""
        ...
