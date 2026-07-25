"""Future orthology boundary; no implementation exists in MVP-0."""

from typing import Protocol

from protein_interaction_hunter.models.evidence import OrthologRecord


class OrthologyProvider(Protocol):
    def get_orthologs(self, protein_ids: list[str]) -> list[OrthologRecord]:
        """Return validated orthology evidence from a future adapter."""
        ...
