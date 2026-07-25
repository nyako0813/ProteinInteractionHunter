"""Future known-interaction boundary; no implementation exists in MVP-0."""

from typing import Protocol

from protein_interaction_hunter.models.evidence import KnownInteractionEvidence


class KnownInteractionProvider(Protocol):
    def get_interactions(
        self,
        query_id: str,
        candidate_ids: list[str],
    ) -> list[KnownInteractionEvidence]:
        """Return validated evidence from a future provider."""
        ...
