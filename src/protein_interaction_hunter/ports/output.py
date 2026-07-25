"""Serialization boundaries for canonical and derived artifacts."""

from collections.abc import Sequence
from pathlib import Path
from typing import Protocol

from protein_interaction_hunter.models.evidence import CandidateEvidenceBundle
from protein_interaction_hunter.models.run import RunManifest
from protein_interaction_hunter.models.structure_queue import StructurePredictionQueueEntry


class EvidenceBundleWriter(Protocol):
    def write(self, records: Sequence[CandidateEvidenceBundle], path: Path) -> Path:
        """Write canonical evidence records."""
        ...


class RunManifestWriter(Protocol):
    def write(self, manifest: RunManifest, path: Path) -> Path:
        """Write a run manifest."""
        ...


class StructureQueueWriter(Protocol):
    def write(self, entries: Sequence[StructurePredictionQueueEntry], path: Path) -> Path:
        """Write the flat manual-review queue."""
        ...
