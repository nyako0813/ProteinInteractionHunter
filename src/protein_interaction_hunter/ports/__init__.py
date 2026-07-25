"""Application port interfaces."""

from protein_interaction_hunter.ports.annotation import AnnotationProvider
from protein_interaction_hunter.ports.genome import GenomeRepository
from protein_interaction_hunter.ports.interaction import KnownInteractionProvider
from protein_interaction_hunter.ports.orthology import OrthologyProvider
from protein_interaction_hunter.ports.output import (
    EvidenceBundleWriter,
    RunManifestWriter,
    StructureQueueWriter,
)
from protein_interaction_hunter.ports.proteome import ProteomeRepository

__all__ = [
    "AnnotationProvider",
    "EvidenceBundleWriter",
    "GenomeRepository",
    "KnownInteractionProvider",
    "OrthologyProvider",
    "ProteomeRepository",
    "RunManifestWriter",
    "StructureQueueWriter",
]
