"""Public domain model API."""

from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    ContradictionSeverity,
    EvidenceOrigin,
    EvidenceStatus,
    EvidenceTier,
    IdentifierMatchStatus,
    ManualReviewStatus,
    ManualStructurePriority,
    PredictedRelationshipType,
    RunStatus,
)
from protein_interaction_hunter.models.evidence import (
    CandidateEvidenceBundle,
    ContradictionEvidence,
    DomainEvidence,
    EvidenceProvenance,
    FunctionalEvidence,
    FusionEvidence,
    GenomeContextEvidence,
    KnownInteractionEvidence,
    LocalizationEvidence,
    OperonEvidence,
    OrthologRecord,
    PhylogeneticProfileEvidence,
)
from protein_interaction_hunter.models.genome import GeneCoordinate
from protein_interaction_hunter.models.identity import IdentifierAlias, IdentifierResolution
from protein_interaction_hunter.models.protein import (
    CandidateProtein,
    ProteinRecord,
    QueryProtein,
)
from protein_interaction_hunter.models.run import InputFileManifest, RunManifest
from protein_interaction_hunter.models.scoring import CandidateScore
from protein_interaction_hunter.models.structure_queue import StructurePredictionQueueEntry

__all__ = [
    "AnnotationRecord",
    "CandidateDisposition",
    "CandidateEvidenceBundle",
    "CandidateProtein",
    "CandidateScore",
    "ContradictionEvidence",
    "ContradictionSeverity",
    "DomainEvidence",
    "EvidenceOrigin",
    "EvidenceProvenance",
    "EvidenceStatus",
    "EvidenceTier",
    "IdentifierAlias",
    "IdentifierMatchStatus",
    "IdentifierResolution",
    "FunctionalEvidence",
    "FusionEvidence",
    "GeneCoordinate",
    "GenomeContextEvidence",
    "InputFileManifest",
    "KnownInteractionEvidence",
    "LocalizationEvidence",
    "ManualReviewStatus",
    "ManualStructurePriority",
    "OperonEvidence",
    "OrthologRecord",
    "PhylogeneticProfileEvidence",
    "PredictedRelationshipType",
    "ProteinRecord",
    "QueryProtein",
    "RunManifest",
    "RunStatus",
    "StructurePredictionQueueEntry",
]
