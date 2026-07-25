"""Manual structure-prediction queue schema; no submission behavior exists."""

from pathlib import Path
from typing import Annotated

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    EvidenceTier,
    ManualReviewStatus,
    ManualStructurePriority,
    PredictedRelationshipType,
)
from protein_interaction_hunter.schemas.versions import SchemaName, schema_version

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class StructurePredictionQueueEntry(StrictModel):
    """A proposal for manual review; component flags are suggestions, not assertions."""

    schema_version: str = Field(
        default_factory=lambda: schema_version(SchemaName.STRUCTURE_PREDICTION_QUEUE)
    )
    rank: int = Field(ge=1)
    query_id: NonEmptyStr
    candidate_id: NonEmptyStr
    candidate_name: NonEmptyStr
    physical_interaction_score: float | None = Field(default=None, ge=0.0, le=1.0)
    functional_association_score: float | None = Field(default=None, ge=0.0, le=1.0)
    evidence_tier: EvidenceTier | None = None
    predicted_relationship_type: PredictedRelationshipType
    manual_structure_priority: ManualStructurePriority
    suggested_stoichiometry: NonEmptyStr | None = None
    include_rna: bool | None = None
    include_dna: bool | None = None
    include_cofactor: bool | None = None
    include_metal: bool | None = None
    include_ligand: bool | None = None
    query_fasta_path: Path
    candidate_fasta_path: Path
    pair_fasta_path: Path
    reason_for_structural_test: NonEmptyStr
    primary_supporting_evidence: NonEmptyStr
    main_contradiction: str = ""
    manual_review_status: ManualReviewStatus = ManualReviewStatus.NOT_REVIEWED
    manual_notes: str = ""
