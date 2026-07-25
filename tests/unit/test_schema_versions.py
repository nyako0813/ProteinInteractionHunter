"""Central schema version registry."""

from protein_interaction_hunter.models import (
    CandidateDisposition,
    CandidateEvidenceBundle,
    PredictedRelationshipType,
    RunManifest,
    StructurePredictionQueueEntry,
)
from protein_interaction_hunter.schemas.versions import SCHEMA_VERSIONS, SchemaName


def test_schema_versions_are_central_and_applied() -> None:
    bundle = CandidateEvidenceBundle(
        run_id="r",
        query_id="q",
        candidate_id="c",
        candidate_disposition=CandidateDisposition.INCLUDED,
        predicted_relationship_type=PredictedRelationshipType.INSUFFICIENT_EVIDENCE,
    )
    assert bundle.schema_version == SCHEMA_VERSIONS[SchemaName.CANDIDATE_EVIDENCE_BUNDLE]
    assert RunManifest.model_fields["schema_version"].default_factory is not None
    assert StructurePredictionQueueEntry.model_fields["schema_version"].default_factory is not None
