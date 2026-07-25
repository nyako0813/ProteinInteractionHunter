"""Single source of truth for public artifact schema versions."""

from enum import StrEnum


class SchemaName(StrEnum):
    """Versioned public schemas."""

    CANDIDATE_EVIDENCE_BUNDLE = "candidate_evidence_bundle"
    RUN_MANIFEST = "run_manifest"
    STRUCTURE_PREDICTION_QUEUE = "structure_prediction_queue"


SCHEMA_VERSIONS: dict[SchemaName, str] = {
    SchemaName.CANDIDATE_EVIDENCE_BUNDLE: "1.2",
    SchemaName.RUN_MANIFEST: "1.2",
    SchemaName.STRUCTURE_PREDICTION_QUEUE: "1.0",
}


def schema_version(name: SchemaName) -> str:
    """Return the registered version for one public schema."""
    return SCHEMA_VERSIONS[name]
