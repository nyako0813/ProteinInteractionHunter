"""Run manifest and input fingerprint models."""

from datetime import datetime
from pathlib import Path
from typing import Annotated, Any

from pydantic import Field, StringConstraints, field_validator, model_validator

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import RunStatus
from protein_interaction_hunter.schemas.versions import SchemaName, schema_version

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]
Sha256 = Annotated[str, StringConstraints(pattern=r"^[0-9a-f]{64}$")]


class InputFileManifest(StrictModel):
    logical_name: NonEmptyStr
    path: Path
    sha256: Sha256 | None = None
    size_bytes: int | None = Field(default=None, ge=0)
    modified_time: datetime | None = None
    exists: bool
    required: bool

    @field_validator("modified_time")
    @classmethod
    def require_aware_timestamp(cls, value: datetime | None) -> datetime | None:
        if value is not None and (value.tzinfo is None or value.utcoffset() is None):
            raise ValueError("modified_time must be timezone-aware")
        return value


class RunManifest(StrictModel):
    schema_version: str = Field(default_factory=lambda: schema_version(SchemaName.RUN_MANIFEST))
    run_id: NonEmptyStr
    run_name: NonEmptyStr
    started_at: datetime
    completed_at: datetime | None = None
    status: RunStatus
    command_line: list[str] = Field(default_factory=list)
    working_directory: Path
    config_path: Path
    config_sha256: Sha256
    config_snapshot_path: Path | None = None
    input_files: list[InputFileManifest] = Field(default_factory=list)
    python_version: NonEmptyStr
    platform: NonEmptyStr
    package_version: NonEmptyStr
    dependency_versions: dict[str, str] = Field(default_factory=dict)
    git_commit: NonEmptyStr | None = None
    random_seed: int
    warnings: list[str] = Field(default_factory=list)
    incomplete_evidence_flags: list[str] = Field(default_factory=list)
    normalization_rule_version: NonEmptyStr | None = None
    gene_context_rule_version: NonEmptyStr | None = None
    orthology_rule_version: NonEmptyStr | None = None
    phylogenetic_profile_rule_version: NonEmptyStr | None = None
    fusion_rule_version: NonEmptyStr | None = None
    known_interactions_rule_version: NonEmptyStr | None = None
    scoring_rule_version: NonEmptyStr | None = None
    scoring_config_snapshot: dict[str, Any] = Field(default_factory=dict)
    policy_settings: dict[str, str | int | bool] = Field(default_factory=dict)
    parser_warnings: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_timestamps(self) -> "RunManifest":
        for field_name, value in (
            ("started_at", self.started_at),
            ("completed_at", self.completed_at),
        ):
            if value is not None and (value.tzinfo is None or value.utcoffset() is None):
                raise ValueError(f"{field_name} must be timezone-aware")
        if self.completed_at is not None and self.completed_at < self.started_at:
            raise ValueError("completed_at cannot be earlier than started_at")
        return self
