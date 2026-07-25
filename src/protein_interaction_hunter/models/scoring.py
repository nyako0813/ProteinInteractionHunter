"""Score container only; no score or tier calculation is performed."""

from pydantic import Field

from protein_interaction_hunter.models.base import StrictModel

OptionalScore = float | None


class CandidateScore(StrictModel):
    physical_interaction_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    functional_association_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    gene_context_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    evolutionary_coupling_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    annotation_confidence_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    contradiction_penalty: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    evidence_completeness: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    total_ranking_score: OptionalScore = Field(default=None, ge=0.0, le=1.0)
    calculation_trace: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
