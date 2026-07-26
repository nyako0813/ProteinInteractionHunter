"""Evidence records capable of representing unavailable and unrun engines."""

from typing import Annotated, Any

from pydantic import Field, StringConstraints

from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    ContextCompleteness,
    ContradictionSeverity,
    CoordinatePosition,
    EvidenceOrigin,
    EvidenceStatus,
    EvidenceTier,
    OperonProxyStatus,
    PredictedRelationshipType,
    RelativePosition,
    StrandRelationship,
    TranscriptionalOrder,
)
from protein_interaction_hunter.models.protein import CandidateProtein
from protein_interaction_hunter.models.scoring import CandidateScore
from protein_interaction_hunter.schemas.versions import SchemaName, schema_version

NonEmptyStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class EvidenceProvenance(StrictModel):
    source_name: NonEmptyStr
    source_version: NonEmptyStr | None = None
    source_record_id: NonEmptyStr | None = None
    method: NonEmptyStr | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)


class BaseEvidence(StrictModel):
    status: EvidenceStatus = EvidenceStatus.NOT_RUN
    origin: EvidenceOrigin = EvidenceOrigin.INFERRED
    quality: float | None = Field(default=None, ge=0.0, le=1.0)
    provenance: list[EvidenceProvenance] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


class GenomeContextEvidence(BaseEvidence):
    calculation_rule_version: NonEmptyStr | None = None
    same_contig: bool | None = None
    same_seqid: bool | None = None
    query_contig: NonEmptyStr | None = None
    candidate_contig: NonEmptyStr | None = None
    query_start: int | None = Field(default=None, ge=1)
    query_end: int | None = Field(default=None, ge=1)
    query_strand: str | None = Field(default=None, pattern=r"^[+\-?]$")
    candidate_start: int | None = Field(default=None, ge=1)
    candidate_end: int | None = Field(default=None, ge=1)
    candidate_strand: str | None = Field(default=None, pattern=r"^[+\-?]$")
    strand_relationship: StrandRelationship | None = None
    distance_bp: int | None = Field(default=None, ge=0)
    edge_to_edge_distance_bp: int | None = Field(default=None, ge=0)
    overlap_bp: int | None = Field(default=None, ge=0)
    relative_position: RelativePosition | None = None
    coordinate_position: CoordinatePosition | None = None
    intervening_feature_count: int | None = Field(default=None, ge=0)
    intervening_gene_count: int | None = Field(default=None, ge=0)
    query_feature_index: int | None = Field(default=None, ge=0)
    candidate_feature_index: int | None = Field(default=None, ge=0)
    feature_index_delta: int | None = Field(default=None, ge=0)
    within_neighborhood_window: bool | None = None
    within_neighborhood_gene_count: bool | None = None
    query_left_edge_distance_bp: int | None = Field(default=None, ge=0)
    query_right_edge_distance_bp: int | None = Field(default=None, ge=0)
    candidate_left_edge_distance_bp: int | None = Field(default=None, ge=0)
    candidate_right_edge_distance_bp: int | None = Field(default=None, ge=0)
    query_distance_to_contig_left_edge: int | None = Field(default=None, ge=0)
    query_distance_to_contig_right_edge: int | None = Field(default=None, ge=0)
    candidate_distance_to_contig_left_edge: int | None = Field(default=None, ge=0)
    candidate_distance_to_contig_right_edge: int | None = Field(default=None, ge=0)
    context_completeness: ContextCompleteness | None = None
    strand_relation: str | None = None
    boundary_flags: list[str] = Field(default_factory=list)


class OperonEvidence(BaseEvidence):
    calculation_rule_version: NonEmptyStr | None = None
    same_contig: bool | None = None
    same_strand: bool | None = None
    is_adjacent: bool | None = None
    intergenic_distance_bp: int | None = None
    overlap_bp: int | None = Field(default=None, ge=0)
    intervening_gene_count: int | None = Field(default=None, ge=0)
    transcriptional_order: TranscriptionalOrder | None = None
    maximum_intergenic_distance_bp: int | None = Field(default=None, ge=0)
    passes_distance_threshold: bool | None = None
    proxy_status: OperonProxyStatus | None = None
    proxy_rule_id: NonEmptyStr | None = None
    supporting_conditions: list[str] = Field(default_factory=list)
    conflicting_conditions: list[str] = Field(default_factory=list)
    support: float | None = Field(default=None, ge=0.0, le=1.0)


class OrthologRecord(BaseEvidence):
    protein_id: NonEmptyStr
    reference_id: NonEmptyStr
    ortholog_id: NonEmptyStr | None = None
    identity: float | None = Field(default=None, ge=0.0, le=1.0)
    query_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    subject_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    evalue: float | None = Field(default=None, ge=0.0)
    orthogroup: NonEmptyStr | None = None
    paralog_ambiguity: bool = False


class PhylogeneticProfileEvidence(BaseEvidence):
    informative_taxa: int | None = Field(default=None, ge=0)
    missing_rate: float | None = Field(default=None, ge=0.0, le=1.0)
    raw_similarity: float | None = Field(default=None, ge=-1.0, le=1.0)
    corrected_similarity: float | None = Field(default=None, ge=-1.0, le=1.0)


class DomainEvidence(BaseEvidence):
    calculation_rule_version: NonEmptyStr | None = None
    protein_id: NonEmptyStr
    source: NonEmptyStr | None = None
    accession: NonEmptyStr | None = None
    name: NonEmptyStr | None = None
    start: int | None = Field(default=None, ge=1)
    end: int | None = Field(default=None, ge=1)
    architecture_index: int | None = Field(default=None, ge=0)
    role: NonEmptyStr | None = None
    pair_rule_id: NonEmptyStr | None = None
    paired_protein_id: NonEmptyStr | None = None
    paired_accession: NonEmptyStr | None = None
    is_shared: bool | None = None
    pair_matched: bool | None = None
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    ruleset_path: NonEmptyStr | None = None


class FunctionalEvidence(BaseEvidence):
    calculation_rule_version: NonEmptyStr | None = None
    query_role: NonEmptyStr | None = None
    candidate_role: NonEmptyStr | None = None
    relationship_hint: PredictedRelationshipType | None = None
    rule_id: NonEmptyStr | None = None
    query_matched_terms: list[str] = Field(default_factory=list)
    candidate_matched_terms: list[str] = Field(default_factory=list)
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    query_annotation_text: NonEmptyStr | None = None
    candidate_annotation_text: NonEmptyStr | None = None
    ruleset_path: NonEmptyStr | None = None
    matched: bool | None = None


class LocalizationEvidence(BaseEvidence):
    calculation_rule_version: NonEmptyStr | None = None
    protein_id: NonEmptyStr
    compartment: NonEmptyStr | None = None
    signal_peptide: bool | None = None
    transmembrane_helices: int | None = Field(default=None, ge=0)
    topology: NonEmptyStr | None = None
    compatibility: bool | None = None
    query_compartment: NonEmptyStr | None = None
    candidate_compartment: NonEmptyStr | None = None
    localization_annotation: NonEmptyStr | None = None
    transmembrane_annotation: NonEmptyStr | None = None
    matched_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    rule_id: NonEmptyStr | None = None
    annotation_source: NonEmptyStr | None = None
    annotation_confidence: float | None = Field(
        default=None,
        ge=0.0,
        le=1.0,
    )


class FusionEvidence(BaseEvidence):
    fused_protein_id: NonEmptyStr | None = None
    reference_id: NonEmptyStr | None = None
    query_region: tuple[int, int] | None = None
    candidate_region: tuple[int, int] | None = None
    taxonomic_support_count: int | None = Field(default=None, ge=0)


class KnownInteractionEvidence(BaseEvidence):
    protein_a: NonEmptyStr
    protein_b: NonEmptyStr
    interaction_type: NonEmptyStr | None = None
    detection_method: NonEmptyStr | None = None
    publication_id: NonEmptyStr | None = None
    taxonomic_distance: NonEmptyStr | None = None


class ContradictionEvidence(BaseEvidence):
    contradiction_type: NonEmptyStr
    severity: ContradictionSeverity
    penalty: float | None = Field(default=None, ge=0.0, le=1.0)
    hard_exclusion: bool = False
    explanation: NonEmptyStr


class CandidateEvidenceBundle(StrictModel):
    schema_version: str = Field(
        default_factory=lambda: schema_version(SchemaName.CANDIDATE_EVIDENCE_BUNDLE)
    )
    run_id: NonEmptyStr
    query_id: NonEmptyStr
    candidate_id: NonEmptyStr
    candidate: CandidateProtein | None = None
    candidate_disposition: CandidateDisposition
    predicted_relationship_type: PredictedRelationshipType
    evidence_tier: EvidenceTier | None = None
    genome_context: list[GenomeContextEvidence] = Field(default_factory=list)
    operon: list[OperonEvidence] = Field(default_factory=list)
    orthology: list[OrthologRecord] = Field(default_factory=list)
    phylogenetic_profile: list[PhylogeneticProfileEvidence] = Field(default_factory=list)
    domains: list[DomainEvidence] = Field(default_factory=list)
    functional: list[FunctionalEvidence] = Field(default_factory=list)
    localization: list[LocalizationEvidence] = Field(default_factory=list)
    fusion: list[FusionEvidence] = Field(default_factory=list)
    known_interactions: list[KnownInteractionEvidence] = Field(default_factory=list)
    contradictions: list[ContradictionEvidence] = Field(default_factory=list)
    score: CandidateScore = Field(default_factory=CandidateScore)
    engine_statuses: dict[str, EvidenceStatus] = Field(default_factory=dict)
    provenance: list[EvidenceProvenance] = Field(default_factory=list)
    policy_settings: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
