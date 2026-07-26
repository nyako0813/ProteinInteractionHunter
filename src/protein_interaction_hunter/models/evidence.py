"""Evidence records capable of representing unavailable and unrun engines."""

from typing import Annotated, Any, Literal

from pydantic import Field, StringConstraints, model_validator

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
    calculation_rule_version: NonEmptyStr | None = None
    protein_id: NonEmptyStr
    reference_id: NonEmptyStr
    ortholog_id: NonEmptyStr | None = None
    identity: float | None = Field(default=None, ge=0.0, le=1.0)
    query_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    subject_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    evalue: float | None = Field(default=None, ge=0.0)
    orthogroup: NonEmptyStr | None = None
    paralog_ambiguity: bool = False
    reference_organism: NonEmptyStr | None = None
    relationship: NonEmptyStr | None = None
    paired_protein_id: NonEmptyStr | None = None
    paired_reference_id: NonEmptyStr | None = None
    paired_ortholog_id: NonEmptyStr | None = None
    paired_orthogroup: NonEmptyStr | None = None
    shared_orthogroup: bool | None = None
    pair_supported: bool | None = None
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    source: NonEmptyStr | None = None
    source_record_id: NonEmptyStr | None = None


class PhylogeneticProfileEvidence(BaseEvidence):
    query_protein_id: NonEmptyStr
    candidate_protein_id: NonEmptyStr
    informative_species_count: int | None = Field(default=None, ge=0)
    shared_presence_count: int | None = Field(default=None, ge=0)
    shared_absence_count: int | None = Field(default=None, ge=0)
    discordant_count: int | None = Field(default=None, ge=0)
    unknown_count: int | None = Field(default=None, ge=0)
    profile_similarity: float | None = Field(default=None, ge=0.0, le=1.0)
    pair_supported: bool | None = None
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    calculation_rule_version: NonEmptyStr | None = None
    source: NonEmptyStr | None = None
    source_record_id: NonEmptyStr | None = None


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


class FusionObservation(BaseEvidence):
    """One fusion observation using 1-based inclusive component coordinates."""

    query_protein_id: NonEmptyStr
    candidate_protein_id: NonEmptyStr
    fusion_protein_id: NonEmptyStr
    reference_organism: NonEmptyStr
    query_component_reference_id: NonEmptyStr | None = None
    candidate_component_reference_id: NonEmptyStr | None = None
    query_component_start: int = Field(ge=1)
    query_component_end: int = Field(ge=1)
    candidate_component_start: int = Field(ge=1)
    candidate_component_end: int = Field(ge=1)
    fusion_protein_length: int = Field(ge=1)
    query_component_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    candidate_component_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    query_component_identity: float | None = Field(default=None, ge=0.0, le=1.0)
    candidate_component_identity: float | None = Field(default=None, ge=0.0, le=1.0)
    evalue_query: float | None = Field(default=None, ge=0.0)
    evalue_candidate: float | None = Field(default=None, ge=0.0)
    component_overlap_length: int | None = Field(default=None, ge=0)
    component_overlap_fraction: float | None = Field(default=None, ge=0.0, le=1.0)
    source: NonEmptyStr | None = None
    source_record_id: NonEmptyStr | None = None

    @model_validator(mode="after")
    def validate_coordinates_and_overlap(self) -> "FusionObservation":
        if self.query_protein_id == self.candidate_protein_id:
            raise ValueError("query and candidate protein IDs must differ")
        if self.query_component_end < self.query_component_start:
            raise ValueError("query component end must be >= start")
        if self.candidate_component_end < self.candidate_component_start:
            raise ValueError("candidate component end must be >= start")
        if self.query_component_end > self.fusion_protein_length:
            raise ValueError("query component exceeds fusion protein length")
        if self.candidate_component_end > self.fusion_protein_length:
            raise ValueError("candidate component exceeds fusion protein length")
        query_length = self.query_component_end - self.query_component_start + 1
        candidate_length = self.candidate_component_end - self.candidate_component_start + 1
        overlap = max(
            0,
            min(self.query_component_end, self.candidate_component_end)
            - max(self.query_component_start, self.candidate_component_start)
            + 1,
        )
        overlap_fraction = overlap / min(query_length, candidate_length)
        if self.component_overlap_length is not None and self.component_overlap_length != overlap:
            raise ValueError("component overlap length is inconsistent with coordinates")
        if (
            self.component_overlap_fraction is not None
            and abs(self.component_overlap_fraction - overlap_fraction) > 1e-12
        ):
            raise ValueError("component overlap fraction is inconsistent with coordinates")
        object.__setattr__(self, "component_overlap_length", overlap)
        object.__setattr__(self, "component_overlap_fraction", overlap_fraction)
        return self


class FusionEvidence(BaseEvidence):
    query_protein_id: NonEmptyStr
    candidate_protein_id: NonEmptyStr
    supporting_record_count: int = Field(ge=0)
    qualifying_record_count: int = Field(ge=0)
    reference_organisms: list[NonEmptyStr] = Field(default_factory=list)
    fusion_protein_ids: list[NonEmptyStr] = Field(default_factory=list)
    best_query_component_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    best_candidate_component_coverage: float | None = Field(default=None, ge=0.0, le=1.0)
    minimum_component_overlap_fraction: float | None = Field(default=None, ge=0.0, le=1.0)
    pair_supported: bool | None = None
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    calculation_rule_version: NonEmptyStr | None = None
    source: NonEmptyStr | None = None
    source_record_ids: list[NonEmptyStr] = Field(default_factory=list)


KnownInteractionType = Literal[
    "physical",
    "direct",
    "genetic",
    "functional_association",
    "co_complex",
    "co_expression",
    "predicted",
    "other",
]


IdentifierMappingStatus = Literal["mapped", "uncertain", "unmapped"]


class KnownInteractionObservation(BaseEvidence):
    protein_a_id: NonEmptyStr
    protein_b_id: NonEmptyStr
    interaction_type: KnownInteractionType
    reference_organism: NonEmptyStr
    detection_method: NonEmptyStr | None = None
    normalized_detection_method: NonEmptyStr | None = None
    publication_id: NonEmptyStr | None = None
    confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    is_direct: bool | None = None
    is_physical: bool | None = None
    is_biological: bool | None = None
    database_version: NonEmptyStr | None = None
    protein_a_reference_id: NonEmptyStr | None = None
    protein_b_reference_id: NonEmptyStr | None = None
    source: NonEmptyStr
    source_record_id: NonEmptyStr
    notes: NonEmptyStr | None = None
    identifier_mapping_status: IdentifierMappingStatus = "mapped"

    @model_validator(mode="after")
    def reject_self_interaction(self) -> "KnownInteractionObservation":
        if self.protein_a_id == self.protein_b_id:
            raise ValueError("interaction protein IDs must differ")
        return self


class KnownInteractionEvidence(BaseEvidence):
    query_protein_id: NonEmptyStr
    candidate_protein_id: NonEmptyStr
    supporting_record_count: int = Field(ge=0)
    qualifying_record_count: int = Field(ge=0)
    direct_record_count: int = Field(ge=0)
    physical_record_count: int = Field(ge=0)
    biological_record_count: int = Field(ge=0)
    independent_publication_count: int = Field(ge=0)
    independent_source_count: int = Field(ge=0)
    interaction_types: list[NonEmptyStr] = Field(default_factory=list)
    detection_methods: list[NonEmptyStr] = Field(default_factory=list)
    publication_ids: list[NonEmptyStr] = Field(default_factory=list)
    reference_organisms: list[NonEmptyStr] = Field(default_factory=list)
    sources: list[NonEmptyStr] = Field(default_factory=list)
    source_record_ids: list[NonEmptyStr] = Field(default_factory=list)
    best_confidence: float | None = Field(default=None, ge=0.0, le=1.0)
    pair_supported: bool | None = None
    direct_interaction_supported: bool | None = None
    physical_interaction_supported: bool | None = None
    functional_association_supported: bool | None = None
    support_terms: list[str] = Field(default_factory=list)
    conflicting_terms: list[str] = Field(default_factory=list)
    calculation_rule_version: NonEmptyStr | None = None


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
