"""Domain enumerations with stable serialized values."""

from enum import StrEnum


class PredictedRelationshipType(StrEnum):
    PHYSICAL_COMPLEX = "physical_complex"
    TRANSIENT_INTERACTION = "transient_interaction"
    ENZYME_SUBSTRATE = "enzyme_substrate"
    ACCESSORY_FACTOR = "accessory_factor"
    PATHWAY_ASSOCIATION = "pathway_association"
    GENE_CONTEXT_ONLY = "gene_context_only"
    FUNCTIONAL_SIMILARITY = "functional_similarity"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"


class EvidenceTier(StrEnum):
    TIER_A_STRONG_MULTI_EVIDENCE = "tier_a_strong_multi_evidence"
    TIER_B_MODERATE = "tier_b_moderate"
    TIER_C_EXPLORATORY = "tier_c_exploratory"
    TIER_D_CONTEXT_ONLY = "tier_d_context_only"
    TIER_E_INSUFFICIENT = "tier_e_insufficient"
    TIER_X_CONFLICTING = "tier_x_conflicting"


class CandidateDisposition(StrEnum):
    INCLUDED = "included"
    DOWN_RANKED = "down_ranked"
    FLAGGED = "flagged"
    EXCLUDED = "excluded"


class IdentifierMatchStatus(StrEnum):
    EXACT_MATCH = "exact_match"
    UNIQUE_ALIAS_MATCH = "unique_alias_match"
    AMBIGUOUS_MATCH = "ambiguous_match"
    NO_MATCH = "no_match"


class EvidenceOrigin(StrEnum):
    EXACT_PROTEIN = "exact_protein"
    EXACT_PAIR = "exact_pair"
    ORTHOLOG_TRANSFERRED = "ortholog_transferred"
    ANNOTATION = "annotation"
    LOCAL_PREDICTION = "local_prediction"
    INFERRED = "inferred"


class EvidenceStatus(StrEnum):
    AVAILABLE = "available"
    MISSING = "missing"
    NOT_APPLICABLE = "not_applicable"
    FAILED = "failed"
    NOT_RUN = "not_run"


class ContradictionSeverity(StrEnum):
    INFO = "info"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class ManualReviewStatus(StrEnum):
    NOT_REVIEWED = "not_reviewed"
    QUEUED = "queued"
    APPROVED = "approved"
    DEFERRED = "deferred"
    REJECTED = "rejected"
    COMPLETED = "completed"


class ManualStructurePriority(StrEnum):
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    URGENT_REVIEW = "urgent_review"


class RunStatus(StrEnum):
    INITIALIZED = "initialized"
    VALIDATED = "validated"
    COMPLETED = "completed"
    PARTIAL = "partial"
    FAILED = "failed"
