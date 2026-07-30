"""Deterministic candidate evidence orchestration and optional integrated scoring."""

from __future__ import annotations

import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import yaml
from pydantic import Field

from protein_interaction_hunter.adapters.local.annotation import LocalAnnotationTsvLoader
from protein_interaction_hunter.adapters.local.domain_rules import (
    LocalDomainRulesLoader,
)
from protein_interaction_hunter.adapters.local.domains import (
    LocalDomainTsvLoader,
)
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.functional_rules import (
    LocalFunctionalRulesLoader,
)
from protein_interaction_hunter.adapters.local.fusion import LocalFusionTsvLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.adapters.local.known_interactions import (
    LocalKnownInteractionsTsvLoader,
)
from protein_interaction_hunter.adapters.local.orthology import (
    LocalOrthologyTsvLoader,
)
from protein_interaction_hunter.adapters.local.phylogenetic_profile import (
    LocalPhylogeneticProfileTsvLoader,
)
from protein_interaction_hunter.application.candidates import generate_candidates
from protein_interaction_hunter.application.domain_pairs import (
    DOMAIN_PAIR_ENGINE_VERSION,
    build_domain_index,
    evaluate_domain_pairs,
)
from protein_interaction_hunter.application.evidence_tiers import (
    EVIDENCE_TIER_RULE_VERSION,
    HIGH_SPECIFICITY_DEFINITION_VERSION,
    evaluate_evidence_tier,
)
from protein_interaction_hunter.application.functional_complementarity import (
    FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION,
    evaluate_functional_complementarity,
)
from protein_interaction_hunter.application.fusion import (
    FUSION_ENGINE_VERSION,
    build_fusion_index,
    evaluate_fusion_pair,
)
from protein_interaction_hunter.application.gene_context import (
    GENE_CONTEXT_RULE_VERSION,
    build_coordinate_index,
    calculate_gene_context,
)
from protein_interaction_hunter.application.identifiers import NORMALIZATION_RULE_VERSION
from protein_interaction_hunter.application.known_interactions import (
    KNOWN_INTERACTIONS_ENGINE_VERSION,
    build_known_interaction_index,
    evaluate_known_interaction_pair,
)
from protein_interaction_hunter.application.localization import (
    LOCALIZATION_ENGINE_VERSION,
    evaluate_localization,
)
from protein_interaction_hunter.application.operon_proxy import (
    OPERON_PROXY_RULE_VERSION,
    calculate_operon_proxy,
)
from protein_interaction_hunter.application.orthology import (
    ORTHOLOGY_ENGINE_VERSION,
    build_orthology_index,
    evaluate_orthology_pair,
)
from protein_interaction_hunter.application.phylogenetic_profile import (
    PHYLOGENETIC_PROFILE_ENGINE_VERSION,
    build_phylogenetic_profile_index,
    evaluate_phylogenetic_profile_pair,
)
from protein_interaction_hunter.application.scoring import (
    SCORING_RULE_VERSION,
    candidate_score_from_integrated,
    rank_scores,
    score_pair,
)
from protein_interaction_hunter.config import AppConfig, load_config
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.manifest import build_input_file_manifest, build_run_manifest
from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.base import StrictModel
from protein_interaction_hunter.models.enums import (
    CandidateDisposition,
    ContextCompleteness,
    EvidenceStatus,
    PredictedRelationshipType,
    RunStatus,
)
from protein_interaction_hunter.models.evidence import (
    CandidateEvidenceBundle,
    DomainEvidence,
    EvidenceProvenance,
    EvidenceTierResult,
    FunctionalEvidence,
    FusionEvidence,
    GenomeContextEvidence,
    IntegratedScore,
    KnownInteractionEvidence,
    LocalizationEvidence,
    OperonEvidence,
    OrthologRecord,
    PhylogeneticProfileEvidence,
)
from protein_interaction_hunter.models.protein import CandidateProtein
from protein_interaction_hunter.models.scoring import CandidateScore
from protein_interaction_hunter.outputs.candidates import (
    CandidateTableTsvWriter,
    WarningSummaryTsvWriter,
)
from protein_interaction_hunter.outputs.excel import ExcelSchemaWriter
from protein_interaction_hunter.outputs.jsonl import (
    JsonlEvidenceBundleWriter,
    JsonRunManifestWriter,
)

_UNIMPLEMENTED_ENGINES: tuple[str, ...] = ()


class CandidateGenerationSummary(StrictModel):
    query_count: int = Field(ge=0)
    protein_count: int = Field(ge=0)
    pair_count: int = Field(ge=0)
    included_count: int = Field(ge=0)
    flagged_count: int = Field(ge=0)
    excluded_count: int = Field(ge=0)
    duplicate_group_count: int = Field(ge=0)
    fragment_candidate_count: int = Field(ge=0)
    hypothetical_protein_count: int = Field(ge=0)
    missing_coordinate_count: int = Field(ge=0)
    missing_annotation_count: int = Field(ge=0)
    ambiguous_mapping_count: int = Field(ge=0)
    same_contig_pair_count: int = Field(default=0, ge=0)
    different_contig_pair_count: int = Field(default=0, ge=0)
    overlapping_pair_count: int = Field(default=0, ge=0)
    missing_context_pair_count: int = Field(default=0, ge=0)
    ambiguous_context_pair_count: int = Field(default=0, ge=0)
    neighborhood_pair_count: int = Field(default=0, ge=0)
    incomplete_context_pair_count: int = Field(default=0, ge=0)
    output_path: Path


class PipelineResult(StrictModel):
    summary: CandidateGenerationSummary
    bundles: list[CandidateEvidenceBundle]
    evidence_path: Path
    candidate_table_path: Path
    manifest_path: Path
    config_snapshot_path: Path
    warning_summary_path: Path
    excel_path: Path | None = None


def _git_commit(repository: Path) -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def _functional_annotation_text(
    protein_id: str,
    annotation_by_protein: dict[str, AnnotationRecord],
    description_by_protein: dict[str, str],
) -> str | None:
    annotation = annotation_by_protein.get(protein_id)

    values = [
        annotation.product if annotation else None,
        annotation.functional_category if annotation else None,
        description_by_protein.get(protein_id),
    ]

    text = " | ".join(value.strip() for value in values if value is not None and value.strip())

    return text or None


def _policy_settings(config: AppConfig) -> dict[str, str | int | bool]:
    policy = config.candidate_generation
    return {
        "include_hypothetical_proteins": policy.include_hypothetical_proteins,
        "include_missing_coordinates": policy.include_missing_coordinates,
        "minimum_length_aa": policy.minimum_length_aa,
        "duplicate_sequence_policy": policy.duplicate_sequence_policy,
        "fragment_policy": policy.fragment_policy,
        "self_candidate_policy": policy.self_candidate_policy,
        "gene_context_enabled": config.gene_context.enabled,
        "neighborhood_gene_count": config.gene_context.neighborhood_gene_count,
        "require_query_coordinates": config.gene_context.require_query_coordinates,
        "operon_proxy_max_intergenic_bp": (config.gene_context.operon_proxy_max_intergenic_bp),
    }


def _write_snapshot(config: AppConfig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(config.model_dump(mode="json"), sort_keys=True, allow_unicode=True),
        encoding="utf-8",
        newline="\n",
    )
    return path


def _excel_context_rows(
    run_id: str,
    contexts: dict[tuple[str, str], GenomeContextEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), context in sorted(contexts.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Same_Contig": context.same_contig,
                "Query_Contig": context.query_contig,
                "Candidate_Contig": context.candidate_contig,
                "Query_Start": context.query_start,
                "Query_End": context.query_end,
                "Query_Strand": context.query_strand,
                "Candidate_Start": context.candidate_start,
                "Candidate_End": context.candidate_end,
                "Candidate_Strand": context.candidate_strand,
                "Strand_Relationship": context.strand_relationship,
                "Relative_Position": context.relative_position,
                "Coordinate_Position": context.coordinate_position,
                "Distance_BP": context.distance_bp,
                "Overlap_BP": context.overlap_bp,
                "Intervening_Gene_Count": context.intervening_gene_count,
                "Intervening_Feature_Count": context.intervening_feature_count,
                "Feature_Index_Delta": context.feature_index_delta,
                "Within_Neighborhood_Window": context.within_neighborhood_window,
                "Context_Completeness": context.context_completeness,
                "Status": context.status,
                "Warnings": "|".join(context.warnings),
                "Source": GENE_CONTEXT_RULE_VERSION,
                "Provenance": "|".join(
                    f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                    for item in context.provenance
                ),
            }
        )
    return rows


def _excel_operon_rows(
    run_id: str,
    operons: dict[tuple[str, str], OperonEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for (query_id, candidate_id), operon in sorted(operons.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Status": operon.status,
                "Proxy_Status": operon.proxy_status,
                "Same_Contig": operon.same_contig,
                "Same_Strand": operon.same_strand,
                "Is_Adjacent": operon.is_adjacent,
                "Intergenic_Distance_BP": operon.intergenic_distance_bp,
                "Overlap_BP": operon.overlap_bp,
                "Intervening_Gene_Count": operon.intervening_gene_count,
                "Transcriptional_Order": operon.transcriptional_order,
                "Maximum_Intergenic_Distance_BP": (operon.maximum_intergenic_distance_bp),
                "Passes_Distance_Threshold": operon.passes_distance_threshold,
                "Supporting_Conditions": "|".join(operon.supporting_conditions),
                "Conflicting_Conditions": "|".join(operon.conflicting_conditions),
                "Rule_Version": operon.calculation_rule_version,
                "Rule_ID": operon.proxy_rule_id,
                "Warnings": "|".join(operon.warnings),
                "Provenance": "|".join(
                    f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                    for item in operon.provenance
                ),
            }
        )

    return rows


def _excel_domain_rows(
    run_id: str,
    domain_evidence: dict[
        tuple[str, str],
        list[DomainEvidence],
    ],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for (query_id, candidate_id), evidence_records in sorted(domain_evidence.items()):
        for evidence in evidence_records:
            rows.append(
                {
                    "Run_ID": run_id,
                    "Query_ID": query_id,
                    "Candidate_ID": candidate_id,
                    "Status": evidence.status,
                    "Pair_Matched": evidence.pair_matched,
                    "Candidate_Protein_ID": evidence.protein_id,
                    "Candidate_Source": evidence.source,
                    "Candidate_Accession": evidence.accession,
                    "Candidate_Domain_Name": evidence.name,
                    "Candidate_Start": evidence.start,
                    "Candidate_End": evidence.end,
                    "Architecture_Index": evidence.architecture_index,
                    "Candidate_Role": evidence.role,
                    "Pair_Rule_ID": evidence.pair_rule_id,
                    "Query_Protein_ID": evidence.paired_protein_id,
                    "Query_Accession": evidence.paired_accession,
                    "Is_Shared": evidence.is_shared,
                    "Support_Terms": "|".join(evidence.support_terms),
                    "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                    "Rule_Version": (evidence.calculation_rule_version),
                    "Ruleset_Path": evidence.ruleset_path,
                    "Warnings": "|".join(evidence.warnings),
                    "Provenance": "|".join(
                        (f"{item.source_name}:{item.source_version or ''}:{item.method or ''}")
                        for item in evidence.provenance
                    ),
                }
            )

    return rows


def _excel_functional_rows(
    run_id: str,
    functional_evidence: dict[
        tuple[str, str],
        list[FunctionalEvidence],
    ],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for (query_id, candidate_id), evidence_records in sorted(functional_evidence.items()):
        for evidence in evidence_records:
            rows.append(
                {
                    "Run_ID": run_id,
                    "Query_ID": query_id,
                    "Candidate_ID": candidate_id,
                    "Status": evidence.status,
                    "Matched": evidence.matched,
                    "Query_Role": evidence.query_role,
                    "Candidate_Role": evidence.candidate_role,
                    "Relationship_Hint": evidence.relationship_hint,
                    "Rule_ID": evidence.rule_id,
                    "Query_Matched_Terms": "|".join(evidence.query_matched_terms),
                    "Candidate_Matched_Terms": "|".join(evidence.candidate_matched_terms),
                    "Support_Terms": "|".join(evidence.support_terms),
                    "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                    "Query_Annotation_Text": (evidence.query_annotation_text),
                    "Candidate_Annotation_Text": (evidence.candidate_annotation_text),
                    "Rule_Version": (evidence.calculation_rule_version),
                    "Ruleset_Path": evidence.ruleset_path,
                    "Warnings": "|".join(evidence.warnings),
                    "Provenance": "|".join(
                        (f"{item.source_name}:{item.source_version or ''}:{item.method or ''}")
                        for item in evidence.provenance
                    ),
                }
            )

    return rows


def _excel_localization_rows(
    run_id: str,
    localization_evidence: dict[
        tuple[str, str],
        LocalizationEvidence,
    ],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for (query_id, candidate_id), evidence in sorted(localization_evidence.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Status": evidence.status,
                "Candidate_Protein_ID": evidence.protein_id,
                "Query_Compartment": evidence.query_compartment,
                "Candidate_Compartment": evidence.candidate_compartment,
                "Compartment": evidence.compartment,
                "Compatibility": evidence.compatibility,
                "Signal_Peptide": evidence.signal_peptide,
                "Transmembrane_Helices": (evidence.transmembrane_helices),
                "Topology": evidence.topology,
                "Localization_Annotation": (evidence.localization_annotation),
                "Transmembrane_Annotation": (evidence.transmembrane_annotation),
                "Matched_Terms": "|".join(evidence.matched_terms),
                "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                "Rule_ID": evidence.rule_id,
                "Rule_Version": (evidence.calculation_rule_version),
                "Annotation_Source": evidence.annotation_source,
                "Annotation_Confidence": (evidence.annotation_confidence),
                "Warnings": "|".join(evidence.warnings),
                "Provenance": "|".join(
                    (f"{item.source_name}:{item.source_version or ''}:{item.method or ''}")
                    for item in evidence.provenance
                ),
            }
        )

    return rows


def _excel_orthology_rows(
    run_id: str,
    orthology_evidence: dict[
        tuple[str, str],
        list[OrthologRecord],
    ],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for (query_id, candidate_id), evidence_records in sorted(orthology_evidence.items()):
        for evidence in evidence_records:
            rows.append(
                {
                    "Run_ID": run_id,
                    "Query_ID": query_id,
                    "Candidate_ID": candidate_id,
                    "Status": evidence.status,
                    "Protein_ID": evidence.protein_id,
                    "Reference_ID": evidence.reference_id,
                    "Ortholog_ID": evidence.ortholog_id,
                    "Identity": evidence.identity,
                    "Query_Coverage": evidence.query_coverage,
                    "Subject_Coverage": evidence.subject_coverage,
                    "Evalue": evidence.evalue,
                    "Orthogroup": evidence.orthogroup,
                    "Paralog_Ambiguity": evidence.paralog_ambiguity,
                    "Reference_Organism": evidence.reference_organism,
                    "Relationship": evidence.relationship,
                    "Paired_Protein_ID": evidence.paired_protein_id,
                    "Paired_Reference_ID": evidence.paired_reference_id,
                    "Paired_Ortholog_ID": evidence.paired_ortholog_id,
                    "Paired_Orthogroup": evidence.paired_orthogroup,
                    "Shared_Orthogroup": evidence.shared_orthogroup,
                    "Pair_Supported": evidence.pair_supported,
                    "Support_Terms": "|".join(evidence.support_terms),
                    "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                    "Rule_Version": evidence.calculation_rule_version,
                    "Source": evidence.source,
                    "Source_Record_ID": evidence.source_record_id,
                    "Warnings": "|".join(evidence.warnings),
                    "Provenance": "|".join(
                        f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                        for item in evidence.provenance
                    ),
                }
            )

    return rows


def _excel_phylogenetic_profile_rows(
    run_id: str,
    profile_evidence: dict[tuple[str, str], PhylogeneticProfileEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), evidence in sorted(profile_evidence.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Status": evidence.status,
                "Informative_Species": evidence.informative_species_count,
                "Shared_Presence": evidence.shared_presence_count,
                "Shared_Absence": evidence.shared_absence_count,
                "Discordant": evidence.discordant_count,
                "Unknown": evidence.unknown_count,
                "Profile_Similarity": evidence.profile_similarity,
                "Pair_Supported": evidence.pair_supported,
                "Support_Terms": "|".join(evidence.support_terms),
                "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                "Rule_Version": evidence.calculation_rule_version,
                "Source": evidence.source,
                "Source_Record_ID": evidence.source_record_id,
                "Warnings": "|".join(evidence.warnings),
                "Provenance": "|".join(
                    "{}:{}:{}".format(
                        item.source_name,
                        item.source_version or "",
                        item.method or "",
                    )
                    for item in evidence.provenance
                ),
            }
        )
    return rows


def _excel_fusion_rows(
    run_id: str,
    fusion_evidence: dict[tuple[str, str], FusionEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), evidence in sorted(fusion_evidence.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Status": evidence.status,
                "Supporting_Record_Count": evidence.supporting_record_count,
                "Qualifying_Record_Count": evidence.qualifying_record_count,
                "Reference_Organisms": "|".join(evidence.reference_organisms),
                "Fusion_Protein_IDs": "|".join(evidence.fusion_protein_ids),
                "Best_Query_Component_Coverage": (evidence.best_query_component_coverage),
                "Best_Candidate_Component_Coverage": (evidence.best_candidate_component_coverage),
                "Minimum_Component_Overlap_Fraction": (evidence.minimum_component_overlap_fraction),
                "Pair_Supported": evidence.pair_supported,
                "Support_Terms": "|".join(evidence.support_terms),
                "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                "Rule_Version": evidence.calculation_rule_version,
                "Source": evidence.source,
                "Source_Record_IDs": "|".join(evidence.source_record_ids),
                "Warnings": "|".join(evidence.warnings),
                "Provenance": "|".join(
                    "{}:{}:{}".format(
                        item.source_name,
                        item.source_version or "",
                        item.method or "",
                    )
                    for item in evidence.provenance
                ),
            }
        )
    return rows


def _excel_known_interaction_rows(
    run_id: str,
    evidence_by_pair: dict[tuple[str, str], KnownInteractionEvidence],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), evidence in sorted(evidence_by_pair.items()):
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Status": evidence.status,
                "Supporting_Record_Count": evidence.supporting_record_count,
                "Qualifying_Record_Count": evidence.qualifying_record_count,
                "Direct_Record_Count": evidence.direct_record_count,
                "Physical_Record_Count": evidence.physical_record_count,
                "Biological_Record_Count": evidence.biological_record_count,
                "Independent_Publication_Count": evidence.independent_publication_count,
                "Independent_Source_Count": evidence.independent_source_count,
                "Interaction_Types": "|".join(evidence.interaction_types),
                "Detection_Methods": "|".join(evidence.detection_methods),
                "Publication_IDs": "|".join(evidence.publication_ids),
                "Reference_Organisms": "|".join(evidence.reference_organisms),
                "Sources": "|".join(evidence.sources),
                "Source_Record_IDs": "|".join(evidence.source_record_ids),
                "Best_Confidence": evidence.best_confidence,
                "Pair_Supported": evidence.pair_supported,
                "Direct_Interaction_Supported": evidence.direct_interaction_supported,
                "Physical_Interaction_Supported": evidence.physical_interaction_supported,
                "Functional_Association_Supported": evidence.functional_association_supported,
                "Support_Terms": "|".join(evidence.support_terms),
                "Conflicting_Terms": "|".join(evidence.conflicting_terms),
                "Rule_Version": evidence.calculation_rule_version,
                "Warnings": "|".join(evidence.warnings),
                "Provenance": "|".join(
                    "{}:{}:{}".format(
                        item.source_name,
                        item.source_version or "",
                        item.method or "",
                    )
                    for item in evidence.provenance
                ),
            }
        )
    return rows


def _excel_integrated_scoring_rows(
    run_id: str,
    scores: dict[tuple[str, str], IntegratedScore],
) -> list[dict[str, object]]:
    return [
        {
            "Run_ID": run_id,
            "Query_ID": query_id,
            "Candidate_ID": candidate_id,
            "Status": score.status,
            "Provisional_Score": score.provisional_score,
            "Integrated_Score": score.output_score,
            "Normalized_Score": score.normalized_score,
            "Rank": score.rank,
            "Tied_Rank": score.tied_rank,
            "Sufficient_Evidence": score.sufficient_evidence,
            "Available_Weight": score.available_weight,
            "Evidence_Category_Count": score.evidence_category_count,
            "Evidence_Component_Count": score.evidence_component_count,
            "Positive_Component_Count": score.positive_component_count,
            "Neutral_Component_Count": score.neutral_component_count,
            "Negative_Component_Count": score.negative_component_count,
            "Rule_Version": score.calculation_rule_version,
            "Support_Terms": "|".join(score.support_terms),
            "Conflicting_Terms": "|".join(score.conflicting_terms),
            "Warnings": "|".join(score.warnings),
            "Provenance": "|".join(
                f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                for item in score.provenance
            ),
        }
        for (query_id, candidate_id), score in sorted(scores.items())
    ]


def _excel_scoring_component_rows(
    run_id: str,
    scores: dict[tuple[str, str], IntegratedScore],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (query_id, candidate_id), score in sorted(scores.items()):
        for component in score.component_scores:
            rows.append(
                {
                    "Run_ID": run_id,
                    "Query_ID": query_id,
                    "Candidate_ID": candidate_id,
                    "Component_Name": component.component_name,
                    "Category_Name": component.category_name,
                    "Evidence_Status": component.evidence_status,
                    "Raw_Value": component.raw_value,
                    "Normalized_Value": component.normalized_value,
                    "Configured_Weight": component.configured_weight,
                    "Effective_Weight": component.effective_weight,
                    "Weighted_Contribution": component.weighted_contribution,
                    "Direction": component.direction,
                    "Applied": component.applied,
                    "Exclusion_Reason": component.exclusion_reason,
                    "Source_Rule_Version": component.source_rule_version,
                    "Support_Terms": "|".join(component.support_terms),
                    "Conflicting_Terms": "|".join(component.conflicting_terms),
                    "Warnings": "|".join(component.warnings),
                }
            )
    return rows


def _excel_evidence_tier_rows(
    run_id: str,
    tiers: dict[tuple[str, str], EvidenceTierResult],
) -> list[dict[str, object]]:
    return [
        {
            "Run_ID": run_id,
            "Query_ID": query_id,
            "Candidate_ID": candidate_id,
            "Status": result.status,
            "Assigned_Tier": result.assigned_tier,
            "Base_Tier": result.base_tier,
            "Tier_Eligible": result.tier_eligible,
            "Formal_Score": result.formal_score,
            "Rank": result.rank,
            "Sufficient_Evidence": result.sufficient_evidence,
            "Evidence_Category_Count": result.evidence_category_count,
            "Evidence_Component_Count": result.evidence_component_count,
            "Available_Weight": result.available_weight,
            "Positive_Component_Count": result.positive_component_count,
            "Neutral_Component_Count": result.neutral_component_count,
            "Negative_Component_Count": result.negative_component_count,
            "High_Specificity_Component_Count": result.high_specificity_component_count,
            "High_Specificity_Components": "|".join(result.high_specificity_components),
            "Direct_Interaction_Supported": result.direct_interaction_supported,
            "Physical_Interaction_Supported": result.physical_interaction_supported,
            "Functional_Association_Supported": result.functional_association_supported,
            "Fusion_Supported": result.fusion_supported,
            "Orthology_Supported": result.orthology_supported,
            "Phylogenetic_Profile_Supported": result.phylogenetic_profile_supported,
            "Explicit_Conflict_Present": result.explicit_conflict_present,
            "Predicted_Only": result.predicted_only,
            "Functional_Association_Only": result.functional_association_only,
            "Applied_Tier_Caps": "|".join(result.applied_tier_caps),
            "Satisfied_Requirements": "|".join(result.satisfied_requirements),
            "Failed_Requirements": "|".join(result.failed_requirements),
            "Rule_Version": result.calculation_rule_version,
            "Support_Terms": "|".join(result.support_terms),
            "Conflicting_Terms": "|".join(result.conflicting_terms),
            "Warnings": "|".join(result.warnings),
            "Provenance": "|".join(
                f"{item.source_name}:{item.source_version or ''}:{item.method or ''}"
                for item in result.provenance
            ),
        }
        for (query_id, candidate_id), result in sorted(tiers.items())
    ]


def _excel_tier_summary_rows(
    run_id: str,
    candidates: list[CandidateProtein],
    tiers: dict[tuple[str, str], EvidenceTierResult],
) -> list[dict[str, object]]:
    if not tiers:
        return []
    rows: list[dict[str, object]] = []
    for query_id in sorted({candidate.query_id for candidate in candidates}):
        query_candidates = [item for item in candidates if item.query_id == query_id]
        results = [
            tiers[(item.query_id, item.protein_id)]
            for item in query_candidates
            if (item.query_id, item.protein_id) in tiers
        ]
        counts = {
            name: sum(result.assigned_tier.value == name for result in results)
            for name in ("tier_1", "tier_2", "tier_3", "tier_4", "unclassified")
        }
        rows.append(
            {
                "Run_ID": run_id,
                "Query_ID": query_id,
                "Tier_1_Count": counts["tier_1"],
                "Tier_2_Count": counts["tier_2"],
                "Tier_3_Count": counts["tier_3"],
                "Tier_4_Count": counts["tier_4"],
                "Unclassified_Count": counts["unclassified"],
                "Eligible_Count": sum(result.tier_eligible for result in results),
                "Total_Candidate_Count": len(query_candidates),
            }
        )
    return rows


def _excel_candidate_ranking_rows(
    candidates: list[CandidateProtein],
    scores: dict[tuple[str, str], IntegratedScore],
    tiers: dict[tuple[str, str], EvidenceTierResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        query_id = candidate.query_id
        candidate_id = candidate.protein_id
        pair = (query_id, candidate_id)
        score = scores.get(pair)
        tier = tiers.get(pair)
        rows.append(
            {
                "Rank": score.rank if score else None,
                "Query_ID": query_id,
                "Candidate_ID": candidate_id,
                "Disposition": candidate.disposition,
                "Relationship_Type": PredictedRelationshipType.INSUFFICIENT_EVIDENCE,
                "Evidence_Tier": tier.assigned_tier if tier else None,
                "Total_Ranking_Score": score.output_score if score else None,
            }
        )
    return rows


class InteractionCandidatePipeline:
    """Generate auditable candidates, evidence, and optional integrated scores."""

    def run(
        self, config_path: Path | None = None, command_line: list[str] | None = None
    ) -> PipelineResult:
        if config_path is None:
            raise ValueError("A config path is required to run the candidate pipeline.")
        resolved_config_path = config_path.expanduser().resolve()
        config = load_config(resolved_config_path)
        proteins = LocalFastaLoader().load(config.input.proteome_fasta)
        document = LocalGff3Loader().load_document(config.input.genome_gff)
        coordinates = document.features
        annotations: list[AnnotationRecord] = []
        if config.input.annotation_table is not None:
            annotations = LocalAnnotationTsvLoader().load(config.input.annotation_table)
        annotation_by_protein = {annotation.protein_id: annotation for annotation in annotations}
        generated = generate_candidates(
            proteins=proteins,
            coordinates=coordinates,
            annotations=annotations,
            query_ids=config.query.protein_ids,
            policy=config.candidate_generation,
        )
        canonical_query_ids = {query.query_id: query.protein_id for query in generated.queries}
        contexts: dict[tuple[str, str], GenomeContextEvidence] = {}
        if config.gene_context.enabled:
            coordinate_index = build_coordinate_index(proteins, document)
            if config.gene_context.require_query_coordinates:
                for query in generated.queries:
                    if query.protein_id in coordinate_index.ambiguous_proteins:
                        raise InputValidationError(
                            f"Query coordinate is ambiguous: {query.protein_id}"
                        )
                    if query.protein_id not in coordinate_index.by_protein:
                        raise InputValidationError(
                            f"Query coordinate is missing: {query.protein_id}"
                        )
            for candidate in generated.candidates:
                contexts[(candidate.query_id, candidate.protein_id)] = calculate_gene_context(
                    canonical_query_ids[candidate.query_id],
                    candidate.protein_id,
                    coordinate_index,
                    config.gene_context.neighborhood_gene_count,
                )
        operons: dict[tuple[str, str], OperonEvidence] = {}
        if config.gene_context.enabled:
            for pair, pair_context in contexts.items():
                operons[pair] = calculate_operon_proxy(
                    pair_context,
                    config.gene_context.operon_proxy_max_intergenic_bp,
                )

        orthology_evidence: dict[
            tuple[str, str],
            list[OrthologRecord],
        ] = {}

        if config.orthology.enabled:
            if config.orthology.source != "local_table":
                raise InputValidationError(
                    "orthology.source must be 'local_table' "
                    "when orthology.enabled is true in MVP-1G"
                )

            orthology_table_path = config.orthology.local_table

            if orthology_table_path is None:
                raise InputValidationError(
                    "orthology.local_table is required when orthology.enabled is true"
                )

            orthology_records = LocalOrthologyTsvLoader().load(orthology_table_path)
            orthology_index = build_orthology_index(orthology_records)

            for candidate in generated.candidates:
                pair = (
                    candidate.query_id,
                    candidate.protein_id,
                )
                query_protein_id = canonical_query_ids[candidate.query_id]

                orthology_evidence[pair] = evaluate_orthology_pair(
                    query_protein_id,
                    candidate.protein_id,
                    orthology_index,
                )

        phylogenetic_profile_evidence: dict[
            tuple[str, str],
            PhylogeneticProfileEvidence,
        ] = {}

        if config.phylogenetic_profile.enabled:
            profile_table_path = config.phylogenetic_profile.local_table
            if profile_table_path is None:
                raise InputValidationError(
                    "phylogenetic_profile.local_table is required when "
                    "phylogenetic_profile.enabled is true"
                )
            profile_records = LocalPhylogeneticProfileTsvLoader().load(profile_table_path)
            profile_index = build_phylogenetic_profile_index(profile_records)
            for candidate in generated.candidates:
                pair = (candidate.query_id, candidate.protein_id)
                phylogenetic_profile_evidence[pair] = evaluate_phylogenetic_profile_pair(
                    canonical_query_ids[candidate.query_id],
                    candidate.protein_id,
                    profile_index,
                    minimum_shared_species=(config.phylogenetic_profile.minimum_shared_species),
                    minimum_informative_species=(
                        config.phylogenetic_profile.minimum_informative_species
                    ),
                    minimum_profile_similarity=(
                        config.phylogenetic_profile.minimum_profile_similarity
                    ),
                )

        fusion_evidence: dict[tuple[str, str], FusionEvidence] = {}
        if config.fusion.enabled:
            fusion_table_path = config.fusion.local_table
            if fusion_table_path is None:
                raise InputValidationError(
                    "fusion.local_table is required when fusion.enabled is true"
                )
            fusion_records = LocalFusionTsvLoader().load(fusion_table_path)
            fusion_index = build_fusion_index(fusion_records)
            for candidate in generated.candidates:
                pair = (candidate.query_id, candidate.protein_id)
                fusion_evidence[pair] = evaluate_fusion_pair(
                    canonical_query_ids[candidate.query_id],
                    candidate.protein_id,
                    fusion_index,
                    minimum_supporting_records=(config.fusion.minimum_supporting_records),
                    minimum_component_coverage=(config.fusion.minimum_component_coverage),
                    maximum_component_overlap_fraction=(
                        config.fusion.maximum_component_overlap_fraction
                    ),
                )

        known_interaction_evidence: dict[tuple[str, str], KnownInteractionEvidence] = {}
        if config.known_interactions.enabled:
            interaction_table_path = config.known_interactions.local_table
            if interaction_table_path is None:
                raise InputValidationError(
                    "known_interactions.local_table is required when "
                    "known_interactions.enabled is true"
                )
            interaction_records = LocalKnownInteractionsTsvLoader().load(interaction_table_path)
            interaction_index = build_known_interaction_index(interaction_records)
            for candidate in generated.candidates:
                pair = (candidate.query_id, candidate.protein_id)
                known_interaction_evidence[pair] = evaluate_known_interaction_pair(
                    canonical_query_ids[candidate.query_id],
                    candidate.protein_id,
                    interaction_index,
                    minimum_supporting_records=(
                        config.known_interactions.minimum_supporting_records
                    ),
                    minimum_direct_records=(config.known_interactions.minimum_direct_records),
                    accepted_interaction_types=(
                        config.known_interactions.accepted_interaction_types
                    ),
                    accepted_evidence_methods=(config.known_interactions.accepted_evidence_methods),
                    excluded_evidence_methods=(config.known_interactions.excluded_evidence_methods),
                    minimum_confidence=config.known_interactions.minimum_confidence,
                )

        domain_evidence: dict[
            tuple[str, str],
            list[DomainEvidence],
        ] = {}

        if config.domains.enabled:
            if config.domains.source != "local_table":
                raise InputValidationError(
                    "domains.source must be 'local_table' when domains.enabled is true"
                )

            domain_table_path = config.domains.local_table
            domain_rules_path = config.domains.rules_path

            if domain_table_path is None:
                raise InputValidationError(
                    "domains.local_table is required when domains.enabled is true"
                )

            if domain_rules_path is None:
                raise InputValidationError(
                    "domains.rules_path is required when domains.enabled is true"
                )

            domain_records = LocalDomainTsvLoader().load(domain_table_path)
            domain_rules = LocalDomainRulesLoader().load(domain_rules_path)
            domain_index = build_domain_index(domain_records)

            for candidate in generated.candidates:
                pair = (
                    candidate.query_id,
                    candidate.protein_id,
                )
                query_protein_id = canonical_query_ids[candidate.query_id]

                domain_evidence[pair] = evaluate_domain_pairs(
                    query_protein_id,
                    candidate.protein_id,
                    domain_index,
                    domain_rules,
                    domain_rules_path,
                )

        functional_evidence: dict[
            tuple[str, str],
            list[FunctionalEvidence],
        ] = {}

        if config.functional_complementarity.enabled:
            rules_path = config.functional_complementarity.rules_path
            if rules_path is None:
                raise InputValidationError(
                    "functional_complementarity.rules_path is required "
                    "when functional_complementarity.enabled is true"
                )

            functional_rules = LocalFunctionalRulesLoader().load(rules_path)
            description_by_protein = {
                protein.protein_id: protein.description for protein in proteins
            }

            for candidate in generated.candidates:
                pair = (
                    candidate.query_id,
                    candidate.protein_id,
                )
                query_protein_id = canonical_query_ids[candidate.query_id]

                query_text = _functional_annotation_text(
                    query_protein_id,
                    annotation_by_protein,
                    description_by_protein,
                )
                candidate_text = _functional_annotation_text(
                    candidate.protein_id,
                    annotation_by_protein,
                    description_by_protein,
                )

                functional_evidence[pair] = evaluate_functional_complementarity(
                    query_text,
                    candidate_text,
                    functional_rules,
                    rules_path,
                )

        localization_evidence: dict[
            tuple[str, str],
            LocalizationEvidence,
        ] = {}

        if config.localization.enabled:
            if config.localization.source != "annotation_only":
                raise InputValidationError(
                    "localization.source must be 'annotation_only' "
                    "when localization.enabled is true"
                )

            for candidate in generated.candidates:
                pair = (
                    candidate.query_id,
                    candidate.protein_id,
                )
                query_protein_id = canonical_query_ids[candidate.query_id]

                localization_evidence[pair] = evaluate_localization(
                    query_protein_id,
                    candidate.protein_id,
                    annotation_by_protein,
                )

        integrated_scores: dict[tuple[str, str], IntegratedScore] = {}
        if config.scoring.enabled:
            for candidate in generated.candidates:
                pair = (candidate.query_id, candidate.protein_id)
                integrated_scores[pair] = score_pair(
                    canonical_query_ids[candidate.query_id],
                    candidate,
                    config.scoring,
                    genome_context=contexts.get(pair),
                    operon=operons.get(pair),
                    domains=domain_evidence.get(pair, []),
                    functional=functional_evidence.get(pair, []),
                    localization=localization_evidence.get(pair),
                    orthology=orthology_evidence.get(pair, []),
                    phylogenetic_profile=phylogenetic_profile_evidence.get(pair),
                    fusion=fusion_evidence.get(pair),
                    known_interactions=known_interaction_evidence.get(pair),
                )
            integrated_scores = rank_scores(
                integrated_scores,
                generated.candidates,
                config.scoring.tie_precision,
            )

        evidence_tiers: dict[tuple[str, str], EvidenceTierResult] = {}
        if config.evidence_tiers.enabled:
            for candidate in generated.candidates:
                pair = (candidate.query_id, candidate.protein_id)
                evidence_tiers[pair] = evaluate_evidence_tier(
                    canonical_query_ids[candidate.query_id],
                    candidate,
                    integrated_scores.get(pair),
                    config.evidence_tiers,
                    known_interactions=known_interaction_evidence.get(pair),
                    fusion=fusion_evidence.get(pair),
                    orthology=orthology_evidence.get(pair, []),
                    phylogenetic_profile=phylogenetic_profile_evidence.get(pair),
                )

        config_hash = build_input_file_manifest(
            "config", resolved_config_path, required=True
        ).sha256
        assert config_hash is not None
        run_id = f"{config.project.run_name}-{config_hash[:12]}"
        output_path = config.output.directory
        output_path.mkdir(parents=True, exist_ok=True)
        snapshot_path = _write_snapshot(config, output_path / "config.snapshot.yaml")
        policies = _policy_settings(config)
        provenance = EvidenceProvenance(
            source_name="candidate_generation",
            source_version=NORMALIZATION_RULE_VERSION,
            method="deterministic_all_proteome_enumeration",
            metadata={"policy_settings": policies},
        )
        bundles: list[CandidateEvidenceBundle] = []
        for candidate in generated.candidates:
            pair = (candidate.query_id, candidate.protein_id)
            context = contexts.get(pair)
            operon = operons.get(pair)
            domains = domain_evidence.get(pair, [])
            functional = functional_evidence.get(pair, [])
            localization = localization_evidence.get(pair)
            orthology = orthology_evidence.get(pair, [])
            profile = phylogenetic_profile_evidence.get(pair)
            fusion = fusion_evidence.get(pair)
            known_interaction = known_interaction_evidence.get(pair)
            integrated_score = integrated_scores.get(pair)
            tier_result = evidence_tiers.get(pair)

            statuses = {engine: EvidenceStatus.NOT_RUN for engine in _UNIMPLEMENTED_ENGINES}
            statuses["gene_context"] = context.status if context else EvidenceStatus.NOT_RUN
            statuses["operon"] = operon.status if operon else EvidenceStatus.NOT_RUN
            statuses["domains"] = domains[0].status if domains else EvidenceStatus.NOT_RUN
            statuses["functional_complementarity"] = (
                functional[0].status if functional else EvidenceStatus.NOT_RUN
            )
            statuses["localization"] = (
                localization.status if localization else EvidenceStatus.NOT_RUN
            )
            statuses["orthology"] = orthology[0].status if orthology else EvidenceStatus.NOT_RUN
            statuses["phylogenetic_profile"] = profile.status if profile else EvidenceStatus.NOT_RUN
            statuses["fusion"] = fusion.status if fusion else EvidenceStatus.NOT_RUN
            statuses["known_interactions"] = (
                known_interaction.status if known_interaction else EvidenceStatus.NOT_RUN
            )
            statuses["scoring"] = (
                integrated_score.status if integrated_score else EvidenceStatus.NOT_RUN
            )
            statuses["evidence_tiers"] = (
                tier_result.status if tier_result else EvidenceStatus.NOT_RUN
            )

            context_warnings = context.warnings if context else []
            operon_warnings = operon.warnings if operon else []
            domain_warnings = [warning for evidence in domains for warning in evidence.warnings]
            functional_warnings = [
                warning for evidence in functional for warning in evidence.warnings
            ]
            localization_warnings = localization.warnings if localization else []
            orthology_warnings = [
                warning for evidence in orthology for warning in evidence.warnings
            ]
            profile_warnings = profile.warnings if profile else []
            fusion_warnings = fusion.warnings if fusion else []
            known_interaction_warnings = known_interaction.warnings if known_interaction else []
            scoring_warnings = integrated_score.warnings if integrated_score else []
            tier_warnings = tier_result.warnings if tier_result else []
            bundles.append(
                CandidateEvidenceBundle(
                    domains=domains,
                    functional=functional,
                    localization=([localization] if localization else []),
                    orthology=orthology,
                    phylogenetic_profile=[profile] if profile else [],
                    fusion=[fusion] if fusion else [],
                    known_interactions=([known_interaction] if known_interaction else []),
                    integrated_scoring=([integrated_score] if integrated_score else []),
                    evidence_tiers=([tier_result] if tier_result else []),
                    evidence_tier=(tier_result.assigned_tier if tier_result else None),
                    score=(
                        candidate_score_from_integrated(integrated_score)
                        if integrated_score
                        else CandidateScore()
                    ),
                    run_id=run_id,
                    query_id=candidate.query_id,
                    candidate_id=candidate.protein_id,
                    candidate=candidate,
                    candidate_disposition=candidate.disposition,
                    predicted_relationship_type=(PredictedRelationshipType.INSUFFICIENT_EVIDENCE),
                    genome_context=[context] if context else [],
                    operon=[operon] if operon else [],
                    engine_statuses=statuses,
                    provenance=[provenance],
                    policy_settings=policies,
                    warnings=sorted(
                        set(
                            candidate.warnings
                            + context_warnings
                            + operon_warnings
                            + domain_warnings
                            + functional_warnings
                            + localization_warnings
                            + orthology_warnings
                            + profile_warnings
                            + fusion_warnings
                            + known_interaction_warnings
                            + scoring_warnings
                            + tier_warnings
                        )
                    ),
                )
            )
        evidence_path = JsonlEvidenceBundleWriter().write(
            bundles, output_path / "candidate_evidence_bundle.jsonl"
        )
        candidate_table_path = CandidateTableTsvWriter().write(
            run_id,
            generated.candidates,
            output_path / "candidate_table.tsv",
            contexts=contexts,
            operons=operons,
            domains=domain_evidence,
            localization=localization_evidence,
            functional=functional_evidence,
            orthology=orthology_evidence,
            phylogenetic_profile=phylogenetic_profile_evidence,
            fusion=fusion_evidence,
            known_interactions=known_interaction_evidence,
            integrated_scores=integrated_scores,
            evidence_tiers=evidence_tiers,
        )
        all_warnings = [warning for bundle in bundles for warning in bundle.warnings]
        warning_summary_path = WarningSummaryTsvWriter().write(
            all_warnings, output_path / "warning_summary.tsv"
        )
        excel_path = None
        if config.output.write_excel:
            excel_path = ExcelSchemaWriter().write(
                output_path / "ProteinInteractionHunter.xlsx",
                rows_by_sheet={
                    "Gene_Context": _excel_context_rows(run_id, contexts),
                    "Operon_Proxy": _excel_operon_rows(run_id, operons),
                    "Domain_Complementarity": _excel_domain_rows(
                        run_id,
                        domain_evidence,
                    ),
                    "Functional_Complementarity": _excel_functional_rows(
                        run_id,
                        functional_evidence,
                    ),
                    "Localization_Evidence": _excel_localization_rows(
                        run_id,
                        localization_evidence,
                    ),
                    "Orthology_Evidence": _excel_orthology_rows(
                        run_id,
                        orthology_evidence,
                    ),
                    "Phylogenetic_Profile_Evidence": (
                        _excel_phylogenetic_profile_rows(
                            run_id,
                            phylogenetic_profile_evidence,
                        )
                    ),
                    "Fusion_Evidence": _excel_fusion_rows(run_id, fusion_evidence),
                    "Known_Interactions_Evidence": (
                        _excel_known_interaction_rows(run_id, known_interaction_evidence)
                    ),
                    "Integrated_Scoring": _excel_integrated_scoring_rows(run_id, integrated_scores),
                    "Scoring_Components": _excel_scoring_component_rows(run_id, integrated_scores),
                    "Evidence_Tiers": _excel_evidence_tier_rows(run_id, evidence_tiers),
                    "Tier_Summary": _excel_tier_summary_rows(
                        run_id, generated.candidates, evidence_tiers
                    ),
                    "Candidate_Ranking": _excel_candidate_ranking_rows(
                        generated.candidates, integrated_scores, evidence_tiers
                    ),
                },
            )
        input_files = [
            build_input_file_manifest("proteome_fasta", config.input.proteome_fasta, required=True),
            build_input_file_manifest("genome_gff", config.input.genome_gff, required=True),
        ]
        if config.input.annotation_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "annotation_table", config.input.annotation_table, required=False
                )
            )
        if config.domains.local_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "domain_local_table",
                    config.domains.local_table,
                    required=config.domains.enabled,
                )
            )
        if config.domains.rules_path is not None:
            input_files.append(
                build_input_file_manifest(
                    "domain_pair_rules",
                    config.domains.rules_path,
                    required=config.domains.enabled,
                )
            )
        if config.orthology.local_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "orthology_local_table",
                    config.orthology.local_table,
                    required=config.orthology.enabled,
                )
            )
        if config.phylogenetic_profile.local_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "phylogenetic_profile_local_table",
                    config.phylogenetic_profile.local_table,
                    required=config.phylogenetic_profile.enabled,
                )
            )
        if config.fusion.enabled and config.fusion.local_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "fusion_local_table",
                    config.fusion.local_table,
                    required=True,
                )
            )
        if config.known_interactions.enabled and config.known_interactions.local_table is not None:
            input_files.append(
                build_input_file_manifest(
                    "known_interactions_local_table",
                    config.known_interactions.local_table,
                    required=True,
                )
            )
        manifest = build_run_manifest(
            run_id=run_id,
            run_name=config.project.run_name,
            config_path=resolved_config_path,
            config_snapshot_path=snapshot_path,
            input_files=input_files,
            random_seed=config.project.random_seed,
            command_line=command_line or sys.argv,
        )
        manifest.status = RunStatus.COMPLETED
        manifest.completed_at = datetime.now(UTC)
        manifest.git_commit = _git_commit(Path(__file__).resolve().parents[3])
        manifest.normalization_rule_version = NORMALIZATION_RULE_VERSION
        manifest.gene_context_rule_version = (
            GENE_CONTEXT_RULE_VERSION if config.gene_context.enabled else None
        )
        manifest.orthology_rule_version = (
            ORTHOLOGY_ENGINE_VERSION if config.orthology.enabled else None
        )
        manifest.phylogenetic_profile_rule_version = (
            PHYLOGENETIC_PROFILE_ENGINE_VERSION if config.phylogenetic_profile.enabled else None
        )
        manifest.fusion_rule_version = FUSION_ENGINE_VERSION if config.fusion.enabled else None
        manifest.known_interactions_rule_version = (
            KNOWN_INTERACTIONS_ENGINE_VERSION if config.known_interactions.enabled else None
        )
        manifest.scoring_rule_version = SCORING_RULE_VERSION if config.scoring.enabled else None
        manifest.scoring_config_snapshot = (
            config.scoring.model_dump(mode="json") if config.scoring.enabled else {}
        )
        manifest.evidence_tier_rule_version = (
            EVIDENCE_TIER_RULE_VERSION if config.evidence_tiers.enabled else None
        )
        manifest.evidence_tier_config_snapshot = (
            config.evidence_tiers.model_dump(mode="json") if config.evidence_tiers.enabled else {}
        )
        manifest.evidence_tier_thresholds = (
            {
                name: getattr(config.evidence_tiers, name).model_dump(mode="json")
                for name in ("tier_1", "tier_2", "tier_3", "tier_4")
            }
            if config.evidence_tiers.enabled
            else {}
        )
        manifest.evidence_tier_caps = (
            {
                "explicit_conflict": config.evidence_tiers.explicit_conflict_tier_cap,
                "predicted_only": config.evidence_tiers.predicted_only_tier_cap,
                "functional_association_only": (
                    config.evidence_tiers.functional_association_only_tier_cap
                ),
                "insufficient_evidence": config.evidence_tiers.insufficient_evidence_tier,
            }
            if config.evidence_tiers.enabled
            else {}
        )
        manifest.evidence_tier_high_specificity_definition_version = (
            HIGH_SPECIFICITY_DEFINITION_VERSION if config.evidence_tiers.enabled else None
        )
        manifest.evidence_tier_scoring_rule_version_dependency = (
            SCORING_RULE_VERSION if config.evidence_tiers.enabled else None
        )

        if config.gene_context.enabled:
            manifest.policy_settings["operon_proxy_rule_version"] = OPERON_PROXY_RULE_VERSION

        if config.functional_complementarity.enabled:
            manifest.policy_settings["functional_complementarity_rule_version"] = (
                FUNCTIONAL_COMPLEMENTARITY_ENGINE_VERSION
            )

        if config.domains.enabled:
            manifest.policy_settings["domain_pair_rule_version"] = DOMAIN_PAIR_ENGINE_VERSION

        if config.localization.enabled:
            manifest.policy_settings["localization_rule_version"] = LOCALIZATION_ENGINE_VERSION

        manifest.policy_settings = policies | manifest.policy_settings
        manifest.warnings = sorted(set(all_warnings))
        manifest.parser_warnings = sorted(
            set(document.warnings)
            | {warning for protein in proteins for warning in protein.warnings}
            | {warning for coordinate in coordinates for warning in coordinate.warnings}
            | {warning for annotation in annotations for warning in annotation.warnings}
        )
        manifest.incomplete_evidence_flags = [
            f"{engine}_not_run" for engine in _UNIMPLEMENTED_ENGINES
        ]
        if not config.gene_context.enabled:
            manifest.incomplete_evidence_flags.extend(["gene_context_not_run", "operon_not_run"])
        if not config.phylogenetic_profile.enabled:
            manifest.incomplete_evidence_flags.append("phylogenetic_profile_not_run")
        if not config.fusion.enabled:
            manifest.incomplete_evidence_flags.append("fusion_not_run")
        if not config.known_interactions.enabled:
            manifest.incomplete_evidence_flags.append("known_interactions_not_run")
        if not config.scoring.enabled:
            manifest.incomplete_evidence_flags.append("scoring_not_run")
        if not config.evidence_tiers.enabled:
            manifest.incomplete_evidence_flags.append("evidence_tiers_not_run")
        manifest_path = JsonRunManifestWriter().write(manifest, output_path / "run_manifest.json")
        disposition_counts = {
            disposition: sum(
                candidate.disposition is disposition for candidate in generated.candidates
            )
            for disposition in CandidateDisposition
        }
        context_values = list(contexts.values())
        summary = CandidateGenerationSummary(
            query_count=len(generated.queries),
            protein_count=len(proteins),
            pair_count=len(generated.candidates),
            included_count=disposition_counts[CandidateDisposition.INCLUDED],
            flagged_count=disposition_counts[CandidateDisposition.FLAGGED],
            excluded_count=disposition_counts[CandidateDisposition.EXCLUDED],
            duplicate_group_count=len(generated.duplicate_groups),
            fragment_candidate_count=sum(
                candidate.is_fragment_candidate for candidate in generated.candidates
            ),
            hypothetical_protein_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if candidate.is_hypothetical
                }
            ),
            missing_coordinate_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if not candidate.has_coordinate
                }
            ),
            missing_annotation_count=len(
                {
                    candidate.protein_id
                    for candidate in generated.candidates
                    if not candidate.has_annotation
                }
            ),
            ambiguous_mapping_count=generated.ambiguous_mapping_count,
            same_contig_pair_count=sum(item.same_contig is True for item in context_values),
            different_contig_pair_count=sum(
                item.status is EvidenceStatus.NOT_APPLICABLE for item in context_values
            ),
            overlapping_pair_count=sum(
                item.status is EvidenceStatus.AVAILABLE and bool(item.overlap_bp)
                for item in context_values
            ),
            missing_context_pair_count=sum(
                item.status is EvidenceStatus.MISSING for item in context_values
            ),
            ambiguous_context_pair_count=sum(
                item.status is EvidenceStatus.FAILED for item in context_values
            ),
            neighborhood_pair_count=sum(
                item.within_neighborhood_window is True for item in context_values
            ),
            incomplete_context_pair_count=sum(
                item.status is EvidenceStatus.AVAILABLE
                and item.context_completeness is not ContextCompleteness.COMPLETE
                for item in context_values
            ),
            output_path=output_path.resolve(),
        )
        return PipelineResult(
            summary=summary,
            bundles=bundles,
            evidence_path=evidence_path,
            candidate_table_path=candidate_table_path,
            manifest_path=manifest_path,
            config_snapshot_path=snapshot_path,
            warning_summary_path=warning_summary_path,
            excel_path=excel_path,
        )
