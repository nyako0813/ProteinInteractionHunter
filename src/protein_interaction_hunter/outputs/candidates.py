"""Deterministic flat candidate and warning summary writers."""

import csv
import io
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from protein_interaction_hunter.models.evidence import (
    DomainEvidence,
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

CANDIDATE_COLUMNS = (
    "run_id",
    "query_id",
    "candidate_id",
    "candidate_description",
    "candidate_disposition",
    "disposition_reasons",
    "sequence_length",
    "gene_id",
    "locus_tag",
    "old_locus_tag",
    "contig",
    "strand",
    "has_coordinate",
    "has_annotation",
    "same_contig_as_query",
    "is_duplicate_sequence",
    "duplicate_group_id",
    "is_fragment_candidate",
    "fragment_reasons",
    "is_hypothetical",
    "identifier_match_status",
    "warnings",
    "same_contig",
    "query_start",
    "query_end",
    "query_strand",
    "candidate_start",
    "candidate_end",
    "candidate_strand",
    "strand_relationship",
    "relative_position",
    "coordinate_position",
    "distance_bp",
    "overlap_bp",
    "intervening_gene_count",
    "intervening_feature_count",
    "feature_index_delta",
    "within_neighborhood_window",
    "context_completeness",
    "gene_context_status",
    "edge_to_edge_distance_bp",
    "within_neighborhood_gene_count",
    "query_distance_to_contig_left_edge",
    "query_distance_to_contig_right_edge",
    "candidate_distance_to_contig_left_edge",
    "candidate_distance_to_contig_right_edge",
    "gene_context_rule_version",
    "gene_context_warnings",
    "operon_status",
    "operon_proxy_status",
    "operon_same_contig",
    "operon_same_strand",
    "operon_is_adjacent",
    "operon_intergenic_distance_bp",
    "operon_overlap_bp",
    "operon_intervening_gene_count",
    "operon_transcriptional_order",
    "operon_maximum_intergenic_distance_bp",
    "operon_passes_distance_threshold",
    "operon_supporting_conditions",
    "operon_conflicting_conditions",
    "operon_rule_version",
    "operon_rule_id",
    "operon_warnings",
    "localization_status",
    "localization_compartment",
    "localization_query_compartment",
    "localization_compatibility",
    "localization_signal_peptide",
    "localization_tm_helices",
    "localization_topology",
    "localization_matched_terms",
    "localization_conflicting_terms",
    "localization_rule_version",
    "localization_annotation_source",
    "localization_annotation_confidence",
    "localization_warnings",
    "functional_status",
    "functional_matched",
    "functional_relationship_hints",
    "functional_rule_ids",
    "functional_query_roles",
    "functional_candidate_roles",
    "functional_query_matched_terms",
    "functional_candidate_matched_terms",
    "functional_support_terms",
    "functional_conflicting_terms",
    "functional_rule_versions",
    "functional_warnings",
    "domain_status",
    "domain_pair_matched",
    "domain_pair_rule_ids",
    "domain_candidate_accessions",
    "domain_query_accessions",
    "domain_candidate_roles",
    "domain_is_shared",
    "domain_support_terms",
    "domain_conflicting_terms",
    "domain_rule_versions",
    "domain_warnings",
    "orthology_status",
    "orthology_reference_organism",
    "orthology_relationship",
    "orthology_paired_protein_id",
    "orthology_paired_reference_id",
    "orthology_paired_ortholog_id",
    "orthology_paired_orthogroup",
    "orthology_shared_orthogroup",
    "orthology_pair_supported",
    "orthology_support_terms",
    "orthology_conflicting_terms",
    "orthology_rule_version",
    "orthology_source",
    "orthology_source_record_id",
    "orthology_warnings",
    "phylogenetic_profile_status",
    "phylogenetic_profile_informative_species",
    "phylogenetic_profile_shared_presence",
    "phylogenetic_profile_shared_absence",
    "phylogenetic_profile_discordant",
    "phylogenetic_profile_unknown",
    "phylogenetic_profile_similarity",
    "phylogenetic_profile_pair_supported",
    "phylogenetic_profile_support_terms",
    "phylogenetic_profile_conflicting_terms",
    "phylogenetic_profile_rule_version",
    "phylogenetic_profile_source",
    "phylogenetic_profile_warnings",
    "fusion_status",
    "fusion_supporting_record_count",
    "fusion_qualifying_record_count",
    "fusion_reference_organisms",
    "fusion_protein_ids",
    "fusion_best_query_component_coverage",
    "fusion_best_candidate_component_coverage",
    "fusion_minimum_component_overlap_fraction",
    "fusion_pair_supported",
    "fusion_support_terms",
    "fusion_conflicting_terms",
    "fusion_rule_version",
    "fusion_source",
    "fusion_source_record_ids",
    "fusion_warnings",
    "known_interactions_status",
    "known_interactions_supporting_record_count",
    "known_interactions_qualifying_record_count",
    "known_interactions_direct_record_count",
    "known_interactions_physical_record_count",
    "known_interactions_biological_record_count",
    "known_interactions_interaction_types",
    "known_interactions_detection_methods",
    "known_interactions_publication_ids",
    "known_interactions_reference_organisms",
    "known_interactions_sources",
    "known_interactions_source_record_ids",
    "known_interactions_best_confidence",
    "known_interactions_pair_supported",
    "known_interactions_direct_supported",
    "known_interactions_physical_supported",
    "known_interactions_functional_association_supported",
    "known_interactions_support_terms",
    "known_interactions_conflicting_terms",
    "known_interactions_rule_version",
    "known_interactions_warnings",
    "scoring_status",
    "integrated_score",
    "provisional_score",
    "normalized_score",
    "rank",
    "tied_rank",
    "sufficient_evidence",
    "available_weight",
    "evidence_category_count",
    "evidence_component_count",
    "positive_component_count",
    "neutral_component_count",
    "negative_component_count",
    "scoring_rule_version",
    "scoring_support_terms",
    "scoring_conflicting_terms",
    "scoring_warnings",
    "score_component_genome_context",
    "score_component_operon_proxy",
    "score_component_domain_pair",
    "score_component_functional_complementarity",
    "score_component_localization",
    "score_component_orthology",
    "score_component_phylogenetic_profile",
    "score_component_fusion",
    "score_component_known_interactions",
)


def _list(values: list[str]) -> str:
    return "|".join(sorted(values))


def _value(value: Any) -> Any:
    if value is None:
        return ""
    return value.value if hasattr(value, "value") else value


def _evidence_values(
    evidence: Sequence[DomainEvidence | FunctionalEvidence | OrthologRecord],
    attribute: str,
) -> str:
    values: set[str] = set()

    for item in evidence:
        value = getattr(item, attribute)

        if value is None:
            continue

        if isinstance(value, list):
            values.update(str(_value(entry)) for entry in value)
        else:
            values.add(str(_value(value)))

    return "|".join(sorted(values))


class CandidateTableTsvWriter:
    def write(
        self,
        run_id: str,
        candidates: Sequence[CandidateProtein],
        path: Path,
        contexts: Mapping[
            tuple[str, str],
            GenomeContextEvidence,
        ]
        | None = None,
        operons: Mapping[
            tuple[str, str],
            OperonEvidence,
        ]
        | None = None,
        domains: Mapping[
            tuple[str, str],
            Sequence[DomainEvidence],
        ]
        | None = None,
        localization: Mapping[
            tuple[str, str],
            LocalizationEvidence,
        ]
        | None = None,
        functional: Mapping[
            tuple[str, str],
            Sequence[FunctionalEvidence],
        ]
        | None = None,
        orthology: Mapping[
            tuple[str, str],
            Sequence[OrthologRecord],
        ]
        | None = None,
        phylogenetic_profile: Mapping[
            tuple[str, str],
            PhylogeneticProfileEvidence,
        ]
        | None = None,
        fusion: Mapping[tuple[str, str], FusionEvidence] | None = None,
        known_interactions: Mapping[tuple[str, str], KnownInteractionEvidence] | None = None,
        integrated_scores: Mapping[tuple[str, str], IntegratedScore] | None = None,
    ) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(
            buffer, fieldnames=CANDIDATE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        context_map = contexts or {}
        operon_map = operons or {}
        domain_map = domains or {}
        localization_map = localization or {}
        functional_map = functional or {}
        orthology_map = orthology or {}
        phylogenetic_profile_map = phylogenetic_profile or {}
        fusion_map = fusion or {}
        known_interactions_map = known_interactions or {}
        integrated_score_map = integrated_scores or {}
        for candidate in sorted(candidates, key=lambda item: (item.query_id, item.protein_id)):
            pair = (candidate.query_id, candidate.protein_id)
            context = context_map.get(pair)
            operon = operon_map.get(pair)
            domain_records = domain_map.get(pair, ())
            localization_record = localization_map.get(pair)
            functional_records = functional_map.get(pair, ())
            orthology_records = orthology_map.get(pair, ())
            profile = phylogenetic_profile_map.get(pair)
            fusion_evidence = fusion_map.get(pair)
            known_interaction = known_interactions_map.get(pair)
            integrated_score = integrated_score_map.get(pair)
            component_values = (
                {
                    item.component_name: item.normalized_value
                    for item in integrated_score.component_scores
                }
                if integrated_score
                else {}
            )
            writer.writerow(
                {
                    "run_id": run_id,
                    "query_id": candidate.query_id,
                    "candidate_id": candidate.protein_id,
                    "candidate_description": candidate.description,
                    "candidate_disposition": candidate.disposition.value,
                    "disposition_reasons": _list(candidate.disposition_reasons),
                    "sequence_length": candidate.sequence_length,
                    "gene_id": _value(candidate.gene_id),
                    "locus_tag": _value(candidate.locus_tag),
                    "old_locus_tag": _value(candidate.old_locus_tag),
                    "contig": _value(candidate.contig),
                    "strand": _value(candidate.strand),
                    "has_coordinate": candidate.has_coordinate,
                    "has_annotation": candidate.has_annotation,
                    "same_contig_as_query": _value(candidate.same_contig_as_query),
                    "is_duplicate_sequence": candidate.is_duplicate_sequence,
                    "duplicate_group_id": _value(candidate.duplicate_sequence_group),
                    "is_fragment_candidate": candidate.is_fragment_candidate,
                    "fragment_reasons": _list(candidate.fragment_reasons),
                    "is_hypothetical": candidate.is_hypothetical,
                    "identifier_match_status": candidate.identifier_match_status.value,
                    "warnings": _list(candidate.warnings),
                    "same_contig": _value(context.same_contig if context else None),
                    "query_start": _value(context.query_start if context else None),
                    "query_end": _value(context.query_end if context else None),
                    "query_strand": _value(context.query_strand if context else None),
                    "candidate_start": _value(context.candidate_start if context else None),
                    "candidate_end": _value(context.candidate_end if context else None),
                    "candidate_strand": _value(context.candidate_strand if context else None),
                    "strand_relationship": _value(context.strand_relationship if context else None),
                    "relative_position": _value(context.relative_position if context else None),
                    "coordinate_position": _value(context.coordinate_position if context else None),
                    "distance_bp": _value(context.distance_bp if context else None),
                    "overlap_bp": _value(context.overlap_bp if context else None),
                    "intervening_gene_count": _value(
                        context.intervening_gene_count if context else None
                    ),
                    "intervening_feature_count": _value(
                        context.intervening_feature_count if context else None
                    ),
                    "feature_index_delta": _value(context.feature_index_delta if context else None),
                    "within_neighborhood_window": _value(
                        context.within_neighborhood_window if context else None
                    ),
                    "context_completeness": _value(
                        context.context_completeness if context else None
                    ),
                    "gene_context_status": _value(context.status if context else None),
                    "edge_to_edge_distance_bp": _value(
                        context.edge_to_edge_distance_bp if context else None
                    ),
                    "within_neighborhood_gene_count": _value(
                        context.within_neighborhood_gene_count if context else None
                    ),
                    "query_distance_to_contig_left_edge": _value(
                        context.query_distance_to_contig_left_edge if context else None
                    ),
                    "query_distance_to_contig_right_edge": _value(
                        context.query_distance_to_contig_right_edge if context else None
                    ),
                    "candidate_distance_to_contig_left_edge": _value(
                        context.candidate_distance_to_contig_left_edge if context else None
                    ),
                    "candidate_distance_to_contig_right_edge": _value(
                        context.candidate_distance_to_contig_right_edge if context else None
                    ),
                    "gene_context_rule_version": _value(
                        context.calculation_rule_version if context else None
                    ),
                    "gene_context_warnings": _list(context.warnings) if context else "",
                    "operon_status": _value(operon.status if operon else None),
                    "operon_proxy_status": _value(operon.proxy_status if operon else None),
                    "operon_same_contig": _value(operon.same_contig if operon else None),
                    "operon_same_strand": _value(operon.same_strand if operon else None),
                    "operon_is_adjacent": _value(operon.is_adjacent if operon else None),
                    "operon_intergenic_distance_bp": _value(
                        operon.intergenic_distance_bp if operon else None
                    ),
                    "operon_overlap_bp": _value(operon.overlap_bp if operon else None),
                    "operon_intervening_gene_count": _value(
                        operon.intervening_gene_count if operon else None
                    ),
                    "operon_transcriptional_order": _value(
                        operon.transcriptional_order if operon else None
                    ),
                    "operon_maximum_intergenic_distance_bp": _value(
                        operon.maximum_intergenic_distance_bp if operon else None
                    ),
                    "operon_passes_distance_threshold": _value(
                        operon.passes_distance_threshold if operon else None
                    ),
                    "operon_supporting_conditions": (
                        _list(operon.supporting_conditions) if operon else ""
                    ),
                    "operon_conflicting_conditions": (
                        _list(operon.conflicting_conditions) if operon else ""
                    ),
                    "operon_rule_version": _value(
                        operon.calculation_rule_version if operon else None
                    ),
                    "operon_rule_id": _value(operon.proxy_rule_id if operon else None),
                    "operon_warnings": _list(operon.warnings) if operon else "",
                    "domain_status": _evidence_values(
                        domain_records,
                        "status",
                    ),
                    "localization_status": _value(
                        localization_record.status if localization_record else None
                    ),
                    "localization_compartment": _value(
                        localization_record.compartment if localization_record else None
                    ),
                    "localization_query_compartment": _value(
                        localization_record.query_compartment if localization_record else None
                    ),
                    "localization_compatibility": _value(
                        localization_record.compatibility if localization_record else None
                    ),
                    "localization_signal_peptide": _value(
                        localization_record.signal_peptide if localization_record else None
                    ),
                    "localization_tm_helices": _value(
                        localization_record.transmembrane_helices if localization_record else None
                    ),
                    "localization_topology": _value(
                        localization_record.topology if localization_record else None
                    ),
                    "localization_matched_terms": (
                        _list(localization_record.matched_terms) if localization_record else ""
                    ),
                    "localization_conflicting_terms": (
                        _list(localization_record.conflicting_terms) if localization_record else ""
                    ),
                    "localization_rule_version": _value(
                        localization_record.calculation_rule_version
                        if localization_record
                        else None
                    ),
                    "localization_annotation_source": _value(
                        localization_record.annotation_source if localization_record else None
                    ),
                    "localization_annotation_confidence": _value(
                        localization_record.annotation_confidence if localization_record else None
                    ),
                    "localization_warnings": (
                        _list(localization_record.warnings) if localization_record else ""
                    ),
                    "domain_pair_matched": _evidence_values(
                        domain_records,
                        "pair_matched",
                    ),
                    "domain_pair_rule_ids": _evidence_values(
                        domain_records,
                        "pair_rule_id",
                    ),
                    "domain_candidate_accessions": _evidence_values(
                        domain_records,
                        "accession",
                    ),
                    "domain_query_accessions": _evidence_values(
                        domain_records,
                        "paired_accession",
                    ),
                    "domain_candidate_roles": _evidence_values(
                        domain_records,
                        "role",
                    ),
                    "domain_is_shared": _evidence_values(
                        domain_records,
                        "is_shared",
                    ),
                    "domain_support_terms": _evidence_values(
                        domain_records,
                        "support_terms",
                    ),
                    "domain_conflicting_terms": _evidence_values(
                        domain_records,
                        "conflicting_terms",
                    ),
                    "domain_rule_versions": _evidence_values(
                        domain_records,
                        "calculation_rule_version",
                    ),
                    "domain_warnings": _evidence_values(
                        domain_records,
                        "warnings",
                    ),
                    "functional_status": _evidence_values(
                        functional_records,
                        "status",
                    ),
                    "functional_matched": _evidence_values(
                        functional_records,
                        "matched",
                    ),
                    "functional_relationship_hints": _evidence_values(
                        functional_records,
                        "relationship_hint",
                    ),
                    "functional_rule_ids": _evidence_values(
                        functional_records,
                        "rule_id",
                    ),
                    "functional_query_roles": _evidence_values(
                        functional_records,
                        "query_role",
                    ),
                    "functional_candidate_roles": _evidence_values(
                        functional_records,
                        "candidate_role",
                    ),
                    "functional_query_matched_terms": _evidence_values(
                        functional_records,
                        "query_matched_terms",
                    ),
                    "functional_candidate_matched_terms": _evidence_values(
                        functional_records,
                        "candidate_matched_terms",
                    ),
                    "functional_support_terms": _evidence_values(
                        functional_records,
                        "support_terms",
                    ),
                    "functional_conflicting_terms": _evidence_values(
                        functional_records,
                        "conflicting_terms",
                    ),
                    "functional_rule_versions": _evidence_values(
                        functional_records,
                        "calculation_rule_version",
                    ),
                    "functional_warnings": _evidence_values(
                        functional_records,
                        "warnings",
                    ),
                    "orthology_status": _evidence_values(
                        orthology_records,
                        "status",
                    ),
                    "orthology_reference_organism": _evidence_values(
                        orthology_records,
                        "reference_organism",
                    ),
                    "orthology_relationship": _evidence_values(
                        orthology_records,
                        "relationship",
                    ),
                    "orthology_paired_protein_id": _evidence_values(
                        orthology_records,
                        "paired_protein_id",
                    ),
                    "orthology_paired_reference_id": _evidence_values(
                        orthology_records,
                        "paired_reference_id",
                    ),
                    "orthology_paired_ortholog_id": _evidence_values(
                        orthology_records,
                        "paired_ortholog_id",
                    ),
                    "orthology_paired_orthogroup": _evidence_values(
                        orthology_records,
                        "paired_orthogroup",
                    ),
                    "orthology_shared_orthogroup": _evidence_values(
                        orthology_records,
                        "shared_orthogroup",
                    ),
                    "orthology_pair_supported": _evidence_values(
                        orthology_records,
                        "pair_supported",
                    ),
                    "orthology_support_terms": _evidence_values(
                        orthology_records,
                        "support_terms",
                    ),
                    "orthology_conflicting_terms": _evidence_values(
                        orthology_records,
                        "conflicting_terms",
                    ),
                    "orthology_rule_version": _evidence_values(
                        orthology_records,
                        "calculation_rule_version",
                    ),
                    "orthology_source": _evidence_values(
                        orthology_records,
                        "source",
                    ),
                    "orthology_source_record_id": _evidence_values(
                        orthology_records,
                        "source_record_id",
                    ),
                    "orthology_warnings": _evidence_values(
                        orthology_records,
                        "warnings",
                    ),
                    "phylogenetic_profile_status": _value(profile.status if profile else None),
                    "phylogenetic_profile_informative_species": _value(
                        profile.informative_species_count if profile else None
                    ),
                    "phylogenetic_profile_shared_presence": _value(
                        profile.shared_presence_count if profile else None
                    ),
                    "phylogenetic_profile_shared_absence": _value(
                        profile.shared_absence_count if profile else None
                    ),
                    "phylogenetic_profile_discordant": _value(
                        profile.discordant_count if profile else None
                    ),
                    "phylogenetic_profile_unknown": _value(
                        profile.unknown_count if profile else None
                    ),
                    "phylogenetic_profile_similarity": _value(
                        profile.profile_similarity if profile else None
                    ),
                    "phylogenetic_profile_pair_supported": _value(
                        profile.pair_supported if profile else None
                    ),
                    "phylogenetic_profile_support_terms": (
                        _list(profile.support_terms) if profile else ""
                    ),
                    "phylogenetic_profile_conflicting_terms": (
                        _list(profile.conflicting_terms) if profile else ""
                    ),
                    "phylogenetic_profile_rule_version": _value(
                        profile.calculation_rule_version if profile else None
                    ),
                    "phylogenetic_profile_source": _value(profile.source if profile else None),
                    "phylogenetic_profile_warnings": (_list(profile.warnings) if profile else ""),
                    "fusion_status": _value(fusion_evidence.status if fusion_evidence else None),
                    "fusion_supporting_record_count": _value(
                        fusion_evidence.supporting_record_count if fusion_evidence else None
                    ),
                    "fusion_qualifying_record_count": _value(
                        fusion_evidence.qualifying_record_count if fusion_evidence else None
                    ),
                    "fusion_reference_organisms": (
                        _list(fusion_evidence.reference_organisms) if fusion_evidence else ""
                    ),
                    "fusion_protein_ids": (
                        _list(fusion_evidence.fusion_protein_ids) if fusion_evidence else ""
                    ),
                    "fusion_best_query_component_coverage": _value(
                        fusion_evidence.best_query_component_coverage if fusion_evidence else None
                    ),
                    "fusion_best_candidate_component_coverage": _value(
                        fusion_evidence.best_candidate_component_coverage
                        if fusion_evidence
                        else None
                    ),
                    "fusion_minimum_component_overlap_fraction": _value(
                        fusion_evidence.minimum_component_overlap_fraction
                        if fusion_evidence
                        else None
                    ),
                    "fusion_pair_supported": _value(
                        fusion_evidence.pair_supported if fusion_evidence else None
                    ),
                    "fusion_support_terms": (
                        _list(fusion_evidence.support_terms) if fusion_evidence else ""
                    ),
                    "fusion_conflicting_terms": (
                        _list(fusion_evidence.conflicting_terms) if fusion_evidence else ""
                    ),
                    "fusion_rule_version": _value(
                        fusion_evidence.calculation_rule_version if fusion_evidence else None
                    ),
                    "fusion_source": _value(fusion_evidence.source if fusion_evidence else None),
                    "fusion_source_record_ids": (
                        _list(fusion_evidence.source_record_ids) if fusion_evidence else ""
                    ),
                    "fusion_warnings": (_list(fusion_evidence.warnings) if fusion_evidence else ""),
                    "known_interactions_status": _value(
                        known_interaction.status if known_interaction else None
                    ),
                    "known_interactions_supporting_record_count": _value(
                        known_interaction.supporting_record_count if known_interaction else None
                    ),
                    "known_interactions_qualifying_record_count": _value(
                        known_interaction.qualifying_record_count if known_interaction else None
                    ),
                    "known_interactions_direct_record_count": _value(
                        known_interaction.direct_record_count if known_interaction else None
                    ),
                    "known_interactions_physical_record_count": _value(
                        known_interaction.physical_record_count if known_interaction else None
                    ),
                    "known_interactions_biological_record_count": _value(
                        known_interaction.biological_record_count if known_interaction else None
                    ),
                    "known_interactions_interaction_types": (
                        _list(known_interaction.interaction_types) if known_interaction else ""
                    ),
                    "known_interactions_detection_methods": (
                        _list(known_interaction.detection_methods) if known_interaction else ""
                    ),
                    "known_interactions_publication_ids": (
                        _list(known_interaction.publication_ids) if known_interaction else ""
                    ),
                    "known_interactions_reference_organisms": (
                        _list(known_interaction.reference_organisms) if known_interaction else ""
                    ),
                    "known_interactions_sources": (
                        _list(known_interaction.sources) if known_interaction else ""
                    ),
                    "known_interactions_source_record_ids": (
                        _list(known_interaction.source_record_ids) if known_interaction else ""
                    ),
                    "known_interactions_best_confidence": _value(
                        known_interaction.best_confidence if known_interaction else None
                    ),
                    "known_interactions_pair_supported": _value(
                        known_interaction.pair_supported if known_interaction else None
                    ),
                    "known_interactions_direct_supported": _value(
                        known_interaction.direct_interaction_supported
                        if known_interaction
                        else None
                    ),
                    "known_interactions_physical_supported": _value(
                        known_interaction.physical_interaction_supported
                        if known_interaction
                        else None
                    ),
                    "known_interactions_functional_association_supported": _value(
                        known_interaction.functional_association_supported
                        if known_interaction
                        else None
                    ),
                    "known_interactions_support_terms": (
                        _list(known_interaction.support_terms) if known_interaction else ""
                    ),
                    "known_interactions_conflicting_terms": (
                        _list(known_interaction.conflicting_terms) if known_interaction else ""
                    ),
                    "known_interactions_rule_version": _value(
                        known_interaction.calculation_rule_version if known_interaction else None
                    ),
                    "known_interactions_warnings": (
                        _list(known_interaction.warnings) if known_interaction else ""
                    ),
                    "scoring_status": _value(integrated_score.status if integrated_score else None),
                    "integrated_score": _value(
                        integrated_score.output_score if integrated_score else None
                    ),
                    "provisional_score": _value(
                        integrated_score.provisional_score if integrated_score else None
                    ),
                    "normalized_score": _value(
                        integrated_score.normalized_score if integrated_score else None
                    ),
                    "rank": _value(integrated_score.rank if integrated_score else None),
                    "tied_rank": _value(integrated_score.tied_rank if integrated_score else None),
                    "sufficient_evidence": _value(
                        integrated_score.sufficient_evidence if integrated_score else None
                    ),
                    "available_weight": _value(
                        integrated_score.available_weight if integrated_score else None
                    ),
                    "evidence_category_count": _value(
                        integrated_score.evidence_category_count if integrated_score else None
                    ),
                    "evidence_component_count": _value(
                        integrated_score.evidence_component_count if integrated_score else None
                    ),
                    "positive_component_count": _value(
                        integrated_score.positive_component_count if integrated_score else None
                    ),
                    "neutral_component_count": _value(
                        integrated_score.neutral_component_count if integrated_score else None
                    ),
                    "negative_component_count": _value(
                        integrated_score.negative_component_count if integrated_score else None
                    ),
                    "scoring_rule_version": _value(
                        integrated_score.calculation_rule_version if integrated_score else None
                    ),
                    "scoring_support_terms": (
                        _list(integrated_score.support_terms) if integrated_score else ""
                    ),
                    "scoring_conflicting_terms": (
                        _list(integrated_score.conflicting_terms) if integrated_score else ""
                    ),
                    "scoring_warnings": (
                        _list(integrated_score.warnings) if integrated_score else ""
                    ),
                    "score_component_genome_context": _value(
                        component_values.get("genome_context")
                    ),
                    "score_component_operon_proxy": _value(component_values.get("operon_proxy")),
                    "score_component_domain_pair": _value(component_values.get("domain_pair")),
                    "score_component_functional_complementarity": _value(
                        component_values.get("functional_complementarity")
                    ),
                    "score_component_localization": _value(component_values.get("localization")),
                    "score_component_orthology": _value(component_values.get("orthology")),
                    "score_component_phylogenetic_profile": _value(
                        component_values.get("phylogenetic_profile")
                    ),
                    "score_component_fusion": _value(component_values.get("fusion")),
                    "score_component_known_interactions": _value(
                        component_values.get("known_interactions")
                    ),
                }
            )
        output_path.write_text(buffer.getvalue(), encoding="utf-8", newline="\n")
        return output_path


class WarningSummaryTsvWriter:
    def write(self, warnings: Sequence[str], path: Path) -> Path:
        output_path = path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(("warning", "count"))
        writer.writerows(sorted(Counter(warnings).items()))
        output_path.write_text(buffer.getvalue(), encoding="utf-8", newline="\n")
        return output_path
