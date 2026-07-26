"""Deterministic local known-interaction pair evaluation."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Sequence

from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    KnownInteractionEvidence,
    KnownInteractionObservation,
)

KNOWN_INTERACTIONS_ENGINE_VERSION = "mvp1j-known-interactions-v1"
KnownInteractionIndex = dict[tuple[str, str], list[KnownInteractionObservation]]

_METHOD_ALIASES = {
    "y2h": "yeast_two_hybrid",
    "yeast_2_hybrid": "yeast_two_hybrid",
    "yeast_two_hybrid": "yeast_two_hybrid",
    "pulldown": "pull_down",
    "pull_down": "pull_down",
    "co_ip": "co_immunoprecipitation",
    "co_immunoprecipitation": "co_immunoprecipitation",
    "ap_ms": "affinity_purification_mass_spectrometry",
    "affinity_purification_mass_spectrometry": ("affinity_purification_mass_spectrometry"),
    "xl_ms": "crosslinking_mass_spectrometry",
    "crosslinking_mass_spectrometry": "crosslinking_mass_spectrometry",
    "spr": "surface_plasmon_resonance",
    "surface_plasmon_resonance": "surface_plasmon_resonance",
    "itc": "isothermal_titration_calorimetry",
    "isothermal_titration_calorimetry": "isothermal_titration_calorimetry",
    "bli": "biolayer_interferometry",
    "biolayer_interferometry": "biolayer_interferometry",
    "biotin_proximity_labeling": "proximity_labeling",
    "proximity_labeling": "proximity_labeling",
    "genetic_interaction": "genetic_interaction",
    "co_expression": "co_expression",
    "database_inference": "database_inference",
}


def normalize_detection_method(value: str | None) -> str | None:
    if value is None or not value.strip():
        return None
    slug = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
    return _METHOD_ALIASES.get(slug, "other")


def canonical_interaction_pair(protein_a: str, protein_b: str) -> tuple[str, str]:
    return (protein_a, protein_b) if protein_a <= protein_b else (protein_b, protein_a)


def build_known_interaction_index(
    records: list[KnownInteractionObservation],
) -> KnownInteractionIndex:
    index: dict[tuple[str, str], list[KnownInteractionObservation]] = defaultdict(list)
    for record in records:
        index[canonical_interaction_pair(record.protein_a_id, record.protein_b_id)].append(record)
    return {
        pair: sorted(
            pair_records,
            key=lambda record: (
                record.source,
                record.source_record_id,
                record.reference_organism,
                record.interaction_type,
                record.normalized_detection_method or "",
            ),
        )
        for pair, pair_records in sorted(index.items())
    }


def _is_direct(record: KnownInteractionObservation) -> bool:
    return record.interaction_type in {"direct", "physical"} or record.is_direct is True


def _is_physical(record: KnownInteractionObservation) -> bool:
    return (
        _is_direct(record)
        or record.interaction_type in {"physical", "co_complex"}
        or record.is_physical is True
    )


def _is_functional(record: KnownInteractionObservation) -> bool:
    return record.interaction_type in {
        "genetic",
        "functional_association",
        "co_expression",
        "predicted",
    }


def _category_support(
    qualifying_count: int,
    undetermined_count: int,
) -> bool | None:
    if qualifying_count:
        return True
    if undetermined_count:
        return None
    return False


def evaluate_known_interaction_pair(
    query_protein_id: str,
    candidate_protein_id: str,
    interaction_index: KnownInteractionIndex,
    *,
    minimum_supporting_records: int,
    minimum_direct_records: int,
    accepted_interaction_types: Sequence[str],
    accepted_evidence_methods: Sequence[str],
    excluded_evidence_methods: Sequence[str],
    minimum_confidence: float | None,
) -> KnownInteractionEvidence:
    records = interaction_index.get(
        canonical_interaction_pair(query_protein_id, candidate_protein_id), []
    )
    accepted_types = set(accepted_interaction_types)
    accepted_methods = {
        normalized
        for value in accepted_evidence_methods
        if (normalized := normalize_detection_method(value)) is not None
    }
    excluded_methods = {
        normalized
        for value in excluded_evidence_methods
        if (normalized := normalize_detection_method(value)) is not None
    }
    provenance = [
        EvidenceProvenance(
            source_name="local_known_interactions_table",
            source_version=KNOWN_INTERACTIONS_ENGINE_VERSION,
            method="canonical_pair_type_method_confidence_evaluation",
            metadata={
                "minimum_supporting_records": minimum_supporting_records,
                "minimum_direct_records": minimum_direct_records,
                "accepted_interaction_types": sorted(accepted_types),
                "accepted_evidence_methods": sorted(accepted_methods),
                "excluded_evidence_methods": sorted(excluded_methods),
                "minimum_confidence": minimum_confidence,
                "ranking_propagation": False,
                "score_propagation": False,
                "relationship_propagation": False,
                "evidence_tier_propagation": False,
            },
        )
    ]
    if not records:
        return KnownInteractionEvidence(
            status=EvidenceStatus.AVAILABLE,
            origin=EvidenceOrigin.EXACT_PAIR,
            query_protein_id=query_protein_id,
            candidate_protein_id=candidate_protein_id,
            supporting_record_count=0,
            qualifying_record_count=0,
            direct_record_count=0,
            physical_record_count=0,
            biological_record_count=0,
            independent_publication_count=0,
            independent_source_count=0,
            pair_supported=False,
            direct_interaction_supported=False,
            physical_interaction_supported=False,
            functional_association_supported=False,
            conflicting_terms=["no_known_interaction_record"],
            calculation_rule_version=KNOWN_INTERACTIONS_ENGINE_VERSION,
            provenance=provenance,
        )

    qualifying: list[KnownInteractionObservation] = []
    undetermined: list[KnownInteractionObservation] = []
    conflicting_terms: set[str] = set()
    missing_confidence_count = 0
    missing_method_count = 0
    uncertain_mapping_count = 0
    for record in records:
        known_failure = False
        unknown = False
        method = record.normalized_detection_method
        if record.interaction_type not in accepted_types:
            conflicting_terms.add("unsupported_interaction_type")
            known_failure = True
        if method is not None and method in excluded_methods:
            conflicting_terms.add("excluded_detection_method")
            known_failure = True
        if accepted_methods:
            if method is None:
                missing_method_count += 1
                unknown = True
            elif method not in accepted_methods:
                conflicting_terms.add("unsupported_detection_method")
                known_failure = True
        if minimum_confidence is not None:
            if record.confidence is None:
                conflicting_terms.add("missing_confidence")
                missing_confidence_count += 1
                unknown = True
            elif record.confidence < minimum_confidence:
                conflicting_terms.add("low_confidence_record")
                known_failure = True
        if record.identifier_mapping_status != "mapped":
            conflicting_terms.add("identifier_mapping_uncertain")
            uncertain_mapping_count += 1
            unknown = True
        if not known_failure and unknown:
            undetermined.append(record)
        elif not known_failure:
            qualifying.append(record)

    direct = [record for record in qualifying if _is_direct(record)]
    physical = [record for record in qualifying if _is_physical(record)]
    biological = [record for record in qualifying if record.is_biological is True]
    functional = [record for record in qualifying if _is_functional(record)]
    undetermined_direct = [record for record in undetermined if _is_direct(record)]
    undetermined_physical = [record for record in undetermined if _is_physical(record)]
    undetermined_functional = [record for record in undetermined if _is_functional(record)]

    supporting_met = len(qualifying) >= minimum_supporting_records
    direct_met = minimum_direct_records == 0 or len(direct) >= minimum_direct_records
    if supporting_met and direct_met:
        status = EvidenceStatus.AVAILABLE
        pair_supported: bool | None = True
    else:
        supporting_possible = len(qualifying) + len(undetermined) >= minimum_supporting_records
        direct_possible = (
            minimum_direct_records == 0
            or len(direct) + len(undetermined_direct) >= minimum_direct_records
        )
        if supporting_possible and direct_possible:
            status = EvidenceStatus.MISSING
            pair_supported = None
        else:
            status = EvidenceStatus.AVAILABLE
            pair_supported = False

    direct_supported = _category_support(len(direct), len(undetermined_direct))
    physical_supported = _category_support(len(physical), len(undetermined_physical))
    functional_supported = _category_support(len(functional), len(undetermined_functional))
    support_terms: set[str] = set()
    if qualifying:
        support_terms.add("known_interaction_record")
    if direct:
        support_terms.add("direct_interaction_record")
    if physical:
        support_terms.add("physical_interaction_record")
    if functional:
        support_terms.add("functional_association_record")
    if len(qualifying) > 1:
        support_terms.add("multiple_supporting_records")
    qualifying_sources = {record.source for record in qualifying}
    qualifying_publications = {
        record.publication_id for record in qualifying if record.publication_id
    }
    if len(qualifying_sources) > 1:
        support_terms.add("multiple_independent_sources")
    if len(qualifying_publications) > 1:
        support_terms.add("multiple_publications")
    if any(
        record.normalized_detection_method not in {None, "database_inference", "other"}
        for record in qualifying
    ):
        support_terms.add("experimentally_supported_interaction")
    if not supporting_met:
        conflicting_terms.add("insufficient_supporting_records")
    if not direct_met:
        conflicting_terms.add("insufficient_direct_records")
    if functional_supported is True and physical_supported is not True:
        conflicting_terms.add("functional_association_only")
    if qualifying and all(record.interaction_type == "predicted" for record in qualifying):
        conflicting_terms.add("predicted_interaction_only")

    confidences = [record.confidence for record in records if record.confidence is not None]
    methods = sorted(
        {
            record.normalized_detection_method
            for record in records
            if record.normalized_detection_method
        }
    )
    warnings: list[str] = []
    if missing_confidence_count:
        warnings.append(
            f"known_interaction_records_with_missing_confidence:{missing_confidence_count}"
        )
    if missing_method_count:
        warnings.append(
            f"known_interaction_records_with_missing_detection_method:{missing_method_count}"
        )
    if uncertain_mapping_count:
        warnings.append(
            f"known_interaction_records_with_uncertain_mapping:{uncertain_mapping_count}"
        )
    return KnownInteractionEvidence(
        status=status,
        origin=EvidenceOrigin.EXACT_PAIR,
        query_protein_id=query_protein_id,
        candidate_protein_id=candidate_protein_id,
        supporting_record_count=len(records),
        qualifying_record_count=len(qualifying),
        direct_record_count=len(direct),
        physical_record_count=len(physical),
        biological_record_count=len(biological),
        independent_publication_count=len(qualifying_publications),
        independent_source_count=len(qualifying_sources),
        interaction_types=sorted({record.interaction_type for record in records}),
        detection_methods=methods,
        publication_ids=sorted(
            {record.publication_id for record in records if record.publication_id}
        ),
        reference_organisms=sorted({record.reference_organism for record in records}),
        sources=sorted({record.source for record in records}),
        source_record_ids=sorted({record.source_record_id for record in records}),
        best_confidence=max(confidences) if confidences else None,
        pair_supported=pair_supported,
        direct_interaction_supported=direct_supported,
        physical_interaction_supported=physical_supported,
        functional_association_supported=functional_supported,
        support_terms=sorted(support_terms),
        conflicting_terms=sorted(conflicting_terms),
        calculation_rule_version=KNOWN_INTERACTIONS_ENGINE_VERSION,
        warnings=sorted(warnings),
        provenance=provenance,
    )
