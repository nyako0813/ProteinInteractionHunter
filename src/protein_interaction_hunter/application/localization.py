"""Deterministic annotation-based localization evidence."""

from __future__ import annotations

import re
import unicodedata

from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    LocalizationEvidence,
)

LOCALIZATION_ENGINE_VERSION = "mvp1f-localization-v1"
LOCALIZATION_RULE_ID = "annotation-localization-compatibility-v1"


def normalize_localization_text(value: str | None) -> str:
    if not value:
        return ""

    normalized = unicodedata.normalize("NFKC", value).casefold()
    normalized = re.sub(r"[^0-9a-z]+", " ", normalized)
    return " ".join(normalized.split())


def _classify_compartment(text: str) -> tuple[str | None, list[str]]:
    rules = (
        ("s_layer", ("s layer", "surface layer")),
        ("secreted", ("secreted", "extracellular", "exported")),
        ("membrane", ("membrane", "membrane bound", "membrane associated")),
        ("cytosolic", ("cytosolic", "cytoplasmic", "cytoplasm")),
    )

    for compartment, terms in rules:
        matched = [term for term in terms if term in text]
        if matched:
            return compartment, matched

    return None, []


def _classify_transmembrane(
    text: str,
) -> tuple[int | None, str | None, list[str]]:
    if not text:
        return None, None, []

    if any(term in text for term in ("none", "no tm", "no transmembrane")):
        return 0, "none", ["no_transmembrane"]

    if any(
        term in text
        for term in ("multi pass", "multipass", "multi transmembrane")
    ):
        return 2, "multi_pass", ["multi_pass"]

    if any(
        term in text
        for term in ("single pass", "single transmembrane")
    ):
        return 1, "single_pass", ["single_pass"]

    if "transmembrane" in text or "tm helix" in text:
        return 1, "transmembrane", ["transmembrane"]

    return None, None, []


def _signal_peptide(
    localization_text: str,
    transmembrane_text: str,
) -> tuple[bool | None, list[str]]:
    combined = f"{localization_text} {transmembrane_text}".strip()

    if "no signal peptide" in combined:
        return False, ["no_signal_peptide"]
    
    if "signal peptide" in combined or "secretion signal" in combined:
        return True, ["signal_peptide"]

    return None, []


def evaluate_localization(
    query_protein_id: str,
    candidate_protein_id: str,
    annotation_index: dict[str, AnnotationRecord],
) -> LocalizationEvidence:
    """Evaluate annotation-only localization without score propagation."""

    query = annotation_index.get(query_protein_id)
    candidate = annotation_index.get(candidate_protein_id)

    provenance = [
        EvidenceProvenance(
            source_name="annotation_localization",
            source_version=LOCALIZATION_ENGINE_VERSION,
            method="deterministic_annotation_term_mapping",
            metadata={
                "ranking_propagation": False,
                "score_propagation": False,
                "evidence_tier_propagation": False,
                "contradiction_propagation": False,
            },
        )
    ]

    if candidate is None:
        return LocalizationEvidence(
            status=EvidenceStatus.MISSING,
            origin=EvidenceOrigin.ANNOTATION,
            calculation_rule_version=LOCALIZATION_ENGINE_VERSION,
            protein_id=candidate_protein_id,
            rule_id=LOCALIZATION_RULE_ID,
            conflicting_terms=["missing_candidate_annotation"],
            warnings=["missing_localization_annotation"],
            provenance=provenance,
        )

    query_text = normalize_localization_text(
        query.localization_annotation if query else None
    )
    candidate_text = normalize_localization_text(
        candidate.localization_annotation
    )
    transmembrane_text = normalize_localization_text(
        candidate.transmembrane_annotation
    )

    query_compartment, query_terms = _classify_compartment(query_text)
    candidate_compartment, candidate_terms = _classify_compartment(
        candidate_text
    )
    helices, topology, topology_terms = _classify_transmembrane(
        transmembrane_text
    )
    signal_peptide, signal_terms = _signal_peptide(
        candidate_text,
        transmembrane_text,
    )

    matched_terms = sorted(
        set(
            [f"query:{term}" for term in query_terms]
            + [f"candidate:{term}" for term in candidate_terms]
            + [f"candidate:{term}" for term in topology_terms]
            + [f"candidate:{term}" for term in signal_terms]
        )
    )

    conflicting_terms: list[str] = []

    if query_compartment and candidate_compartment:
        compatibility = query_compartment == candidate_compartment
        if not compatibility:
            conflicting_terms.append("different_compartment")
    else:
        compatibility = None

    if not candidate_text and not transmembrane_text:
        return LocalizationEvidence(
            status=EvidenceStatus.MISSING,
            origin=EvidenceOrigin.ANNOTATION,
            calculation_rule_version=LOCALIZATION_ENGINE_VERSION,
            protein_id=candidate_protein_id,
            query_compartment=query_compartment,
            compatibility=None,
            rule_id=LOCALIZATION_RULE_ID,
            annotation_source=candidate.annotation_source,
            annotation_confidence=candidate.annotation_confidence,
            conflicting_terms=["missing_localization_annotation"],
            warnings=["missing_localization_annotation"],
            provenance=provenance,
        )

    return LocalizationEvidence(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.ANNOTATION,
        calculation_rule_version=LOCALIZATION_ENGINE_VERSION,
        protein_id=candidate_protein_id,
        compartment=candidate_compartment,
        signal_peptide=signal_peptide,
        transmembrane_helices=helices,
        topology=topology,
        compatibility=compatibility,
        query_compartment=query_compartment,
        candidate_compartment=candidate_compartment,
        localization_annotation=candidate.localization_annotation,
        transmembrane_annotation=candidate.transmembrane_annotation,
        matched_terms=matched_terms,
        conflicting_terms=conflicting_terms,
        rule_id=LOCALIZATION_RULE_ID,
        annotation_source=candidate.annotation_source,
        annotation_confidence=candidate.annotation_confidence,
        provenance=provenance,
    )