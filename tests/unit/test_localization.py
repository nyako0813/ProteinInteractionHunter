from protein_interaction_hunter.application.localization import (
    evaluate_localization,
    normalize_localization_text,
)
from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.enums import EvidenceStatus


def annotation(
    protein_id: str,
    localization: str | None,
    transmembrane: str | None = None,
) -> AnnotationRecord:
    return AnnotationRecord(
        protein_id=protein_id,
        localization_annotation=localization,
        transmembrane_annotation=transmembrane,
        annotation_source="synthetic_curated",
        annotation_confidence=0.8,
        status=EvidenceStatus.AVAILABLE,
    )


def test_normalize_localization_text() -> None:
    assert normalize_localization_text("  Multi-pass / MEMBRANE ") == (
        "multi pass membrane"
    )


def test_same_cytosolic_compartment_is_compatible() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic", "none"),
        "CANDIDATE": annotation("CANDIDATE", "cytoplasmic", "none"),
    }

    evidence = evaluate_localization(
        "QUERY",
        "CANDIDATE",
        index,
    )

    assert evidence.status is EvidenceStatus.AVAILABLE
    assert evidence.query_compartment == "cytosolic"
    assert evidence.candidate_compartment == "cytosolic"
    assert evidence.compatibility is True
    assert evidence.transmembrane_helices == 0


def test_membrane_candidate_differs_from_cytosolic_query() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic", "none"),
        "MEM": annotation("MEM", "membrane", "multi-pass"),
    }

    evidence = evaluate_localization("QUERY", "MEM", index)

    assert evidence.compartment == "membrane"
    assert evidence.topology == "multi_pass"
    assert evidence.transmembrane_helices == 2
    assert evidence.compatibility is False
    assert evidence.conflicting_terms == ["different_compartment"]


def test_signal_peptide_is_detected() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic"),
        "SECRETED": annotation(
            "SECRETED",
            "secreted protein with signal peptide",
        ),
    }

    evidence = evaluate_localization(
        "QUERY",
        "SECRETED",
        index,
    )

    assert evidence.compartment == "secreted"
    assert evidence.signal_peptide is True


def test_missing_candidate_annotation_is_missing() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic"),
    }

    evidence = evaluate_localization(
        "QUERY",
        "UNKNOWN",
        index,
    )

    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.compatibility is None


def test_present_annotation_without_localization_is_missing() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic"),
        "UNKNOWN": annotation("UNKNOWN", None, None),
    }

    evidence = evaluate_localization(
        "QUERY",
        "UNKNOWN",
        index,
    )

    assert evidence.status is EvidenceStatus.MISSING
    assert evidence.conflicting_terms == [
        "missing_localization_annotation"
    ]

def test_no_signal_peptide_is_not_positive() -> None:
    index = {
        "QUERY": annotation("QUERY", "cytosolic"),
        "CANDIDATE": annotation(
            "CANDIDATE",
            "cytosolic",
            "no signal peptide",
        ),
    }

    evidence = evaluate_localization(
        "QUERY",
        "CANDIDATE",
        index,
    )

    assert evidence.signal_peptide is False
    assert "candidate:no_signal_peptide" in evidence.matched_terms