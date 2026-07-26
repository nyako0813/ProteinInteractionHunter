from protein_interaction_hunter.models.enums import EvidenceStatus
from protein_interaction_hunter.models.evidence import (
    LocalizationEvidence,
)


def test_localization_evidence_defaults_are_shadow_only() -> None:
    evidence = LocalizationEvidence(protein_id="TEST_001")

    assert evidence.status is EvidenceStatus.NOT_RUN
    assert evidence.compatibility is None
    assert evidence.signal_peptide is None
    assert evidence.transmembrane_helices is None
    assert evidence.matched_terms == []
    assert evidence.conflicting_terms == []


def test_localization_evidence_records_annotation_result() -> None:
    evidence = LocalizationEvidence(
        status=EvidenceStatus.AVAILABLE,
        calculation_rule_version="mvp1f-localization-v1",
        protein_id="MEM_001",
        compartment="membrane",
        signal_peptide=False,
        transmembrane_helices=2,
        topology="multi_pass",
        compatibility=False,
        query_compartment="cytosolic",
        candidate_compartment="membrane",
        localization_annotation="membrane",
        transmembrane_annotation="multi-pass",
        matched_terms=[
            "candidate:membrane",
            "candidate:multi-pass",
        ],
        conflicting_terms=[
            "different_compartment",
        ],
        rule_id="annotation-localization-compatibility-v1",
        annotation_source="synthetic_curated",
        annotation_confidence=0.8,
    )

    assert evidence.compartment == "membrane"
    assert evidence.transmembrane_helices == 2
    assert evidence.compatibility is False
    assert evidence.quality is None


def test_localization_list_defaults_are_independent() -> None:
    first = LocalizationEvidence(protein_id="A")
    second = LocalizationEvidence(protein_id="B")

    first.matched_terms.append("cytosolic")

    assert second.matched_terms == []