import pytest
from pydantic import ValidationError

from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import OrthologRecord


def test_ortholog_record_defaults_are_shadow_only() -> None:
    record = OrthologRecord(
        protein_id="CANDIDATE_001",
        reference_id="reference_archaea",
    )

    assert record.status is EvidenceStatus.NOT_RUN
    assert record.pair_supported is None
    assert record.shared_orthogroup is None
    assert record.paralog_ambiguity is False
    assert record.support_terms == []
    assert record.conflicting_terms == []


def test_ortholog_record_represents_supported_pair() -> None:
    record = OrthologRecord(
        status=EvidenceStatus.AVAILABLE,
        origin=EvidenceOrigin.ORTHOLOG_TRANSFERRED,
        calculation_rule_version="mvp1g-orthology-v1",
        protein_id="CANDIDATE_001",
        reference_id="reference_archaea",
        ortholog_id="REF_CANDIDATE",
        identity=0.72,
        query_coverage=0.91,
        subject_coverage=0.88,
        evalue=1e-40,
        orthogroup="OG0001",
        reference_organism="Methanosarcina mazei",
        relationship="reciprocal_best_hit",
        paired_protein_id="QUERY_001",
        paired_reference_id="reference_archaea",
        paired_ortholog_id="REF_QUERY",
        paired_orthogroup="OG0001",
        shared_orthogroup=True,
        pair_supported=True,
        support_terms=["shared_orthogroup"],
        source="synthetic_local_table",
        source_record_id="row-1",
    )

    assert record.pair_supported is True
    assert record.shared_orthogroup is True
    assert record.identity == 0.72
    assert record.quality is None


def test_ortholog_record_rejects_fraction_above_one() -> None:
    with pytest.raises(ValidationError):
        OrthologRecord(
            protein_id="CANDIDATE_001",
            reference_id="reference_archaea",
            identity=1.01,
        )


def test_ortholog_list_defaults_are_independent() -> None:
    first = OrthologRecord(
        protein_id="A",
        reference_id="REF",
    )
    second = OrthologRecord(
        protein_id="B",
        reference_id="REF",
    )

    first.support_terms.append("shared_orthogroup")

    assert second.support_terms == []
