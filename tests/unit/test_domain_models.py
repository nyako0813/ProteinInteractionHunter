from protein_interaction_hunter.models.enums import EvidenceStatus
from protein_interaction_hunter.models.evidence import DomainEvidence


def test_domain_evidence_defaults_are_observation_only() -> None:
    evidence = DomainEvidence(protein_id="TEST_001")

    assert evidence.status is EvidenceStatus.NOT_RUN
    assert evidence.pair_matched is None
    assert evidence.is_shared is None
    assert evidence.support_terms == []
    assert evidence.conflicting_terms == []


def test_domain_evidence_records_domain_pair_without_score() -> None:
    evidence = DomainEvidence(
        status=EvidenceStatus.AVAILABLE,
        calculation_rule_version="mvp1e-domain-pair-v1",
        protein_id="CANDIDATE_001",
        source="local_table",
        accession="PF00001",
        name="Accessory domain",
        start=5,
        end=100,
        architecture_index=0,
        role="accessory_domain",
        pair_rule_id="enzyme-accessory-domain-v1",
        paired_protein_id="QUERY_001",
        paired_accession="PF00002",
        is_shared=False,
        pair_matched=True,
        support_terms=[
            "query:PF00002",
            "candidate:PF00001",
        ],
        ruleset_path="rules/domain_pairs.v1.yaml",
    )

    assert evidence.pair_matched is True
    assert evidence.pair_rule_id == "enzyme-accessory-domain-v1"
    assert evidence.quality is None


def test_domain_evidence_list_defaults_are_independent() -> None:
    first = DomainEvidence(protein_id="A")
    second = DomainEvidence(protein_id="B")

    first.support_terms.append("PF00001")

    assert second.support_terms == []