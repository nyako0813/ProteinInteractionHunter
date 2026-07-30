from pathlib import Path

from protein_interaction_hunter.config import load_config

ROOT = Path(__file__).resolve().parents[2]
CONFIG = (
    ROOT
    / "data"
    / "pilot"
    / "methanosarcina_acetivorans_MA_4115"
    / "config"
    / "pilot_formal_scoring.yaml"
)


def test_ma4115_formal_config_preserves_engine_and_scoring_contract() -> None:
    config = load_config(CONFIG)

    assert config.gene_context.enabled is True
    assert config.localization.enabled is True
    assert config.orthology.enabled is True
    assert config.phylogenetic_profile.enabled is True
    assert config.functional_complementarity.enabled is False
    assert config.domains.enabled is False
    assert config.known_interactions.enabled is False
    assert config.fusion.enabled is False
    assert config.structure_prediction_queue.enabled is False
    assert config.structure_prediction_queue.automatic_structure_prediction is False

    assert config.scoring.enabled is True
    assert config.scoring.rule_version == "mvp1k-integrated-scoring-v1"
    assert config.scoring.minimum_evidence_weight == 1.0
    assert config.scoring.minimum_evidence_categories == 2
    assert config.scoring.missing_policy == "exclude_from_denominator"
    assert config.scoring.tie_precision == 8
    assert config.scoring.weights.model_dump() == {
        "genome_context": 1.0,
        "operon_proxy": 1.0,
        "domain_pair": 1.0,
        "functional_complementarity": 1.0,
        "localization": 0.5,
        "orthology": 0.75,
        "phylogenetic_profile": 1.0,
        "fusion": 1.5,
        "known_interactions": 1.5,
    }
    assert config.scoring.category_caps.model_dump() == {
        "genomic_context": 1.5,
        "functional_annotation": 1.5,
        "cellular_compatibility": 0.5,
        "evolutionary": 2.0,
        "direct_interaction": 2.0,
    }
    assert config.scoring.penalties.model_dump() == {
        "contradictory_evidence": 0.25,
        "ambiguous_mapping": 0.10,
    }

    assert config.evidence_tiers.enabled is True
    assert config.evidence_tiers.rule_version == "mvp1l-evidence-tiers-v1"
    assert config.evidence_tiers.predicted_only_tier_cap == "tier_3"
    assert config.evidence_tiers.explicit_conflict_tier_cap == "tier_3"
    assert config.evidence_tiers.functional_association_only_tier_cap == "tier_3"
    assert config.evidence_tiers.insufficient_evidence_tier == "unclassified"
