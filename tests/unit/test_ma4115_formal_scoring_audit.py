from __future__ import annotations

import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.audit_ma4115_formal_scoring import independently_score  # noqa: E402


def _bundle() -> dict[str, Any]:
    return {
        "candidate": {"identifier_match_status": "exact_match"},
        "genome_context": [
            {
                "status": "available",
                "distance_bp": 5,
                "coordinate_position": "left_of_query",
                "overlap_bp": 0,
                "same_contig": True,
                "within_neighborhood_window": True,
                "intervening_gene_count": 0,
            }
        ],
        "operon": [{"status": "available", "proxy_status": "supported"}],
        "localization": [{"status": "available", "compatibility": True}],
        "orthology": [
            {
                "status": "available",
                "pair_supported": True,
                "paralog_ambiguity": False,
            }
        ],
        "phylogenetic_profile": [
            {
                "status": "available",
                "profile_similarity": 0.8,
                "shared_presence_count": 8,
                "shared_absence_count": 2,
                "informative_species_count": 10,
                "discordant_count": 0,
                "conflicting_terms": [],
            }
        ],
    }


def test_independent_audit_applies_category_grouping_and_exact_cap() -> None:
    score = independently_score(_bundle())

    assert score["category_count"] == 3
    assert score["component_count"] == 5
    assert score["components"]["genome_context"]["effective_weight"] == 0.75
    assert score["components"]["operon_proxy"]["effective_weight"] == 0.75
    assert score["category_scores"]["genomic_context"]["available_weight"] == 1.5
    assert score["category_scores"]["evolutionary"]["available_weight"] == 1.75
    assert score["sufficient"] is True


def test_independent_audit_excludes_missing_and_invalid_values() -> None:
    bundle = _bundle()
    bundle["genome_context"] = [{"status": "failed"}]
    bundle["operon"] = [{"status": "failed", "proxy_status": "unknown"}]
    bundle["localization"] = [{"status": "missing", "compatibility": None}]
    bundle["orthology"] = [{"status": "missing", "pair_supported": None}]
    bundle["phylogenetic_profile"] = [
        {
            "status": "available",
            "profile_similarity": None,
            "shared_presence_count": 0,
            "shared_absence_count": 0,
            "informative_species_count": 0,
            "discordant_count": 0,
            "conflicting_terms": ["insufficient_informative_species"],
        }
    ]

    score = independently_score(bundle)

    assert score["available_weight"] == 0.0
    assert score["category_count"] == 0
    assert score["output_score"] is None


def test_independent_audit_preserves_negative_evidence_without_missing_penalty() -> None:
    bundle = deepcopy(_bundle())
    bundle["localization"][0]["compatibility"] = False
    bundle["phylogenetic_profile"][0].update(
        {
            "profile_similarity": 0.1,
            "shared_presence_count": 1,
            "shared_absence_count": 0,
            "informative_species_count": 10,
            "discordant_count": 9,
        }
    )

    score = independently_score(bundle)

    assert score["negative_count"] == 2
    assert score["components"]["localization"]["value"] == -0.25
    assert score["components"]["phylogenetic_profile"]["value"] == -0.25
    assert score["normalized_score"] >= 0.0
