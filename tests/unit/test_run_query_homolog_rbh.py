from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.run_query_homolog_rbh import Hit, classify_species  # noqa: E402

POLICY = {
    "search": {
        "evalue_maximum": 1e-5,
        "minimum_bit_score": 40,
        "minimum_percent_identity": 20,
        "minimum_alignment_length": 60,
        "minimum_query_coverage": 0.5,
        "minimum_subject_coverage": 0.5,
        "minimum_length_ratio": 0.5,
        "maximum_length_ratio": 2.0,
        "best_hit_tie_relative_bit_score": 0.01,
        "secondary_hit_relative_bit_score": 0.9,
    },
    "classification": {
        "fragment_length_ratio_maximum": 0.65,
        "multi_copy_minimum_secondary_hits": 1,
    },
}


def _hit(subject: str, bits: float = 100, *, length: int = 180, slen: int = 200) -> Hit:
    return Hit("Q", subject, 40.0, length, 200, slen, 1e-20, bits)


def test_no_hit_is_not_biological_absence() -> None:
    result = classify_species([], {}, query_normalized_id="Q", policy=POLICY)
    assert result[0] == "no_detectable_homolog"


def test_fragment_only() -> None:
    result = classify_species(
        [_hit("T", length=50, slen=70)], {}, query_normalized_id="Q", policy=POLICY
    )
    assert result[0] == "fragment_only"


def test_weak_homolog() -> None:
    weak = Hit("Q", "T", 10.0, 180, 200, 200, 1e-20, 100)
    assert classify_species([weak], {}, query_normalized_id="Q", policy=POLICY)[0] == "weak_homolog"


def test_unique_rbh() -> None:
    forward = _hit("T")
    reciprocal = Hit("T", "Q", 40.0, 180, 200, 200, 1e-20, 100)
    result = classify_species(
        [forward], {"T": [reciprocal]}, query_normalized_id="Q", policy=POLICY
    )
    assert result[:3] == (
        "unique_RBH",
        "conservative_thresholds_and_unique_reciprocal_best_hit",
        True,
    )


def test_forward_tie_is_ambiguous() -> None:
    hits = [_hit("T1", 100), _hit("T2", 99.5)]
    result = classify_species(hits, {}, query_normalized_id="Q", policy=POLICY)
    assert result[0] == "ambiguous_RBH"


def test_near_top_secondary_is_multi_copy() -> None:
    hits = [_hit("T1", 100), _hit("T2", 92)]
    reciprocal = {"T1": [Hit("T1", "Q", 40.0, 180, 200, 200, 1e-20, 100)]}
    result = classify_species(hits, reciprocal, query_normalized_id="Q", policy=POLICY)
    assert result[0] == "multi_copy_homolog_family"
    assert result[4] == 1


def test_candidate_without_reciprocal_support() -> None:
    reciprocal = {"T": [_hit("OTHER")]}
    result = classify_species([_hit("T")], reciprocal, query_normalized_id="Q", policy=POLICY)
    assert result[0] == "candidate_homolog"
