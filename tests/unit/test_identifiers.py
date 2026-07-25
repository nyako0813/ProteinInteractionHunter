"""Identifier normalization and query resolution behavior."""

import pytest

from protein_interaction_hunter.application.candidates import resolve_queries
from protein_interaction_hunter.application.identifiers import (
    IdentifierIndex,
    normalize_identifier,
)
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models import IdentifierMatchStatus, ProteinRecord


def proteins() -> list[ProteinRecord]:
    return [
        ProteinRecord(
            protein_id="Prot_A.1",
            sequence="MSTK",
            gene_id="Gene-A",
            locus_tag="LOC_A",
            aliases=["shared"],
        ),
        ProteinRecord(
            protein_id="Prot_B.1",
            sequence="MAAA",
            gene_id="Gene-B",
            locus_tag="LOC_B",
            aliases=["shared"],
        ),
    ]


def test_normalization_is_deterministic_for_case_space_and_known_prefix() -> None:
    assert normalize_identifier("  protein_id:Prot_A.1 ") == "prot_a.1"
    assert normalize_identifier("LOCUS_TAG=Loc_A") == "loc_a"


def test_exact_identifier_match_preserves_original() -> None:
    resolution = IdentifierIndex(proteins()).resolve("Prot_A.1")
    assert resolution.status is IdentifierMatchStatus.EXACT_MATCH
    assert resolution.input_identifier == "Prot_A.1"
    assert resolution.normalized_identifier == "prot_a.1"
    assert resolution.canonical_protein_id == "Prot_A.1"


def test_unique_alias_match() -> None:
    resolution = IdentifierIndex(proteins()).resolve(" gene:GENE-A ")
    assert resolution.status is IdentifierMatchStatus.UNIQUE_ALIAS_MATCH
    assert resolution.canonical_protein_id == "Prot_A.1"


def test_ambiguous_alias_is_not_guessed() -> None:
    resolution = IdentifierIndex(proteins()).resolve("SHARED")
    assert resolution.status is IdentifierMatchStatus.AMBIGUOUS_MATCH
    assert resolution.canonical_protein_id is None
    assert resolution.candidate_protein_ids == ["Prot_A.1", "Prot_B.1"]


def test_missing_and_ambiguous_queries_are_fatal() -> None:
    index = IdentifierIndex(proteins())
    assert index.resolve("absent").status is IdentifierMatchStatus.NO_MATCH
    with pytest.raises(InputValidationError, match="not found"):
        resolve_queries(["absent"], index)
    with pytest.raises(InputValidationError, match="ambiguous"):
        resolve_queries(["shared"], index)


def test_multiple_queries_and_duplicate_canonical_resolution() -> None:
    index = IdentifierIndex(proteins())
    resolved = resolve_queries(["Prot_A.1", "Gene-B"], index)
    assert [query.protein_id for query in resolved] == ["Prot_A.1", "Prot_B.1"]
    with pytest.raises(InputValidationError, match="same protein"):
        resolve_queries(["Prot_A.1", "Gene-A"], index)
