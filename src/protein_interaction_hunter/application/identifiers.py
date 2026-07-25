"""Deterministic cross-source identifier normalization and resolution."""

import re
import unicodedata
from collections import defaultdict
from collections.abc import Iterable

from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.enums import IdentifierMatchStatus
from protein_interaction_hunter.models.genome import GeneCoordinate
from protein_interaction_hunter.models.identity import IdentifierAlias, IdentifierResolution
from protein_interaction_hunter.models.protein import ProteinRecord

NORMALIZATION_RULE_VERSION = "mvp1a-id-v1"
_KNOWN_PREFIX = re.compile(
    r"^(?:protein_id|protein|gene_id|gene|locus_tag|locus|old_locus_tag|cds|id|parent)\s*[:=]\s*",
    flags=re.IGNORECASE,
)


def normalize_identifier(value: str) -> str:
    """Normalize for matching while preserving the original identifier separately."""
    normalized = unicodedata.normalize("NFKC", value).strip()
    return _KNOWN_PREFIX.sub("", normalized).casefold()


class IdentifierIndex:
    """Many-to-many alias index that never guesses through an ambiguity."""

    def __init__(self, proteins: list[ProteinRecord]) -> None:
        self._exact = {protein.protein_id: protein.protein_id for protein in proteins}
        self._aliases: dict[str, set[str]] = defaultdict(set)
        self._records: dict[str, list[IdentifierAlias]] = defaultdict(list)
        for protein in proteins:
            self._bind(protein.protein_id, self.protein_identifiers(protein))

    @staticmethod
    def protein_identifiers(protein: ProteinRecord) -> list[tuple[str, str, str]]:
        values = [(protein.protein_id, "protein_id", "fasta")]
        if protein.gene_id:
            values.append((protein.gene_id, "gene_id", "fasta"))
        if protein.locus_tag:
            values.append((protein.locus_tag, "locus_tag", "fasta"))
        values.extend((alias, "alias", "fasta") for alias in protein.aliases)
        return values

    @staticmethod
    def gff_identifiers(record: GeneCoordinate) -> list[tuple[str, str, str]]:
        return [
            (value, kind, "gff3")
            for value, kind in (
                (record.protein_id, "protein_id"),
                (record.feature_id, "gff_id"),
                (record.parent_id, "gff_parent"),
                (record.locus_tag, "locus_tag"),
                (record.old_locus_tag, "old_locus_tag"),
            )
            if value
        ]

    @staticmethod
    def annotation_identifiers(record: AnnotationRecord) -> list[tuple[str, str, str]]:
        values = [(record.protein_id, "protein_id", "annotation")]
        if record.gene_name:
            values.append((record.gene_name, "gene_id", "annotation"))
        if record.locus_tag:
            values.append((record.locus_tag, "locus_tag", "annotation"))
        return values

    def _bind(self, canonical_id: str, values: Iterable[tuple[str, str, str]]) -> None:
        for original, kind, source in values:
            normalized = normalize_identifier(original)
            if not normalized:
                continue
            self._aliases[normalized].add(canonical_id)
            alias = IdentifierAlias(
                original=original, normalized=normalized, kind=kind, source=source
            )
            if alias not in self._records[normalized]:
                self._records[normalized].append(alias)

    def _candidates(self, values: Iterable[tuple[str, str, str]]) -> set[str]:
        candidates: set[str] = set()
        for original, _, _ in values:
            candidates.update(self._aliases.get(normalize_identifier(original), set()))
        return candidates

    def _add_external(self, groups: list[list[tuple[str, str, str]]]) -> None:
        for _ in range(2):
            for values in groups:
                candidates = self._candidates(values)
                if len(candidates) == 1:
                    self._bind(next(iter(candidates)), values)
                elif len(candidates) > 1:
                    for canonical_id in candidates:
                        self._bind(canonical_id, values)

    def add_gff_records(self, records: list[GeneCoordinate]) -> None:
        self._add_external([self.gff_identifiers(record) for record in records])

    def add_annotations(self, records: list[AnnotationRecord]) -> None:
        self._add_external([self.annotation_identifiers(record) for record in records])

    def resolve(self, identifier: str) -> IdentifierResolution:
        stripped = identifier.strip()
        normalized = normalize_identifier(stripped)
        if stripped in self._exact:
            canonical = self._exact[stripped]
            return IdentifierResolution(
                input_identifier=identifier,
                normalized_identifier=normalized,
                status=IdentifierMatchStatus.EXACT_MATCH,
                canonical_protein_id=canonical,
                candidate_protein_ids=[canonical],
                matched_aliases=self.aliases_for(normalized),
            )
        candidates = sorted(self._aliases.get(normalized, set()))
        if len(candidates) == 1:
            return IdentifierResolution(
                input_identifier=identifier,
                normalized_identifier=normalized,
                status=IdentifierMatchStatus.UNIQUE_ALIAS_MATCH,
                canonical_protein_id=candidates[0],
                candidate_protein_ids=candidates,
                matched_aliases=self.aliases_for(normalized),
            )
        return IdentifierResolution(
            input_identifier=identifier,
            normalized_identifier=normalized,
            status=(
                IdentifierMatchStatus.AMBIGUOUS_MATCH
                if candidates
                else IdentifierMatchStatus.NO_MATCH
            ),
            candidate_protein_ids=candidates,
            matched_aliases=self.aliases_for(normalized),
        )

    def aliases_for(self, normalized: str) -> list[IdentifierAlias]:
        return sorted(
            self._records.get(normalized, []),
            key=lambda item: (item.source, item.kind, item.original),
        )

    def resolve_values(self, values: Iterable[tuple[str, str, str]]) -> IdentifierResolution | None:
        materialized = list(values)
        if not materialized:
            return None
        direct = next((value for value, kind, _ in materialized if kind == "protein_id"), None)
        if direct is not None:
            resolution = self.resolve(direct)
            if resolution.status is IdentifierMatchStatus.EXACT_MATCH:
                return resolution
        candidates = sorted(self._candidates(materialized))
        joined = "|".join(value for value, _, _ in materialized)
        normalized = "|".join(normalize_identifier(value) for value, _, _ in materialized)
        return IdentifierResolution(
            input_identifier=joined,
            normalized_identifier=normalized,
            status=(
                IdentifierMatchStatus.UNIQUE_ALIAS_MATCH
                if len(candidates) == 1
                else IdentifierMatchStatus.AMBIGUOUS_MATCH
                if candidates
                else IdentifierMatchStatus.NO_MATCH
            ),
            canonical_protein_id=candidates[0] if len(candidates) == 1 else None,
            candidate_protein_ids=candidates,
        )

    def resolve_gff(self, record: GeneCoordinate) -> IdentifierResolution | None:
        return self.resolve_values(self.gff_identifiers(record))

    def resolve_annotation(self, record: AnnotationRecord) -> IdentifierResolution:
        resolution = self.resolve_values(self.annotation_identifiers(record))
        assert resolution is not None
        return resolution
