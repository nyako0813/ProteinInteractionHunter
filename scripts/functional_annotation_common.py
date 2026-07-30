"""Shared, dependency-light helpers for the MA_4115 functional workflow."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from protein_interaction_hunter.exceptions import InputValidationError

GO_ROOTS = frozenset({"GO:0003674", "GO:0008150", "GO:0005575"})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class GoTerm:
    go_id: str
    name: str
    namespace: str
    parents: tuple[tuple[str, str], ...]
    obsolete: bool
    replaced_by: tuple[str, ...]
    consider: tuple[str, ...]


@dataclass(frozen=True)
class GoOntology:
    data_version: str
    terms: dict[str, GoTerm]


def parse_go_obo(path: Path) -> GoOntology:
    """Parse the fields needed for deterministic mapping from a GO OBO file."""
    if not path.is_file():
        raise InputValidationError(f"GO OBO file not found: {path}")
    data_version = ""
    records: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "[Term]":
                if current is not None:
                    records.append(current)
                current = {}
                continue
            if line.startswith("[") and line.endswith("]"):
                if current is not None:
                    records.append(current)
                    current = None
                continue
            if not line or line.startswith("!"):
                continue
            if ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if current is None:
                if key == "data-version":
                    data_version = value
                continue
            current.setdefault(key, []).append(value)
    if current is not None:
        records.append(current)

    terms: dict[str, GoTerm] = {}
    for record in records:
        ids = record.get("id", [])
        if len(ids) != 1:
            raise InputValidationError("GO term must contain exactly one id")
        go_id = ids[0]
        if go_id in terms:
            raise InputValidationError(f"Duplicate GO term: {go_id}")
        parents: set[tuple[str, str]] = set()
        for value in record.get("is_a", []):
            parents.add(("is_a", value.split(" ! ", 1)[0]))
        for value in record.get("relationship", []):
            parts = value.split()
            if len(parts) < 2:
                raise InputValidationError(f"Malformed GO relationship for {go_id}")
            parents.add((parts[0], parts[1]))
        terms[go_id] = GoTerm(
            go_id=go_id,
            name=(record.get("name") or [""])[0],
            namespace=(record.get("namespace") or [""])[0],
            parents=tuple(sorted(parents)),
            obsolete=(record.get("is_obsolete") or ["false"])[0].casefold() == "true",
            replaced_by=tuple(sorted(record.get("replaced_by", []))),
            consider=tuple(sorted(record.get("consider", []))),
        )
    if not terms:
        raise InputValidationError("GO OBO file contains no [Term] records")
    return GoOntology(data_version=data_version, terms=terms)


def ancestor_paths(
    ontology: GoOntology,
    start_id: str,
    allowed_relations: frozenset[str],
) -> dict[str, tuple[tuple[str, ...], ...]]:
    """Return every deterministic acyclic path from a term to each ancestor."""
    if start_id not in ontology.terms:
        raise InputValidationError(f"Unknown GO term: {start_id}")
    paths: dict[str, set[tuple[str, ...]]] = {start_id: {(start_id,)}}

    def visit(term_id: str, path: tuple[str, ...]) -> None:
        for relation, parent_id in ontology.terms[term_id].parents:
            if relation not in allowed_relations:
                continue
            if parent_id not in ontology.terms:
                raise InputValidationError(
                    f"GO parent {parent_id} referenced by {term_id} is absent"
                )
            if parent_id in path:
                raise InputValidationError("GO cycle detected: " + " -> ".join((*path, parent_id)))
            parent_path = (*path, parent_id)
            paths.setdefault(parent_id, set()).add(parent_path)
            visit(parent_id, parent_path)

    visit(start_id, (start_id,))
    return {term_id: tuple(sorted(term_paths)) for term_id, term_paths in sorted(paths.items())}
