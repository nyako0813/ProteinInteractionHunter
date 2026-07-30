#!/usr/bin/env python3
"""Audit eggNOG/InterPro GO and PFAM overlap without creating evidence."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path

from protein_interaction_hunter.exceptions import InputValidationError

if __package__:
    from scripts.convert_eggnog_annotations import AUDIT_COLUMNS as EGGNOG_COLUMNS
    from scripts.functional_annotation_common import ancestor_paths, parse_go_obo
else:
    from convert_eggnog_annotations import (  # type: ignore[import-not-found,no-redef]
        AUDIT_COLUMNS as EGGNOG_COLUMNS,
    )
    from functional_annotation_common import (  # type: ignore[import-not-found,no-redef]
        ancestor_paths,
        parse_go_obo,
    )

UNION_COLUMNS = (
    "protein_id",
    "eggnog_only_GO",
    "interpro_only_GO",
    "shared_GO",
    "GO_aspects",
    "exact_agreement_count",
    "ancestor_descendant_agreement",
    "conflict_like_difference",
    "unresolved_GO",
    "eggnog_only_PFAM",
    "interpro_only_PFAM",
    "shared_PFAM",
)


def _read_eggnog(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != EGGNOG_COLUMNS:
            raise InputValidationError("Unexpected eggNOG audit header")
        rows = {row["protein_id"]: dict(row) for row in reader}
    return rows


def _split(value: str, delimiter: str = ",") -> set[str]:
    return {item.strip() for item in value.split(delimiter) if item.strip()}


def audit_union(
    *,
    eggnog_audit: Path,
    interpro_audit: Path,
    go_obo: Path,
) -> tuple[tuple[dict[str, object], ...], tuple[tuple[str, int], ...]]:
    eggnog = _read_eggnog(eggnog_audit)
    ontology = parse_go_obo(go_obo)
    interpro_go: dict[str, set[str]] = defaultdict(set)
    interpro_pfam: dict[str, set[str]] = defaultdict(set)
    with interpro_audit.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "analysis", "signature_accession", "go_terms"}
        if not required <= set(reader.fieldnames or ()):
            raise InputValidationError("InterPro audit lacks union fields")
        for row in reader:
            protein_id = row["protein_id"]
            interpro_go[protein_id].update(
                raw.split("(", 1)[0] for raw in _split(row["go_terms"], "|")
            )
            if row["analysis"].casefold() == "pfam" and row["signature_accession"]:
                interpro_pfam[protein_id].add(row["signature_accession"])
    unknown_ids = (set(interpro_go) | set(interpro_pfam)) - set(eggnog)
    if unknown_ids:
        raise InputValidationError(
            "InterPro union contains unknown IDs: " + ", ".join(sorted(unknown_ids))
        )

    rows: list[dict[str, object]] = []
    for protein_id in sorted(eggnog):
        egg_go = _split(eggnog[protein_id]["GO_terms"])
        ipr_go = interpro_go.get(protein_id, set())
        egg_pfam = _split(eggnog[protein_id]["PFAMs"])
        ipr_pfam = interpro_pfam.get(protein_id, set())
        egg_only_go = egg_go - ipr_go
        ipr_only_go = ipr_go - egg_go
        unresolved = sorted(
            go_id
            for go_id in egg_go | ipr_go
            if go_id not in ontology.terms or ontology.terms[go_id].obsolete
        )
        ancestor_pairs: set[str] = set()
        for egg_id in sorted(egg_only_go - set(unresolved)):
            egg_ancestors = ancestor_paths(ontology, egg_id, frozenset({"is_a"}))
            for ipr_id in sorted(ipr_only_go - set(unresolved)):
                ipr_ancestors = ancestor_paths(ontology, ipr_id, frozenset({"is_a"}))
                if ipr_id in egg_ancestors:
                    ancestor_pairs.add(f"{egg_id}>{ipr_id}")
                elif egg_id in ipr_ancestors:
                    ancestor_pairs.add(f"{ipr_id}>{egg_id}")
        aspects = sorted(
            {
                ontology.terms[go_id].namespace
                for go_id in egg_go | ipr_go
                if go_id in ontology.terms
            }
        )
        rows.append(
            {
                "protein_id": protein_id,
                "eggnog_only_GO": "|".join(sorted(egg_only_go)),
                "interpro_only_GO": "|".join(sorted(ipr_only_go)),
                "shared_GO": "|".join(sorted(egg_go & ipr_go)),
                "GO_aspects": "|".join(aspects),
                "exact_agreement_count": len(egg_go & ipr_go),
                "ancestor_descendant_agreement": "|".join(sorted(ancestor_pairs)),
                "conflict_like_difference": False,
                "unresolved_GO": "|".join(unresolved),
                "eggnog_only_PFAM": "|".join(sorted(egg_pfam - ipr_pfam)),
                "interpro_only_PFAM": "|".join(sorted(ipr_pfam - egg_pfam)),
                "shared_PFAM": "|".join(sorted(egg_pfam & ipr_pfam)),
            }
        )
    metrics = Counter(
        {
            "proteome_total": len(rows),
            "proteins_with_eggnog_go": sum(
                bool(row["eggnog_only_GO"] or row["shared_GO"]) for row in rows
            ),
            "proteins_with_interpro_go": sum(
                bool(row["interpro_only_GO"] or row["shared_GO"]) for row in rows
            ),
            "proteins_with_shared_go": sum(bool(row["shared_GO"]) for row in rows),
            "proteins_with_ancestor_descendant_go": sum(
                bool(row["ancestor_descendant_agreement"]) for row in rows
            ),
            "proteins_with_unresolved_go": sum(bool(row["unresolved_GO"]) for row in rows),
            "proteins_with_eggnog_pfam": sum(
                bool(row["eggnog_only_PFAM"] or row["shared_PFAM"]) for row in rows
            ),
            "proteins_with_interpro_pfam": sum(
                bool(row["interpro_only_PFAM"] or row["shared_PFAM"]) for row in rows
            ),
            "proteins_with_shared_pfam": sum(bool(row["shared_PFAM"]) for row in rows),
            "conflict_like_differences": 0,
        }
    )
    return tuple(rows), tuple(sorted(metrics.items()))


def write_union(
    rows: Sequence[dict[str, object]],
    metrics: Sequence[tuple[str, int]],
    *,
    audit_output: Path,
    summary_output: Path,
) -> None:
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    with audit_output.open("w", encoding="utf-8", newline="\n") as handle:
        audit_writer = csv.DictWriter(
            handle,
            fieldnames=UNION_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        audit_writer.writeheader()
        audit_writer.writerows(rows)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    with summary_output.open("w", encoding="utf-8", newline="\n") as handle:
        summary_writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        summary_writer.writerow(("metric", "count"))
        summary_writer.writerows(metrics)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eggnog-audit", required=True, type=Path)
    parser.add_argument("--interpro-audit", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    rows, metrics = audit_union(
        eggnog_audit=args.eggnog_audit,
        interpro_audit=args.interpro_audit,
        go_obo=args.go_obo,
    )
    write_union(
        rows,
        metrics,
        audit_output=args.audit_output,
        summary_output=args.summary_output,
    )
    print(f"Functional source union rows: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
