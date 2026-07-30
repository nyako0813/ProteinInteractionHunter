#!/usr/bin/env python3
"""Build a formal annotation table from versioned, source-specific mappings.

Raw source terms and mapping decisions remain in a separate audit. Description
keywords never create formal categories, GO Cellular Component is excluded,
and only explicitly accepted mappings can populate ``functional_category``.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from protein_interaction_hunter.adapters.local.annotation import (
    ANNOTATION_COLUMNS,
    LocalAnnotationTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError

if __package__:
    from scripts.convert_eggnog_annotations import AUDIT_COLUMNS as EGGNOG_COLUMNS
    from scripts.functional_annotation_common import (
        GO_ROOTS,
        GoOntology,
        ancestor_paths,
        parse_go_obo,
        sha256_file,
    )
else:
    from convert_eggnog_annotations import (  # type: ignore[import-not-found,no-redef]
        AUDIT_COLUMNS as EGGNOG_COLUMNS,
    )
    from functional_annotation_common import (  # type: ignore[import-not-found,no-redef]
        GO_ROOTS,
        GoOntology,
        ancestor_paths,
        parse_go_obo,
        sha256_file,
    )

MAPPING_AUDIT_COLUMNS = (
    "protein_id",
    "source",
    "source_identifier",
    "source_label",
    "source_aspect",
    "mapping_id",
    "formal_category",
    "match_type",
    "ancestor_path",
    "mapping_status",
    "mapping_confidence",
    "conflict",
    "exclusion_reason",
    "rule_version",
)
CATEGORY_FIELDS = {
    "category_id",
    "label",
    "definition",
    "inclusion_criteria",
    "exclusion_criteria",
    "allowed_source_types",
    "accepted_go_aspects",
    "accepted_go_ancestor_relations",
    "taxonomic_applicability",
    "evidence_quality",
    "ambiguity_policy",
    "conflict_policy",
}
MAPPING_FIELDS = {
    "mapping_id",
    "formal_category",
    "source_type",
    "source_identifier",
    "source_aspect",
    "match_type",
    "allowed_relations",
    "minimum_specificity",
    "priority",
    "mapping_status",
    "rationale",
    "reference",
    "rule_version",
}
ALLOWED_RELATIONS = frozenset({"is_a", "part_of"})
ALLOWED_MAPPING_STATUSES = frozenset({"accepted", "manual_review", "rejected"})
ALLOWED_MATCH_TYPES = frozenset({"exact", "ancestor"})


@dataclass(frozen=True)
class MappingRule:
    mapping_id: str
    formal_category: str
    source_type: str
    source_identifier: str
    source_aspect: str
    match_type: str
    allowed_relations: frozenset[str]
    priority: int
    mapping_status: str
    rationale: str
    references: tuple[str, ...]
    rule_version: str


@dataclass(frozen=True)
class MappingPolicy:
    rule_version: str
    categories: frozenset[str]
    rules: tuple[MappingRule, ...]


@dataclass(frozen=True)
class SourceTerm:
    protein_id: str
    source: str
    source_type: str
    identifier: str
    label: str
    aspect: str


@dataclass(frozen=True)
class FunctionalBuild:
    annotation_rows: tuple[dict[str, str], ...]
    audit_rows: tuple[dict[str, object], ...]
    coverage_rows: tuple[tuple[str, int, str], ...]
    metadata: dict[str, object]


def _require_text(value: object, field: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise InputValidationError(f"Mapping policy has empty {field}")
    return text


def load_mapping_policy(path: Path, ontology: GoOntology) -> MappingPolicy:
    if not path.is_file():
        raise InputValidationError(f"Functional mapping policy not found: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise InputValidationError("Functional mapping policy must be a mapping")
    rule_version = _require_text(raw.get("mapping_rule_version"), "mapping_rule_version")
    raw_categories = raw.get("formal_categories")
    raw_mappings = raw.get("mappings")
    if not isinstance(raw_categories, list) or not raw_categories:
        raise InputValidationError("Mapping policy needs formal_categories")
    if not isinstance(raw_mappings, list):
        raise InputValidationError("Mapping policy mappings must be a list")

    categories: set[str] = set()
    for item in raw_categories:
        if not isinstance(item, dict) or set(item) != CATEGORY_FIELDS:
            raise InputValidationError("Formal category fields do not match the schema")
        category_id = _require_text(item["category_id"], "category_id")
        if category_id in categories:
            raise InputValidationError(f"Duplicate formal category: {category_id}")
        if not item["allowed_source_types"] or not item["accepted_go_aspects"]:
            raise InputValidationError(f"Formal category {category_id} has empty policy lists")
        categories.add(category_id)

    rules: list[MappingRule] = []
    seen_mapping_ids: set[str] = set()
    source_targets: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    source_priorities: set[tuple[str, str, int]] = set()
    for item in raw_mappings:
        if not isinstance(item, dict) or set(item) != MAPPING_FIELDS:
            raise InputValidationError("Mapping record fields do not match the schema")
        mapping_id = _require_text(item["mapping_id"], "mapping_id")
        if mapping_id in seen_mapping_ids:
            raise InputValidationError(f"Duplicate mapping_id: {mapping_id}")
        seen_mapping_ids.add(mapping_id)
        category = _require_text(item["formal_category"], "formal_category")
        if category not in categories:
            raise InputValidationError(f"Unknown formal category: {category}")
        source_type = _require_text(item["source_type"], "source_type")
        source_identifier = _require_text(item["source_identifier"], "source_identifier")
        source_aspect = _require_text(item["source_aspect"], "source_aspect")
        match_type = _require_text(item["match_type"], "match_type")
        status = _require_text(item["mapping_status"], "mapping_status")
        if match_type not in ALLOWED_MATCH_TYPES:
            raise InputValidationError(f"Invalid match_type: {match_type}")
        if status not in ALLOWED_MAPPING_STATUSES:
            raise InputValidationError(f"Invalid mapping_status: {status}")
        relation_values = frozenset(str(value) for value in item["allowed_relations"])
        if not relation_values <= ALLOWED_RELATIONS:
            raise InputValidationError(f"Invalid GO ancestor relation in mapping {mapping_id}")
        if match_type == "ancestor" and not relation_values:
            raise InputValidationError(f"Ancestor mapping {mapping_id} has no relation")
        if match_type == "exact" and relation_values:
            raise InputValidationError(f"Exact mapping {mapping_id} has ancestor relations")
        if source_type == "go":
            if source_identifier not in ontology.terms:
                raise InputValidationError(f"Unreachable GO ID: {source_identifier}")
            term = ontology.terms[source_identifier]
            if term.obsolete:
                raise InputValidationError(f"Obsolete GO mapping: {source_identifier}")
            if source_identifier in GO_ROOTS:
                raise InputValidationError(f"Overbroad GO root mapping: {source_identifier}")
            if source_aspect != term.namespace:
                raise InputValidationError(
                    f"GO aspect mismatch for {source_identifier}: {source_aspect}"
                )
        references = item["reference"]
        if (
            not isinstance(references, list)
            or not references
            or not all(str(value).strip() for value in references)
        ):
            raise InputValidationError(f"Mapping {mapping_id} has no reference")
        rationale = _require_text(item["rationale"], "rationale")
        if _require_text(item["rule_version"], "rule_version") != rule_version:
            raise InputValidationError(f"Mapping {mapping_id} has wrong rule version")
        try:
            priority = int(item["priority"])
        except (TypeError, ValueError) as exc:
            raise InputValidationError(f"Mapping {mapping_id} has invalid priority") from exc
        priority_key = (source_type, source_identifier, priority)
        if priority_key in source_priorities:
            raise InputValidationError(
                f"Priority tie for {source_type}:{source_identifier}:{priority}"
            )
        source_priorities.add(priority_key)
        if status == "accepted":
            source_targets[(source_type, source_identifier, match_type)].add(category)
        rules.append(
            MappingRule(
                mapping_id=mapping_id,
                formal_category=category,
                source_type=source_type,
                source_identifier=source_identifier,
                source_aspect=source_aspect,
                match_type=match_type,
                allowed_relations=relation_values,
                priority=priority,
                mapping_status=status,
                rationale=rationale,
                references=tuple(str(value) for value in references),
                rule_version=rule_version,
            )
        )
    conflicts = {key: targets for key, targets in source_targets.items() if len(targets) > 1}
    if conflicts:
        raise InputValidationError(f"Conflicting accepted mappings: {conflicts}")
    return MappingPolicy(
        rule_version=rule_version,
        categories=frozenset(categories),
        rules=tuple(sorted(rules, key=lambda rule: (-rule.priority, rule.mapping_id))),
    )


def _read_exact_tsv(path: Path, columns: Sequence[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise InputValidationError(f"TSV not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise InputValidationError(f"Unexpected TSV header in {path}")
        return [dict(row) for row in reader]


def _split(value: str, delimiter: str = ",") -> tuple[str, ...]:
    return tuple(sorted({item.strip() for item in value.split(delimiter) if item.strip()}))


def collect_source_terms(
    eggnog_rows: Sequence[dict[str, str]],
    interpro_audit: Path | None,
    ontology: GoOntology,
) -> tuple[SourceTerm, ...]:
    terms: set[SourceTerm] = set()
    for row in eggnog_rows:
        protein_id = row["protein_id"]
        for go_id in _split(row["GO_terms"]):
            term = ontology.terms.get(go_id)
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="go",
                    identifier=go_id,
                    label=term.name if term else "",
                    aspect=term.namespace if term else "unknown",
                )
            )
        for ko_id in _split(row["KEGG_ko"]):
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="kegg_ko",
                    identifier=ko_id,
                    label="",
                    aspect="not_applicable",
                )
            )
        for ec_id in _split(row["EC"]):
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="ec",
                    identifier=ec_id,
                    label="",
                    aspect="not_applicable",
                )
            )
        for og_id in _split(row["eggNOG_OGs"]):
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="eggnog_og",
                    identifier=og_id,
                    label="",
                    aspect="not_applicable",
                )
            )
        if row["COG_category"]:
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="cog_category",
                    identifier=row["COG_category"],
                    label="",
                    aspect="not_applicable",
                )
            )
        if row["description"]:
            terms.add(
                SourceTerm(
                    protein_id=protein_id,
                    source="eggNOG-mapper",
                    source_type="description",
                    identifier=row["description"],
                    label=row["description"],
                    aspect="not_applicable",
                )
            )
    if interpro_audit:
        with interpro_audit.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"protein_id", "go_terms", "interpro_accession", "interpro_description"}
            if not required <= set(reader.fieldnames or ()):
                raise InputValidationError("InterPro audit lacks required union columns")
            for row in reader:
                for raw_go in _split(row.get("go_terms") or "", "|"):
                    go_id = raw_go.split("(", 1)[0]
                    term = ontology.terms.get(go_id)
                    terms.add(
                        SourceTerm(
                            protein_id=row["protein_id"],
                            source="InterProScan",
                            source_type="go",
                            identifier=go_id,
                            label=term.name if term else "",
                            aspect=term.namespace if term else "unknown",
                        )
                    )
    return tuple(
        sorted(
            terms,
            key=lambda term: (
                term.protein_id,
                term.source,
                term.source_type,
                term.identifier,
            ),
        )
    )


def _audit_base(term: SourceTerm, policy: MappingPolicy) -> dict[str, object]:
    return {
        "protein_id": term.protein_id,
        "source": term.source,
        "source_identifier": term.identifier,
        "source_label": term.label,
        "source_aspect": term.aspect,
        "mapping_id": "",
        "formal_category": "",
        "match_type": "",
        "ancestor_path": "",
        "mapping_status": "unmapped",
        "mapping_confidence": "",
        "conflict": False,
        "exclusion_reason": "no_accepted_mapping",
        "rule_version": policy.rule_version,
    }


def map_source_term(
    term: SourceTerm,
    policy: MappingPolicy,
    ontology: GoOntology,
) -> list[dict[str, object]]:
    base = _audit_base(term, policy)
    if term.source_type == "description":
        base.update(
            {
                "mapping_status": "manual_review",
                "exclusion_reason": "description_not_mapping_evidence",
            }
        )
        return [base]
    if term.source_type == "go":
        go_term = ontology.terms.get(term.identifier)
        if go_term is None:
            base["exclusion_reason"] = "unknown_go_id"
            return [base]
        if go_term.obsolete:
            replacement = ",".join(go_term.replaced_by)
            consider = ",".join(go_term.consider)
            reason = "obsolete_go_id"
            if replacement:
                reason += f";replaced_by={replacement}"
            if consider:
                reason += f";consider={consider}"
            base["exclusion_reason"] = reason
            return [base]
        if go_term.namespace == "cellular_component":
            base.update(
                {
                    "mapping_status": "excluded",
                    "exclusion_reason": "go_cellular_component_excluded",
                }
            )
            return [base]
    matches: list[tuple[MappingRule, str, str]] = []
    for rule in policy.rules:
        if rule.mapping_status != "accepted" or rule.source_type != term.source_type:
            continue
        if rule.match_type == "exact" and rule.source_identifier == term.identifier:
            matches.append((rule, "exact", term.identifier))
        elif rule.match_type == "ancestor" and term.source_type == "go":
            paths = ancestor_paths(ontology, term.identifier, rule.allowed_relations)
            if rule.source_identifier in paths:
                encoded = ";".join(">".join(path) for path in paths[rule.source_identifier])
                matches.append((rule, "ancestor", encoded))
    if not matches:
        return [base]
    highest_priority = max(rule.priority for rule, _, _ in matches)
    selected = [item for item in matches if item[0].priority == highest_priority]
    categories = {rule.formal_category for rule, _, _ in selected}
    conflict = len(categories) > 1
    rows: list[dict[str, object]] = []
    for rule, match_type, path in selected:
        row = dict(base)
        row.update(
            {
                "mapping_id": rule.mapping_id,
                "formal_category": rule.formal_category,
                "match_type": match_type,
                "ancestor_path": path,
                "mapping_status": "conflict" if conflict else "accepted",
                "mapping_confidence": "high" if match_type == "exact" else "medium",
                "conflict": conflict,
                "exclusion_reason": "conflicting_top_priority_mapping" if conflict else "",
            }
        )
        rows.append(row)
    return rows


def build_functional_table(
    *,
    annotation_table: Path,
    eggnog_audit: Path,
    go_obo: Path,
    mapping_policy: Path,
    query_id: str,
    interpro_audit: Path | None = None,
) -> FunctionalBuild:
    LocalAnnotationTsvLoader().load(annotation_table)
    annotation_rows = _read_exact_tsv(annotation_table, ANNOTATION_COLUMNS)
    eggnog_rows = _read_exact_tsv(eggnog_audit, EGGNOG_COLUMNS)
    annotation_ids = {row["protein_id"] for row in annotation_rows}
    eggnog_ids = {row["protein_id"] for row in eggnog_rows}
    unknown_annotation_ids = sorted(annotation_ids - eggnog_ids)
    if unknown_annotation_ids:
        raise InputValidationError(
            "Annotation IDs absent from eggNOG audit: " + ", ".join(unknown_annotation_ids)
        )
    if query_id not in annotation_ids:
        raise InputValidationError(f"Query {query_id} is absent from annotations")
    ontology = parse_go_obo(go_obo)
    policy = load_mapping_policy(mapping_policy, ontology)
    source_terms = collect_source_terms(eggnog_rows, interpro_audit, ontology)
    unknown_term_ids = sorted({term.protein_id for term in source_terms} - eggnog_ids)
    if unknown_term_ids:
        raise InputValidationError(
            "Functional sources contain unknown protein IDs: " + ", ".join(unknown_term_ids)
        )
    audit_rows = tuple(
        row for term in source_terms for row in map_source_term(term, policy, ontology)
    )
    categories_by_protein: dict[str, set[str]] = defaultdict(set)
    for audit_row in audit_rows:
        if audit_row["mapping_status"] == "accepted":
            categories_by_protein[str(audit_row["protein_id"])].add(
                str(audit_row["formal_category"])
            )
    enhanced_rows: list[dict[str, str]] = []
    for annotation_row in annotation_rows:
        enhanced = dict(annotation_row)
        enhanced["functional_category"] = "|".join(
            sorted(categories_by_protein.get(annotation_row["protein_id"], set()))
        )
        enhanced_rows.append(enhanced)

    total = len(eggnog_ids)
    source_ids = {term.protein_id for term in source_terms}
    accepted_ids = set(categories_by_protein)
    status_by_protein: dict[str, set[str]] = defaultdict(set)
    for audit_row in audit_rows:
        status_by_protein[str(audit_row["protein_id"])].add(str(audit_row["mapping_status"]))
    metrics: list[tuple[str, int]] = [
        ("proteome_total", total),
        ("proteins_with_formal_category", len(accepted_ids)),
        (
            "proteins_with_multiple_formal_categories",
            sum(len(values) > 1 for values in categories_by_protein.values()),
        ),
        (
            "proteins_with_exact_mapping",
            len(
                {
                    str(row["protein_id"])
                    for row in audit_rows
                    if row["mapping_status"] == "accepted" and row["match_type"] == "exact"
                }
            ),
        ),
        (
            "proteins_with_ancestor_mapping",
            len(
                {
                    str(row["protein_id"])
                    for row in audit_rows
                    if row["mapping_status"] == "accepted" and row["match_type"] == "ancestor"
                }
            ),
        ),
        (
            "proteins_with_ambiguous_mapping",
            sum("ambiguous" in statuses for statuses in status_by_protein.values()),
        ),
        (
            "proteins_with_conflicting_mapping",
            sum("conflict" in statuses for statuses in status_by_protein.values()),
        ),
        (
            "proteins_with_only_unmapped_source_terms",
            sum(
                bool(statuses) and not ({"accepted", "conflict", "ambiguous"} & statuses)
                for statuses in status_by_protein.values()
            ),
        ),
        (
            "proteins_with_source_but_no_formal_category",
            len(source_ids - accepted_ids),
        ),
        ("proteins_with_no_source_annotation", len(eggnog_ids - source_ids)),
        ("formal_annotation_table_rows", len(enhanced_rows)),
        ("proteome_ids_without_formal_annotation_row", len(eggnog_ids - annotation_ids)),
    ]
    for category, count in sorted(
        Counter(
            category for values in categories_by_protein.values() for category in values
        ).items()
    ):
        metrics.append((f"category:{category}", count))
    coverage_rows = tuple(
        (
            metric,
            count,
            f"{(100.0 * count / total) if total else 0.0:.6f}",
        )
        for metric, count in metrics
    )
    query_categories = sorted(categories_by_protein.get(query_id, set()))
    metadata: dict[str, object] = {
        "mapping_rule_version": policy.rule_version,
        "annotation_table_sha256": sha256_file(annotation_table),
        "eggnog_audit_sha256": sha256_file(eggnog_audit),
        "go_data_version": ontology.data_version,
        "go_obo_sha256": sha256_file(go_obo),
        "mapping_policy_sha256": sha256_file(mapping_policy),
        "interpro_audit_sha256": (sha256_file(interpro_audit) if interpro_audit else None),
        "query_id": query_id,
        "query_formal_categories": query_categories,
        "query_mapping_status": "ready" if query_categories else "unresolved",
        "query_unresolved_reason": (
            "" if query_categories else "no_accepted_source_to_category_mapping"
        ),
        "formal_category_count": len(policy.categories),
        "accepted_mapping_rule_count": sum(
            rule.mapping_status == "accepted" for rule in policy.rules
        ),
        "complementarity_rules_status": "not_ready",
        "functional_evidence_status": "not_ready",
        "scoring_status": "not_ready",
    }
    return FunctionalBuild(
        annotation_rows=tuple(enhanced_rows),
        audit_rows=audit_rows,
        coverage_rows=coverage_rows,
        metadata=metadata,
    )


def _write_dict_rows(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[dict[str, object] | dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_functional_build(
    build: FunctionalBuild,
    *,
    annotation_output: Path,
    mapping_audit_output: Path,
    coverage_output: Path,
    metadata_output: Path,
) -> None:
    _write_dict_rows(annotation_output, ANNOTATION_COLUMNS, build.annotation_rows)
    LocalAnnotationTsvLoader().load(annotation_output)
    _write_dict_rows(mapping_audit_output, MAPPING_AUDIT_COLUMNS, build.audit_rows)
    coverage_output.parent.mkdir(parents=True, exist_ok=True)
    with coverage_output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "count", "percent_of_proteome"))
        writer.writerows(build.coverage_rows)
    metadata = dict(build.metadata)
    metadata.update(
        {
            "annotation_output_sha256": sha256_file(annotation_output),
            "mapping_audit_output_sha256": sha256_file(mapping_audit_output),
            "coverage_output_sha256": sha256_file(coverage_output),
        }
    )
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-table", required=True, type=Path)
    parser.add_argument("--eggnog-audit", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--mapping-policy", required=True, type=Path)
    parser.add_argument("--query-id", required=True)
    parser.add_argument("--interpro-audit", type=Path)
    parser.add_argument("--annotation-output", required=True, type=Path)
    parser.add_argument("--mapping-audit-output", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    build = build_functional_table(
        annotation_table=args.annotation_table,
        eggnog_audit=args.eggnog_audit,
        go_obo=args.go_obo,
        mapping_policy=args.mapping_policy,
        query_id=args.query_id,
        interpro_audit=args.interpro_audit,
    )
    write_functional_build(
        build,
        annotation_output=args.annotation_output,
        mapping_audit_output=args.mapping_audit_output,
        coverage_output=args.coverage_output,
        metadata_output=args.metadata_output,
    )
    print(f"Functional annotation rows: {len(build.annotation_rows)}")
    print(f"Mapping audit rows: {len(build.audit_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
