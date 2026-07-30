"""Convert OrthoFinder 3 pairwise output to the formal local orthology contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

FORMAL_COLUMNS = (
    "protein_id",
    "reference_id",
    "ortholog_id",
    "reference_organism",
    "identity",
    "query_coverage",
    "subject_coverage",
    "evalue",
    "orthogroup",
    "relationship",
    "paralog_ambiguity",
    "source",
    "source_record_id",
)

AUDIT_COLUMNS = (
    "protein_id",
    "normalized_protein_id",
    "reference_id",
    "reference_assembly",
    "ortholog_id",
    "normalized_ortholog_id",
    "orthogroup",
    "relationship",
    "query_copy_count",
    "reference_copy_count",
    "paralog_ambiguity",
    "decision",
    "decision_reason",
    "source_row",
    "source_name",
    "source_version",
    "source_command",
    "rule_version",
    "source_path",
    "source_sha256",
    "formal_sha256",
)

COVERAGE_COLUMNS = ("metric", "count", "percent_of_query_proteome")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _split_ids(value: str | None) -> list[str]:
    return sorted({item.strip() for item in (value or "").split(",") if item.strip()})


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def _relationship(query_count: int, reference_count: int) -> str:
    if query_count == 1 and reference_count == 1:
        return "one_to_one"
    if query_count == 1:
        return "one_to_many"
    if reference_count == 1:
        return "many_to_one"
    return "many_to_many"


def convert_orthofinder(
    *,
    orthologues_path: Path,
    orthogroups_path: Path,
    mapping_path: Path,
    panel_path: Path,
    formal_output: Path,
    audit_output: Path,
    coverage_output: Path,
    metadata_output: Path,
    query_species_id: str,
    query_assembly: str,
    query_protein_id: str,
    source_version: str,
    source_command: str,
) -> dict[str, Any]:
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        mapping_rows = list(csv.DictReader(handle, delimiter="\t"))
    normalized_mapping = {row["normalized_protein_id"]: row for row in mapping_rows}
    query_mapping = {
        row["normalized_protein_id"]: row
        for row in mapping_rows
        if row["species_id"] == query_species_id and row["assembly_accession"] == query_assembly
    }
    query_raw_ids = {row["raw_protein_id"] for row in query_mapping.values()}
    if not query_mapping:
        raise ValueError("No query-species IDs in normalization mapping")
    with panel_path.open(encoding="utf-8", newline="") as handle:
        panel_rows = list(csv.DictReader(handle, delimiter="\t"))
    panel_by_species = {
        row["species_id"]: row for row in panel_rows if row.get("selection_status") == "selected"
    }
    reference_species = sorted(set(panel_by_species) - {query_species_id})
    source_sha = _sha256(orthologues_path)
    rule_version = "ma4115-orthofinder-orthology-v1"
    formal_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    seen: dict[tuple[str, str, str], tuple[str, str, bool]] = {}
    exact_duplicates = 0
    unknown_query_ids: set[str] = set()
    unknown_reference_ids: set[str] = set()
    malformed = 0
    with orthologues_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"Orthogroup", "Species", query_species_id, "Orthologs"}
        missing = expected - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "OrthoFinder orthologues missing columns: " + ", ".join(sorted(missing))
            )
        for source_row, row in enumerate(reader, 2):
            orthogroup = (row.get("Orthogroup") or "").strip()
            species_id = (row.get("Species") or "").strip()
            query_ids = _split_ids(row.get(query_species_id))
            ortholog_ids = _split_ids(row.get("Orthologs"))
            if not orthogroup or not species_id or not query_ids or not ortholog_ids:
                malformed += 1
                raise ValueError(f"Malformed OrthoFinder orthologue row {source_row}")
            if species_id not in panel_by_species or species_id == query_species_id:
                raise ValueError(
                    f"Unknown or invalid reference species on row {source_row}: {species_id}"
                )
            relationship = _relationship(len(query_ids), len(ortholog_ids))
            ambiguity = len(query_ids) > 1 or len(ortholog_ids) > 1
            for normalized_query in query_ids:
                query_record = query_mapping.get(normalized_query)
                if query_record is None:
                    unknown_query_ids.add(normalized_query)
                    continue
                for normalized_ortholog in ortholog_ids:
                    ortholog_record = normalized_mapping.get(normalized_ortholog)
                    if ortholog_record is None or ortholog_record["species_id"] != species_id:
                        unknown_reference_ids.add(normalized_ortholog)
                        continue
                    identity = (
                        query_record["raw_protein_id"],
                        species_id,
                        ortholog_record["raw_protein_id"],
                    )
                    signature = (orthogroup, relationship, ambiguity)
                    previous = seen.get(identity)
                    if previous == signature:
                        exact_duplicates += 1
                        continue
                    if previous is not None:
                        raise ValueError(f"Conflicting duplicate orthology relation: {identity}")
                    seen[identity] = signature
                    source_record_id = (
                        f"{orthogroup}:{query_record['raw_protein_id']}:"
                        f"{species_id}:{ortholog_record['raw_protein_id']}"
                    )
                    formal_rows.append(
                        {
                            "protein_id": query_record["raw_protein_id"],
                            "reference_id": species_id,
                            "ortholog_id": ortholog_record["raw_protein_id"],
                            "reference_organism": panel_by_species[species_id]["organism_name"],
                            "identity": "",
                            "query_coverage": "",
                            "subject_coverage": "",
                            "evalue": "",
                            "orthogroup": orthogroup,
                            "relationship": relationship,
                            "paralog_ambiguity": str(ambiguity).lower(),
                            "source": f"OrthoFinder {source_version}",
                            "source_record_id": source_record_id,
                        }
                    )
                    audit_rows.append(
                        {
                            "protein_id": query_record["raw_protein_id"],
                            "normalized_protein_id": normalized_query,
                            "reference_id": species_id,
                            "reference_assembly": panel_by_species[species_id][
                                "assembly_accession"
                            ],
                            "ortholog_id": ortholog_record["raw_protein_id"],
                            "normalized_ortholog_id": normalized_ortholog,
                            "orthogroup": orthogroup,
                            "relationship": relationship,
                            "query_copy_count": len(query_ids),
                            "reference_copy_count": len(ortholog_ids),
                            "paralog_ambiguity": str(ambiguity).lower(),
                            "decision": "accepted",
                            "decision_reason": "official_pairwise_orthologue_output",
                            "source_row": source_row,
                            "source_name": "OrthoFinder",
                            "source_version": source_version,
                            "source_command": source_command,
                            "rule_version": rule_version,
                            "source_path": str(orthologues_path.resolve()),
                            "source_sha256": source_sha,
                            "formal_sha256": "",
                        }
                    )
    if unknown_query_ids or unknown_reference_ids:
        raise ValueError(
            "ID round-trip failed: "
            f"query={len(unknown_query_ids)}, reference={len(unknown_reference_ids)}"
        )
    formal_rows.sort(
        key=lambda row: (
            row["protein_id"],
            row["reference_id"],
            row["orthogroup"],
            row["ortholog_id"],
        )
    )
    _write_tsv(formal_output, FORMAL_COLUMNS, formal_rows)
    formal_sha = _sha256(formal_output)
    for row in audit_rows:
        row["formal_sha256"] = formal_sha
    audit_rows.sort(
        key=lambda row: (
            row["protein_id"],
            row["reference_id"],
            row["orthogroup"],
            row["ortholog_id"],
        )
    )
    _write_tsv(audit_output, AUDIT_COLUMNS, audit_rows)
    orthogroup_by_query: dict[str, str] = {}
    with orthogroups_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if query_species_id not in (reader.fieldnames or ()):
            raise ValueError(f"Orthogroups table missing query species: {query_species_id}")
        for row in reader:
            for normalized_id in _split_ids(row.get(query_species_id)):
                mapping = query_mapping.get(normalized_id)
                if mapping is None:
                    raise ValueError(f"Unknown query ID in orthogroups: {normalized_id}")
                raw_id = mapping["raw_protein_id"]
                if raw_id in orthogroup_by_query:
                    raise ValueError(f"Query protein assigned to multiple orthogroups: {raw_id}")
                orthogroup_by_query[raw_id] = row["Orthogroup"]
    relationships: dict[str, set[str]] = defaultdict(set)
    external_species: dict[str, set[str]] = defaultdict(set)
    ambiguous: set[str] = set()
    for row in formal_rows:
        protein_id = row["protein_id"]
        relationships[protein_id].add(row["relationship"])
        external_species[protein_id].add(row["reference_id"])
        if row["paralog_ambiguity"] == "true":
            ambiguous.add(protein_id)
    total = len(query_raw_ids)
    metrics = {
        "query_proteins_total": total,
        "proteins_assigned_to_orthogroup": len(orthogroup_by_query),
        "proteins_unassigned": total - len(orthogroup_by_query),
        "proteins_with_at_least_one_external_ortholog": len(external_species),
        "proteins_with_unique_one_to_one_ortholog": sum(
            "one_to_one" in values for values in relationships.values()
        ),
        "proteins_with_one_to_many_relation": sum(
            "one_to_many" in values for values in relationships.values()
        ),
        "proteins_with_many_to_one_relation": sum(
            "many_to_one" in values for values in relationships.values()
        ),
        "proteins_with_many_to_many_relation": sum(
            "many_to_many" in values for values in relationships.values()
        ),
        "proteins_with_paralog_ambiguity": len(ambiguous),
        "proteins_with_fragment_only_support": 0,
        "proteins_with_no_external_ortholog": total - len(external_species),
        "unknown_query_ids": 0,
        "unknown_comparison_ids": 0,
        "malformed_records": malformed,
        "exact_duplicate_records": exact_duplicates,
        "formal_records": len(formal_rows),
        "comparison_species_total": len(reference_species),
    }
    coverage_rows = [
        {
            "metric": metric,
            "count": count,
            "percent_of_query_proteome": f"{100 * count / total:.6f}"
            if metric.startswith("proteins_") or metric.startswith("query_proteins")
            else "",
        }
        for metric, count in metrics.items()
    ]
    _write_tsv(coverage_output, COVERAGE_COLUMNS, coverage_rows)
    if query_protein_id not in query_raw_ids:
        raise ValueError(f"Query protein absent from query assembly: {query_protein_id}")
    query_records = [row for row in formal_rows if row["protein_id"] == query_protein_id]
    query_relationships = Counter(row["relationship"] for row in query_records)
    metadata = {
        "schema_version": "1.0",
        "rule_version": rule_version,
        "source_name": "OrthoFinder",
        "source_version": source_version,
        "source_command": source_command,
        "source_path": str(orthologues_path.resolve()),
        "source_sha256": source_sha,
        "orthogroups_sha256": _sha256(orthogroups_path),
        "mapping_sha256": _sha256(mapping_path),
        "panel_sha256": _sha256(panel_path),
        "formal_sha256": formal_sha,
        "audit_sha256": _sha256(audit_output),
        "coverage_sha256": _sha256(coverage_output),
        "metrics": metrics,
        "query": {
            "protein_id": query_protein_id,
            "orthogroup": orthogroup_by_query.get(query_protein_id),
            "ortholog_record_count": len(query_records),
            "ortholog_species_count": len({row["reference_id"] for row in query_records}),
            "relationships": dict(sorted(query_relationships.items())),
            "paralog_ambiguity": any(row["paralog_ambiguity"] == "true" for row in query_records),
        },
    }
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthologues", required=True, type=Path)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--formal-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--query-species-id", required=True)
    parser.add_argument("--query-assembly", required=True)
    parser.add_argument("--query-protein-id", required=True)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--source-command", required=True)
    args = parser.parse_args()
    metadata = convert_orthofinder(
        orthologues_path=args.orthologues,
        orthogroups_path=args.orthogroups,
        mapping_path=args.mapping,
        panel_path=args.panel,
        formal_output=args.formal_output,
        audit_output=args.audit_output,
        coverage_output=args.coverage_output,
        metadata_output=args.metadata_output,
        query_species_id=args.query_species_id,
        query_assembly=args.query_assembly,
        query_protein_id=args.query_protein_id,
        source_version=args.source_version,
        source_command=args.source_command,
    )
    print(json.dumps(metadata, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
