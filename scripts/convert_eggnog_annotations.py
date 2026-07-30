#!/usr/bin/env python3
"""Convert eggNOG-mapper v2.1 annotations into a lossless audit table.

No source annotation is promoted to a formal interaction category here.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.exceptions import InputValidationError

if __package__:
    from scripts.functional_annotation_common import parse_go_obo, sha256_file
else:
    from functional_annotation_common import (  # type: ignore[import-not-found,no-redef]
        parse_go_obo,
        sha256_file,
    )

RAW_COLUMNS = (
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)
AUDIT_COLUMNS = (
    "protein_id",
    "seed_ortholog",
    "seed_ortholog_evalue",
    "seed_ortholog_score",
    "eggNOG_OGs",
    "max_annot_level",
    "COG_category",
    "description",
    "preferred_name",
    "GO_terms",
    "EC",
    "KEGG_ko",
    "KEGG_pathway",
    "KEGG_module",
    "KEGG_reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_reaction",
    "PFAMs",
    "source",
    "source_version",
    "database_version",
    "command",
    "raw_row_number",
    "parse_status",
    "exclusion_reason",
)
RAW_TO_AUDIT = {
    "query": "protein_id",
    "seed_ortholog": "seed_ortholog",
    "evalue": "seed_ortholog_evalue",
    "score": "seed_ortholog_score",
    "eggNOG_OGs": "eggNOG_OGs",
    "max_annot_lvl": "max_annot_level",
    "COG_category": "COG_category",
    "Description": "description",
    "Preferred_name": "preferred_name",
    "GOs": "GO_terms",
    "EC": "EC",
    "KEGG_ko": "KEGG_ko",
    "KEGG_Pathway": "KEGG_pathway",
    "KEGG_Module": "KEGG_module",
    "KEGG_Reaction": "KEGG_reaction",
    "KEGG_rclass": "KEGG_rclass",
    "BRITE": "BRITE",
    "KEGG_TC": "KEGG_TC",
    "CAZy": "CAZy",
    "BiGG_Reaction": "BiGG_reaction",
    "PFAMs": "PFAMs",
}
UNKNOWN_DESCRIPTIONS = ("unknown function", "function unknown", "uncharacterized")
BROAD_DESCRIPTIONS = ("protein of unknown function", "hypothetical protein")


@dataclass(frozen=True)
class EggnogRow:
    values: tuple[str, ...]
    raw_row_number: int

    @property
    def protein_id(self) -> str:
        return self.values[0]


@dataclass(frozen=True)
class EggnogConversion:
    audit_rows: tuple[dict[str, object], ...]
    coverage_rows: tuple[tuple[str, int, str], ...]
    metadata: dict[str, object]


def optional(value: str) -> str:
    stripped = value.strip()
    return "" if stripped == "-" else stripped


def percent(count: int, total: int) -> str:
    return f"{(100.0 * count / total) if total else 0.0:.6f}"


def parse_eggnog_annotations(path: Path) -> tuple[list[EggnogRow], list[str], int]:
    """Validate the exact v2.1 schema and collapse only identical duplicates."""
    if not path.is_file():
        raise InputValidationError(f"eggNOG annotations not found: {path}")
    metadata_lines: list[str] = []
    header_seen = False
    rows_by_id: dict[str, EggnogRow] = {}
    identical_duplicates = 0
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##"):
                metadata_lines.append(line)
                continue
            if line.startswith("#"):
                parsed_header = tuple(line[1:].split("\t"))
                if header_seen:
                    raise InputValidationError("eggNOG output contains multiple headers")
                if parsed_header != RAW_COLUMNS:
                    raise InputValidationError(
                        f"eggNOG header mismatch: expected {RAW_COLUMNS!r}, found {parsed_header!r}"
                    )
                header_seen = True
                continue
            if not header_seen:
                raise InputValidationError(f"eggNOG data precedes the header on line {line_number}")
            values = tuple(optional(value) for value in line.split("\t"))
            if len(values) != len(RAW_COLUMNS):
                raise InputValidationError(
                    f"eggNOG row {line_number} has {len(values)} columns; "
                    f"expected {len(RAW_COLUMNS)}"
                )
            protein_id = values[0]
            if not protein_id:
                raise InputValidationError(f"eggNOG row {line_number} has an empty query ID")
            row = EggnogRow(values=values, raw_row_number=line_number)
            prior = rows_by_id.get(protein_id)
            if prior is not None:
                if prior.values == row.values:
                    identical_duplicates += 1
                    continue
                raise InputValidationError(
                    f"Conflicting duplicate eggNOG rows for query {protein_id}"
                )
            rows_by_id[protein_id] = row
    if not header_seen:
        raise InputValidationError("eggNOG output is empty or lacks its exact header")
    return list(rows_by_id.values()), metadata_lines, identical_duplicates


def go_aspects(go_terms: str, namespaces: dict[str, str]) -> set[str]:
    aspects: set[str] = set()
    for go_id in (value for value in go_terms.split(",") if value):
        namespace = namespaces.get(go_id)
        if namespace is None:
            aspects.add("unknown")
        elif namespace == "molecular_function":
            aspects.add("MF")
        elif namespace == "biological_process":
            aspects.add("BP")
        elif namespace == "cellular_component":
            aspects.add("CC")
        else:
            aspects.add("unknown")
    return aspects


def read_query_ids(path: Path | None, label: str) -> set[str]:
    if path is None:
        return set()
    if not path.is_file():
        raise InputValidationError(f"eggNOG {label} file not found: {path}")
    result: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            protein_id = line.split("	", 1)[0].strip()
            if not protein_id:
                raise InputValidationError(
                    f"eggNOG {label} has an empty query on line {line_number}"
                )
            result.add(protein_id)
    return result


def convert_eggnog(
    *,
    input_path: Path,
    proteome_fasta: Path,
    query_id: str,
    source_version: str,
    database_version: str,
    command: str,
    go_obo: Path | None = None,
    hits_path: Path | None = None,
    seed_orthologs_path: Path | None = None,
) -> EggnogConversion:
    parsed_rows, metadata_lines, identical_duplicates = parse_eggnog_annotations(input_path)
    if command == "auto":
        recorded_commands = [
            line[3:] for line in metadata_lines if "emapper.py" in line and line.startswith("## ")
        ]
        if len(recorded_commands) != 1:
            raise InputValidationError("Could not resolve exactly one command from eggNOG metadata")
        command = recorded_commands[0]
    proteins = LocalFastaLoader().load(proteome_fasta)
    protein_ids = {str(record.protein_id) for record in proteins}
    hit_ids = read_query_ids(hits_path, "hits")
    seed_ids = read_query_ids(seed_orthologs_path, "seed orthologs")
    if query_id not in protein_ids:
        raise InputValidationError(f"Query {query_id} is absent from the proteome")
    unknown_ids = sorted(
        ({row.protein_id for row in parsed_rows} | hit_ids | seed_ids) - protein_ids
    )
    if unknown_ids:
        raise InputValidationError(
            "eggNOG output contains proteome-unknown IDs: " + ", ".join(unknown_ids)
        )
    by_id = {row.protein_id: row for row in parsed_rows}
    namespaces = (
        {go_id: term.namespace for go_id, term in parse_go_obo(go_obo).terms.items()}
        if go_obo
        else {}
    )
    audit_rows: list[dict[str, object]] = []
    for protein_id in sorted(protein_ids):
        parsed = by_id.get(protein_id)
        row: dict[str, object] = dict.fromkeys(AUDIT_COLUMNS, "")
        row.update(
            {
                "protein_id": protein_id,
                "source": "eggNOG-mapper",
                "source_version": source_version,
                "database_version": database_version,
                "command": command,
            }
        )
        if parsed is None:
            if seed_orthologs_path and protein_id in seed_ids:
                row.update(
                    {
                        "parse_status": "hit_no_annotation",
                        "exclusion_reason": "seed_assigned_but_no_transfer_annotation",
                    }
                )
            elif hits_path and protein_id not in hit_ids:
                row.update(
                    {
                        "parse_status": "no_hit",
                        "exclusion_reason": "no_significant_search_hit",
                    }
                )
            else:
                row.update(
                    {
                        "parse_status": "missing",
                        "exclusion_reason": "missing_raw_annotation_row",
                    }
                )
        else:
            for raw_name, value in zip(RAW_COLUMNS, parsed.values, strict=True):
                row[RAW_TO_AUDIT[raw_name]] = value
            row["raw_row_number"] = parsed.raw_row_number
            if parsed.values[1]:
                row["parse_status"] = "annotated"
            else:
                row.update(
                    {
                        "parse_status": "no_hit",
                        "exclusion_reason": "no_seed_ortholog",
                    }
                )
        audit_rows.append(row)
    if sum(row["protein_id"] == query_id for row in audit_rows) != 1:
        raise InputValidationError(f"Query {query_id} was not preserved exactly once")

    total = len(protein_ids)
    annotated = [row for row in audit_rows if row["parse_status"] == "annotated"]
    metrics = [
        ("proteome_total", total),
        ("raw_rows", len(parsed_rows)),
        (
            "search_hits",
            len(hit_ids) if hits_path else sum(bool(row["seed_ortholog"]) for row in audit_rows),
        ),
        (
            "seed_ortholog_assigned",
            len(seed_ids)
            if seed_orthologs_path
            else sum(bool(row["seed_ortholog"]) for row in audit_rows),
        ),
        ("description_present", sum(bool(row["description"]) for row in audit_rows)),
        ("preferred_name_present", sum(bool(row["preferred_name"]) for row in audit_rows)),
        ("eggnog_og_present", sum(bool(row["eggNOG_OGs"]) for row in audit_rows)),
        ("cog_present", sum(bool(row["COG_category"]) for row in audit_rows)),
        (
            "go_mf_present",
            sum("MF" in go_aspects(str(row["GO_terms"]), namespaces) for row in audit_rows),
        ),
        (
            "go_bp_present",
            sum("BP" in go_aspects(str(row["GO_terms"]), namespaces) for row in audit_rows),
        ),
        (
            "go_cc_present",
            sum("CC" in go_aspects(str(row["GO_terms"]), namespaces) for row in audit_rows),
        ),
        (
            "go_unknown_aspect_present",
            sum("unknown" in go_aspects(str(row["GO_terms"]), namespaces) for row in audit_rows),
        ),
        (
            "go_annotation_count",
            sum(
                len([value for value in str(row["GO_terms"]).split(",") if value])
                for row in audit_rows
            ),
        ),
        ("ec_present", sum(bool(row["EC"]) for row in audit_rows)),
        ("kegg_ko_present", sum(bool(row["KEGG_ko"]) for row in audit_rows)),
        ("kegg_pathway_present", sum(bool(row["KEGG_pathway"]) for row in audit_rows)),
        ("kegg_module_present", sum(bool(row["KEGG_module"]) for row in audit_rows)),
        ("kegg_reaction_present", sum(bool(row["KEGG_reaction"]) for row in audit_rows)),
        ("pfam_present", sum(bool(row["PFAMs"]) for row in audit_rows)),
        ("no_hit", sum(row["parse_status"] == "no_hit" for row in audit_rows)),
        (
            "hit_no_annotation",
            sum(row["parse_status"] == "hit_no_annotation" for row in audit_rows),
        ),
        ("missing_raw_row", sum(row["parse_status"] == "missing" for row in audit_rows)),
        (
            "unknown_function_description",
            sum(
                any(term in str(row["description"]).casefold() for term in UNKNOWN_DESCRIPTIONS)
                for row in annotated
            ),
        ),
        (
            "broad_only_annotation",
            sum(
                any(term in str(row["description"]).casefold() for term in BROAD_DESCRIPTIONS)
                for row in annotated
            ),
        ),
        ("identical_duplicate_rows", identical_duplicates),
        ("conflicting_duplicate_rows", 0),
        ("malformed_rows", 0),
        ("unknown_ids", 0),
    ]
    coverage_rows = tuple((name, count, percent(count, total)) for name, count in metrics)
    metadata: dict[str, object] = {
        "schema": "eggnog-mapper-v2.1-annotations-21-columns",
        "source_version": source_version,
        "database_version": database_version,
        "command": command,
        "raw_source_path": str(input_path.expanduser().resolve()),
        "raw_source_sha256": sha256_file(input_path),
        "proteome_path": str(proteome_fasta.expanduser().resolve()),
        "proteome_sha256": sha256_file(proteome_fasta),
        "query_id": query_id,
        "metadata_lines": metadata_lines,
        "identical_duplicate_rows": identical_duplicates,
        "raw_row_count": len(parsed_rows),
        "audit_row_count": len(audit_rows),
        "parse_status_counts": dict(
            sorted(Counter(str(row["parse_status"]) for row in audit_rows).items())
        ),
    }
    if go_obo:
        metadata.update(
            {
                "go_obo_path": str(go_obo.expanduser().resolve()),
                "go_obo_sha256": sha256_file(go_obo),
            }
        )
    return EggnogConversion(
        audit_rows=tuple(audit_rows),
        coverage_rows=coverage_rows,
        metadata=metadata,
    )


def write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_conversion(
    result: EggnogConversion,
    *,
    audit_output: Path,
    coverage_output: Path,
    metadata_output: Path,
) -> None:
    write_tsv(
        audit_output,
        AUDIT_COLUMNS,
        tuple(tuple(row[column] for column in AUDIT_COLUMNS) for row in result.audit_rows),
    )
    write_tsv(
        coverage_output,
        ("metric", "count", "percent_of_proteome"),
        result.coverage_rows,
    )
    metadata = dict(result.metadata)
    metadata["audit_output_sha256"] = sha256_file(audit_output)
    metadata["coverage_output_sha256"] = sha256_file(coverage_output)
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input", required=True, type=Path)
    result.add_argument("--proteome-fasta", required=True, type=Path)
    result.add_argument("--query-id", required=True)
    result.add_argument("--source-version", required=True)
    result.add_argument("--database-version", required=True)
    result.add_argument("--command", required=True)
    result.add_argument("--go-obo", type=Path)
    result.add_argument("--hits", type=Path)
    result.add_argument("--seed-orthologs", type=Path)
    result.add_argument("--audit-output", required=True, type=Path)
    result.add_argument("--coverage-output", required=True, type=Path)
    result.add_argument("--metadata-output", required=True, type=Path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    result = convert_eggnog(
        input_path=args.input,
        proteome_fasta=args.proteome_fasta,
        query_id=args.query_id,
        source_version=args.source_version,
        database_version=args.database_version,
        command=args.command,
        go_obo=args.go_obo,
        hits_path=args.hits,
        seed_orthologs_path=args.seed_orthologs,
    )
    write_conversion(
        result,
        audit_output=args.audit_output,
        coverage_output=args.coverage_output,
        metadata_output=args.metadata_output,
    )
    print(f"eggNOG audit rows: {len(result.audit_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
