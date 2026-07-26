#!/usr/bin/env python3
"""Convert InterProScan 5 TSV output into ProteinInteractionHunter domain TSV.

Only fields supported by the formal domain loader are emitted there. InterPro,
GO, pathway, score, status, and run metadata are retained in a separate audit.
No biological category or interaction rule is inferred.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from protein_interaction_hunter.adapters.local.domains import DOMAIN_COLUMNS
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.exceptions import InputValidationError

ALLOWED_COLUMN_COUNTS = {11, 13, 14, 15}
AUDIT_COLUMNS = (
    "protein_id",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score_or_evalue",
    "status",
    "run_date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathway_terms",
    "architecture_index",
    "known_protein_id",
    "included_in_domain_table",
    "interproscan_version",
    "source_file",
)
_HYPOTHETICAL_TERMS = ("hypothetical", "uncharacterized")


@dataclass(frozen=True)
class InterProHit:
    protein_id: str
    sequence_md5: str
    sequence_length: int
    analysis: str
    signature_accession: str
    signature_description: str
    start: int
    end: int
    score_or_evalue: str
    status: str
    run_date: str
    interpro_accession: str = ""
    interpro_description: str = ""
    go_terms: str = ""
    pathway_terms: str = ""


@dataclass(frozen=True)
class DomainRow:
    protein_id: str
    source: str
    accession: str
    name: str
    start: int
    end: int
    architecture_index: int


@dataclass(frozen=True)
class AuditRow:
    protein_id: str
    sequence_md5: str
    sequence_length: int
    analysis: str
    signature_accession: str
    signature_description: str
    start: int
    end: int
    score_or_evalue: str
    status: str
    run_date: str
    interpro_accession: str
    interpro_description: str
    go_terms: str
    pathway_terms: str
    architecture_index: int
    known_protein_id: bool
    included_in_domain_table: bool
    interproscan_version: str
    source_file: str


@dataclass(frozen=True)
class ConversionResult:
    domains: tuple[DomainRow, ...]
    audit_rows: tuple[AuditRow, ...]
    coverage: tuple[tuple[str, str, str], ...]
    exact_duplicates_excluded: int
    unknown_rows_excluded: int


def _optional(value: str) -> str:
    stripped = value.strip()
    return "" if stripped == "-" else stripped


def _require(value: str, field: str, line_number: int) -> str:
    normalized = _optional(value)
    if not normalized:
        raise InputValidationError(f"InterProScan TSV has empty {field} on line {line_number}")
    return normalized


def parse_interproscan_tsv(path: Path) -> tuple[list[InterProHit], int]:
    """Parse InterProScan 5 TSV, removing only exact duplicate rows."""
    source_path = path.expanduser().resolve()
    if not source_path.is_file():
        raise InputValidationError(f"InterProScan TSV not found: {source_path}")
    hits: list[InterProHit] = []
    seen_rows: set[tuple[str, ...]] = set()
    duplicate_count = 0
    with source_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].lstrip().startswith("#"):
                continue
            if len(row) not in ALLOWED_COLUMN_COUNTS:
                raise InputValidationError(
                    "InterProScan TSV must contain 11, 13, 14, or 15 columns; "
                    f"found {len(row)} on line {line_number}"
                )
            raw_identity = tuple(row)
            if raw_identity in seen_rows:
                duplicate_count += 1
                continue
            seen_rows.add(raw_identity)
            try:
                length = int(_require(row[2], "sequence length", line_number))
                start = int(_require(row[6], "start", line_number))
                end = int(_require(row[7], "end", line_number))
            except ValueError as exc:
                raise InputValidationError(
                    f"Invalid numeric InterProScan value on line {line_number}: {exc}"
                ) from exc
            if length < 1 or start < 1 or end < start or end > length:
                raise InputValidationError(
                    f"Invalid InterProScan coordinates on line {line_number}: "
                    f"length={length}, start={start}, end={end}"
                )
            padded = row + [""] * (15 - len(row))
            hits.append(
                InterProHit(
                    protein_id=_require(padded[0], "protein accession", line_number),
                    sequence_md5=_require(padded[1], "sequence MD5", line_number),
                    sequence_length=length,
                    analysis=_require(padded[3], "analysis", line_number),
                    signature_accession=_require(padded[4], "signature accession", line_number),
                    signature_description=_optional(padded[5]),
                    start=start,
                    end=end,
                    score_or_evalue=_optional(padded[8]),
                    status=_require(padded[9], "status", line_number),
                    run_date=_optional(padded[10]),
                    interpro_accession=_optional(padded[11]),
                    interpro_description=_optional(padded[12]),
                    go_terms=_optional(padded[13]),
                    pathway_terms=_optional(padded[14]),
                )
            )
    return hits, duplicate_count


def _hit_sort_key(hit: InterProHit) -> tuple[object, ...]:
    return (
        hit.protein_id,
        hit.start,
        hit.end,
        hit.analysis.casefold(),
        hit.signature_accession,
        hit.interpro_accession,
        hit.signature_description,
        hit.score_or_evalue,
    )


def _percent(count: int, total: int) -> str:
    return f"{(100.0 * count / total) if total else 0.0:.6f}"


def _coverage_rows(
    hits: Sequence[InterProHit],
    descriptions: dict[str, str],
    query_id: str,
    exact_duplicates: int,
    unknown_rows: int,
) -> tuple[tuple[str, str, str], ...]:
    protein_ids = set(descriptions) or {hit.protein_id for hit in hits}
    hit_ids = {hit.protein_id for hit in hits} & protein_ids
    hypothetical_ids = {
        protein_id
        for protein_id, description in descriptions.items()
        if any(term in description.casefold() for term in _HYPOTHETICAL_TERMS)
    }
    non_hypothetical_ids = protein_ids - hypothetical_ids
    by_protein: dict[str, list[InterProHit]] = defaultdict(list)
    for hit in hits:
        by_protein[hit.protein_id].append(hit)
    repeated = sum(
        any(
            count > 1
            for count in Counter(
                (hit.analysis, hit.signature_accession) for hit in protein_hits
            ).values()
        )
        for protein_id, protein_hits in by_protein.items()
        if protein_id in protein_ids
    )
    metrics: list[tuple[str, int]] = [
        ("total_proteins", len(protein_ids)),
        ("proteins_with_domain", len(hit_ids)),
        ("proteins_without_domain", len(protein_ids - hit_ids)),
        ("total_domain_hits", len(hits)),
        ("unique_domain_accessions", len({hit.signature_accession for hit in hits})),
        ("pfam_hit_count", sum(hit.analysis.casefold() == "pfam" for hit in hits)),
        (
            "pfam_protein_coverage",
            len(
                {hit.protein_id for hit in hits if hit.analysis.casefold() == "pfam"} & protein_ids
            ),
        ),
        ("interpro_accession_hit_count", sum(bool(hit.interpro_accession) for hit in hits)),
        (
            "interpro_protein_coverage",
            len({hit.protein_id for hit in hits if hit.interpro_accession} & protein_ids),
        ),
        ("go_hit_count", sum(bool(hit.go_terms) for hit in hits)),
        (
            "go_protein_coverage",
            len({hit.protein_id for hit in hits if hit.go_terms} & protein_ids),
        ),
        ("pathway_hit_count", sum(bool(hit.pathway_terms) for hit in hits)),
        (
            "pathway_protein_coverage",
            len({hit.protein_id for hit in hits if hit.pathway_terms} & protein_ids),
        ),
        ("hypothetical_proteins", len(hypothetical_ids)),
        ("hypothetical_proteins_with_domain", len(hypothetical_ids & hit_ids)),
        ("non_hypothetical_proteins", len(non_hypothetical_ids)),
        ("non_hypothetical_proteins_with_domain", len(non_hypothetical_ids & hit_ids)),
        ("query_domain_count", len(by_protein.get(query_id, []))),
        (
            "multi_domain_protein_count",
            sum(
                len(protein_hits) > 1
                for protein_id, protein_hits in by_protein.items()
                if protein_id in protein_ids
            ),
        ),
        ("repeated_domain_protein_count", repeated),
        ("exact_duplicate_rows_excluded", exact_duplicates),
        ("unknown_protein_rows_excluded", unknown_rows),
        ("malformed_rows", 0),
    ]
    for analysis, count in sorted(Counter(hit.analysis for hit in hits).items()):
        metrics.append((f"analysis:{analysis}:hit_count", count))
    return tuple((name, str(count), _percent(count, len(protein_ids))) for name, count in metrics)


def convert_interproscan(
    *,
    input_path: Path,
    proteome_fasta: Path | None = None,
    query_id: str = "",
    interproscan_version: str = "",
) -> ConversionResult:
    """Convert hits deterministically and audit IDs against a proteome."""
    parsed_hits, exact_duplicates = parse_interproscan_tsv(input_path)
    proteins = LocalFastaLoader().load(proteome_fasta) if proteome_fasta else []
    descriptions = {record.protein_id: record.description for record in proteins}
    known_ids = set(descriptions)
    if not proteome_fasta:
        known_ids = {hit.protein_id for hit in parsed_hits}
        descriptions = dict.fromkeys(known_ids, "")
    hits = sorted(
        (hit for hit in parsed_hits if hit.protein_id in known_ids),
        key=_hit_sort_key,
    )
    unknown_rows = len(parsed_hits) - len(hits)
    by_protein: dict[str, list[InterProHit]] = defaultdict(list)
    for hit in hits:
        by_protein[hit.protein_id].append(hit)
    indexes = {
        hit: index for protein_hits in by_protein.values() for index, hit in enumerate(protein_hits)
    }
    formal_identity: set[tuple[str, str, str, int, int]] = set()
    domains: list[DomainRow] = []
    for hit in hits:
        identity = (
            hit.protein_id,
            hit.analysis,
            hit.signature_accession,
            hit.start,
            hit.end,
        )
        if identity in formal_identity:
            raise InputValidationError(
                f"InterProScan rows collapse to a conflicting duplicate domain identity: {identity}"
            )
        formal_identity.add(identity)
        domains.append(
            DomainRow(
                protein_id=hit.protein_id,
                source=hit.analysis,
                accession=hit.signature_accession,
                name=hit.signature_description,
                start=hit.start,
                end=hit.end,
                architecture_index=indexes[hit],
            )
        )
    audit_rows = tuple(
        AuditRow(
            **asdict(hit),
            architecture_index=indexes.get(hit, -1),
            known_protein_id=hit.protein_id in known_ids,
            included_in_domain_table=hit.protein_id in known_ids,
            interproscan_version=interproscan_version,
            source_file=str(input_path.expanduser().resolve()),
        )
        for hit in sorted(parsed_hits, key=_hit_sort_key)
    )
    return ConversionResult(
        domains=tuple(domains),
        audit_rows=audit_rows,
        coverage=_coverage_rows(hits, descriptions, query_id, exact_duplicates, unknown_rows),
        exact_duplicates_excluded=exact_duplicates,
        unknown_rows_excluded=unknown_rows,
    )


def _write_dicts(path: Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_conversion_outputs(
    result: ConversionResult,
    *,
    domain_output: Path,
    audit_output: Path,
    coverage_output: Path,
) -> None:
    """Write formal, audit, and coverage TSVs as UTF-8 with LF endings."""
    _write_dicts(domain_output, DOMAIN_COLUMNS, (asdict(row) for row in result.domains))
    _write_dicts(audit_output, AUDIT_COLUMNS, (asdict(row) for row in result.audit_rows))
    _write_dicts(
        coverage_output,
        ("metric", "count", "percent_of_scope"),
        (
            {"metric": name, "count": count, "percent_of_scope": percent}
            for name, count, percent in result.coverage
        ),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--proteome-fasta", type=Path)
    parser.add_argument("--query-id", default="")
    parser.add_argument("--interproscan-version", default="")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = convert_interproscan(
        input_path=args.input,
        proteome_fasta=args.proteome_fasta,
        query_id=args.query_id,
        interproscan_version=args.interproscan_version,
    )
    write_conversion_outputs(
        result,
        domain_output=args.output,
        audit_output=args.audit_output,
        coverage_output=args.coverage_output,
    )
    print(f"Formal domain rows: {len(result.domains)}")
    print(f"Audited InterProScan rows: {len(result.audit_rows)}")
    print(f"Exact duplicates excluded: {result.exact_duplicates_excluded}")
    print(f"Unknown protein rows excluded: {result.unknown_rows_excluded}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
