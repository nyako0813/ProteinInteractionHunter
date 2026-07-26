#!/usr/bin/env python3
"""Build a deterministic ProteinInteractionHunter annotation TSV from NCBI FASTA/GFF.

The workflow uses only explicit source annotations. It never infers functional category,
localization, transmembrane topology, or annotation confidence from product text.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from protein_interaction_hunter.adapters.local.annotation import ANNOTATION_COLUMNS
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.adapters.local.gff import LocalGff3Loader
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.genome import GeneCoordinate

AUDIT_COLUMNS = (
    "protein_id",
    "fasta_description",
    "gff_product",
    "gene_feature_id",
    "gene_name",
    "locus_tag",
    "old_locus_tags",
    "ncbi_gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "cds_record_count",
    "duplicate_cds_count",
    "mapping_status",
    "mapping_warning",
    "product_match",
)

_TAXON_SUFFIX = re.compile(r"\s+\[([^=\[\]]+)\]\s*$")
_HYPOTHETICAL_TERMS = ("hypothetical", "uncharacterized")


@dataclass(frozen=True)
class AnnotationRow:
    protein_id: str
    gene_name: str = ""
    locus_tag: str = ""
    product: str = ""
    functional_category: str = ""
    localization_annotation: str = ""
    transmembrane_annotation: str = ""
    annotation_source: str = ""
    annotation_confidence: str = ""


@dataclass(frozen=True)
class AuditRow:
    protein_id: str
    fasta_description: str
    gff_product: str
    gene_feature_id: str
    gene_name: str
    locus_tag: str
    old_locus_tags: str
    ncbi_gene_id: str
    contig: str
    start: int | str
    end: int | str
    strand: str
    cds_record_count: int
    duplicate_cds_count: int
    mapping_status: str
    mapping_warning: str
    product_match: str


@dataclass(frozen=True)
class BuildSummary:
    total_proteins: int
    gff_cds_records: int
    gff_distinct_protein_ids: int
    fasta_gff_intersection: int
    fasta_only: int
    gff_only: int
    unique_one_to_one: int
    multi_locus_ambiguous: int
    formal_annotation_rows: int
    gene_name_coverage: int
    locus_tag_coverage: int
    product_coverage: int
    functional_category_coverage: int
    localization_coverage: int
    transmembrane_coverage: int
    annotation_confidence_coverage: int
    missing_product: int
    missing_locus_tag: int
    missing_gene_name: int
    hypothetical_or_uncharacterized_product: int
    duplicate_cds_records: int
    duplicate_annotation_rows: int
    malformed_rows: int
    product_mismatch_accessions: int
    query_id: str
    query_covered: bool


@dataclass(frozen=True)
class BuildResult:
    annotations: tuple[AnnotationRow, ...]
    audit_rows: tuple[AuditRow, ...]
    summary: BuildSummary


def clean_fasta_description(value: str) -> str:
    """Remove only a final non-key/value taxon suffix such as ``[Organism name]``."""
    return _TAXON_SUFFIX.sub("", value).strip()


def _unique(values: Iterable[str | None]) -> list[str]:
    return sorted({value.strip() for value in values if value and value.strip()})


def _attribute_values(
    records: Sequence[GeneCoordinate],
    parent_genes: dict[str, GeneCoordinate],
    key: str,
) -> list[str]:
    values: list[str | None] = []
    for record in records:
        values.extend(record.attributes.get(key, []))
        for parent_id in record.parent_ids:
            parent = parent_genes.get(parent_id)
            if parent is not None:
                values.extend(parent.attributes.get(key, []))
    return _unique(values)


def _locus_identity(record: GeneCoordinate) -> tuple[str, ...]:
    if record.parent_ids:
        return ("parent", *sorted(record.parent_ids))
    if record.locus_tag:
        return ("locus_tag", record.locus_tag)
    return (
        "coordinate",
        record.seqid,
        str(record.start),
        str(record.end),
        record.strand or "?",
    )


def _cds_signature(record: GeneCoordinate) -> tuple[object, ...]:
    return (
        record.seqid,
        record.start,
        record.end,
        record.strand,
        tuple(sorted(record.parent_ids)),
        record.protein_id,
        record.locus_tag,
        tuple(record.attributes.get("product", [])),
    )


def _ncbi_gene_ids(
    records: Sequence[GeneCoordinate], parent_genes: dict[str, GeneCoordinate]
) -> list[str]:
    return _unique(
        value.removeprefix("GeneID:")
        for value in _attribute_values(records, parent_genes, "Dbxref")
        if value.startswith("GeneID:")
    )


def _single_or_blank(values: Sequence[str], warning: str, warnings: list[str]) -> str:
    if len(values) == 1:
        return values[0]
    if len(values) > 1:
        warnings.append(warning)
    return ""


def _audit_locus(
    *,
    protein_id: str,
    fasta_description: str,
    records: Sequence[GeneCoordinate],
    parent_genes: dict[str, GeneCoordinate],
    mapping_status: str,
    extra_warnings: Sequence[str] = (),
) -> AuditRow:
    warnings = list(extra_warnings)
    products = _attribute_values(records, parent_genes, "product")
    gene_names = _attribute_values(records, parent_genes, "gene")
    locus_tags = _attribute_values(records, parent_genes, "locus_tag")
    gene_ids = _ncbi_gene_ids(records, parent_genes)
    contigs = _unique(record.seqid for record in records)
    strands = _unique(record.strand for record in records)
    parent_ids = _unique(parent for record in records for parent in record.parent_ids)
    old_locus_tags = _attribute_values(records, parent_genes, "old_locus_tag")
    signatures = Counter(_cds_signature(record) for record in records)
    duplicate_count = sum(count - 1 for count in signatures.values())
    if duplicate_count:
        warnings.append("duplicate_cds_record")
    gff_product = _single_or_blank(products, "conflicting_gff_product", warnings)
    gene_name = _single_or_blank(gene_names, "conflicting_gene_name", warnings)
    locus_tag = _single_or_blank(locus_tags, "conflicting_locus_tag", warnings)
    ncbi_gene_id = _single_or_blank(gene_ids, "conflicting_ncbi_gene_id", warnings)
    contig = _single_or_blank(contigs, "conflicting_contig", warnings)
    strand = _single_or_blank(strands, "conflicting_strand", warnings) or "?"
    cleaned_fasta = clean_fasta_description(fasta_description)
    product_match = ""
    if cleaned_fasta and gff_product:
        product_match = str(cleaned_fasta.casefold() == gff_product.casefold()).lower()
        if product_match == "false":
            warnings.append("fasta_gff_product_mismatch")
    return AuditRow(
        protein_id=protein_id,
        fasta_description=cleaned_fasta,
        gff_product=gff_product,
        gene_feature_id="|".join(parent_ids),
        gene_name=gene_name,
        locus_tag=locus_tag,
        old_locus_tags="|".join(old_locus_tags),
        ncbi_gene_id=ncbi_gene_id,
        contig=contig,
        start=min(record.start for record in records),
        end=max(record.end for record in records),
        strand=strand,
        cds_record_count=len(records),
        duplicate_cds_count=duplicate_count,
        mapping_status=mapping_status,
        mapping_warning="|".join(sorted(set(warnings))),
        product_match=product_match,
    )


def build_annotation_tables(
    *,
    fasta_path: Path,
    gff_path: Path,
    annotation_source: str,
    query_id: str = "",
) -> BuildResult:
    """Build formal and audit records without any biological inference."""
    proteins = LocalFastaLoader().load(fasta_path)
    document = LocalGff3Loader().load_document(gff_path)
    fasta_by_id = {record.protein_id: record for record in proteins}
    parent_genes = {
        record.feature_id: record
        for record in document.features
        if record.feature_type.casefold() == "gene" and record.feature_id
    }
    cds_records = [
        record
        for record in document.features
        if record.feature_type.casefold() == "cds" and record.protein_id
    ]
    cds_by_protein: dict[str, list[GeneCoordinate]] = defaultdict(list)
    for record in cds_records:
        assert record.protein_id is not None
        cds_by_protein[record.protein_id].append(record)

    fasta_ids = set(fasta_by_id)
    gff_ids = set(cds_by_protein)
    annotations: list[AnnotationRow] = []
    audit_rows: list[AuditRow] = []
    ambiguous_ids: set[str] = set()

    for protein_id in sorted(fasta_ids | gff_ids):
        fasta_record = fasta_by_id.get(protein_id)
        fasta_description = fasta_record.description if fasta_record else ""
        records = cds_by_protein.get(protein_id, [])
        if not records:
            cleaned = clean_fasta_description(fasta_description)
            annotations.append(
                AnnotationRow(
                    protein_id=protein_id,
                    product=cleaned,
                    annotation_source=annotation_source,
                )
            )
            audit_rows.append(
                AuditRow(
                    protein_id=protein_id,
                    fasta_description=cleaned,
                    gff_product="",
                    gene_feature_id="",
                    gene_name="",
                    locus_tag="",
                    old_locus_tags="",
                    ncbi_gene_id="",
                    contig="",
                    start="",
                    end="",
                    strand="",
                    cds_record_count=0,
                    duplicate_cds_count=0,
                    mapping_status="fasta_only",
                    mapping_warning="missing_gff_cds",
                    product_match="",
                )
            )
            continue

        by_locus: dict[tuple[str, ...], list[GeneCoordinate]] = defaultdict(list)
        for record in records:
            by_locus[_locus_identity(record)].append(record)
        is_ambiguous = len(by_locus) > 1
        if is_ambiguous:
            ambiguous_ids.add(protein_id)
        for locus_key in sorted(by_locus):
            audit_rows.append(
                _audit_locus(
                    protein_id=protein_id,
                    fasta_description=fasta_description,
                    records=by_locus[locus_key],
                    parent_genes=parent_genes,
                    mapping_status=(
                        "ambiguous_multi_locus"
                        if is_ambiguous
                        else ("gff_only" if fasta_record is None else "one_to_one")
                    ),
                    extra_warnings=(
                        ("excluded_from_formal_annotation:multi_locus",)
                        if is_ambiguous
                        else (("absent_from_proteome",) if fasta_record is None else ())
                    ),
                )
            )
        if is_ambiguous or fasta_record is None:
            continue

        locus_audit = audit_rows[-1]
        gff_product = locus_audit.gff_product
        product = gff_product or clean_fasta_description(fasta_description)
        annotations.append(
            AnnotationRow(
                protein_id=protein_id,
                gene_name=locus_audit.gene_name,
                locus_tag=locus_audit.locus_tag,
                product=product,
                annotation_source=annotation_source,
            )
        )

    annotations.sort(key=lambda row: row.protein_id)
    audit_rows.sort(
        key=lambda row: (
            row.protein_id,
            str(row.contig),
            int(row.start) if isinstance(row.start, int) else -1,
            int(row.end) if isinstance(row.end, int) else -1,
        )
    )
    duplicate_cds_records = sum(row.duplicate_cds_count for row in audit_rows)
    mismatches = {row.protein_id for row in audit_rows if row.product_match == "false"}
    summary = BuildSummary(
        total_proteins=len(proteins),
        gff_cds_records=len(cds_records),
        gff_distinct_protein_ids=len(gff_ids),
        fasta_gff_intersection=len(fasta_ids & gff_ids),
        fasta_only=len(fasta_ids - gff_ids),
        gff_only=len(gff_ids - fasta_ids),
        unique_one_to_one=len((fasta_ids & gff_ids) - ambiguous_ids),
        multi_locus_ambiguous=len(ambiguous_ids & fasta_ids),
        formal_annotation_rows=len(annotations),
        gene_name_coverage=sum(bool(row.gene_name) for row in annotations),
        locus_tag_coverage=sum(bool(row.locus_tag) for row in annotations),
        product_coverage=sum(bool(row.product) for row in annotations),
        functional_category_coverage=sum(bool(row.functional_category) for row in annotations),
        localization_coverage=sum(bool(row.localization_annotation) for row in annotations),
        transmembrane_coverage=sum(bool(row.transmembrane_annotation) for row in annotations),
        annotation_confidence_coverage=sum(bool(row.annotation_confidence) for row in annotations),
        missing_product=sum(not row.product for row in annotations),
        missing_locus_tag=sum(not row.locus_tag for row in annotations),
        missing_gene_name=sum(not row.gene_name for row in annotations),
        hypothetical_or_uncharacterized_product=sum(
            any(term in row.product.casefold() for term in _HYPOTHETICAL_TERMS)
            for row in annotations
        ),
        duplicate_cds_records=duplicate_cds_records,
        duplicate_annotation_rows=(len(annotations) - len({row.protein_id for row in annotations})),
        malformed_rows=0,
        product_mismatch_accessions=len(mismatches),
        query_id=query_id,
        query_covered=bool(query_id and any(row.protein_id == query_id for row in annotations)),
    )
    return BuildResult(tuple(annotations), tuple(audit_rows), summary)


def _write_dataclass_tsv(
    path: Path,
    rows: Sequence[AnnotationRow | AuditRow],
    column_names: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(column_names), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_build_outputs(
    result: BuildResult,
    *,
    annotation_output: Path,
    audit_output: Path,
    coverage_output: Path | None = None,
) -> None:
    """Write UTF-8/LF formal annotation, mapping audit, and optional coverage TSV."""
    _write_dataclass_tsv(annotation_output, result.annotations, ANNOTATION_COLUMNS)
    _write_dataclass_tsv(audit_output, result.audit_rows, AUDIT_COLUMNS)
    if coverage_output is None:
        return
    coverage_output.parent.mkdir(parents=True, exist_ok=True)
    total = result.summary.total_proteins
    with coverage_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "count", "percent_of_proteome"])
        for field in fields(result.summary):
            value = getattr(result.summary, field.name)
            if isinstance(value, (bool, str)):
                writer.writerow([field.name, str(value).lower(), ""])
            else:
                percent = (100.0 * value / total) if total else 0.0
                writer.writerow([field.name, value, f"{percent:.6f}"])


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a deterministic NCBI-derived annotation table and mapping audit."
    )
    parser.add_argument("--fasta", type=Path, required=True, help="NCBI protein FASTA")
    parser.add_argument("--gff", type=Path, required=True, help="Matching NCBI GFF3")
    parser.add_argument(
        "--annotation-output", type=Path, required=True, help="Formal annotation TSV"
    )
    parser.add_argument("--audit-output", type=Path, required=True, help="Mapping audit TSV")
    parser.add_argument("--coverage-output", type=Path, help="Optional coverage TSV")
    parser.add_argument(
        "--annotation-source",
        required=True,
        help="Literal source label written to each formal row",
    )
    parser.add_argument("--query-id", default="", help="Optional query coverage check")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        result = build_annotation_tables(
            fasta_path=args.fasta,
            gff_path=args.gff,
            annotation_source=args.annotation_source,
            query_id=args.query_id,
        )
        write_build_outputs(
            result,
            annotation_output=args.annotation_output,
            audit_output=args.audit_output,
            coverage_output=args.coverage_output,
        )
    except (InputValidationError, OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    for field in fields(result.summary):
        print(f"{field.name}: {getattr(result.summary, field.name)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
