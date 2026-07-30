#!/usr/bin/env python3
"""Convert PSORTb 3 terse output into the formal annotation-table contract.

PSORTb-specific score, provenance, mapping, and row-level decisions are kept in
an audit TSV. The formal output retains the existing nine-column annotation
loader contract. No product/domain text, signal peptide, or transmembrane
topology is inferred.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from protein_interaction_hunter.adapters.local.annotation import (
    ANNOTATION_COLUMNS,
    LocalAnnotationTsvLoader,
)
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.exceptions import InputValidationError

PSORTB_MAPPING_RULE_VERSION = "psortb-3.0-archaea-terse-mapping-v1"
PSORTB_SOURCE_ROLE = "predicted_subcellular_localization"

AUDIT_COLUMNS = (
    "protein_id",
    "raw_prediction",
    "raw_score",
    "formal_localization",
    "mapping_status",
    "mapping_rule_version",
    "known_protein_id",
    "included_in_formal_table",
    "duplicate_kind",
    "raw_row_number",
    "normalization_decision",
    "exclusion_reason",
    "source_name",
    "source_version",
    "source_role",
    "organism_mode",
    "source_command",
    "source_file",
    "raw_input_sha256",
    "formal_output_sha256",
)


@dataclass(frozen=True)
class MappingRule:
    raw_label: str
    formal_label: str
    status: str


PSORTB_ARCHAEA_MAPPING = (
    MappingRule("Cytoplasmic", "cytosolic", "accepted"),
    MappingRule("CytoplasmicMembrane", "membrane", "accepted"),
    MappingRule("Cellwall", "", "unsupported"),
    MappingRule("Extracellular", "secreted", "accepted"),
    MappingRule("Unknown", "", "unknown"),
)


def _mapping_key(value: str) -> str:
    return "".join(character for character in value.casefold() if character.isalnum())


_MAPPING_BY_KEY = {_mapping_key(rule.raw_label): rule for rule in PSORTB_ARCHAEA_MAPPING}


@dataclass(frozen=True)
class PsortbRow:
    protein_id: str
    raw_prediction: str
    score: float
    score_text: str
    line_number: int
    duplicate_kind: str = ""


@dataclass(frozen=True)
class AuditRow:
    protein_id: str
    raw_prediction: str
    raw_score: str
    formal_localization: str
    mapping_status: str
    mapping_rule_version: str
    known_protein_id: bool
    included_in_formal_table: bool
    duplicate_kind: str
    raw_row_number: int
    normalization_decision: str
    exclusion_reason: str
    source_name: str
    source_version: str
    source_role: str
    organism_mode: str
    source_command: str
    source_file: str
    raw_input_sha256: str


@dataclass(frozen=True)
class ConversionResult:
    formal_rows: tuple[dict[str, object], ...]
    audit_rows: tuple[AuditRow, ...]
    coverage: tuple[tuple[str, str, str], ...]
    raw_input_sha256: str
    exact_duplicate_count: int
    unknown_protein_id_count: int
    missing_protein_id_count: int
    query_prediction: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _score_text(value: float) -> str:
    return format(value, ".15g")


def parse_psortb_terse(path: Path) -> tuple[list[PsortbRow], int]:
    """Parse headerless or headered official three-column PSORTb terse output."""
    source_path = path.expanduser().resolve()
    if not source_path.is_file():
        raise InputValidationError(f"PSORTb terse output not found: {source_path}")

    rows: list[PsortbRow] = []
    exact_duplicates = 0
    seen_by_id: dict[str, tuple[str, float]] = {}
    with source_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, values in enumerate(reader, start=1):
            if not values or (len(values) == 1 and not values[0].strip()):
                continue
            if values[0].lstrip().startswith("#"):
                continue
            if [value.strip().casefold() for value in values] == [
                "seqid",
                "localization",
                "score",
            ]:
                continue
            if len(values) != 3:
                raise InputValidationError(
                    "PSORTb terse output must contain exactly 3 columns "
                    f"(SeqId, Localization, Score); found {len(values)} on line {line_number}"
                )
            protein_id, prediction, score_raw = (value.strip() for value in values)
            if not protein_id:
                raise InputValidationError(
                    f"PSORTb terse output has empty SeqId on line {line_number}"
                )
            if not prediction:
                raise InputValidationError(
                    f"PSORTb terse output has empty Localization on line {line_number}"
                )
            try:
                score = float(score_raw)
            except ValueError as exc:
                raise InputValidationError(
                    f"Invalid PSORTb score on line {line_number}: {score_raw!r}"
                ) from exc
            if not math.isfinite(score) or not 0.0 <= score <= 10.0:
                raise InputValidationError(
                    f"PSORTb score must be finite and between 0 and 10 on line {line_number}"
                )

            identity = (prediction, score)
            prior = seen_by_id.get(protein_id)
            if prior is not None:
                if prior != identity:
                    raise InputValidationError(
                        "Conflicting duplicate PSORTb SeqId "
                        f"{protein_id!r} on line {line_number}"
                    )
                exact_duplicates += 1
                rows.append(
                    PsortbRow(
                        protein_id=protein_id,
                        raw_prediction=prediction,
                        score=score,
                        score_text=_score_text(score),
                        line_number=line_number,
                        duplicate_kind="exact_duplicate",
                    )
                )
                continue
            seen_by_id[protein_id] = identity
            rows.append(
                PsortbRow(
                    protein_id=protein_id,
                    raw_prediction=prediction,
                    score=score,
                    score_text=_score_text(score),
                    line_number=line_number,
                )
            )
    return rows, exact_duplicates


def _percent(count: int, total: int) -> str:
    return f"{(100.0 * count / total) if total else 0.0:.6f}"


def _source_value(existing: str | None, source_version: str) -> str:
    psortb = f"PSORTb {source_version}" if source_version else "PSORTb"
    values = [value for value in (existing, psortb) if value]
    return "; ".join(dict.fromkeys(values))


def _formal_annotation_row(record: object | None, protein_id: str) -> dict[str, object]:
    confidence = getattr(record, "annotation_confidence", None)
    return {
        "protein_id": protein_id,
        "gene_name": getattr(record, "gene_name", None) or "",
        "locus_tag": getattr(record, "locus_tag", None) or "",
        "product": getattr(record, "product", None) or "",
        "functional_category": getattr(record, "functional_category", None) or "",
        "localization_annotation": getattr(record, "localization_annotation", None) or "",
        "transmembrane_annotation": getattr(record, "transmembrane_annotation", None) or "",
        "annotation_source": getattr(record, "annotation_source", None) or "",
        "annotation_confidence": confidence if confidence is not None else "",
    }


def _rule_for(raw_prediction: str) -> MappingRule:
    return _MAPPING_BY_KEY.get(
        _mapping_key(raw_prediction),
        MappingRule(raw_prediction, "", "unsupported"),
    )


def convert_psortb(
    *,
    input_path: Path,
    proteome_fasta: Path | None = None,
    base_annotation_table: Path | None = None,
    query_id: str = "",
    source_version: str = "",
    source_command: str = "",
    organism_mode: str = "archaea",
) -> ConversionResult:
    """Convert PSORTb rows, preserving unknown/unsupported results only in audit."""
    if organism_mode != "archaea":
        raise InputValidationError("This workflow requires organism_mode='archaea'")

    parsed_rows, exact_duplicates = parse_psortb_terse(input_path)
    unique_rows = [row for row in parsed_rows if not row.duplicate_kind]
    prediction_by_id = {row.protein_id: row for row in unique_rows}
    if query_id and query_id not in prediction_by_id:
        raise InputValidationError(f"Query ID missing from PSORTb output: {query_id}")

    proteins = LocalFastaLoader().load(proteome_fasta) if proteome_fasta else []
    proteome_ids = {protein.protein_id for protein in proteins}
    base_records = (
        LocalAnnotationTsvLoader().load(base_annotation_table)
        if base_annotation_table
        else []
    )
    base_by_id = {record.protein_id: record for record in base_records}

    if base_records:
        # Preserve the base table's one-row-per-ID scope. Adding otherwise
        # absent IDs would change candidate annotation presence even while
        # localization evidence is disabled in a coverage-only pilot.
        formal_scope = sorted(base_by_id)
    elif proteins:
        formal_scope = sorted(proteome_ids)
    else:
        formal_scope = sorted(prediction_by_id)
        proteome_ids = set(formal_scope)
    formal_ids = set(formal_scope)

    unknown_ids = set(prediction_by_id) - proteome_ids
    missing_ids = proteome_ids - set(prediction_by_id)
    raw_hash = _sha256(input_path.expanduser().resolve())

    formal_rows: list[dict[str, object]] = []
    for protein_id in formal_scope:
        base = base_by_id.get(protein_id)
        formal = _formal_annotation_row(base, protein_id)
        prediction = prediction_by_id.get(protein_id)
        if prediction is not None:
            rule = _rule_for(prediction.raw_prediction)
            formal_label = rule.formal_label if rule.status == "accepted" else ""
            existing = str(formal["localization_annotation"])
            if formal_label and existing and existing != formal_label:
                raise InputValidationError(
                    "Base annotation conflicts with PSORTb formal localization for "
                    f"{protein_id}: {existing!r} versus {formal_label!r}"
                )
            if formal_label:
                formal["localization_annotation"] = formal_label
            formal["annotation_source"] = _source_value(
                str(formal["annotation_source"]) or None,
                source_version,
            )
        formal_rows.append(formal)

    audit_rows: list[AuditRow] = []
    for row in sorted(parsed_rows, key=lambda item: (item.protein_id, item.line_number)):
        rule = _rule_for(row.raw_prediction)
        known = row.protein_id in proteome_ids
        included = (
            known
            and row.protein_id in formal_ids
            and not row.duplicate_kind
            and rule.status == "accepted"
            and bool(rule.formal_label)
        )
        if row.duplicate_kind:
            exclusion = "exact_duplicate"
        elif not known:
            exclusion = "unknown_protein_id"
        elif row.protein_id not in formal_ids:
            exclusion = "missing_base_annotation"
        elif rule.status == "unknown":
            exclusion = "prediction_unknown"
        elif rule.status == "unsupported":
            exclusion = "unsupported_prediction"
        else:
            exclusion = ""
        if rule.status == "accepted":
            decision = f"explicit_mapping:{rule.raw_label}->{rule.formal_label}"
        elif rule.status == "unknown":
            decision = "preserved_as_missing"
        else:
            decision = "excluded_from_formal_localization"
        audit_rows.append(
            AuditRow(
                protein_id=row.protein_id,
                raw_prediction=row.raw_prediction,
                raw_score=row.score_text,
                formal_localization=rule.formal_label,
                mapping_status=rule.status,
                mapping_rule_version=PSORTB_MAPPING_RULE_VERSION,
                known_protein_id=known,
                included_in_formal_table=included,
                duplicate_kind=row.duplicate_kind,
                raw_row_number=row.line_number,
                normalization_decision=decision,
                exclusion_reason=exclusion,
                source_name="PSORTb",
                source_version=source_version,
                source_role=PSORTB_SOURCE_ROLE,
                organism_mode=organism_mode,
                source_command=source_command,
                source_file=str(input_path.expanduser().resolve()),
                raw_input_sha256=raw_hash,
            )
        )

    known_unique = [row for row in unique_rows if row.protein_id in proteome_ids]
    mapped = [row for row in known_unique if _rule_for(row.raw_prediction).status == "accepted"]
    unknown = [row for row in known_unique if _rule_for(row.raw_prediction).status == "unknown"]
    unsupported = len(known_unique) - len(mapped) - len(unknown)
    distribution = Counter(_rule_for(row.raw_prediction).formal_label for row in mapped)
    metrics: list[tuple[str, int]] = [
        ("proteome_total", len(proteome_ids)),
        ("raw_data_rows", len(parsed_rows)),
        ("unique_prediction_rows", len(unique_rows)),
        ("known_proteins_represented", len({row.protein_id for row in known_unique})),
        ("formal_table_rows", len(formal_scope)),
        (
            "mapped_predictions_in_formal_table",
            sum(row.protein_id in formal_ids for row in mapped),
        ),
        ("non_unknown_predictions", len(mapped)),
        ("unknown_predictions", len(unknown)),
        ("proteins_absent_from_raw_output", len(missing_ids)),
        ("malformed_rows", 0),
        ("exact_duplicate_rows", exact_duplicates),
        ("conflicting_duplicate_rows", 0),
        ("unknown_protein_ids", len(unknown_ids)),
        ("unsupported_predictions", unsupported),
        (
            "missing_base_annotation_predictions",
            sum(row.protein_id not in formal_ids for row in known_unique),
        ),
        ("excluded_rows", sum(bool(row.exclusion_reason) for row in audit_rows)),
    ]
    for label, count in sorted(distribution.items()):
        metrics.append((f"localization:{label}", count))
    query_prediction = (
        prediction_by_id[query_id].raw_prediction
        if query_id and query_id in prediction_by_id
        else ""
    )
    return ConversionResult(
        formal_rows=tuple(formal_rows),
        audit_rows=tuple(audit_rows),
        coverage=tuple(
            (name, str(count), _percent(count, len(proteome_ids)))
            for name, count in metrics
        ),
        raw_input_sha256=raw_hash,
        exact_duplicate_count=exact_duplicates,
        unknown_protein_id_count=len(unknown_ids),
        missing_protein_id_count=len(missing_ids),
        query_prediction=query_prediction,
    )


def _write_dicts(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[dict[str, object]],
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


def write_conversion_outputs(
    result: ConversionResult,
    *,
    annotation_output: Path,
    audit_output: Path,
    coverage_output: Path,
    metadata_output: Path,
) -> None:
    """Write formal annotation, audit, coverage, and metadata as UTF-8/LF."""
    _write_dicts(annotation_output, ANNOTATION_COLUMNS, result.formal_rows)
    formal_hash = _sha256(annotation_output)
    _write_dicts(
        audit_output,
        AUDIT_COLUMNS,
        (
            asdict(row) | {"formal_output_sha256": formal_hash}
            for row in result.audit_rows
        ),
    )
    _write_dicts(
        coverage_output,
        ("metric", "count", "percent_of_proteome"),
        (
            {
                "metric": name,
                "count": count,
                "percent_of_proteome": percent,
            }
            for name, count, percent in result.coverage
        ),
    )
    metadata = (
        ("mapping_rule_version", PSORTB_MAPPING_RULE_VERSION),
        ("source_role", PSORTB_SOURCE_ROLE),
        ("raw_input_sha256", result.raw_input_sha256),
        ("formal_output_sha256", formal_hash),
        ("audit_output_sha256", _sha256(audit_output)),
        ("coverage_output_sha256", _sha256(coverage_output)),
        ("query_prediction", result.query_prediction),
    )
    _write_dicts(
        metadata_output,
        ("key", "value"),
        ({"key": key, "value": value} for key, value in metadata),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--proteome-fasta", type=Path)
    parser.add_argument("--base-annotation-table", type=Path)
    parser.add_argument("--query-id", default="")
    parser.add_argument("--source-version", default="")
    parser.add_argument("--source-command", default="")
    parser.add_argument("--organism-mode", default="archaea")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    result = convert_psortb(
        input_path=args.input,
        proteome_fasta=args.proteome_fasta,
        base_annotation_table=args.base_annotation_table,
        query_id=args.query_id,
        source_version=args.source_version,
        source_command=args.source_command,
        organism_mode=args.organism_mode,
    )
    write_conversion_outputs(
        result,
        annotation_output=args.output,
        audit_output=args.audit_output,
        coverage_output=args.coverage_output,
        metadata_output=args.metadata_output,
    )
    print(f"Formal annotation rows: {len(result.formal_rows)}")
    print(f"Audited PSORTb rows: {len(result.audit_rows)}")
    print(f"Exact duplicate rows: {result.exact_duplicate_count}")
    print(f"Unknown protein IDs: {result.unknown_protein_id_count}")
    print(f"Missing protein IDs: {result.missing_protein_id_count}")
    print(f"Query prediction: {result.query_prediction or 'not requested'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
