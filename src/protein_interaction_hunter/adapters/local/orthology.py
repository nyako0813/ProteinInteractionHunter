"""Load and validate local orthology annotation TSV files."""

import csv
from pathlib import Path

from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import (
    EvidenceOrigin,
    EvidenceStatus,
)
from protein_interaction_hunter.models.evidence import OrthologRecord

ORTHOLOGY_COLUMNS = (
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


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


def _optional_float(value: str | None) -> float | None:
    stripped = (value or "").strip()
    return float(stripped) if stripped else None


def _parse_bool(value: str | None) -> bool:
    normalized = (value or "").strip().casefold()

    if normalized in {"", "false", "0", "no"}:
        return False
    if normalized in {"true", "1", "yes"}:
        return True

    raise ValueError("paralog_ambiguity must be true/false, 1/0, yes/no, or blank")


class LocalOrthologyTsvLoader:
    def load(self, path: Path) -> list[OrthologRecord]:
        orthology_path = path.expanduser().resolve()

        if not orthology_path.is_file():
            raise InputValidationError(f"Orthology TSV not found: {orthology_path}")

        records: list[OrthologRecord] = []
        seen: set[tuple[str, str, str]] = set()

        with orthology_path.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(ORTHOLOGY_COLUMNS) - set(reader.fieldnames or [])

            if missing:
                raise InputValidationError(
                    "Orthology TSV missing columns: " + ", ".join(sorted(missing))
                )

            for line_number, row in enumerate(reader, start=2):
                try:
                    record = OrthologRecord(
                        status=EvidenceStatus.AVAILABLE,
                        origin=EvidenceOrigin.ORTHOLOG_TRANSFERRED,
                        protein_id=(row.get("protein_id") or "").strip(),
                        reference_id=(row.get("reference_id") or "").strip(),
                        ortholog_id=_optional(row.get("ortholog_id")),
                        reference_organism=_optional(row.get("reference_organism")),
                        identity=_optional_float(row.get("identity")),
                        query_coverage=_optional_float(row.get("query_coverage")),
                        subject_coverage=_optional_float(row.get("subject_coverage")),
                        evalue=_optional_float(row.get("evalue")),
                        orthogroup=_optional(row.get("orthogroup")),
                        relationship=_optional(row.get("relationship")),
                        paralog_ambiguity=_parse_bool(row.get("paralog_ambiguity")),
                        source=_optional(row.get("source")),
                        source_record_id=_optional(row.get("source_record_id")),
                    )
                except (ValueError, ValidationError) as exc:
                    raise InputValidationError(
                        f"Invalid orthology annotation on line {line_number}: {exc}"
                    ) from exc

                identity = (
                    record.protein_id,
                    record.reference_id,
                    record.ortholog_id or "",
                )

                if identity in seen:
                    raise InputValidationError(
                        f"Duplicate orthology annotation on line {line_number}: {identity}"
                    )

                seen.add(identity)
                records.append(record)

        return records
