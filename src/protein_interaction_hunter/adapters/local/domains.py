"""Load and validate local domain annotation TSV files."""

import csv
from pathlib import Path

from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.domain import DomainAnnotationRecord

DOMAIN_COLUMNS = (
    "protein_id",
    "source",
    "accession",
    "name",
    "start",
    "end",
    "architecture_index",
)


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


class LocalDomainTsvLoader:
    def load(self, path: Path) -> list[DomainAnnotationRecord]:
        domain_path = path.expanduser().resolve()

        if not domain_path.is_file():
            raise InputValidationError(
                f"Domain annotation TSV not found: {domain_path}"
            )

        records: list[DomainAnnotationRecord] = []
        seen: set[tuple[str, str, str, int, int]] = set()

        with domain_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(DOMAIN_COLUMNS) - set(reader.fieldnames or [])

            if missing:
                raise InputValidationError(
                    "Domain annotation TSV missing columns: "
                    + ", ".join(sorted(missing))
                )

            for line_number, row in enumerate(reader, start=2):
                try:
                    protein_id = (row.get("protein_id") or "").strip()
                    source = (row.get("source") or "").strip()
                    accession = (row.get("accession") or "").strip()
                    start = int((row.get("start") or "").strip())
                    end = int((row.get("end") or "").strip())
                    architecture_index = int(
                        (row.get("architecture_index") or "").strip()
                    )

                    record = DomainAnnotationRecord(
                        protein_id=protein_id,
                        source=source,
                        accession=accession,
                        name=_optional(row.get("name")),
                        start=start,
                        end=end,
                        architecture_index=architecture_index,
                    )
                except (ValueError, ValidationError) as exc:
                    raise InputValidationError(
                        f"Invalid domain annotation on line "
                        f"{line_number}: {exc}"
                    ) from exc

                identity = (
                    record.protein_id,
                    record.source,
                    record.accession,
                    record.start,
                    record.end,
                )

                if identity in seen:
                    raise InputValidationError(
                        "Duplicate domain annotation on line "
                        f"{line_number}: {identity}"
                    )

                seen.add(identity)
                records.append(record)

        return records