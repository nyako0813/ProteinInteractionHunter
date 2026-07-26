"""Load and validate deterministic local phylogenetic profile TSV files."""

import csv
from pathlib import Path

from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.phylogenetic_profile import (
    PhylogeneticProfileObservation,
)

REQUIRED_PHYLOGENETIC_PROFILE_COLUMNS = ("protein_id", "species_id", "presence")


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


def _presence(value: str | None) -> bool | None:
    normalized = (value or "").strip().casefold()
    if not normalized:
        return None
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError("presence must be true/false, 1/0, yes/no, or blank")


class LocalPhylogeneticProfileTsvLoader:
    def load(self, path: Path) -> list[PhylogeneticProfileObservation]:
        profile_path = path.expanduser().resolve()
        if not profile_path.is_file():
            raise InputValidationError(f"Phylogenetic profile TSV not found: {profile_path}")

        records: list[PhylogeneticProfileObservation] = []
        seen: set[tuple[str, str]] = set()

        with profile_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(REQUIRED_PHYLOGENETIC_PROFILE_COLUMNS) - set(reader.fieldnames or [])
            if missing:
                raise InputValidationError(
                    "Phylogenetic profile TSV missing columns: " + ", ".join(sorted(missing))
                )

            for line_number, row in enumerate(reader, start=2):
                try:
                    record = PhylogeneticProfileObservation(
                        protein_id=(row.get("protein_id") or "").strip(),
                        species_id=(row.get("species_id") or "").strip(),
                        presence=_presence(row.get("presence")),
                        taxonomic_group=_optional(row.get("taxonomic_group")),
                        source=_optional(row.get("source")),
                        source_record_id=_optional(row.get("source_record_id")),
                    )
                except (ValueError, ValidationError) as exc:
                    raise InputValidationError(
                        f"Invalid phylogenetic profile observation on line {line_number}: {exc}"
                    ) from exc

                identity = (record.protein_id, record.species_id)
                if identity in seen:
                    raise InputValidationError(
                        "Duplicate phylogenetic profile observation on line "
                        f"{line_number}: {identity}"
                    )
                seen.add(identity)
                records.append(record)

        return records
