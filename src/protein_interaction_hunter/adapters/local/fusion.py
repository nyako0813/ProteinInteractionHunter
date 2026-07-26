"""Load validated local gene-fusion observations from TSV."""

import csv
from pathlib import Path

from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import EvidenceProvenance, FusionObservation

FUSION_REQUIRED_COLUMNS = (
    "query_protein_id",
    "candidate_protein_id",
    "fusion_protein_id",
    "reference_organism",
    "query_component_start",
    "query_component_end",
    "candidate_component_start",
    "candidate_component_end",
    "fusion_protein_length",
)


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


def _integer(value: str | None, field: str) -> int:
    stripped = (value or "").strip()
    if not stripped:
        raise ValueError(f"{field} must not be blank")
    return int(stripped)


def _optional_float(value: str | None) -> float | None:
    stripped = (value or "").strip()
    return float(stripped) if stripped else None


def _duplicate_identity(record: FusionObservation) -> tuple[object, ...]:
    first, second = sorted((record.query_protein_id, record.candidate_protein_id))
    if record.query_protein_id == first:
        first_region = (record.query_component_start, record.query_component_end)
        second_region = (record.candidate_component_start, record.candidate_component_end)
    else:
        first_region = (record.candidate_component_start, record.candidate_component_end)
        second_region = (record.query_component_start, record.query_component_end)
    return (
        first,
        second,
        record.reference_organism,
        record.fusion_protein_id,
        first_region,
        second_region,
    )


class LocalFusionTsvLoader:
    """Load 1-based inclusive fusion-component coordinates."""

    def load(self, path: Path) -> list[FusionObservation]:
        fusion_path = path.expanduser().resolve()
        if not fusion_path.is_file():
            raise InputValidationError(f"Fusion TSV not found: {fusion_path}")

        records: list[FusionObservation] = []
        seen: set[tuple[object, ...]] = set()
        with fusion_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            missing = set(FUSION_REQUIRED_COLUMNS) - set(reader.fieldnames or [])
            if missing:
                raise InputValidationError(
                    "Fusion TSV missing columns: " + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                source_record_id = _optional(row.get("source_record_id"))
                try:
                    record = FusionObservation(
                        status=EvidenceStatus.AVAILABLE,
                        origin=EvidenceOrigin.INFERRED,
                        query_protein_id=(row.get("query_protein_id") or "").strip(),
                        candidate_protein_id=(row.get("candidate_protein_id") or "").strip(),
                        fusion_protein_id=(row.get("fusion_protein_id") or "").strip(),
                        reference_organism=(row.get("reference_organism") or "").strip(),
                        query_component_reference_id=_optional(
                            row.get("query_component_reference_id")
                        ),
                        candidate_component_reference_id=_optional(
                            row.get("candidate_component_reference_id")
                        ),
                        query_component_start=_integer(
                            row.get("query_component_start"), "query_component_start"
                        ),
                        query_component_end=_integer(
                            row.get("query_component_end"), "query_component_end"
                        ),
                        candidate_component_start=_integer(
                            row.get("candidate_component_start"),
                            "candidate_component_start",
                        ),
                        candidate_component_end=_integer(
                            row.get("candidate_component_end"), "candidate_component_end"
                        ),
                        fusion_protein_length=_integer(
                            row.get("fusion_protein_length"), "fusion_protein_length"
                        ),
                        query_component_coverage=_optional_float(
                            row.get("query_component_coverage")
                        ),
                        candidate_component_coverage=_optional_float(
                            row.get("candidate_component_coverage")
                        ),
                        query_component_identity=_optional_float(
                            row.get("query_component_identity")
                        ),
                        candidate_component_identity=_optional_float(
                            row.get("candidate_component_identity")
                        ),
                        evalue_query=_optional_float(row.get("evalue_query")),
                        evalue_candidate=_optional_float(row.get("evalue_candidate")),
                        source=_optional(row.get("source")),
                        source_record_id=source_record_id,
                        provenance=[
                            EvidenceProvenance(
                                source_name="local_fusion_table",
                                source_record_id=source_record_id,
                                method="validated_1_based_inclusive_tsv_row",
                            )
                        ],
                    )
                except (ValueError, ValidationError) as exc:
                    raise InputValidationError(
                        f"Invalid fusion observation on line {line_number}: {exc}"
                    ) from exc
                identity = _duplicate_identity(record)
                if identity in seen:
                    raise InputValidationError(
                        f"Duplicate fusion observation on line {line_number}: {identity}"
                    )
                seen.add(identity)
                records.append(record)
        return records
