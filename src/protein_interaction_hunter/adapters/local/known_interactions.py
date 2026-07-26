"""Load validated local known-interaction observations from TSV."""

import csv
import re
from pathlib import Path
from typing import cast

from pydantic import ValidationError

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.enums import EvidenceOrigin, EvidenceStatus
from protein_interaction_hunter.models.evidence import (
    EvidenceProvenance,
    IdentifierMappingStatus,
    KnownInteractionObservation,
    KnownInteractionType,
)

KNOWN_INTERACTION_REQUIRED_COLUMNS = (
    "protein_a_id",
    "protein_b_id",
    "interaction_type",
    "reference_organism",
    "source",
    "source_record_id",
)

_METHOD_ALIASES = {
    "y2h": "yeast_two_hybrid",
    "yeast_2_hybrid": "yeast_two_hybrid",
    "yeast_two_hybrid": "yeast_two_hybrid",
    "pulldown": "pull_down",
    "pull_down": "pull_down",
    "co_ip": "co_immunoprecipitation",
    "co_immunoprecipitation": "co_immunoprecipitation",
    "ap_ms": "affinity_purification_mass_spectrometry",
    "affinity_purification_mass_spectrometry": ("affinity_purification_mass_spectrometry"),
    "xl_ms": "crosslinking_mass_spectrometry",
    "crosslinking_mass_spectrometry": "crosslinking_mass_spectrometry",
    "spr": "surface_plasmon_resonance",
    "surface_plasmon_resonance": "surface_plasmon_resonance",
    "itc": "isothermal_titration_calorimetry",
    "isothermal_titration_calorimetry": "isothermal_titration_calorimetry",
    "bli": "biolayer_interferometry",
    "biolayer_interferometry": "biolayer_interferometry",
    "biotin_proximity_labeling": "proximity_labeling",
    "proximity_labeling": "proximity_labeling",
    "genetic_interaction": "genetic_interaction",
    "co_expression": "co_expression",
    "database_inference": "database_inference",
}


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


def _optional_float(value: str | None) -> float | None:
    stripped = (value or "").strip()
    return float(stripped) if stripped else None


def _optional_bool(value: str | None) -> bool | None:
    normalized = (value or "").strip().casefold()
    if not normalized:
        return None
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError("boolean must be true/false, 1/0, yes/no, or blank")


def normalize_detection_method(value: str | None) -> str | None:
    original = _optional(value)
    if original is None:
        return None
    slug = re.sub(r"[^a-z0-9]+", "_", original.casefold()).strip("_")
    return _METHOD_ALIASES.get(slug, "other")


class LocalKnownInteractionsTsvLoader:
    def load(self, path: Path) -> list[KnownInteractionObservation]:
        interaction_path = path.expanduser().resolve()
        if not interaction_path.is_file():
            raise InputValidationError(f"Known interactions TSV not found: {interaction_path}")

        records: list[KnownInteractionObservation] = []
        seen_source_records: set[tuple[str, str]] = set()
        with interaction_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="	")
            missing = set(KNOWN_INTERACTION_REQUIRED_COLUMNS) - set(reader.fieldnames or [])
            if missing:
                raise InputValidationError(
                    "Known interactions TSV missing columns: " + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                source = (row.get("source") or "").strip()
                source_record_id = (row.get("source_record_id") or "").strip()
                try:
                    record = KnownInteractionObservation(
                        status=EvidenceStatus.AVAILABLE,
                        origin=EvidenceOrigin.EXACT_PAIR,
                        protein_a_id=(row.get("protein_a_id") or "").strip(),
                        protein_b_id=(row.get("protein_b_id") or "").strip(),
                        interaction_type=cast(
                            KnownInteractionType,
                            (row.get("interaction_type") or "").strip(),
                        ),
                        reference_organism=(row.get("reference_organism") or "").strip(),
                        detection_method=_optional(row.get("detection_method")),
                        normalized_detection_method=normalize_detection_method(
                            row.get("detection_method")
                        ),
                        publication_id=_optional(row.get("publication_id")),
                        confidence=_optional_float(row.get("confidence")),
                        is_direct=_optional_bool(row.get("is_direct")),
                        is_physical=_optional_bool(row.get("is_physical")),
                        is_biological=_optional_bool(row.get("is_biological")),
                        database_version=_optional(row.get("database_version")),
                        protein_a_reference_id=_optional(row.get("protein_a_reference_id")),
                        protein_b_reference_id=_optional(row.get("protein_b_reference_id")),
                        source=source,
                        source_record_id=source_record_id,
                        notes=_optional(row.get("notes")),
                        identifier_mapping_status=cast(
                            IdentifierMappingStatus,
                            _optional(row.get("identifier_mapping_status")) or "mapped",
                        ),
                        provenance=[
                            EvidenceProvenance(
                                source_name=source or "local_known_interactions_table",
                                source_version=_optional(row.get("database_version")),
                                source_record_id=source_record_id or None,
                                method="validated_local_interaction_tsv_row",
                            )
                        ],
                    )
                except (ValueError, ValidationError) as exc:
                    raise InputValidationError(
                        f"Invalid known interaction observation on line {line_number}: {exc}"
                    ) from exc
                identity = (record.source, record.source_record_id)
                if identity in seen_source_records:
                    raise InputValidationError(
                        f"Duplicate known interaction source record on line "
                        f"{line_number}: {identity}"
                    )
                seen_source_records.add(identity)
                records.append(record)
        return records
