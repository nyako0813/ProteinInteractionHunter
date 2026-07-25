"""Local annotation TSV validation."""

import csv
from pathlib import Path

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.annotation import AnnotationRecord
from protein_interaction_hunter.models.enums import EvidenceStatus

ANNOTATION_COLUMNS = (
    "protein_id",
    "gene_name",
    "locus_tag",
    "product",
    "functional_category",
    "localization_annotation",
    "transmembrane_annotation",
    "annotation_source",
    "annotation_confidence",
)


def _optional(value: str | None) -> str | None:
    stripped = (value or "").strip()
    return stripped or None


class LocalAnnotationTsvLoader:
    def load(self, path: Path) -> list[AnnotationRecord]:
        annotation_path = path.expanduser().resolve()
        if not annotation_path.is_file():
            raise InputValidationError(f"Annotation TSV not found: {annotation_path}")
        records: list[AnnotationRecord] = []
        seen: set[str] = set()
        with annotation_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(ANNOTATION_COLUMNS) - set(reader.fieldnames or [])
            if missing:
                raise InputValidationError(
                    "Annotation TSV missing columns: " + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                protein_id = (row.get("protein_id") or "").strip()
                if not protein_id:
                    raise InputValidationError(
                        f"Annotation TSV has empty protein_id on line {line_number}"
                    )
                if protein_id in seen:
                    raise InputValidationError(
                        f"Duplicate annotation protein_id: {protein_id}"
                    )
                seen.add(protein_id)
                confidence_text = _optional(row.get("annotation_confidence"))
                try:
                    confidence = (
                        float(confidence_text) if confidence_text is not None else None
                    )
                    has_annotation = any(
                        _optional(row.get(column)) for column in ANNOTATION_COLUMNS[1:]
                    )
                    records.append(
                        AnnotationRecord(
                            protein_id=protein_id,
                            gene_name=_optional(row.get("gene_name")),
                            locus_tag=_optional(row.get("locus_tag")),
                            product=_optional(row.get("product")),
                            functional_category=_optional(
                                row.get("functional_category")
                            ),
                            localization_annotation=_optional(
                                row.get("localization_annotation")
                            ),
                            transmembrane_annotation=_optional(
                                row.get("transmembrane_annotation")
                            ),
                            annotation_source=_optional(row.get("annotation_source")),
                            annotation_confidence=confidence,
                            status=(
                                EvidenceStatus.AVAILABLE
                                if has_annotation
                                else EvidenceStatus.MISSING
                            ),
                        )
                    )
                except ValueError as exc:
                    raise InputValidationError(
                        f"Invalid annotation on line {line_number}: {exc}"
                    ) from exc
        return records
