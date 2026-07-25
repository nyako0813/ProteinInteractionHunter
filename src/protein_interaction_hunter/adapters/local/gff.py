"""Minimal GFF3 parser that stops before neighborhood analysis."""

from pathlib import Path
from urllib.parse import unquote

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.genome import GeneCoordinate


def parse_attributes(text: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = {}
    for part in text.split(";"):
        if not part.strip():
            continue
        key, separator, raw_value = part.partition("=")
        value = unquote(raw_value if separator else "")
        attributes.setdefault(key.strip(), []).extend(
            item.strip() for item in value.split(",") if item.strip()
        )
    return attributes


def _first(attributes: dict[str, list[str]], key: str) -> str | None:
    values = attributes.get(key, [])
    return values[0] if values else None


class LocalGff3Loader:
    """Load gene/CDS coordinates and identifier attributes."""

    def load(self, path: Path) -> list[GeneCoordinate]:
        gff_path = path.expanduser().resolve()
        if not gff_path.is_file():
            raise InputValidationError(f"Genome GFF not found: {gff_path}")
        records: list[GeneCoordinate] = []
        for line_number, raw_line in enumerate(
            gff_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise InputValidationError(
                    f"GFF line {line_number} must contain exactly 9 columns"
                )
            attributes = parse_attributes(columns[8])
            try:
                start = int(columns[3])
                end = int(columns[4])
            except ValueError as exc:
                raise InputValidationError(
                    f"GFF line {line_number} has non-integer coordinates"
                ) from exc
            strand = columns[6] if columns[6] in {"+", "-"} else None
            try:
                record = GeneCoordinate(
                    seqid=columns[0],
                    feature_type=columns[2],
                    start=start,
                    end=end,
                    strand=strand,
                    feature_id=_first(attributes, "ID"),
                    parent_id=_first(attributes, "Parent"),
                    protein_id=_first(attributes, "protein_id"),
                    locus_tag=_first(attributes, "locus_tag"),
                    old_locus_tag=_first(attributes, "old_locus_tag"),
                    attributes=attributes,
                )
            except ValueError as exc:
                raise InputValidationError(f"Invalid GFF line {line_number}: {exc}") from exc
            records.append(record)
        if not records:
            raise InputValidationError(f"No GFF features found: {gff_path}")
        return records


def coordinates_by_protein(records: list[GeneCoordinate]) -> dict[str, GeneCoordinate]:
    """Index direct protein_id attributes; ambiguous duplicate IDs are invalid."""
    index: dict[str, GeneCoordinate] = {}
    for record in records:
        if record.protein_id is None:
            continue
        if record.protein_id in index:
            raise InputValidationError(
                f"Duplicate GFF protein_id coordinate: {record.protein_id}"
            )
        index[record.protein_id] = record
    return index
