"""Strict GFF3 coordinate and sequence-region loader."""

from pathlib import Path
from urllib.parse import unquote

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.genome import GeneCoordinate, GffDocument, SequenceRegion


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
    """Load GFF features plus optional contig boundary metadata."""

    def load(self, path: Path) -> list[GeneCoordinate]:
        return self.load_document(path).features

    def load_document(self, path: Path) -> GffDocument:
        gff_path = path.expanduser().resolve()
        if not gff_path.is_file():
            raise InputValidationError(f"Genome GFF not found: {gff_path}")
        records: list[GeneCoordinate] = []
        sequence_regions: dict[str, SequenceRegion] = {}
        warnings: list[str] = []
        for line_number, raw_line in enumerate(
            gff_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if raw_line.startswith("##sequence-region"):
                parts = raw_line.split()
                if len(parts) != 4:
                    raise InputValidationError(
                        f"GFF line {line_number} has malformed ##sequence-region"
                    )
                try:
                    region = SequenceRegion(seqid=parts[1], start=int(parts[2]), end=int(parts[3]))
                except ValueError as exc:
                    raise InputValidationError(
                        f"Invalid GFF sequence-region on line {line_number}: {exc}"
                    ) from exc
                if region.seqid in sequence_regions and sequence_regions[region.seqid] != region:
                    raise InputValidationError(f"Conflicting ##sequence-region for {region.seqid}")
                sequence_regions[region.seqid] = region
                continue
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise InputValidationError(f"GFF line {line_number} must contain exactly 9 columns")
            attributes = parse_attributes(columns[8])
            try:
                start = int(columns[3])
                end = int(columns[4])
            except ValueError as exc:
                raise InputValidationError(
                    f"GFF line {line_number} has non-integer coordinates"
                ) from exc
            strand = columns[6] if columns[6] in {"+", "-"} else "?"
            parent_ids = attributes.get("Parent", [])
            try:
                record = GeneCoordinate(
                    seqid=columns[0],
                    feature_type=columns[2],
                    source=columns[1] if columns[1] != "." else None,
                    start=start,
                    end=end,
                    strand=strand,
                    feature_id=_first(attributes, "ID"),
                    parent_id=parent_ids[0] if parent_ids else None,
                    parent_ids=parent_ids,
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
        for record in records:
            feature_region = sequence_regions.get(record.seqid)
            if feature_region and (
                record.start < feature_region.start or record.end > feature_region.end
            ):
                raise InputValidationError(
                    f"GFF feature {record.feature_id or record.seqid} lies outside sequence-region"
                )
        for seqid in sorted({record.seqid for record in records} - sequence_regions.keys()):
            warnings.append(f"missing_sequence_region:{seqid}")
        return GffDocument(
            features=records,
            sequence_regions=sequence_regions,
            warnings=warnings,
        )


def coordinates_by_protein(records: list[GeneCoordinate]) -> dict[str, GeneCoordinate]:
    """Index direct protein_id attributes; ambiguous duplicate IDs are invalid."""
    index: dict[str, GeneCoordinate] = {}
    for record in records:
        if record.protein_id is None:
            continue
        if record.protein_id in index:
            raise InputValidationError(f"Duplicate GFF protein_id coordinate: {record.protein_id}")
        index[record.protein_id] = record
    return index
