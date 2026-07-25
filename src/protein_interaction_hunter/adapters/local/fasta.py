"""Small standard-library FASTA loader for validated MVP-0 fixtures."""

import re
from pathlib import Path

from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.protein import ProteinRecord

_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
_METADATA_PATTERN = re.compile(r"\[(gene_id|locus_tag)=([^\]]+)\]")


class LocalFastaLoader:
    """Load protein FASTA without biological analysis."""

    def load(self, path: Path) -> list[ProteinRecord]:
        fasta_path = path.expanduser().resolve()
        if not fasta_path.is_file():
            raise InputValidationError(f"Proteome FASTA not found: {fasta_path}")

        records: list[ProteinRecord] = []
        seen_ids: set[str] = set()
        header: str | None = None
        sequence_parts: list[str] = []

        def finish_record() -> None:
            if header is None:
                return
            record = self._build_record(header, sequence_parts)
            if record.protein_id in seen_ids:
                raise InputValidationError(f"Duplicate FASTA ID: {record.protein_id}")
            seen_ids.add(record.protein_id)
            records.append(record)

        for line_number, raw_line in enumerate(
            fasta_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                sequence_parts = []
                if not header:
                    raise InputValidationError(f"Empty FASTA header on line {line_number}")
            else:
                if header is None:
                    raise InputValidationError(
                        f"FASTA sequence appears before a header on line {line_number}"
                    )
                sequence_parts.append("".join(line.split()).upper())
        finish_record()
        if not records:
            raise InputValidationError(f"No FASTA records found: {fasta_path}")
        return records

    @staticmethod
    def _build_record(header: str, sequence_parts: list[str]) -> ProteinRecord:
        first, _, remainder = header.partition(" ")
        protein_id = first.strip()
        if not protein_id:
            raise InputValidationError("FASTA protein ID must not be empty")
        sequence = "".join(sequence_parts).upper()
        if not sequence:
            raise InputValidationError(f"Empty FASTA sequence for {protein_id}")
        invalid = sorted(set(sequence) - _AMINO_ACIDS)
        if invalid:
            raise InputValidationError(
                f"Invalid amino-acid characters for {protein_id}: {''.join(invalid)}"
            )
        metadata = dict(_METADATA_PATTERN.findall(remainder))
        description = _METADATA_PATTERN.sub("", remainder).strip()
        return ProteinRecord(
            protein_id=protein_id,
            description=description,
            sequence=sequence,
            gene_id=metadata.get("gene_id"),
            locus_tag=metadata.get("locus_tag"),
        )


def duplicate_sequence_groups(records: list[ProteinRecord]) -> list[list[str]]:
    """Return deterministic groups of IDs sharing an exact sequence."""
    by_sequence: dict[str, list[str]] = {}
    for record in records:
        by_sequence.setdefault(record.sequence, []).append(record.protein_id)
    return [
        sorted(ids)
        for _, ids in sorted(by_sequence.items(), key=lambda item: item[0])
        if len(ids) > 1
    ]
