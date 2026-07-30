"""Validate and normalize panel proteomes for deterministic orthology analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ALLOWED_RESIDUES = set("ACDEFGHIKLMNPQRSTVWYBJZXUO")
STANDARD_RESIDUES = set("ACDEFGHIKLMNPQRSTVWY")
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")

MAPPING_COLUMNS = (
    "species_id",
    "assembly_accession",
    "normalized_protein_id",
    "raw_protein_id",
    "original_header",
    "sequence_length",
    "sequence_sha256",
    "nonstandard_residues",
    "duplicate_sequence_group",
    "input_fasta_sha256",
    "normalized_fasta_sha256",
)

AUDIT_COLUMNS = (
    "species_id",
    "assembly_accession",
    "input_path",
    "input_sha256",
    "output_path",
    "output_sha256",
    "protein_count",
    "duplicate_id_count",
    "duplicate_sequence_group_count",
    "duplicate_sequence_protein_count",
    "empty_sequence_count",
    "invalid_residue_sequence_count",
    "nonstandard_residue_sequence_count",
    "internal_stop_sequence_count",
    "extremely_short_sequence_count",
    "extremely_long_sequence_count",
    "minimum_length",
    "maximum_length",
    "query_protein_count",
    "status",
    "warnings",
)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header.split()[0], header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
            if not header:
                raise ValueError(f"Empty FASTA header on line {line_number}: {path}")
        elif line:
            if header is None:
                raise ValueError(f"Sequence before FASTA header on line {line_number}: {path}")
            chunks.append(line)
    if header is not None:
        records.append((header.split()[0], header, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records: {path}")
    return records


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def prepare_proteomes(
    *,
    panel_path: Path,
    data_root: Path,
    normalized_directory: Path,
    mapping_output: Path,
    audit_output: Path,
    manifest_output: Path,
    query_assembly: str,
    query_protein_id: str,
    expected_query_fasta_sha256: str,
) -> dict[str, Any]:
    with panel_path.open(encoding="utf-8", newline="") as handle:
        panel_rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"species_id", "assembly_accession", "selection_status"}
    missing = required - set(panel_rows[0] if panel_rows else ())
    if missing:
        raise ValueError("Panel table missing columns: " + ", ".join(sorted(missing)))
    selected = [row for row in panel_rows if row["selection_status"] == "selected"]
    if not selected:
        raise ValueError("Panel contains no selected assemblies")
    if len({row["assembly_accession"] for row in selected}) != len(selected):
        raise ValueError("Panel contains duplicate selected assemblies")
    normalized_directory.mkdir(parents=True, exist_ok=True)
    mapping_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    total_proteins = 0
    query_count = 0
    for panel_row in sorted(selected, key=lambda row: row["species_id"]):
        species_id = panel_row["species_id"]
        assembly = panel_row["assembly_accession"]
        if not SAFE_ID_RE.fullmatch(species_id):
            raise ValueError(f"Unsafe species_id: {species_id!r}")
        input_path = data_root / assembly / "protein.faa"
        if not input_path.is_file():
            raise ValueError(f"Missing protein FASTA: {input_path}")
        input_sha = _sha256_file(input_path)
        if assembly == query_assembly and input_sha != expected_query_fasta_sha256:
            raise ValueError(
                f"Query proteome checksum mismatch: {input_sha} != {expected_query_fasta_sha256}"
            )
        records = _read_fasta(input_path)
        id_counts = Counter(record[0] for record in records)
        duplicate_ids = sorted(protein_id for protein_id, count in id_counts.items() if count > 1)
        if duplicate_ids:
            raise ValueError(f"Duplicate FASTA IDs in {assembly}: {', '.join(duplicate_ids[:10])}")
        empty = [protein_id for protein_id, _, sequence in records if not sequence]
        if empty:
            raise ValueError(f"Empty FASTA sequences in {assembly}: {', '.join(empty[:10])}")
        invalid = [
            protein_id
            for protein_id, _, sequence in records
            if set(sequence) - ALLOWED_RESIDUES - {"*"}
        ]
        if invalid:
            raise ValueError(f"Invalid FASTA residues in {assembly}: {', '.join(invalid[:10])}")
        internal_stop = [
            protein_id for protein_id, _, sequence in records if "*" in sequence.rstrip("*")
        ]
        if internal_stop:
            raise ValueError(
                f"Internal stop residues in {assembly}: {', '.join(internal_stop[:10])}"
            )
        terminal_stop = [
            protein_id for protein_id, _, sequence in records if sequence.endswith("*")
        ]
        if terminal_stop:
            raise ValueError(
                f"Terminal stop residues require explicit policy in {assembly}: "
                + ", ".join(terminal_stop[:10])
            )
        by_sequence: dict[str, list[str]] = defaultdict(list)
        for protein_id, _, sequence in records:
            by_sequence[sequence].append(protein_id)
        duplicate_groups = {
            sequence: f"{species_id}_dupseq_{index:06d}"
            for index, sequence in enumerate(
                sorted(sequence for sequence, ids in by_sequence.items() if len(ids) > 1), 1
            )
        }
        output_path = normalized_directory / f"{species_id}.faa"
        output_parts: list[str] = []
        species_mapping: list[dict[str, Any]] = []
        nonstandard_count = 0
        short_count = 0
        long_count = 0
        species_query_count = 0
        for raw_id, original_header, sequence in sorted(records, key=lambda record: record[0]):
            if not SAFE_ID_RE.fullmatch(raw_id):
                raise ValueError(f"Unsafe raw protein ID in {assembly}: {raw_id!r}")
            normalized_id = f"{species_id}__{raw_id}"
            nonstandard = "".join(sorted(set(sequence) - STANDARD_RESIDUES))
            nonstandard_count += bool(nonstandard)
            short_count += len(sequence) < 20
            long_count += len(sequence) > 10_000
            species_query_count += raw_id == query_protein_id
            output_parts.append(f">{normalized_id}\n")
            output_parts.extend(
                sequence[index : index + 60] + "\n" for index in range(0, len(sequence), 60)
            )
            species_mapping.append(
                {
                    "species_id": species_id,
                    "assembly_accession": assembly,
                    "normalized_protein_id": normalized_id,
                    "raw_protein_id": raw_id,
                    "original_header": original_header,
                    "sequence_length": len(sequence),
                    "sequence_sha256": _sha256_bytes(sequence.encode("ascii")),
                    "nonstandard_residues": nonstandard,
                    "duplicate_sequence_group": duplicate_groups.get(sequence, ""),
                    "input_fasta_sha256": input_sha,
                    "normalized_fasta_sha256": "",
                }
            )
        output_path.write_text("".join(output_parts), encoding="utf-8", newline="\n")
        output_sha = _sha256_file(output_path)
        for row in species_mapping:
            row["normalized_fasta_sha256"] = output_sha
        mapping_rows.extend(species_mapping)
        lengths = [len(record[2]) for record in records]
        duplicate_proteins = sum(len(ids) for ids in by_sequence.values() if len(ids) > 1)
        warnings = []
        if duplicate_groups:
            warnings.append("exact_duplicate_sequences")
        if nonstandard_count:
            warnings.append("nonstandard_residues")
        if short_count:
            warnings.append("extremely_short_proteins")
        if long_count:
            warnings.append("extremely_long_proteins")
        audit_rows.append(
            {
                "species_id": species_id,
                "assembly_accession": assembly,
                "input_path": str(input_path.resolve()),
                "input_sha256": input_sha,
                "output_path": str(output_path.resolve()),
                "output_sha256": output_sha,
                "protein_count": len(records),
                "duplicate_id_count": 0,
                "duplicate_sequence_group_count": len(duplicate_groups),
                "duplicate_sequence_protein_count": duplicate_proteins,
                "empty_sequence_count": 0,
                "invalid_residue_sequence_count": 0,
                "nonstandard_residue_sequence_count": nonstandard_count,
                "internal_stop_sequence_count": 0,
                "extremely_short_sequence_count": short_count,
                "extremely_long_sequence_count": long_count,
                "minimum_length": min(lengths),
                "maximum_length": max(lengths),
                "query_protein_count": species_query_count,
                "status": "valid",
                "warnings": "|".join(warnings),
            }
        )
        total_proteins += len(records)
        query_count += species_query_count if assembly == query_assembly else 0
    if query_count != 1:
        raise ValueError(f"Query protein must occur once in query assembly; observed {query_count}")
    mapping_rows.sort(key=lambda row: (row["species_id"], row["raw_protein_id"]))
    audit_rows.sort(key=lambda row: row["species_id"])
    _write_tsv(mapping_output, MAPPING_COLUMNS, mapping_rows)
    _write_tsv(audit_output, AUDIT_COLUMNS, audit_rows)
    manifest = {
        "schema_version": "1.0",
        "normalization_rule_version": "ma4115-orthology-fasta-v1",
        "query_assembly": query_assembly,
        "query_protein_id": query_protein_id,
        "expected_query_fasta_sha256": expected_query_fasta_sha256,
        "panel_path": str(panel_path.resolve()),
        "panel_sha256": _sha256_file(panel_path),
        "species_count": len(audit_rows),
        "total_protein_count": total_proteins,
        "mapping_path": str(mapping_output.resolve()),
        "mapping_sha256": _sha256_file(mapping_output),
        "audit_path": str(audit_output.resolve()),
        "audit_sha256": _sha256_file(audit_output),
        "normalized_directory": str(normalized_directory.resolve()),
        "normalized_fasta_checksums": {
            row["species_id"]: row["output_sha256"] for row in audit_rows
        },
    }
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--normalized-directory", required=True, type=Path)
    parser.add_argument("--mapping-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--manifest-output", required=True, type=Path)
    parser.add_argument("--query-assembly", required=True)
    parser.add_argument("--query-protein-id", required=True)
    parser.add_argument("--expected-query-fasta-sha256", required=True)
    args = parser.parse_args()
    manifest = prepare_proteomes(
        panel_path=args.panel,
        data_root=args.data_root,
        normalized_directory=args.normalized_directory,
        mapping_output=args.mapping_output,
        audit_output=args.audit_output,
        manifest_output=args.manifest_output,
        query_assembly=args.query_assembly,
        query_protein_id=args.query_protein_id,
        expected_query_fasta_sha256=args.expected_query_fasta_sha256,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
