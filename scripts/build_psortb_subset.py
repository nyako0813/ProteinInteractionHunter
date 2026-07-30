#!/usr/bin/env python3
"""Build a deterministic accession-only FASTA subset for a PSORTb pilot."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

from protein_interaction_hunter.adapters.local.domains import LocalDomainTsvLoader
from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.protein import ProteinRecord

AUDIT_COLUMNS = ("protein_id", "length_aa", "description", "selection_reasons")
_STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
_HYDROPHOBIC = frozenset("AILMFWVY")


def _read_seed(path: Path) -> dict[str, set[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"protein_id", "selection_reasons"}
        if not required <= set(reader.fieldnames or ()):
            raise InputValidationError("Seed audit lacks protein_id/selection_reasons")
        reasons: dict[str, set[str]] = {}
        for row in reader:
            protein_id = (row.get("protein_id") or "").strip()
            if not protein_id:
                raise InputValidationError("Seed audit contains an empty protein_id")
            reasons.setdefault(protein_id, set()).update(
                value
                for value in (row.get("selection_reasons") or "").split("|")
                if value
            )
    return reasons


def _hydrophobic_fraction(sequence: str) -> float:
    prefix = sequence[: min(25, len(sequence))]
    return sum(residue in _HYDROPHOBIC for residue in prefix) / len(prefix)


def _maximum_window_fraction(sequence: str, width: int = 19) -> float:
    if len(sequence) <= width:
        return sum(residue in _HYDROPHOBIC for residue in sequence) / len(sequence)
    return max(
        sum(residue in _HYDROPHOBIC for residue in sequence[index : index + width]) / width
        for index in range(len(sequence) - width + 1)
    )


def build_subset(
    *,
    proteome_fasta: Path,
    seed_audit: Path,
    domain_table: Path,
    sample_size: int = 8,
) -> tuple[list[ProteinRecord], dict[str, set[str]]]:
    proteins = LocalFastaLoader().load(proteome_fasta)
    by_id = {protein.protein_id: protein for protein in proteins}
    reasons = _read_seed(seed_audit)
    unknown_seed = sorted(set(reasons) - set(by_id))
    if unknown_seed:
        raise InputValidationError(
            "Seed subset IDs absent from proteome: " + ", ".join(unknown_seed)
        )

    domains = LocalDomainTsvLoader().load(domain_table)
    domain_counts = Counter(record.protein_id for record in domains)
    unknown_domains = sorted(set(domain_counts) - set(by_id))
    if unknown_domains:
        raise InputValidationError(
            "Domain IDs absent from proteome: " + ", ".join(unknown_domains)
        )

    def add(ids: Sequence[str], reason: str) -> None:
        for protein_id in ids:
            reasons.setdefault(protein_id, set()).add(reason)

    add(
        sorted(
            (protein_id for protein_id, count in domain_counts.items() if count > 1)
        )[:sample_size],
        "multiple_pfam_domains",
    )
    add(
        sorted(set(by_id) - set(domain_counts))[:sample_size],
        "no_pfam_hit",
    )
    add(
        sorted(
            protein.protein_id
            for protein in proteins
            if set(protein.sequence) - _STANDARD
        )[:sample_size],
        "non_standard_residue",
    )
    add(
        [
            protein.protein_id
            for protein in sorted(
                proteins,
                key=lambda item: (-_hydrophobic_fraction(item.sequence), item.protein_id),
            )[:sample_size]
        ],
        "n_terminal_hydrophobic",
    )
    add(
        [
            protein.protein_id
            for protein in sorted(
                proteins,
                key=lambda item: (_maximum_window_fraction(item.sequence), item.protein_id),
            )[:sample_size]
        ],
        "no_obvious_hydrophobic_segment",
    )
    selected = [by_id[protein_id] for protein_id in sorted(reasons)]
    return selected, reasons


def write_subset(
    proteins: Sequence[ProteinRecord],
    reasons: dict[str, set[str]],
    *,
    fasta_output: Path,
    audit_output: Path,
) -> None:
    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    with fasta_output.open("w", encoding="utf-8", newline="\n") as handle:
        for protein in proteins:
            protein_id = str(protein.protein_id)
            sequence = str(protein.sequence)
            handle.write(f">{protein_id}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")

    audit_output.parent.mkdir(parents=True, exist_ok=True)
    with audit_output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=AUDIT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for protein in proteins:
            protein_id = str(protein.protein_id)
            writer.writerow(
                {
                    "protein_id": protein_id,
                    "length_aa": len(str(protein.sequence)),
                    "description": str(protein.description),
                    "selection_reasons": "|".join(sorted(reasons[protein_id])),
                }
            )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteome-fasta", required=True, type=Path)
    parser.add_argument("--seed-audit", required=True, type=Path)
    parser.add_argument("--domain-table", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--sample-size", type=int, default=8)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.sample_size < 1:
        raise InputValidationError("sample-size must be at least 1")
    proteins, reasons = build_subset(
        proteome_fasta=args.proteome_fasta,
        seed_audit=args.seed_audit,
        domain_table=args.domain_table,
        sample_size=args.sample_size,
    )
    write_subset(
        proteins,
        reasons,
        fasta_output=args.fasta_output,
        audit_output=args.audit_output,
    )
    print(f"PSORTb subset proteins: {len(proteins)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
