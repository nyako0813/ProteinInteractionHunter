#!/usr/bin/env python3
"""Build the deterministic MA_4115 eggNOG readiness subset."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path

from protein_interaction_hunter.adapters.local.fasta import LocalFastaLoader
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models.protein import ProteinRecord

AUDIT_COLUMNS = ("protein_id", "length_aa", "description", "selection_reasons")
KEYWORD_POLICIES = {
    "enzyme_product_sample": re.compile(r"\b(?:enzyme|ase)\b", re.IGNORECASE),
    "rna_product_sample": re.compile(r"\b(?:rna|ribosom|trna|rrna)\b", re.IGNORECASE),
    "sulfur_product_sample": re.compile(
        r"\b(?:sulfur|sulphur|thiol|thio|sulfide|sulfate)\b",
        re.IGNORECASE,
    ),
    "atp_product_sample": re.compile(r"\b(?:atp|atpase)\b", re.IGNORECASE),
    "redox_product_sample": re.compile(
        r"\b(?:redox|oxid|reduct|ferredoxin|thioredoxin)\b",
        re.IGNORECASE,
    ),
    "transporter_product_sample": re.compile(
        r"\b(?:transport|permease|symporter|antiporter)\b",
        re.IGNORECASE,
    ),
}


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise InputValidationError(f"Subset source not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required <= set(reader.fieldnames or ()):
            raise InputValidationError(
                f"Subset source {path} lacks: "
                + ", ".join(sorted(required - set(reader.fieldnames or ())))
            )
        return [dict(row) for row in reader]


def _read_seed(path: Path) -> dict[str, set[str]]:
    rows = _read_tsv(path, {"protein_id", "selection_reasons"})
    result: dict[str, set[str]] = {}
    for row in rows:
        protein_id = row["protein_id"].strip()
        if not protein_id:
            raise InputValidationError("Subset seed contains an empty protein_id")
        result.setdefault(protein_id, set()).update(
            reason for reason in row["selection_reasons"].split("|") if reason
        )
    return result


def _add(
    reasons: dict[str, set[str]],
    protein_ids: Iterable[str],
    reason: str,
) -> None:
    for protein_id in protein_ids:
        reasons.setdefault(protein_id, set()).add(reason)


def build_subset(
    *,
    proteome_fasta: Path,
    seed_audit: Path,
    annotation_table: Path,
    profile_pair_audit: Path,
    orthology_table: Path,
    interpro_audit: Path,
    sample_size: int = 4,
    ranked_sample_size: int = 6,
) -> tuple[list[ProteinRecord], dict[str, set[str]]]:
    proteins = LocalFastaLoader().load(proteome_fasta)
    by_id = {str(protein.protein_id): protein for protein in proteins}
    reasons = _read_seed(seed_audit)
    unknown_seed = sorted(set(reasons) - set(by_id))
    if unknown_seed:
        raise InputValidationError(
            "Subset seed IDs absent from proteome: " + ", ".join(unknown_seed)
        )

    annotations = _read_tsv(
        annotation_table,
        {"protein_id", "product", "localization_annotation"},
    )
    annotation_by_id = {row["protein_id"]: row for row in annotations}
    unknown_annotation_ids = sorted(set(annotation_by_id) - set(by_id))
    if unknown_annotation_ids:
        raise InputValidationError(
            "Annotation IDs absent from proteome: " + ", ".join(unknown_annotation_ids)
        )
    for reason, pattern in KEYWORD_POLICIES.items():
        matches = sorted(
            protein_id
            for protein_id, row in annotation_by_id.items()
            if pattern.search(row["product"])
        )
        _add(reasons, matches[:sample_size], reason)

    localization_groups: dict[str, list[str]] = {}
    for protein_id, row in annotation_by_id.items():
        localization = row["localization_annotation"].strip()
        if localization:
            localization_groups.setdefault(localization, []).append(protein_id)
    for localization, protein_ids in sorted(localization_groups.items()):
        safe = re.sub(r"[^a-z0-9]+", "_", localization.casefold()).strip("_")
        _add(
            reasons,
            sorted(protein_ids)[:sample_size],
            f"psortb_localization_{safe}",
        )

    profile_rows = _read_tsv(
        profile_pair_audit,
        {"candidate_protein_id", "profile_similarity", "threshold_result"},
    )
    profile_sorted = sorted(
        profile_rows,
        key=lambda row: (
            -float(row["profile_similarity"] or -1.0),
            row["candidate_protein_id"],
        ),
    )
    threshold_ids = [
        row["candidate_protein_id"]
        for row in profile_sorted
        if row["threshold_result"].casefold() == "true"
    ][:ranked_sample_size]
    _add(reasons, threshold_ids, "profile_threshold_pass_top")
    _add(
        reasons,
        (row["candidate_protein_id"] for row in profile_sorted[:ranked_sample_size]),
        "profile_high_similarity_top",
    )

    orthology_rows = _read_tsv(
        orthology_table,
        {"protein_id", "relationship", "paralog_ambiguity"},
    )
    one_to_one_counts = Counter(
        row["protein_id"] for row in orthology_rows if row["relationship"] == "one_to_one"
    )
    _add(
        reasons,
        (
            protein_id
            for protein_id, _ in sorted(
                one_to_one_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )[:ranked_sample_size]
        ),
        "one_to_one_orthology_count_top",
    )
    ambiguous_ids = sorted(
        {
            row["protein_id"]
            for row in orthology_rows
            if row["paralog_ambiguity"].casefold() == "true"
        }
    )
    _add(
        reasons,
        ambiguous_ids[:ranked_sample_size],
        "paralog_ambiguity_sample",
    )
    assigned_ids = {row["protein_id"] for row in orthology_rows}
    _add(
        reasons,
        sorted(set(by_id) - assigned_ids)[:ranked_sample_size],
        "orthogroup_unassigned_sample",
    )

    interpro_rows = _read_tsv(interpro_audit, {"protein_id", "go_terms"})
    interpro_go_ids = {row["protein_id"] for row in interpro_rows if row["go_terms"].strip()}
    interpro_seen_ids = {row["protein_id"] for row in interpro_rows}
    _add(
        reasons,
        sorted(interpro_go_ids)[:sample_size],
        "interpro_go_present_sample",
    )
    _add(
        reasons,
        sorted(set(by_id) - interpro_go_ids)[:sample_size],
        "interpro_go_absent_sample",
    )
    _add(
        reasons,
        sorted(interpro_seen_ids)[:sample_size],
        "pfam_present_sample",
    )
    _add(
        reasons,
        sorted(set(by_id) - interpro_seen_ids)[:sample_size],
        "pfam_absent_sample",
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
            handle.write(f">{protein.protein_id}\n")
            sequence = str(protein.sequence)
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
            writer.writerow(
                {
                    "protein_id": str(protein.protein_id),
                    "length_aa": len(str(protein.sequence)),
                    "description": str(protein.description),
                    "selection_reasons": "|".join(sorted(reasons[str(protein.protein_id)])),
                }
            )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteome-fasta", required=True, type=Path)
    parser.add_argument("--seed-audit", required=True, type=Path)
    parser.add_argument("--annotation-table", required=True, type=Path)
    parser.add_argument("--profile-pair-audit", required=True, type=Path)
    parser.add_argument("--orthology-table", required=True, type=Path)
    parser.add_argument("--interpro-audit", required=True, type=Path)
    parser.add_argument("--sample-size", type=int, default=4)
    parser.add_argument("--ranked-sample-size", type=int, default=6)
    parser.add_argument("--fasta-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.sample_size < 1 or args.ranked_sample_size < 1:
        raise InputValidationError("Subset sample sizes must be positive")
    proteins, reasons = build_subset(
        proteome_fasta=args.proteome_fasta,
        seed_audit=args.seed_audit,
        annotation_table=args.annotation_table,
        profile_pair_audit=args.profile_pair_audit,
        orthology_table=args.orthology_table,
        interpro_audit=args.interpro_audit,
        sample_size=args.sample_size,
        ranked_sample_size=args.ranked_sample_size,
    )
    write_subset(
        proteins,
        reasons,
        fasta_output=args.fasta_output,
        audit_output=args.audit_output,
    )
    print(f"Functional subset proteins: {len(proteins)}")
    print(f"Functional subset amino acids: {sum(len(str(p.sequence)) for p in proteins)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
