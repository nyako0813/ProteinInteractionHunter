"""Build formal three-state phylogenetic profiles from orthology observations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

FORMAL_COLUMNS = (
    "protein_id",
    "species_id",
    "presence",
    "taxonomic_group",
    "source",
    "source_record_id",
)

AUDIT_COLUMNS = (
    "protein_id",
    "species_id",
    "state",
    "binary_presence",
    "ortholog_count",
    "relationships",
    "uncertain",
    "decision_reason",
    "rule_version",
)

PAIR_COLUMNS = (
    "query_protein_id",
    "candidate_protein_id",
    "query_informative_species",
    "candidate_informative_species",
    "jointly_informative_species",
    "shared_presence",
    "shared_absence",
    "discordant_presence",
    "query_only_presence",
    "candidate_only_presence",
    "missing_species",
    "similarity_numerator",
    "similarity_denominator",
    "profile_similarity",
    "minimum_shared_species",
    "minimum_informative_species",
    "minimum_profile_similarity",
    "threshold_result",
    "shared_absence_fraction_of_numerator",
    "calculation_rule_version",
)

COVERAGE_COLUMNS = ("metric", "count", "percent")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def _state(records: list[dict[str, str]]) -> tuple[str, str, str]:
    if not records:
        return "not_detected", "false", "valid_proteome_no_orthofinder_orthologue"
    usable = [row for row in records if row["relationship"] != "fragment_only"]
    if not usable:
        return "fragment_only", "", "fragment_only_support"
    relationships = {row["relationship"] for row in usable}
    if relationships & {"many_to_one", "many_to_many"}:
        return "present_ambiguous", "", "query_side_or_bilateral_paralogy"
    if "one_to_many" in relationships or len({row["ortholog_id"] for row in records}) > 1:
        return "present_multi_copy", "true", "multiple_reference_orthologues"
    return "present_unique", "true", "single_pairwise_orthologue"


def _availability(row: dict[str, str]) -> str:
    value = (row.get("proteome_status") or "valid").strip()
    allowed = {"valid", "species_missing", "proteome_invalid", "not_evaluated"}
    if value not in allowed:
        raise ValueError(f"Unknown proteome_status: {value}")
    return value


def build_profiles(
    *,
    orthology_path: Path,
    mapping_path: Path,
    panel_path: Path,
    formal_output: Path,
    audit_output: Path,
    pair_audit_output: Path,
    coverage_output: Path,
    metadata_output: Path,
    query_species_id: str,
    query_protein_id: str,
    source_version: str,
    minimum_shared_species: int,
    minimum_informative_species: int,
    minimum_profile_similarity: float,
) -> dict[str, Any]:
    with orthology_path.open(encoding="utf-8", newline="") as handle:
        orthology = list(csv.DictReader(handle, delimiter="\t"))
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        mappings = list(csv.DictReader(handle, delimiter="\t"))
    proteins = sorted(
        {row["raw_protein_id"] for row in mappings if row["species_id"] == query_species_id}
    )
    if query_protein_id not in proteins:
        raise ValueError(f"Query protein absent from mapping: {query_protein_id}")
    with panel_path.open(encoding="utf-8", newline="") as handle:
        panel = list(csv.DictReader(handle, delimiter="\t"))
    species_rows = sorted(
        (
            row
            for row in panel
            if row.get("selection_status") == "selected" and row["species_id"] != query_species_id
        ),
        key=lambda row: row["species_id"],
    )
    species_ids = [row["species_id"] for row in species_rows]
    group_by_species = {row["species_id"]: row["taxonomic_group"] for row in species_rows}
    availability = {row["species_id"]: _availability(row) for row in species_rows}
    observations: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    protein_set = set(proteins)
    species_set = set(species_ids)
    unknown_proteins: set[str] = set()
    unknown_species: set[str] = set()
    for row in orthology:
        if row["protein_id"] not in protein_set:
            unknown_proteins.add(row["protein_id"])
            continue
        if row["reference_id"] not in species_set:
            unknown_species.add(row["reference_id"])
            continue
        observations[(row["protein_id"], row["reference_id"])].append(row)
    if unknown_proteins or unknown_species:
        raise ValueError(
            f"Unknown orthology keys: proteins={len(unknown_proteins)}, "
            f"species={len(unknown_species)}"
        )
    rule_version = "ma4115-orthology-profile-v1"
    formal_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    profiles: dict[str, dict[str, bool | None]] = defaultdict(dict)
    state_counts: Counter[str] = Counter()
    for protein_id in proteins:
        for species_id in species_ids:
            records = observations.get((protein_id, species_id), [])
            if availability[species_id] == "valid":
                state, presence, reason = _state(records)
            else:
                state = availability[species_id]
                presence = ""
                reason = f"comparison_{state}"
            value = True if presence == "true" else False if presence == "false" else None
            profiles[protein_id][species_id] = value
            state_counts[state] += 1
            relationships = sorted({row["relationship"] for row in records})
            formal_rows.append(
                {
                    "protein_id": protein_id,
                    "species_id": species_id,
                    "presence": presence,
                    "taxonomic_group": group_by_species[species_id],
                    "source": f"OrthoFinder {source_version}",
                    "source_record_id": (
                        "|".join(sorted({row["source_record_id"] for row in records}))
                        if records
                        else f"{protein_id}:{species_id}:not_detected"
                    ),
                }
            )
            audit_rows.append(
                {
                    "protein_id": protein_id,
                    "species_id": species_id,
                    "state": state,
                    "binary_presence": presence,
                    "ortholog_count": len({row["ortholog_id"] for row in records}),
                    "relationships": "|".join(relationships),
                    "uncertain": str(value is None).lower(),
                    "decision_reason": reason,
                    "rule_version": rule_version,
                }
            )
    _write_tsv(formal_output, FORMAL_COLUMNS, formal_rows)
    _write_tsv(audit_output, AUDIT_COLUMNS, audit_rows)
    query_profile = profiles[query_protein_id]
    pair_rows: list[dict[str, Any]] = []
    comparable = 0
    for candidate in proteins:
        candidate_profile = profiles[candidate]
        shared_presence = shared_absence = query_only = candidate_only = missing = 0
        query_informative = sum(value is not None for value in query_profile.values())
        candidate_informative = sum(value is not None for value in candidate_profile.values())
        for species_id in species_ids:
            query_value = query_profile[species_id]
            candidate_value = candidate_profile[species_id]
            if query_value is None or candidate_value is None:
                missing += 1
            elif query_value and candidate_value:
                shared_presence += 1
            elif not query_value and not candidate_value:
                shared_absence += 1
            elif query_value:
                query_only += 1
            else:
                candidate_only += 1
        discordant = query_only + candidate_only
        jointly_informative = shared_presence + shared_absence + discordant
        numerator = shared_presence + shared_absence
        similarity = numerator / jointly_informative if jointly_informative else None
        threshold_result = (
            jointly_informative >= minimum_informative_species
            and shared_presence >= minimum_shared_species
            and similarity is not None
            and similarity >= minimum_profile_similarity
        )
        comparable += jointly_informative >= minimum_informative_species
        pair_rows.append(
            {
                "query_protein_id": query_protein_id,
                "candidate_protein_id": candidate,
                "query_informative_species": query_informative,
                "candidate_informative_species": candidate_informative,
                "jointly_informative_species": jointly_informative,
                "shared_presence": shared_presence,
                "shared_absence": shared_absence,
                "discordant_presence": discordant,
                "query_only_presence": query_only,
                "candidate_only_presence": candidate_only,
                "missing_species": missing,
                "similarity_numerator": numerator,
                "similarity_denominator": jointly_informative,
                "profile_similarity": f"{similarity:.12f}" if similarity is not None else "",
                "minimum_shared_species": minimum_shared_species,
                "minimum_informative_species": minimum_informative_species,
                "minimum_profile_similarity": f"{minimum_profile_similarity:.12f}",
                "threshold_result": str(threshold_result).lower(),
                "shared_absence_fraction_of_numerator": (
                    f"{shared_absence / numerator:.12f}" if numerator else ""
                ),
                "calculation_rule_version": "mvp1h-phylogenetic-profile-v1",
            }
        )
    _write_tsv(pair_audit_output, PAIR_COLUMNS, pair_rows)
    total_cells = len(proteins) * len(species_ids)
    query_values = list(query_profile.values())
    metrics = {
        "query_proteins_total": len(proteins),
        "profiles_generated": len(proteins),
        "profiles_with_minimum_informative_species": sum(
            sum(value is not None for value in profile.values()) >= minimum_informative_species
            for profile in profiles.values()
        ),
        "profiles_below_minimum_informative_species": sum(
            sum(value is not None for value in profile.values()) < minimum_informative_species
            for profile in profiles.values()
        ),
        "query_profile_informative_species": sum(value is not None for value in query_values),
        "query_profile_present_species": sum(value is True for value in query_values),
        "query_profile_absent_species": sum(value is False for value in query_values),
        "query_profile_uncertain_species": sum(value is None for value in query_values),
        "candidate_profiles_comparable_to_query": comparable,
        "candidate_profiles_not_comparable": len(proteins) - comparable,
        "species_total": len(species_ids),
        "species_valid": sum(value == "valid" for value in availability.values()),
        "species_excluded": sum(value != "valid" for value in availability.values()),
        "missing_profile_cells": state_counts["species_missing"]
        + state_counts["proteome_invalid"]
        + state_counts["not_evaluated"],
        "uncertain_profile_cells": state_counts["present_ambiguous"]
        + state_counts["fragment_only"],
        "total_profile_cells": total_cells,
    }
    coverage_rows = [
        {
            "metric": metric,
            "count": count,
            "percent": f"{100 * count / total_cells:.6f}" if metric.endswith("_cells") else "",
        }
        for metric, count in metrics.items()
    ]
    _write_tsv(coverage_output, COVERAGE_COLUMNS, coverage_rows)
    similarities = sorted(
        float(row["profile_similarity"]) for row in pair_rows if row["profile_similarity"]
    )
    metadata = {
        "schema_version": "1.0",
        "rule_version": rule_version,
        "binary_mapping": {
            "present_unique": True,
            "present_multi_copy": True,
            "present_ambiguous": None,
            "fragment_only": None,
            "not_detected": False,
            "species_missing": None,
            "proteome_invalid": None,
            "not_evaluated": None,
        },
        "missing_is_absence": False,
        "not_detected_requires_valid_proteome": True,
        "source_version": source_version,
        "source_orthology_sha256": _sha256(orthology_path),
        "formal_sha256": _sha256(formal_output),
        "audit_sha256": _sha256(audit_output),
        "pair_audit_sha256": _sha256(pair_audit_output),
        "coverage_sha256": _sha256(coverage_output),
        "state_counts": dict(sorted(state_counts.items())),
        "metrics": metrics,
        "similarity": {
            "minimum": min(similarities) if similarities else None,
            "maximum": max(similarities) if similarities else None,
            "threshold_met_count": sum(row["threshold_result"] == "true" for row in pair_rows),
        },
    }
    metadata_output.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthology", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--formal-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--pair-audit-output", required=True, type=Path)
    parser.add_argument("--coverage-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--query-species-id", required=True)
    parser.add_argument("--query-protein-id", required=True)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--minimum-shared-species", type=int, default=2)
    parser.add_argument("--minimum-informative-species", type=int, default=3)
    parser.add_argument("--minimum-profile-similarity", type=float, default=0.8)
    args = parser.parse_args()
    metadata = build_profiles(
        orthology_path=args.orthology,
        mapping_path=args.mapping,
        panel_path=args.panel,
        formal_output=args.formal_output,
        audit_output=args.audit_output,
        pair_audit_output=args.pair_audit_output,
        coverage_output=args.coverage_output,
        metadata_output=args.metadata_output,
        query_species_id=args.query_species_id,
        query_protein_id=args.query_protein_id,
        source_version=args.source_version,
        minimum_shared_species=args.minimum_shared_species,
        minimum_informative_species=args.minimum_informative_species,
        minimum_profile_similarity=args.minimum_profile_similarity,
    )
    print(json.dumps(metadata, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
