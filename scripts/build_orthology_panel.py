"""Build a deterministic, auditable comparative-proteome panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

ACCESSION_RE = re.compile(r"^(GC[AF]_\d+)\.(\d+)$")

FORMAL_COLUMNS = (
    "species_id",
    "assembly_accession",
    "organism_name",
    "tax_id",
    "species_tax_id",
    "taxonomic_group",
    "panel_layer",
    "assembly_level",
    "refseq_category",
    "selection_status",
    "selection_reason",
)

AUDIT_COLUMNS = (
    *FORMAL_COLUMNS,
    "target_taxon",
    "assembly_name",
    "release_date",
    "annotation_name",
    "annotation_provider",
    "protein_coding_gene_count",
    "contig_count",
    "scaffold_count",
    "total_sequence_length",
    "gc_percent",
    "checkm_completeness",
    "checkm_contamination",
    "strain",
    "isolate",
    "biosample_accession",
    "bioproject_accession",
    "assembly_status",
    "atypical",
    "mag",
    "duplicate_species",
    "strain_redundancy",
    "missing_metadata",
    "taxonomy_conflict",
    "manual_review",
    "exclusion_reasons",
    "policy_version",
    "source_name",
    "source_version",
    "dataformat_version",
    "source_command",
    "query_date_utc",
    "source_path",
    "source_sha256",
    "raw_metadata",
)


@dataclass(frozen=True)
class Target:
    target_taxon: str
    panel_layer: str
    taxonomic_group: str
    target_reason: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _nested(record: dict[str, Any], *keys: str) -> Any:
    value: Any = record
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _species_id(accession: str) -> str:
    match = ACCESSION_RE.fullmatch(accession)
    if match is None:
        raise ValueError(f"Versioned assembly accession required: {accession!r}")
    return f"{match.group(1)}_{match.group(2)}"


def _load_targets(path: Path) -> list[Target]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_taxon", "panel_layer", "taxonomic_group", "target_reason"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Target table missing columns: " + ", ".join(sorted(missing)))
        targets = [
            Target(**{key: (row.get(key) or "").strip() for key in required}) for row in reader
        ]
    if any(not value for target in targets for value in target.__dict__.values()):
        raise ValueError("Target table contains an empty required value")
    names = [target.target_taxon.casefold() for target in targets]
    if len(names) != len(set(names)):
        raise ValueError("Target table contains duplicate taxa")
    return targets


def _match_target(organism_name: str, targets: list[Target]) -> Target | None:
    normalized = organism_name.casefold()
    matches = [
        target
        for target in targets
        if normalized == target.target_taxon.casefold()
        or normalized.startswith(target.target_taxon.casefold() + " ")
    ]
    if not matches:
        return None
    return max(matches, key=lambda target: len(target.target_taxon))


def _candidate_row(
    record: dict[str, Any],
    target: Target,
    *,
    policy_version: str,
    source_name: str,
    source_version: str,
    dataformat_version: str,
    source_command: str,
    query_date_utc: str,
    source_path: Path,
    source_sha256: str,
) -> dict[str, str]:
    accession = _text(record.get("accession"))
    species_id = _species_id(accession)
    assembly_info = record.get("assembly_info") or {}
    annotation_info = record.get("annotation_info") or {}
    stats = record.get("assembly_stats") or {}
    checkm = record.get("checkm_info") or {}
    organism = record.get("organism") or {}
    biosample = assembly_info.get("biosample") or {}
    infraspecific = organism.get("infraspecific_names") or {}
    gene_counts = _nested(annotation_info, "stats", "gene_counts") or {}
    atypical_value = _nested(assembly_info, "atypical", "is_atypical")
    mag_value = _nested(assembly_info, "biosample", "models")
    missing_fields = [
        name
        for name, value in (
            ("organism_name", organism.get("organism_name")),
            ("tax_id", organism.get("tax_id")),
            ("assembly_level", assembly_info.get("assembly_level")),
            ("assembly_status", assembly_info.get("assembly_status")),
            ("annotation_name", annotation_info.get("name")),
            ("protein_coding_gene_count", gene_counts.get("protein_coding")),
        )
        if value in (None, "")
    ]
    taxonomy_status = _nested(record, "average_nucleotide_identity", "taxonomy_check_status")
    taxonomy_conflict = bool(taxonomy_status and taxonomy_status != "OK")
    return {
        "species_id": species_id,
        "assembly_accession": accession,
        "organism_name": _text(organism.get("organism_name")),
        "tax_id": _text(organism.get("tax_id")),
        "species_tax_id": _text(checkm.get("checkm_species_tax_id")),
        "taxonomic_group": target.taxonomic_group,
        "panel_layer": target.panel_layer,
        "assembly_level": _text(assembly_info.get("assembly_level")),
        "refseq_category": _text(assembly_info.get("refseq_category")),
        "selection_status": "candidate",
        "selection_reason": "",
        "target_taxon": target.target_taxon,
        "assembly_name": _text(assembly_info.get("assembly_name")),
        "release_date": _text(assembly_info.get("release_date")),
        "annotation_name": _text(annotation_info.get("name")),
        "annotation_provider": _text(annotation_info.get("provider")),
        "protein_coding_gene_count": _text(gene_counts.get("protein_coding")),
        "contig_count": _text(stats.get("number_of_contigs")),
        "scaffold_count": _text(stats.get("number_of_scaffolds")),
        "total_sequence_length": _text(stats.get("total_sequence_length")),
        "gc_percent": _text(stats.get("gc_percent")),
        "checkm_completeness": _text(checkm.get("completeness")),
        "checkm_contamination": _text(checkm.get("contamination")),
        "strain": _text(infraspecific.get("strain") or biosample.get("strain")),
        "isolate": _text(infraspecific.get("isolate") or biosample.get("isolate")),
        "biosample_accession": _text(biosample.get("accession")),
        "bioproject_accession": _text(assembly_info.get("bioproject_accession")),
        "assembly_status": _text(assembly_info.get("assembly_status")),
        "atypical": _text(atypical_value),
        "mag": "true"
        if isinstance(mag_value, list) and any("metagenome" in _text(x).lower() for x in mag_value)
        else "false",
        "duplicate_species": "false",
        "strain_redundancy": "false",
        "missing_metadata": "|".join(missing_fields),
        "taxonomy_conflict": str(taxonomy_conflict).lower(),
        "manual_review": str(bool(missing_fields or taxonomy_conflict)).lower(),
        "exclusion_reasons": "",
        "policy_version": policy_version,
        "source_name": source_name,
        "source_version": source_version,
        "dataformat_version": dataformat_version,
        "source_command": source_command,
        "query_date_utc": query_date_utc,
        "source_path": str(source_path.resolve()),
        "source_sha256": source_sha256,
        "raw_metadata": json.dumps(
            record, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ),
    }


def _rank(row: dict[str, str], query_accession: str) -> tuple[Any, ...]:
    category = {"reference genome": 0, "representative genome": 1, "": 2}
    level = {"Complete Genome": 0, "Chromosome": 1}
    completeness = float(row["checkm_completeness"] or "-1")
    contamination = float(row["checkm_contamination"] or "999")
    return (
        0 if row["assembly_accession"] == query_accession else 1,
        category.get(row["refseq_category"], 3),
        level.get(row["assembly_level"], 9),
        0 if row["annotation_name"] else 1,
        -completeness,
        contamination,
        row["assembly_accession"],
    )


def _policy_exclusions(row: dict[str, str], policy: dict[str, Any]) -> list[str]:
    selection = policy.get("selection") or {}
    reasons: list[str] = []
    if not ACCESSION_RE.fullmatch(row["assembly_accession"]):
        reasons.append("versionless_or_invalid_accession")
    if not row["assembly_accession"].startswith("GCF_"):
        reasons.append("not_refseq")
    allowed_levels = set(selection.get("allowed_assembly_levels") or ())
    if row["assembly_level"] not in allowed_levels:
        reasons.append("assembly_level_not_allowed")
    if row["assembly_status"].casefold() not in {"current", ""}:
        reasons.append("obsolete_or_suppressed")
    if row["atypical"].casefold() in {"true", "1", "yes"}:
        reasons.append("atypical_assembly")
    if row["mag"] == "true":
        reasons.append("metagenome_assembled_genome")
    if selection.get("require_annotation") and not row["annotation_name"]:
        reasons.append("annotation_missing")
    count_text = row["protein_coding_gene_count"]
    if count_text:
        count = int(count_text)
        if count < int(selection.get("minimum_protein_count", 0)):
            reasons.append("protein_count_below_minimum")
        if count > int(selection.get("maximum_protein_count", 10**12)):
            reasons.append("protein_count_above_maximum")
    else:
        reasons.append("protein_count_missing")
    return reasons


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def build_panel(
    *,
    metadata_path: Path,
    targets_path: Path,
    policy_path: Path,
    formal_output: Path,
    audit_output: Path,
    accessions_output: Path,
    manifest_output: Path,
    source_version: str,
    dataformat_version: str,
    source_command: str,
    query_date_utc: str,
) -> dict[str, Any]:
    policy = yaml.safe_load(policy_path.read_text(encoding="utf-8"))
    if not isinstance(policy, dict):
        raise ValueError("Panel policy root must be a mapping")
    policy_version = _text(policy.get("policy_version"))
    query_accession = _text(policy.get("query_assembly"))
    _species_id(query_accession)
    targets = _load_targets(targets_path)
    source_sha256 = _sha256(metadata_path)
    rows: list[dict[str, str]] = []
    seen_accessions: set[str] = set()
    for line_number, line in enumerate(metadata_path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        accession = _text(record.get("accession"))
        if accession in seen_accessions:
            raise ValueError(f"Duplicate assembly accession on line {line_number}: {accession}")
        seen_accessions.add(accession)
        organism_name = _text(_nested(record, "organism", "organism_name"))
        target = _match_target(organism_name, targets)
        if target is None:
            continue
        rows.append(
            _candidate_row(
                record,
                target,
                policy_version=policy_version,
                source_name="NCBI Datasets",
                source_version=source_version,
                dataformat_version=dataformat_version,
                source_command=source_command,
                query_date_utc=query_date_utc,
                source_path=metadata_path,
                source_sha256=source_sha256,
            )
        )
    by_target: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        reasons = _policy_exclusions(row, policy)
        if reasons and row["assembly_accession"] != query_accession:
            row["selection_status"] = "excluded"
            row["selection_reason"] = "policy_exclusion"
            row["exclusion_reasons"] = "|".join(sorted(reasons))
            row["manual_review"] = "true"
            continue
        if reasons:
            row["manual_review"] = "true"
            row["exclusion_reasons"] = "|".join(sorted(reasons))
        by_target.setdefault(row["target_taxon"], []).append(row)
    selected: list[dict[str, str]] = []
    for target in targets:
        candidates = sorted(
            by_target.get(target.target_taxon, []), key=lambda row: _rank(row, query_accession)
        )
        if not candidates:
            rows.append(
                {
                    **{column: "" for column in AUDIT_COLUMNS},
                    "target_taxon": target.target_taxon,
                    "taxonomic_group": target.taxonomic_group,
                    "panel_layer": target.panel_layer,
                    "selection_status": "excluded",
                    "selection_reason": "no_policy_eligible_assembly_returned",
                    "exclusion_reasons": "no_policy_eligible_assembly_returned",
                    "manual_review": "true",
                    "policy_version": policy_version,
                    "source_name": "NCBI Datasets",
                    "source_version": source_version,
                    "dataformat_version": dataformat_version,
                    "source_command": source_command,
                    "query_date_utc": query_date_utc,
                    "source_path": str(metadata_path.resolve()),
                    "source_sha256": source_sha256,
                }
            )
            continue
        chosen = candidates[0]
        chosen["selection_status"] = "selected"
        chosen["selection_reason"] = (
            "fixed_query_assembly"
            if chosen["assembly_accession"] == query_accession
            else "highest_policy_ranked_representative"
        )
        selected.append(chosen)
        for rejected in candidates[1:]:
            rejected["selection_status"] = "excluded"
            rejected["selection_reason"] = "strain_redundancy"
            rejected["strain_redundancy"] = "true"
            rejected["duplicate_species"] = "true"
            rejected["exclusion_reasons"] = "lower_ranked_same_target_taxon"
    if query_accession not in {row["assembly_accession"] for row in selected}:
        raise ValueError(f"Fixed query assembly was not selected: {query_accession}")

    def ordering(row: dict[str, str]) -> tuple[str, str, str, str]:
        return (
            row["panel_layer"],
            row["taxonomic_group"],
            row["organism_name"],
            row["assembly_accession"],
        )

    selected.sort(key=ordering)
    rows.sort(key=lambda row: (row["panel_layer"], row["target_taxon"], row["assembly_accession"]))
    _write_tsv(formal_output, FORMAL_COLUMNS, selected)
    _write_tsv(audit_output, AUDIT_COLUMNS, rows)
    accessions_output.parent.mkdir(parents=True, exist_ok=True)
    accessions_output.write_text(
        "".join(f"{row['assembly_accession']}\n" for row in selected),
        encoding="utf-8",
        newline="\n",
    )
    layer_counts: dict[str, int] = {}
    for row in selected:
        layer_counts[row["panel_layer"]] = layer_counts.get(row["panel_layer"], 0) + 1
    manifest = {
        "schema_version": "1.0",
        "policy_version": policy_version,
        "query_assembly": query_accession,
        "query_date_utc": query_date_utc,
        "source_name": "NCBI Datasets",
        "source_version": source_version,
        "dataformat_version": dataformat_version,
        "source_command": source_command,
        "inputs": {
            "metadata": {"path": str(metadata_path.resolve()), "sha256": source_sha256},
            "targets": {"path": str(targets_path.resolve()), "sha256": _sha256(targets_path)},
            "policy": {"path": str(policy_path.resolve()), "sha256": _sha256(policy_path)},
        },
        "selected_assembly_count": len(selected),
        "target_taxon_count": len(targets),
        "unresolved_target_count": sum(
            row["selection_reason"] == "no_policy_eligible_assembly_returned" for row in rows
        ),
        "layer_counts": dict(sorted(layer_counts.items())),
    }
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--policy", required=True, type=Path)
    parser.add_argument("--formal-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--accessions-output", required=True, type=Path)
    parser.add_argument("--manifest-output", required=True, type=Path)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--dataformat-version", required=True)
    parser.add_argument("--source-command", required=True)
    parser.add_argument(
        "--query-date-utc",
        default=None,
    )
    args = parser.parse_args()
    query_date = args.query_date_utc or datetime.now(UTC).isoformat().replace("+00:00", "Z")
    manifest = build_panel(
        metadata_path=args.metadata,
        targets_path=args.targets,
        policy_path=args.policy,
        formal_output=args.formal_output,
        audit_output=args.audit_output,
        accessions_output=args.accessions_output,
        manifest_output=args.manifest_output,
        source_version=args.source_version,
        dataformat_version=args.dataformat_version,
        source_command=args.source_command,
        query_date_utc=query_date,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
