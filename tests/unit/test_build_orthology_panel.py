from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.build_orthology_panel import (  # noqa: E402
    AUDIT_COLUMNS,
    FORMAL_COLUMNS,
    build_panel,
)


def _record(
    accession: str,
    organism: str,
    *,
    level: str = "Complete Genome",
    category: str = "",
    status: str = "current",
    annotation: bool = True,
    proteins: int = 2000,
    completeness: float | None = 98.0,
    contamination: float | None = 1.0,
    atypical: bool = False,
    models: list[str] | None = None,
) -> dict[str, Any]:
    record = {
        "accession": accession,
        "organism": {
            "organism_name": organism,
            "tax_id": 1,
            "infraspecific_names": {"strain": "strain-a"},
        },
        "assembly_info": {
            "assembly_level": level,
            "assembly_name": "ASM",
            "assembly_status": status,
            "refseq_category": category,
            "release_date": "2020-01-01",
            "bioproject_accession": "PRJNA1",
            "biosample": {
                "accession": "SAMN1",
                "models": models or ["Microbe"],
                "strain": "strain-a",
            },
            "atypical": {"is_atypical": atypical},
        },
        "assembly_stats": {
            "number_of_contigs": 1,
            "number_of_scaffolds": 1,
            "total_sequence_length": 2_000_000,
            "gc_percent": 50.0,
        },
        "annotation_info": {
            "name": "PGAP" if annotation else "",
            "provider": "NCBI RefSeq",
            "stats": {"gene_counts": {"protein_coding": proteins}},
        },
        "average_nucleotide_identity": {"taxonomy_check_status": "OK"},
    }
    if completeness is not None:
        record["checkm_info"] = {
            "completeness": completeness,
            "contamination": contamination,
            "checkm_species_tax_id": 2,
        }
    return record


def _run(
    tmp_path: Path,
    records: list[dict[str, Any]],
    *,
    query: str = "GCF_000000001.1",
) -> tuple[
    dict[str, Any],
    list[dict[str, str]],
    list[dict[str, str]],
    Path,
    Path,
    Path,
    Path,
]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    metadata = tmp_path / "metadata.jsonl"
    metadata.write_text(
        "".join(json.dumps(record) + "\n" for record in records),
        encoding="utf-8",
        newline="\n",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "target_taxon\tpanel_layer\ttaxonomic_group\ttarget_reason\n"
        "Species alpha\tlayer_1\tgroup_a\tquery\n"
        "Species beta\tlayer_2\tgroup_b\tcomparison\n",
        encoding="utf-8",
        newline="\n",
    )
    policy = tmp_path / "policy.yaml"
    policy.write_text(
        yaml.safe_dump(
            {
                "policy_version": "test-panel-v1",
                "query_assembly": query,
                "selection": {
                    "allowed_assembly_levels": ["Complete Genome", "Chromosome"],
                    "require_annotation": True,
                    "minimum_protein_count": 500,
                    "maximum_protein_count": 12000,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
        newline="\n",
    )
    formal = tmp_path / "panel.tsv"
    audit = tmp_path / "audit.tsv"
    accessions = tmp_path / "accessions.txt"
    manifest = tmp_path / "manifest.json"
    result = build_panel(
        metadata_path=metadata,
        targets_path=targets,
        policy_path=policy,
        formal_output=formal,
        audit_output=audit,
        accessions_output=accessions,
        manifest_output=manifest,
        source_version="18.33.1",
        dataformat_version="18.33.1",
        source_command="datasets summary genome taxon",
        query_date_utc="2026-07-30T00:00:00Z",
    )
    with formal.open(encoding="utf-8", newline="") as handle:
        formal_rows = list(csv.DictReader(handle, delimiter="\t"))
    with audit.open(encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    return result, formal_rows, audit_rows, formal, audit, accessions, manifest


def test_selects_query_and_best_representative_deterministically(tmp_path: Path) -> None:
    records = [
        _record("GCF_000000003.1", "Species beta strain 2", completeness=99.0),
        _record(
            "GCF_000000002.1",
            "Species beta strain 1",
            category="reference genome",
            completeness=95.0,
        ),
        _record("GCF_000000001.1", "Species alpha strain q", completeness=None),
    ]
    result, formal_rows, audit_rows, formal, audit, accessions, manifest = _run(tmp_path, records)
    assert result["selected_assembly_count"] == 2
    assert {row["assembly_accession"] for row in formal_rows} == {
        "GCF_000000001.1",
        "GCF_000000002.1",
    }
    assert any(row["selection_reason"] == "strain_redundancy" for row in audit_rows)
    assert tuple(formal_rows[0]) == FORMAL_COLUMNS
    assert tuple(audit_rows[0]) == AUDIT_COLUMNS
    assert accessions.read_text().endswith("\n")
    assert formal.read_bytes().count(b"\r") == 0
    assert audit.read_bytes().count(b"\r") == 0
    assert json.loads(manifest.read_text())["policy_version"] == "test-panel-v1"


@pytest.mark.parametrize(
    ("record", "reason"),
    [
        (
            _record("GCF_000000010.1", "Species beta x", level="Scaffold"),
            "assembly_level_not_allowed",
        ),
        (
            _record("GCF_000000010.1", "Species beta x", status="suppressed"),
            "obsolete_or_suppressed",
        ),
        (_record("GCF_000000010.1", "Species beta x", atypical=True), "atypical_assembly"),
        (
            _record(
                "GCF_000000010.1",
                "Species beta x",
                models=["Metagenome-assembled Genome"],
            ),
            "metagenome_assembled_genome",
        ),
        (_record("GCF_000000010.1", "Species beta x", annotation=False), "annotation_missing"),
        (_record("GCF_000000010.1", "Species beta x", proteins=100), "protein_count_below_minimum"),
    ],
)
def test_explicit_policy_exclusions(tmp_path: Path, record: dict[str, Any], reason: str) -> None:
    query = _record("GCF_000000001.1", "Species alpha query")
    _, formal_rows, audit_rows, *_ = _run(tmp_path, [query, record])
    assert len(formal_rows) == 1
    rejected = next(row for row in audit_rows if row["assembly_accession"] == record["accession"])
    assert reason in rejected["exclusion_reasons"]
    assert rejected["selection_status"] == "excluded"


def test_missing_metadata_is_retained_and_flagged_for_query(tmp_path: Path) -> None:
    query = _record("GCF_000000001.1", "Species alpha query", annotation=False)
    _, formal_rows, audit_rows, *_ = _run(tmp_path, [query])
    assert len(formal_rows) == 1
    query_audit = next(row for row in audit_rows if row["assembly_accession"])
    assert query_audit["manual_review"] == "true"
    assert "annotation_name" in query_audit["missing_metadata"]


def test_versionless_accession_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Versioned assembly accession"):
        _run(tmp_path, [_record("GCF_000000001", "Species alpha query")])


def test_duplicate_assembly_rejected(tmp_path: Path) -> None:
    record = _record("GCF_000000001.1", "Species alpha query")
    with pytest.raises(ValueError, match="Duplicate assembly accession"):
        _run(tmp_path, [record, record])


def test_order_and_checksums_are_deterministic(tmp_path: Path) -> None:
    records = [
        _record("GCF_000000001.1", "Species alpha query"),
        _record("GCF_000000002.1", "Species beta comparison"),
    ]
    first = _run(tmp_path / "a", records)
    second = _run(tmp_path / "b", list(reversed(records)))
    assert first[1] == second[1]
    assert first[5].read_bytes() == second[5].read_bytes()
    first_manifest = json.loads(first[6].read_text())
    assert len(first_manifest["inputs"]["metadata"]["sha256"]) == 64
