from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.prepare_orthology_proteomes import (  # noqa: E402
    AUDIT_COLUMNS,
    MAPPING_COLUMNS,
    prepare_proteomes,
)


def _write_panel(path: Path) -> Path:
    path.write_text(
        "species_id\tassembly_accession\tselection_status\n"
        "GCF_000000001_1\tGCF_000000001.1\tselected\n"
        "GCF_000000002_1\tGCF_000000002.1\tselected\n",
        encoding="utf-8",
        newline="\n",
    )
    return path


def _setup(
    tmp_path: Path, query_fasta: str, comparison_fasta: str
) -> tuple[Path, Path, str, Path, Path, Path, Path]:
    data_root = tmp_path / "data"
    for accession, fasta in (
        ("GCF_000000001.1", query_fasta),
        ("GCF_000000002.1", comparison_fasta),
    ):
        directory = data_root / accession
        directory.mkdir(parents=True)
        (directory / "protein.faa").write_text(fasta, encoding="utf-8", newline="\n")
    panel = _write_panel(tmp_path / "panel.tsv")
    expected = hashlib.sha256(
        (data_root / "GCF_000000001.1" / "protein.faa").read_bytes()
    ).hexdigest()
    normalized = tmp_path / "normalized"
    mapping = tmp_path / "mapping.tsv"
    audit = tmp_path / "audit.tsv"
    manifest = tmp_path / "manifest.json"
    return data_root, panel, expected, normalized, mapping, audit, manifest


def _run(
    tmp_path: Path, query_fasta: str, comparison_fasta: str
) -> tuple[dict[str, Any], Path, Path, Path, Path]:
    args = _setup(tmp_path, query_fasta, comparison_fasta)
    data_root, panel, expected, normalized, mapping, audit, manifest = args
    result = prepare_proteomes(
        panel_path=panel,
        data_root=data_root,
        normalized_directory=normalized,
        mapping_output=mapping,
        audit_output=audit,
        manifest_output=manifest,
        query_assembly="GCF_000000001.1",
        query_protein_id="QUERY",
        expected_query_fasta_sha256=expected,
    )
    return result, normalized, mapping, audit, manifest


def test_valid_fasta_normalizes_and_round_trips(tmp_path: Path) -> None:
    result, normalized, mapping, audit, manifest = _run(
        tmp_path,
        ">QUERY query header\nACDEFGHIKLMNPQRSTVWY\n>A other header\nAAAA\n",
        ">B comparison header\nCCCC\n>C duplicate\nCCCC\n",
    )
    assert result["species_count"] == 2
    assert result["total_protein_count"] == 4
    with mapping.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(rows[0]) == MAPPING_COLUMNS
    query = next(row for row in rows if row["raw_protein_id"] == "QUERY")
    assert query["normalized_protein_id"] == "GCF_000000001_1__QUERY"
    assert query["original_header"] == "QUERY query header"
    assert {row["raw_protein_id"] for row in rows} == {"QUERY", "A", "B", "C"}
    with audit.open(encoding="utf-8", newline="") as handle:
        audits = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(audits[0]) == AUDIT_COLUMNS
    duplicate_audit = next(row for row in audits if row["assembly_accession"].endswith("2.1"))
    assert duplicate_audit["duplicate_sequence_group_count"] == "1"
    assert mapping.read_bytes().count(b"\r") == 0
    assert all(path.read_bytes().count(b"\r") == 0 for path in normalized.glob("*.faa"))
    assert len(json.loads(manifest.read_text())["mapping_sha256"]) == 64


@pytest.mark.parametrize(
    ("fasta", "message"),
    [
        (">QUERY\nAAAA\n>QUERY\nCCCC\n", "Duplicate FASTA IDs"),
        (">QUERY\n\n>A\nCCCC\n", "Empty FASTA sequences"),
        (">QUERY\nAAA?AAA\n", "Invalid FASTA residues"),
        (">QUERY\nAAA*AAA\n", "Internal stop residues"),
        (">QUERY\nAAAA*\n", "Terminal stop residues"),
    ],
)
def test_rejects_invalid_query_fasta(tmp_path: Path, fasta: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _run(tmp_path, fasta, ">B\nCCCC\n")


def test_checksum_mismatch_rejected(tmp_path: Path) -> None:
    data_root, panel, _, normalized, mapping, audit, manifest = _setup(
        tmp_path, ">QUERY\nAAAA\n", ">B\nCCCC\n"
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        prepare_proteomes(
            panel_path=panel,
            data_root=data_root,
            normalized_directory=normalized,
            mapping_output=mapping,
            audit_output=audit,
            manifest_output=manifest,
            query_assembly="GCF_000000001.1",
            query_protein_id="QUERY",
            expected_query_fasta_sha256="0" * 64,
        )


def test_query_must_round_trip_once(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="observed 0"):
        _run(tmp_path, ">NOT_QUERY\nAAAA\n", ">B\nCCCC\n")


def test_nonstandard_residue_is_audited_not_silently_removed(tmp_path: Path) -> None:
    _, normalized, mapping, audit, _ = _run(tmp_path, ">QUERY\nAXAA\n", ">B\nCCCC\n")
    assert "\tX\t" in mapping.read_text()
    assert "nonstandard_residues" in audit.read_text()
    assert "AXAA" in next(normalized.glob("GCF_000000001_1.faa")).read_text()


def test_output_is_deterministic_for_input_order(tmp_path: Path) -> None:
    first = _run(tmp_path / "a", ">QUERY\nAAAA\n>A\nCCCC\n", ">B\nDDDD\n")
    second = _run(tmp_path / "b", ">A\nCCCC\n>QUERY\nAAAA\n", ">B\nDDDD\n")
    assert (first[1] / "GCF_000000001_1.faa").read_bytes() == (
        second[1] / "GCF_000000001_1.faa"
    ).read_bytes()
    with first[2].open(encoding="utf-8", newline="") as handle:
        first_rows = list(csv.DictReader(handle, delimiter="	"))
    with second[2].open(encoding="utf-8", newline="") as handle:
        second_rows = list(csv.DictReader(handle, delimiter="	"))
    for rows in (first_rows, second_rows):
        for row in rows:
            row.pop("input_fasta_sha256")
    assert first_rows == second_rows
