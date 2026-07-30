from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.exceptions import InputValidationError  # noqa: E402
from scripts.convert_eggnog_annotations import (  # noqa: E402
    RAW_COLUMNS,
    convert_eggnog,
    parse_eggnog_annotations,
    write_conversion,
)


def _fasta(path: Path, ids: tuple[str, ...] = ("Q", "C")) -> Path:
    path.write_text("".join(f">{protein_id}\nMAAA\n" for protein_id in ids))
    return path


def _row(protein_id: str, **values: str) -> str:
    row = dict.fromkeys(RAW_COLUMNS, "-")
    row.update(
        {
            "query": protein_id,
            "seed_ortholog": "1.seed",
            "evalue": "1e-20",
            "score": "100",
            "eggNOG_OGs": "COG1@1|root",
        }
    )
    row.update(values)
    return "\t".join(row[column] for column in RAW_COLUMNS)


def _raw(path: Path, rows: tuple[str, ...], header: tuple[str, ...] = RAW_COLUMNS) -> Path:
    path.write_text(
        "## emapper-v2.1.15\n"
        "## /tool/emapper.py -i input.faa --data_dir db\n"
        + "#"
        + "\t".join(header)
        + "\n"
        + "\n".join(rows)
        + ("\n" if rows else "")
    )
    return path


def test_valid_metadata_empty_and_comma_fields(tmp_path: Path) -> None:
    raw = _raw(
        tmp_path / "raw.tsv",
        (
            _row(
                "Q",
                Description="specific function",
                GOs="GO:0001,GO:0002",
                KEGG_ko="ko:K1,ko:K2",
            ),
        ),
    )
    rows, metadata, duplicates = parse_eggnog_annotations(raw)
    assert rows[0].protein_id == "Q"
    assert rows[0].values[9] == "GO:0001,GO:0002"
    assert metadata[0] == "## emapper-v2.1.15"
    assert duplicates == 0


@pytest.mark.parametrize(
    "header",
    (
        RAW_COLUMNS[:-1],
        (*RAW_COLUMNS, "extra"),
        ("wrong", *RAW_COLUMNS[1:]),
    ),
)
def test_exact_header_required(tmp_path: Path, header: tuple[str, ...]) -> None:
    raw = _raw(tmp_path / "raw.tsv", (), header)
    with pytest.raises(InputValidationError, match="header mismatch"):
        parse_eggnog_annotations(raw)


def test_empty_file_rejected_but_header_only_accepted(tmp_path: Path) -> None:
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    with pytest.raises(InputValidationError, match="empty or lacks"):
        parse_eggnog_annotations(empty)
    header_only = _raw(tmp_path / "header.tsv", ())
    rows, _, _ = parse_eggnog_annotations(header_only)
    assert rows == []


def test_identical_duplicate_collapses_and_conflict_fails(tmp_path: Path) -> None:
    same = _row("Q")
    raw = _raw(tmp_path / "same.tsv", (same, same))
    rows, _, duplicates = parse_eggnog_annotations(raw)
    assert len(rows) == 1
    assert duplicates == 1
    conflict = _raw(
        tmp_path / "conflict.tsv",
        (_row("Q"), _row("Q", score="101")),
    )
    with pytest.raises(InputValidationError, match="Conflicting duplicate"):
        parse_eggnog_annotations(conflict)


def test_malformed_and_empty_query_fail(tmp_path: Path) -> None:
    malformed = _raw(tmp_path / "bad.tsv", ("\t".join(["Q"] * 20),))
    with pytest.raises(InputValidationError, match="20 columns"):
        parse_eggnog_annotations(malformed)
    empty_query = _raw(tmp_path / "empty-query.tsv", (_row(""),))
    with pytest.raises(InputValidationError, match="empty query ID"):
        parse_eggnog_annotations(empty_query)


def test_unknown_id_fails_and_query_is_preserved(tmp_path: Path) -> None:
    raw = _raw(tmp_path / "raw.tsv", (_row("UNKNOWN"),))
    with pytest.raises(InputValidationError, match="proteome-unknown"):
        convert_eggnog(
            input_path=raw,
            proteome_fasta=_fasta(tmp_path / "proteome.faa"),
            query_id="Q",
            source_version="2.1.15",
            database_version="5.0.2",
            command="auto",
        )


def test_hits_seed_no_hit_and_hit_no_annotation(tmp_path: Path) -> None:
    raw = _raw(tmp_path / "raw.tsv", (_row("Q"),))
    hits = tmp_path / "hits.tsv"
    hits.write_text("Q\tseed\nC\tseed\n")
    seeds = tmp_path / "seeds.tsv"
    seeds.write_text("Q\tseed\t1e-20\t100\nC\tseed\t1e-10\t80\n")
    result = convert_eggnog(
        input_path=raw,
        proteome_fasta=_fasta(tmp_path / "proteome.faa", ("Q", "C", "N")),
        query_id="Q",
        source_version="2.1.15",
        database_version="5.0.2",
        command="auto",
        hits_path=hits,
        seed_orthologs_path=seeds,
    )
    statuses = {row["protein_id"]: row["parse_status"] for row in result.audit_rows}
    assert statuses == {"C": "hit_no_annotation", "N": "no_hit", "Q": "annotated"}
    assert dict((name, count) for name, count, _ in result.coverage_rows)["no_hit"] == 1


def test_deterministic_lf_outputs_versions_command_and_checksums(tmp_path: Path) -> None:
    raw = _raw(tmp_path / "raw.tsv", (_row("Q"),))
    kwargs: dict[str, Any] = {
        "input_path": raw,
        "proteome_fasta": _fasta(tmp_path / "proteome.faa", ("Q",)),
        "query_id": "Q",
        "source_version": "2.1.15",
        "database_version": "5.0.2",
        "command": "auto",
    }
    result = convert_eggnog(**kwargs)
    audit = tmp_path / "audit.tsv"
    coverage = tmp_path / "coverage.tsv"
    metadata = tmp_path / "metadata.json"
    write_conversion(
        result,
        audit_output=audit,
        coverage_output=coverage,
        metadata_output=metadata,
    )
    first = (audit.read_bytes(), coverage.read_bytes(), metadata.read_bytes())
    write_conversion(
        convert_eggnog(**kwargs),
        audit_output=audit,
        coverage_output=coverage,
        metadata_output=metadata,
    )
    assert first == (audit.read_bytes(), coverage.read_bytes(), metadata.read_bytes())
    assert b"\r" not in audit.read_bytes()
    parsed = json.loads(metadata.read_text())
    assert parsed["source_version"] == "2.1.15"
    assert parsed["database_version"] == "5.0.2"
    assert parsed["command"].startswith("/tool/emapper.py")
    assert len(parsed["raw_source_sha256"]) == 64
