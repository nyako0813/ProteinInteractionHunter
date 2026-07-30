from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest
from openpyxl import Workbook

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.exceptions import InputValidationError  # noqa: E402
from scripts.audit_coverage_ab import audit_ab  # noqa: E402


def _write_fixture(path: Path, *, run_id: str, rank: str = "") -> None:
    path.mkdir()
    fields = (
        "run_id",
        "protein_id",
        "functional_status",
        "scoring_status",
        "evidence_tier_status",
        "integrated_score",
        "rank",
        "evidence_tier",
    )
    with (path / "candidate_table.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "run_id": run_id,
                "protein_id": "P1",
                "functional_status": "",
                "scoring_status": "",
                "evidence_tier_status": "",
                "integrated_score": "",
                "rank": rank,
                "evidence_tier": "",
            }
        )
    bundle = {
        "run_id": run_id,
        "protein_id": "P1",
        "engine_statuses": {
            "functional_complementarity": "not_run",
            "scoring": "not_run",
            "evidence_tiers": "not_run",
        },
        "functional": None,
        "score": {"total_ranking_score": None},
        "evidence_tier": None,
    }
    (path / "candidate_evidence_bundle.jsonl").write_text(
        json.dumps(bundle) + "\n", encoding="utf-8"
    )
    (path / "warning_summary.tsv").write_text("warning\tcount\nnone\t0\n", encoding="utf-8")
    workbook = Workbook()
    worksheet = workbook.active
    assert worksheet is not None
    worksheet.title = "Candidates"
    worksheet.append(("protein_id",))
    worksheet.append(("P1",))
    workbook.save(path / "ProteinInteractionHunter.xlsx")
    (path / "run_manifest.json").write_text(
        json.dumps({"run_id": run_id, "stable": "same"}), encoding="utf-8"
    )


def test_coverage_ab_ignores_run_metadata_and_requires_disabled_outputs(
    tmp_path: Path,
) -> None:
    baseline = tmp_path / "baseline"
    comparison = tmp_path / "comparison"
    _write_fixture(baseline, run_id="old")
    _write_fixture(comparison, run_id="new")
    rows = audit_ab(baseline, comparison)
    assert all(value == "True" for check, value, _ in rows if check != "manifest_changed_keys")
    details = {check: detail for check, _, detail in rows}
    assert details["manifest_changed_keys"] == "run_id"


def test_coverage_ab_rejects_ranked_comparison(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    comparison = tmp_path / "comparison"
    _write_fixture(baseline, run_id="old")
    _write_fixture(comparison, run_id="new", rank="1")
    with pytest.raises(InputValidationError, match="A/B comparison failed"):
        audit_ab(baseline, comparison)
