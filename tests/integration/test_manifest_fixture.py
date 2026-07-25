"""Expected run manifest fixture validation."""

from pathlib import Path

from protein_interaction_hunter.models import RunManifest


def test_expected_manifest_parses_with_aware_timestamp(fixture_dir: Path) -> None:
    manifest = RunManifest.model_validate_json(
        (fixture_dir / "expected" / "run_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest.started_at.utcoffset() is not None
    assert manifest.incomplete_evidence_flags == ["all_analysis_engines_not_run"]
