"""Candidate pipeline artifacts, provenance, and CLI end-to-end behavior."""

import csv
import json
from pathlib import Path

import yaml
from typer.testing import CliRunner

from protein_interaction_hunter.application.gene_context import GENE_CONTEXT_RULE_VERSION
from protein_interaction_hunter.application.identifiers import NORMALIZATION_RULE_VERSION
from protein_interaction_hunter.application.pipeline import InteractionCandidatePipeline
from protein_interaction_hunter.cli import app
from protein_interaction_hunter.models import CandidateEvidenceBundle, EvidenceStatus, RunManifest
from protein_interaction_hunter.outputs.candidates import CANDIDATE_COLUMNS
from protein_interaction_hunter.outputs.jsonl import JsonlEvidenceBundleWriter


def e2e_config(valid_config_path: Path, tmp_path: Path) -> Path:
    raw = yaml.safe_load(valid_config_path.read_text(encoding="utf-8"))
    fixture_dir = valid_config_path.parent
    raw["input"]["proteome_fasta"] = str(fixture_dir / "synthetic_proteome.fasta")
    raw["input"]["genome_gff"] = str(fixture_dir / "synthetic_genome.gff3")
    raw["input"]["annotation_table"] = str(fixture_dir / "synthetic_annotations.tsv")
    raw["output"]["directory"] = str(tmp_path / "generated")
    path = tmp_path / "e2e.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=True), encoding="utf-8")
    return path


def test_pipeline_gene_context_scores_and_jsonl_round_trip(
    valid_config_path: Path, tmp_path: Path
) -> None:
    result = InteractionCandidatePipeline().run(e2e_config(valid_config_path, tmp_path))
    assert len(result.bundles) == 12
    for bundle in result.bundles:
        assert len(bundle.genome_context) == 1
        assert bundle.engine_statuses["gene_context"] is bundle.genome_context[0].status
        assert all(
            status is EvidenceStatus.NOT_RUN
            for engine, status in bundle.engine_statuses.items()
            if engine != "gene_context"
        )
        scores = bundle.score.model_dump()
        assert all(value is None for key, value in scores.items() if key.endswith("_score"))
        assert bundle.score.contradiction_penalty is None
        assert bundle.score.evidence_completeness is None
        assert bundle.evidence_tier is None
    assert JsonlEvidenceBundleWriter().read(result.evidence_path) == result.bundles
    assert all(
        CandidateEvidenceBundle.model_validate_json(line)
        for line in result.evidence_path.read_text(encoding="utf-8").splitlines()
    )


def test_candidate_tsv_header_and_manifest_provenance(
    valid_config_path: Path, tmp_path: Path
) -> None:
    config = e2e_config(valid_config_path, tmp_path)
    result = InteractionCandidatePipeline().run(config)
    with result.candidate_table_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        assert tuple(next(reader)) == CANDIDATE_COLUMNS
        assert len(list(reader)) == 12
    with result.candidate_table_path.open(encoding="utf-8", newline="") as handle:
        rows_by_id = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows_by_id["NEAR_001"]["distance_bp"] == "29"
    assert rows_by_id["CONTIG2_001"]["distance_bp"] == ""
    assert rows_by_id["FRAG_001"]["gene_context_status"] == "missing"
    manifest = RunManifest.model_validate_json(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest.config_sha256
    assert all(item.sha256 for item in manifest.input_files)
    assert manifest.package_version == "0.1.0"
    assert manifest.git_commit
    assert manifest.normalization_rule_version == NORMALIZATION_RULE_VERSION
    assert manifest.gene_context_rule_version == GENE_CONTEXT_RULE_VERSION
    assert manifest.policy_settings["fragment_policy"] == "flag"
    assert "scoring_not_run" in manifest.incomplete_evidence_flags
    assert "gene_context_not_run" not in manifest.incomplete_evidence_flags
    assert result.config_snapshot_path.is_file()
    assert result.warning_summary_path.is_file()
    assert result.excel_path is not None
    assert result.excel_path.is_file()


def test_generate_candidates_cli_e2e_has_no_ranking_output(
    valid_config_path: Path, tmp_path: Path
) -> None:
    config = e2e_config(valid_config_path, tmp_path)
    result = CliRunner().invoke(app, ["generate-candidates", "--config", str(config)])
    assert result.exit_code == 0, result.stdout
    assert "query_count: 1" in result.stdout
    assert "protein_count: 12" in result.stdout
    assert "query_candidate_pair_count: 12" in result.stdout
    assert "duplicate_group_count: 1" in result.stdout
    assert "same_contig_pair_count:" in result.stdout
    assert "different_contig_pair_count: 1" in result.stdout
    assert "overlapping_pair_count:" in result.stdout
    assert "output_path:" in result.stdout
    assert "ranking" not in result.stdout.casefold()
    output = tmp_path / "generated"
    assert {path.name for path in output.iterdir()} == {
        "ProteinInteractionHunter.xlsx",
        "candidate_evidence_bundle.jsonl",
        "candidate_table.tsv",
        "config.snapshot.yaml",
        "run_manifest.json",
        "warning_summary.tsv",
    }
    first = (output / "candidate_table.tsv").read_text(encoding="utf-8")
    rerun = CliRunner().invoke(app, ["generate-candidates", "--config", str(config)])
    assert rerun.exit_code == 0
    assert (output / "candidate_table.tsv").read_text(encoding="utf-8") == first
    payload = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["normalization_rule_version"] == NORMALIZATION_RULE_VERSION
    assert payload["gene_context_rule_version"] == GENE_CONTEXT_RULE_VERSION
