from __future__ import annotations

import json
from pathlib import Path

import yaml

from protein_interaction_hunter.application.pipeline import InteractionCandidatePipeline
from protein_interaction_hunter.application.validation import validate_local_inputs
from protein_interaction_hunter.config import load_config

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def test_disabled_local_tables_are_validated_and_fingerprinted_without_propagation(
    tmp_path: Path,
) -> None:
    raw = yaml.safe_load((FIXTURES / "config.valid.yaml").read_text())
    for key in ("proteome_fasta", "genome_gff", "annotation_table"):
        raw["input"][key] = str((FIXTURES / raw["input"][key]).resolve())
    raw["domains"]["local_table"] = str((FIXTURES / raw["domains"]["local_table"]).resolve())
    raw["domains"]["rules_path"] = str((FIXTURES / raw["domains"]["rules_path"]).resolve())
    raw["functional_complementarity"]["rules_path"] = str(
        (FIXTURES / raw["functional_complementarity"]["rules_path"]).resolve()
    )
    raw["fusion"]["local_table"] = str((FIXTURES / raw["fusion"]["local_table"]).resolve())
    raw["known_interactions"]["local_table"] = str(
        (FIXTURES / raw["known_interactions"]["local_table"]).resolve()
    )
    raw["orthology"]["enabled"] = False
    raw["orthology"]["local_table"] = str((FIXTURES / "synthetic_orthology.tsv").resolve())
    raw["phylogenetic_profile"]["enabled"] = False
    raw["phylogenetic_profile"]["local_table"] = str(
        (FIXTURES / "synthetic_phylogenetic_profiles.tsv").resolve()
    )
    raw["scoring"]["enabled"] = False
    raw["evidence_tiers"]["enabled"] = False
    raw["output"]["directory"] = str(tmp_path / "output")
    raw["cache"]["directory"] = str(tmp_path / "cache")
    raw["logging"]["directory"] = str(tmp_path / "logs")
    config_path = tmp_path / "coverage.yaml"
    config_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")

    summary = validate_local_inputs(load_config(config_path))
    assert summary.orthology_annotation_count > 0
    assert summary.orthology_protein_count > 0
    assert summary.unknown_orthology_id_count == 0
    assert summary.phylogenetic_profile_observation_count > 0
    assert summary.phylogenetic_profile_protein_count > 0
    assert summary.phylogenetic_profile_species_count > 0
    assert summary.unknown_phylogenetic_profile_id_count == 0

    result = InteractionCandidatePipeline().run(config_path)
    manifest = json.loads(result.manifest_path.read_text())
    logical_names = {row["logical_name"]: row for row in manifest["input_files"]}
    assert logical_names["orthology_local_table"]["required"] is False
    assert logical_names["phylogenetic_profile_local_table"]["required"] is False
    bundles = [json.loads(line) for line in result.evidence_path.read_text().splitlines()]
    assert all(bundle["orthology"] == [] for bundle in bundles)
    assert all(bundle["phylogenetic_profile"] == [] for bundle in bundles)
    assert all(bundle["engine_statuses"]["orthology"] == "not_run" for bundle in bundles)
    assert all(bundle["engine_statuses"]["phylogenetic_profile"] == "not_run" for bundle in bundles)
    assert all(bundle["score"]["total_ranking_score"] is None for bundle in bundles)
    assert all(bundle["evidence_tier"] is None for bundle in bundles)
