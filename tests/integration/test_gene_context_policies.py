"""MVP-1B pipeline policy and alias behavior."""

from pathlib import Path
from typing import Any

import pytest
import yaml

from protein_interaction_hunter.application.pipeline import InteractionCandidatePipeline
from protein_interaction_hunter.exceptions import InputValidationError
from protein_interaction_hunter.models import EvidenceStatus


def config_copy(
    source: Path,
    tmp_path: Path,
    name: str,
) -> tuple[dict[str, Any], Path]:
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    fixture_dir = source.parent

    raw["input"]["proteome_fasta"] = str((fixture_dir / raw["input"]["proteome_fasta"]).resolve())
    raw["input"]["genome_gff"] = str((fixture_dir / raw["input"]["genome_gff"]).resolve())

    annotation_table = raw["input"].get("annotation_table")
    if annotation_table is not None:
        raw["input"]["annotation_table"] = str((fixture_dir / annotation_table).resolve())

    domain_table = raw["domains"].get("local_table")
    if domain_table is not None:
        raw["domains"]["local_table"] = str((fixture_dir / domain_table).resolve())

    orthology_table = raw["orthology"].get("local_table")
    if orthology_table is not None:
        raw["orthology"]["local_table"] = str((fixture_dir / orthology_table).resolve())

    profile_table = raw["phylogenetic_profile"].get("local_table")
    if profile_table is not None:
        raw["phylogenetic_profile"]["local_table"] = str((fixture_dir / profile_table).resolve())

    fusion_table = raw["fusion"].get("local_table")
    if fusion_table is not None:
        raw["fusion"]["local_table"] = str((fixture_dir / fusion_table).resolve())

    domain_rules_path = raw["domains"].get("rules_path")
    if domain_rules_path is not None:
        raw["domains"]["rules_path"] = str((fixture_dir / domain_rules_path).resolve())

    functional_rules_path = raw["functional_complementarity"].get("rules_path")
    if functional_rules_path is not None:
        raw["functional_complementarity"]["rules_path"] = str(
            (fixture_dir / functional_rules_path).resolve()
        )

    return raw, tmp_path / f"{name}.yaml"


def write_config(raw: dict[str, Any], path: Path) -> Path:
    path.write_text(yaml.safe_dump(raw, sort_keys=True), encoding="utf-8")
    return path


def test_disabled_gene_context_remains_not_run(valid_config_path: Path, tmp_path: Path) -> None:
    raw, path = config_copy(valid_config_path, tmp_path, "disabled")
    raw["gene_context"]["enabled"] = False
    result = InteractionCandidatePipeline().run(write_config(raw, path))
    assert all(not bundle.genome_context for bundle in result.bundles)
    assert all(
        bundle.engine_statuses["gene_context"] is EvidenceStatus.NOT_RUN
        for bundle in result.bundles
    )
    assert result.summary.same_contig_pair_count == 0


def test_missing_query_coordinate_is_fatal_when_required(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw, path = config_copy(valid_config_path, tmp_path, "required_missing")
    raw["query"]["protein_ids"] = ["FRAG_001"]
    raw["gene_context"]["require_query_coordinates"] = True
    with pytest.raises(InputValidationError, match="Query coordinate is missing"):
        InteractionCandidatePipeline().run(write_config(raw, path))


def test_missing_query_coordinate_is_represented_when_allowed(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw, path = config_copy(valid_config_path, tmp_path, "allowed_missing")
    raw["query"]["protein_ids"] = ["FRAG_001"]
    result = InteractionCandidatePipeline().run(write_config(raw, path))
    assert all(
        bundle.engine_statuses["gene_context"] is EvidenceStatus.MISSING
        for bundle in result.bundles
    )


def test_ambiguous_query_coordinate_is_fatal_when_required(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw, path = config_copy(valid_config_path, tmp_path, "required_ambiguous")
    source = Path(raw["input"]["genome_gff"])
    ambiguous_gff = tmp_path / "ambiguous.gff3"
    ambiguous_gff.write_text(
        source.read_text(encoding="utf-8")
        + "contig1\tsynthetic\tCDS\t800\t850\t.\t+\t0\t"
        + "ID=cds_query_second;protein_id=QUERY_001\n",
        encoding="utf-8",
    )
    raw["input"]["genome_gff"] = str(ambiguous_gff)
    raw["gene_context"]["require_query_coordinates"] = True
    with pytest.raises(InputValidationError, match="Query coordinate is ambiguous"):
        InteractionCandidatePipeline().run(write_config(raw, path))


def test_query_alias_uses_canonical_query_coordinate(
    valid_config_path: Path, tmp_path: Path
) -> None:
    raw, path = config_copy(valid_config_path, tmp_path, "query_alias")
    raw["query"]["protein_ids"] = ["LT0001"]
    result = InteractionCandidatePipeline().run(write_config(raw, path))
    near = next(bundle for bundle in result.bundles if bundle.candidate_id == "NEAR_001")
    assert near.query_id == "LT0001"
    assert near.genome_context[0].status is EvidenceStatus.AVAILABLE
    assert near.genome_context[0].query_start == 50
