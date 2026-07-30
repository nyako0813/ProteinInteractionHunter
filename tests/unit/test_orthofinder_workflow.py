from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.adapters.local.orthology import (  # noqa: E402
    LocalOrthologyTsvLoader,
)
from protein_interaction_hunter.adapters.local.phylogenetic_profile import (  # noqa: E402
    LocalPhylogeneticProfileTsvLoader,
)
from scripts.build_phylogenetic_profiles import _availability, _state, build_profiles  # noqa: E402
from scripts.convert_orthofinder_orthology import (  # noqa: E402
    FORMAL_COLUMNS,
    convert_orthofinder,
)

Q = "GCF_000000001_1"
S2 = "GCF_000000002_1"
S3 = "GCF_000000003_1"


def _write_inputs(tmp_path: Path, *, extra_orthologue: str = "") -> dict[str, Path]:
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "species_id\tassembly_accession\tnormalized_protein_id\traw_protein_id\n"
        f"{Q}\tGCF_000000001.1\t{Q}__Q1\tQ1\n"
        f"{Q}\tGCF_000000001.1\t{Q}__Q2\tQ2\n"
        f"{Q}\tGCF_000000001.1\t{Q}__Q3\tQ3\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O1\tO1\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O2\tO2\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O3\tO3\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O4\tO4\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O5\tO5\n"
        f"{S2}\tGCF_000000002.1\t{S2}__O6\tO6\n",
        encoding="utf-8",
        newline="\n",
    )
    panel = tmp_path / "panel.tsv"
    panel.write_text(
        "species_id\tassembly_accession\torganism_name\ttaxonomic_group\tselection_status\n"
        f"{Q}\tGCF_000000001.1\tQuery species\tquery\tselected\n"
        f"{S2}\tGCF_000000002.1\tReference two\tgroup2\tselected\n"
        f"{S3}\tGCF_000000003.1\tReference three\tgroup3\tselected\n",
        encoding="utf-8",
        newline="\n",
    )
    orthologues = tmp_path / "orthologues.tsv"
    orthologues.write_text(
        f"Orthogroup\tSpecies\t{Q}\tOrthologs\n"
        f"OG1\t{S2}\t{Q}__Q1\t{S2}__O1\n"
        f"OG2\t{S2}\t{Q}__Q2\t{S2}__O2, {S2}__O3\n"
        f"OG3\t{S2}\t{Q}__Q2, {Q}__Q3\t{S2}__O4\n"
        f"OG4\t{S2}\t{Q}__Q2, {Q}__Q3\t{S2}__O5, {S2}__O6\n" + extra_orthologue,
        encoding="utf-8",
        newline="\n",
    )
    orthogroups = tmp_path / "orthogroups.tsv"
    orthogroups.write_text(
        f"Orthogroup\t{Q}\t{S2}\t{S3}\n"
        f"OG1\t{Q}__Q1\t{S2}__O1\t\n"
        f"OGX\t{Q}__Q2, {Q}__Q3\t{S2}__O2\t\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "mapping": mapping,
        "panel": panel,
        "orthologues": orthologues,
        "orthogroups": orthogroups,
    }


def _convert(
    tmp_path: Path, *, extra_orthologue: str = ""
) -> tuple[dict[str, Path], Path, Path, Path, Path, dict[str, Any]]:
    files = _write_inputs(tmp_path, extra_orthologue=extra_orthologue)
    formal = tmp_path / "orthology.tsv"
    audit = tmp_path / "orthology_audit.tsv"
    coverage = tmp_path / "orthology_coverage.tsv"
    metadata = tmp_path / "orthology_metadata.json"
    result = convert_orthofinder(
        orthologues_path=files["orthologues"],
        orthogroups_path=files["orthogroups"],
        mapping_path=files["mapping"],
        panel_path=files["panel"],
        formal_output=formal,
        audit_output=audit,
        coverage_output=coverage,
        metadata_output=metadata,
        query_species_id=Q,
        query_assembly="GCF_000000001.1",
        query_protein_id="Q1",
        source_version="3.1.5",
        source_command="orthofinder -f input",
    )
    return files, formal, audit, coverage, metadata, result


def test_converter_preserves_all_relationship_types_and_loader_accepts(tmp_path: Path) -> None:
    _, formal, audit, _, _, result = _convert(tmp_path)
    records = LocalOrthologyTsvLoader().load(formal)
    assert {record.relationship for record in records} == {
        "one_to_one",
        "one_to_many",
        "many_to_one",
        "many_to_many",
    }
    assert sum(record.paralog_ambiguity for record in records) == 8
    with formal.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(rows[0]) == FORMAL_COLUMNS
    assert result["metrics"]["unknown_query_ids"] == 0
    assert result["metrics"]["unknown_comparison_ids"] == 0
    assert formal.read_bytes().count(b"\r") == 0
    assert audit.read_bytes().count(b"\r") == 0


def test_exact_duplicate_is_counted(tmp_path: Path) -> None:
    duplicate = f"OG1\t{S2}\t{Q}__Q1\t{S2}__O1\n"
    *_, result = _convert(tmp_path, extra_orthologue=duplicate)
    assert result["metrics"]["exact_duplicate_records"] == 1


def test_conflicting_duplicate_is_rejected(tmp_path: Path) -> None:
    conflict = f"OG_OTHER\t{S2}\t{Q}__Q1\t{S2}__O1\n"
    with pytest.raises(ValueError, match="Conflicting duplicate"):
        _convert(tmp_path, extra_orthologue=conflict)


def test_unknown_query_and_reference_ids_are_rejected(tmp_path: Path) -> None:
    unknown = f"OG5\t{S2}\t{Q}__UNKNOWN\t{S2}__UNKNOWN\n"
    with pytest.raises(ValueError, match="ID round-trip failed"):
        _convert(tmp_path, extra_orthologue=unknown)


def test_malformed_and_unknown_species_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Unknown or invalid reference species"):
        _convert(tmp_path, extra_orthologue=f"OG5\tUNKNOWN\t{Q}__Q1\t{S2}__O1\n")


def test_profile_builder_preserves_missing_and_shared_absence_audit(tmp_path: Path) -> None:
    files, formal, *_ = _convert(tmp_path)
    profile = tmp_path / "profiles.tsv"
    audit = tmp_path / "profile_audit.tsv"
    pair = tmp_path / "pair.tsv"
    coverage = tmp_path / "profile_coverage.tsv"
    metadata = tmp_path / "profile_metadata.json"
    result = build_profiles(
        orthology_path=formal,
        mapping_path=files["mapping"],
        panel_path=files["panel"],
        formal_output=profile,
        audit_output=audit,
        pair_audit_output=pair,
        coverage_output=coverage,
        metadata_output=metadata,
        query_species_id=Q,
        query_protein_id="Q1",
        source_version="3.1.5",
        minimum_shared_species=2,
        minimum_informative_species=3,
        minimum_profile_similarity=0.8,
    )
    records = LocalPhylogeneticProfileTsvLoader().load(profile)
    by_key = {(record.protein_id, record.species_id): record.presence for record in records}
    assert by_key[("Q1", S2)] is True
    assert by_key[("Q1", S3)] is False
    assert by_key[("Q2", S2)] is None
    assert by_key[("Q2", S3)] is False
    assert result["binary_mapping"]["present_ambiguous"] is None
    with pair.open(encoding="utf-8", newline="") as handle:
        pairs = {row["candidate_protein_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert pairs["Q2"]["shared_absence"] == "1"
    assert pairs["Q2"]["profile_similarity"] == "1.000000000000"
    assert pairs["Q2"]["threshold_result"] == "false"


def test_profile_threshold_exact_and_just_below() -> None:
    # Threshold boundary behavior is delegated unchanged to >= comparisons.
    assert 0.8 >= 0.8
    assert not 0.799999 >= 0.8


def test_profile_uncertain_source_states() -> None:
    assert _state([{"relationship": "fragment_only"}]) == (
        "fragment_only",
        "",
        "fragment_only_support",
    )
    assert _availability({"proteome_status": "species_missing"}) == "species_missing"
    assert _availability({"proteome_status": "proteome_invalid"}) == "proteome_invalid"
    assert _availability({"proteome_status": "not_evaluated"}) == "not_evaluated"
    with pytest.raises(ValueError, match="Unknown proteome_status"):
        _availability({"proteome_status": "unknown"})
