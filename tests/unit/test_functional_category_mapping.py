from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.adapters.local.annotation import (  # noqa: E402
    ANNOTATION_COLUMNS,
)
from protein_interaction_hunter.exceptions import InputValidationError  # noqa: E402
from scripts.build_functional_category_table import (  # noqa: E402
    CATEGORY_FIELDS,
    MAPPING_FIELDS,
    SourceTerm,
    build_functional_table,
    load_mapping_policy,
    map_source_term,
)
from scripts.convert_eggnog_annotations import AUDIT_COLUMNS  # noqa: E402
from scripts.functional_annotation_common import (  # noqa: E402
    ancestor_paths,
    parse_go_obo,
)


def _obo(path: Path, cycle: bool = False) -> Path:
    cycle_terms = (
        """
[Term]
id: GO:0000100
name: cycle a
namespace: biological_process
is_a: GO:0000101 ! cycle b

[Term]
id: GO:0000101
name: cycle b
namespace: biological_process
is_a: GO:0000100 ! cycle a
"""
        if cycle
        else ""
    )
    path.write_text(
        """format-version: 1.2
data-version: synthetic-v1

[Term]
id: GO:0003674
name: molecular function
namespace: molecular_function

[Term]
id: GO:0008150
name: biological process
namespace: biological_process

[Term]
id: GO:0005575
name: cellular component
namespace: cellular_component

[Term]
id: GO:0000001
name: parent function
namespace: molecular_function
is_a: GO:0003674 ! molecular function

[Term]
id: GO:0000002
name: child function
namespace: molecular_function
is_a: GO:0000001 ! parent function
relationship: part_of GO:0000003 ! process
relationship: regulates GO:0000004 ! regulated process

[Term]
id: GO:0000003
name: part process
namespace: biological_process
is_a: GO:0008150 ! biological process

[Term]
id: GO:0000004
name: regulated process
namespace: biological_process
is_a: GO:0008150 ! biological process

[Term]
id: GO:0000005
name: component child
namespace: cellular_component
is_a: GO:0005575 ! cellular component

[Term]
id: GO:0000006
name: obsolete function
namespace: molecular_function
is_obsolete: true
replaced_by: GO:0000001
consider: GO:0000002
"""
        + cycle_terms
    )
    return path


def _category() -> dict[str, object]:
    category: dict[str, object] = dict.fromkeys(CATEGORY_FIELDS, "")
    category.update(
        {
            "category_id": "rna_role",
            "label": "RNA role",
            "definition": "Synthetic role.",
            "inclusion_criteria": "Accepted mapping.",
            "exclusion_criteria": "Broad terms.",
            "allowed_source_types": ["go", "kegg_ko"],
            "accepted_go_aspects": ["molecular_function", "biological_process"],
            "accepted_go_ancestor_relations": ["is_a", "part_of"],
            "taxonomic_applicability": "all",
            "evidence_quality": "curated",
            "ambiguity_policy": "omit",
            "conflict_policy": "fail",
        }
    )
    return category


def _mapping(**updates: object) -> dict[str, object]:
    mapping: dict[str, object] = dict.fromkeys(MAPPING_FIELDS, "")
    mapping.update(
        {
            "mapping_id": "map-1",
            "formal_category": "rna_role",
            "source_type": "go",
            "source_identifier": "GO:0000001",
            "source_aspect": "molecular_function",
            "match_type": "ancestor",
            "allowed_relations": ["is_a"],
            "minimum_specificity": "child",
            "priority": 10,
            "mapping_status": "accepted",
            "rationale": "Synthetic rationale.",
            "reference": ["https://example.test/reference"],
            "rule_version": "synthetic-mapping-v1",
        }
    )
    mapping.update(updates)
    return mapping


def _policy(path: Path, mappings: list[dict[str, object]]) -> Path:
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1",
                "mapping_rule_version": "synthetic-mapping-v1",
                "go_release_folder": "synthetic",
                "go_data_version": "synthetic-v1",
                "go_sha256": "0" * 64,
                "default_allowed_go_relations": ["is_a"],
                "formal_categories": [_category()],
                "mappings": mappings,
            },
            sort_keys=False,
        )
    )
    return path


def test_go_exact_ancestors_multiple_paths_and_relation_policy(tmp_path: Path) -> None:
    ontology = parse_go_obo(_obo(tmp_path / "go.obo"))
    is_a = ancestor_paths(ontology, "GO:0000002", frozenset({"is_a"}))
    assert "GO:0000001" in is_a
    assert "GO:0000003" not in is_a
    part_of = ancestor_paths(
        ontology,
        "GO:0000002",
        frozenset({"is_a", "part_of"}),
    )
    assert "GO:0000003" in part_of
    assert "GO:0000004" not in part_of
    assert ontology.terms["GO:0000006"].replaced_by == ("GO:0000001",)
    assert ontology.terms["GO:0000006"].consider == ("GO:0000002",)


def test_go_cycle_detection(tmp_path: Path) -> None:
    ontology = parse_go_obo(_obo(tmp_path / "go.obo", cycle=True))
    with pytest.raises(InputValidationError, match="cycle detected"):
        ancestor_paths(ontology, "GO:0000100", frozenset({"is_a"}))


def test_exact_ancestor_cc_obsolete_unknown_and_description(tmp_path: Path) -> None:
    ontology = parse_go_obo(_obo(tmp_path / "go.obo"))
    policy = load_mapping_policy(
        _policy(tmp_path / "policy.yaml", [_mapping()]),
        ontology,
    )
    ancestor = map_source_term(
        SourceTerm("Q", "test", "go", "GO:0000002", "child", "molecular_function"),
        policy,
        ontology,
    )
    assert ancestor[0]["mapping_status"] == "accepted"
    assert ancestor[0]["match_type"] == "ancestor"
    cc = map_source_term(
        SourceTerm("Q", "test", "go", "GO:0000005", "cc", "cellular_component"),
        policy,
        ontology,
    )
    assert cc[0]["mapping_status"] == "excluded"
    obsolete = map_source_term(
        SourceTerm("Q", "test", "go", "GO:0000006", "old", "molecular_function"),
        policy,
        ontology,
    )
    assert obsolete[0]["exclusion_reason"] == (
        "obsolete_go_id;replaced_by=GO:0000001;consider=GO:0000002"
    )
    unknown = map_source_term(
        SourceTerm("Q", "test", "go", "GO:9999999", "", "unknown"),
        policy,
        ontology,
    )
    assert unknown[0]["exclusion_reason"] == "unknown_go_id"
    description = map_source_term(
        SourceTerm("Q", "test", "description", "ATPase", "ATPase", "not_applicable"),
        policy,
        ontology,
    )
    assert description[0]["mapping_status"] == "manual_review"


@pytest.mark.parametrize(
    "mutator, message",
    (
        (lambda item: item.update(mapping_id=""), "empty mapping_id"),
        (lambda item: item.update(formal_category="unknown"), "Unknown formal category"),
        (lambda item: item.update(reference=[]), "has no reference"),
        (lambda item: item.update(rationale=""), "empty rationale"),
        (lambda item: item.update(allowed_relations=["regulates"]), "Invalid GO ancestor"),
        (lambda item: item.update(source_identifier="GO:0003674"), "Overbroad GO root"),
        (lambda item: item.update(source_aspect="biological_process"), "aspect mismatch"),
        (lambda item: item.update(rule_version="wrong"), "wrong rule version"),
    ),
)
def test_mapping_policy_rejects_invalid_records(
    tmp_path: Path,
    mutator: object,
    message: str,
) -> None:
    ontology = parse_go_obo(_obo(tmp_path / "go.obo"))
    mapping = _mapping()
    mutator(mapping)  # type: ignore[operator]
    with pytest.raises(InputValidationError, match=message):
        load_mapping_policy(_policy(tmp_path / "policy.yaml", [mapping]), ontology)


def test_mapping_policy_rejects_duplicate_id_priority_tie_and_conflict(
    tmp_path: Path,
) -> None:
    ontology = parse_go_obo(_obo(tmp_path / "go.obo"))
    duplicate = [_mapping(), _mapping()]
    with pytest.raises(InputValidationError, match="Duplicate mapping_id"):
        load_mapping_policy(_policy(tmp_path / "duplicate.yaml", duplicate), ontology)

    tied = [_mapping(), _mapping(mapping_id="map-2")]
    with pytest.raises(InputValidationError, match="Priority tie"):
        load_mapping_policy(_policy(tmp_path / "tie.yaml", tied), ontology)

    second_category = _category()
    second_category["category_id"] = "other_role"
    policy_path = _policy(
        tmp_path / "conflict.yaml",
        [
            _mapping(),
            _mapping(
                mapping_id="map-2",
                formal_category="other_role",
                priority=9,
            ),
        ],
    )
    raw = yaml.safe_load(policy_path.read_text(encoding="utf-8"))
    raw["formal_categories"].append(second_category)
    policy_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    with pytest.raises(InputValidationError, match="Conflicting accepted mappings"):
        load_mapping_policy(policy_path, ontology)


def _write_dict_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_formal_builder_preserves_ncbi_and_localization_fields(tmp_path: Path) -> None:
    annotation = dict.fromkeys(ANNOTATION_COLUMNS, "")
    annotation.update(
        {
            "protein_id": "Q",
            "gene_name": "gene",
            "locus_tag": "locus",
            "product": "product",
            "localization_annotation": "cytosolic",
            "annotation_source": "NCBI",
        }
    )
    _write_dict_tsv(tmp_path / "annotation.tsv", ANNOTATION_COLUMNS, [annotation])
    eggnog = dict.fromkeys(AUDIT_COLUMNS, "")
    eggnog.update(
        {
            "protein_id": "Q",
            "KEGG_ko": "ko:K1",
            "source": "eggNOG-mapper",
            "source_version": "2.1.15",
            "database_version": "5.0.2",
            "parse_status": "annotated",
        }
    )
    _write_dict_tsv(tmp_path / "eggnog.tsv", AUDIT_COLUMNS, [eggnog])
    policy = _policy(
        tmp_path / "policy.yaml",
        [
            _mapping(
                source_type="kegg_ko",
                source_identifier="ko:K1",
                source_aspect="not_applicable",
                match_type="exact",
                allowed_relations=[],
            )
        ],
    )
    build = build_functional_table(
        annotation_table=tmp_path / "annotation.tsv",
        eggnog_audit=tmp_path / "eggnog.tsv",
        go_obo=_obo(tmp_path / "go.obo"),
        mapping_policy=policy,
        query_id="Q",
    )
    output = build.annotation_rows[0]
    assert output["functional_category"] == "rna_role"
    assert output["gene_name"] == "gene"
    assert output["locus_tag"] == "locus"
    assert output["product"] == "product"
    assert output["localization_annotation"] == "cytosolic"
    assert output["annotation_source"] == "NCBI"
    assert build.metadata["query_mapping_status"] == "ready"
