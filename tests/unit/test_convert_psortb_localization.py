import csv
import hashlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.adapters.local.annotation import (  # noqa: E402
    ANNOTATION_COLUMNS,
    LocalAnnotationTsvLoader,
)
from protein_interaction_hunter.exceptions import InputValidationError  # noqa: E402
from scripts.convert_psortb_localization import (  # noqa: E402
    AUDIT_COLUMNS,
    PSORTB_ARCHAEA_MAPPING,
    PSORTB_MAPPING_RULE_VERSION,
    convert_psortb,
    parse_psortb_terse,
    write_conversion_outputs,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def _fasta(path: Path, ids: tuple[str, ...] = ("Q", "A", "B")) -> Path:
    return _write(
        path,
        "".join(f">{protein_id} product text\n{'M' * 20}\n" for protein_id in ids),
    )


def _base_annotation(path: Path, ids: tuple[str, ...] = ("Q", "A")) -> Path:
    rows = [
        {
            "protein_id": protein_id,
            "gene_name": "",
            "locus_tag": "",
            "product": "membrane Pfam product words",
            "functional_category": "",
            "localization_annotation": "",
            "transmembrane_annotation": "",
            "annotation_source": "NCBI RefSeq",
            "annotation_confidence": "",
        }
        for protein_id in ids
    ]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ANNOTATION_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_valid_archaeal_terse_maps_explicit_vocabulary_and_preserves_id(
    tmp_path: Path,
) -> None:
    raw = _write(
        tmp_path / "psortb.tsv",
        "SeqID\tLocalization\tScore\n"
        "Q\tCytoplasmic\t9.50\n"
        "A\tCytoplasmicMembrane\t8\n"
        "B\tCellwall\t7.5\n",
    )
    result = convert_psortb(
        input_path=raw,
        proteome_fasta=_fasta(tmp_path / "proteome.faa"),
        query_id="Q",
        source_version="3.0.6",
        source_command="psort -a -o terse input.faa",
    )

    formal = {row["protein_id"]: row for row in result.formal_rows}
    assert formal["Q"]["localization_annotation"] == "cytosolic"
    assert formal["A"]["localization_annotation"] == "membrane"
    assert formal["B"]["localization_annotation"] == ""
    assert next(row for row in result.audit_rows if row.protein_id == "B").mapping_status == (
        "unsupported"
    )
    assert formal["Q"]["transmembrane_annotation"] == ""
    assert result.query_prediction == "Cytoplasmic"
    assert result.audit_rows[0].source_command == "psort -a -o terse input.faa"


def test_empty_header_comments_and_blank_lines_are_accepted(tmp_path: Path) -> None:
    raw = _write(
        tmp_path / "empty.tsv",
        "# PSORTb 3.0.6\n\nSeqId\tLocalization\tScore\n",
    )
    parsed, duplicates = parse_psortb_terse(raw)
    assert parsed == []
    assert duplicates == 0


def test_unknown_is_missing_not_negative_and_product_is_not_inferred(
    tmp_path: Path,
) -> None:
    raw = _write(tmp_path / "unknown.tsv", "Q\tUnknown\t0\n")
    result = convert_psortb(
        input_path=raw,
        proteome_fasta=_fasta(tmp_path / "proteome.faa", ("Q",)),
        base_annotation_table=_base_annotation(tmp_path / "base.tsv", ("Q",)),
        query_id="Q",
    )
    row = result.formal_rows[0]
    assert row["product"] == "membrane Pfam product words"
    assert row["localization_annotation"] == ""
    assert row["transmembrane_annotation"] == ""
    assert result.audit_rows[0].mapping_status == "unknown"
    assert result.audit_rows[0].exclusion_reason == "prediction_unknown"


def test_missing_query_is_fatal(tmp_path: Path) -> None:
    raw = _write(tmp_path / "missing.tsv", "A\tCytoplasmic\t9\n")
    with pytest.raises(InputValidationError, match="Query ID missing"):
        convert_psortb(input_path=raw, query_id="Q")


def test_exact_duplicate_is_audited_separately(tmp_path: Path) -> None:
    raw = _write(
        tmp_path / "duplicate.tsv",
        "Q\tCytoplasmic\t9\nQ\tCytoplasmic\t9.0\n",
    )
    result = convert_psortb(input_path=raw)
    assert result.exact_duplicate_count == 1
    assert len(result.formal_rows) == 1
    assert result.audit_rows[1].duplicate_kind == "exact_duplicate"
    assert result.audit_rows[1].exclusion_reason == "exact_duplicate"


def test_conflicting_duplicate_is_fatal(tmp_path: Path) -> None:
    raw = _write(
        tmp_path / "conflict.tsv",
        "Q\tCytoplasmic\t9\nQ\tExtracellular\t8\n",
    )
    with pytest.raises(InputValidationError, match="Conflicting duplicate"):
        parse_psortb_terse(raw)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("Q\tCytoplasmic\n", "exactly 3 columns"),
        ("Q\tCytoplasmic\tnot-a-score\n", "Invalid PSORTb score"),
        ("Q\tCytoplasmic\t11\n", "between 0 and 10"),
        ("\tCytoplasmic\t9\n", "empty SeqId"),
    ],
)
def test_malformed_rows_are_fatal(tmp_path: Path, text: str, message: str) -> None:
    raw = _write(tmp_path / "malformed.tsv", text)
    with pytest.raises(InputValidationError, match=message):
        parse_psortb_terse(raw)


def test_unknown_and_missing_proteome_ids_are_audited(tmp_path: Path) -> None:
    raw = _write(
        tmp_path / "coverage.tsv",
        "Q\tCytoplasmic\t9\nUNKNOWN\tExtracellular\t8\n",
    )
    result = convert_psortb(
        input_path=raw,
        proteome_fasta=_fasta(tmp_path / "proteome.faa", ("Q", "MISSING")),
    )
    assert result.unknown_protein_id_count == 1
    assert result.missing_protein_id_count == 1
    audit = {row.protein_id: row for row in result.audit_rows}
    assert not audit["UNKNOWN"].known_protein_id
    assert audit["UNKNOWN"].exclusion_reason == "unknown_protein_id"
    metrics = {name: count for name, count, _ in result.coverage}
    assert metrics["known_proteins_represented"] == "1"
    assert metrics["proteins_absent_from_raw_output"] == "1"


def test_unsupported_prediction_is_not_silently_mapped(tmp_path: Path) -> None:
    raw = _write(tmp_path / "unsupported.tsv", "Q\tPeriplasmic\t8\n")
    result = convert_psortb(input_path=raw)
    assert result.formal_rows[0]["localization_annotation"] == ""
    assert result.audit_rows[0].mapping_status == "unsupported"
    assert result.audit_rows[0].exclusion_reason == "unsupported_prediction"


def test_outputs_are_deterministic_loadable_exact_header_lf_and_checksummed(
    tmp_path: Path,
) -> None:
    raw = _write(
        tmp_path / "raw.tsv",
        "B\tExtracellular\t8\nA\tCytoplasmic\t9\nQ\tUnknown\t0\n",
    )
    fasta = _fasta(tmp_path / "proteome.faa", ("Q", "B", "A"))
    base = _base_annotation(tmp_path / "base.tsv", ("Q", "A"))
    result = convert_psortb(
        input_path=raw,
        proteome_fasta=fasta,
        base_annotation_table=base,
        query_id="Q",
        source_version="3.0.6",
    )
    annotation = tmp_path / "annotation.tsv"
    audit = tmp_path / "audit.tsv"
    coverage = tmp_path / "coverage.tsv"
    metadata = tmp_path / "metadata.tsv"
    write_conversion_outputs(
        result,
        annotation_output=annotation,
        audit_output=audit,
        coverage_output=coverage,
        metadata_output=metadata,
    )

    assert [record.protein_id for record in LocalAnnotationTsvLoader().load(annotation)] == [
        "A",
        "Q",
    ]
    assert next(row for row in result.audit_rows if row.protein_id == "B").exclusion_reason == (
        "missing_base_annotation"
    )
    metrics = {name: count for name, count, _ in result.coverage}
    assert metrics["formal_table_rows"] == "2"
    assert metrics["missing_base_annotation_predictions"] == "1"
    assert annotation.read_text(encoding="utf-8").splitlines()[0].split("\t") == list(
        ANNOTATION_COLUMNS
    )
    assert audit.read_text(encoding="utf-8").splitlines()[0].split("\t") == list(
        AUDIT_COLUMNS
    )
    for output in (annotation, audit, coverage, metadata):
        assert b"\r" not in output.read_bytes()
    formal_hash = hashlib.sha256(annotation.read_bytes()).hexdigest()
    assert formal_hash in audit.read_text(encoding="utf-8")
    assert formal_hash in metadata.read_text(encoding="utf-8")


def test_mapping_table_is_versioned_and_complete_for_archaea() -> None:
    assert PSORTB_MAPPING_RULE_VERSION.endswith("-v1")
    mapping = [
        (rule.raw_label, rule.formal_label, rule.status)
        for rule in PSORTB_ARCHAEA_MAPPING
    ]
    assert mapping == [
        ("Cytoplasmic", "cytosolic", "accepted"),
        ("CytoplasmicMembrane", "membrane", "accepted"),
        ("Cellwall", "", "unsupported"),
        ("Extracellular", "secreted", "accepted"),
        ("Unknown", "", "unknown"),
    ]


def test_non_archaeal_mode_is_rejected(tmp_path: Path) -> None:
    raw = _write(tmp_path / "raw.tsv", "Q\tCytoplasmic\t9\n")
    with pytest.raises(InputValidationError, match="organism_mode='archaea'"):
        convert_psortb(input_path=raw, organism_mode="positive")


def test_conflicting_existing_localization_is_not_overwritten(tmp_path: Path) -> None:
    raw = _write(tmp_path / "raw.tsv", "Q\tCytoplasmic\t9\n")
    base = _base_annotation(tmp_path / "base.tsv", ("Q",))
    lines = base.read_text(encoding="utf-8").splitlines()
    values = lines[1].split("\t")
    values[list(ANNOTATION_COLUMNS).index("localization_annotation")] = "membrane"
    _write(base, lines[0] + "\n" + "\t".join(values) + "\n")
    with pytest.raises(InputValidationError, match="Base annotation conflicts"):
        convert_psortb(
            input_path=raw,
            proteome_fasta=_fasta(tmp_path / "proteome.faa", ("Q",)),
            base_annotation_table=base,
        )
