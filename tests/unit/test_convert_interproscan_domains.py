import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from protein_interaction_hunter.adapters.local.domains import LocalDomainTsvLoader  # noqa: E402
from protein_interaction_hunter.exceptions import InputValidationError  # noqa: E402
from scripts.convert_interproscan_domains import (  # noqa: E402
    convert_interproscan,
    parse_interproscan_tsv,
    write_conversion_outputs,
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def _hit(
    protein_id: str = "Q",
    accession: str = "PF00001",
    description: str = "Alpha hydrolase family",
    start: int = 10,
    end: int = 150,
    extra: tuple[str, ...] = (
        "IPR000001",
        "Hydrolase superfamily",
        "GO:0016787(InterPro)",
        "MetaCyc:PWY-1",
    ),
) -> str:
    columns = (
        protein_id,
        "0123456789abcdef0123456789abcdef",
        "200",
        "Pfam",
        accession,
        description,
        str(start),
        str(end),
        "1.2E-20",
        "T",
        "26-07-2026",
        *extra,
    )
    return "\t".join(columns) + "\n"


def test_single_pfam_hit_preserves_query_and_audit_fields(tmp_path: Path) -> None:
    raw = _write(tmp_path / "query.tsv", "# InterProScan 5.78-109.0\n" + _hit())
    fasta = _write(tmp_path / "proteome.faa", ">Q alpha hydrolase\n" + "M" * 200 + "\n")

    result = convert_interproscan(
        input_path=raw,
        proteome_fasta=fasta,
        query_id="Q",
        interproscan_version="5.78-109.0",
    )

    assert len(result.domains) == 1
    domain = result.domains[0]
    assert (domain.protein_id, domain.source, domain.accession) == (
        "Q",
        "Pfam",
        "PF00001",
    )
    assert (domain.start, domain.end, domain.architecture_index) == (10, 150, 0)
    audit = result.audit_rows[0]
    assert audit.interpro_accession == "IPR000001"
    assert audit.go_terms == "GO:0016787(InterPro)"
    assert audit.pathway_terms == "MetaCyc:PWY-1"
    assert audit.interproscan_version == "5.78-109.0"
    assert dict((name, count) for name, count, _ in result.coverage)["query_domain_count"] == "1"


def test_nonoverlapping_repeat_is_retained_and_exact_duplicate_removed(
    tmp_path: Path,
) -> None:
    first = _hit(start=5, end=50)
    second = _hit(start=100, end=145)
    raw = _write(tmp_path / "repeats.tsv", second + first + first)

    result = convert_interproscan(input_path=raw)

    assert [(row.start, row.end) for row in result.domains] == [(5, 50), (100, 145)]
    assert [row.architecture_index for row in result.domains] == [0, 1]
    assert result.exact_duplicates_excluded == 1
    assert (
        dict((name, count) for name, count, _ in result.coverage)["repeated_domain_protein_count"]
        == "1"
    )


def test_empty_input_and_no_hit_protein_write_loadable_table(tmp_path: Path) -> None:
    raw = _write(tmp_path / "empty.tsv", "# no matches\n")
    fasta = _write(
        tmp_path / "proteome.faa",
        ">Q query\n" + "M" * 200 + "\n>N no hit\n" + "A" * 80 + "\n",
    )
    result = convert_interproscan(input_path=raw, proteome_fasta=fasta, query_id="Q")
    domain_output = tmp_path / "domains.tsv"
    audit_output = tmp_path / "audit.tsv"
    coverage_output = tmp_path / "coverage.tsv"

    write_conversion_outputs(
        result,
        domain_output=domain_output,
        audit_output=audit_output,
        coverage_output=coverage_output,
    )

    assert LocalDomainTsvLoader().load(domain_output) == []
    metrics = dict((name, count) for name, count, _ in result.coverage)
    assert metrics["total_proteins"] == "2"
    assert metrics["proteins_without_domain"] == "2"
    assert b"\r" not in domain_output.read_bytes()
    assert b"\r" not in audit_output.read_bytes()


def test_missing_interpro_fields_and_special_characters(tmp_path: Path) -> None:
    raw = _write(
        tmp_path / "base.tsv",
        _hit(
            description="Beta-propeller / alpha & omega",
            extra=(),
        ),
    )

    result = convert_interproscan(input_path=raw)

    assert result.domains[0].name == "Beta-propeller / alpha & omega"
    assert result.audit_rows[0].interpro_accession == ""
    assert result.audit_rows[0].go_terms == ""
    assert result.audit_rows[0].pathway_terms == ""


def test_unknown_protein_is_audited_and_excluded(tmp_path: Path) -> None:
    raw = _write(tmp_path / "unknown.tsv", _hit("Q") + _hit("UNKNOWN"))
    fasta = _write(tmp_path / "proteome.faa", ">Q query\n" + "M" * 200 + "\n")

    result = convert_interproscan(input_path=raw, proteome_fasta=fasta, query_id="Q")

    assert [row.protein_id for row in result.domains] == ["Q"]
    unknown = next(row for row in result.audit_rows if row.protein_id == "UNKNOWN")
    assert unknown.known_protein_id is False
    assert unknown.included_in_domain_table is False
    assert unknown.architecture_index == -1
    assert result.unknown_rows_excluded == 1


def test_order_independent_deterministic_output(tmp_path: Path) -> None:
    first = _write(
        tmp_path / "first.tsv",
        _hit("B", "PF2", start=90, end=120) + _hit("A", "PF1", start=1, end=50),
    )
    second = _write(
        tmp_path / "second.tsv",
        _hit("A", "PF1", start=1, end=50) + _hit("B", "PF2", start=90, end=120),
    )

    result_one = convert_interproscan(input_path=first)
    result_two = convert_interproscan(input_path=second)

    assert result_one.domains == result_two.domains
    assert [row.protein_id for row in result_one.domains] == ["A", "B"]


@pytest.mark.parametrize(
    "text, message",
    [
        ("\t".join(["x"] * 12) + "\n", "11, 13, 14, or 15 columns"),
        (_hit(start=151, end=150), "Invalid InterProScan coordinates"),
    ],
)
def test_malformed_input_is_rejected(tmp_path: Path, text: str, message: str) -> None:
    raw = _write(tmp_path / "malformed.tsv", text)

    with pytest.raises(InputValidationError, match=message):
        parse_interproscan_tsv(raw)
